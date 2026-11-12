"""Merge the Dl / Vv / POA monoamine concentrations (July 2026 HPLC runs) into a
copy of the tidy monoamine dataset, so the Rosner monoamine script can read all
five regions from one file.

Reads (never writes to either):
  <datasets>/GOT_monoamine_data_AUG2026.csv       tidy dataset, Dl/Vv/POA analytes empty
  <concentrations.csv>                            HPLC peak areas -> concentrations

Writes:
  data/GOT_monoamine_data_with_Dl_Vv_POA.csv   (input of 05_statistics/plasma_monoamine_models.R)
  data/merge_log_Dl_Vv_POA.md

Usage: python 04_datasets/merge_dl_vv_poa_concentrations.py [<tidy csv> <concentrations csv>]

Rules (numbers only: status, protein notes and LOW_SN flags are carried along in
the log but do not exclude anything):
  * the numeric columns 5-HT, 5-HIAA, 5-HIAA/5-HT, DA, DOPAC, DOPAC/DA are taken
    as given; an empty cell stays missing (peak not quantified);
  * chromatogram file "127 Dl" is sample 122 Dl (no fish 127; noted in the csv);
  * "81 POA ASL TEST" is not one of the 48 fish and is left out;
  * 92 Vv was injected twice: concentrations are the mean of the runs that have a
    value, and the ratios are recomputed from those means.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "data"
DEFAULT_TIDY = OUT_DIR / "GOT_monoamine_data_AUG2026.csv"
DEFAULT_CONC = OUT_DIR / "raw" / "Dl_Vv_POA_concentrations.csv"   # HPLC peak areas -> concentrations

AREA_TO_REGION = {"DL": "Dl", "VV": "Vv", "POA": "POA"}
CONC_COLS = {"5-HT": "HT5", "5-HIAA": "HIAA5", "DA": "DA", "DOPAC": "DOPAC"}
RATIO_COLS = {"5-HIAA/5-HT": "HIAA5_HT5", "DOPAC/DA": "DOPAC_DA"}
ID_REMAP = {("127", "DL"): "122"}


def main(tidy_path=DEFAULT_TIDY, conc_path=DEFAULT_CONC):
    tidy = pd.read_csv(tidy_path)
    conc = pd.read_csv(conc_path, dtype={"SAMPLE_ID": str})
    log = []

    conc["run2"] = conc["SAMPLE_ID"].str.contains("run 2", na=False)
    conc["sid"] = conc["SAMPLE_ID"].str.replace(r"\s*\(run 2\)", "", regex=True).str.strip()
    for (old, area), new in ID_REMAP.items():
        m = (conc["sid"] == old) & (conc["AREA"] == area)
        if m.any():
            conc.loc[m, "sid"] = new
            log.append(f"- {old} {area} taken as sample {new} (as noted in concentrations.csv).")
    conc["sid"] = conc["sid"].astype(int)

    fish = set(tidy["sample_id"].unique())
    extra = conc[~conc["sid"].isin(fish)]
    for _, r in extra.iterrows():
        log.append(f"- {r['SAMPLE_ID']} {r['AREA']} is not one of the 48 fish: left out ({r['file']}).")
    conc = conc[conc["sid"].isin(fish)].copy()

    rows = []
    for (sid, area), g in conc.groupby(["sid", "AREA"]):
        rec = {"sample_id": sid, "brain_region": AREA_TO_REGION[area]}
        if len(g) > 1:
            for src, dst in CONC_COLS.items():
                rec[dst] = g[src].mean(skipna=True)  # NaN if no run has a value
            rec["HIAA5_HT5"] = rec["HIAA5"] / rec["HT5"]
            rec["DOPAC_DA"] = rec["DOPAC"] / rec["DA"]
            log.append(
                f"- {sid} {area}: {len(g)} injections; concentrations averaged over the runs "
                f"with a value, ratios recomputed from the means."
            )
        else:
            r = g.iloc[0]
            for src, dst in {**CONC_COLS, **RATIO_COLS}.items():
                rec[dst] = r[src]
        rows.append(rec)
    new = pd.DataFrame(rows)

    out = tidy.copy()
    key = ["sample_id", "brain_region"]
    out = out.merge(new, on=key, how="left", suffixes=("", "_new"))
    is_new_region = out["brain_region"].isin(AREA_TO_REGION.values())
    for dst in list(CONC_COLS.values()) + list(RATIO_COLS.values()):
        out.loc[is_new_region, dst] = out.loc[is_new_region, dst + "_new"]
        out = out.drop(columns=dst + "_new")
    analyte_cols = list(CONC_COLS.values())
    has_any = out[analyte_cols].notna().any(axis=1)
    in_conc = out.set_index(key).index.isin(new.set_index(key).index)
    out.loc[is_new_region & has_any, "analyte_status"] = "OK"
    out.loc[is_new_region & ~has_any & in_conc, "analyte_status"] = "run, no peak quantified"
    out.loc[is_new_region & ~in_conc, "analyte_status"] = "no HPLC record"

    OUT_DIR.mkdir(exist_ok=True)
    out.to_csv(OUT_DIR / "GOT_monoamine_data_with_Dl_Vv_POA.csv", index=False, na_rep="NA")

    counts = (
        out[is_new_region]
        .groupby("brain_region")[analyte_cols + list(RATIO_COLS.values())]
        .count()
    )
    status = out[is_new_region].groupby(["brain_region", "analyte_status"]).size().unstack(fill_value=0)
    missing = out[is_new_region & ~in_conc][key]

    # numbers-only bookkeeping: how many values rest on imputed / extrapolated protein
    c2 = conc.copy()
    c2["imputed"] = c2["status"].str.contains("imputed", case=False, na=False)
    c2["extrap"] = c2["protein_note"].fillna("").str.contains("extrapolated")
    c2["has_value"] = c2[list(CONC_COLS)].notna().any(axis=1)
    prot = c2[c2["has_value"] & ~c2["run2"]].groupby("AREA")[["imputed", "extrap"]].sum()
    prot["samples_with_values"] = c2[c2["has_value"] & ~c2["run2"]].groupby("AREA").size()

    md = ["# Merge of Dl / Vv / POA concentrations into the tidy monoamine dataset", "",
          f"Source: `{conc_path}`", "", "## Decisions", ""] + log + [
          "", "## Non-missing values per region (of 48 fish)", "", counts.to_markdown(), "",
          "## Analyte status", "", status.to_markdown(), "",
          "## Fish with no chromatogram in the file", "",
          missing.to_markdown(index=False) if len(missing) else "none", "",
          "## Protein basis of the samples that have at least one value", "",
          "Not used to exclude anything; given so the concentrations can be read with it in mind.",
          "", prot.to_markdown(), ""]
    (OUT_DIR / "merge_log_Dl_Vv_POA.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main(*[Path(a) for a in sys.argv[1:3]])
