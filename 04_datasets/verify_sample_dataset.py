# -*- coding: utf-8 -*-
"""Independent verification of GOT_monoamine_data_AUG2026.csv.

Re-reads every upstream source and checks the tidy sheet against it.  This does not
import build_sample_dataset.py - it re-derives the expected values so a bug in the build
cannot hide behind a shared helper.
"""
import csv
import os
import openpyxl
from collections import defaultdict

# Paths are relative to the repository root: the source workbooks are expected in
# data/raw/, the metadata workbooks in data/, and the tidy dataset is written to data/.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASETS = os.path.join(REPO_ROOT, "data")
RAW = os.path.join(DATASETS, "raw")
OUT_DIR = DATASETS
HEDDA = RAW + "/Monoamine concentrations Hedda GOT.xlsx"
ASLAUG = RAW + "/Monoamine concentrations_GOT_aslaug.xlsx"

CHOICE = sorted(list(range(31, 47)) + list(range(71, 79)) +
                list(range(87, 95)) + list(range(111, 127)))
REGIONS = ["Dl", "Vv", "POA", "DM", "RN"]
MD_REGIONS = {"Dl", "Vv", "POA"}
ANALYTE_COLS = ["HT5", "HIAA5", "HIAA5_HT5", "DA", "DOPAC", "DOPAC_DA"]

fail = []


def chk(cond, msg):
    print(("  OK   " if cond else "  FAIL ") + msg)
    if not cond:
        fail.append(msg)


def near(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= 1e-9 * max(abs(a), 1e-300)


def f(v):
    """CSV cell -> float or None."""
    if v is None or v == "NA" or v == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


rows = list(csv.DictReader(open(os.path.join(OUT_DIR, "GOT_monoamine_data_AUG2026.csv"),
                                encoding="utf-8")))
idx = {(int(r["sample_id"]), r["brain_region"]): r for r in rows}

print("\n[1] Shape")
chk(len(rows) == 240, "240 rows (48 fish x 5 regions), got %d" % len(rows))
chk(len(idx) == 240, "sample_id x brain_region is unique")
chk(sorted({int(r["sample_id"]) for r in rows}) == CHOICE,
    "the 48 fish are exactly 31-46, 71-78, 87-94, 111-126")
per = defaultdict(set)
for r in rows:
    per[int(r["sample_id"])].add(r["brain_region"])
chk(all(v == set(REGIONS) for v in per.values()), "every fish has all five regions")
chk(not any(int(r["sample_id"]) in range(47, 71) or int(r["sample_id"]) in range(79, 87)
            or int(r["sample_id"]) in range(95, 111) for r in rows),
    "no morphine/naloxone fish leaked in")

print("\n[2] DM/RN analytes match the master students' Concentrations sheets")
conc = {}
for path, sheet, reg in ((HEDDA, "Concentrations Dm", "DM"), (HEDDA, "Concentrations Rh", "RN"),
                         (ASLAUG, "Concentrations Dm", "DM"), (ASLAUG, "Concentrations Rn", "RN")):
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    for i in range(2, ws.max_row + 1):
        sid = ws.cell(i, 26).value
        if not isinstance(sid, (int, float)):
            continue
        g = lambda c: ws.cell(i, c).value if isinstance(ws.cell(i, c).value, (int, float)) else None
        conc[(int(sid), reg)] = dict(DOPAC=g(28), HIAA5=g(29), DA=g(30), HT5=g(31),
                                     DOPAC_DA=g(32), HIAA5_HT5=g(33))
errs = []
for sid in CHOICE:
    for reg in ("DM", "RN"):
        r = idx[(sid, reg)]
        src = conc.get((sid, reg))
        if src is None:
            if any(f(r[c]) is not None for c in ANALYTE_COLS):
                errs.append("%s %s has analytes but no source" % (sid, reg))
            continue
        for c in ANALYTE_COLS:
            if not near(f(r[c]), src[c]):
                errs.append("%s %s %s: %r vs source %r" % (sid, reg, c, r[c], src[c]))
chk(not errs, "every DM/RN analyte matches: %s" % (errs[:3] or "all %d values" % (95 * 6)))

print("\n[3] The DM 5-HT repair actually happened")
mono = {}
ws = openpyxl.load_workbook(DATASETS + "/monoamine concentration ALL GOT.xlsx",
                            data_only=True)["sheet1"]
for i in range(2, ws.max_row + 1):
    sid = ws.cell(i, 1).value
    if sid is None:
        continue
    mono[(int(sid), ws.cell(i, 2).value)] = ws.cell(i, 10).value
differ = [sid for sid in CHOICE
          if (sid, "DM") in mono and f(idx[(sid, "DM")]["HT5"]) is not None
          and not near(f(idx[(sid, "DM")]["HT5"]), mono[(sid, "DM")])]
chk(len(differ) == 43, "43 DM rows carry a 5-HT that differs from the old master file, got %d" % len(differ))
bad_ratio = []
for sid in CHOICE:
    for reg in ("DM", "RN"):
        r = idx[(sid, reg)]
        ht, hiaa, rat = f(r["HT5"]), f(r["HIAA5"]), f(r["HIAA5_HT5"])
        if None in (ht, hiaa, rat) or ht == 0:
            continue
        if abs(hiaa / ht - rat) > 1e-9 * max(abs(rat), 1):
            bad_ratio.append((sid, reg))
chk(not bad_ratio, "5-HIAA/5-HT is now internally consistent in every DM/RN row: %s"
    % (bad_ratio[:4] or "yes"))
drag = [sid for sid in CHOICE if f(idx[(sid, "DM")]["HT5"]) is not None
        and abs(round(f(idx[(sid, "DM")]["HT5"]) % 1, 10) - 0.0576637584) < 1e-9]
chk(not drag, "no fill-drag value survives: %s" % (drag or "none"))

print("\n[4] Punches and protein come from the Bradford sheets, not the shifted aggregate")
brad = {}
for path, sheet, reg, cs, cm, cp in ((HEDDA, "Bradford Dm", "DM", 1, 3, 4),
                                     (HEDDA, "Bradford Rh", "RN", 1, 3, 4),
                                     (ASLAUG, "Bradford Dm", "DM", 1, 2, 3),
                                     (ASLAUG, "Bradford Rn", "RN", 1, 2, 3)):
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    for i in range(2, ws.max_row + 1):
        sid = ws.cell(i, cs).value
        if not isinstance(sid, (int, float)):
            continue
        g = lambda c: ws.cell(i, c).value if isinstance(ws.cell(i, c).value, (int, float)) else None
        brad[(int(sid), reg)] = (g(cm), g(cp))
errs = []
for sid in CHOICE:
    for reg in ("DM", "RN"):
        r = idx[(sid, reg)]
        src = brad.get((sid, reg))
        if src is None:
            if f(r["punches"]) is not None or f(r["protein_mg"]) is not None:
                errs.append("%s %s has tissue data but no Bradford source" % (sid, reg))
            continue
        if not near(f(r["protein_mg"]), src[0]):
            errs.append("%s %s protein %r vs %r" % (sid, reg, r["protein_mg"], src[0]))
        if not near(f(r["punches"]), src[1]):
            errs.append("%s %s punches %r vs %r" % (sid, reg, r["punches"], src[1]))
chk(not errs, "every DM/RN punch and protein matches the Bradford sheets: %s" % (errs[:3] or "yes"))

print("\n[5] Dl/Vv/POA punches and protein match the combined workbook")
ws = openpyxl.load_workbook(ONEDRIVE + "/MD_punches_x_Bradford_protein_COMBINED.xlsx",
                            data_only=True)["Combined"]
H = {ws.cell(4, j).value: j for j in range(1, ws.max_column + 1) if ws.cell(4, j).value}
errs = []
for i in range(5, ws.max_row + 1):
    sid, reg = ws.cell(i, H["Sample ID"]).value, ws.cell(i, H["Brain region"]).value
    if sid is None or int(sid) not in CHOICE:
        continue
    r = idx[(int(sid), reg)]
    for label, col in (("punches", "Punches"), ("protein_mg", "Protein mg")):
        want = ws.cell(i, H[col]).value
        want = want if isinstance(want, (int, float)) else None
        if not near(f(r[label]), want):
            errs.append("%s %s %s: %r vs %r" % (sid, reg, label, r[label], want))
chk(not errs, "every Dl/Vv/POA punch and protein matches: %s" % (errs[:3] or "all 144 rows"))

print("\n[6] Analyte block for Dl/Vv/POA is empty, not zero")
leak = [(r["sample_id"], r["brain_region"], c) for r in rows if r["brain_region"] in MD_REGIONS
        for c in ANALYTE_COLS if r[c] != "NA"]
chk(not leak, "all 144 Dl/Vv/POA rows x 6 analytes are NA: %s" % (leak[:3] or "yes"))
chk(all(r["analyte_status"] == "not quantified yet" for r in rows
        if r["brain_region"] in MD_REGIONS), "all Dl/Vv/POA rows say 'not quantified yet'")

print("\n[7] Operator")
chk(all(r["operator"] == "Martina" for r in rows if r["brain_region"] in MD_REGIONS),
    "operator is Martina on every Dl/Vv/POA row")
ops = {r["operator"] for r in rows if r["brain_region"] in ("DM", "RN")}
chk(ops <= {"Hedda", "Aslaug", "NA"}, "DM/RN operators are Hedda/Aslaug: %s" % sorted(ops))

print("\n[8] Per-fish metadata is constant across a fish's five rows")
META = ["group", "treatment", "trial", "sampling_date", "sex", "weight_g", "length_cm", "tank",
        "hematocrit", "hematocrit2", "haemoglobin", "acid_inj_vol_ml", "second_inj_vol_ml",
        "brain_analysis", "plasma_cortisol", "glucose_mmol_L", "lactate_mmol_L", "sampling_notes"]
errs = []
for sid in CHOICE:
    for c in META:
        vals = {idx[(sid, reg)][c] for reg in REGIONS}
        if len(vals) > 1:
            errs.append("%s %s varies: %s" % (sid, c, vals))
chk(not errs, "all 18 metadata fields constant per fish: %s" % (errs[:3] or "yes"))

print("\n[9] Metadata values against sampling data GOT.xlsx")
ws = openpyxl.load_workbook(DATASETS + "/sampling data GOT.xlsx", data_only=True)["sampling data"]
errs = []
for i in range(2, ws.max_row + 1):
    sid = ws.cell(i, 1).value
    if not isinstance(sid, (int, float)) or int(sid) not in CHOICE:
        continue
    r = idx[(int(sid), "Dl")]
    d = ws.cell(i, 2).value
    if r["sampling_date"] != d.date().isoformat():
        errs.append("%s date %r vs %r" % (sid, r["sampling_date"], d))
    if r["sex"] != str(ws.cell(i, 13).value):
        errs.append("%s sex %r vs %r" % (sid, r["sex"], ws.cell(i, 13).value))
    if not near(f(r["weight_g"]), ws.cell(i, 5).value):
        errs.append("%s weight %r" % (sid, r["weight_g"]))
chk(not errs, "date, sex and weight match the master for all 48 fish: %s" % (errs[:3] or "yes"))
chk(all(r["sampling_date"].startswith("2024-11") for r in rows),
    "every sampling date is in Nov 2024 (the 2025 typos are gone)")
for sid, want in ((71, "2024-11-20"), (87, "2024-11-22"), (119, "2024-11-27")):
    chk(idx[(sid, "Dl")]["sampling_date"] == want,
        "sample %d date is %s, not the monoamine file's wrong value" % (sid, want))

print("\n[10] Canonical vocabularies")
chk({r["group"] for r in rows} == {"exe. Choice", "no choice"},
    "group values: %s" % sorted({r["group"] for r in rows}))
chk({r["treatment"] for r in rows} == {"exercise", "control"},
    "treatment values: %s" % sorted({r["treatment"] for r in rows}))
chk({r["trial"] for r in rows} == {"1", "2", "3"}, "trial values: %s" % sorted({r["trial"] for r in rows}))
errs = [sid for sid in CHOICE
        if idx[(sid, "Dl")]["trial"] != str(1 if sid <= 62 else (2 if sid <= 94 else 3))]
chk(not errs, "trial matches the sample-ID block rule for all 48: %s" % (errs or "yes"))

print("\n[11] Derived columns")
errs = []
for r in rows:
    mg, ug, pun, per = f(r["protein_mg"]), f(r["protein_ug"]), f(r["punches"]), f(r["ug_protein_per_punch"])
    if mg is None:
        if ug is not None or per is not None:
            errs.append("%s %s derived without protein" % (r["sample_id"], r["brain_region"]))
        continue
    if not near(ug, mg * 1000):
        errs.append("%s %s ug" % (r["sample_id"], r["brain_region"]))
    if pun:
        if not near(per, ug / pun):
            errs.append("%s %s per-punch" % (r["sample_id"], r["brain_region"]))
    elif per is not None:
        errs.append("%s %s per-punch without punches" % (r["sample_id"], r["brain_region"]))
chk(not errs, "protein_ug and ug_protein_per_punch are correct everywhere: %s" % (errs[:3] or "yes"))

print("\n[12] Known gaps, and only those")
gaps = sorted((int(r["sample_id"]), r["brain_region"]) for r in rows if r["protein_mg"] == "NA")
expected_no_record = [(35, "DM"), (37, "DM"), (36, "POA"), (72, "Dl"), (74, "Vv")]
chk(all(g in gaps for g in expected_no_record),
    "the five known missing combinations are present as gaps")
chk(idx[(37, "DM")]["analyte_status"] == "no HPLC record", "37 DM: no HPLC record")
chk(idx[(35, "DM")]["analyte_status"] == "not computable (protein #NUM!)"
    and idx[(35, "DM")]["punches"] == "4",
    "35 DM: Bradford failed (#NUM!) so analytes are not computable, but 4 punches are recorded")
chk(sum(1 for r in rows if r["protein_mg"] == "NA") == 27,
    "27 of 240 rows have no protein, got %d" % sum(1 for r in rows if r["protein_mg"] == "NA"))

print("\n[13] CSV is clean for R")
hdr = list(rows[0].keys())
chk(len(hdr) == len(set(hdr)), "no duplicate column names")
chk(all(h == h.strip() and " " not in h for h in hdr), "no spaces or padding in header names")
chk(all(len(r) == len(hdr) for r in rows), "every row has the full column count")
chk(not any((r[c] or "").strip() != (r[c] or "") for r in rows for c in hdr),
    "no leading/trailing whitespace in any cell")

print("\n" + ("ALL CHECKS PASSED" if not fail else "FAILURES: %d" % len(fail)))
for x in fail:
    print("  -", x)
raise SystemExit(1 if fail else 0)
