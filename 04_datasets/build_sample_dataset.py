# -*- coding: utf-8 -*-
"""Build the tidy GOT choice-experiment datasheet: one row per sample x brain region.

Scope: the 48 exercise/no-choice fish (31-46, 71-78, 87-94, 111-126) x 5 brain regions
(Dl, Vv, POA, DM, RN) = 240 rows.  The morphine/naloxone arm is deliberately excluded.

Analyte values for Dl, Vv and POA are left EMPTY here; they are added by
merge_dl_vv_poa_concentrations.py.
Punches and protein for those three regions are already known and are filled in.

SOURCE PRECEDENCE:

  DM / RN analytes + protein + punches
      the two master students' workbooks, sheets 'Concentrations *' and 'Bradford *'.
      These are the authoritative chain (peak areas -> standards -> concentrations).
      NOT 'monoamine concentration ALL GOT.xlsx', whose DM 5-HT column is wrong in 88 of
      92 rows, and NOT its aggregated 'punches' column, which is shifted by one row.
  Dl / Vv / POA punches + protein
      MD_punches_x_Bradford_protein_COMBINED.xlsx, which carries the seven
      paper-confirmed plate-label corrections.
  per-fish metadata
      sampling data GOT.xlsx (master: real datetimes, sex, weight, length, tank, notes)
      + GOT glu_lac_cort updated march 2025.xlsx (cortisol, glucose, lactate)
      + trial from monoamine concentration ALL GOT.xlsx (the only source that has it)

Reads only; writes CSV + XLSX to data/.  Deterministic and re-runnable.
"""
import csv
import os
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# Paths are relative to the repository root: the source workbooks are expected in
# data/raw/, the metadata workbooks in data/, and the tidy dataset is written to data/.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASETS = os.path.join(REPO_ROOT, "data")
RAW = os.path.join(DATASETS, "raw")
OUT_DIR = DATASETS

SRC = {
    "hedda":    RAW + "/Monoamine concentrations Hedda GOT.xlsx",
    "aslaug":   RAW + "/Monoamine concentrations_GOT_aslaug.xlsx",
    "combined": RAW + "/MD_punches_x_Bradford_protein_COMBINED.xlsx",
    "sampling": DATASETS + "/sampling data GOT.xlsx",
    "glc":      DATASETS + "/GOT glu_lac_cort updated march 2025.xlsx",
    "mono":     DATASETS + "/monoamine concentration ALL GOT.xlsx",
}

CHOICE = sorted(list(range(31, 47)) + list(range(71, 79)) +
                list(range(87, 95)) + list(range(111, 127)))
REGIONS = ["Dl", "Vv", "POA", "DM", "RN"]
MD_REGIONS = {"Dl", "Vv", "POA"}
HPLC_REGIONS = {"DM", "RN"}

# canonical vocabularies
GROUP_MAP = {"exe. choice": "exe. Choice", "exercise choice": "exe. Choice",
             "no choice": "no choice"}


def load(key, sheet):
    return openpyxl.load_workbook(SRC[key], data_only=True)[sheet]


def num(v):
    return v if isinstance(v, (int, float)) else None


# ---------------------------------------------------------------------------
# 1. DM / RN  -  analytes from 'Concentrations *', punches+protein from 'Bradford *'
# ---------------------------------------------------------------------------
# 'Concentrations' block starts at column Z: Sample, NE, DOPAC, 5-HIAA, DA, 5-HT,
# DOPAC/DA, 5HIAA/5HT  (columns 26..33)
CONC_SHEETS = [("hedda", "Concentrations Dm", "DM", "Hedda"),
               ("hedda", "Concentrations Rh", "RN", "Hedda"),
               ("aslaug", "Concentrations Dm", "DM", "Aslaug"),
               ("aslaug", "Concentrations Rn", "RN", "Aslaug")]

analytes = {}       # (sid, region) -> dict
for key, sheet, region, operator in CONC_SHEETS:
    ws = load(key, sheet)
    for i in range(2, ws.max_row + 1):
        sid = ws.cell(i, 26).value
        if not isinstance(sid, (int, float)):
            continue
        analytes[(int(sid), region)] = dict(
            operator=operator,
            DOPAC=num(ws.cell(i, 28).value), HIAA5=num(ws.cell(i, 29).value),
            DA=num(ws.cell(i, 30).value), HT5=num(ws.cell(i, 31).value),
            DOPAC_DA=num(ws.cell(i, 32).value), HIAA5_HT5=num(ws.cell(i, 33).value))

# Bradford sheets differ in layout between the two students
BRAD_SHEETS = [("hedda", "Bradford Dm", "DM", 1, 3, 4, 5),
               ("hedda", "Bradford Rh", "RN", 1, 3, 4, 5),
               ("aslaug", "Bradford Dm", "DM", 1, 2, 3, None),
               ("aslaug", "Bradford Rn", "RN", 1, 2, 3, None)]

tissue = {}         # (sid, region) -> dict
for key, sheet, region, c_sid, c_mg, c_pun, c_note in BRAD_SHEETS:
    ws = load(key, sheet)
    for i in range(2, ws.max_row + 1):
        sid = ws.cell(i, c_sid).value
        if not isinstance(sid, (int, float)):
            continue
        tissue[(int(sid), region)] = dict(
            protein_mg=num(ws.cell(i, c_mg).value),
            punches=num(ws.cell(i, c_pun).value),
            notes=ws.cell(i, c_note).value if c_note else None)

# ---------------------------------------------------------------------------
# 2. Dl / Vv / POA  -  punches, protein, plate provenance from the combined workbook
# ---------------------------------------------------------------------------
ws = load("combined", "Combined")
HDR = {ws.cell(4, j).value: j for j in range(1, ws.max_column + 1) if ws.cell(4, j).value}
md_rows = {}
for i in range(5, ws.max_row + 1):
    sid = ws.cell(i, HDR["Sample ID"]).value
    reg = ws.cell(i, HDR["Brain region"]).value
    if sid is None:
        continue
    md_rows[(int(sid), reg)] = dict(
        punches=num(ws.cell(i, HDR["Punches"]).value),
        protein_mg=num(ws.cell(i, HDR["Protein mg"]).value),
        status=ws.cell(i, HDR["Measurement status"]).value,
        plate=ws.cell(i, HDR["Plate"]).value,
        well=ws.cell(i, HDR["Well"]).value,
        notes=ws.cell(i, HDR["Bradford note"]).value,
        md_date=ws.cell(i, HDR["MD date"]).value,
        md_notes=ws.cell(i, HDR["MD notes"]).value,
        quality=ws.cell(i, HDR["Sample quality flag"]).value)

# ---------------------------------------------------------------------------
# 3. per-fish metadata
# ---------------------------------------------------------------------------
ws = load("sampling", "sampling data")
meta = {}
for i in range(2, ws.max_row + 1):
    sid = ws.cell(i, 1).value
    if not isinstance(sid, (int, float)):
        continue
    d = ws.cell(i, 2).value
    meta[int(sid)] = dict(
        sampling_date=d.date().isoformat() if hasattr(d, "date") else d,
        tank=ws.cell(i, 3).value, weight_g=num(ws.cell(i, 5).value),
        length_cm=num(ws.cell(i, 6).value),
        group=GROUP_MAP.get(str(ws.cell(i, 7).value).strip().lower(), ws.cell(i, 7).value),
        acid_inj_vol_ml=num(ws.cell(i, 8).value), second_inj_vol_ml=num(ws.cell(i, 9).value),
        hematocrit=num(ws.cell(i, 10).value), hematocrit2=num(ws.cell(i, 11).value),
        haemoglobin=num(ws.cell(i, 12).value), sex=ws.cell(i, 13).value,
        sampling_notes=ws.cell(i, 14).value, brain_analysis=ws.cell(i, 15).value)

ws = load("glc", "results")
for i in range(2, ws.max_row + 1):
    sid = ws.cell(i, 1).value
    if not isinstance(sid, (int, float)) or int(sid) not in meta:
        continue
    meta[int(sid)].update(plasma_cortisol=num(ws.cell(i, 5).value),
                          glucose_mmol_L=num(ws.cell(i, 6).value),
                          lactate_mmol_L=num(ws.cell(i, 8).value))

# treatment and trial: the monoamine file is the only source of trial
ws = load("mono", "sheet1")
for i in range(2, ws.max_row + 1):
    sid = ws.cell(i, 1).value
    if not isinstance(sid, (int, float)) or int(sid) not in meta:
        continue
    m = meta[int(sid)]
    m.setdefault("treatment", str(ws.cell(i, 5).value).strip().lower())
    m.setdefault("trial", num(ws.cell(i, 7).value))

# trial is also implied by the sampling blocks - verify rather than trust
for sid, m in meta.items():
    if sid not in CHOICE:
        continue
    implied = 1 if sid <= 62 else (2 if sid <= 94 else 3)
    if m.get("trial") != implied:
        raise SystemExit("trial mismatch for sample %s: file says %r, sample-ID block implies %r"
                         % (sid, m.get("trial"), implied))

# ---------------------------------------------------------------------------
# 4. assemble
# ---------------------------------------------------------------------------
COLUMNS = ["sample_id", "brain_region", "operator",
           "punches", "protein_mg", "protein_ug", "ug_protein_per_punch",
           "HT5", "HIAA5", "HIAA5_HT5", "DA", "DOPAC", "DOPAC_DA",
           "analyte_status", "tissue_status", "bradford_plate", "bradford_well", "tissue_notes",
           "group", "treatment", "trial", "sampling_date", "sex", "weight_g", "length_cm",
           "tank", "hematocrit", "hematocrit2", "haemoglobin",
           "acid_inj_vol_ml", "second_inj_vol_ml", "brain_analysis",
           "plasma_cortisol", "glucose_mmol_L", "lactate_mmol_L", "sampling_notes",
           "md_date", "md_notes", "sample_quality_flag"]

ANALYTE_COLS = ["HT5", "HIAA5", "HIAA5_HT5", "DA", "DOPAC", "DOPAC_DA"]

rows = []
for sid in CHOICE:
    m = meta[sid]
    for reg in REGIONS:
        r = {c: None for c in COLUMNS}
        r["sample_id"] = sid
        r["brain_region"] = reg
        for k in ("group", "treatment", "trial", "sampling_date", "sex", "weight_g", "length_cm",
                  "tank", "hematocrit", "hematocrit2", "haemoglobin", "acid_inj_vol_ml",
                  "second_inj_vol_ml", "brain_analysis", "plasma_cortisol", "glucose_mmol_L",
                  "lactate_mmol_L", "sampling_notes"):
            r[k] = m.get(k)

        if reg in MD_REGIONS:
            r["operator"] = "Martina"
            t = md_rows.get((sid, reg))
            if t:
                r["punches"] = t["punches"]
                r["protein_mg"] = t["protein_mg"]
                r["bradford_plate"] = t["plate"]
                r["bradford_well"] = t["well"]
                r["tissue_notes"] = t["notes"]
                r["tissue_status"] = t["status"]
                r["md_date"] = t["md_date"]
                r["md_notes"] = t["md_notes"]
                r["sample_quality_flag"] = t["quality"]
            r["analyte_status"] = "not quantified yet"
        else:
            a = analytes.get((sid, reg))
            t = tissue.get((sid, reg))
            if a:
                r["operator"] = a["operator"]
                for c in ANALYTE_COLS:
                    r[c] = a[c]
                # a row can exist yet hold #NUM! throughout when its Bradford protein failed,
                # because every concentration is divided by protein
                r["analyte_status"] = ("OK" if any(a[c] is not None for c in ANALYTE_COLS)
                                       else "not computable (protein #NUM!)")
            else:
                r["analyte_status"] = "no HPLC record"
            if t:
                r["punches"] = t["punches"]
                r["protein_mg"] = t["protein_mg"]
                r["tissue_notes"] = t["notes"]
                r["tissue_status"] = ("OK" if t["protein_mg"] is not None
                                      else "protein not computable (#NUM!)")
            else:
                r["tissue_status"] = "no Bradford record"

        if r["protein_mg"] is not None:
            r["protein_ug"] = r["protein_mg"] * 1000
            if isinstance(r["punches"], (int, float)) and r["punches"]:
                r["ug_protein_per_punch"] = r["protein_ug"] / r["punches"]
        rows.append(r)

assert len(rows) == len(CHOICE) * len(REGIONS) == 240, len(rows)

# ---------------------------------------------------------------------------
# 5. write CSV
# ---------------------------------------------------------------------------
csv_path = os.path.join(OUT_DIR, "GOT_monoamine_data_AUG2026.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(COLUMNS)
    for r in rows:
        w.writerow(["NA" if r[c] is None else r[c] for c in COLUMNS])

# ---------------------------------------------------------------------------
# 6. write XLSX
# ---------------------------------------------------------------------------
AR = "Arial"
HDRF = Font(name=AR, sz=10, bold=True, color="FFFFFFFF")
HFIL = PatternFill("solid", fgColor="FF4472C4")
AFIL = PatternFill("solid", fgColor="FF7030A0")      # analyte header
INPUT = PatternFill("solid", fgColor="FFFFF2CC")     # cells you type into
BODY = Font(name=AR, sz=10)
THIN = Side(style="thin", color="FFD9D9D9")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

wb = openpyxl.Workbook()
sh = wb.active
sh.title = "tidy"
sh.cell(1, 1, "GOT choice experiment - one row per sample x brain region").font = Font(name=AR, sz=12, bold=True)
sh.cell(2, 1, "48 fish x 5 regions = 240 rows. Yellow cells are the analyte block for Dl/Vv/POA, "
              "still to be quantified.").font = Font(name=AR, sz=10, italic=True, color="FF595959")
sh.cell(3, 1, "CAUTION: DM/RN punch counts (4-54) are on a different scale from Dl/Vv/POA (3-14), so "
              "ug_protein_per_punch is comparable within a region but NOT across the two sets. "
              "See 04_datasets/build_sample_dataset.py.").font = Font(name=AR, sz=10, bold=True, color="FFC00000")

for j, c in enumerate(COLUMNS, 1):
    cell = sh.cell(4, j, c)
    cell.font, cell.border = HDRF, BOX
    cell.fill = AFIL if c in ANALYTE_COLS else HFIL
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    sh.column_dimensions[get_column_letter(j)].width = max(9, min(20, len(c) + 3))
sh.row_dimensions[4].height = 30

for i, r in enumerate(rows, 5):
    pending = r["brain_region"] in MD_REGIONS
    for j, c in enumerate(COLUMNS, 1):
        cell = sh.cell(i, j, r[c])
        cell.font, cell.border = BODY, BOX
        if c in ("sample_id", "brain_region", "punches", "trial", "bradford_plate", "bradford_well"):
            cell.alignment = Alignment(horizontal="center")
        if c == "protein_mg":
            cell.number_format = "0.000000"
        if c in ("protein_ug", "ug_protein_per_punch"):
            cell.number_format = "0.000"
        if c in ANALYTE_COLS:
            cell.number_format = "0.0000"
            if pending:
                cell.fill = INPUT
sh.freeze_panes = "D5"
sh.auto_filter.ref = "A4:%s%d" % (get_column_letter(len(COLUMNS)), 4 + len(rows))
sh.sheet_view.showGridLines = False
xlsx_path = os.path.join(OUT_DIR, "GOT_monoamine_data_AUG2026.xlsx")
wb.save(xlsx_path)

# ---------------------------------------------------------------------------
print("wrote %s" % csv_path)
print("wrote %s" % xlsx_path)
filled = sum(1 for r in rows if r["HT5"] is not None)
pun = sum(1 for r in rows if r["punches"] is not None)
pro = sum(1 for r in rows if r["protein_mg"] is not None)
print("rows=%d  punches filled=%d  protein filled=%d  analytes filled=%d (DM/RN) "
      "pending=%d (Dl/Vv/POA)" % (len(rows), pun, pro, filled, 48 * 3))
gaps = [(r["sample_id"], r["brain_region"], r["tissue_status"], r["analyte_status"])
        for r in rows if r["protein_mg"] is None]
print("rows with no protein: %d" % len(gaps))
for g in gaps:
    print("   %-4s %-4s %-26s %s" % g)
