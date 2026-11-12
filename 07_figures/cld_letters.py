"""Compact letter display (CLD) helpers for the manuscript figures.
Two conditions must hold for a panel to carry letters:
 1. the FULL interaction among the factors shown in the panel is significant (treatment x group
    for a four-cell panel, treatment x group x interval for an eight-point interval panel,
    treatment x group x zone for panels split by zone) -- a significant main effect or
    lower-order interaction alone gives no letters (see statement_p / anova_table_p);
 2. the cells differ: every letter function returns None when all cells share the same letters."""
import os, re, string
import pandas as pd

ALPHA = 0.05


def statement_p(path, key, key_col="variable"):
    """p-value of one pre-formatted pipeline statement ('... p = 0.131 ...' / '... p < 0.001 ...')
    for `key` in an anova_stmt_*.csv. Returns None if the file, key or p-value is missing."""
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, encoding="utf-8-sig")
    if key_col not in d.columns or "statement" not in d.columns:
        return None
    s = d.loc[d[key_col] == key, "statement"]
    m = re.search(r"p\s*[<=]\s*([0-9.]+)", str(s.iloc[0])) if len(s) else None
    return float(m.group(1)) if m else None


def anova_table_p(path, term, model=None, p_col="Pr(>Chisq)"):
    """p-value of `term` in a raw ANOVA table (columns term, Pr(>Chisq)[, model])."""
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path)
    if model is not None and "model" in d.columns:
        d = d[d.model == model]
    p = d.loc[d.term == term, p_col]
    return float(p.iloc[0]) if len(p) == 1 else None


def significant(p):
    """True only for a known p-value below ALPHA (missing information = no letters)."""
    return p is not None and p < ALPHA


def canon_group(x):
    return "No choice" if str(x).strip().lower().startswith("no") else "Exercise choice"


def letters_from_pvalues(cells, pvals, alpha=0.05):
    """cells: ordered list of cell keys (plot order). pvals: {(cell_i, cell_j): adjusted p}.
    Insert-and-absorb algorithm; cells sharing a letter do not differ. Letters are assigned in
    plot order, so the leftmost cell always carries "a". Returns {cell: letters} or None."""
    idx = {c: i for i, c in enumerate(cells)}
    sig = [(idx[a], idx[b]) for (a, b), p in pvals.items() if p < alpha]
    if not sig:
        return None
    sets = [set(range(len(cells)))]
    for i, j in sig:
        new = []
        for s in sets:
            if i in s and j in s:
                new += [s - {i}, s - {j}]
            else:
                new.append(s)
        new = [s for s in new if s and not any(s < t for t in new)]
        sets = [s for k, s in enumerate(new) if s not in new[:k]]
    sets.sort(key=lambda s: sorted(s))
    out = {c: "".join(string.ascii_lowercase[k] for k, s in enumerate(sets) if i in s) for c, i in idx.items()}
    return out if len(set(out.values())) > 1 else None


def relabel_in_plot_order(letters, order):
    """Rename the letters so that, sweeping the cells in plot order (left to right), they appear
    as a, b, c ...: the leftmost cell carries "a" whenever it has a single letter. Only the names
    change; which cells share a letter is untouched."""
    m = {}
    for c in order:
        for ch in sorted(letters.get(c, "")):
            m.setdefault(ch, string.ascii_lowercase[len(m)])
    return {c: "".join(sorted(m.get(ch, ch) for ch in l)) for c, l in letters.items()}


def letters_from_cld_csv(path, by=("Group", "Treatment"), order=None):
    """Read a pipeline cld/<key>.csv (columns Treatment, Group[, interval], cld).
    Returns {(Group, Treatment[, interval]): letters}, or None if the file is missing or all
    cells share the same letters. order: cell keys in plot order (left to right), used to
    rename the letters so the leftmost cell gets "a" (see relabel_in_plot_order)."""
    if not path or not os.path.exists(path):
        return None
    d = pd.read_csv(path, encoding="utf-8-sig")
    if "cld" not in d.columns or not set(by) <= set(d.columns):
        return None
    d["Group"] = d.Group.map(canon_group)
    d["cld"] = d.cld.astype(str).str.strip()
    if d.cld.nunique() <= 1:
        return None
    out = {tuple(int(r[c]) if c == "interval" else r[c] for c in by): r["cld"] for _, r in d.iterrows()}
    return relabel_in_plot_order(out, order) if order else out


def draw_letters(ax, xs, tops, labels, frac=0.04, fontsize=11):
    """Write one letter label above each (x, top); frac = gap as a fraction of the y-axis span.
    Extends the upper y-limit if a label would otherwise be clipped."""
    lo, hi = ax.get_ylim()
    gap = frac * (hi - lo)
    need = hi
    ax.set_xlim(ax.get_xlim()); ax.set_ylim(lo, hi)      # freeze the limits used to test overlaps
    renderer = ax.figure.canvas.get_renderer()
    placed = []
    for x, t, l in sorted(zip(xs, tops, labels), key=lambda r: r[1]):   # lowest label first
        if not l:
            continue
        # translucent white box: a line of another series passing behind stays visible but
        # does not make the letters unreadable
        txt = ax.text(x, t + gap, l, ha="center", va="bottom", fontweight="bold", fontsize=fontsize, zorder=6,
                      bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=0.6))
        # two labels of neighbouring cells can land on the same spot: raise this one until it is
        # clear of every label already placed
        for _ in range(20):
            bb = txt.get_window_extent(renderer).expanded(1.08, 1.08)
            hit = [p for p in placed if bb.overlaps(p)]
            if not hit:
                break
            y_px = max(p.y1 for p in hit) + 2
            txt.set_y(ax.transData.inverted().transform((0, y_px))[1])
        placed.append(txt.get_window_extent(renderer))
        need = max(need, txt.get_position()[1] + 0.07 * (hi - lo))
    if need > hi:
        ax.set_ylim(lo, need)
