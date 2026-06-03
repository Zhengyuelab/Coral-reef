#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 01. MAG-level KO presence/absence analysis.

Purpose
-------
This script assigns MAGs to ecological/health states using a genus-to-state
mapping table, converts MAG KO count profiles into a KO presence/absence matrix,
tests KO prevalence differences among states, and extracts state-associated
characteristic KOs and the MAGs carrying them.

Required input files
--------------------
Default input files are placed under ./data/:

1. data/genus_state_mapping.tsv
   Required columns:
   - Genus
   - State

2. data/MAG_taxonomy.tsv
   Required columns:
   - MAG_ID
   - genus or classification
   Optional columns:
   - MAG_set, domain, phylum, class, order, family, species

3. data/MAG_KO_count_matrix.tsv
   Required columns:
   - MAG_ID
   - one or more KO columns, e.g., K00001, K00002, ...

Outputs
-------
Results are written to:
MAG_KO_presence_absence_results/
"""

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

warnings.filterwarnings("ignore")

# =========================================================
# Input files
# =========================================================

DATA_DIR = Path("data")

GENUS_STATE_FILE = DATA_DIR / "genus_state_mapping.tsv"
MAG_TAXONOMY_FILE = DATA_DIR / "MAG_taxonomy.tsv"
MAG_KO_COUNT_FILE = DATA_DIR / "MAG_KO_count_matrix.tsv"

OUTDIR = Path("MAG_KO_presence_absence_results")
OUTDIR.mkdir(exist_ok=True)

STATE_ORDER = ["Good", "Moderate", "Poor"]

# =========================================================
# Main thresholds
# =========================================================

MIN_MAGS_PER_STATE = 3

# A KO must occur in at least this many MAGs overall before statistical testing.
MIN_TOTAL_PRESENT_MAGS = 2

# A state-associated KO must occur in at least this many MAGs in the target state.
MIN_TARGET_PRESENT_MAGS = 3

# Minimum target-state prevalence.
MIN_TARGET_PREVALENCE = 0.10

# Minimum prevalence difference between the target state and the mean of other states.
MIN_PREVALENCE_DIFF = 0.20

# Minimum prevalence ratio between the target state and the mean of other states.
MIN_PREVALENCE_RATIO = 2.0

# FDR cutoff for the overall chi-square test.
FDR_CUTOFF = 0.05

# Maximum number of characteristic KOs retained per state.
# Use None to keep all qualified KOs.
TOP_N_PER_STATE = 30

PSEUDO = 1e-6

# =========================================================
# Helper functions
# =========================================================

def clean_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def normalize_genus(x):
    x = clean_text(x)
    if x in ["", "-", "NA", "NaN", "nan", "None"]:
        return ""
    x = x.replace(";", "")
    if not x.startswith("g__"):
        x = "g__" + x
    if x == "g__":
        return ""
    return x

def normalize_state(x):
    x = clean_text(x)
    lx = x.lower()
    if lx == "good":
        return "Good"
    if lx == "moderate":
        return "Moderate"
    if lx == "poor":
        return "Poor"
    return x

def extract_genus_from_classification(classification):
    classification = clean_text(classification)
    for p in classification.split(";"):
        p = p.strip()
        if p.startswith("g__"):
            return normalize_genus(p)
    return ""

def bh_fdr(pvalues):
    p = np.array(pvalues, dtype=float)
    q = np.full(len(p), np.nan)

    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return q

    order = np.argsort(pv)
    ranked = pv[order]
    q_ranked = ranked * n / (np.arange(n) + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.minimum(q_ranked, 1.0)

    q_valid = np.empty(n)
    q_valid[order] = q_ranked
    q[valid] = q_valid
    return q

def read_count_matrix(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]

    if "MAG_ID" not in df.columns:
        raise SystemExit(f"ERROR: MAG_ID column not found in {path}")

    ko_cols = [c for c in df.columns if c != "MAG_ID"]

    if len(ko_cols) == 0:
        raise SystemExit(f"ERROR: no KO columns found in {path}")

    df[ko_cols] = df[ko_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    return df

def load_taxonomy(path, genus_to_state):
    tax = pd.read_csv(path, sep="\t", dtype=str)
    tax.columns = [str(c).strip() for c in tax.columns]

    if "MAG_ID" not in tax.columns:
        raise SystemExit(f"ERROR: MAG_ID column not found in {path}")

    if "genus" in tax.columns:
        tax["Genus_norm"] = tax["genus"].map(normalize_genus)
    elif "classification" in tax.columns:
        tax["Genus_norm"] = tax["classification"].map(extract_genus_from_classification)
    else:
        raise SystemExit(f"ERROR: no genus/classification column found in {path}")

    if "MAG_set" not in tax.columns:
        tax["MAG_set"] = "MAG_set_1"

    tax["Matched_state"] = tax["Genus_norm"].map(genus_to_state).fillna("Unmatched")

    keep = ["MAG_set", "MAG_ID", "Genus_norm", "Matched_state"]
    for c in ["classification", "domain", "phylum", "class", "order", "family", "genus", "species"]:
        if c in tax.columns and c not in keep:
            keep.append(c)

    return tax[keep].copy()

# =========================================================
# 1. Genus -> state mapping
# =========================================================

print("===== Step 1: Read genus-state mapping =====")

mapping = pd.read_csv(GENUS_STATE_FILE, sep="\t", dtype=str)
mapping.columns = [str(c).strip() for c in mapping.columns]

if "Genus" not in mapping.columns:
    raise SystemExit("ERROR: 'Genus' column not found in genus-state mapping table.")

state_col = None
for candidate in ["State", "Health_status", "health status", "Group"]:
    if candidate in mapping.columns:
        state_col = candidate
        break

if state_col is None:
    raise SystemExit("ERROR: no state column found. Use one of: State, Health_status, health status, Group.")

mapping["Genus_norm"] = mapping["Genus"].map(normalize_genus)
mapping["State"] = mapping[state_col].map(normalize_state)

mapping = mapping[
    (mapping["Genus_norm"] != "") &
    (mapping["State"].isin(STATE_ORDER))
].copy()

# If a genus maps to multiple states, it is treated as ambiguous and excluded.
genus_state = (
    mapping.groupby("Genus_norm")["State"]
    .apply(lambda x: sorted(set(x)))
    .reset_index()
)
genus_state["n_states"] = genus_state["State"].apply(len)

conflict = genus_state[genus_state["n_states"] > 1].copy()
conflict.to_csv(OUTDIR / "00_conflicting_genus_state_mapping.tsv", sep="\t", index=False)

unique_genus = genus_state[genus_state["n_states"] == 1].copy()
unique_genus["State"] = unique_genus["State"].apply(lambda x: x[0])

genus_to_state = dict(zip(unique_genus["Genus_norm"], unique_genus["State"]))

unique_genus[["Genus_norm", "State"]].to_csv(
    OUTDIR / "01_Genus_to_state_mapping.tsv",
    sep="\t",
    index=False
)

print("Raw genus-state rows:", len(mapping))
print("Unambiguous genus:", len(genus_to_state))
print("Conflicting genus:", len(conflict))

# =========================================================
# 2. Assign state to MAGs
# =========================================================

print("===== Step 2: Assign MAG state =====")

mag_state = load_taxonomy(MAG_TAXONOMY_FILE, genus_to_state)

mag_state.to_csv(OUTDIR / "02_MAG_state_assignment.tsv", sep="\t", index=False)

matched_mag_state = mag_state[mag_state["Matched_state"].isin(STATE_ORDER)].copy()
matched_mag_state.to_csv(OUTDIR / "02_MAG_state_assignment_matched_only.tsv", sep="\t", index=False)

print("Total MAG taxonomy rows:", len(mag_state))
print("Matched MAGs:", len(matched_mag_state))
print(matched_mag_state["Matched_state"].value_counts().reindex(STATE_ORDER).fillna(0).astype(int))

# =========================================================
# 3. Read KO count matrix and convert to presence/absence
# =========================================================

print("===== Step 3: Read KO count matrix and convert to presence/absence =====")

combined_count = read_count_matrix(MAG_KO_COUNT_FILE)

if "MAG_set" not in combined_count.columns:
    combined_count.insert(0, "MAG_set", "MAG_set_1")

ko_cols = [c for c in combined_count.columns if c not in ["MAG_set", "MAG_ID"]]
ko_cols = sorted(ko_cols)

combined_count = combined_count[["MAG_set", "MAG_ID"] + ko_cols]
combined_count.to_csv(OUTDIR / "03_MAG_KO_count_matrix.tsv", sep="\t", index=False)

combined_pa = combined_count.copy()
combined_pa[ko_cols] = (combined_pa[ko_cols] > 0).astype(int)
combined_pa.to_csv(OUTDIR / "04_MAG_KO_presence_absence_matrix.tsv", sep="\t", index=False)

matched_count = combined_count.merge(
    matched_mag_state[["MAG_ID", "Genus_norm", "Matched_state"]],
    on="MAG_ID",
    how="inner"
)

matched_pa = combined_pa.merge(
    matched_mag_state[["MAG_ID", "Genus_norm", "Matched_state"]],
    on="MAG_ID",
    how="inner"
)

meta_cols = ["MAG_set", "MAG_ID", "Genus_norm", "Matched_state"]
matched_count = matched_count[meta_cols + ko_cols]
matched_pa = matched_pa[meta_cols + ko_cols]

matched_count.to_csv(OUTDIR / "05_MAG_KO_count_matrix_matched.tsv", sep="\t", index=False)
matched_pa.to_csv(OUTDIR / "06_MAG_KO_presence_absence_matrix_matched.tsv", sep="\t", index=False)

print("Combined MAGs:", len(combined_count))
print("Matched MAGs:", len(matched_pa))
print("KO columns:", len(ko_cols))
print("Matched MAGs per state:")
print(matched_pa["Matched_state"].value_counts().reindex(STATE_ORDER).fillna(0).astype(int))

# =========================================================
# 4. KO prevalence summary
# =========================================================

print("===== Step 4: KO prevalence summary =====")

summary_rows = []

for ko in ko_cols:
    row = {"KEGG_KO": ko}

    total_present = int(matched_pa[ko].sum())
    row["Total_present_MAGs"] = total_present
    row["Total_prevalence"] = float(total_present / len(matched_pa)) if len(matched_pa) else np.nan
    row["Total_copy_number"] = int(matched_count[ko].sum())

    for state in STATE_ORDER:
        idx = matched_pa["Matched_state"] == state
        n = int(idx.sum())
        present = int(matched_pa.loc[idx, ko].sum())
        copy_sum = int(matched_count.loc[idx, ko].sum())

        row[f"{state}_n_MAGs"] = n
        row[f"{state}_present_MAGs"] = present
        row[f"{state}_prevalence"] = float(present / n) if n > 0 else np.nan
        row[f"{state}_copy_sum"] = copy_sum
        row[f"{state}_mean_copy_number"] = float(matched_count.loc[idx, ko].mean()) if n > 0 else np.nan

    prevs = {state: row[f"{state}_prevalence"] for state in STATE_ORDER}
    max_state = max(prevs, key=lambda s: prevs[s] if np.isfinite(prevs[s]) else -1)
    other_states = [s for s in STATE_ORDER if s != max_state]
    other_prev_mean = np.mean([prevs[s] for s in other_states])

    row["Max_prevalence_state"] = max_state
    row["Max_prevalence"] = prevs[max_state]
    row["Other_states_prevalence_mean"] = other_prev_mean
    row["Prevalence_diff_max_vs_others"] = row["Max_prevalence"] - other_prev_mean
    row["Prevalence_ratio_max_vs_others"] = (row["Max_prevalence"] + PSEUDO) / (other_prev_mean + PSEUDO)

    summary_rows.append(row)

prevalence_df = pd.DataFrame(summary_rows)
prevalence_df.to_csv(OUTDIR / "07_KO_prevalence_by_state.tsv", sep="\t", index=False)

# =========================================================
# 5. Overall chi-square test
# =========================================================

print("===== Step 5: Overall chi-square test =====")

test_rows = []

for ko in ko_cols:
    total_present = int(matched_pa[ko].sum())
    if total_present < MIN_TOTAL_PRESENT_MAGS:
        continue

    table = []
    sufficient = True

    for state in STATE_ORDER:
        idx = matched_pa["Matched_state"] == state
        n = int(idx.sum())
        present = int(matched_pa.loc[idx, ko].sum())
        absent = n - present

        if n < MIN_MAGS_PER_STATE:
            sufficient = False

        table.append([present, absent])

    row = {"KEGG_KO": ko}

    if not sufficient:
        row["Chi2_stat"] = np.nan
        row["Chi2_p"] = np.nan
        row["Test_status"] = "Skipped_insufficient_MAGs"
    else:
        try:
            if np.sum(table) == 0:
                row["Chi2_stat"] = np.nan
                row["Chi2_p"] = np.nan
                row["Test_status"] = "No_data"
            else:
                stat, p, dof, expected = chi2_contingency(np.array(table), correction=False)
                row["Chi2_stat"] = stat
                row["Chi2_p"] = p
                row["Chi2_dof"] = dof
                row["Test_status"] = "OK"
        except Exception as e:
            row["Chi2_stat"] = np.nan
            row["Chi2_p"] = np.nan
            row["Test_status"] = f"Error:{e}"

    test_rows.append(row)

chi_df = pd.DataFrame(test_rows)
chi_df["Chi2_FDR"] = bh_fdr(chi_df["Chi2_p"].values)

chi_df = chi_df.merge(prevalence_df, on="KEGG_KO", how="left")
chi_df.to_csv(OUTDIR / "08_KO_presence_absence_overall_chisq.tsv", sep="\t", index=False)

# =========================================================
# 6. Pairwise Fisher exact tests
# =========================================================

print("===== Step 6: Pairwise Fisher exact tests =====")

pairs = [("Good", "Moderate"), ("Good", "Poor"), ("Moderate", "Poor")]
pair_rows = []

for ko in ko_cols:
    total_present = int(matched_pa[ko].sum())
    if total_present < MIN_TOTAL_PRESENT_MAGS:
        continue

    for a, b in pairs:
        ia = matched_pa["Matched_state"] == a
        ib = matched_pa["Matched_state"] == b

        na = int(ia.sum())
        nb = int(ib.sum())

        pa = int(matched_pa.loc[ia, ko].sum())
        pb = int(matched_pa.loc[ib, ko].sum())

        aa = na - pa
        ab = nb - pb

        row = {
            "KEGG_KO": ko,
            "Group1": a,
            "Group2": b,
            "Group1_n_MAGs": na,
            "Group2_n_MAGs": nb,
            "Group1_present_MAGs": pa,
            "Group2_present_MAGs": pb,
            "Group1_prevalence": pa / na if na > 0 else np.nan,
            "Group2_prevalence": pb / nb if nb > 0 else np.nan,
            "Prevalence_diff_Group1_minus_Group2": (pa / na if na > 0 else np.nan) - (pb / nb if nb > 0 else np.nan),
            "Prevalence_ratio_Group1_vs_Group2": ((pa / na if na > 0 else 0) + PSEUDO) / ((pb / nb if nb > 0 else 0) + PSEUDO),
        }

        if na < MIN_MAGS_PER_STATE or nb < MIN_MAGS_PER_STATE:
            row["Odds_ratio"] = np.nan
            row["Fisher_p"] = np.nan
            row["Test_status"] = "Skipped_insufficient_MAGs"
        else:
            try:
                odds, p = fisher_exact([[pa, aa], [pb, ab]], alternative="two-sided")
                row["Odds_ratio"] = odds
                row["Fisher_p"] = p
                row["Test_status"] = "OK"
            except Exception as e:
                row["Odds_ratio"] = np.nan
                row["Fisher_p"] = np.nan
                row["Test_status"] = f"Error:{e}"

        pair_rows.append(row)

pair_df = pd.DataFrame(pair_rows)
pair_df["Fisher_FDR"] = bh_fdr(pair_df["Fisher_p"].values)
pair_df.to_csv(OUTDIR / "09_KO_presence_absence_pairwise_fisher.tsv", sep="\t", index=False)

# =========================================================
# 7. Main characteristic KOs by presence/absence
# =========================================================

print("===== Step 7: Select main characteristic KOs =====")

main_rows = []

for _, r in chi_df.iterrows():
    ko = r["KEGG_KO"]

    if not np.isfinite(r["Chi2_FDR"]):
        continue

    if r["Chi2_FDR"] > FDR_CUTOFF:
        continue

    for state in STATE_ORDER:
        other_states = [s for s in STATE_ORDER if s != state]

        target_prev = float(r[f"{state}_prevalence"])
        target_present = int(r[f"{state}_present_MAGs"])
        target_n = int(r[f"{state}_n_MAGs"])
        target_copy_sum = int(r[f"{state}_copy_sum"])
        target_mean_copy = float(r[f"{state}_mean_copy_number"])

        other_prev_mean = np.mean([float(r[f"{s}_prevalence"]) for s in other_states])
        other_present_sum = int(sum([int(r[f"{s}_present_MAGs"]) for s in other_states]))
        other_copy_sum = int(sum([int(r[f"{s}_copy_sum"]) for s in other_states]))

        prev_diff = target_prev - other_prev_mean
        prev_ratio = (target_prev + PSEUDO) / (other_prev_mean + PSEUDO)

        if r["Max_prevalence_state"] != state:
            continue
        if target_present < MIN_TARGET_PRESENT_MAGS:
            continue
        if target_prev < MIN_TARGET_PREVALENCE:
            continue
        if prev_diff < MIN_PREVALENCE_DIFF:
            continue
        if prev_ratio < MIN_PREVALENCE_RATIO:
            continue

        score = (
            -math.log10(float(r["Chi2_FDR"]) + 1e-300)
            * prev_diff
            * math.log1p(target_present)
            * prev_ratio
        )

        main_rows.append({
            "State": state,
            "KEGG_KO": ko,
            "Chi2_p": r["Chi2_p"],
            "Chi2_FDR": r["Chi2_FDR"],
            "Target_n_MAGs": target_n,
            "Target_present_MAGs": target_present,
            "Target_prevalence": target_prev,
            "Other_states_present_MAGs": other_present_sum,
            "Other_states_prevalence_mean": other_prev_mean,
            "Prevalence_diff_state_vs_others": prev_diff,
            "Prevalence_ratio_state_vs_others": prev_ratio,
            "Target_copy_sum": target_copy_sum,
            "Other_states_copy_sum": other_copy_sum,
            "Target_mean_copy_number": target_mean_copy,
            "Characteristic_score": score,
            "Good_present_MAGs": r["Good_present_MAGs"],
            "Moderate_present_MAGs": r["Moderate_present_MAGs"],
            "Poor_present_MAGs": r["Poor_present_MAGs"],
            "Good_prevalence": r["Good_prevalence"],
            "Moderate_prevalence": r["Moderate_prevalence"],
            "Poor_prevalence": r["Poor_prevalence"],
            "Good_copy_sum": r["Good_copy_sum"],
            "Moderate_copy_sum": r["Moderate_copy_sum"],
            "Poor_copy_sum": r["Poor_copy_sum"],
            "Good_mean_copy_number": r["Good_mean_copy_number"],
            "Moderate_mean_copy_number": r["Moderate_mean_copy_number"],
            "Poor_mean_copy_number": r["Poor_mean_copy_number"],
        })

main_df = pd.DataFrame(main_rows)

if len(main_df) > 0:
    main_df = main_df.sort_values(
        ["State", "Characteristic_score", "Prevalence_diff_state_vs_others", "Target_prevalence"],
        ascending=[True, False, False, False]
    )

    if TOP_N_PER_STATE is not None:
        main_df = (
            main_df.groupby("State", group_keys=False)
            .head(TOP_N_PER_STATE)
            .copy()
        )

main_df.to_csv(OUTDIR / "10_characteristic_KOs_by_presence_absence.tsv", sep="\t", index=False)

print("Main characteristic KOs:")
if len(main_df) > 0:
    print(main_df["State"].value_counts().reindex(STATE_ORDER).fillna(0).astype(int))
else:
    print("None")

# =========================================================
# 8. Extract feature MAGs and KO copy number
# =========================================================

print("===== Step 8: Extract MAGs carrying characteristic KOs =====")

detail_rows = []

if len(main_df) > 0:
    for _, fr in main_df.iterrows():
        state = fr["State"]
        ko = fr["KEGG_KO"]

        sub = matched_pa[
            (matched_pa["Matched_state"] == state) &
            (matched_pa[ko] > 0)
        ].copy()

        for _, mr in sub.iterrows():
            count_row = matched_count.loc[matched_count["MAG_ID"] == mr["MAG_ID"], ko]
            copy_number = int(count_row.iloc[0]) if len(count_row) > 0 else 0

            detail_rows.append({
                "State": state,
                "MAG_set": mr["MAG_set"],
                "MAG_ID": mr["MAG_ID"],
                "Genus": mr["Genus_norm"],
                "KEGG_KO": ko,
                "KO_presence": 1,
                "KO_copy_number_in_MAG": copy_number,
                "Chi2_FDR": fr["Chi2_FDR"],
                "Target_prevalence": fr["Target_prevalence"],
                "Other_states_prevalence_mean": fr["Other_states_prevalence_mean"],
                "Prevalence_diff_state_vs_others": fr["Prevalence_diff_state_vs_others"],
                "Prevalence_ratio_state_vs_others": fr["Prevalence_ratio_state_vs_others"],
                "Characteristic_score": fr["Characteristic_score"],
            })

detail_df = pd.DataFrame(detail_rows)
detail_df.to_csv(OUTDIR / "11_MAG_characteristic_KO_detail.tsv", sep="\t", index=False)

if len(detail_df) > 0:
    def join_unique(x):
        return ";".join(sorted(set(map(str, x))))

    def join_ko_copy(df):
        vals = []
        for _, rr in df.iterrows():
            vals.append(f"{rr['KEGG_KO']}:{rr['KO_copy_number_in_MAG']}")
        return ";".join(sorted(vals))

    mag_summary = (
        detail_df.groupby(["State", "MAG_set", "MAG_ID", "Genus"], as_index=False)
        .agg(
            Characteristic_KO_number=("KEGG_KO", lambda x: len(set(x))),
            Total_characteristic_KO_copy_number=("KO_copy_number_in_MAG", "sum"),
            Characteristic_KO_list=("KEGG_KO", join_unique),
        )
    )

    copy_detail = (
        detail_df.groupby(["State", "MAG_set", "MAG_ID", "Genus"])
        .apply(join_ko_copy)
        .reset_index()
        .rename(columns={0: "Characteristic_KO_copy_detail"})
    )

    mag_summary = mag_summary.merge(
        copy_detail,
        on=["State", "MAG_set", "MAG_ID", "Genus"],
        how="left"
    )

    mag_summary = mag_summary.sort_values(
        ["State", "Characteristic_KO_number", "Total_characteristic_KO_copy_number"],
        ascending=[True, False, False]
    )
else:
    mag_summary = pd.DataFrame()

mag_summary.to_csv(OUTDIR / "12_MAG_characteristic_KO_summary.tsv", sep="\t", index=False)

# =========================================================
# 9. Heatmap PDF for characteristic KO prevalence
# =========================================================

print("===== Step 9: Plot KO prevalence heatmap =====")

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    heatmap_pdf = OUTDIR / "13_characteristic_KO_prevalence_heatmap.pdf"
    heatmap_matrix = OUTDIR / "13_characteristic_KO_prevalence_heatmap_matrix.tsv"

    if len(main_df) > 0:
        hm = main_df.copy()
        hm["State"] = pd.Categorical(hm["State"], categories=STATE_ORDER, ordered=True)
        hm = hm.sort_values(["State", "Characteristic_score"], ascending=[True, False])

        mat = hm[["Good_prevalence", "Moderate_prevalence", "Poor_prevalence"]].astype(float)
        mat.columns = ["Good", "Moderate", "Poor"]

        # Row Z-score highlights state specificity.
        mat_z = mat.copy()
        row_mean = mat_z.mean(axis=1)
        row_std = mat_z.std(axis=1, ddof=0).replace(0, np.nan)
        mat_z = mat_z.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0)

        heatmap_out = pd.concat(
            [hm[["State", "KEGG_KO"]].reset_index(drop=True), mat_z.reset_index(drop=True)],
            axis=1
        )
        heatmap_out.to_csv(heatmap_matrix, sep="\t", index=False)

        n_rows = mat_z.shape[0]
        fig_height = max(6, 0.22 * n_rows + 2.5)
        fig_width = 7

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        vmax = max(2, np.nanmax(np.abs(mat_z.values)))
        im = ax.imshow(
            mat_z.values,
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax
        )

        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(["Good", "Moderate", "Poor"], fontsize=11)

        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(hm["KEGG_KO"].astype(str).tolist(), fontsize=6 if n_rows > 60 else 8)

        state_list = hm["State"].astype(str).tolist()
        for i in range(1, len(state_list)):
            if state_list[i] != state_list[i - 1]:
                ax.axhline(i - 0.5, color="black", linewidth=0.8)

        group_positions = {}
        for state in STATE_ORDER:
            idx = [i for i, x in enumerate(state_list) if x == state]
            if idx:
                group_positions[state] = (min(idx) + max(idx)) / 2

        for state, ypos in group_positions.items():
            ax.text(-0.9, ypos, state, va="center", ha="right", fontsize=10, fontweight="bold", rotation=90)

        ax.set_title("Characteristic KOs based on MAG-level presence/absence prevalence", fontsize=12, pad=12)
        ax.set_xlabel("State")
        ax.set_ylabel("Characteristic KOs")

        cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
        cbar.set_label("Row Z-score of KO prevalence", fontsize=9)

        plt.tight_layout()

        with PdfPages(heatmap_pdf) as pdf:
            pdf.savefig(fig, bbox_inches="tight")

        plt.close(fig)

        print("Saved heatmap:", heatmap_pdf)
        print("Saved heatmap matrix:", heatmap_matrix)

except Exception as e:
    print("WARNING: heatmap plotting failed:", e)
    print("The statistical results were still generated.")

# =========================================================
# 10. Summary
# =========================================================

summary_file = OUTDIR / "00_analysis_summary.txt"

with open(summary_file, "w", encoding="utf-8") as f:
    f.write("===== Analysis type =====\n")
    f.write("Main analysis based on KO presence/absence at the MAG level.\n")
    f.write("KO copy number is retained only for descriptive summaries of MAGs carrying characteristic KOs.\n\n")

    f.write("===== Genus-state mapping =====\n")
    f.write(f"Raw genus-state rows\t{len(mapping)}\n")
    f.write(f"Unambiguous genus\t{len(genus_to_state)}\n")
    f.write(f"Conflicting genus\t{len(conflict)}\n\n")

    f.write("===== MAG matching =====\n")
    f.write(f"Total MAG taxonomy rows\t{len(mag_state)}\n")
    f.write(f"Matched MAGs\t{len(matched_pa)}\n")
    for state in STATE_ORDER:
        f.write(f"Matched MAGs in {state}\t{int((matched_pa['Matched_state'] == state).sum())}\n")

    f.write("\n===== KO matrix =====\n")
    f.write(f"Total KO columns\t{len(ko_cols)}\n")
    f.write(f"KOs tested by chi-square\t{len(chi_df)}\n\n")

    f.write("===== Thresholds =====\n")
    f.write(f"FDR_CUTOFF\t{FDR_CUTOFF}\n")
    f.write(f"MIN_TARGET_PRESENT_MAGS\t{MIN_TARGET_PRESENT_MAGS}\n")
    f.write(f"MIN_TARGET_PREVALENCE\t{MIN_TARGET_PREVALENCE}\n")
    f.write(f"MIN_PREVALENCE_DIFF\t{MIN_PREVALENCE_DIFF}\n")
    f.write(f"MIN_PREVALENCE_RATIO\t{MIN_PREVALENCE_RATIO}\n")
    f.write(f"TOP_N_PER_STATE\t{TOP_N_PER_STATE}\n\n")

    f.write("===== Characteristic KOs =====\n")
    if len(main_df) > 0:
        for state in STATE_ORDER:
            f.write(f"{state}\t{int((main_df['State'] == state).sum())}\n")
    else:
        for state in STATE_ORDER:
            f.write(f"{state}\t0\n")

    f.write("\n===== MAGs carrying characteristic KOs =====\n")
    if len(mag_summary) > 0:
        for state in STATE_ORDER:
            f.write(f"{state}\t{mag_summary.loc[mag_summary['State'] == state, 'MAG_ID'].nunique()}\n")
    else:
        for state in STATE_ORDER:
            f.write(f"{state}\t0\n")

print("===== Presence/absence analysis finished =====")
print("Output directory:", OUTDIR)
print("Summary:", summary_file)
print("Characteristic KOs:", OUTDIR / "10_characteristic_KOs_by_presence_absence.tsv")
print("MAG-KO detail:", OUTDIR / "11_MAG_characteristic_KO_detail.tsv")
print("MAG summary:", OUTDIR / "12_MAG_characteristic_KO_summary.tsv")
