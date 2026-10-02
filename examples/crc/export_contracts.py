#!/usr/bin/env python3
"""Export the frozen CRC panel results to general v1 workflow contracts."""
import argparse
import csv
import json
from pathlib import Path


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        yield from csv.DictReader(stream)


def write(path, fieldnames, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def source(study):
    return "geo" if study == "GSE236581" else "scctdb"


def endpoint(study):
    return "RECIST" if study == "GSE236581" else "pCR_non_pCR"


def main(root):
    panel = root / "panel"
    output = root / "processed" / "contracts_v1"
    obs = []
    for r in rows(panel / "03_tumor_treatment_trajectories" / "source" / "trajectory_plot_data.csv"):
        obs.append({"source": source(r["cohort"]), "study_id": r["cohort"],
                    "patient_id": r["patient"], "timepoint": r["stage"],
                    "compartment": "malignant_epithelial", "feature": r["gene"],
                    "value": r["expression"], "unit": "log2(CPM+1)",
                    "n_cells": r["n_malignant_cells"], "response": r["response"],
                    "response_endpoint": endpoint(r["cohort"]),
                    "annotation_source": r["annotation_source"],
                    "raw_count": r["raw_count"], "library_size": r["full_gene_library"]})
    obs_cols = ("source", "study_id", "patient_id", "timepoint", "compartment",
                "feature", "value", "unit", "n_cells", "response",
                "response_endpoint", "annotation_source", "raw_count", "library_size")
    write(output / "patient_observations.csv", obs_cols, obs)
    statistics = []
    for panel_name, analysis_id in (("01_MHC_I_correlation", "mhc_i_correlation"),
                                    ("02_T_cell_state_correlation", "t_cell_state_correlation")):
        for r in rows(panel / panel_name / "source" / "correlation_statistics.csv"):
            statistics.append({"analysis_id": analysis_id, "source": source(r["cohort"]),
                               "study_id": r["cohort"], "comparison": r["outcome"],
                               "feature": r["target"], "method": r["p_method"],
                               "n": r["n"], "estimate": r["rho"], "p_value": r["p"],
                               "adjusted_p": r["q"], "adjustment_family": r["family"],
                               "status": r["status"], "ci_low": r["ci_low"],
                               "ci_high": r["ci_high"]})
    for r in rows(panel / "03_tumor_treatment_trajectories" / "source" / "paired_t_test_statistics.csv"):
        statistics.append({"analysis_id": "treatment_change", "source": source(r["cohort"]),
                           "study_id": r["cohort"],
                           "comparison": f"{r['stage']}_vs_Pre_{r['group']}",
                           "feature": r["gene"], "method": "paired_t_two_sided",
                           "n": r["n_pairs"], "estimate": r["mean_difference"],
                           "p_value": r["p"], "adjusted_p": r["q"],
                           "adjustment_family": f"{r['cohort']}_{r['group']}",
                           "status": "tested" if r["p"] else "insufficient_pairs",
                           "ci_low": r["ci_low"], "ci_high": r["ci_high"]})
    stat_cols = ("analysis_id", "source", "study_id", "comparison", "feature",
                 "method", "n", "estimate", "p_value", "adjusted_p",
                 "adjustment_family", "status", "ci_low", "ci_high")
    write(output / "statistical_results.csv", stat_cols, statistics)
    (output / "manifest.json").write_text(
        json.dumps({"contract_version": "1.0", "expression_unit": "log2(CPM+1)",
                    "files": {"patient_observations.csv": len(obs),
                              "statistical_results.csv": len(statistics)}}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Exported {len(obs)} observations and {len(statistics)} results to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    main(args.work_root.resolve())
