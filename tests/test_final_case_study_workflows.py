import importlib.util
from pathlib import Path


MODULE = Path("engine/scripts/generate_case_study_workflows.py").resolve()


def load_module():
    spec = importlib.util.spec_from_file_location("case_study_workflows", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_final_case_set_excludes_invalid_phix_reads_and_labels_claim_scope():
    mod = load_module()

    assert [case[0] for case in mod.FINAL_CASES] == [
        "final_01_pancreas_scanpy_e2e",
        "final_02_zymo_metaphlan_e2e",
        "final_03_gse110004_chri_smoke",
    ]

    for _, steps, title, benchmark in mod.FINAL_CASES:
        graph = mod.build(steps, title, benchmark)
        mod.validate(graph)
        assert graph["extra"]["benchmark"] == benchmark
        assert "/Users/" not in str(graph)
        assert "paper_phix174/reads_" not in str(graph)

    by_name = {case[0]: case for case in mod.FINAL_CASES}
    pancreas = by_name["final_01_pancreas_scanpy_e2e"][3]
    assert pancreas["claim_level"] == "technical_e2e"
    assert pancreas["dimensions"] == {"cells": 3696, "genes": 27998}

    zymo = by_name["final_02_zymo_metaphlan_e2e"][3]
    assert zymo["runtime_database_required"] is True
    assert zymo["read_pairs"] == 100000
    assert set(zymo["dataset_sha256"]) == {"reads_R1.fastq.gz", "reads_R2.fastq.gz"}

    gse = by_name["final_03_gse110004_chri_smoke"][3]
    assert gse["claim_level"] == "smoke"
    assert gse["read_pairs_per_run"] == 50000
    assert gse["runs"] == ["SRR6357070", "SRR6357072", "SRR6357076", "SRR6357077"]

    phix = next(item for item in mod.FINAL_STATUS if item["case_id"] == "CS-E2")
    assert phix["status"] == "blocked"
    assert "SRR5458066" in phix["reason"]
