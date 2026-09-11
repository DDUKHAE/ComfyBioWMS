import importlib.util
from pathlib import Path
import pandas as pd
import pytest

from nodes.class_2.common import discover_samples
from nodes.class_1.tximport import Tximport
from nodes.class_2.salmon import SalmonQuantReads
from nodes.class_2.trim_galore import TrimGalore
from nodes.class_2.fastqc import FastQC
from nodes.class_2.fastp import Fastp


def test_discover_samples_single_file(tmp_path):
    fwd = tmp_path / "test_R1.fastq.gz"
    rev = tmp_path / "test_R2.fastq.gz"
    fwd.touch()
    rev.touch()

    res = discover_samples(str(fwd), str(rev))
    assert len(res) == 1
    assert res[0][0] == "test"
    assert res[0][1] == fwd.resolve()
    assert res[0][2] == rev.resolve()


def test_discover_samples_directory_paired_and_single(tmp_path):
    (tmp_path / "sampleA_1.fq.gz").touch()
    (tmp_path / "sampleA_2.fq.gz").touch()
    (tmp_path / "sampleB_R1.fastq.gz").touch()
    (tmp_path / "sampleB_R2.fastq.gz").touch()
    (tmp_path / "singleC.fastq.gz").touch()

    samples = discover_samples(str(tmp_path))
    assert len(samples) == 3
    sample_dict = {s[0]: (s[1], s[2]) for s in samples}

    assert "sampleA" in sample_dict
    assert sample_dict["sampleA"][1] is not None
    assert "sampleB" in sample_dict
    assert sample_dict["sampleB"][1] is not None
    assert "singleC" in sample_dict
    assert sample_dict["singleC"][1] is None


def test_discover_samples_comma_separated(tmp_path):
    f1 = tmp_path / "s1_1.fq.gz"
    f2 = tmp_path / "s1_2.fq.gz"
    f3 = tmp_path / "s2_1.fq.gz"
    f4 = tmp_path / "s2_2.fq.gz"
    for f in (f1, f2, f3, f4):
        f.touch()

    samples = discover_samples(f"{f1},{f3}", f"{f2},{f4}")
    assert len(samples) == 2
    assert samples[0][1] == f1.resolve()
    assert samples[0][2] == f2.resolve()
    assert samples[1][1] == f3.resolve()
    assert samples[1][2] == f4.resolve()


def test_tximport_directory_discovery(tmp_path):
    # Simulate Salmon output directory containing multiple sample subdirectories
    salmon_base = tmp_path / "SalmonQuantReads"
    s1_dir = salmon_base / "sample_WT1"
    s2_dir = salmon_base / "sample_KO1"
    s1_dir.mkdir(parents=True)
    s2_dir.mkdir(parents=True)

    header = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
    (s1_dir / "quant.sf").write_text(header + "tx1\t1000\t800\t50.0\t100\n" + "tx2\t2000\t1800\t20.0\t50\n")
    (s2_dir / "quant.sf").write_text(header + "tx1\t1000\t800\t80.0\t160\n" + "tx2\t2000\t1800\t10.0\t25\n")

    # Pass the directory path directly to Tximport
    tximport_node = Tximport()
    counts_tsv, tpm_tsv, summary_json = tximport_node.run(
        quant_files=str(salmon_base),
        output_dir=str(tmp_path / "tximport_out"),
    )

    assert Path(counts_tsv).is_file()
    assert Path(tpm_tsv).is_file()

    counts_df = pd.read_csv(counts_tsv, sep="\t")
    # Verify columns were inferred as sample folder names!
    assert "sample_WT1" in counts_df.columns
    assert "sample_KO1" in counts_df.columns
    assert len(counts_df) == 2


def test_salmon_quant_reads_batch_directory(tmp_path, monkeypatch):
    # Create input reads directory with 2 samples
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "sample1_1.fq.gz").touch()
    (reads_dir / "sample1_2.fq.gz").touch()
    (reads_dir / "sample2_1.fq.gz").touch()
    (reads_dir / "sample2_2.fq.gz").touch()

    fake_idx = tmp_path / "index"
    fake_idx.mkdir()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/salmon" if cmd == "salmon" else None)

    import nodes.class_2.salmon as salmon_mod
    def mock_run(argv, cwd):
        out_p = Path(cwd)
        (out_p / "quant.sf").write_text("dummy quant.sf")
    monkeypatch.setattr(salmon_mod, "_run", mock_run)

    node = SalmonQuantReads()
    quant_sfs, base_out = node.run(
        salmon_index_dir=str(fake_idx),
        reads_fwd=str(reads_dir),
        output_dir=str(tmp_path / "salmon_results"),
    )

    sfs = quant_sfs.split(",")
    assert len(sfs) == 2
    assert all(Path(s).is_file() for s in sfs)
    assert "sample1" in sfs[0]
    assert "sample2" in sfs[1]


def test_trim_galore_batch_directory(tmp_path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "sample1_1.fq.gz").touch()
    (reads_dir / "sample1_2.fq.gz").touch()
    (reads_dir / "sample2_1.fq.gz").touch()
    (reads_dir / "sample2_2.fq.gz").touch()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/trim_galore" if cmd == "trim_galore" else None)

    import nodes.class_2.trim_galore as tg_mod
    def mock_run(argv, cwd):
        out_p = Path(cwd)
        (out_p / "sample_val_1.fq.gz").write_text("trimmed R1")
        (out_p / "sample_val_2.fq.gz").write_text("trimmed R2")
        (out_p / "sample_trimming_report.txt").write_text("report")
    monkeypatch.setattr(tg_mod, "_run", mock_run)

    node = TrimGalore()
    fwd_out, rev_out, rep_out = node.run(
        reads_fwd=str(reads_dir),
        output_dir=str(tmp_path / "tg_results"),
    )

    fwds = fwd_out.split(",")
    revs = rev_out.split(",")
    assert len(fwds) == 2
    assert len(revs) == 2
    assert all(Path(f).is_file() for f in fwds)
    assert all(Path(r).is_file() for r in revs)


def test_fastqc_batch_directory(tmp_path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "s1.fastq.gz").touch()
    (reads_dir / "s2.fastq.gz").touch()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/fastqc" if cmd == "fastqc" else None)

    import nodes.class_2.fastqc as fqc_mod
    def mock_run(argv, cwd):
        out_p = Path(cwd)
        for target in argv:
            if "s1" in target:
                (out_p / "s1_fastqc.html").write_text("<html>s1</html>")
                (out_p / "s1_fastqc.zip").write_text("dummy")
            elif "s2" in target:
                (out_p / "s2_fastqc.html").write_text("<html>s2</html>")
                (out_p / "s2_fastqc.zip").write_text("dummy")
    monkeypatch.setattr(fqc_mod, "_run", mock_run)

    node = FastQC()
    htmls, zips = node.run(
        input_file=str(reads_dir),
        output_dir=str(tmp_path / "fastqc_results"),
    )

    html_list = htmls.split(",")
    assert len(html_list) == 2
    assert all(Path(h).is_file() for h in html_list)


def test_fastp_batch_directory(tmp_path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "s1_1.fq.gz").touch()
    (reads_dir / "s1_2.fq.gz").touch()
    (reads_dir / "s2_1.fq.gz").touch()
    (reads_dir / "s2_2.fq.gz").touch()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/fastp" if cmd == "fastp" else None)

    import nodes.class_2.fastp as fastp_mod
    def mock_run(argv, cwd):
        out_p = Path(cwd)
        for i, a in enumerate(argv):
            if a == "-o":
                Path(argv[i+1]).write_text("r1")
            elif a == "-O":
                Path(argv[i+1]).write_text("r2")
            elif a == "-j":
                Path(argv[i+1]).write_text('{"summary": {}}')
            elif a == "-h":
                Path(argv[i+1]).write_text("<html></html>")
    monkeypatch.setattr(fastp_mod, "_run", mock_run)

    node = fastp_mod.Fastp()
    r1, r2, js, ht = node.run(read1=str(reads_dir), output_dir=str(tmp_path / "fastp_out"))

    assert len(r1.split(",")) == 2
    assert len(r2.split(",")) == 2


def test_kallisto_batch_directory(tmp_path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "s1_1.fq.gz").touch()
    (reads_dir / "s1_2.fq.gz").touch()
    (reads_dir / "s2_1.fq.gz").touch()
    (reads_dir / "s2_2.fq.gz").touch()

    fake_idx = tmp_path / "index.idx"
    fake_idx.touch()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/kallisto" if cmd == "kallisto" else None)

    import subprocess
    def mock_subp_run(argv, capture_output=True, text=True):
        out_dir = Path(argv[argv.index("-o") + 1])
        (out_dir / "abundance.tsv").write_text("target_id\tlength\teff_length\test_counts\ttpm\n")
        (out_dir / "abundance.h5").write_text("h5")
        class Result:
            returncode = 0
            stderr = ""
        return Result()
    monkeypatch.setattr(subprocess, "run", mock_subp_run)

    import nodes.class_2.kallisto as kal_mod
    node = kal_mod.KallistoQuant()
    tsv_out, h5_out = node.run(
        kallisto_index=str(fake_idx),
        reads_fwd=str(reads_dir),
        output_dir=str(tmp_path / "kallisto_out"),
    )

    tsvs = tsv_out.split(",")
    assert len(tsvs) == 2
    assert all(Path(t).is_file() for t in tsvs)


def test_star_batch_directory(tmp_path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "s1_1.fq.gz").touch()
    (reads_dir / "s1_2.fq.gz").touch()
    (reads_dir / "s2_1.fq.gz").touch()
    (reads_dir / "s2_2.fq.gz").touch()

    fake_idx = tmp_path / "star_index"
    fake_idx.mkdir()

    import shutil
    monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/STAR" if cmd == "STAR" else None)

    import nodes.class_2.star as star_mod
    def mock_run(argv, cwd):
        out_p = Path(cwd)
        (out_p / "Aligned.sortedByCoord.out.bam").write_text("BAM")
        (out_p / "Aligned.toTranscriptome.out.bam").write_text("TxBAM")
        (out_p / "ReadsPerGene.out.tab").write_text("GENE")
        (out_p / "SJ.out.tab").write_text("SJ")
        (out_p / "Log.final.out").write_text("LOG")
    monkeypatch.setattr(star_mod, "_run", mock_run)

    node = star_mod.STARAlignReads()
    bams, tx_bams, counts, sjs, logs, unmapped = node.run(
        star_index_dir=str(fake_idx),
        reads_fwd=str(reads_dir),
        output_dir=str(tmp_path / "star_out"),
    )

    bam_list = bams.split(",")
    assert len(bam_list) == 2
    assert all(Path(b).is_file() for b in bam_list)
