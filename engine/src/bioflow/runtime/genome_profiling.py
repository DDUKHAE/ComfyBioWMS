from pathlib import Path

from bioflow.runtime.command_runner import CondaCommandRunner, conda_command

ENV_NAME = "genome_assembly"


def estimate_genome_size_mb(fastq_dir: Path, *, kmer_size: int = 21, runner=None) -> float:
    """Rough haploid genome size estimate from a k-mer coverage histogram peak.

    Runs `jellyfish count` then `jellyfish histo` against every FASTQ file in `fastq_dir`,
    finds the histogram bin with the most distinct k-mers (the coverage peak — the first
    few low-depth bins are sequencing-error noise and are skipped), and estimates size as
    (k-mers in the peak bin) / (peak depth) / 1,000,000, in megabases. This is a coarse
    scale check (microbial ~1-15 Mb vs. eukaryotic 100s-1000s of Mb), not a precise
    genome-size estimate — it does not need GenomeScope2's full model fit for that purpose.
    """
    fastq_dir = Path(fastq_dir)
    if not fastq_dir.is_dir():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")
    fastq_files = sorted(str(p) for p in fastq_dir.glob("*.fastq")) + sorted(str(p) for p in fastq_dir.glob("*.fastq.gz"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in: {fastq_dir}")

    runner = runner or CondaCommandRunner()
    jf_index = fastq_dir / "reads.jf"
    runner.run(
        conda_command(
            ENV_NAME,
            "jellyfish",
            "count",
            "-m",
            str(kmer_size),
            "-s",
            "100M",
            "-t",
            "4",
            "-o",
            str(jf_index),
            *fastq_files,
        ),
        cwd=fastq_dir,
    )
    histo_record = runner.run(
        conda_command(ENV_NAME, "jellyfish", "histo", str(jf_index)),
        cwd=fastq_dir,
    )

    peak_depth = 0
    peak_count = 0
    for line in histo_record.stdout.strip().splitlines():
        depth_str, count_str = line.split()
        depth, count = int(depth_str), int(count_str)
        if depth <= 3:
            continue  # skip the sequencing-error spike at very low depth
        if count > peak_count:
            peak_depth, peak_count = depth, count

    if peak_depth == 0:
        raise ValueError(f"Could not find a coverage peak in jellyfish histogram for: {fastq_dir}")
    return (peak_count / peak_depth) / 1_000_000
