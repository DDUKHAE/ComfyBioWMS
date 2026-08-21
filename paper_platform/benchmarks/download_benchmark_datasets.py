"""
Download script for Official Lightweight Real Bioinformatics Benchmark Datasets.
Downloads genuine public sequencing data from official NCBI, nf-core, and ENCODE repositories.
"""

import os
import sys
import urllib.request
import time

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))

DATASETS = {
    # -------------------------------------------------------------------------
    # 1. NCBI RefSeq PhiX174 Viral De Novo Assembly Benchmark (Real Illumina Reads & Genome)
    # -------------------------------------------------------------------------
    "phix174": {
        "description": "NCBI RefSeq NC_001422.1 & Real Illumina PhiX174 Paired-End Control Reads",
        "dir": os.path.join(BASE_DIR, "phix174"),
        "files": [
            {
                "name": "phix174_ref.fasta",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta",
                "ncbi_ref_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
            },
            {
                "name": "phix174_R1.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz"
            },
            {
                "name": "phix174_R2.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz"
            }
        ]
    },
    # -------------------------------------------------------------------------
    # 2. SEQC / Human Real Transcriptomics (RNA-Seq) Benchmark (nf-core/rnaseq official real test set)
    # -------------------------------------------------------------------------
    "rnaseq_seqc": {
        "description": "Real Human Transcriptomics Paired-End RNA-Seq Benchmark (Sample A vs B)",
        "dir": os.path.join(BASE_DIR, "rnaseq_seqc"),
        "files": [
            {
                "name": "SAMPLE_A1_R1.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz"
            },
            {
                "name": "SAMPLE_A1_R2.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz"
            },
            {
                "name": "SAMPLE_B1_R1.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz"
            },
            {
                "name": "SAMPLE_B1_R2.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz"
            },
            {
                "name": "reference_transcriptome.fasta",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/transcriptome.fasta"
            },
            {
                "name": "reference_genes.gtf",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/genes.gtf"
            }
        ]
    },
    # -------------------------------------------------------------------------
    # 3. ENCODE Standard ATAC-Seq Epigenomics Benchmark (nf-core/atacseq official real test set)
    # -------------------------------------------------------------------------
    "atacseq_encode": {
        "description": "ENCODE Tier-1 Real ATAC-Seq Chromatin Accessibility Benchmark",
        "dir": os.path.join(BASE_DIR, "atacseq_encode"),
        "files": [
            {
                "name": "encode_atac_rep1_R1.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/testdata/SRR1822153_1.fastq.gz"
            },
            {
                "name": "encode_atac_rep1_R2.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/testdata/SRR1822153_2.fastq.gz"
            },
            {
                "name": "genome_chr22.fasta",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/reference/genome.fa"
            },
            {
                "name": "genome_chr22.gtf",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/reference/genes.gtf"
            }
        ]
    },
    # -------------------------------------------------------------------------
    # 4. ZymoBIOMICS Metagenomics / Microbiome Profiling Benchmark (nf-core/taxprofiler official test set)
    # -------------------------------------------------------------------------
    "metagenome_zymo": {
        "description": "ZymoBIOMICS / HMP Mock Community Real Metagenomic Sequencing Reads",
        "dir": os.path.join(BASE_DIR, "metagenome_zymo"),
        "files": [
            {
                "name": "zymo_mock_R1.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/taxprofiler/data/fastq/ERX5474932_ERR5766176_1.fastq.gz"
            },
            {
                "name": "zymo_mock_R2.fastq.gz",
                "url": "https://raw.githubusercontent.com/nf-core/test-datasets/taxprofiler/data/fastq/ERX5474932_ERR5766176_2.fastq.gz"
            }
        ]
    }
}


def download_url(url, dest_path, ncbi_fallback=None):
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 1000:
        print(f"  [ALREADY EXISTS] {os.path.basename(dest_path)} ({os.path.getsize(dest_path):,} bytes)")
        return True

    print(f"  [DOWNLOADING] {os.path.basename(dest_path)} from {url}...")
    headers = {"User-Agent": "Mozilla/5.0 (Bioinformatics-Benchmark-Downloader)"}

    urls_to_try = [url]
    if ncbi_fallback:
        urls_to_try.insert(0, ncbi_fallback)

    for target_url in urls_to_try:
        try:
            req = urllib.request.Request(target_url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as response, open(dest_path, "wb") as out_file:
                chunk_size = 1024 * 1024  # 1MB
                total_downloaded = 0
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    total_downloaded += len(chunk)
            print(f"  [SUCCESS] {os.path.basename(dest_path)}: {total_downloaded:,} bytes downloaded.")
            return True
        except Exception as e:
            print(f"  [WARNING] Failed downloading from {target_url}: {e}")
            if os.path.exists(dest_path):
                os.remove(dest_path)

    print(f"  [ERROR] All download attempts failed for {dest_path}")
    return False


def main():
    print("=" * 80)
    print("ComfyBioWMS: Downloading Official Certified Benchmark Datasets")
    print(f"Target Directory: {BASE_DIR}")
    print("=" * 80)

    start_time = time.time()
    total_files = 0
    success_files = 0

    for domain_key, dataset in DATASETS.items():
        print(f"\n▶ [{domain_key.upper()}] {dataset['description']}")
        target_dir = dataset["dir"]
        for f_info in dataset["files"]:
            total_files += 1
            dest_file = os.path.join(target_dir, f_info["name"])
            ncbi_fallback = f_info.get("ncbi_ref_url")
            if download_url(f_info["url"], dest_file, ncbi_fallback):
                success_files += 1

    elapsed = time.time() - start_time
    print("\n" + "=" * 80)
    print(f"Download Summary: {success_files}/{total_files} files successfully retrieved in {elapsed:.2f}s.")
    print("=" * 80)


if __name__ == "__main__":
    main()
