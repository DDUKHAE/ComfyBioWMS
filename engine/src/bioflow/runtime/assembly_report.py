import argparse
from pathlib import Path

from bioflow.runtime.report import domain_report_markdown


def write_assembly_report(qc_dir: Path, plot_dir: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        domain_report_markdown(
            "ComfyBIO Genome Assembly Report",
            "QUAST Assembly Quality Metrics",
            [("QUAST output directory", qc_dir), ("Assembly summary plot", Path(plot_dir) / "assembly_summary.png")],
        ),
        encoding="utf-8",
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a ComfyBIO genome assembly markdown report.")
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_assembly_report(args.qc_dir, args.plot_dir, args.output)


if __name__ == "__main__":
    main()
