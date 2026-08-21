import argparse
from pathlib import Path

from bioflow.runtime.report import domain_report_markdown


def write_variant_report(vcf_dir: Path, plot_dir: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        domain_report_markdown(
            "ComfyBIO Variant Analysis Report",
            "Filtered Variants",
            [("Filtered VCF directory", vcf_dir), ("Variant summary plot", Path(plot_dir) / "variant_summary.png")],
        ),
        encoding="utf-8",
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a ComfyBIO variant analysis markdown report.")
    parser.add_argument("--vcf-dir", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_variant_report(args.vcf_dir, args.plot_dir, args.output)


if __name__ == "__main__":
    main()
