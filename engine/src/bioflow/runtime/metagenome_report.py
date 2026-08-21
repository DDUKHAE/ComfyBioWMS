import argparse
import sys
from pathlib import Path

# Ensure engine/src is on sys.path for standalone execution
_ENGINE_SRC = Path(__file__).resolve().parents[2]
if str(_ENGINE_SRC) not in sys.path:
    sys.path.insert(0, str(_ENGINE_SRC))

from bioflow.runtime.report import domain_report_markdown


def write_metagenome_report(bracken_dir: Path, plot_dir: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        domain_report_markdown(
            "ComfyBIO Metagenome Taxonomic Profiling Report",
            "Bracken Abundance Estimates",
            [("Bracken output directory", bracken_dir), ("Taxonomic summary plot", Path(plot_dir) / "metagenome_summary.png")],
        ),
        encoding="utf-8",
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a ComfyBIO metagenome taxonomic profiling markdown report.")
    parser.add_argument("--bracken-dir", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_metagenome_report(args.bracken_dir, args.plot_dir, args.output)


if __name__ == "__main__":
    main()
