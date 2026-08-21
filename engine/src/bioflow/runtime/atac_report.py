import argparse
import sys
from pathlib import Path

# Ensure engine/src is on sys.path for standalone execution
_ENGINE_SRC = Path(__file__).resolve().parents[2]
if str(_ENGINE_SRC) not in sys.path:
    sys.path.insert(0, str(_ENGINE_SRC))

from bioflow.runtime.report import domain_report_markdown


def write_atac_report(peaks_dir: Path, plot_dir: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        domain_report_markdown(
            "ComfyBIO ATAC-seq Report",
            "Called Peaks",
            [("Peaks directory", peaks_dir), ("Peak summary plot", Path(plot_dir) / "atac_summary.png")],
        ),
        encoding="utf-8",
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a ComfyBIO ATAC-seq markdown report.")
    parser.add_argument("--peaks-dir", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_atac_report(args.peaks_dir, args.plot_dir, args.output)


if __name__ == "__main__":
    main()
