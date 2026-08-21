import argparse
from pathlib import Path


def domain_report_markdown(title: str, section_heading: str, entries: list[tuple[str, Path]]) -> str:
    lines = [f"# {title}", "", f"## {section_heading}", ""]
    lines += [f"- {label}: `{path}`" for label, path in entries]
    return "\n".join(lines) + "\n"


def report_markdown(results_csv: Path, plot_dir: Path) -> str:
    return (
        "# ComfyBIO Report\n\n"
        "## DESeq2 Results\n\n"
        f"- Results table: `{results_csv}`\n"
        f"- PCA: `{plot_dir / 'pca.png'}`\n"
        f"- MA plot: `{plot_dir / 'ma.png'}`\n"
        f"- Volcano: `{plot_dir / 'volcano.png'}`\n"
        f"- Heatmap: `{plot_dir / 'heatmap.png'}`\n"
    )


def write_report(results_csv: Path, plot_dir: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report_markdown(results_csv, plot_dir), encoding="utf-8")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a ComfyBIO markdown report.")
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_report(args.results, args.plot_dir, args.output)


if __name__ == "__main__":
    main()
