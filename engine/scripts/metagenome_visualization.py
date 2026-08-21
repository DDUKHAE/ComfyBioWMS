import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_bracken_report(report_path: Path, top_n: int = 5) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    with report_path.open(encoding="utf-8") as handle:
        first_line = handle.readline().rstrip("\n")
        if not first_line:
            return [("No data", 1.0)]
        
        fields = first_line.split("\t")
        if "name" in fields and "fraction_total_reads" in fields:
            # Standard Bracken TSV output
            name_idx = fields.index("name")
            fraction_idx = fields.index("fraction_total_reads")
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(name_idx, fraction_idx):
                    try:
                        rows.append((f[name_idx].strip(), float(f[fraction_idx])))
                    except ValueError:
                        continue
        else:
            # Kraken2 / Kraken report format (6 columns: %reads, count, direct_count, rank, taxid, name)
            def parse_kraken_line(line_str):
                f = line_str.rstrip("\n").split("\t")
                if len(f) >= 6:
                    try:
                        pct = float(f[0].strip()) / 100.0
                        name = f[5].strip()
                        return (name, pct)
                    except ValueError:
                        pass
                return None

            first_parsed = parse_kraken_line(first_line)
            if first_parsed:
                rows.append(first_parsed)
            for line in handle:
                p = parse_kraken_line(line)
                if p:
                    rows.append(p)

    if not rows:
        return [("Unclassified", 1.0)]

    # Filter out pure 'root' or 'cellular organisms' if specific species present
    filtered = [r for r in rows if r[0] not in ("root", "cellular organisms")]
    if filtered:
        rows = filtered

    rows.sort(key=lambda row: row[1], reverse=True)
    return rows[:top_n]


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot top taxa by abundance fraction per sample from Bracken/Kraken reports.")
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()

    report_files = sorted(args.reports_dir.glob("*/bracken_output.txt"))
    if not report_files:
        report_files = sorted(args.reports_dir.glob("*/kraken2_report.txt"))
    if not report_files and (args.reports_dir / "bracken_output.txt").exists():
        report_files = [args.reports_dir / "bracken_output.txt"]
    if not report_files and (args.reports_dir / "kraken2_report.txt").exists():
        report_files = [args.reports_dir / "kraken2_report.txt"]
    if not report_files:
        raise SystemExit(f"No bracken_output.txt or kraken2_report.txt files found under {args.reports_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(report_files), figsize=(max(6 * len(report_files), 6), 4), squeeze=False)
    for ax, report_path in zip(axes[0], report_files):
        sample_name = report_path.parent.name
        top_taxa = parse_bracken_report(report_path, args.top_n)
        names = [name for name, _ in top_taxa]
        fractions = [fraction for _, fraction in top_taxa]
        ax.barh(names, fractions, color="#3b82f6")
        ax.set_title(f"Sample: {sample_name}")
        ax.set_xlabel("Relative Abundance")
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Saved metagenome summary plot to {args.output}")


if __name__ == "__main__":
    main()
