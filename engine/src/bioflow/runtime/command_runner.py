import subprocess
import shlex
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class CommandRecord:
    argv: list[str]
    cwd: Path
    dry_run: bool
    returncode: int | None = None
    stdout: str = ""
    stderr: str = ""


@dataclass
class DryRunCommandRunner:
    commands: list[CommandRecord] = field(default_factory=list)

    def run(self, argv: list[str], cwd: Path) -> CommandRecord:
        record = CommandRecord(argv=argv, cwd=cwd, dry_run=True)
        self.commands.append(record)
        return record


@dataclass
class CondaCommandRunner:
    commands: list[CommandRecord] = field(default_factory=list)

    def run(self, argv: list[str], cwd: Path) -> CommandRecord:
        completed = subprocess.run(
            argv,
            check=False,
            cwd=cwd,
            capture_output=True,
            text=True,
        )
        record = CommandRecord(
            argv=argv,
            cwd=cwd,
            dry_run=False,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
        self.commands.append(record)
        if completed.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {completed.returncode}: {' '.join(argv)}\n{completed.stderr}")
        return record


def conda_command(env_name: str, executable: str, *args: str) -> list[str]:
    return ["conda", "run", "-n", env_name, executable, *args]


def parse_extra_command_tokens(extra_command: str) -> list[str]:
    tokens: list[str] = []
    for raw_line in extra_command.splitlines() or [extra_command]:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        tokens.extend(shlex.split(line))
    return tokens


def fastp_trim_argv(env_name: str, sample, sample_output_dir, threads, extra_command="") -> list[str]:
    out = Path(sample_output_dir)
    args = ["-i", str(sample.fastq_1)]
    if sample.fastq_2 is not None:
        args += ["-I", str(sample.fastq_2)]
    args += ["--out1", str(out / "R1.fastq")]
    if sample.fastq_2 is not None:
        args += ["--out2", str(out / "R2.fastq")]
    args += ["-w", str(threads), "-j", str(out / "fastp.json"), "-h", str(out / "fastp.html")]
    return conda_command(env_name, "fastp", *args, *parse_extra_command_tokens(extra_command))
