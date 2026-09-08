import datetime
import json
import os
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_TOOL_ENV_MAPPINGS: dict[str, str] = {
    "fastp": "bulk_rna_seq",
    "salmon": "bulk_rna_seq",
    "Rscript": "bulk_rna_seq",
    "STAR": "bulk_rna_seq",
    "trimmomatic": "bulk_rna_seq",
    "fastqc": "bulk_rna_seq",
    "samtools": "variant_analysis",
    "bcftools": "variant_analysis",
    "bwa": "genome_assembly",
    "bwa-mem2": "variant_analysis",
    "spades": "genome_assembly",
    "spades.py": "genome_assembly",
    "quast": "genome_assembly",
    "quast.py": "genome_assembly",
    "kraken2": "metagenome",
    "bracken": "metagenome",
    "macs3": "epigenomics",
    "bedtools": "epigenomics",
    "flye": "genome_assembly",
}


def get_conda_base() -> Path | None:
    """Find the base directory for conda / miniconda installation."""
    if "CONDA_PREFIX" in os.environ:
        base = Path(os.environ["CONDA_PREFIX"])
        if (base / "envs").exists():
            return base
        if base.parent.name == "envs":
            return base.parent.parent

    for candidate in [
        Path("/opt/miniconda3"),
        Path("/opt/anaconda3"),
        Path.home() / "miniconda3",
        Path.home() / "anaconda3",
    ]:
        if candidate.exists() and (candidate / "envs").exists():
            return candidate

    conda_which = shutil.which("conda")
    if conda_which:
        p = Path(conda_which).resolve()
        # usually .../bin/conda -> base is parent of bin
        base = p.parent.parent
        if (base / "envs").exists():
            return base

    return None


def resolve_tool_environment(tool: str, requested_env: str | None = None) -> tuple[str, str]:
    """Resolve the execution environment and executable path for a bioinformatics tool.

    Returns:
        (env_name, executable_path)

    Raises:
        EnvironmentError: If requested_env is specified but does not exist.
        FileNotFoundError: If the tool cannot be found in the specified or canonical environments.
    """
    conda_base = get_conda_base()

    if requested_env:
        # User/Workflow explicitly requested a specific conda environment
        if not conda_base or not (conda_base / "envs" / requested_env).exists():
            raise EnvironmentError(
                f"Requested conda environment '{requested_env}' does not exist. "
                f"Silent fallback to PATH is prohibited for reproducibility."
            )
        env_dir = conda_base / "envs" / requested_env
        bin_path = env_dir / "bin" / tool
        if bin_path.exists() and os.access(bin_path, os.X_OK):
            return requested_env, str(bin_path)
        # Check alternative common extensions or names (e.g. spades -> spades.py)
        if (env_dir / "bin" / f"{tool}.py").exists():
            return requested_env, str(env_dir / "bin" / f"{tool}.py")
        raise FileNotFoundError(
            f"Tool '{tool}' not found in conda environment '{requested_env}' at {env_dir / 'bin'}."
        )

    # If not explicitly requested, check canonical domain mapping
    canonical_env = DEFAULT_TOOL_ENV_MAPPINGS.get(tool)
    if canonical_env and conda_base and (conda_base / "envs" / canonical_env).exists():
        env_dir = conda_base / "envs" / canonical_env
        bin_path = env_dir / "bin" / tool
        if bin_path.exists() and os.access(bin_path, os.X_OK):
            return canonical_env, str(bin_path)
        if (env_dir / "bin" / f"{tool}.py").exists():
            return canonical_env, str(env_dir / "bin" / f"{tool}.py")

    # Search other available conda environments if tool is not in canonical env
    if conda_base and (conda_base / "envs").exists():
        for env_path in (conda_base / "envs").iterdir():
            if env_path.is_dir():
                bin_path = env_path / "bin" / tool
                if bin_path.exists() and os.access(bin_path, os.X_OK):
                    return env_path.name, str(bin_path)

    # Finally check system PATH as a fallback if no conda environment provided it
    system_tool = shutil.which(tool)
    if system_tool:
        return "system_path", system_tool

    raise FileNotFoundError(
        f"Tool '{tool}' could not be resolved in any registered conda environment or system PATH."
    )


def validate_extra_command_conflicts(extra_command: str, managed_flags: set[str] | list[str]) -> list[str]:
    """Parse extra_command tokens and check for conflicts with UI-managed flags.

    Raises ValueError if a conflicting option is present.
    """
    tokens = parse_extra_command_tokens(extra_command)
    managed_set = set(managed_flags)
    for token in tokens:
        if token.startswith("-"):
            flag_name = token.split("=")[0]
            if flag_name in managed_set:
                raise ValueError(
                    f"Conflict detected: option '{flag_name}' in extra_command conflicts with UI-managed parameter."
                )
    return tokens


@dataclass(frozen=True)
class CommandRecord:
    argv: list[str]
    cwd: Path
    dry_run: bool
    returncode: int | None = None
    stdout: str = ""
    stderr: str = ""
    environment: str = ""
    duration_seconds: float = 0.0
    stdout_path: Path | None = None
    stderr_path: Path | None = None
    manifest_path: Path | None = None


@dataclass
class DryRunCommandRunner:
    commands: list[CommandRecord] = field(default_factory=list)

    def run(self, argv: list[str], cwd: Path, **kwargs: Any) -> CommandRecord:
        record = CommandRecord(argv=argv, cwd=cwd, dry_run=True)
        self.commands.append(record)
        try:
            self.export_manifest(cwd)
        except Exception:
            pass
        return record

    def export_manifest(self, output_dir: Path) -> tuple[Path, Path]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        sh_path = output_dir / "run_manifest.sh"
        json_path = output_dir / "run_manifest.json"

        sh_lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
        manifest_data = []
        for i, cmd in enumerate(self.commands, 1):
            argv_str = " ".join(shlex.quote(arg) for arg in cmd.argv)
            sh_lines.append(f"# Step {i}: CWD={cmd.cwd}")
            sh_lines.append(f"(cd {shlex.quote(str(cmd.cwd))} && {argv_str})")
            sh_lines.append("")

            manifest_data.append({
                "step": i,
                "argv": cmd.argv,
                "cwd": str(cmd.cwd),
                "dry_run": cmd.dry_run,
                "returncode": cmd.returncode,
            })

        sh_path.write_text("\n".join(sh_lines) + "\n", encoding="utf-8")
        sh_path.chmod(0o755)
        json_path.write_text(json.dumps(manifest_data, indent=2) + "\n", encoding="utf-8")
        return sh_path, json_path


@dataclass
class CondaCommandRunner:
    commands: list[CommandRecord] = field(default_factory=list)

    def run(self, argv: list[str], cwd: Path) -> CommandRecord:
        sub_env = os.environ.copy()
        sub_env["CONDA_NO_PLUGINS"] = "true"
        completed = subprocess.run(
            argv,
            check=False,
            cwd=cwd,
            capture_output=True,
            text=True,
            env=sub_env,
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
            raise RuntimeError(
                f"Command failed with exit code {completed.returncode}: {' '.join(argv)}\n{completed.stderr}"
            )
        try:
            self.export_manifest(cwd)
        except Exception:
            pass
        return record

    def export_manifest(self, output_dir: Path) -> tuple[Path, Path]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        sh_path = output_dir / "run_manifest.sh"
        json_path = output_dir / "run_manifest.json"

        sh_lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
        manifest_data = []
        for i, cmd in enumerate(self.commands, 1):
            argv_str = " ".join(shlex.quote(arg) for arg in cmd.argv)
            sh_lines.append(f"# Step {i}: CWD={cmd.cwd}")
            sh_lines.append(f"(cd {shlex.quote(str(cmd.cwd))} && {argv_str})")
            sh_lines.append("")

            manifest_data.append({
                "step": i,
                "argv": cmd.argv,
                "cwd": str(cmd.cwd),
                "dry_run": cmd.dry_run,
                "returncode": cmd.returncode,
            })

        sh_path.write_text("\n".join(sh_lines) + "\n", encoding="utf-8")
        sh_path.chmod(0o755)
        json_path.write_text(json.dumps(manifest_data, indent=2) + "\n", encoding="utf-8")
        return sh_path, json_path


class BioCommandRunner:
    """High-reliability command runner for bioinformatics CLI tools.

    Provides:
    - Dynamic Conda environment resolution and isolated execution
    - Streaming execution logs directly to disk (preventing memory exhaustion)
    - Comprehensive audit manifest generation (run_manifest.json & run_manifest.sh)
    - Input & output artifact tracking
    """

    @classmethod
    def run(
        cls,
        argv: list[str],
        cwd: Path | str,
        env_name: str | None = None,
        node_type: str = "",
        node_id: str = "",
        workflow_id: str = "",
        inputs: list[str | Path] | None = None,
        outputs: list[str | Path] | None = None,
        extra_command: str = "",
        managed_flags: set[str] | list[str] | None = None,
    ) -> CommandRecord:
        cwd_path = Path(cwd).resolve()
        cwd_path.mkdir(parents=True, exist_ok=True)

        # Validate extra_command against managed flags if provided
        extra_tokens: list[str] = []
        if extra_command.strip():
            extra_tokens = validate_extra_command_conflicts(extra_command, managed_flags or set())

        full_argv = list(argv) + extra_tokens

        # Resolve executable and environment
        tool_name = Path(full_argv[0]).name
        resolved_env, exec_path = resolve_tool_environment(tool_name, env_name)
        full_argv[0] = exec_path

        # Setup environment variables for execution
        exec_env = os.environ.copy()
        exec_bin_dir = str(Path(exec_path).parent)
        exec_env["PATH"] = f"{exec_bin_dir}:{exec_env.get('PATH', '')}"
        exec_env["CONDA_NO_PLUGINS"] = "true"

        # Prepare log files
        stdout_log = cwd_path / "stdout.log"
        stderr_log = cwd_path / "stderr.log"

        # Record inputs metadata
        input_records: list[dict[str, Any]] = []
        if inputs:
            for inp in inputs:
                p = Path(inp).resolve()
                if p.exists():
                    st = p.stat()
                    input_records.append({
                        "path": str(p),
                        "exists": True,
                        "size": st.st_size,
                        "mtime_ns": st.st_mtime_ns,
                    })
                else:
                    input_records.append({"path": str(p), "exists": False})

        start_time_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()
        t0 = time.perf_counter()

        with open(stdout_log, "w", encoding="utf-8") as out_f, open(stderr_log, "w", encoding="utf-8") as err_f:
            proc = subprocess.run(
                full_argv,
                cwd=str(cwd_path),
                stdout=out_f,
                stderr=err_f,
                check=False,
                env=exec_env,
            )

        duration = time.perf_counter() - t0
        end_time_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()
        returncode = proc.returncode

        # Record outputs metadata
        output_records: list[dict[str, Any]] = []
        if outputs:
            for out in outputs:
                p = Path(out).resolve()
                if p.exists():
                    st = p.stat()
                    output_records.append({
                        "path": str(p),
                        "exists": True,
                        "size": st.st_size,
                        "mtime_ns": st.st_mtime_ns,
                    })
                else:
                    output_records.append({"path": str(p), "exists": False})

        # Generate run_manifest.json and run_manifest.sh
        manifest_data = {
            "node_type": node_type or tool_name,
            "node_id": node_id,
            "workflow_id": workflow_id,
            "environment": resolved_env,
            "executable": exec_path,
            "argv": full_argv,
            "cwd": str(cwd_path),
            "start_time": start_time_iso,
            "end_time": end_time_iso,
            "duration_seconds": round(duration, 4),
            "returncode": returncode,
            "status": "SUCCESS" if returncode == 0 else "FAILED",
            "stdout_log": str(stdout_log),
            "stderr_log": str(stderr_log),
            "inputs": input_records,
            "outputs": output_records,
        }

        manifest_json_path = cwd_path / "run_manifest.json"
        manifest_json_path.write_text(json.dumps(manifest_data, indent=2) + "\n", encoding="utf-8")

        manifest_sh_path = cwd_path / "run_manifest.sh"
        sh_lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"# Execution manifest for {node_type or tool_name}",
            f"# Environment: {resolved_env}",
            f"export PATH={shlex.quote(exec_bin_dir)}:$PATH",
            f"cd {shlex.quote(str(cwd_path))}",
            " ".join(shlex.quote(a) for a in full_argv),
            "",
        ]
        manifest_sh_path.write_text("\n".join(sh_lines), encoding="utf-8")
        manifest_sh_path.chmod(0o755)

        record = CommandRecord(
            argv=full_argv,
            cwd=cwd_path,
            dry_run=False,
            returncode=returncode,
            environment=resolved_env,
            duration_seconds=duration,
            stdout_path=stdout_log,
            stderr_path=stderr_log,
            manifest_path=manifest_json_path,
        )

        if returncode != 0:
            # Read last 15 lines of stderr for concise and diagnostic error message
            stderr_tail = ""
            if stderr_log.exists():
                lines = stderr_log.read_text(encoding="utf-8", errors="replace").splitlines()
                stderr_tail = "\n".join(lines[-15:])
            raise RuntimeError(
                f"[{node_type or tool_name}] failed with returncode {returncode}.\n"
                f"Command: {' '.join(full_argv)}\n"
                f"Stderr tail:\n{stderr_tail}\n"
                f"Audit logs saved to: {cwd_path}"
            )

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

