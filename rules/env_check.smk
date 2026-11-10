# Software environment check: run each selected tool's version command through
# the same Conda environment / container used by the analysis rules.
import shutil

# tool: (conda env, containers key, version command)
_ENV_CHECK_REGISTRY = {
    "bcftools": ("../envs/basic.yaml", "bcftools", "bcftools --version"),
    "samtools": ("../envs/basic.yaml", "samtools", "samtools --version"),
    "fastp": ("../envs/ngs/fastp.yaml", "fastp", "fastp --version"),
    "fastqc": ("../envs/ngs/fastqc.yaml", "fastqc", "fastqc --version"),
    "gatk4": ("../envs/ngs/gatk4.yaml", "gatk4", "gatk --version"),
    "multiqc": ("../envs/multiqc.yaml", "multiqc", "multiqc --version"),
    "minimap2": ("../envs/long_read/ont_align.yaml", "minimap2", "minimap2 --version"),
    "sniffles": ("../envs/long_read/ont_sniffles.yaml", "sniffles", "sniffles --version"),
    "cutesv": ("../envs/long_read/ont_cutesv.yaml", "cutesv", "cuteSV --version"),
    "truvari": ("../envs/long_read/ont_truvari.yaml", "truvari", "truvari version"),
}

ENV_CHECK_DIR = f"{config['output']['root']}/env_check"
ENV_CHECK_TOOLS = list(config.get("env_check", {}).get("tools", ["bcftools", "fastp", "minimap2"]))
_unknown = sorted(set(ENV_CHECK_TOOLS) - set(_ENV_CHECK_REGISTRY))
if not ENV_CHECK_TOOLS or _unknown:
    raise ValueError("env_check.tools must be a nonempty subset of: "
                     + ", ".join(_ENV_CHECK_REGISTRY) + (f"; unknown: {', '.join(_unknown)}" if _unknown else ""))

# Requested deployment methods, e.g. "conda", "apptainer" or "apptainer conda".
ENV_CHECK_METHODS = " ".join(sorted(m.name.lower() for m in workflow.deployment_settings.deployment_method))

for _tool in ENV_CHECK_TOOLS:
    _env, _container_key, _command = _ENV_CHECK_REGISTRY[_tool]

    rule:
        name: f"env_check_{_tool}"
        output:
            f"{ENV_CHECK_DIR}/tools/{_tool}.tsv"
        params:
            tool = _tool,
            command = _command,
            binary = _command.split()[0],
            methods = ENV_CHECK_METHODS
        threads: 1
        resources:
            mem_mb = 2048,
            time = "00:30:00"
        conda:
            _env
        container:
            config.get("containers", {}).get(_container_key, "")
        log:
            f"{LOGS_DIR}/env_check/{_tool}.log"
        shell:
            r"""
            exec 2> {log:q}
            if ! version=$({params.command} 2>&1); then
                printf '%s\n' "$version" >&2
                echo "ERROR: '{params.command}' failed" >&2
                exit 1
            fi
            version=$(printf '%s\n' "$version" | sed -n 1p)
            binary=$(command -v {params.binary})
            container=${{APPTAINER_CONTAINER:-${{SINGULARITY_CONTAINER:--}}}}
            case " {params.methods} " in
                *" apptainer "*)
                    if [ "$container" = "-" ]; then
                        echo "ERROR: {params.tool} did not run in a container; check containers.* and apptainer/singularity" >&2
                        exit 1
                    fi ;;
            esac
            case " {params.methods} " in
                *" conda "*)
                    case "$binary" in
                        "${{CONDA_PREFIX:-/nonexistent}}"/*) ;;
                        *) echo "ERROR: $binary is not from the job Conda environment" >&2; exit 1 ;;
                    esac ;;
            esac
            printf '%s\t%s\t%s\t%s\t%s\n' {params.tool} "$version" "$binary" \
                "${{CONDA_PREFIX:--}}" "$container" > {output:q}
            """


rule env_check_report:
    input:
        expand(f"{ENV_CHECK_DIR}/tools/{{tool}}.tsv", tool=ENV_CHECK_TOOLS)
    output:
        f"{ENV_CHECK_DIR}/report.tsv"
    localrule: True
    run:
        from snakemake import __version__ as snakemake_version

        # Host details describe the machine running the main Snakemake process.
        lines = [f"# snakemake\t{snakemake_version}",
                 f"# deployment\t{ENV_CHECK_METHODS or 'none'}"]
        lines += [f"# {cmd}\t{shutil.which(cmd) or '-'}"
                  for cmd in ("conda", "mamba", "apptainer", "singularity", "sbatch")]
        lines.append("tool\tversion\tbinary\tconda_prefix\tcontainer")
        for path in input:
            with open(path) as fh:
                lines.append(fh.read().rstrip("\n"))
        with open(output[0], "w") as fh:
            fh.write("\n".join(lines) + "\n")
