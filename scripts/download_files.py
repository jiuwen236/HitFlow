#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Set
sys.path.insert(0, str(Path(__file__).parent))
from snakemake.merge_yaml_duplicates import load_yaml_with_merged_duplicates

data_dir_name = "ref_data"

def parse_args():
    parser = argparse.ArgumentParser(description="Download reference files based on config.yaml settings")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Show files to download without downloading them")
    parser.add_argument("-r", "--resume", action="store_true", help="Resume interrupted downloads using curl -C -")
    parser.add_argument("--all", action="store_true", help="Download all reference files")
    parser.add_argument("--proxyHost", type=str, help="Proxy host")
    parser.add_argument("--proxyPort", type=str, help="Proxy port")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config.yaml")
    parser.add_argument("--cnvpytor", action="store_true", dest="cnvpytor_only", help="Only download CNVpytor resources")
    return parser.parse_args()


def setup_proxy(proxy_host: str, proxy_port: str):
    if proxy_port and not proxy_host:
        proxy_host = "127.0.0.1"
    if proxy_host and proxy_port:
        proxy_url = f"http://{proxy_host}:{proxy_port}"
        os.environ["http_proxy"] = proxy_url
        os.environ["https_proxy"] = proxy_url
        os.environ["HTTP_PROXY"] = proxy_url
        os.environ["HTTPS_PROXY"] = proxy_url


def get_required_files(config: Dict, download_all: bool) -> Set[str]:
    required = set()
    
    workflows = config.get("workflows", {})
    enable_ngs = workflows.get("enable_basic_ngs_workflows", False)
    enable_ngs_single = workflows.get("enable_ngs_single_sample_filtering", False)
    enable_ngs_sv = workflows.get("enable_ngs_sv_detection", False)
    enable_ngs_joint = workflows.get("enable_ngs_joint_calling", False)
    enable_tgs_sv = workflows.get("enable_tgs_sv_detection", False)
    enable_tgs_var = workflows.get("enable_tgs_variant_calling", False)
    enable_tgs_anno = workflows.get("enable_tgs_variant_annotation", False)
    enable_tgs_joint = workflows.get("enable_tgs_joint_calling", False)
    ngs_sv_config = config.get("ngs_sv", {})
    cnvpytor_enabled = ngs_sv_config.get("cnvpytor", {}).get("enabled", False)
    
    if download_all or enable_ngs or enable_ngs_single or enable_ngs_sv or enable_ngs_joint or enable_tgs_sv or enable_tgs_var or enable_tgs_anno or enable_tgs_joint:
        required.add("reference")
    
    if download_all or enable_ngs or enable_ngs_joint or enable_tgs_joint:
        required.add("gatk_dbsnp")
    
    if download_all or enable_ngs:
        required.add("verifybamid")
    
    if download_all or (enable_ngs_joint and config.get("vqsr", {}).get("enabled", False)):
        required.add("vqsr")
    
    if download_all or (enable_ngs_joint and config.get("phasing", {}).get("enabled", False)):
        required.add("phasing")
    
    if download_all or enable_ngs_sv:
        required.add("annotsv")
    if download_all or (enable_ngs_sv and cnvpytor_enabled):
        required.add("cnvpytor")
    annotation_cfg = config.get("ngs_single_sample_variants", {}).get("annotation", {})
    annotation_tool = str(annotation_cfg.get("tool", "")).strip().lower()
    if download_all or (annotation_tool == "snpeff" and (enable_ngs or enable_ngs_single)) or enable_tgs_anno:
        required.add("snpeff")
    if download_all or (annotation_tool == "funcotator" and (enable_ngs or enable_ngs_single)):
        required.add("funcotator")
    
    return required


def download_file(url: str, path: str, resume: bool = False) -> bool:
    parent_dir = os.path.dirname(path)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)

    local_size = os.path.getsize(path) if os.path.exists(path) else 0

    if not url:
        print(f"manual: {path}")
        return False

    check_cmd = ["curl", "-s", "-I", "-L", url]
    result = subprocess.run(check_cmd, capture_output=True, text=True)
    status_line = result.stdout.splitlines()[0] if result.stdout else ""
    # Accept common successful/redirect/resume status codes.
    if not any(code in status_line for code in ("200", "206", "302", "303", "307", "308")):
        print(f"badurl: {path}, {url}")
        return False

    remote_size = None
    for line in result.stdout.splitlines():
        if line.lower().startswith("content-length"):
            try:
                remote_size = int(line.split(":", 1)[1].strip())
            except (IndexError, ValueError):
                remote_size = None
            break

    if local_size:
        if remote_size and local_size >= remote_size:
            print(f"skip: {path}")
            return True
        if not resume:
            if remote_size:
                print(f"incomplete: {path} ({local_size}/{remote_size} bytes), re-downloading")
            else:
                print(f"skip: {path}")
                return True
        else:
            display_total = remote_size if remote_size is not None else "?"
            print(f"resume: {path} ({local_size}/{display_total} bytes)")

    download_cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "5",
        "--connect-timeout",
        "30",
        "--max-time",
        "0",
    ]
    if resume:
        download_cmd.extend(["-C", "-"])
    download_cmd.extend(["-o", path, url])

    result = subprocess.run(download_cmd)
    if result.returncode != 0:
        print(f"fail: {path}, {url}")
        return False

    final_size = os.path.getsize(path) if os.path.exists(path) else 0
    if remote_size and final_size < remote_size:
        print(f"fail: {path} (downloaded {final_size}/{remote_size} bytes)")
        return False

    print(f"done: {path}")
    return True


def download_reference_files(required: Set[str], data_dir: str, dry_run: bool = False, resume: bool = False):
    if not dry_run:
        os.makedirs(data_dir, exist_ok=True)
    prev_cwd = os.getcwd()
    try:
        if not dry_run:
            os.chdir(data_dir)
        base_ftp = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"
        base_verifybamid = "https://raw.githubusercontent.com/Griffan/VerifyBamID/master/resource"
        base_gatk = "http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/"
        base_hg38 = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
        base_gmap = "https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/maps/b38"

        check_cmd = ["curl", "-s", "--head", base_hg38]
        result = subprocess.run(check_cmd, capture_output=True, text=True)
        status_line = result.stdout.splitlines()[0] if result.stdout else ""
        if "200" not in status_line:
            base_hg38 = base_gatk

        files = {}

        if "reference" in required:
            files.update({
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt",
                "ref/grch38/GRCh38_full_analysis_set_plus_decoy_hla.dict": f"{base_ftp}/GRCh38_full_analysis_set_plus_decoy_hla.dict",
            })

        if "verifybamid" in required:
            files.update({
                "ref/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.bed": f"{base_verifybamid}/1000g.phase3.100k.b38.vcf.gz.dat.bed",
                "ref/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu": f"{base_verifybamid}/1000g.phase3.100k.b38.vcf.gz.dat.mu",
                "ref/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD": f"{base_verifybamid}/1000g.phase3.100k.b38.vcf.gz.dat.UD",
                "ref/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.V": f"{base_verifybamid}/1000g.phase3.100k.b38.vcf.gz.dat.V",
            })

        if "gatk_dbsnp" in required:
            files.update({
                "hg38_GATK/dbsnp_138.hg38.vcf.gz": f"{base_gatk}/dbsnp_138.hg38.vcf.gz",
                "hg38_GATK/dbsnp_138.hg38.vcf.gz.tbi": f"{base_gatk}/dbsnp_138.hg38.vcf.gz.tbi",
            })

        if "vqsr" in required:
            files.update({
                "vqsr_resources/hapmap_3.3.hg38.vcf.gz": f"{base_hg38}/hapmap_3.3.hg38.vcf.gz",
                "vqsr_resources/hapmap_3.3.hg38.vcf.gz.tbi": f"{base_hg38}/hapmap_3.3.hg38.vcf.gz.tbi",
                "vqsr_resources/1000G_omni2.5.hg38.vcf.gz": f"{base_hg38}/1000G_omni2.5.hg38.vcf.gz",
                "vqsr_resources/1000G_omni2.5.hg38.vcf.gz.tbi": f"{base_hg38}/1000G_omni2.5.hg38.vcf.gz.tbi",
                "vqsr_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz": f"{base_hg38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "vqsr_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi": f"{base_hg38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
                "vqsr_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz": f"{base_hg38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "vqsr_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi": f"{base_hg38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
            })

        if "phasing" in required:
            for chrom in list(range(1, 23)) + ["X"]:
                files[f"phasing_resources/chr{chrom}.b38.gmap.gz"] = f"{base_gmap}/chr{chrom}.b38.gmap.gz"

        if dry_run:
            if files:
                print("Dry run - the following files would be downloaded:")
            else:
                print("Dry run - no files need to be downloaded.")
            for path, url in sorted(files.items()):
                full_path = os.path.join(data_dir, path)
                if not os.path.exists(full_path):
                    print(f"missing: {path} <- {url}")
            return

        for path, url in files.items():
            download_file(url, path, resume=resume)
    finally:
        os.chdir(prev_cwd)


def download_snpeff_database(project_root: str, dry_run: bool = False, resume: bool = False):
    root = os.fspath(project_root)
    snpeff_dir = os.path.join(root, data_dir_name, "snpEff")
    zip_name = "snpEff_v5_0_GRCh38.99.zip"
    zip_path = os.path.join(snpeff_dir, zip_name)
    genome_dir = os.path.join(snpeff_dir, "data", "GRCh38.99")
    url = "https://snpeff.odsp.astrazeneca.com/databases/v5_0/snpEff_v5_0_GRCh38.99.zip"

    rel_zip = os.path.relpath(zip_path, root)
    rel_genome_dir = os.path.relpath(genome_dir, root)

    if dry_run:
        if not os.path.isdir(genome_dir):
            print(f"missing: {rel_zip} <- {url}")
        return

    os.makedirs(snpeff_dir, exist_ok=True)

    if not download_file(url, zip_path, resume=resume):
        return

    if not os.path.isdir(genome_dir):
        try:
            shutil.unpack_archive(zip_path, snpeff_dir)
            print(f"unzipped: {rel_genome_dir}")
        except (shutil.ReadError, ValueError) as exc:
            print(f"fail_unzip: {rel_zip} ({exc})")
    else:
        print(f"skip_unzip: {rel_genome_dir}")


def download_funcotator_datasources(project_root: str, config: Dict, dry_run: bool = False):
    """Download Funcotator germline data sources using GATK."""
    root = os.fspath(project_root)
    annotation_cfg = config.get("ngs_single_sample_variants", {}).get("annotation", {})
    funcotator_cfg = annotation_cfg.get("funcotator", {})
    data_sources_path = funcotator_cfg.get("data_sources", "ref_data/funcotator_dataSources/funcotator_dataSources.v1.8.hg38.20230908g")
    
    # Convert relative path to absolute
    if not os.path.isabs(data_sources_path):
        data_sources_path = os.path.join(root, data_sources_path)
    
    # Check if data sources already exist
    if os.path.isdir(data_sources_path):
        print(f"skip: Funcotator data sources already exist at {os.path.relpath(data_sources_path, root)}")
        return
    
    if dry_run:
        print(f"missing: Funcotator data sources at {os.path.relpath(data_sources_path, root)}")
        print("Would download using GATK FuncotatorDataSourceDownloader")
        return
    
    print("Downloading Funcotator germline data sources...")
    
    # Prepare output directory
    output_dir = os.path.dirname(data_sources_path)
    os.makedirs(output_dir, exist_ok=True)
    
    # Check for proxy settings
    proxy_host = os.environ.get("http_proxy") or os.environ.get("HTTP_PROXY", "")
    proxy_opts = ""
    if proxy_host:
        # Parse proxy from http://host:port format
        import re
        match = re.search(r'http://([^:]+):(\d+)', proxy_host)
        if match:
            host, port = match.groups()
            proxy_opts = f"-Dhttp.proxyHost={host} -Dhttp.proxyPort={port} -Dhttps.proxyHost={host} -Dhttps.proxyPort={port}"
    
    # Build GATK command
    log_file = os.path.join(root, "FuncotatorDataSourceDownloader.log")
    java_opts = f'"{proxy_opts}"' if proxy_opts else '""'
    
    cmd = f"""
    conda run -n gatk4-tmp gatk --java-options {java_opts} \
      FuncotatorDataSourceDownloader \
      --germline \
      --hg38 \
      --validate-integrity \
      --extract-after-download \
      -O {output_dir}/data \
      > {log_file} 2>&1
    """
    
    # Check if gatk4-tmp conda environment exists, if not create it
    check_env = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True,
        text=True
    )
    
    if "gatk4-tmp" not in check_env.stdout:
        print("Creating gatk4-tmp conda environment...")
        create_env = subprocess.run(
            ["conda", "create", "-n", "gatk4-tmp", "-c", "bioconda", "-c", "conda-forge", "gatk4=4.6.2.0", "-y"],
            cwd=root
        )
        if create_env.returncode != 0:
            print("fail: Could not create gatk4-tmp conda environment")
            return
    
    print(f"Running GATK FuncotatorDataSourceDownloader (check {os.path.relpath(log_file, root)} for progress)...")
    result = subprocess.run(cmd, shell=True, executable="/bin/bash", cwd=root)
    
    if result.returncode == 0 and os.path.isdir(data_sources_path):
        print(f"done: Funcotator data sources downloaded to {os.path.relpath(data_sources_path, root)}")
    else:
        print(f"fail: Funcotator data source download failed. Check {os.path.relpath(log_file, root)} for details")


# Matches the default AnnotSV container (quay.io/biocontainers/annotsv:3.5.3).
ANNOTSV_VERSION = "3.5.3"


def run_annotsv_install(data_dir: str, dry_run: bool = False):
    annotsv_path = os.path.join(data_dir, "AnnotSV_annotations/Annotations_Human")
    if not os.path.isdir(annotsv_path):
        script_path = os.path.join(os.path.dirname(data_dir), "scripts/pipeline/install_annotsv_annotations.sh")
        if dry_run:
            print(f"Dry run - would run AnnotSV installation script: {script_path}")
            return
        subprocess.run(["bash", script_path, ANNOTSV_VERSION, data_dir], check=True)


def run_cnvpytor_download(project_root: str):
    print("Downloading CNVpytor resources...")

    def _run() -> subprocess.CompletedProcess:
        if shutil.which("snakemake"):
            return subprocess.run(["snakemake", "cnvpytor_download_resources", "--use-conda"], cwd=project_root)
        command = "source ~/.bashrc && conda activate snakemake && snakemake cnvpytor_download_resources --use-conda --cores 4"
        return subprocess.run(command, shell=True, executable="/bin/bash", cwd=project_root)

    result = _run()
    if result.returncode != 0:
        print("CNVpytor resource download failed (1). Try to fix PosixPath bug.")

    cnvpytor_log = os.path.join(project_root, "logs/ngs_sv/cnvpytor_download.log")
    if os.path.exists(cnvpytor_log):
        with open(cnvpytor_log) as f:
            log_content = f.read()
        if "'PosixPath' object has no attribute 'split'" in log_content:
            import re
            match = re.search(r'File "(.+genome\.py)"', log_content)
            if match:
                genome_py = match.group(1)
                if os.path.exists(genome_py):
                    print(f"Fixing CNVpytor PosixPath bug in: {genome_py}")
                    with open(genome_py) as f:
                        content = f.read()
                    content = re.sub(r'fn = res\.split\("/"\)[^]]*\]', 'fn = res.name', content)
                    with open(genome_py, "w") as f:
                        f.write(content)
                    _run()


def hitsv_enabled(config: Dict) -> bool:
    # Mirrors the SRS/LRS/Hybrid switches in rules/hitsv.smk.
    workflows = config.get("workflows", {})
    srs = workflows.get("enable_ngs_sv_detection", False) and config.get("ngs_sv", {}).get("hitsv", {}).get("enabled", False)
    lrs = workflows.get("enable_tgs_sv_detection", False) and config.get("tgs_config", {}).get("sv_detection", {}).get("hitsv", {}).get("enabled", False)
    return bool(srs or lrs or config.get("hitsv", {}).get("hybrid", {}).get("enabled", False))


def run_import_script(script_name: str, project_root: str):
    script_path = os.path.join(project_root, "scripts", script_name)
    if os.path.exists(script_path):
        print(f"Running {script_name}...")
        subprocess.run(["bash", script_path], cwd=project_root)


def download_1000g_phased_panel(project_root: str, dry_run: bool = False, resume: bool = False):
    """
    Download 1000G high coverage phased reference panel for SHAPEIT5.
    Files are downloaded from:
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
    """
    base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
    output_dir = os.path.join(project_root, data_dir_name, "phasing_resources", "1000g_phased")
    
    def _ensure_symlink(link_path: str, target_path: str) -> None:
        """Create/refresh a symlink at link_path -> target_path (relative when possible)."""
        if link_path == target_path:
            return
        try:
            # If existing and correct, keep it.
            if os.path.islink(link_path):
                if os.readlink(link_path) == os.path.relpath(target_path, os.path.dirname(link_path)):
                    return
                os.unlink(link_path)
            elif os.path.exists(link_path):
                # Don't overwrite real files.
                return

            rel_target = os.path.relpath(target_path, os.path.dirname(link_path))
            os.symlink(rel_target, link_path)
        except OSError:
            # Best-effort; symlinks may be unsupported on some FS.
            return

    # List of files to download (non-recursive, main directory only)
    # chr1-22: 1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    # chrX:    1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
    chromosomes = list(range(1, 23)) + ["X"]
    
    files = {}
    for chrom in chromosomes:
        if str(chrom).upper() == "X":
            vcf_name = "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
        else:
            vcf_name = f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        tbi_name = f"{vcf_name}.tbi"
        files[vcf_name] = f"{base_url}/{vcf_name}"
        files[tbi_name] = f"{base_url}/{tbi_name}"
    
    # Additional files in the directory
    additional_files = [
        "README_1kGP_phased_panel_110722.pdf",
    ]
    for fname in additional_files:
        files[fname] = f"{base_url}/{fname}"
    
    if dry_run:
        print(f"Dry run - would download {len(files)} files to {os.path.relpath(output_dir, project_root)}")
        for fname, url in sorted(files.items()):
            full_path = os.path.join(output_dir, fname)
            if not os.path.exists(full_path):
                print(f"missing: {fname}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    prev_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        print(f"Downloading 1000G phased reference panel to {output_dir}...")
        success_count = 0
        downloaded: Set[str] = set()
        for fname, url in sorted(files.items()):
            if download_file(url, fname, resume=resume):
                success_count += 1
                downloaded.add(fname)

        # Create compatibility symlinks for chrX v2 file (without ".v2").
        v2_vcf = "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
        v2_tbi = f"{v2_vcf}.tbi"
        legacy_vcf = "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        legacy_tbi = f"{legacy_vcf}.tbi"
        if os.path.exists(v2_vcf):
            _ensure_symlink(legacy_vcf, v2_vcf)
        if os.path.exists(v2_tbi):
            _ensure_symlink(legacy_tbi, v2_tbi)

        print(f"Downloaded {success_count}/{len(files)} files.")
    finally:
        os.chdir(prev_cwd)


def main():
    args = parse_args()
    
    setup_proxy(args.proxyHost, args.proxyPort)
    
    project_root = Path(__file__).parent.parent.absolute()

    if args.cnvpytor_only:
        if args.dry_run:
            print("Dry run - can not judge CNVpytor download. Please run without --dry-run.")
        else:
            run_cnvpytor_download(project_root)
        return
    
    config_path = os.path.join(project_root, args.config)
    
    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)
    
    config = load_yaml_with_merged_duplicates(config_path)
    required_files = get_required_files(config, args.all)
    phasing_cfg = config.get("phasing", {})
    phasing_enabled = phasing_cfg.get("enabled", False)
    use_reference_panel = phasing_cfg.get("use_reference", False)
    download_1000g_panel = args.all or (phasing_enabled and use_reference_panel)
    
    data_dir = os.path.join(project_root, data_dir_name)
    download_reference_files(required_files, data_dir, dry_run=args.dry_run, resume=args.resume)
    if "snpeff" in required_files:
        download_snpeff_database(project_root, dry_run=args.dry_run, resume=args.resume)
    if "funcotator" in required_files:
        download_funcotator_datasources(project_root, config, dry_run=args.dry_run)
    
    if "annotsv" in required_files:
        run_annotsv_install(data_dir, dry_run=args.dry_run)

    if download_1000g_panel:
        download_1000g_phased_panel(project_root, dry_run=args.dry_run, resume=args.resume)
    
    if args.dry_run:
        print("Dry run completed.")
        return

    workflows = config.get("workflows", {})
    enable_ngs_sv = workflows.get("enable_ngs_sv_detection", False)
    enable_ngs_joint = workflows.get("enable_ngs_joint_calling", False)
    
    if args.all or enable_ngs_sv:
        ngs_sv = config.get("ngs_sv", {})
        if ngs_sv.get("manta", {}).get("enabled", False):
            run_import_script("pipeline/import_manta.sh", project_root)
    
    if args.all or (enable_ngs_joint and phasing_enabled):
        run_import_script("pipeline/import_shapeit5.sh", project_root)

    # Installs the pinned binary and the Hybrid -k error profile.
    if args.all or hitsv_enabled(config):
        run_import_script("pipeline/install_hitsv.sh", project_root)
    
    print("Download process completed.")


if __name__ == "__main__":
    main()
