from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def crisp_align(id, sample, fq1, fq2, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-standard-16", work_disk="auto", work_disk_type="pd-balanced", preemptible=False):
    """Start CRISP ALIGN pipeline

    - id (``str``) unique run identifier
    - sample (``str``) SAMPLE in read-group
    - fq1 (``str``) fastq1 GS location
    - fq2 (``str``) fastq2 GS location
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "crisp-align-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        fq_size_gb = size_gcp(fq1, env) + size_gcp(fq2, env)
        work_disk = "%sG" % max(200, round(fq_size_gb * 10))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "SAMPLE":sample,
        "FQ1":fq1,
        "FQ2":fq2,
        "ALIGN_UMI":env["CRISP_ALIGN_UMI"],
        "ALIGN_FASTA":env["CRISP_ALIGN_FASTA"],
        "ALIGN_INDEX":env["CRISP_ALIGN_INDEX"],
        "ALIGN_CUTARGS":env["CRISP_ALIGN_CUTARGS"],
        "ALIGN_MRGARGS":env["CRISP_ALIGN_MRGARGS"],
        "ALIGN_RRNA":env["CRISP_ALIGN_RRNA"],
        "ALIGN_GENOTYPE":env["CRISP_ALIGN_GENOTYPE"],
    }
    check_done(id, "crisp-align", env, overwrite, local)
    with open(env["ROOT"] / "pipe/crisp-align/crisp-align_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/crisp-align/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/crisp-align/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/crisp-align/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/crisp-align/crisp-align_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/crisp-align/crisp-align_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
        
@task
def crisp_germline(id, tsample, alnt, name=None, local=False, nowait=False, overwrite=False,
                   machine="n2-highcpu-32", work_disk="250G", work_disk_type="pd-ssd", preemptible=False):
    """Execute CRISP germline pipeline

    - id (``str``) unique run identifier
    - tsample (``str``) tumor SAMPLE in read-group
    - alnt (``str``) crisp-alignment tumor GS location
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "crisp-germline-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "TSAMPLE":tsample,
        "ALNT":alnt,
        "ALIGN_FASTA":env["CRISP_ALIGN_FASTA"],
        "ALIGN_DBSNP":env["CRISP_ALIGN_DBSNP"],
        "ALIGN_INDEL":env["CRISP_ALIGN_INDEL"],
        "GERMLINE_DNASCOPEARGS":env["CRISP_GERMLINE_DNASCOPEARGS"],
        "GERMLINE_DBSNP":env["CRISP_GERMLINE_DBSNP"],
        "GERMLINE_DNASCOPE":env["CRISP_GERMLINE_DNASCOPE"],
    }
    check_done(id, "crisp-germline", env, overwrite, local)
    with open(env["ROOT"] / "pipe/crisp-germline/crisp-germline_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/crisp-germline/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/crisp-germline/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/crisp-germline/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/crisp-germline/crisp-germline_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/crisp-germline/crisp-germline_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def crisp_quant(id, fq1s, fq2s, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-standard-16", work_disk="200G", work_disk_type="pd-balanced", preemptible=False):
    """Execute CRISP QUANT pipeline

    - id (``str``) unique run identifier
    - fq1s (``str``) fastq1 GS locations ';' separated
    - fq2s (``str``) fastq2 GS locations ';' separated
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "crisp-quant-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "FQ1S":fq1s,
        "FQ2S":fq2s,
        "ALIGN_CUTARGS":env["CRISP_ALIGN_CUTARGS"],
        "ALIGN_RRNA":env["CRISP_ALIGN_RRNA"],
        "QUANT_INDEX":env["CRISP_QUANT_INDEX"],
        "QUANT_ARGS":env["CRISP_QUANT_ARGS"],
        "QUANT_GTF":env["CRISP_QUANT_GTF"]
    }
    check_done(id, "crisp-quant", env, overwrite, local)
    with open(env["ROOT"] / "pipe/crisp-quant/crisp-quant_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/crisp-quant/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/crisp-quant/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/crisp-quant/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
            log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/crisp-quant/crisp-quant_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/crisp-quant/crisp-quant_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
        
@task
def crisp_codac(id, aln, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-standard-8", work_disk="400G", work_disk_type="pd-standard", preemptible=False):
    """Execute CRISP CODAC pipeline

    - id (``str``) unique run identifier
    - aln (``str``) CRISP alignment GS locations ';' separated
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "crisp-codac-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "ALN":aln,
        "ALIGN_FASTA":env["CRISP_ALIGN_FASTA"],
        "CODAC_GTF":env["CRISP_CODAC_GTF"],
        "CODAC_CONFIG":env["CRISP_CODAC_CONFIG"],
        "CODAC_INDEX":env["CRISP_CODAC_INDEX"],
        "CODAC_NEO":env["CRISP_CODAC_NEO"]
    }
    check_done(id, "crisp-codac", env, overwrite, local)
    with open(env["ROOT"] / "pipe/crisp-codac/crisp-codac_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/crisp-codac/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/crisp-codac/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/crisp-codac/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/crisp-codac/crisp-codac_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/crisp-codac/crisp-codac_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def crisp_quasr(id, aln, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-highcpu-16", work_disk="200G", work_disk_type="pd-standard", preemptible=False):
    """Execute CRISP QUASR pipeline

    - id (``str``) unique run identifier
    - aln (``str``) CRISP alignment GS locations ';' separated
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "crisp-quasr-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "ALN":aln,
        "ALIGN_FASTA":env["CRISP_ALIGN_FASTA"],
        "QUASR_COUNTGTF":env["CRISP_QUASR_COUNTGTF"],
        "QUASR_COUNTARGS":env["CRISP_QUASR_COUNTARGS"],
        "QUASR_MIXCRARGS":env["CRISP_QUASR_MIXCRARGS"],
        "QUASR_KEEP":env["CRISP_QUASR_KEEP"],
        "QUASR_COUNT":env["CRISP_QUASR_COUNT"],
        "QUASR_MIXCR":env["CRISP_QUASR_MIXCR"]
    }
    check_done(id, "crisp-quasr", env, overwrite, local)
    with open(env["ROOT"] / "pipe/crisp-quasr/crisp-quasr_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/crisp-quasr/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/crisp-quasr/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/crisp-quasr/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/crisp-quasr/crisp-quasr_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/crisp-quasr/crisp-quasr_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
