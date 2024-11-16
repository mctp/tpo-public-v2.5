#!/usr/bin/env python
from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def cords_unalign(id, aln, name=None, local=False, nowait=False, overwrite=False,
                  machine="n2-highcpu-8", work_disk="auto", work_disk_type="pd-balanced", preemptible=False):
    """Execute CORDS UNALIGN pipeline

    - id (``str``) unique run identifier
    - aln (``str``) Alignment file GS location
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-unalign-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        bam_size_gb = size_gcp(aln, env)
        work_disk = "%sG" % max(200, round(bam_size_gb * 10))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "ALN":aln,
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"]
    }
    check_done(id, "cords-unalign", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-unalign/cords-unalign_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-unalign/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-unalign/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-unalign/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-unalign/cords-unalign_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/cords-unalign/cords-unalign_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_align(id, sample, fq1, fq2, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-highcpu-32", work_disk="auto", work_disk_type="pd-balanced", preemptible=False):
    """Execute CORDS ALIGN pipeline

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
    env["GCP_NAME"] = name if name else "cords-align-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        fq_size_gb = size_gcp(fq1, env) + size_gcp(fq2, env)
        work_disk = "%sG" % max(200, round(fq_size_gb * 15))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "SAMPLE":sample,
        "FQ1":fq1,
        "FQ2":fq2,
        "ALIGN_BWA":env["CORDS_ALIGN_BWA"],
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_INDEX":env["CORDS_ALIGN_INDEX"],
        "ALIGN_DBSNP":env["CORDS_ALIGN_DBSNP"],
        "ALIGN_INDEL":env["CORDS_ALIGN_INDEL"],
        "ALIGN_TRIM":env["CORDS_ALIGN_TRIM"],
        "ALIGN_QUAL":env["CORDS_ALIGN_QUAL"],
        "ALIGN_GENOTYPE":env["CORDS_ALIGN_GENOTYPE"],
    }
    check_done(id, "cords-align", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-align/cords-align_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-align/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-align/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-align/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-align/cords-align_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/cords-align/cords-align_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_umialign(id, sample, fq1s, fq2s, name=None, local=False, nowait=False, overwrite=False,
                   machine="n2-highcpu-32", work_disk="200G", work_disk_type="pd-balanced", preemptible=False):
    """Execute CORDS UMIALIGN pipeline

    - id (``str``) unique run identifier
    - sample (``str``) SAMPLE in read-group
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
    env["GCP_NAME"] = name if name else "cords-umialign-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "SAMPLE":sample,
        "FQ1S":fq1s,
        "FQ2S":fq2s,
        "ALIGN_BWA":env["CORDS_ALIGN_BWA"],
        "ALIGN_UMI":env["CORDS_ALIGN_UMI"],
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_INDEX":env["CORDS_ALIGN_INDEX"],
        "ALIGN_DBSNP":env["CORDS_ALIGN_DBSNP"],
        "ALIGN_INDEL":env["CORDS_ALIGN_INDEL"],
        "ALIGN_TRIM":env["CORDS_ALIGN_TRIM"],
        "ALIGN_GENOTYPE":env["CORDS_ALIGN_GENOTYPE"],
    }
    check_done(id, "cords-umialign", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-umialign/cords-umialign_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-umialign/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-umialign/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-umialign/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-umialign/cords-umialign_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/cords-umialign/cords-umialign_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_misc(id, paln, name=None, local=False, nowait=False, overwrite=False,
               machine="n2-standard-8", work_disk="auto", work_disk_type="pd-ssd", preemptible=False):
    """Execute CORDS Misc pipeline

    - id (``str``) unique run identifier
    - paln (``str``) post-alignment GS location
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-misc-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        bam_size_gb=0
        for bam_loc in paln.split(';'):
            bam_size_gb += size_gcp(os.path.join(bam_loc, ''), env) * 2
        work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "PALN":paln,
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_GENOTYPE":env["CORDS_ALIGN_GENOTYPE"],
        "MISC_TARGETS":env["CORDS_MISC_TARGETS"],
        "MISC_VIRUSFASTA":env["CORDS_MISC_VIRUSFASTA"],
        "MISC_QCMETRICS": env["CORDS_MISC_QCMETRICS"],
        "MISC_GENOTYPE": env["CORDS_MISC_GENOTYPE"],
        "MISC_COVMETRICS":env["CORDS_MISC_COVMETRICS"],
        "MISC_COVPERBASE":env["CORDS_MISC_COVPERBASE"],
        "MISC_VIRMER":env["CORDS_MISC_VIRMER"],
        "MISC_VIRMERARGS":env["CORDS_MISC_VIRMERARGS"],
        "MISC_MERGE":env["CORDS_MISC_MERGE"],
    }
    check_done(id, "cords-misc", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-misc/cords-misc_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-misc/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-misc/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-misc/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-misc/cords-misc_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/cords-misc/cords-misc_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
        
@task
def cords_postalign(id, aln, name=None, local=False, nowait=False, overwrite=False,
                    machine="n2-highcpu-32", work_disk="auto", work_disk_type="pd-ssd", preemptible=False):
    """Execute CORDS POSTALIGN pipeline

    - id (``str``) unique run identifier
    - aln (``str``) alignment GS locations
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-postalign-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
      bam_size_gb=0
      for bam_loc in aln.split(';'):
        bam_size_gb += size_gcp(os.path.join(bam_loc, ''), env)*3
      work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "ALN":aln,
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_DBSNP":env["CORDS_ALIGN_DBSNP"],
        "ALIGN_INDEL":env["CORDS_ALIGN_INDEL"],
        "ALIGN_FIXRG":env["CORDS_ALIGN_FIXRG"],
    }
    check_done(id, "cords-postalign", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-postalign/cords-postalign_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-postalign/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-postalign/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-postalign/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-postalign/cords-postalign_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/cords-postalign/cords-postalign_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_somatic(id, alnt=None, alnn=None, regions=None, name=None, local=False, nowait=False, overwrite=False,
                  machine="n2-highcpu-32", work_disk="auto", work_disk_type="pd-ssd", preemptible=False):
    """Execute CORDS SOMATIC pipeline

    - id (``str``) unique run identifier
    - alnt (``str``) post-alignment tumor GS location
    - alnn (``str``) post-alignment normal GS location
    - regions (``str``) only process select regions of input [bed-file or chr-start:end]
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-somatic-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
      bam_size_gb=0
      if alnt is not None:
        bam_size_gb += size_gcp(os.path.join(alnt, ''), env) * 2.5
      if alnn is not None:
        bam_size_gb += size_gcp(os.path.join(alnn, ''), env) * 2.5
      work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    if regions:
        regions_gcp(regions, "cords-somatic", env)
    vals = {
        "ID":id,
        "ALNT":alnt.rstrip("/") if alnt is not None else None,
        "ALNN":alnn.rstrip("/") if alnn is not None else None,
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_INDEL":env["CORDS_ALIGN_INDEL"],
        "ALIGN_INDEX":env["CORDS_ALIGN_INDEX"],
        "SOMATIC_TNSCOPEARGS":env["CORDS_SOMATIC_TNSCOPEARGS"],
        "SOMATIC_TINSCOPEARGS":env["CORDS_SOMATIC_TINSCOPEARGS"],
        "SOMATIC_DNASCOPEARGS":env["CORDS_SOMATIC_DNASCOPEARGS"],
        "SOMATIC_DBSNP":env["CORDS_SOMATIC_DBSNP"],
        "SOMATIC_REALIGN":env["CORDS_SOMATIC_REALIGN"],
        "SOMATIC_TNSCOPE":env["CORDS_SOMATIC_TNSCOPE"],
        "SOMATIC_TINSCOPE":env["CORDS_SOMATIC_TINSCOPE"],
        "SOMATIC_DNASCOPE":env["CORDS_SOMATIC_DNASCOPE"],
        "SOMATIC_GERMLINEDB":env["CORDS_SOMATIC_GERMLINEDB"],
    }
    check_done(id, "cords-somatic", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-somatic/cords-somatic_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-somatic/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-somatic/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-somatic/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        if alnn and alnt:
            ret = run_app(env["ROOT"] / "pipe/cords-somatic/cords-somatic_tn_local.sh", [], env=env)
            log_run_outerr(ret)
        else:
            ret = run_app(env["ROOT"] / "pipe/cords-somatic/cords-somatic_to_local.sh", [], env=env)
            log_run_outerr(ret)
    else:
        if alnn and alnt:
            ret = run_app(env["ROOT"] / 'pipe/cords-somatic/cords-somatic_tn_gcp.sh', env=env)
            log_run_outerr(ret)
        else:
            ret = run_app(env["ROOT"] / 'pipe/cords-somatic/cords-somatic_to_gcp.sh', env=env)
            log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_germline(id, aln=None, sex=None, regions=None, name=None, local=False, nowait=False, overwrite=False,
                   machine="n2-highcpu-32", work_disk="250G", work_disk_type="pd-ssd", preemptible=False):
    """Execute CORDS germline pipeline

    - id (``str``) unique run identifier
    - aln (``str``) post-alignment GS location
    - sex (``str``) chromosomal sex (XX, XY) of the sample (blank for auto-detection)
    - regions (``str``) only process select regions of input [bed-file or chr-start:end]
    - name (``str``) name of GCP instance [default from `id`] 
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-germline-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    if regions:
        regions_gcp(regions, "cords-germline", env)
    vals = {
        "ID":id,
        "ALN":aln.rstrip("/"),
        "SEX":sex.upper() if sex else "",
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_INDEL":env["CORDS_ALIGN_INDEL"],
        "GERMLINE_DNASCOPEARGS":env["CORDS_GERMLINE_DNASCOPEARGS"],
        "GERMLINE_DBSNP":env["CORDS_GERMLINE_DBSNP"],
        "GERMLINE_SEXSNP":env["CORDS_GERMLINE_SEXSNP"],
        "GERMLINE_INTERVALS":env["CORDS_GERMLINE_INTERVALS"],
        "GERMLINE_REALIGN":env["CORDS_GERMLINE_REALIGN"],
        "GERMLINE_DNASCOPE":env["CORDS_GERMLINE_DNASCOPE"],
    }
    check_done(id, "cords-germline", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-germline/cords-germline_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-germline/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-germline/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-germline/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-germline/cords-germline_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/cords-germline/cords-germline_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)

@task
def cords_structural(id, alnt=None, alnn=None, regions=None, name=None, local=False, nowait=False, overwrite=False,
                     machine="n2-custom-8-32768", work_disk="auto", work_disk_type="pd-ssd", preemptible=False):
    """Execute CORDS Structural pipeline

    - id (``str``) unique run identifier
    - alnt (``str``) post-alignment tumor GS location
    - alnn (``str``) post-alignment normal GS location
    - regions (``str``) only process select regions of input [bed-file or chr-start:end]
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-structural-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        bam_size_gb=0
        if alnt is not None:
            bam_size_gb += size_gcp(os.path.join(alnt, ''), env) * 2.5
        if alnn is not None:
            bam_size_gb += size_gcp(os.path.join(alnn, ''), env) * 2.5
        work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    if regions:
        regions_gcp(regions, "cords-structural", env)
    vals = {
        "ID":id,
        "ALNT":alnt.rstrip("/"),
        "ALNN":alnn.rstrip("/"),
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "ALIGN_INDEX":env["CORDS_ALIGN_INDEX"],
        "STRUCTURAL_DEBUG":env["CORDS_STRUCTURAL_DEBUG"],
        "STRUCTURAL_TNSCOPE":env["CORDS_STRUCTURAL_TNSCOPE"],
        "STRUCTURAL_TNSCOPEARGS":env["CORDS_STRUCTURAL_TNSCOPEARGS"],
        "STRUCTURAL_MANTA":env["CORDS_STRUCTURAL_MANTA"],
        "STRUCTURAL_MANTAARGS":env["CORDS_STRUCTURAL_MANTAARGS"],
        "STRUCTURAL_SVABA":env["CORDS_STRUCTURAL_SVABA"],
        "STRUCTURAL_SVABAARGS":env["CORDS_STRUCTURAL_SVABAARGS"],
        "STRUCTURAL_GRIDSS":env["CORDS_STRUCTURAL_GRIDSS"],
        "STRUCTURAL_GRIDSSARGS":env["CORDS_STRUCTURAL_GRIDSSARGS"],
        "STRUCTURAL_GRIPSSARGS":env["CORDS_STRUCTURAL_GRIPSSARGS"]
    }
    check_done(id, "cords-structural", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-structural/cords-structural_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-structural/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-structural/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cords-structural/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-structural/cords-structural_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/cords-structural/cords-structural_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
        
@task
def cords_cnvex(id, alnt=None, alnn=None, tvar=None, gvar=None, cnvx=None, sex=None, picks=None, pickg=None,
                plots=None, plotg=None, name=None, local=False, nowait=False, overwrite=False,
                machine="n2-standard-16", work_disk="auto", work_disk_type="pd-balanced", preemptible=False):
    """Execute CORDS CNVEX pipeline

    Provide a combination of `alnt` (required) and optinally `alnn`, `tvar`, `gvar` 
    Alternatively provied `cnvx` 

    - id (``str``) unique run identifier
    - alnt (``str``) CORDS tumor post-alignment GS location
    - alnn (``str``) CORDS normal post-alignment GS location
    - tvar (``str``) CORDS exome variant GS location
    - gvar (``str``) CORDS genome variant GS location
    - cnvx (``str``) CORDS CNVEX output folder (optional)
    - sex (``str``) [male or female] (optional)
    - picks (``str``) Pick model(s) for somatic disgest(s) [topN or custom 'purity:ploidy']
    - pickg (``str``) Pick model(s) for germline disgest(s) [topN or custom 'purity:ploidy']
    - plots (``str``) Plot model(s) for somatic disgest(s) ['arranged-plot']
    - plotg (``str``) Plot model(s) for germline disgest(s) ['arranged-plot']
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cords-cnvex-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
      bam_size_gb=0
      if alnt is not None:
        bam_size_gb += size_gcp(os.path.join(alnt, ''), env) * 2.5
      if alnn is not None:
        bam_size_gb += size_gcp(os.path.join(alnn, ''), env) * 2.5
      work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "ALNT":alnt.rstrip("/") if alnt else "",
        "ALNN":alnn.rstrip("/") if alnn else "",
        "TVAR":tvar if tvar else "",
        "GVAR":gvar if gvar else "",
        "CNVX":cnvx if cnvx else "",
        "CNVSEX":sex if sex else "",
        "ALIGN_FASTA":env["CORDS_ALIGN_FASTA"],
        "CNVEX_ASSEMBLY":env["CORDS_CNVEX_ASSEMBLY"],
        "CNVEX_SETTINGS":env["CORDS_CNVEX_SETTINGS"],
        "CNVEX_CAPTURE":env["CORDS_CNVEX_CAPTURE"],
        "CNVEX_POOL":env["CORDS_CNVEX_POOL"],
        "CNVEX_POPAF":env["CORDS_CNVEX_POPAF"],
        "CNVEX_PROCESS":env["CORDS_CNVEX_PROCESS"],
        "CNVEX_SOMATIC":env["CORDS_CNVEX_SOMATIC"],
        "CNVEX_SOMATICSEARCH":env["CORDS_CNVEX_SOMATICSEARCH"],
        "CNVEX_SOMATICPICK":picks if picks else env["CORDS_CNVEX_SOMATICPICK"], ## override config
        "CNVEX_SOMATICPLOT":plots if plots else env["CORDS_CNVEX_SOMATICPLOT"], ## override config
        "CNVEX_GERMLINE":env["CORDS_CNVEX_GERMLINE"],
        "CNVEX_GERMLINESEARCH":env["CORDS_CNVEX_GERMLINESEARCH"],
        "CNVEX_GERMLINEPICK":pickg if pickg else env["CORDS_CNVEX_GERMLINEPICK"], ## override config
        "CNVEX_GERMLINEPLOT":plotg if plotg else env["CORDS_CNVEX_GERMLINEPLOT"], ## override config
    }
    check_done(id, "cords-cnvex", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cords-cnvex/cords-cnvex_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cords-cnvex/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cords-cnvex/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = "gsutil -m cp %s gs://%s/pipe/cords-cnvex/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cords-cnvex/cords-cnvex_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / 'pipe/cords-cnvex/cords-cnvex_gcp.sh', env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
