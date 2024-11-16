from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def bcl(id, tar, lib, name=None, reverse=False, local=False, nowait=False, overwrite=False,
        machine="n2-highcpu-32", work_disk="500G", work_disk_type="pd-balanced", preemptible=False):
    """Execute BCL pipeline

    - id (``str``) unique run identifier
    - tar (``str``) Flowcell tar file GS location
    - lib (``str``) Flowcel library sheet GS location
    - name (``str``) instance name (overrides default)
    - reverse (``bool``) reverse complement index2 (1.5 chemistry)
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "bcl-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "TAR":tar,
        "LIB":lib,
        "REVERSE":reverse,
        "BCL2FASTQ_ARGS":env["BCL_BCL2FASTQ_ARGS"],
        "BCL2FASTQ_SHEET":env["BCL_BCL2FASTQ_SHEET"],
        "OUTPUT_FORMAT":env["BCL_OUTPUT_FORMAT"]
    }
    check_done(id, "bcl", env, overwrite, local)
    with open(env["ROOT"] / "pipe/bcl/bcl_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/bcl/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/bcl/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/bcl/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/bcl/bcl_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/bcl/bcl_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
