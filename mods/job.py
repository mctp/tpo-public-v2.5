from moke import *
from . import *
from string import Template

@task
def job(id, script, params, name=None, local=False, nowait=False, overwrite=False,
             machine="n2-highcpu-8", work_disk="200G", work_disk_type="pd-balanced", preemptible=False):
    """Execute JOB pipeline

    - id (``str``) unique run identifier
    - script (``str``) script to be executed (local file)
    - params (``str``) configuration file for shell script (local file)
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "job-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id
    }
    job_vals = {"job_params":params, "job_script":script}
    check_done(id, "job", env, overwrite, local)
    with open(env["ROOT"] / "pipe/job/job_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/job/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/job/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                for k, v in job_vals.items():
                    cmd = "cp %s %s/job/%s/%s" % (v, env["RUNTIME_RUNS"],  env["GCP_NAME"], k)
                    ret = run_app(cmd, env=env)
                    log_run_outerr(ret)
            else:
                cmd = GSUTIL + " cp %s gs://%s/pipe/job/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                for k, v in job_vals.items():
                    cmd = GSUTIL + " cp %s gs://%s/pipe/job/%s/%s"
                    cmd = cmd % (v, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"], k)
                    ret = run_app(cmd, env=env)
                    log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/job/job_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/job/job_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
