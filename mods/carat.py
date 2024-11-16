from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def carat_anno(id, somatic=None, structural=None, name=None, local=False, nowait=False, overwrite=False,
               ncores=None, machine="n2-standard-8", work_disk="auto", work_disk_type="boot-disk", preemptible=False):
    """Execute CARAT ANNO pipeline

    - id (``str``) unique run identifier
    - somatic (``str``) cords-somatic GS location(s)
    - structural (``str``) cords-structural GS location(s)
    - name (``str``) name of GCP instance [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - ncores (``str``) limit the number of cores used
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "carat-anno-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    if work_disk == "auto":
        bam_size_gb=0
        for bam_loc in somatic.split(';'):
            bam_size_gb += size_gcp(os.path.join(bam_loc, ''), env) * 1.5
        work_disk = "%sG" % max(200, round(bam_size_gb))
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    env["GCP_MACHINE_NCORES"] = ncores if ncores else ""
    vals = {
        "ID":id,
        "SOMATIC":somatic if somatic else "",
        "STRUCTURAL":structural if structural else "",
        "ANNO_ASSEMBLY":env["CARAT_ANNO_ASSEMBLY"],
        "ANNO_REGIONS":env["CARAT_ANNO_REGIONS"],
        "ANNO_VEPREFS":env["CARAT_ANNO_VEPREFS"],
        "ANNO_VEPARGS":env["CARAT_ANNO_VEPARGS"],
        "ANNO_CONFIG":env["CARAT_ANNO_CONFIG"],
        "ANNO_REPORTARGS":env["CARAT_ANNO_REPORTARGS"],
        "ANNO_FILTERSOMATIC":env["CARAT_ANNO_FILTERSOMATIC"],
        "ANNO_FILTERGERMLINE":env["CARAT_ANNO_FILTERGERMLINE"],
    }
    check_done(id, "carat-anno", env, overwrite, local)
    with open(env["ROOT"] / "pipe/carat-anno/carat-anno_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/carat-anno/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/carat-anno/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/carat-anno/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/carat-anno/carat-anno_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/carat-anno/carat-anno_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
