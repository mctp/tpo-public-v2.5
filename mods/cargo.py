from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def cargo_vault(id, case=None, cohort=None, anno=None, misc_t=None, misc_n=None, cnvex=None, tquasr=None, nquasr=None, codac=None, name=None, local=False, nowait=False,
                overwrite=False, machine="n2-custom-8-32768", work_disk="200G", work_disk_type="boot-disk", preemptible=False):
    """Execute CARGO VAULT pipeline

    - id (``str``) unique run identifier
    - case (``str``) case id
    - cohort (``str``) cohort id
    - anno (``str``) CARAT ANNO output GS location
    - misc_t (``str``) CORDS MISC tumor output GS location
    - misc_n (``str``) CORDS MISC normal output GS location
    - cnvex (``str``) CORDS CNVEX output GS location
    - tquasr (``str``) CRISP QUASR tumor output GS location
    - nquasr (``str``) CRISP QUASR normal output GS location
    - codac (``str``) CRISP CODAC output GS location
    - name (``str``) instance name [default from `id`]
    - local (``bool``) execute locally
    - nowait (``bool``) do not wait for instance to shut down
    - overwrite (``bool``) run even if output is present
    - machine (``str``) machine type
    - work_disk (``str``) work disk size
    - work_disk_type (``str``) work disk type
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name if name else "cargo-vault-" + id.lower().replace("_", "-").replace(".", "-")
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    vals = {
        "ID":id,
        "CASE":case if case else "",
        "COHORT":cohort if cohort else "",
        "CARAT_ANNO":anno if anno else "",
        "CORDS_CNVEX":cnvex if cnvex else "",
        "CORDS_MISC_TUMOR":misc_t if misc_t else "",
        "CORDS_MISC_NORMAL":misc_n if misc_n else "",
        "CRISP_QUASR_TUMOR":tquasr if tquasr else "",
        "CRISP_QUASR_NORMAL":nquasr if nquasr else "",
        "CRISP_CODAC":codac if codac else "",
        "VAULT_GENEGTF":env["CARGO_VAULT_GENEGTF"],
        "VAULT_HOMOPOLYMERS":env["CARGO_VAULT_HOMOPOLYMERS"],
        "VAULT_TARGETS":env["CARGO_VAULT_TARGETS"],
        "VAULT_SETTINGS":env["CARGO_VAULT_SETTINGS"],
        "VAULT_DUMP":env["CARGO_VAULT_DUMP"]
    }
    check_done(id, "cargo-vault", env, overwrite, local)
    with open(env["ROOT"] / "pipe/cargo-vault/cargo-vault_config.template") as fh:
        src = Template(fh.read())
        src = src.substitute(vals)
        with tmp_file() as ofh:
            ofh.write(src)
            ofh.flush()
            if local:
                cmd = "mkdir -p %s/cargo-vault/%s" % (env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
                cmd = "cp %s %s/cargo-vault/%s/config.txt" % (ofh.name, env["RUNTIME_RUNS"],  env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
            else:
                cmd = GSUTIL + "cp %s gs://%s/pipe/cargo-vault/%s/config.txt"
                cmd = cmd % (ofh.name, env["RUNTIME_WORK_BUCKET"], env["GCP_NAME"])
                ret = run_app(cmd, env=env)
                log_run_outerr(ret)
    if local:
        ret = run_app(env["ROOT"] / "pipe/cargo-vault/cargo-vault_local.sh", [], env=env)
        log_run_outerr(ret)
    else:
        ret = run_app(env["ROOT"] / "pipe/cargo-vault/cargo-vault_gcp.sh", [], env=env)
        log_run_outerr(ret)
    if nowait == False and not local:
        wait_gcp(env)
