#!/usr/bin/env python3
import os
import sys
import logging
from string import Template
from moke import *
from mods import *
from mods.bcl import *
from mods.cords import *
from mods.crisp import *
from mods.carat import *
from mods.cargo import *
from mods.job import *
from mods.unit_test import *
import mods

## REFS

@task
def refs_pull():
    """Pull references from GCP
    """
    out = env["RUNTIME_REFS"]
    if not os.path.exists(out):
        os.makedirs(out)
    cmd = GSUTIL + "rsync -r gs://%s %s"
    cmd = cmd % (env["TPO_REFS_VER"], out)
    ret = run_app(cmd, env=env)
    log_run_outerr(ret)

@task
def refs_push():
    """Push references to GCP
    """
    inp = env["RUNTIME_REFS"]
    if not os.path.exists(inp):
        log("directory '%s' is missing" % inp, logging.ERROR)
        sys.exit(1)
    cmd = GSUTIL + "rsync -d -r %s gs://%s"
    cmd = cmd % (inp, env["TPO_REFS_VER"])
    ret = run_app(cmd, env=env)
    log_run_outerr(ret)

## CODE

@task
def code_build(dirty=False, nocache=False, remove=False):
    """Build TPO Code Images

    - dirty (``bool``) build code including uncommited files
    - nocache (``bool`) docker build ignore cache
    - remove (``bool``) remove temporary built R libraries
    """
    env["BUILD_COMMIT"] = "dirty" if dirty else "commit"
    env["BUILD_DOCKER_NOCACHE"] = "--no-cache" if nocache else ""
    env["BUILD_REMOVE"] = "remove" if remove else "keep"
    
    cmd = "%s/base/images/build_tpocode.sh"
    cmd = cmd % (env["ROOT"],)
    ret = run_app(cmd, env=env)
    log_run_outerr(ret)

@task
def code_push():
    """Push code image to GCR
    """
    cmd = "%s/base/gcp/gcp-code-push.sh"
    cmd = cmd % (env["ROOT"],)
    ret = run_app(cmd, [], env=env)
    log_run_outerr(ret)

@task
def code_pull():
    """Pull code image from GCR
    """
    cmd = "%s/base/gcp/gcp-code-pull.sh"
    cmd = cmd % (env["ROOT"],)
    ret = run_app(cmd, [], env=env)
    log_run_outerr(ret)
    
## BOOT

@task
def boot_build(image='all', nocache=False):
    """Build TPO Boot Images

    - image (``str``) which image to build, other options: base, refs, cords, crisp, carat, dev
    - nocache (``bool`) docker build ignore cache
    """
    env["BUILD_DOCKER_NOCACHE"] = "--no-cache" if nocache else ""

    if image=="base" or image == 'all': 
        cmd = "%s/base/images/build_tpobase.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

    if image=="refs" or image == 'all':
        cmd = "%s/base/images/build_tporefs.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

    if image=="cords" or image == 'all':
        cmd = "%s/base/images/build_tpocords.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

    if image=="crisp" or image == 'all':
        cmd = "%s/base/images/build_tpocrisp.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

    if image=="carat" or image == 'all':
        cmd = "%s/base/images/build_tpocarat.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

    if image=="dev":
        cmd = "%s/base/images/build_tpodev.sh"
        cmd = cmd % (env["ROOT"],)
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)

@task
def boot_push():
    """Push TPO Boot images to GCR
    """
    cmd = "%s/base/gcp/gcp-boot-push.sh"
    cmd = cmd % (env["ROOT"],)
    ret = run_app(cmd, [], env=env)
    log_run_outerr(ret)

@task
def boot_pull():
    """Pull TPO Boot images from GCR
    """
    cmd = "%s/base/gcp/gcp-boot-pull.sh"
    cmd = cmd % (env["ROOT"],)
    ret = run_app(cmd, [], env=env)
    log_run_outerr(ret)

## GCP

@task
def gcp_configure(nfs_share=False, nfs_name="tpo-nfs", nfs_tier="STANDARD", nfs_size="1T"):
    """Configure needed GCP resources

    TODO: open-ssh network tags
    """
    if nfs_share:
        env["GCP_NFS_SHARE"] = nfs_name
        env["GCP_NFS_TIER"] = nfs_tier
        env["GCP_NFS_SIZE"] = nfs_size
        ret = run_app(env["ROOT"] / "base/gcp/gcp-nfs-instance.sh", [], env=env)
        log_run_outerr(ret)

@task
def gcp_boot_build():
    """Build boot GCP image/disk
    """
    env["GCP_NAME"] = "boot-" + env["TPO_BOOT_VER"].replace(".", "")
    ret = run_app(env["ROOT"] / "base/gcp/gcp-boot-build.sh", [], env=env)
    log_run_outerr(ret)

@task
def gcp_root_build():
    """Build root GCP image/disk
    """
    env["GCP_NAME"] = "root-" + env["TPO_ROOT_VER"].replace(".", "")
    cmd = GSUTIL + "cp %s gs://%s/root"
    cmd = cmd % (env["ROOT"] / "base/root/*", env["RUNTIME_WORK_BUCKET"])
    ret = run_app(cmd, env=env)
    log_run_outerr(ret)
    ret = run_app(env["ROOT"] / "base/gcp/gcp-root-build.sh", [], env=env)
    log_run_outerr(ret)

@task
def gcp_root_push(dest=None, work_disk="250G", work_disk_type="pd-balanced"):
    """Push root/tpo from image to GS
    """
    env["GCP_NAME"] = "push-" + env["TPO_ROOT_VER"].replace(".", "")
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["RUNTIME_WORK_BUCKET"] = dest if dest else env["RUNTIME_WORK_BUCKET"]
    ret = run_app(env["ROOT"] / "base/gcp/gcp-root-push.sh", [], env=env)
    log_run_outerr(ret)

@task
def gcp_root_pull(dest=None):
    """Pull root/tpo from GS to local
    """
    env["TEMP_BUCKET"] = dest if dest else env["RUNTIME_WORK_BUCKET"]
    ret = run_app(env["ROOT"] / "base/gcp/gcp-root-pull.sh", [], env=env)
    log_run_outerr(ret)

@task
def gcp_work_build():
    """Build work GCP image/disk
    """
    env["GCP_NAME"] = "work-" + env["GCP_WORK_VOLUME"].replace(".", "")
    ret = run_app(env["ROOT"] / "base/gcp/gcp-work-build.sh", [], env=env)
    log_run_outerr(ret)

@task
def gcp_instance(name, job=None, tags="open-ssh", nfs_share=None,
                 machine="n2-standard-8", boot_disk="100G", work_disk="200G", work_disk_type="local-ssd",
                 preemptible=False):
    """Start GCP instance to develop pipelines

    - name (``str``) instance name
    - job (``str``) custom script (gcp_script.sh,docker_script.sh)
    - tags (``str``) additional network tags [open-ssh]
    - nfs_share (``str``) name of nfs share
    - machine (``str``) machine type
    - boot_disk (``str``) TPO boot machine disk size
    - work_disk (``str``) TPO work machine disk size
    - work_disk_type (``str``) work disk size
    - preemptible (``bool``) start pre-emptible
    """
    env["GCP_NAME"] = name
    env["GCP_JOB"] = job if job else ""
    if tags:
        env["GCP_TAG"] = tags
    if not nfs_share is None:
        env["GCP_NFS_SHARE"] =  nfs_share
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_BOOT_DISK"] = boot_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"] + (" --preemptible" if preemptible else "")
    ret = run_app(env["ROOT"] / "base/gcp/gcp-instance.sh", [], env=env)
    log_run_outerr(ret)

@task
def curve_db_init(name, tags="open-ssh",
                 machine="n2-standard-8", work_disk="500G", work_disk_type="pd-ssd"):
    """Start GCP instance and initialize a curve DB

    - name (``str``) instance name
    - tags (``str``) additional network tags [open-ssh]
    - machine (``str``) machine type
    - boot_disk (``str``) TPO boot machine disk size
    - root_disk (``str``) TPO root machine disk size
    - work_disk (``str``) TPO work machine disk size
    - work_disk_type (``str``) work disk size
    """
    env["GCP_NAME"] = name
    if tags:
        env["GCP_TAG"] = tags
    env["GCP_MACHINE"] = machine
    env["GCP_WORK_DISK"] = work_disk
    env["GCP_WORK_DISK_TYPE"] = work_disk_type
    env["GCP_MACHINE_ARGS"] = env["GCP_MACHINE_ARGS"]
    ret = run_app(env["ROOT"] / "base/gcp/curve-db-init.sh", [], env=env)
    log_run_outerr(ret)

@task
def curve_db_add(db_name,tpo_path="./"):
    """Add a new database to an existing GCP Curve instance

    - db_name (``str``) name of database to add
    - tpo_path (``str``) path to TPO [./]
    """
    env["CURVE_DB"]=os.environ["CURVE_DB"]
    ret = run_app(env["ROOT"] / "base/gcp/curve-db-add.sh" ,[db_name, tpo_path], env=env)
    log_run_outerr(ret)

@task
def curve_init_db(instance=None, local=False):
    """Initialize a curve DB on local or remote instance

    - instance (``str``) instance name
    - local (``bool``) execute locally
    """
    if local:
        ret = run_app(env["ROOT"] / "rlibs/curve/exec/curve-init-db.sh", [], env=env)
        log_run_outerr(ret)
    elif instance:
        env["CURVE_INSTANCE"] = instance
        ret = run_app(env["ROOT"] / "base/gcp/gcp-curve-init-db.sh", [], env=env)
        log_run_outerr(ret)
    else:
        log("Instance name or local missing.", ERROR)
        sys.exit(1)       

@task
def curve_add_db(dbname=None, instance=None, local=False):
    """Add a new database to an existing CURVE SQL server

    Parses the CURVE_DB environmental variable

    The `dbname` argument takes precedence over CURVE_DB and configuration file.
    GCP instance name require if not local.

    - dbname (``str``) new database name
    - instance (``str``) new gcp instance
    - local (``bool``) execute locally
    """
    ## parse CURVE_DB env variable
    CURVE_DB = os.environ.get("CURVE_DB", None)
    if CURVE_DB:
        host, port, name, user, pswd = CURVE_DB.split(":")
        env["CURVE_DB_HOST"] = host
        env["CURVE_DB_PORT"] = port
        env["CURVE_DB_NAME"] = dbname if dbname else name
        env["CURVE_DB_USER"] = user
        env["CURVE_DB_PASS"] = pswd
    if local:
        ret = run_app(env["ROOT"] / "rlibs/curve/exec/curve-add-db.sh", [], env=env)
        log_run_outerr(ret)
    elif instance:
        env["CURVE_INSTANCE"] = instance
        ret = run_app(env["ROOT"] / "base/gcp/gcp-curve-add-db.sh", [], env=env)
        log_run_outerr(ret)
    else:
        log("Instance name or local missing.", ERROR)
        sys.exit(1)

if __name__ == "__main__":
    ROOT = path(os.path.dirname(os.path.abspath(__file__)))
    env = make_env(ROOT)
    mods.bcl.__dict__["env"] = env
    mods.cords.__dict__["env"] = env
    mods.crisp.__dict__["env"] = env
    mods.carat.__dict__["env"] = env
    mods.cargo.__dict__["env"] = env
    mods.job.__dict__["env"] = env
    mods.unit_test.__dict__["env"] = env
    task()
