import configparser
import os
import sys
import time
import multiprocessing
import logging
from collections import defaultdict
from random import randrange
from string import Template
from moke import * #@UnusedWildImport
env = defaultdict(lambda: "") # suppresses warnings

GSUTIL = "gsutil -q -m -o GSUtil:parallel_composite_upload_threshold=100M " + \
                      "-o GSUtil:parallel_process_count=8 " + \
                      "-o GSUtil:parallel_thread_count=2 "

def log_run_outerr(code_stdout_stderr_cmd, pfx=""):
    code, stdout, stderr, cmd = code_stdout_stderr_cmd
    if code:
        level = logging.ERROR
    else:
        level = logging.DEFAULT
    cmd = "%sshell[%s]: %s" % (pfx, code, cmd)
    log(cmd, level=level)
    outs = stdout.decode("utf-8").splitlines()
    for out in outs:
        out = "%sstdout: %s" % (pfx, out)
        log(out, level=level)
    errs = stderr.decode("utf-8").splitlines()
    for err in errs:
        err = "%sstderr: %s" % (pfx, err)
        log(err, level=level)
    return code

def make_env(root):
    ## this is a hack!
    config_file = dict(task._funcparse(None).parse_args().__dict__).get("config").name
    if config_file.endswith("dist-packages/moke/data/mokefile.ini"):
        sys.stderr.write("Configuration file not provided, failing.\n")
        sys.exit(1)
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform=str
    cp.read(config_file)
    env = defaultdict(lambda: "")
    ## This is needed for gsutil to find a suitable python when call through Popen
    env["PATH"] = os.environ.get("PATH")
    env["ROOT"] = root
    ## Sections without hard-coded defaults
    sections = cp.sections()
    sections.remove("RUNTIME")
    for section in sections:
        env.update([(section + "_" + k, v) for k,v in cp.items(section)])
    ## GCP
    zones = env["GCP_ZONE"].split(",")
    env["GCP_ZONE"] = zones.pop(randrange(len(zones)))
    ## CODE
    if cp.get("TPO", "CODE_VER"):
        env["TPO_CODE_VER"] = cp.get("TPO", "CODE_VER")
    else:
        commit = run_app("git --git-dir %s/.git rev-parse --short HEAD" % env["ROOT"])[1].strip()
        env["TPO_CODE_VER"] = commit
    ## SECRETS
    SENTIEON_LICENSE = os.environ.get("SENTIEON_LICENSE")
    if SENTIEON_LICENSE:
        env["SECRETS_SENTIEON_LICENSE"] = SENTIEON_LICENSE
    AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID")
    env["SECRETS_AWS_ACCESS_KEY_ID"] = AWS_ACCESS_KEY_ID if AWS_ACCESS_KEY_ID else env.get("SECRETS_AWS_ACCESS_KEY_ID", "NOAWSKEY")
    AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY")
    env["SECRETS_AWS_SECRET_ACCESS_KEY"] = AWS_SECRET_ACCESS_KEY if AWS_SECRET_ACCESS_KEY else env.get("SECRETS_AWS_SECRET_ACCESS_KEY", "NOAWSKEY")
    ## RUNTIME
    if not cp.get("RUNTIME", "WORK_BUCKET") or not cp.get("RUNTIME", "ROOT") or not cp.get("RUNTIME", "WORK"):
        raise Exception("Required RUNTIME settings not set.")
    env["RUNTIME_WORK_BUCKET"] = path(cp.get("RUNTIME", "WORK_BUCKET"))
    env["RUNTIME_ROOT"] = path(cp.get("RUNTIME", "ROOT"))
    env["RUNTIME_WORK"] = path(cp.get("RUNTIME", "WORK"))
    DEFAULT_TEMP = env["RUNTIME_WORK"] / "tmp"
    env["RUNTIME_TEMP"] = path(cp.get("RUNTIME", "TEMP") if cp.get("RUNTIME", "TEMP") else DEFAULT_TEMP)
    DEFAULT_RUNS = env["RUNTIME_WORK"] / "runs"
    env["RUNTIME_RUNS"] = path(cp.get("RUNTIME", "RUNS") if cp.get("RUNTIME", "RUNS") else DEFAULT_RUNS)
    DEFAULT_REFS = env["RUNTIME_WORK"] / "refs"
    env["RUNTIME_REFS"] = path(cp.get("RUNTIME", "REFS") if cp.get("RUNTIME", "REFS") else DEFAULT_REFS)
    DEFAULT_WORK_DONE = "skip"
    env["RUNTIME_WORK_DONE"] = path(cp.get("RUNTIME", "WORK_DONE") if cp.get("RUNTIME", "WORK_DONE") else DEFAULT_WORK_DONE)
    ## Values with overrides from -cargs
    config_args = dict(task._funcparse(None).parse_args().__dict__).get("cargs", '')
    if config_args:
        overrides = [kv.split("::") for kv in config_args.split(";")]
        env.update(overrides)
    return(env)

def wait_gcp(env, delay=60):
    cmd = 'gcloud compute instances list --project=%s --filter="NAME=%s" --format="value(NAME)"'
    cmd = cmd % (env["GCP_PROJECT"], env["GCP_NAME"])
    while True:
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)
        code, stdout, stderr, cmd = ret
        if not code and not stdout:
            break
        time.sleep(delay)

def done_gcp(id, pipe, env):
    cmd = GSUTIL + "ls -d gs://%s/repo/%s/%s"
    cmd = cmd % (env["RUNTIME_WORK_BUCKET"], pipe, id)
    ret = run_app(cmd, env=env)
    exists = ret[0] == 0
    return exists

def size_gcp(gs, env, unit=1024*1024*1024):
    cmd = GSUTIL + "du -s %s"
    cmd = cmd % (gs,)
    ret = run_app(cmd, env=env)
    if ret[0] != 0:
        log("Input file missing.")
        sys.exit(1)
    size_bytes = int(ret[1].split(b" ")[0])
    size_unit = size_bytes / unit
    return size_unit

def regions_gcp(regions, pipe, env):
    ## region or file
    with tmp_file() as ofh:
        if not os.path.exists(regions):
            ## write region as bed
            ofh.write(regions.replace(":", "\t").replace("-", "\t") + "\n")
            ofh.flush()
            bed = ofh.name
        else:
            bed = regions
        cmd = GSUTIL + "cp %s gs://%s/pipe/%s/%s/regions.bed"
        cmd = cmd % (bed, env["RUNTIME_WORK_BUCKET"], pipe, env["GCP_NAME"])
        ret = run_app(cmd, env=env)
        log_run(ret)
    return None

def done_local(id, pipe, env):
    exists = False
    return exists
        
def check_done(id, pipe, env, overwrite, local):
    if local:
        done = done_local(id, pipe, env)
    else:
        done = done_gcp(id, pipe, env)
    skip = env["RUNTIME_WORK_DONE"] == ""
    if not done:
        log("Output does not exists.")
    elif done and overwrite:
        log("Overwrite enabled, contiune.")
    elif done and skip:
        log("Overwrite disabled, skipping.")
    else:
        log("Overwrite disabled, skip disabled, failing.")
        sys.exit(1)
    
def copy_file(f, pipe, env, local):
    if local:
        cmd = "mkdir -p %s/%s/%s" % (env["RUNTIME_RUNS"], pipe, env["GCP_NAME"])
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)
        cmd = "cp %s %s/%s/%s/%s" % (f, env["RUNTIME_RUNS"], pipe, env["GCP_NAME"], f.basename())
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)
    else:
        cmd = GSUTIL + "cp %s gs://%s/pipe/%s/%s/%s"
        cmd = cmd % (f, env["RUNTIME_WORK_BUCKET"], pipe, env["GCP_NAME"], f.basename())
        ret = run_app(cmd, env=env)
        log_run_outerr(ret)
    return(f.basename())
