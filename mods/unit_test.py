from moke import * #@UnusedWildImport
from . import * #@UnusedWildImport
from string import Template

@task
def unit_test(test='HCC2218', start_from='align', skip_r_tests=False, skip_cleanup=False, dry_run=False):
    """unit test TPO

    - test (``str``) what unit test sample to test
    - start_from (``str``) which pipeline to start from
    - dry_run (``bool``) save rather than run commands
    - skip_r_tests (``bool``) Skip the R test
    - skip_cleanup (``bool``) Skip removing alignments
    """
    env["CONFIG_FILE"] = dict(task._funcparse(None).parse_args().__dict__).get("config").name
    env["SKIP_R_TESTS"] = 'true' if skip_r_tests else 'false'
    env["SKIP_CLEANUP"] =  'true' if skip_cleanup else 'false'
    env["DRYRUN"] = 'true' if dry_run else 'false'
    env["START_FROM"] = start_from
    env["TEST_CFG"] = test
    ## config file variables
    env["REFS_VER"] = env["TPO_REFS_VER"]
    env["UNIT_TEST_BUCKET"] = env["TEST_UNIT_TEST_BUCKET"]
    env["UNIT_TEST_DIR"] = env["TEST_UNIT_TEST"]
    
    ret = run_app(env["ROOT"] / "pipe/unit-test/unit_test.sh", [], env=env)
    log_run_outerr(ret)
