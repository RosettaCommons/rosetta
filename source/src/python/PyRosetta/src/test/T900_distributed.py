from __future__ import print_function
import shutil
import sys
import time


if tuple(sys.version_info) < (3, 6):
    print("Unsupported python version for pyrosetta.distributed: %s" % str( sys.version_info ) )
    sys.exit(0)

else:
    print("sys.version_info: {0}".format(sys.version_info))


import pyrosetta.rosetta
import subprocess

if not hasattr(pyrosetta.rosetta, "cereal"):
    print("Unsupported non-serialization build for pyrosetta.distributed.")
    sys.exit(0)


# check if required packages was installed
try:
    import attrs
    import cloudpickle
    import dask
    import pandas
except ImportError as e:
    print(f"{e}\nSome packages required for pyrosetta.distributed is missing, skipping the tests...")
    sys.exit(0)


if shutil.which("conda"):
    try:
        export = subprocess.check_output(
            "conda env export --prefix $(conda env list | grep '*' | awk '{print $NF}')",
            shell=True,
            text=True,
        )
        print("Current conda environment:", export, sep="\n")
    except subprocess.CalledProcessError as ex:
        print("Printing conda environment failed with return code: {0}.".format(ex.returncode))
else:
    try:
        freeze = subprocess.check_output(
            "{0} -m pip freeze".format(sys.executable),
            shell=True,
            text=True,
        )
        print("Current pip environment:", freeze, sep="\n")
    except subprocess.CalledProcessError as ex:
        print("Printing pip environment failed with return code: {0}.".format(ex.returncode))


def e(cmd, sleep=1):
    """Run command getting return code and output with subsequent sleep step to provide extra time to flush loggers."""
    print("Executing:\n{0}".format(cmd))
    status, output = subprocess.getstatusoutput("{0} && sleep {1}".format(cmd, sleep))
    print("Output:\n{0}".format(output))
    if status != 0:
        print(
            "Encountered error(s) with exit code {0} while running: {1}\nTerminating...".format(
                status, cmd
            )
        )
        sys.exit(1)

distributed_cluster_test_cases = [
    "pyrosetta.tests.distributed.cluster.test_smoke.SmokeTest.test_smoke",
    "pyrosetta.tests.distributed.cluster.test_smoke.SmokeTest.test_ignore_errors",
    "pyrosetta.tests.distributed.cluster.test_smoke.SmokeTestMulti.test_smoke_multi",
    "pyrosetta.tests.distributed.cluster.test_smoke.SmokeTestMulti.test_smoke_multi_from_instance",
    "pyrosetta.tests.distributed.cluster.test_smoke.SaveAllTest.test_save_all",
    "pyrosetta.tests.distributed.cluster.test_smoke.SerializationTest.test_serialization",
    "pyrosetta.tests.distributed.cluster.test_smoke.SaveAllTest.test_save_all_dry_run",
    "pyrosetta.tests.distributed.cluster.test_smoke.ScoresTest.test_detached_scores",
    "pyrosetta.tests.distributed.cluster.test_smoke.ScoresTest.test_detached_scores_in_protocol",
    "pyrosetta.tests.distributed.cluster.test_smoke.ScoresTest.test_detached_scores_with_reserve_scores",
    "pyrosetta.tests.distributed.cluster.test_smoke.MultipleClientsTest.test_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.ResourcesTest.test_resources",
    "pyrosetta.tests.distributed.cluster.test_smoke.ResourcesTest.test_resources_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_builtin_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_multi_user_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_partition_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_user_client",
    "pyrosetta.tests.distributed.cluster.test_smoke.RuntimeTest.test_timing_multi_instance",
    "pyrosetta.tests.distributed.cluster.test_smoke.RuntimeTest.test_timing_single_instance",
    "pyrosetta.tests.distributed.cluster.test_logging.LoggingTest.test_logging",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_minimizer_nstruct",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_nstruct",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_separate",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_from_reproduce",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_decoy_ids",
]
distributed_test_suites = [
    "pyrosetta.tests.bindings.init.test_init_files",
    "pyrosetta.tests.bindings.core.test_pose",
    "pyrosetta.tests.distributed.test_concurrency",
    "pyrosetta.tests.distributed.test_dask",
    "pyrosetta.tests.distributed.test_gil",
    "pyrosetta.tests.distributed.test_smoke",
    "pyrosetta.tests.distributed.test_viewer",
    "pyrosetta.tests.numeric.test_alignment",
]
tests = distributed_cluster_test_cases + distributed_test_suites

for test in tests:
    sleep = 5 if test in distributed_cluster_test_cases else 1
    t0 = time.time()
    e("{python} -m unittest {test}".format(python=sys.executable, test=test), sleep=sleep)
    t1 = time.time()
    dt = t1 - t0
    print("Finished running test in {0} seconds: {1}\n".format(round(dt, 6), test))
