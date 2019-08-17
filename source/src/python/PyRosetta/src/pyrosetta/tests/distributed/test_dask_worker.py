# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import contextlib
import os
import subprocess
import sys
import tempfile
import time
import unittest


class TestDaskArgs(unittest.TestCase):

    def test_worker_extra(self):
        """worker_extra worker plugin controls rosetta flags and local dir.
        pyrosetta.distributed.dask.worker_extra can be used to specify a set of
        PyRosetta initialization flags and local directory information that will be used to
        pre-initialize all worker processes in a dask cluster.
        This provides support for protocols which may require command line
        flags including (e.g. additional residue types, logging flags, ...).
        """

        # Import locally to allow test discovery if dependencies aren't present
        from dask import delayed
        from dask.distributed import Client, LocalCluster

        import pyrosetta
        import pyrosetta.distributed.dask
        import pyrosetta.distributed.io as io
        import pyrosetta.distributed.tasks.score as score

        # Setup cluster scheduler and working directory
        with tempfile.TemporaryDirectory() as workdir, LocalCluster(
            n_workers=0, diagnostics_port=None
        ) as cluster:

            # Context manager controls launch & teardown of test worker
            @contextlib.contextmanager
            def one_worker(init_flags):
                init = pyrosetta.distributed.dask.worker_extra(
                    init_flags=init_flags, local_directory=workdir
                )

                worker_command = (
                    "%s -m distributed.cli.dask_worker %s --nthreads 1 %s"
                    % (sys.executable, cluster.scheduler_address, " ".join(init)),
                )
                worker = subprocess.Popen(worker_command, shell=True)

                try:
                    wait = 10.0
                    while len(cluster.scheduler.workers) == 0:
                        worker.poll()
                        assert worker.returncode is None
                        time.sleep(.1)
                        wait -= .1
                        assert wait > 0, "Timeout waiting for worker launch."

                    yield None
                finally:
                    worker.terminate()
                    try:
                        worker.wait(1)
                    except TimeoutError:
                        worker.kill()
                        worker.wait()

            HBI_fa_params_path = os.path.join(workdir, "HBI.fa.params")

            with open(HBI_fa_params_path, "w") as of:
                of.write(_HBI_fa_params)

            # Unformatted PyRosetta command line flags
            flags = """
               -extra_res_fa   {0}   #Test flag comment
            -ignore_unrecognized_res 1   -use_input_sc  1
             -ex4 ### Test flag comment
            """.format(
                HBI_fa_params_path
            )

            # Initialize PyRosetta
            pyrosetta.distributed.dask.init_notebook(flags)

            # Setup protocol
            protocol = score.ScorePoseTask("ref2015_cart")

            # Setup pose
            pose = io.pose_from_sequence("TESTING [HBI]")

            # Score locally
            local_result = protocol(pose)

            # Score on dask-worker with correct command-line flags
            with one_worker(flags), Client(cluster):
                cluster_result = delayed(protocol)(pose).compute()

            self.assertEqual(
                cluster_result.scores["total_score"],
                local_result.scores["total_score"]
            )


# Full-atom ligand .params file not enabled by default in the Rosetta
# database residue types set
_HBI_fa_params = """NAME HBI
IO_STRING HBI Z
TYPE LIGAND
AA UNK
ATOM  C4  aroC  X   0.04
ATOM  C2  aroC  X   0.03
ATOM  N1  Nhis  X   -0.79
ATOM  C11 aroC  X   0.70
ATOM  C1  CH3   X   -0.53
ATOM  H2  Hapo  X   0.22
ATOM  H3  Hapo  X   0.19
ATOM  H1  Hapo  X   0.19
ATOM  N2  Npro  X   -1.10
ATOM  C3  CH3   X   -0.17
ATOM  H5  Hapo  X   0.17
ATOM  H6  Hapo  X   0.17
ATOM  H4  Hapo  X   0.22
ATOM  C12 CNH2  X   1.08
ATOM  O2  ONH2  X   -0.83
ATOM  C5  aroC  X   -0.26
ATOM  C6  aroC  X   -0.15
ATOM  C8  aroC  X   0.35
ATOM  F1  F     X   -0.46
ATOM  C10 aroC  X   0.38
ATOM  O1  OOC   X   -0.83
ATOM  C9  aroC  X   0.35
ATOM  C7  aroC  X   -0.16
ATOM  H7  Haro  X   0.18
ATOM  F2  F     X   -0.46
ATOM  H8  Haro  X   0.24
ATOM  H9  Haro  X   0.21
BOND_TYPE  C1   H2  1   
BOND_TYPE  C1   H3  1   
BOND_TYPE  C1   H1  1   
BOND_TYPE  C1   C11 1   
BOND_TYPE  N1   C11 2   
BOND_TYPE  N1   C2  1   
BOND_TYPE  O1   C10 1   
BOND_TYPE  C2   C4  2   
BOND_TYPE  C2   C12 1   
BOND_TYPE  N2   C11 1   
BOND_TYPE  N2   C3  1   
BOND_TYPE  N2   C12 4   
BOND_TYPE  O2   C12 2   
BOND_TYPE  C3   H5  1   
BOND_TYPE  C3   H6  1   
BOND_TYPE  C3   H4  1   
BOND_TYPE  C4   C5  1   
BOND_TYPE  C4   H9  1   
BOND_TYPE  C5   C6  4   
BOND_TYPE  C5   C7  4   
BOND_TYPE  C6   C8  4   
BOND_TYPE  C6   H8  1   
BOND_TYPE  C7   C9  4   
BOND_TYPE  C7   H7  1   
BOND_TYPE  C8   F1  1   
BOND_TYPE  C8   C10 4   
BOND_TYPE  C9   C10 4   
BOND_TYPE  C9   F2  1   
NBR_ATOM  C4 
NBR_RADIUS 5.662663
ICOOR_INTERNAL    C4     0.000000    0.000000    0.000000   C4    C2    N1 
ICOOR_INTERNAL    C2     0.000000  180.000000    1.349492   C4    C2    N1 
ICOOR_INTERNAL    N1     0.000000   52.312127    1.416411   C2    C4    N1 
ICOOR_INTERNAL    C11 -179.876986   71.981668    1.292864   N1    C2    C4 
ICOOR_INTERNAL    C1  -177.823063   57.141926    1.521803   C11   N1    C2 
ICOOR_INTERNAL    H2    -0.677983   56.211953    1.209650   C1    C11   N1 
ICOOR_INTERNAL    H3  -116.300960   67.773410    1.140605   C1    C11   H2 
ICOOR_INTERNAL    H1  -125.935649   65.293694    1.136383   C1    C11   H3 
ICOOR_INTERNAL    N2   177.824660   68.608152    1.418583   C11   N1    C1 
ICOOR_INTERNAL    C3   179.909099   54.484510    1.445087   N2    C11   N1 
ICOOR_INTERNAL    H5   -59.800186   70.467375    1.110117   C3    N2    C11
ICOOR_INTERNAL    H6   119.675359   70.450405    1.110181   C3    N2    H5 
ICOOR_INTERNAL    H4   120.174433   68.846291    1.110631   C3    N2    H6 
ICOOR_INTERNAL    C12 -179.933962   72.728853    1.340776   N2    C11   C3 
ICOOR_INTERNAL    O2   179.964406   54.938317    1.222126   C12   N2    C11
ICOOR_INTERNAL    C5     0.003093   52.791761    1.497704   C4    C2    N1 
ICOOR_INTERNAL    C6    -0.352261   55.946592    1.413198   C5    C4    C2 
ICOOR_INTERNAL    C8   179.949796   59.531212    1.399355   C6    C5    C4 
ICOOR_INTERNAL    F1  -179.990868   60.282211    1.356833   C8    C6    C5 
ICOOR_INTERNAL    C10  179.985408   59.719418    1.396930   C8    C6    F1 
ICOOR_INTERNAL    O1  -179.980291   59.936638    1.344611   C10   C8    C6 
ICOOR_INTERNAL    C9  -179.990699   60.126809    1.396219   C10   C8    O1 
ICOOR_INTERNAL    C7    -0.021432   59.905469    1.397765   C9    C10   C8 
ICOOR_INTERNAL    H7  -179.979744   61.128745    1.083448   C7    C9    C10
ICOOR_INTERNAL    F2  -179.975061   59.878212    1.356616   C9    C10   C7 
ICOOR_INTERNAL    H8   179.987853   57.836467    1.072807   C6    C5    C8 
ICOOR_INTERNAL    H9   179.889455   63.618049    1.087136   C4    C2    C5
"""
