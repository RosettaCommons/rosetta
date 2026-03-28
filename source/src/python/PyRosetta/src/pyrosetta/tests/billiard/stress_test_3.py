import billiard
import functools
import os
import random
import tempfile
import time
import traceback
import unittest

from distributed import Client, LocalCluster, as_completed


def propagate_errors(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        errq = kwargs.pop("_errq", None)
        try:
            return func(*args, **kwargs)
        except Exception:
            if errq is not None:
                errq.put(traceback.format_exc())
            raise

    return wrapper


def get_remodel_options():
    return {
        "-remodel:num_trajectory": "1",
        "-remodel:quick_and_dirty": "1",
        "-remodel:dr_cycles": str(random.randint(1, 5)),
        "-remodel:use_clusters": str(random.choice([0, 1])),
        "-remodel:core_cutoff": str(random.randint(10, 20)),
        "-remodel:boundary_cutoff": str(random.randint(5, 15)),
        "-remodel:generic_aa": random.choice(list("ACDEFGHIKLMNPQRSTVWY")),
        "-remodel:hb_srbb": str(random.choice([0.0, 1.0])),
        "-remodel:vdw": str(random.random()),
        "-remodel:rama": str(random.random()),
        "-remodel:cbeta": str(random.random()),
        "-remodel:cenpack": str(random.random()),
        "-remodel:rg": str(random.random()),
        "-ex1": str(random.choice([0, 1])),
        "-ex2": str(random.choice([0, 1])),
        "-ex3": str(random.choice([0, 1])),
        "-ex4": str(random.choice([0, 1])),
        "-ex1:level": str(random.choice([0, 1])),
        "-ex2:level": str(random.choice([0, 1])),
        "-ex3:level": str(random.choice([0, 1])),
        "-ex4:level": str(random.choice([0, 1])),
        "-extrachi_cutoff": str(random.choice([0, 18])),
        "-use_input_sc": str(random.choice([0, 1])),
    }


@propagate_errors
def make_pose_and_serialize(q):
    import pyrosetta
    import pyrosetta.distributed.io as io

    pid = os.getpid()
    print(f"[PID {pid}] initializing 'beta_nov16'...", flush=True)
    pyrosetta.init("-beta_nov16 1 -out:level 0 -multithreading:total_threads 1", silent=True)

    pose = pyrosetta.pose_from_sequence("AAAA")

    # Score
    sfxn = pyrosetta.get_score_function()
    print(f"[PID {pid}] score (beta_nov16):", sfxn(pose), flush=True)

    # Serialize Pose
    packed = io.to_packed(pose)

    print(f"[PID {pid}] sending PackedPose", flush=True)
    q.put(packed)

    time.sleep(0.5)


@propagate_errors
def receive_and_remodel(q):
    import pyrosetta
    import pyrosetta.distributed
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

    pid = os.getpid()
    random.seed(pid)

    print(f"[PID {pid}] initializing 'beta_jan25'...", flush=True)
    flags = get_remodel_options()
    flags["-beta_jan25"] = "1"
    flags["-out:level"] = "200"
    flags["-multithreading:total_threads"] = "1"
    options = pyrosetta.distributed._normflags(flags)
    pyrosetta.init(options=options, silent=True)

    print(f"[PID {pid}] Initialized:", options, flush=True)

    packed = q.get()
    print(f"[PID {pid}] received PackedPose", flush=True)

    # Deserialize Pose
    pose = packed.pose

    # Score before Remodel
    sfxn = pyrosetta.get_score_function()
    print(f"[PID {pid}] Score (beta_jan25, pre-Remodel):", sfxn(pose), flush=True)

    # Create temporary blueprint file
    with tempfile.TemporaryDirectory() as tmp_dir:
        bp = os.path.join(tmp_dir, "tmp.bp")

        with open(bp, "w") as f:
            for i in range(1, pose.size() + 1):
                ss = "H" if i > (pose.size() - 2) else "."
                f.write(f"{i} A {ss} NATAA\n")
            f.write("0 X H PIKAA A")  # Simple extension, no new line at end

        xml = XmlObjects.create_from_string(
            f"""<MOVERS>
                <RemodelMover name="remodel" blueprint="{bp}"/>
            </MOVERS>"""
        ).get_mover("remodel")

        print(f"[PID {pid}] running Remodel...", flush=True)
        xml.apply(pose)

    # Score after Remodel
    print(f"[PID {pid}] Score (beta_jan25, post-Remodel):", sfxn(pose), flush=True)


def run_once(r):
    print(f"--- Round {r} ---", flush=True)

    ctx = billiard.get_context("spawn")
    q = ctx.Queue()
    errq = ctx.Queue()

    p1 = ctx.Process(
        target=make_pose_and_serialize,
        args=(q,),
        kwargs={"_errq": errq},
    )
    p2 = ctx.Process(
        target=receive_and_remodel,
        args=(q,),
        kwargs={"_errq": errq},
    )

    p1.start()
    p2.start()

    p1.join()
    p2.join()

    errors = []
    while not errq.empty():
        errors.append(errq.get())

    if errors:
        raise RuntimeError("Subprocess error:\n" + "\n".join(errors))
    if p1.exitcode != 0:
        raise RuntimeError(f"Function `make_pose_and_serialize` failed (exitcode={p1.exitcode})")
    if p2.exitcode != 0:
        raise RuntimeError(f"Function `receive_and_remodel` failed (exitcode={p2.exitcode})")

    return r


class StressTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = self.tmpdir.name
        self.cwd = os.getcwd()
        self.remodel_1_pdb_file = os.path.join(self.cwd, "1.pdb")
        self.assertFalse(os.path.isfile(self.remodel_1_pdb_file))
        os.chdir(self.workdir) # Prevent Remodel dumping '1.pdb' in current working directory

    def tearDown(self):
        if os.path.isfile(self.remodel_1_pdb_file):
            try:
                os.remove(self.remodel_1_pdb_file)
            except Exception:
                print(f"Warning: Remodel output file still exists: {self.remodel_1_pdb_file}")
        os.chdir(self.cwd)
        self.tmpdir.cleanup()

    def test_stress(self):
        rounds = 100
        n_workers = 8
        threads_per_worker = 1
        with LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker) as cluster:
            with Client(cluster) as client:
                futures = []
                for r in range(rounds):
                    arg = client.scatter(r, broadcast=False, hash=False)
                    future = client.submit(run_once, arg, pure=False)
                    futures.append(future)
                seq = as_completed(futures)
                for i, future in enumerate(seq, start=1):
                    r = future.result()
                    print(f"Gathered round {r} at step {i}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
