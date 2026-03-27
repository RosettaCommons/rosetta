import billiard
import os
import pyrosetta
import pyrosetta.distributed.io as io
import queue
import random
import tempfile
import time
import traceback
import unittest

from distributed import Client, LocalCluster, as_completed
from pyrosetta.distributed.cluster import Serialization

CONSTANT_FLAGS = "-out:level 100 -multithreading:total_threads 1"
COMPRESSION = "xz"


def score_function_is_available(name):
    if not pyrosetta.rosetta.basic.was_init_called():
        raise RuntimeError("PyRosetta must be initialized.")
    if not isinstance(name, str):
        raise ValueError("Score function name must be a `str` object.")
    if not name.endswith(".wts"):
        name += ".wts"

    db = pyrosetta.rosetta.basic.database.full_name("scoring/weights")

    return os.path.isfile(os.path.join(db, name))


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


def write_remodel_blueprint(pose, workdir, n_res_cterm_ext=4, n_term_threshold=2):
    line_template = "{res} {name1} {ss}{resfile_cmd}"
    lines = []
    for res in range(1, pose.size() + 1):
        name1 = pose.residue(res).name1()
        ss = random.choice(list("HEL")) if res > max((pose.size() - n_term_threshold), 0) else "."
        resfile_cmd = " NATAA" if ss != "." else ""
        line = line_template.format(res=res, name1=name1, ss=ss, resfile_cmd=resfile_cmd)
        lines.append(line)
    for _ in range(n_res_cterm_ext):
        res = 0
        name1 = "X"
        ss = random.choice(list("HEL"))
        resfile_cmd = " PIKAA " + random.choice(list("ACDEFGHIKLMNPQRSTVWY"))
        line = line_template.format(res=res, name1=name1, ss=ss, resfile_cmd=resfile_cmd)
        lines.append(line)
    blueprint_file = os.path.join(workdir, "tmp.bp")
    with open(blueprint_file, "w") as f:
        f.write("\n".join(lines))

    return blueprint_file


def target(r, proc_id, q, sfxn_flags):
    import pyrosetta
    import pyrosetta.distributed
    import pyrosetta.distributed.io as io
    from pyrosetta.distributed.cluster import Serialization
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

    prefix = f"[Process ID: {os.getpid()}, Number: {proc_id}]"

    try:
        random.seed(proc_id)
        sfxn_flag = sfxn_flags[proc_id % len(sfxn_flags)]
        remodel_flags = pyrosetta.distributed._normflags(get_remodel_options())
        options = " ".join([sfxn_flag, remodel_flags, CONSTANT_FLAGS])
        pyrosetta.distributed.maybe_init(options=options, extra_options="", silent=True)
        seed = pyrosetta.rosetta.numeric.random.rg().get_seed()

        # Unpack pose
        serializer = Serialization(compression=COMPRESSION)
        if r is None:
            seq = "".join(random.sample("ACDEFGHIKLMNPQRSTVWY", k=random.randint(2, 7)))
            pose = pyrosetta.pose_from_sequence(seq)
        else:
            pose = serializer.decompress_packed_pose(r).pose

        # Do minimal work
        scorefxn_fa = pyrosetta.get_score_function()
        scorefxn_fa(pose)
        # Make Remodel blueprint file
        tmp_path = tempfile.TemporaryDirectory()
        n_res_cterm_ext = random.randint(3, 5)
        n_term_threshold = random.randint(2, 3)
        blueprint_file = write_remodel_blueprint(
            pose,
            tmp_path.name,
            n_res_cterm_ext=n_res_cterm_ext,
            n_term_threshold=n_term_threshold,
        )
        # Run Remodel
        xml_obj = XmlObjects.create_from_string(
            f"""<MOVERS><RemodelMover name="remodel" blueprint="{blueprint_file}"/></MOVERS>"""
        ).get_mover("remodel")
        xml_obj.apply(pose)
        # Clean up
        tmp_path.cleanup()
        # Score
        total_score = scorefxn_fa(pose)
        print(f"{prefix}, Seed: {seed}, Total Score: {total_score}, Initialized: {options}", flush=True)
        # Add counter
        if "counter" not in pose.cache.all_keys:
            pose.cache["counter"] = 0
        else:
            pose.cache["counter"] += 1
        # Put result
        r = serializer.compress_packed_pose(io.to_packed(pose))
        q.put(r)

        # Keep process alive briefly to overlap teardown
        time.sleep(random.uniform(0, 0.4))

    except Exception as ex:
        print(f"{prefix} Exception: {ex}", flush=True)
        print(traceback.format_exc(), flush=True)


def run_once(r, sfxn_flags, round, n_procs):
    print(f"--- Round {round} ---", flush=True)

    context = billiard.get_context("spawn")
    q = context.Queue()

    procs = []
    for proc_id in range(n_procs):
        p = context.Process(target=target, args=(r, proc_id, q, sfxn_flags))
        p.start()
        procs.append(p)

    # Drain while running
    alive = True
    while alive:
        alive = any(p.is_alive() for p in procs)
        try:
            while True:
                r = q.get_nowait()
        except queue.Empty:
            pass
        time.sleep(0.05)

    for p in procs:
        p.join()

    # Final drain
    try:
        while True:
            r = q.get_nowait()
    except queue.Empty:
        pass

    return r


class StressTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = self.tmpdir.name
        self.cwd = os.getcwd()
        self.remodel_1_pdb_file = os.path.join(self.cwd, "1.pdb")
        self.assertFalse(os.path.isfile(self.remodel_1_pdb_file))
        os.chdir(self.workdir) # Prevent RosettaRemodel dumping '1.pdb' in current working directory

    def tearDown(self):
        if os.path.isfile(self.remodel_1_pdb_file):
            try:
                os.remove(self.remodel_1_pdb_file)
            except Exception:
                print(f"Warning: RosettaRemodel output file still exists: {self.remodel_1_pdb_file}")
        os.chdir(self.cwd)
        self.tmpdir.cleanup()

    def test_rounds(self, rounds=100):

        pyrosetta.init(options=CONSTANT_FLAGS, extra_options="-run:constant_seed 1 -score:weights ref2015", silent=True)

        sfxn_flags = [
            "-beta_nov16 1",
            "-beta_nov16_cart 1",
            "-score:weights ref2015",
            "-restore_talaris_behavior 1",
            "-gen_potential 1",
        ]
        if score_function_is_available("beta_jan25"):
            sfxn_flags.append("-beta_jan25 1")
        print("Available scorefunction flags:", sfxn_flags)

        serializer = Serialization(compression=COMPRESSION)

        n_procs = 8
        n_workers = 4
        threads_per_worker = 1
        with LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker) as cluster:
            with Client(cluster) as client:
                # Submit futures
                futures = []
                for round in range(rounds):
                    args = client.scatter((None, sfxn_flags, round, n_procs), broadcast=False, hash=False)
                    future = client.submit(run_once, *args, pure=False)
                    futures.append(future)
                seq = as_completed(futures)
                # Gather results
                for i, future in enumerate(seq, start=1):
                    r = future.result()
                    # Deserialize
                    pose = serializer.decompress_packed_pose(r).pose
                    if pose.cache["counter"] >= 1:
                        continue
                    # Serialize
                    r = serializer.compress_packed_pose(io.to_packed(pose))
                    # Submit future
                    args = client.scatter((r, sfxn_flags, rounds + i, n_procs), broadcast=False, hash=False)
                    future = client.submit(run_once, *args, pure=False)
                    seq.add(future)


if __name__ == "__main__":
    unittest.main(verbosity=2)
