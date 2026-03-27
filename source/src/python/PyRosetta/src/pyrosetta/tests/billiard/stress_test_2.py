import multiprocessing as mp
import os
import random
import time
import unittest


FLAGS = [
    "-beta_nov16",
    "-beta_nov16_cart",
    "-score:weights ref2015",
    "-restore_talaris_behavior 1",
    # "-beta_july15 1",
    # "-beta_july15_cart 1",
    "-gen_potential 1",
    "-beta_jan25"
]


def worker(proc_id, barrier):
    import pyrosetta

    flag = FLAGS[proc_id % len(FLAGS)]
    opts = f"{flag} -out:level 0 -multithreading:total_threads 1"

    barrier.wait()
    pyrosetta.init(opts, silent=True)
    print(f"[PID {os.getpid()}] init with: {opts}", flush=True)

    # minimal work
    for _ in range(200):
        pose = pyrosetta.pose_from_sequence("AAAA")
        sfxn = pyrosetta.get_score_function()
        sfxn(pose)

    global_store = []
    for _ in range(200):
        pose = pyrosetta.pose_from_sequence("AAAA")
        global_store.append(pose)

    # keep process alive briefly to overlap teardown
    time.sleep(random.uniform(0.5, 1.5))


class StressTest(unittest.TestCase):
    def test_stress(self):
        for round in range(10):  # loop to amplify probability
            print(f"\n--- round {round} ---", flush=True)

            ctx = mp.get_context("spawn")

            n_procs = os.cpu_count() * 3  # increase pressure
            barrier = ctx.Barrier(n_procs)

            procs = []
            for i in range(n_procs):
                p = ctx.Process(target=worker, args=(i, barrier))
                p.start()
                procs.append(p)

            for p in procs:
                p.join()


if __name__ == "__main__":
    unittest.main(verbosity=2)
