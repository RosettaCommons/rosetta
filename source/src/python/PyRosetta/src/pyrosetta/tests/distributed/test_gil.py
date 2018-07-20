import unittest
import threading
import time
import numpy
import logging

import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.score as score


class HeartBeat(threading.Thread):

    @property
    def beat_intervals(self):
        return self.beats[1:] - self.beats[:-1]

    @property
    def beats(self):
        return numpy.array(self._beats)

    def __init__(self, interval):
        self._beats = []
        self.alive = False
        self.interval = interval

        threading.Thread.__init__(self)

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.alive = False
        self.join()

        return False

    def run(self):
        self._beats = []
        self.alive = True

        while self.alive:
            self.tick()
            time.sleep(self.interval)

    def tick(self):
        t = time.time()
        logging.info("HeartBeat.tick: %s", t)
        self._beats.append(t)


class TestGIL(unittest.TestCase):

    def test_gil_score(self):
        with HeartBeat(10e-3) as hb:
            test_pose = io.pose_from_sequence("TESTTESTTEST" * 10)
            score.ScorePoseTask()(test_pose)

        numpy.testing.assert_allclose(hb.beat_intervals, hb.interval, rtol=1)
        self.assertGreater(len(hb.beats), 4)

    def test_gil_sleep(self):
        with HeartBeat(10e-3) as hb:
            time.sleep(60e-3)

        numpy.testing.assert_allclose(hb.beat_intervals, hb.interval, rtol=1)
        self.assertGreater(len(hb.beats), 4)
