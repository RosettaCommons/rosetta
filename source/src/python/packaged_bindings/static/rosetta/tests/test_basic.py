import unittest

class TestInitialization(unittest.TestCase):
    def test_import(self):
        import rosetta

    def test_init(self):
        import rosetta
        rosetta.init()

    def test_pose_creation(self):
        import rosetta
        rosetta.init()
        test_pose = rosetta.pose_from_sequence("TESTTEST")
