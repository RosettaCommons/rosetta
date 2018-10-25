import unittest


def test():
    return unittest.TextTestRunner().run(unittest.defaultTestLoader.discover("pyrosetta.tests"))
