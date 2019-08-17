import glob
import pyrosetta.distributed.io as io
import pyrosetta.distributed.viewer as viewer
import os
import tempfile
import unittest 


class TestViewer(unittest.TestCase):

    with tempfile.TemporaryDirectory() as workdir:

        def setUp(self, local_dir=workdir):

            if not os.path.isdir(local_dir):
                os.mkdir(local_dir)

            poses = [io.pose_from_sequence("TEST" * i) for i in range(1, 4)]
            for i, pose in enumerate(poses, start=1):
                with open(os.path.join(local_dir, "tmp_{0}.pdb".format(i)), "w") as f:
                    f.write(io.to_pdbstring(pose))

        def tearDown(self, local_dir=workdir):

            if os.path.isdir(local_dir):
                pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
                for pdbfile in pdbfiles:
                    os.remove(pdbfile)
                os.rmdir(local_dir)

        def test_viewer_with_pdbfiles(self, local_dir=workdir):

            pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
            viewer.presets.coreBoundarySurface(pdbfiles, continuous_update=True)
            viewer.presets.ligandsAndMetals(pdbfiles, window_size=(200.01, 200.01))
            view = viewer.init(pdbfiles, (1600, 400), delay=0.1234567890) \
                + viewer.setBackgroundColor("black") \
                + viewer.setStyle(style="line", colorscheme="blueCarbon")
            view.show()
            self.assertEqual(view.poses, [None] * len(pdbfiles))
            self.assertEqual(len(view.modules), 2)
            view.reinit()
            self.assertEqual(view.poses, [None] * len(pdbfiles))
            self.assertEqual(len(view.modules), 0)
            view.reset()
            self.assertIsNone(view.poses)
            self.assertIsNone(view.pdbstrings)

        def test_viewer_with_poses(self, local_dir=workdir):

            pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
            packed_poses = [io.pose_from_file(pdbfile) for pdbfile in pdbfiles]
            poses = [io.to_pose(p) for p in packed_poses]
            viewer.presets.coreBoundarySurface(packed_poses, delay=0)
            viewer.presets.ligandsAndMetals(packed_poses, continuous_update=True, window_size=(100., 100.))
            modules = [
                viewer.setBackgroundColor("grey"),
                viewer.setStyle(style="sphere", colorscheme="greenCarbon", radius=1.)
            ]
            view = sum([viewer.init(poses)] + modules)
            view()
            self.assertEqual(view.poses, poses)
            self.assertEqual(len(view.modules), 2)
            view.reinit()
            self.assertEqual(view.poses, poses)
            self.assertEqual(len(view.modules), 0)
            view.modules = modules
            self.assertEqual(len(view.modules), len(modules))
            view.reset()
            self.assertIsNone(view.poses)
            self.assertIsNone(view.pdbstrings)
            view = viewer.init(poses, modules=modules)
            self.assertEqual(len(view.modules), len(modules))
            view.reset()
            self.assertIsNone(view.modules)
