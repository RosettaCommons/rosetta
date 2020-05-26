import sys
try:
    import py3Dmol
except ModuleNotFoundError as ex:
    print("ModuleNotFoundError: {0}. Skipping pyrosetta.distributed.viewer unit test!".format(ex))
    sys.exit(0)

import glob
import os
import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.viewer as viewer
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
            pose = poses[0]
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

            metals_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
                 pyrosetta.rosetta.core.chemical.ResidueProperty(31)
            )
            ligands_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
                 pyrosetta.rosetta.core.chemical.ResidueProperty(2)
            )
            view = viewer.init(poses, window_size=(800, 600)) \
                + viewer.setStyle() \
                + viewer.setStyle(residue_selector=ligands_selector, style="stick", colorscheme="magentaCarbon", radius=0.5) \
                + viewer.setStyle(residue_selector=metals_selector, style="sphere", colorscheme="chainHetatm", radius=1.5)
            view.reset()

            polar_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
                pyrosetta.rosetta.core.chemical.ResidueProperty(52)
            )
            view = viewer.init(packed_poses)
            view.add(viewer.setStyle(radius=0.1))
            view.add(viewer.setStyle(residue_selector=polar_residue_selector, colorscheme="whiteCarbon", radius=0.25, label=False))
            view.add(viewer.setHydrogens(color="white", polar_only=True, radius=0.1))
            view.add(viewer.setHydrogenBonds(color="black"))
            view.add(viewer.setDisulfides(radius=0.1))
            view()
            view.reset()

            view = sum(
                [
                    viewer.init(poses),
                    viewer.setStyle(cartoon=False, style="sphere", radius=1.5, colorscheme="darkgreyCarbon"),
                    viewer.setZoom(factor=0.95)
                ]
            )
            view()
            view.reset()

            command_tuple = {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}}
            command_dict = {"hetflag": True}
            chA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
            chB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")
            view = sum(
                [
                    viewer.init(poses),
                    viewer.setStyle(cartoon_color="lightgrey", radius=0.25),
                    viewer.setSurface(residue_selector=chA, colorscheme="greenCarbon", opacity=0.65, surface_type="VDW"),
                    viewer.setSurface(residue_selector=chB, color="blue", opacity=0.75, surface_type="SAS"),
                    viewer.setDisulfides(radius=0.25),
                    viewer.setZoom(factor=1.5),
                    viewer.setStyle(command=command_tuple),
                    viewer.setStyle(command=command_dict)
                ]
            )
            view()
            view.reset()

            helix_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("H")
            sheet_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("E")
            loop_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("L")
            modules = [
                viewer.setBackgroundColor(color="black"),
                viewer.setStyle(residue_selector=helix_selector, cartoon_color="blue", label=False, radius=0),
                viewer.setStyle(residue_selector=sheet_selector, cartoon_color="red", label=False, radius=0),
                viewer.setStyle(residue_selector=loop_selector, cartoon_color="white", label=False, radius=0),
                viewer.setZoomTo(residue_selector=sheet_selector)
            ]
            viewer.init(poses, window_size=(1200, 600), modules=modules).show()

            view = viewer.init(pose, delay=0.15) \
                + viewer.setStyle(radius=0.1) \
                + viewer.setDisulfides(radius=0.1)
            backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
            minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            for _ in range(3):
                backrub.apply(pose)
                minimize.apply(pose)
                view.show()
            view.reset()

            def myCustomPreset(packed_and_poses_and_pdbs=None, *args, **kwargs):
                """
                Add a description of the preset Viewer here
                """
                # Add custrom ResidueSelectors
                metals_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
                     pyrosetta.rosetta.core.chemical.ResidueProperty(31)
                )
                ligands_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
                     pyrosetta.rosetta.core.chemical.ResidueProperty(2)
                )
                # Add custom Viewer commands
                view = viewer.init(packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs) \
                    + viewer.setBackgroundColor("white") \
                    + viewer.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15) \
                    + viewer.setStyle(residue_selector=ligands_selector, style="stick", colorscheme="brownCarbon", radius=0.5, label=True) \
                    + viewer.setStyle(residue_selector=metals_selector, style="sphere", colorscheme="chainHetatm", radius=1.5, label=True) \
                    + viewer.setHydrogenBonds() \
                    + viewer.setDisulfides(radius=0.15) \
                    + viewer.setHydrogens(color="white", radius=0.033, polar_only=True) \
                    + viewer.setSurface(residue_selector=ligands_selector, surface_type="VDW", opacity=0.5, color="magenta") \
                    + viewer.setSurface(residue_selector=metals_selector, surface_type="VDW", opacity=0.5, color="magenta") \
                    + viewer.setZoomTo(residue_selector=ligands_selector)
                return view()

            myCustomPreset(pose)
