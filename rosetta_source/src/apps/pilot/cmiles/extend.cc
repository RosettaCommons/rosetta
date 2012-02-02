// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

using namespace std;
using namespace core;
using namespace core::kinematics;
using namespace core::pose;
using namespace protocols::loops;

void generate_extended_pose(Pose* extended_pose, const string& sequence) {
  core::pose::make_pose_from_sequence(*extended_pose, sequence, *(chemical::ChemicalManager::get_instance()->residue_type_set(chemical::CENTROID)));

  for (Size i = 1; i <= extended_pose->total_residue(); ++i) {
    extended_pose->set_phi(i, -150);
    extended_pose->set_psi(i, 150);
    extended_pose->set_omega(i, 180);
  }
}

void copy_residues(const Pose& src, Size start, Size stop, Pose* dst) {
  using namespace core::conformation;
  using namespace core::id;

  for (Size i = start; i <= stop; ++i) {
    const Residue& r = src.conformation().residue(i);
    const Residue& s = dst->conformation().residue(i);

    for (Size j = 1; j <= r.natoms(); ++j) {
      AtomID id(j, i);
      dst->set_xyz(id, src.xyz(id));
    }
  }
}

int main(int argc, char* argv[]) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  devel::init(argc, argv);

  // structure to borrow from
  PoseOP reference = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
  core::util::switch_to_residue_type_set(*reference, core::chemical::CENTROID);
  reference->dump_pdb("reference.pdb");

  Pose pose;
  generate_extended_pose(&pose, reference->sequence());
  pose.dump_pdb("pose_extended.pdb");

  // bring structures into same coordinate frame
  core::scoring::calpha_superimpose_pose(pose, *reference);
  pose.dump_pdb("pose_superimposed.pdb");

  // residues to borrow from reference
  Loop f1(48, 60), f2(157, 169);

  // right-to-left propagation
  FoldTree r2l_tree;
  r2l_tree.add_edge(pose.total_residue(), 1, core::kinematics::Edge::PEPTIDE);
  pose.fold_tree(r2l_tree);
  copy_residues(*reference, f1.start(), f1.stop(), &pose);
  core::conformation::idealize_position(f1.start() - 1, pose.conformation());  // Attach residues preceding f1.start() to f1

  // left-to-right propagation
  FoldTree l2r_tree;
  l2r_tree.add_edge(1, pose.total_residue(), core::kinematics::Edge::PEPTIDE);
  pose.fold_tree(l2r_tree);

  core::conformation::idealize_position(f1.stop(), pose.conformation());  // Attach residues following f1.stop() to f1

  copy_residues(*reference, f2.start(), f2.stop(), &pose);
  core::conformation::idealize_position(f2.stop(), pose.conformation());  // Attach residues following f2.stop() to f2

  pose.dump_pdb("pose_final.pdb");
  return 0;
}
