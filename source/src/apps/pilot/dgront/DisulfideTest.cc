// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * ShowInput.cc
 *
 *  Created on: Jan 26, 2009
 *      Author: dgront
 */

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <protocols/Protocol.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer tr("DisulfideTest");

void register_options() {

  OPT( in::file::native );
}

void compute_chi(pose::Pose &pose) {
  //this little method demonstrates access to CHI-angles, atom-coordinates, Atom-Tree, Stubs
  // what we are going to do is to compute the position of SG on CYS outside of the pose
  // we use the same mechanism as deep inside the atomtree, i.e., we get the stub on CB and transform this stub to get SG.
  // we copy the internal coordinates d, theta and the chi-angle from the pose and then compare with the pose-computed position of SG.
  using namespace core;
  using namespace id;
  using namespace numeric;
  using namespace kinematics;
  //      using namespace conformation;
  // ATTENTION no error checking in this code. If you have an invalid AtomID you are screwed.

  Size const seqpos(286);//a CYS

  id::NamedAtomID sulfur_id("SG", seqpos);
  core::kinematics::tree::Atom const& sulfur(pose.conformation().atom_tree().atom(id::AtomID(sulfur_id, pose)));

  id::DOF_Type dof_id_theta(THETA); //the angle SG-CB-CA
  id::DOF_Type dof_id_d(D); //the distance SG-CB

  core::Real const theta( //getting the angle of SG
      sulfur.dof(dof_id_theta));

  //three different ways of getting the S - CB distance.
  core::Real const S_CB_dist(sulfur.distance(*sulfur.parent())); //I know CB is parent in Atom-Tree ( BAD BEHAVIOUR )
  core::Real const S_CB_dist2(sulfur.dof(dof_id_d)); // Thats the better way of doing it. BUT STILL BAD !
  core::Real const S_CB_dist3(distance(pose.xyz(sulfur_id), pose.xyz(NamedAtomID("CB", seqpos)))); //GOOD BEHAVIOUR

  tr.Info << "theta " << theta << " dist1 " << S_CB_dist << " dist2 " << S_CB_dist2 << std::endl;

  //the STUB of CB is where we start
  id::NamedStubID CB_stub_id("CB", "CA", "N", seqpos); //sequence center-atom, parent, grand-parent
  kinematics::Stub manual_stub(pose.stub_from_id(CB_stub_id)); //GOOD WAY OF GETTING IT
  //just to demonstrate how the atom-tree works the BAD WAY OF GETTING THE STUB
  kinematics::Stub atree_stub( //get the same stub directly from the AtmoTree --- asuming N2C folding direction
      pose.conformation().atom_tree().atom(id::AtomID(id::NamedAtomID("CB", seqpos), pose)).get_stub());

  tr.Info << "\natree Stub " << atree_stub << "\n" << "manual     " << manual_stub << std::endl;
  //get the current chi-angle that determines the position of SG ( chi1 ).
  TorsionID chi1_id(seqpos, CHI, 1);
  core::Real const chi(pose.torsion(chi1_id));
  tr.Info << "chi1 " << chi << std::endl;

  //now compute the new STUB for SG. the center of the new stub is the position of SG.
  manual_stub.M *= x_rotation_matrix_degrees(chi);
  Stub new_stub(manual_stub.M * z_rotation_matrix_radians(theta), manual_stub.v);
  new_stub.v += S_CB_dist * new_stub.M.col_x(); //new center of stub -- should coincide with sulfur

  Vector sulfur_xyz(pose.xyz(sulfur_id));
  tr.Info << " compare to pose.xyz " << distance(sulfur_xyz, new_stub.v) << std::endl;

}

int main(int argc, char * argv[]) {
try {
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  register_options();
  devel::init(argc, argv);

  pose::Pose init_pose;

  if (option[in::file::native].user()) {
    core::import_pose::pose_from_pdb(init_pose, option[in::file::native]());
    //    pose::set_ss_from_phipsi(init_pose);

    init_pose.dump_pdb("input.pdb", "");
    std::cout << init_pose.fold_tree() << std::endl;
    utility::vector1<int> cuts = init_pose.fold_tree().cutpoints();
    for (Size i = 1; i <= cuts.size(); i++)
      std::cout << cuts[i] << " ";
    std::cout << std::endl;
  }
} catch ( utility::excn::EXCN_Base const & e ) {
                          std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                              }
    return 0;
}
