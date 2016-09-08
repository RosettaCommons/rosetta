// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file    apps/public/design/zinc1_homodimer_setup.cc
/// @brief   Grafts a 2-residue + zinc match onto the match's scaffold, duplicates the pose, then performs two discrete flips (rollmoves) to the second chain that generates a symmetric pose with a tetrahedral metal binding site.
/// @details Takes a 2-residue + zinc match as -s or -l, requires the scaffold pdb indicated as an option.  The first flip (rollmove) is along the axis that connects p1 with p2 (p1 is an xyzVector one zinc-coordinating atom), and zinc is translated to the origin at the time of the flip.  The second rollmove is along the axis that bisects the p1:p2 and p3:p4.  This protocol preceeds SymMetalInterface_OneZN_design.cc
    /////////////////////parallel///////////////
    //         p3      p1
    // ChainB      Zn      ChainA
    //         p4      p2
    ////////////////////////////////////////////
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/metal_interface/MatchGrafter.hh>
#include <protocols/metal_interface/FindClosestAtom.hh>
#include <protocols/rigid/RollMover.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <protocols/jd2/JobDistributor.hh>


//tracers
using basic::Error;
using basic::Warning;
//using numeric::conversions::degrees;

static THREAD_LOCAL basic::Tracer TR( "apps.public.design.zinc1_homodimer_setup" );

typedef numeric::xyzVector<core::Real> point;
typedef point axis;

basic::options::StringOptionKey const scaffold_pdb("scaffold_pdb");

/// @brief
class zinc1_homodimer_setup : public protocols::moves::Mover {
public:
  zinc1_homodimer_setup()
  {
  }
  virtual ~zinc1_homodimer_setup(){};


  virtual
  void
  apply( core::pose::Pose & match ){

    assert(match.size() == 3); // Two ligands plus zinc

    core::pose::Pose scaffold;
    core::import_pose::pose_from_file( scaffold, basic::options::option[scaffold_pdb].value() , core::import_pose::PDB_file);


    ///////////////////////////Graft match onto scaffold///////////////////////////////////////////////////////////

    protocols::metal_interface::MatchGrafterOP match_grafter = new protocols::metal_interface::MatchGrafter;
    core::pose::Pose scaffold_with_match = match_grafter->graft( match, scaffold );

    core::pose::Pose const homodimer_with_matches = match_grafter->build_combined_pose_with_zinc_overlay( scaffold_with_match, scaffold_with_match );
    core::pose::Pose homodimer_A( homodimer_with_matches );
    core::pose::Pose homodimer_B( homodimer_with_matches );


    ///////////////////////////Need the atom_xyz's to define axes of rotation, which requires the seqpos and atom name ///////////////////////

    utility::vector1< core::Size > metalsite_seqpos( 5, 0 ); //safer when using pushback
    utility::vector1< point > metalsite_atom_xyz( 5 );
    utility::vector1< std::string > metalsite_atom_name( 5, "" );

    core::Size match_res1 = match.pdb_info()->number(1); // gets residue number according to .pdb file, not renumbered starting at 1
    core::Size match_res2 = match.pdb_info()->number(2);
    metalsite_seqpos[1] = match_res1;
    metalsite_seqpos[2] = match_res2;
    metalsite_seqpos[3] = scaffold.size() + 1/*zinc is chain B*/ + match_res1;
    metalsite_seqpos[4] = scaffold.size() + 1/*zinc is chain B*/ + match_res2;
    metalsite_seqpos[5] = scaffold.size() + 1/*zinc is chain B*/;

    metalsite_atom_xyz[5] = homodimer_with_matches.residue( metalsite_seqpos[5] ).atom(1).xyz();

    metalsite_atom_name[1] = protocols::metal_interface::find_closest_atom( homodimer_with_matches.residue(metalsite_seqpos[1]), metalsite_atom_xyz[5] );
    metalsite_atom_name[2] = protocols::metal_interface::find_closest_atom( homodimer_with_matches.residue(metalsite_seqpos[2]), metalsite_atom_xyz[5] );
    metalsite_atom_name[3] = protocols::metal_interface::find_closest_atom( homodimer_with_matches.residue(metalsite_seqpos[3]), metalsite_atom_xyz[5] );
    metalsite_atom_name[4] = protocols::metal_interface::find_closest_atom( homodimer_with_matches.residue(metalsite_seqpos[4]), metalsite_atom_xyz[5] );
    metalsite_atom_name[5] = "ZN";

    metalsite_atom_xyz[1] = homodimer_with_matches.residue( metalsite_seqpos[1] ).atom( metalsite_atom_name[1] ).xyz();
    metalsite_atom_xyz[2] = homodimer_with_matches.residue( metalsite_seqpos[2] ).atom( metalsite_atom_name[2] ).xyz();
    metalsite_atom_xyz[3] = homodimer_with_matches.residue( metalsite_seqpos[3] ).atom( metalsite_atom_name[3] ).xyz();
    metalsite_atom_xyz[4] = homodimer_with_matches.residue( metalsite_seqpos[4] ).atom( metalsite_atom_name[4] ).xyz();

    for ( core::Size i(1); i <= 5; ++i ) {
      TR << "metalsite: " << metalsite_seqpos[i] << " " << metalsite_atom_name[i] << " " << metalsite_atom_xyz[i] << std::endl;
    }

    /////////////////////////////////////Define parallel axis, apply 180 degree flip//////////////////////////////////////

    axis const parallel_axis = metalsite_atom_xyz[2] - metalsite_atom_xyz[1];
    core::Size const chain_begin = scaffold.conformation().chain_begin(1);
    core::Size const chain_end =  scaffold.conformation().chain_end(1);
    protocols::rigid::RollMoverOP parallel_rollmover = new protocols::rigid::RollMover( chain_begin, chain_end, 180, 180, parallel_axis, metalsite_atom_xyz[5] );
    parallel_rollmover->apply( homodimer_A ); // A and B will diverge after 90 degree or -90 degree flip
    parallel_rollmover->apply( homodimer_B );

    homodimer_A.dump_pdb("squareplanar_homodimer.pdb");


    ////////////////////at this point the homodimers have perfect square-planar metal geometry.  To get tetrahedral, 90 or -90 degree flips are made

    /////////////////////parallel///////////////
    //         p3      p1
    // ChainB      Zn      ChainA
    //         p4      p2
    ////////////////////////////////////////////

    // must update p1 and p2 xyz coordinates after rotation (p3 and p4 remained fixed)

    //this is the axis through zinc, parallel to 3 -> 1 vector (point 1 was changed, needs to be updated)
    axis const axis_90_degree = metalsite_atom_xyz[3] - homodimer_A.residue( metalsite_seqpos[1] ).atom( metalsite_atom_name[1] ).xyz();

    protocols::rigid::RollMoverOP rollmover_90_degree = new protocols::rigid::RollMover( chain_begin, chain_end, 90, 90, axis_90_degree, metalsite_atom_xyz[5] );
    rollmover_90_degree->apply( homodimer_A );
    // also want to flip -90 degrees
    rollmover_90_degree->set_min_max_angles( -90.0, -90.0 ); // change rollmove angle to -90 degrees
    rollmover_90_degree->apply( homodimer_B );

    utility::file::FileName match_name = match.pdb_info()->name(); // assumes match will be appropriately named per the scaffold

    std::string homodimer_A_name = match_name.base() + "_homodimer_A.pdb";
    std::string homodimer_B_name = match_name.base() + "_homodimer_B.pdb";

    homodimer_A.dump_pdb( homodimer_A_name );
    homodimer_B.dump_pdb( homodimer_B_name );

    return;
  }


	virtual
	std::string
	get_name() const { return "zinc1_homodimer_setup"; }


private:
  //utility::vector1< core::Size > metalsite_seqpos_;
  //utility::vector1< point > metalsite_atom_xyz_;
  //utility::vector1< std::string > metalsite_atom_name_;
};

typedef utility::pointer::owning_ptr< zinc1_homodimer_setup > zinc1_homodimer_setupOP;

int main( int argc, char* argv[] )
{
	try {
  basic::options::option.add( scaffold_pdb, "protein monomer for metal-mediated dimerization" ).def("3DE8_A.pdb");
  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go(new zinc1_homodimer_setup);

  TR << "************************d**o**n**e**************************************" << std::endl;

  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }

  return 0;
}

