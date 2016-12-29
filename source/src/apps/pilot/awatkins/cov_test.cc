// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   cov_hbs.cc
/// @brief  Sidechain conjugation to acryl amides
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Numeric Headers
#include <numeric/conversions.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;

using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main ( int argc, char* argv[] )
{
	try {
		//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );

		devel::init(argc, argv);

		using namespace core;
		using namespace utility;
		using namespace scoring;
		using namespace pose;
		using namespace core::chemical;
		using namespace conformation;
		using namespace func;
		using namespace constraints;

		using namespace core::id;
		using namespace core::pack;
		using namespace core::pack::task;


		ScoreFunctionOP scorefxn = get_score_function();

		//Get the residue set we are drawing from.
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		Pose pose;

		ResidueType const & cys = residue_set_cap->name_map( "CYS:S-conjugated" );
		ResidueType const & vdp = residue_set_cap->name_map( "VDP:sidechain_electrophile_conjugated" );
		Residue res_cys( cys, true );
		Residue res_vdp( vdp, true );

		std::cout << pose;
		//pose.pdb_info()->attach_to( pose.conformation() );
		//pose.pdb_info()->show(std::cout);
		pose.append_residue_by_jump( res_cys, 1 );

		std::cout << pose;
		//pose.pdb_info()->show(std::cout);

		std::cout << "Cys residue type has " << cys.n_possible_residue_connections() << " and vdp has " << vdp.n_possible_residue_connections() << std::endl;
		runtime_assert( cys.n_possible_residue_connections() == 3 &&
			cys.lower_connect_id() != 3 &&
			cys.upper_connect_id() != 3 );

		runtime_assert( vdp.n_possible_residue_connections() == 3 &&
			vdp.lower_connect_id() != 3 );

		// hardcoded connection numbers
		//pose.append_residue_by_bond( res_vdp, true, 3, 1, 3 );
		pose.append_residue_by_jump( res_vdp, 2 ); //true, 3, 1, 3 );

		kinematics::FoldTree f = pose.fold_tree();
		f.clear();
		f.add_edge( 1, 1, -1 );
		f.add_edge( 1, 2, -2 );


		//pose.pdb_info()->set_resinfo( 1, 'A', 1, ' ' );
		//pose.pdb_info()->set_resinfo( 2, 'B', 1, ' ' );

		std::cout << pose;
		//pose.pdb_info()->show(std::cout);

		core::id::AtomID const C( cys.atom_index( "C" ), 1 );
		core::id::AtomID const CA( cys.atom_index( "CA" ), 1 );
		core::id::AtomID const CB( cys.atom_index( "CB" ), 1 );
		core::id::AtomID const SG( cys.atom_index( "SG" ), 1 );
		core::id::AtomID const CZ( vdp.atom_index( "CZ"  ), 2 );
		core::id::AtomID const CE2( vdp.atom_index( "CE2" ), 2 );
		core::id::AtomID const CD( vdp.atom_index( "CD"  ), 2 );

		pose.conformation().set_torsion_angle( C, CA, CB, SG, numeric::conversions::radians(106.5) );
		pose.conformation().set_torsion_angle( CA, CB, SG, CZ, numeric::conversions::radians(-60.0) );
		pose.conformation().set_torsion_angle( CB, SG, CZ, CE2, numeric::conversions::radians(180.0) );
		pose.conformation().set_torsion_angle( SG, CZ, CE2, CD, numeric::conversions::radians(100.5) );
		//pose.conformation().set_bond_angle( CB, SG, CZ, numeric::conversions::radians(90) );
		pose.conformation().set_bond_length( SG, CZ, 1.8 );

		//pose.pdb_info()->rebuild_pdb2pose();
		std::cout << "about to dump" << std::endl;
		pose.dump_pdb("out.pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
