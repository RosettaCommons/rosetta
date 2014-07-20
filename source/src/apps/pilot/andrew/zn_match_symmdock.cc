// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <devel/init.hh>
#include <devel/znhash/ZnHash.hh>
#include <devel/znhash/SymmZnMoversAndTaskOps.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/Remarks.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>

#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RampingMover.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesMovers.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>

// option keys
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <sstream>


OPT_1GRP_KEY( String, zn_match_symmdock, zn_tworesidue_matches_list )
OPT_1GRP_KEY( String, zn_match_symmdock, reference_pdb )
OPT_1GRP_KEY( String, zn_match_symmdock, two_residue_match_constraint_file )
OPT_1GRP_KEY( String, zn_match_symmdock, four_residue_match_constraint_file )
OPT_1GRP_KEY( Real, zn_match_symmdock, metalhash_constraint_weight )
OPT_1GRP_KEY( Real, zn_match_symmdock, clash_weight )
OPT_1GRP_KEY( Real, zn_match_symmdock, znreach )
OPT_1GRP_KEY( Real, zn_match_symmdock, orbital_dist )
OPT_1GRP_KEY( Real, zn_match_symmdock, orbital_reach )
OPT_1GRP_KEY( Real, zn_match_symmdock, znwelldepth )
OPT_1GRP_KEY( Real, zn_match_symmdock, final_cstscore_limit )
OPT_1GRP_KEY( Boolean, zn_match_symmdock, require_3H )
OPT_1GRP_KEY( Boolean, zn_match_symmdock, preserve_input_virtual_atoms )

static basic::Tracer TR("apps.pilot.andrew.zn_match_symmdock");


void initialize_initalizeZNcst(
	devel::znhash::InitializeZNCoordinationConstraintMoverOP init_zn
)
{
	// Required options
	if ( ! basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::reference_pdb].user() ) {
		utility_exit_with_message("zn_match_symmdock::reference_pdb was not given on the command line" );
	}
	if ( ! basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::two_residue_match_constraint_file].user() ) {
		utility_exit_with_message( "Must provide a file for -zn_match_symmdock::two_residue_match_constraint_file on the comamnd line" );
	}
	if ( ! basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::four_residue_match_constraint_file].user() ) {
		utility_exit_with_message( "Must provide a file for -zn_match_symmdock::four_residue_match_constraint_file on the comamnd line" );
	}
	if ( ! basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::zn_tworesidue_matches_list].user() ) {
		utility_exit_with_message("zn_match_symmdock::zn_tworesidue_matches_list was not given on the command line" );
	}



	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::znreach].user() ) {
		init_zn->set_zn_reach( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::znreach] );
	}
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::orbital_dist].user() ) {
		init_zn->set_orbital_dist( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::orbital_dist] );
	}
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::orbital_reach].user() ) {
		init_zn->set_orbital_reach( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::orbital_reach] );
	}
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::znwelldepth].user() ) {
		init_zn->set_zn_well_depth( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::znwelldepth] );
	}
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::require_3H].user() ) {
		init_zn->require_3H( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::require_3H] );
	}
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::preserve_input_virtual_atoms].user() ) {
		init_zn->set_idealize_input_virtual_atoms(
			! basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::preserve_input_virtual_atoms]
		);
	}
	init_zn->set_matcher_constraint_file_name(
			basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::two_residue_match_constraint_file]
		);

	init_zn->set_reference_pdb( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::reference_pdb] );

	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::clash_weight ].user() ) {
		init_zn->set_clash_weight( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::clash_weight ] );
	}

	/// now open the match-cst-pair list file and read off all the pairs.  Insert them into the znscore object
	init_zn->set_match_pdb_listfilename( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::zn_tworesidue_matches_list]  );


}


int main( int argc, char * argv [] )
{
	try {

	//using namespace core;
	//using namespace core::chemical;
	//using namespace core::id;
	//usng namespace core::io::pdb;
	//using namespace core::pose;

	//using namespace basic::options;
	//using namespace basic::options::OptionKeys::measure;

	NEW_OPT( zn_match_symmdock::zn_tworesidue_matches_list, "REQUIRED: File listing the output-match .pdb files", "");
	NEW_OPT( zn_match_symmdock::reference_pdb, "REQUIRED: Reference PDB file that was used to generate the matches", "");
	NEW_OPT( zn_match_symmdock::two_residue_match_constraint_file, "REQUIRED: Enzyme Design / Matcher / Geometric Constraint file that was used generate the matches.  This should have two residues coordinating zinc", "");
	NEW_OPT( zn_match_symmdock::four_residue_match_constraint_file, "REQUIRED: Enzyme Design / Matcher / Geometric Constraint file similar to the one that was used to generate the matches except that it contains information for how four residues coordinating zinc simultaneously", "");
	NEW_OPT( zn_match_symmdock::metalhash_constraint_weight, "Strength of the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::clash_weight, "Strength of the clash weight in the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::znreach, "Strength of the clash weight in the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::orbital_dist, "Strength of the clash weight in the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::orbital_reach, "Strength of the clash weight in the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::znwelldepth, "Strength of the clash weight in the constraint weight", 1.0 );
	NEW_OPT( zn_match_symmdock::require_3H, "Strength of the clash weight in the constraint weight", false );
	NEW_OPT( zn_match_symmdock::preserve_input_virtual_atoms, "Use the coordinates for zinc virtual atoms in the match PDB instead of idealizing those coordinates", false );
	NEW_OPT( zn_match_symmdock::final_cstscore_limit, "Upper limit for the constraint score after the final refinemen stage, above which a trajectory is said to have failed", 25.0 );



	devel::init( argc, argv );

	/// stolen from SymDockProtocol.cc
	using namespace protocols::simple_moves::symmetry;
	using namespace protocols::symmetric_docking;

	/// Sequence mover will hold three movers for this protocol.
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;

	/// 1. prep for symmetry
	SetupForSymmetryMoverOP setup_mover = new SetupForSymmetryMover;
	seq_mover->add_mover( setup_mover );

	/// 2. initialize the Zn coordination constraint
	devel::znhash::InitializeZNCoordinationConstraintMoverOP init_zn =
		new devel::znhash::InitializeZNCoordinationConstraintMover;
	initialize_initalizeZNcst( init_zn );
	seq_mover->add_mover( init_zn );

	/// 3. Perform low-resolution docking.
	protocols::symmetric_docking::SymDockProtocolOP dock_mover = new protocols::symmetric_docking::SymDockProtocol;
	core::scoring::ScoreFunctionOP docking_sfxn;
	if ( basic::options::option[ basic::options::OptionKeys::docking::low_patch ].user() ) {
		docking_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen",
			basic::options::option[ basic::options::OptionKeys::docking::low_patch ] );
	} else {
		docking_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
	}
	docking_sfxn->set_weight( core::scoring::metalhash_constraint,
		basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::metalhash_constraint_weight] );
	dock_mover->set_lowres_scorefxn( docking_sfxn );
	dock_mover->set_fullatom( false );
	seq_mover->add_mover( dock_mover );


	// 4. Go to full atom
	seq_mover->add_mover( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ));

	// 5. ZNCoordinationConstraintReporterMover to report on which residues are zn coordinated
	seq_mover->add_mover( new devel::znhash::ZNCoordinationConstraintReporterMover( init_zn ) );


	// 6. Place best match on the pose and restore the native sidechains
	devel::znhash::ZNCoordinationConstraintPlacerMoverOP znplacer = new devel::znhash::ZNCoordinationConstraintPlacerMover( init_zn );
	if ( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::final_cstscore_limit ].user() ) {
		znplacer->set_constraint_energy_cutoff( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::final_cstscore_limit ] );
	}
	znplacer->set_four_residue_cst_fname( basic::options::option[ basic::options::OptionKeys::zn_match_symmdock::four_residue_match_constraint_file] );

	seq_mover->add_mover( znplacer );

	// GO!
	protocols::jd2::JobDistributor::get_instance()->go( seq_mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

