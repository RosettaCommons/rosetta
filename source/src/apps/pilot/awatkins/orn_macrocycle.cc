// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyz.functions.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace core::chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("OrnMacrocycle");

void
add_orn_cst(
	core::pose::Pose & pose,
	Size res1,
	Size res2 ) {

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( res1 ).atom_index("CA"), res1 ),
		*new AtomID( pose.residue( res1 ).atom_index("C" ), res1 ),
		*new AtomID( pose.residue( res2 ).atom_index("NE"), res2 ),
		*new AtomID( pose.residue( res2 ).atom_index("CD"), res2 ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( res1 ).atom_index( "O" ), res1 ),
		*new AtomID( pose.residue( res1 ).atom_index( "C" ), res1 ),
		*new AtomID( pose.residue( res2 ).atom_index( "NE"), res2 ),
		*new AtomID( pose.residue( res2 ).atom_index("1HE"), res2 ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( res1 ).atom_index( "O" ), res1 ),
		*new AtomID( pose.residue( res1 ).atom_index( "C" ), res1 ),
		*new AtomID( pose.residue( res1 ).atom_index( "CA"), res1 ),
		*new AtomID( pose.residue( res2 ).atom_index( "NE"), res2 ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );

	pose.add_constraint(
		AngleConstraintOP( new AngleConstraint(
		*new AtomID( pose.residue( res1 ).atom_index( "O" ), res1 ),
		*new AtomID( pose.residue( res1 ).atom_index( "C" ), res1 ),
		*new AtomID( pose.residue( res2 ).atom_index( "NE"), res2 ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi()*2.0/3.0, 0.02 ) ) ) ) );

	pose.add_constraint(
		AtomPairConstraintOP( new AtomPairConstraint(
		*new AtomID( pose.residue( res1 ).atom_index("C" ), res1 ),
		*new AtomID( pose.residue( res2 ).atom_index("NE"), res2 ),
		//TopOutFuncOP( new TopOutFunc( 100, 1.33, 5 ) ) ) ) );
		HarmonicFuncOP( new HarmonicFunc( 1.33, .02 ) ) ) ) );
}


int
main( int argc, char* argv[] )
{
	try {

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		core::pose::Pose pose;

		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

		core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		ResidueType const & ala_type = residue_set_cap->name_map( "ALA" );
		ResidueType const & ala_CT_type = residue_set_cap->name_map( "ALA:C-term_conjugated" );
		ResidueType const & orn_type = residue_set_cap->name_map( "C40:NtermProteinFull:N-conjugated" );
		Residue ala( ala_type, true );
		Residue alaCT( ala_CT_type, true );
		Residue orn( orn_type, true );

		pose.conformation().append_residue_by_jump( orn, 1 );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( alaCT, true );

		pose.conformation().append_residue_by_jump( orn, 6 );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( ala, true );
		pose.conformation().append_residue_by_bond( alaCT, true );

		// Hardcoded connection IDs; never do this! This just illustrates the use case.
		// Note that the NE connections are still #2 because they're ntermproteinfull.

		ResidueOP new_orn = ResidueFactory::create_residue( orn_type, pose.residue( 1 ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( 1 ), *new_orn, pose.conformation() );
		new_orn->residue_connection_partner( 2, 12, 2 );
		pose.conformation().replace_residue( 1, *new_orn, false );
		ResidueOP new_ala = ResidueFactory::create_residue( ala_CT_type, pose.residue( 12 ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( 12 ), *new_ala, pose.conformation() );
		new_ala->residue_connection_partner( 2,  1, 2 );
		pose.conformation().replace_residue( 12, *new_ala, false );
		pose.conformation().declare_chemical_bond( 12, "C", 1, "NE" );


		ResidueOP new_orn2 = ResidueFactory::create_residue( orn_type, pose.residue( 7 ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( 7 ), *new_orn2, pose.conformation() );
		new_orn2->residue_connection_partner( 2,  6, 2 );
		pose.conformation().replace_residue( 7, *new_orn2, false );
		ResidueOP new_ala2 = ResidueFactory::create_residue( ala_CT_type, pose.residue(  6 ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue(  6 ), *new_ala2, pose.conformation() );
		new_ala2->residue_connection_partner( 2,  7, 2 );
		pose.conformation().replace_residue(  6, *new_ala2, false );
		pose.conformation().declare_chemical_bond(  6, "C", 7, "NE" );

		Size nres = pose.size();

		add_orn_cst( pose, nres, 1 );
		add_orn_cst( pose, 6, 7 );

		for ( Size ii = 1; ii <= nres; ++ii ) {
			pose.set_phi( ii, -150 );
			pose.set_psi( ii, 150 );
			pose.set_omega( ii, 180.0 );
		}

		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( 1, true );

		// create minimization mover
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

		/********************************************************\
		Initial MC
		\********************************************************/
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );

		for ( Size i = 1; i <= pose.size(); ++i ) {
			pert_pep_mm->set_bb( i );
		}

		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );
		simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_shear->angle_max( 'H', 2.0 );
		pert_pep_shear->angle_max( 'L', 2.0 );
		pert_pep_shear->angle_max( 'E', 2.0 );

		moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
		pert_pep_random->add_mover( pert_pep_small, 1 );
		pert_pep_random->add_mover( pert_pep_shear, 1 );
		moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 10 ) );

		pert_sequence->add_mover( pert_pep_repeat );
		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		// core::Size hbs_position = 1;
		Real wt = 0.01;
		for ( Size j = 1; j <= 100; ++j ) {

			score_fxn->set_weight( core::scoring::atom_pair_constraint, wt );
			score_fxn->set_weight( core::scoring::angle_constraint, wt );
			score_fxn->set_weight( core::scoring::dihedral_constraint, wt );

			minM->apply( pose );

			pert_trial->apply( pose );

			pert_mc->recover_low( pose );
			//TR<< "pre mc->boltzmann" << std::endl;
			//pert_mc->show_state();
			pert_mc->boltzmann( pose );
			//TR<< "post mc->boltzmann" << std::endl;
			//pert_mc->show_state();

			wt += 0.01;
			TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;
		}
		pert_mc->recover_low( pose );
		minM->apply( pose );
		TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;
		pose.dump_pdb( "out.pdb" );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
