// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   triazolamer_creator
/// @brief  Make a silly triazolamer
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
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

void
perturb_and_rescore( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sfxn ) {

	using namespace core;
	using namespace id;

	Real apc_wt = sfxn->get_weight( scoring::atom_pair_constraint );
	Real ang_wt = sfxn->get_weight( scoring::angle_constraint );
	Real dih_wt = sfxn->get_weight( scoring::dihedral_constraint );

	//sfxn->set_weight( scoring::atom_pair_constraint, 0 );
	//sfxn->set_weight( scoring::angle_constraint, 0 );
	//sfxn->set_weight( scoring::dihedral_constraint, 0 );

	for ( Size j = 1; j <= 10; ++j ) {

		Real initial_score = ( *sfxn )( pose );
		std::cout << "Perturbation " << j << " / 10: initial score " << initial_score << std::endl;

		pose::Pose copy_pose = pose;

		for ( Size i = 1; i <= pose.n_residue()-2; ++i ) {
			Real phi_like = pose.conformation().torsion_angle(
				// replacing w/ equiv that is in atom tree
				AtomID( pose.residue( i+1 ).atom_index( /*"CT2"*/"NT2" ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( "CA"  ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( /*"CT1"*/"C" ), i+1 ) );

			pose.conformation().set_torsion_angle(
				AtomID( pose.residue( i+1 ).atom_index( /*"CT2"*/"NT2" ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( "CA"  ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( /*"CT1"*/"C" ), i+1 ),
				phi_like + numeric::random::rg().gaussian() * 2.0 );
		}
		/*Size i = pose.n_residue()-1;
		Real phi_like = pose.conformation().torsion_angle(
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "CH3" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "CT1" ), i+1 ) );

		pose.conformation().set_torsion_angle(
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "CH3" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "CT1" ), i+1 ),
		phi_like + numeric::random::rg().gaussian() * 20.0 );
		*/

		for ( Size i = 2; i <= pose.n_residue()-1; ++i ) {

			Real psi_like = pose.conformation().torsion_angle(
				AtomID( pose.residue( i   ).atom_index( "NT1" ), i   ),
				AtomID( pose.residue( i   ).atom_index( "CA"  ), i   ),
				AtomID( pose.residue( i   ).atom_index( /*"CT1"*/"C" ), i   ),
				AtomID( pose.residue( i+1 ).atom_index( /*"NT3"*/"N" ), i+1 ) );

			pose.conformation().set_torsion_angle(
				AtomID( pose.residue( i   ).atom_index( "NT1" ), i   ),
				AtomID( pose.residue( i   ).atom_index( "CA"  ), i   ),
				AtomID( pose.residue( i   ).atom_index( /*"CT1"*/"C" ), i   ),
				AtomID( pose.residue( i+1 ).atom_index( /*"NT3"*/"N" ), i+1 ),
				psi_like + numeric::random::rg().gaussian() * 2.0 );
		}
		/*i = 1;
		Real psi_like = pose.conformation().torsion_angle(
		AtomID( pose.residue( i   ).atom_index( "NT1" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CH3" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ) );

		pose.conformation().set_torsion_angle(
		AtomID( pose.residue( i   ).atom_index( "NT1" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CH3" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ),
		psi_like + numeric::random::rg().gaussian() * 20.0 );
		*/

		core::Real new_score = ( *sfxn )( pose );
		std::cout << "Perturbation " << j << " / 10: new score " << new_score << std::endl;

		// Reset dihedrals if rejected
		// temp 1.0, whatever
		if ( new_score < initial_score ) {
			std::cout << "Perturbation " << j << " / 10: accepted" << std::endl;
		} else if ( numeric::random::rg().uniform() > exp(-1.0 * ( new_score - initial_score ) ) ) {
			// reject
			std::cout << "Perturbation " << j << " / 10: rejected" << std::endl;
			pose = copy_pose;
		} else {
			std::cout << "Perturbation " << j << " / 10: accepted thermally" << std::endl;
		}
	}

	sfxn->set_weight( scoring::atom_pair_constraint, apc_wt );
	sfxn->set_weight( scoring::angle_constraint, ang_wt );
	sfxn->set_weight( scoring::dihedral_constraint, dih_wt );

}

class TriazoleCreator : public Mover {

public:

	//default ctor
	TriazoleCreator(): Mover("A3BPeptideBuilder"){}

	//default dtor
	virtual ~TriazoleCreator(){}

	//methods

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "TriazoleCreator"; }

};

typedef utility::pointer::shared_ptr< TriazoleCreator > TriazoleCreatorOP;
typedef utility::pointer::shared_ptr< TriazoleCreator const > TriazoleCreatorCOP;

int main ( int argc, char* argv[] )
{
	try {
		//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );

		devel::init(argc, argv);

		TriazoleCreatorOP builder( new TriazoleCreator() );
		protocols::jd2::JobDistributor::get_instance()->go( builder );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void
TriazoleCreator::apply(
	core::pose::Pose & pose
) {

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

	//first, load the file of residue types to get min energies for.

	//now do initialization stuff.
	TaskFactoryOP task_factory = TaskFactoryOP( new TaskFactory );
	task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
	//need these to keep pack_rotamers from redesigning the residue.
	operation::RestrictToRepackingOP rtrop = operation::RestrictToRepackingOP( new operation::RestrictToRepacking );
	task_factory->push_back( rtrop );

	ScoreFunctionOP scorefxn = get_score_function();




	//Get the residue set we are drawing from.
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose.clear();

	ResidueType const & restype_first = residue_set_cap->name_map( "BZO:triazolamerC" );
	ResidueType const & internal_ALA = residue_set_cap->name_map( "ALA:triazolamerN:triazolamerC" );
	ResidueType const & internal_LYS = residue_set_cap->name_map( "LYS:triazolamerN:triazolamerC" );
	ResidueType const & internal_GLU = residue_set_cap->name_map( "GLU:triazolamerN:triazolamerC" );
	ResidueType const & restype_last = residue_set_cap->name_map( "BZA:triazolamerN" );
	Residue res_first( restype_first, true );
	Residue res_int_ALA( internal_ALA, true );
	Residue res_int_LYS( internal_LYS, true );
	Residue res_int_GLU( internal_GLU, true );
	Residue res_last( restype_last, true );
	pose.append_residue_by_jump( res_first, 1 );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	/*pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	*/pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_last, true );

	//HarmonicFuncOP hf( new HarmonicFunc( 1.347, 0.1 ) );
	//HarmonicFuncOP zf( new HarmonicFunc( 0.000, 0.01 ) );
	//CircularHarmonicFuncOP chf( new CircularHarmonicFunc( 0, 0.01 ) );
	for ( Size i = 1; i <= pose.total_residue() - 1; ++i ) {

		pose.conformation().declare_chemical_bond( i, "CT1", i+1, "NT3" );
		core::pose::ncbb::add_triazole_constraint( pose, i );
	}

	// create movemap for peptide
	kinematics::MoveMapOP pert_mm( new kinematics::MoveMap() );
	pert_mm->set_bb( 1, true );
	pert_mm->set_chi( 1, true );
	pert_mm->set_bb( pose.n_residue(), true );
	pert_mm->set_chi( pose.n_residue(), true );

	for ( Size i = 1+1; i <= pose.n_residue()-1; ++i ) {

		pert_mm->set_bb( i, true );

		pose.set_phi(   i, -135 );
		pose.set_psi(   i,  135 );
		pose.set_omega( i,  180 );
	}

	pose.conformation().detect_bonds();
	//pose.conformation().detect_pseudobonds();
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		pose.conformation().update_polymeric_connection(i);
	}

	protocols::simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( pert_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
	desn_min->cartesian( true );

	std::cout << "Dump initial" << std::endl;
	pose.dump_scored_pdb( "B3A_initial.pdb", *scorefxn);
	Real wt = 0.001;
	for ( Size i = 1; i <= 1000; ++i ) {

		//perturb_and_rescore( pose, scorefxn );

		scorefxn->set_weight( core::scoring::atom_pair_constraint, wt*2 );
		scorefxn->set_weight( core::scoring::angle_constraint, wt );
		scorefxn->set_weight( core::scoring::dihedral_constraint, wt );
		desn_min->apply( pose );
		wt += 0.001;
	}

	for ( Size i = 1; i <= 500; ++i ) {
		wt -= 0.001;
		scorefxn->set_weight( core::scoring::atom_pair_constraint, wt*2 );
		scorefxn->set_weight( core::scoring::angle_constraint, wt );
		scorefxn->set_weight( core::scoring::dihedral_constraint, wt );
		desn_min->apply( pose );
	}

	std::cout << "Dump initial min" << std::endl;
	pose.dump_scored_pdb( "B3A_min.pdb", *scorefxn);
}
