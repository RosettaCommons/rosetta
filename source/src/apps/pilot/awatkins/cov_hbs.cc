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
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
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
#include <numeric/xyz.functions.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
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

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
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

static basic::Tracer TR("CovDD");

class CovalentPeptidomimeticDockDesign : public Mover {

public:

	//default ctor
	CovalentPeptidomimeticDockDesign(): Mover( "CovalentPeptidomimeticDockDesign" ){}

	//default dtor
	virtual ~CovalentPeptidomimeticDockDesign(){}

	//methods

	virtual void apply( core::pose::Pose & pose );
	void update_hydrogens( core::pose::Pose & pose );
	virtual std::string get_name() const { return "CovalentPeptidomimeticDockDesign"; }

};

typedef utility::pointer::shared_ptr< CovalentPeptidomimeticDockDesign > CovalentPeptidomimeticDockDesignOP;
typedef utility::pointer::shared_ptr< CovalentPeptidomimeticDockDesign const > CovalentPeptidomimeticDockDesignCOP;

int main ( int argc, char* argv[] )
{
	try {
		//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );

		devel::init(argc, argv);
		option[ OptionKeys::in::file::read_pdb_link_records ].value( true );
		option[ OptionKeys::in::constraints_from_link_records ].value( true );
		option[ OptionKeys::out::file::write_pdb_link_records ].value( true );

		CovalentPeptidomimeticDockDesignOP builder( new CovalentPeptidomimeticDockDesign() );
		protocols::jd2::JobDistributor::get_instance()->go( builder );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void
CovalentPeptidomimeticDockDesign::apply(
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

	ScoreFunctionOP scorefxn = get_score_function();



	// We assume that there is already a covalent constraint in place if one is possible
	// So we won't start the simulation by looking for one
	// But, if a subsequent one gets created (for example, because two relevant residues become drawn together)
	// I am open to the idea of locking it down

	// Also, just to help out, we will set the constraint weights in case the user hasn't.
	if ( scorefxn->get_weight( atom_pair_constraint ) == 0.0 ) {
		TR << "Setting atom pair constraint weight to 1 manually" << std::endl;
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}
	if ( scorefxn->get_weight( angle_constraint ) == 0.0 ) {
		scorefxn->set_weight( angle_constraint, 1.0 );
	}
	if ( scorefxn->get_weight( dihedral_constraint ) == 0.0 ) {
		scorefxn->set_weight( dihedral_constraint, 1.0 );
	}

	Size vdp_resi = 0;
	Size cys_resi = 0;

	//protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver( pose, true );

	// Simple minimization first
	kinematics::MoveMapOP mm( new kinematics::MoveMap() );
	PackerTaskOP limited_task( TaskFactory::create_packer_task( pose ) );
	utility::vector1_bool packable( pose.total_residue(), false );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		if ( pose.residue( ii ).type().has_property( "ELECTROPHILE" ) ) {
			vdp_resi = ii;
		}
		if ( pose.residue( ii ).type().has_property( "SIDECHAIN_THIOL" ) ) {
			cys_resi = ii;
		}
		if ( pose.residue( ii ).type().has_property( "ELECTROPHILE" ) ||
				pose.residue( ii ).type().has_property( "SIDECHAIN_THIOL" ) ) {
			mm->set_bb(  ii, true );
			mm->set_chi( ii, true );
			packable[ ii ] = true;
			limited_task->nonconst_residue_task( ii ).or_ex1_sample_level( EX_SIX_QUARTER_STEP_STDDEVS);
			limited_task->nonconst_residue_task( ii ).or_ex2_sample_level( EX_SIX_QUARTER_STEP_STDDEVS);
			limited_task->nonconst_residue_task( ii ).or_ex3_sample_level( EX_SIX_QUARTER_STEP_STDDEVS);
			limited_task->nonconst_residue_task( ii ).or_ex4_sample_level( EX_SIX_QUARTER_STEP_STDDEVS);
		} else {
			mm->set_bb(  ii, false );
			mm->set_chi( ii, false );
		}
	}
	//mm->set_jump( 1, true );

	// legit gonna hard code in ids to check
	kinematics::FoldTree f = pose.fold_tree();
	f.slide_jump( 1, cys_resi, vdp_resi );
	pose.conformation().declare_chemical_bond( cys_resi, "SG", vdp_resi, "CZ" );
	pose.fold_tree( f );
	TR << pose.fold_tree() << std::endl;

	pose::set_reasonable_fold_tree( pose );
	mm->set_branches( cys_resi, true );

	limited_task->restrict_to_residues( packable );
	limited_task->restrict_to_repacking();

	simple_moves::MinMoverOP min( new simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true )  );

	// Use cartesian minimization
	//min->cartesian( true );

	min->apply( pose );
	update_hydrogens( pose );

	//Get the residue set we are drawing from.
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	/******************************************************************************
	Rotamer Trials Setup
	*******************************************************************************/

	using core::pack::task::operation::TaskOperationCOP;
	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *scorefxn, 1 ) );
	// create a task factory and task operations
	//TaskFactoryOP pert_tf(new TaskFactory());
	//pert_tf->push_back( operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	//operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	//pert_rrop->filename("two.res");
	//pert_tf->push_back( pert_rrop );

	// create a rotamer trials mover

	//simple_moves::RotamerTrialsMoverOP pert_rt(new simple_moves::EnergyCutRotamerTrialsMover( scorefxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	protocols::simple_moves::sidechain_moves::SidechainMoverOP SCmover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
	SCmover->set_task( limited_task );
	SCmover->set_prob_uniform( 100 );


	/*********************************************************
	Design Setup
	**********************************************************/

	// create a task factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );
	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	/*
	desn_tf->push_back( pert_rtio ); //not bothering to construct it a second time
	*/


	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( scorefxn );

	// Design loop
	for ( Size ii = 1; ii < 10; ++ii ) {

		desn_pr->apply( pose );

		// Inner perturbation loop- Docking by packing the two key residues
		for ( Size jj = 1; jj < 10; ++jj ) {
			//pert_rt->apply( pose );
			//SCmover->apply( pose );
			//update_hydrogens( pose );
			min->apply( pose );
			//update_hydrogens( pose );
		}
	}

	// Before dumping...
	Real preupdate_score = ( *scorefxn )( pose );
	//update_hydrogens( pose );
	Real final_score = ( *scorefxn )( pose );
	TR << "Final scores: " << preupdate_score << " " << final_score << std::endl;

}

void
CovalentPeptidomimeticDockDesign::update_hydrogens(
	core::pose::Pose & pose
) {

	using namespace core;
	using namespace core::chemical;
	using namespace id;

	// Look at any residues with the property ELECTROPHILE
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ! pose.residue( ii ).type().has_property( "ELECTROPHILE" ) ) {
			continue;
		}

		using numeric::conversions::radians;
		using numeric::conversions::degrees;

		AtomID aidVSG( pose.residue( ii ).atom_index( "VSG" ), ii );
		AtomID aidCZ(  pose.residue( ii ).atom_index( "CZ"  ), ii );
		AtomID aid1HZ( pose.residue( ii ).atom_index( "1HZ" ), ii );
		AtomID aid2HZ( pose.residue( ii ).atom_index( "2HZ" ), ii );
		AtomID aidCE2( pose.residue( ii ).atom_index( "CE2" ), ii );
		AtomID aidCD(  pose.residue( ii ).atom_index( "CD"  ), ii );
		AtomID aidNG(  pose.residue( ii ).atom_index( "NG"  ), ii );

		// Assume conjugation is final connection for now
		Size conn_no = pose.residue( ii ).type().n_possible_residue_connections();
		Size jj = pose.residue( ii ).residue_connection_partner( conn_no );

		//Vector const& VSG_xyz( pose.residue( ii ).xyz( "VSG" ) );
		Vector const& CZ_xyz(  pose.residue( ii ).xyz( "CZ"  ) );
		Vector const& SG_xyz(  pose.residue( jj ).xyz( "SG"  ) );
		Vector const& CE2_xyz( pose.residue( ii ).xyz( "CE2" ) );
		Vector const& CD_xyz(  pose.residue( ii ).xyz( "CD"  ) );
		//Vector const& NG_xyz(  pose.residue( ii ).xyz( "NG"  ) );


		Real VSG_torsion_correction = numeric::dihedral_degrees( SG_xyz, CZ_xyz, CE2_xyz, CD_xyz ) - degrees( pose.conformation().torsion_angle( aidVSG, aidCZ, aidCE2, aidCD ) );
		//if ( VSG_torsion_correction > 0.01 || VSG_torsion_correction < 0.01) {
		TR << "Initial 1HZ-CZ-CE2-CD: " << degrees(pose.conformation().torsion_angle(aid1HZ,aidCZ,aidCE2,aidCD)) << std::endl;
		TR << "Initial VSG-CZ-CE2-CD: "<< degrees(pose.conformation().torsion_angle(aidVSG,aidCZ,aidCE2,aidCD)) <<std::endl;
		TR << "Initial SG-CZ-CE2-CD: "<< numeric::dihedral_degrees(SG_xyz,CZ_xyz,CE2_xyz,CD_xyz) <<std::endl;
		TR << " Out of sync by " << VSG_torsion_correction << " degrees " << std::endl;
		Real torsion_1HZ = degrees(pose.conformation().torsion_angle(aid1HZ,aidCZ,aidCE2,aidCD)) + VSG_torsion_correction;
		pose.conformation().set_torsion_angle(aid1HZ,aidCZ,aidCE2,aidCD,radians(torsion_1HZ));
		//pose.dump_pdb( "rosetta_out_oop_pre_mvH.pdb" );
		TR << "Final   VSG-CZ-CE2-CD: "<< degrees(pose.conformation().torsion_angle(aidVSG,aidCZ,aidCE2,aidCD)) <<std::endl;
		TR << "Final    SG-CZ-CE2-CD: "<< numeric::dihedral_degrees(SG_xyz,CZ_xyz,CE2_xyz,CD_xyz) <<std::endl;
		TR << "Final   1HZ-CZ-CE2-CD: "<< degrees(pose.conformation().torsion_angle(aid1HZ,aidCZ,aidCE2,aidCD)) <<std::endl;
		//}

	}
}

