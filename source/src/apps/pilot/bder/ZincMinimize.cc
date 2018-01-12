// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rjha/ZincMinimize.cc
/// @brief This application is designed to help filter through the matches from RosettaMatch.  It calculates some of the geometry about the metal atom (currently the metal-ligand atom distances and the tetrahedral coordinations six ligand-metal-ligand angles (all of which should be 109.5).  It also outputs a sum of squares for those values.
/// @author Bryan Der

// Unit Headers
#include <devel/metal_interface/FindClosestAtom.hh>
#include <devel/metal_interface/ZincSiteFinder.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>
#include <devel/metal_interface/AddZincSiteConstraints.hh>
// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/moves/Mover.hh>
#include <core/id/AtomID.hh>
// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.io.hh>
// Utility Headers
#include <utility/file/FileName.hh>
#include <devel/init.hh>
//#include <core/types.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <core/kinematics/MoveMap.hh>
// C++ headers
#include <iostream>

#include <core/kinematics/FoldTree.hh>
#include <core/scoring/Energies.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/scoring/constraints/ConstraintSet.hh>



//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.bder.ZincMinimize");

typedef numeric::xyzVector<Real> point;
typedef point axis;

using namespace core;


//local options
namespace local {
basic::options::BooleanOptionKey const zn_minimize("zn_minimize");
basic::options::BooleanOptionKey const bb_min("bb_min");
basic::options::BooleanOptionKey const zn_min_pack("zn_min_pack");
basic::options::IntegerOptionKey const zn_min_pack_cycles("zn_min_pack_cycles");
basic::options::RealOptionKey const zn_constraint_weight("zn_constraint_weight");
}//local



/// @brief
class ZincMinimize : public protocols::moves::Mover {
public:
	ZincMinimize()
	{
	}
	virtual ~ZincMinimize(){};



	virtual
	void
	add_zinc_constraints( core::pose::Pose & pose ) {
		TR<< "FoldTree: \n" << pose.fold_tree() << std::endl;
		pdbname_ = pose.pdb_info()->name();
		pose_filename_ = pdbname_; // not sure if this is correct

		for ( Size i(1); i <= pose.size(); ++i ) {
			std::string name3 = pose.residue(i).name3();
			if ( name3 == "ZNX" || name3 == " ZN" || name3 == "ZN " || name3 == "ZN" || name3 == "HIZ" ) {

				TR << "ADDING CONSTRAINTS FOR ZINC RESIDUE: " << i << std::endl;
				//only partly set up for multiple zincs: all constraints should be added, but msr_ is only the most recent msr
				devel::metal_interface::ZincSiteFinderOP find_zinc = new devel::metal_interface::ZincSiteFinder();
				find_zinc->set_expecting_n_ligands(2);
				utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr( find_zinc->find_zinc_site( pose ) );
				msr_ = msr;
				for ( Size i(1); i <= msr_.size(); i++ ) {
					TR << i << " " << msr_[i]->get_seqpos() << " " << msr_[i]->get_ligand_atom_name() << " " << msr_[i]->get_ligand_atom_xyz() << std::endl;
				}

				//AddZincSiteConstraints
				devel::metal_interface::AddZincSiteConstraintsOP zinc_constraints = new devel::metal_interface::AddZincSiteConstraints( msr_ );
				zinc_constraints->add_constraints( pose );
				zinc_constraints->evaluate_constraints( pose ); // prints stuff
			}
		}

		TR << "This is how many constraints the pose has: " << pose.constraint_set()->get_all_constraints().size() << std::endl;

		return;
	}

	virtual
	void
	set_scorefunction() {

		//at the end, we want to score the pose without constraints
		Real zn_cst_weight( basic::options::option[local::zn_constraint_weight] );
		using namespace core::scoring;
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
		ScoreFunctionOP scorefxn_zn = ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
		//ScoreFunctionOP scorefxn_zn = new ScoreFunction;
		scorefxn_ = scorefxn;
		scorefxn_zn_ = scorefxn_zn;
		scorefxn_zn_->set_weight( atom_pair_constraint, zn_cst_weight );
		scorefxn_zn_->set_weight( angle_constraint, zn_cst_weight );
		scorefxn_zn_->set_weight( dihedral_constraint, zn_cst_weight );

		return;
	}


	virtual
	void
	set_min_mover() {

		//movemap
		kinematics::MoveMapOP movemap = new kinematics::MoveMap();
		movemap->set_chi(true);
		movemap->set_bb(false); //make this an option in the future
		movemap->set_jump(true); //will move this HIZ residue
		movemap_ = movemap;

		//minmover
		protocols::minimization_packing::MinMoverOP min_mover = new protocols::minimization_packing::MinMover( movemap_, scorefxn_zn_, "dfpmin_armijo", 0.01, true );
		min_mover_ = min_mover;

		return;
	}




	virtual
	void
	zn_min_pack( core::pose::Pose & pose ) {

		//read resfile
		using namespace core::pack::task;
		using namespace basic::options;
		TaskFactoryOP task_factory = new TaskFactory();
		task_factory->push_back(new operation::InitializeFromCommandline()); //ex1, ex1, minimize sidechains, use_input_sc

		// first, read in user resfile.  Intended to only contain NATAA or NATRO to specify additional residues to not mutation  MUST have ALLAA as default!!
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			task_factory->push_back( new operation::ReadResfile );
			TR << "Reading resfile from user input." << std::endl;
		}

		protocols::minimization_packing::PackRotamersMoverOP packrot_mover = new protocols::minimization_packing::PackRotamersMover;
		packrot_mover->score_function( scorefxn_zn_ );
		packrot_mover->task_factory( task_factory );

		for ( Size i(1); i <= basic::options::option[local::zn_min_pack_cycles]; ++i ) {
			TR << "Iteration " << i << std::endl;
			TR << "Score before pack: " << scorefxn_zn_->score(pose) << std::endl << std::endl;
			packrot_mover->apply( pose );
			TR << "Score before min : " << scorefxn_zn_->score(pose) << std::endl << std::endl;
			min_mover_->apply( pose );
			TR << "Score after  min : " << scorefxn_zn_->score(pose) << std::endl << std::endl;
		}

		//at the end, we want to score the pose without constraints
		scorefxn_->score( pose );

		return;
	}



	virtual
	void
	apply( core::pose::Pose & pose ){

		add_zinc_constraints( pose );

		set_scorefunction();

		set_min_mover();


		if ( basic::options::option[local::zn_minimize] ) {
			min_mover_->apply( pose );
			if ( basic::options::option[local::bb_min] ) {
				movemap_->set_bb(true);
				movemap_->set_jump(true);
				min_mover_->apply( pose );
			}

			std::string name_min = pose_filename_.base() + "_min.pdb";
			pose.dump_scored_pdb(name_min, *scorefxn_);
		}


		if ( basic::options::option[local::zn_min_pack] ) {
			zn_min_pack( pose );
		}

		return;
	}


	virtual
	std::string
	get_name() const { return "ZincMinimize"; }


private:
	std::string pdbname_;
	utility::file::FileName pose_filename_;
	utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_zn_;
	kinematics::MoveMapOP movemap_;
	protocols::minimization_packing::MinMoverOP min_mover_;

};

typedef utility::pointer::owning_ptr< ZincMinimize > ZincMinimizeOP;

int main( int argc, char* argv[] )
{
	using basic::options::option;
	option.add( local::zn_minimize, "minimize structure containing zinc site").def(true);
	option.add( local::bb_min, "bb_min").def(false);
	option.add( local::zn_min_pack, "iterate between minimization and packing for structure containing zinc site").def(false);
	option.add( local::zn_min_pack_cycles, "min pack iterations").def(5);
	option.add( local::zn_constraint_weight, "zinc constraint weight").def(1.0);

	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new ZincMinimize);

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;
}

