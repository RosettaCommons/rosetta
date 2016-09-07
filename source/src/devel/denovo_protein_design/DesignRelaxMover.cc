// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/DenovoProteinDesign/DesignRelaxMover.cc
/// @brief DesignRelaxMover methods implemented
/// @author Grant Murphy g.s.murphy@gmail.com


// Unit Headers
#include <devel/denovo_protein_design/DesignRelaxMover.hh>

// Package Headers


#include <basic/options/option.hh>
#include <basic/options/keys/DenovoProteinDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loopfinder.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <utility>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// C++ Headers
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core::scoring;
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR( "devel.DenovoProteinDesign.DesignRelaxMover" );

namespace devel {
namespace denovo_protein_design {

using namespace core::pack::task;

/// @details
void DesignRelaxMover::apply( core::pose::Pose & pose )
{
	using core::pack::task::operation::TaskOperationCOP;

	bool debug( false );
	if ( basic::options::option[ basic::options::OptionKeys::run::debug ].user() ) {
		debug = true;
	}

	// we don't know how the pose coming in was created but it needs to be a full atom pose
	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set(pose, core::chemical::FA_STANDARD);
		TR << "switching pose to be full atom" << std::endl;
	}

	// design step - taskfactory logic determines which positions are designable, packable and fixed
	// this can be done outside the mover or on the commandline or from resfile
	// the designtaskfactory_ comes from the ctor
	// core::pack::task::TaskFactoryOP designtaskfactory_ = new core::pack::task::TaskFactory;
	// designrelaxmover works with resfile and commandline operations
	designtaskfactory_->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	// check for resfile - if present use it
	if ( basic::options::option [ basic::options::OptionKeys::packing::resfile ].user()  ) {
		designtaskfactory_->push_back( TaskOperationCOP( new operation::ReadResfile ) );
	}


	// make decision about what loops may be rebuilt
	bool loop_file_is_present = option[ OptionKeys::loops::loop_file].user();
	// read in loops from a file if given
	protocols::loops::LoopsOP loops( new protocols::loops::Loops( loop_file_is_present ) );

	// we need this for the rama term
	basic::options::option[ basic::options::OptionKeys::loops::nonpivot_torsion_sampling ].value(true);
	if ( !loop_file_is_present ) {
		protocols::loops::loopfinder( pose, *loops );
	}

	// set up variables to control designrelaxcycles
	core::Size designrelaxcycles(1);
	core::Energy CurrentEnergy = (*relaxfxn_)( pose );
	core::Energy BestEnergy( 0 );
	core::Energy energychange( 0 );
	core::pose::Pose bestpose;

	core::scoring::ScoreFunctionOP centfxn( ( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::CENTROID_WTS ) ));

	core::pose::Pose native_pose;
	if (  basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		native_pose = pose;
	}

	do {

		// by default the designfxn_ will be soft
		// convert taskfactory into PackerTaskOP
		core::pack::task::PackerTaskOP designtask = designtaskfactory_->create_task_and_apply_taskoperations( pose );
		protocols::simple_moves::PackRotamersMover design_step( designfxn_, designtask );
		design_step.apply( pose );

		if ( debug ) {
			std::string design_tag = "design_step" + ObjexxFCL::lead_zero_string_of( designrelaxcycles, 6 );
			core::io::pdb::dump_pdb( pose, design_tag);
		}

		// default is false, so loop rebuilding must be turned off if not desired
		if ( option[ OptionKeys::loops::no_looprebuild ].value() == true  ) { continue; }
		else {

			// changing the number of inner and outer cycles for LoopMover_Refine_KIC
			// giving the user the chance to increase or decrease cycles
			// with the default LoopMover_Refine_KIC cycles a designrelax simulation takes ~ 60 min
			// without much/any decrease in energy (number of outer cycles is set through the option system)
			if ( !basic::options::option[ basic::options::OptionKeys::loops::max_inner_cycles ].user() ) {
				basic::options::option[ basic::options::OptionKeys::loops::max_inner_cycles ].value(3);
			}


			protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine_alc( loops, relaxfxn_ );
			//   refine_alc.set_native_pose( native_pose );
			refine_alc.apply( pose );

			if ( debug ) {
				std::string loop_tag = "loop_step" + ObjexxFCL::lead_zero_string_of( designrelaxcycles, 6 );
				core::io::pdb::dump_pdb( pose, loop_tag);
			}

			/*
			protocols::loops::kinematic_closure::KinematicMoverOP kin_moverOP( new protocols::loops::kinematic_closure::KinematicMover(2.0));
			kin_moverOP->set_vary_bondangles( true );
			for( protocols::loops::Loops::const_iterator it= loops.begin(), it_end=loops.end();
			it != it_end; ++it ){
			protocols::loops::kinematic_closure::KinematicWrapper kinwrapper( kin_moverOP, *it , 100);

			protocols::moves::MoverOP moverop; // need this to pass kinwrapper to TrialMover
			moverop = kinwrapper;

			// put this into a monte carlo mover
			using protocols::moves::MonteCarloOP;
			using protocols::moves::MonteCarlo;
			MonteCarloOP loop_mc( new MonteCarlo( pose, *centfxn , 0.8) );
			protocols::moves::TrialMoverOP trial_mover = new protocols::moves::TrialMover( moverop, loop_mc );

			protocols::moves::SequenceMoverOP loop_rebuild( new protocols::moves::SequenceMover() );

			for( core::Size nn = 1; nn <= 10; nn++){
			loop_rebuild->add_mover( trial_mover );
			}

			loop_rebuild->apply( pose );
			loop_mc->recover_low( pose );

			}
			*/

		}

		// should we use monte carlo to check the quality of the loop remodel?


		// by default the relaxfxn_ will be full
		protocols::relax::ClassicRelax relax_step( relaxfxn_ );
		relax_step.apply( pose );

		if ( debug ) {
			std::string relax_tag = "relax_step" + ObjexxFCL::lead_zero_string_of( designrelaxcycles, 6 );
			core::io::pdb::dump_pdb( pose, relax_tag);
		}

		CurrentEnergy = (*relaxfxn_)( pose ); // update current energy

		// if CurrentEnergy is less than BestEnergy - read better -
		// so energychange is negative and we need to upate
		energychange =  CurrentEnergy - BestEnergy;
		if ( debug ) { TR << " Current Energy " << CurrentEnergy << " BestEnergy " << BestEnergy << " energychange "<< energychange << " cycles " << designrelaxcycles << std::endl; }
		if ( CurrentEnergy < BestEnergy ) {
			bestpose = pose; // update bestpose
			BestEnergy = CurrentEnergy; // update BestEnergy
			TR << "updated bestpose and bestenergy " << std::endl;
		}

		++designrelaxcycles;
		pose = bestpose;


		// continue design relax cycles until energy change is lessthan -1.00
		// because we probably stuck in an energy well by this point
	} while ( energychange <= -1.00 );

	TR << "design relax mover finished in " << designrelaxcycles - 1 << std::endl;


	// i should do this in a monte carlo setting
	// now that we have a sequence we like - we may want to optimize the loop conformation
	if ( option[ OptionKeys::DenovoProteinDesign::optimize_loops ].user() ) {
		protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kin_moverOP( new protocols::loops::loop_closure::kinematic_closure::KinematicMover() );
		kin_moverOP->set_vary_bondangles( true );
		kin_moverOP->set_temperature( 2.0 );
		 for ( auto const & it : *loops ) {
			protocols::loops::loop_closure::kinematic_closure::KinematicWrapper kinwrapper( kin_moverOP, it , 1000);
			kinwrapper.apply(pose);
		}
	}

}

std::string
DesignRelaxMover::get_name() const {
	return "DesignRelaxMover";
}

/// @brief default ctor

DesignRelaxMover::DesignRelaxMover() : Mover()
{
	Mover::type( "DesignRelaxMover" );
	core::pack::task::TaskFactoryOP designtaskfactory( new core::pack::task::TaskFactory );
	core::scoring::ScoreFunctionOP softfxn(
		core::scoring::ScoreFunctionFactory::create_score_function(SOFT_REP_DESIGN_WTS, SCORE12_PATCH));
	core::scoring::ScoreFunctionOP fullfxn(
		core::scoring::get_score_function() );
	designtaskfactory_ = designtaskfactory;
	designfxn_ = softfxn;
	relaxfxn_ = fullfxn;
}

/// @brief ctor with hand made designtaskfactory and default design/relax scorefunctions
DesignRelaxMover::DesignRelaxMover( core::pack::task::TaskFactoryOP designtaskfactory
) : Mover(), designtaskfactory_(std::move( designtaskfactory ))
{
	Mover::type( "DesignRelaxMover" );
	core::scoring::ScoreFunctionOP softfxn(
		core::scoring::ScoreFunctionFactory::create_score_function(SOFT_REP_DESIGN_WTS, SCORE12_PATCH));
	core::scoring::ScoreFunctionOP fullfxn(
		core::scoring::get_score_function());
	designfxn_ = softfxn;
	relaxfxn_ = fullfxn;
}

DesignRelaxMover::DesignRelaxMover(
	core::pack::task::TaskFactoryOP designtaskfactory,
	core::scoring::ScoreFunctionOP designfxn,
	core::scoring::ScoreFunctionOP relaxfxn
) :
	Mover(), designtaskfactory_(std::move( designtaskfactory )),designfxn_(std::move( designfxn )),relaxfxn_(std::move( relaxfxn ))
{
	Mover::type( "DesignRelaxMover" );
}


DesignRelaxMover::~DesignRelaxMover() = default;

} // DenovoProteinDesign
} // devel
