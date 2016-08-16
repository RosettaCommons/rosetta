// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/packer/StepWiseMasterPacker.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/packer/StepWiseMasterPacker.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/screener/O2PrimeScreener.hh>
#include <protocols/stepwise/screener/PhosphateScreener.hh>
#include <protocols/stepwise/screener/PackScreener.hh>
#include <protocols/stepwise/screener/SugarInstantiator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.packer.StepWiseMasterPacker" );

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Unified treatment of protein side-chains, RNA 2'-OH packing/trials for stepwise modeling,
//   including smart visualization if there is not information to place the groups.
//
// Also now including 'packing' of terminal phosphates. This is a separate step in which phosphates
//  go on as backbone variants. I also tried a mode where those phosphates were 'side chains' to
//  enable simultaneous optimization with 2'-OH, but code got really slow due to rotamer explosion.
//
//    -- Rhiju, 2014.
//
// To do (maybe):
//   - test instantiation of RNA sugars (would need to update score term so that a single 2'-OH
//      could actually tip the balance and favor instantiation.
//   - develop code to remove/add peptide termini -- N-acetylation and C-terminal methylamidation;
//      currently those variants are just instantiated before modeler and left on even if they
//      don't contact anything.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;
using namespace protocols::stepwise::screener;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

//Constructor
StepWiseMasterPacker::StepWiseMasterPacker( working_parameters::StepWiseWorkingParametersCOP working_parameters,
	options::StepWiseModelerOptionsCOP options ):
	working_parameters_( working_parameters),
	options_( options )
{}

//Destructor
StepWiseMasterPacker::~StepWiseMasterPacker()
{}

//////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::initialize( pose::Pose const & pose ) {

	if ( options_->sampler_perform_phosphate_pack() ) {
		phosphate_sampler_ = rna::phosphate::MultiPhosphateSamplerOP( new rna::phosphate::MultiPhosphateSampler( pose ) );
		phosphate_sampler_->set_moving_partition_res( working_parameters_->working_moving_partition_res() );
		phosphate_sampler_->set_scorefxn( rna::phosphate::get_phosphate_scorefxn( scorefxn_->energy_method_options() ) );
		phosphate_sampler_->set_force_phosphate_instantiation( options_->force_phosphate_instantiation() );
	}

	if ( options_->o2prime_legacy_mode() ) { // this is the only action that is non-const for the pose... deprecate?
		rna::check_instantiated_O2Prime_hydrogen( pose ); // i.e., instantiate all 2'-OH.
		o2prime_packer_ = rna::o2prime::O2PrimePackerOP( new rna::o2prime::O2PrimePacker( pose, scorefxn_, working_parameters_->working_moving_partition_res() ) );
	}

	initialize_packer();
}

/////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::set_working_pack_res( utility::vector1< core::Size > const & setting ){
	working_pack_res_ = setting;
	packer_->set_working_pack_res( working_pack_res_ );
}

/////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::initialize_packer() {
	packer_ = packer::get_packer( scorefxn_, get_all_working_moving_res( working_parameters_ ), options_ );
	if ( options_->o2prime_legacy_mode() ) {
		packer_->set_pack_o2prime_hydrogens( false ); // in legacy mode, HO2' was handled by separate packer
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::do_prepack( pose::Pose & pose ) {
	if ( phosphate_sampler_ != 0 ) phosphate_sampler_->do_prepack( pose, working_parameters_->working_moving_res_list() ); // must be fixed.
	packer_->do_prepack( pose );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::add_packer_screeners( utility::vector1< screener::StepWiseScreenerOP > & screeners,
	pose::Pose const & pose,
	pose::PoseOP sugar_instantiation_pose ){

	// testing -- OK to move down here? No won't be inherited by the actual pose, right?
	if ( options_->sampler_perform_phosphate_pack() ) {
		screeners.push_back( screener::StepWiseScreenerOP( new PhosphateScreener( phosphate_sampler_ ) ) );
	}

	if ( options_->sampler_try_sugar_instantiation() &&
			working_parameters_->floating_base() ) {
		screeners.push_back( stepwise::screener::StepWiseScreenerOP( new SugarInstantiator( *sugar_instantiation_pose, working_parameters_->working_moving_res() ) ) );
	}

	if ( options_->o2prime_legacy_mode() ) {
		screeners.push_back( screener::StepWiseScreenerOP( new O2PrimeScreener( o2prime_packer_ ) ) );
		if ( !protein::contains_protein( pose ) ) return; // don't put in PackScreener below.
	}

	packer_pose_ = pose.clone();
	screeners.push_back( screener::StepWiseScreenerOP( new PackScreener( *packer_pose_, packer_ ) ) ); // includes HO2' for RNA.
}

///////////////////////////////////////////////////////////////////////////
void
StepWiseMasterPacker::reset( pose::Pose const & pose ){
	if ( phosphate_sampler_ != 0 ) phosphate_sampler_->reset( pose );
	packer_->reset( pose );
}

///////////////////////////////////////////////////////////////////////////
packer::StepWisePackerCOP
StepWiseMasterPacker::packer(){
	return packer_;
}

} //packer
} //modeler
} //stepwise
} //protocols
