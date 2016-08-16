// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseProteinModelerOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

/////////////////////////////////////////////////////////////////////////////////////
//Constructor
StepWiseProteinModelerOptions::StepWiseProteinModelerOptions()
{
	initialize_variables();
}

/////////////////////////////////////////////////////////////////////////////////////
//Destructor
StepWiseProteinModelerOptions::~StepWiseProteinModelerOptions()
{
}

/// @brief copy constructor
StepWiseProteinModelerOptions::StepWiseProteinModelerOptions( StepWiseProteinModelerOptions const & src ) :
	ResourceOptions ( src )
{
	*this = src;
}


/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinModelerOptions::initialize_variables(){

	global_optimize_ = false;
	mapfile_activated_ = false;
	sample_beta_ = false;
	move_jumps_between_chains_ = false;
	disable_sampling_of_loop_takeoff_ = false;
	cart_min_ = false;
	n_sample_ = 18;
	filter_native_big_bins_ = false;
	allow_virtual_side_chains_ = false;
	prepack_ = false;
	centroid_output_ = false;
	centroid_screen_ = false;
	centroid_score_diff_cut_ = 20.0;
	centroid_weights_ = "score3.wts";
	nstruct_centroid_ = 0;
	ghost_loops_ = false;
	ccd_close_ = false; // should this be in here, or in job parameters?
	cluster_by_all_atom_rmsd_ = false;
	pack_weights_ = "";
	expand_loop_takeoff_ = false;
	skip_coord_constraints_ = false;
}

/// @brief clone the options
StepWiseProteinModelerOptionsOP
StepWiseProteinModelerOptions::clone() const
{
	return StepWiseProteinModelerOptionsOP( new StepWiseProteinModelerOptions( *this ) );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinModelerOptions::initialize_from_command_line(){

	global_optimize_ = option[ basic::options::OptionKeys::stepwise::protein::global_optimize ]();
	mapfile_activated_ = option[edensity::mapfile].user();
	sample_beta_ = option[ basic::options::OptionKeys::stepwise::protein::sample_beta ]();
	move_jumps_between_chains_ = option[ basic::options::OptionKeys::stepwise::protein::move_jumps_between_chains ]();
	disable_sampling_of_loop_takeoff_ = option[ basic::options::OptionKeys::stepwise::protein::disable_sampling_of_loop_takeoff ]();
	cart_min_ = option[ basic::options::OptionKeys::stepwise::protein::cart_min ]();
	n_sample_ = option[ basic::options::OptionKeys::stepwise::protein::n_sample ]();
	filter_native_big_bins_ = option[ basic::options::OptionKeys::stepwise::protein::filter_native_big_bins ]();
	allow_virtual_side_chains_ = option[ basic::options::OptionKeys::stepwise::protein::allow_virtual_side_chains ]();
	prepack_ = option[ basic::options::OptionKeys::stepwise::protein::protein_prepack ]();
	centroid_output_ = option[ basic::options::OptionKeys::stepwise::protein::centroid_output ]();
	centroid_screen_ = option[ basic::options::OptionKeys::stepwise::protein::centroid_screen ]();
	centroid_score_diff_cut_ = option[ basic::options::OptionKeys::stepwise::protein::centroid_score_diff_cut ]();
	centroid_weights_ = option[ basic::options::OptionKeys::stepwise::protein::centroid_weights ]();
	nstruct_centroid_ = option[ basic::options::OptionKeys::stepwise::protein::nstruct_centroid ]();
	ghost_loops_ = option[ basic::options::OptionKeys::stepwise::protein::ghost_loops ]();
	ccd_close_ = option[ basic::options::OptionKeys::stepwise::protein::ccd_close ]();
	cluster_by_all_atom_rmsd_ = option[ basic::options::OptionKeys::stepwise::protein::cluster_by_all_atom_rmsd ]();
	expand_loop_takeoff_ = option[ basic::options::OptionKeys::stepwise::protein::expand_loop_takeoff ]() ;
	skip_coord_constraints_ = option[ OptionKeys::stepwise::protein::skip_coord_constraints ]();
	frag_files_ = option[ OptionKeys::in::file::frag_files ]();
	bridge_res_ = option[ OptionKeys::stepwise::protein::bridge_res ]();

	//  pack_weights_ = "stepwise/protein/pack_no_hb_env_dep.wts";
	if ( option[ basic::options::OptionKeys::score::pack_weights ].user() ) pack_weights_ = option[ score::pack_weights ]();
}

} //options
} //modeler
} //stepwise
} //protocols
