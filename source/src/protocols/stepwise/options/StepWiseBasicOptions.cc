// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/options/StepWiseBasicOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/options/StepWiseBasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "protocols.stepwise.options.StepWiseBasicOptions" );

namespace protocols {
namespace stepwise {
namespace options {

//Constructor
StepWiseBasicOptions::StepWiseBasicOptions()
{
	initialize_variables();
}

//Destructor
StepWiseBasicOptions::~StepWiseBasicOptions()
{}

/// @brief copy constructor
StepWiseBasicOptions::StepWiseBasicOptions( StepWiseBasicOptions const & src ) :
	ResourceOptions ( src )
{
	*this = src;
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicOptions::initialize_variables(){
	silent_file_ = "";
	sampler_silent_file_ = "";
	sampler_num_pose_kept_ = 0; // signal to use separate defaults in protein/RNA modeler
	cluster_rmsd_ = 0.0; // signal for clusterers to use their default value
	num_pose_minimize_ = 0; // signal to minimize all.
	num_random_samples_ = 20;
	max_tries_multiplier_for_ccd_ = 10;
	atr_rep_screen_ = true;
	rmsd_screen_ = 0.0;
	output_minimized_pose_list_ = false;
	output_cluster_size_ = false;
	min_type_ = "dfpmin_armijo_nonmonotone"; // used to be dfpmin
	min_tolerance_ = 0.000025;
	vary_rna_bond_geometry_ = false;
	vary_polar_hydrogen_geometry_ = false;
	use_packer_instead_of_rotamer_trials_ = false;
}

/// @brief clone the options
StepWiseBasicOptionsOP
StepWiseBasicOptions::clone() const
{
	return StepWiseBasicOptionsOP( new StepWiseBasicOptions( *this ) );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicOptions::initialize_from_command_line(){
	sampler_silent_file_ = option[ basic::options::OptionKeys::stepwise::sampler_silent_file ]();
	num_random_samples_ = option[ basic::options::OptionKeys::stepwise::num_random_samples ]();
	max_tries_multiplier_for_ccd_ = option[ basic::options::OptionKeys::stepwise::max_tries_multiplier_for_ccd ]();
	sampler_num_pose_kept_ = option[ basic::options::OptionKeys::stepwise::rna::sampler_num_pose_kept ]();
	if ( option[ basic::options::OptionKeys::cluster::radius ].user() ) cluster_rmsd_ = option[ basic::options::OptionKeys::cluster::radius ]();
	atr_rep_screen_ = option[ basic::options::OptionKeys::stepwise::atr_rep_screen ]();
	rmsd_screen_ = option[ basic::options::OptionKeys::stepwise::rmsd_screen ]();
	VDW_rep_screen_info_ = option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ]();
	if ( option[ basic::options::OptionKeys::stepwise::num_pose_minimize ].user() ) num_pose_minimize_ = option[ basic::options::OptionKeys::stepwise::num_pose_minimize ]();
	min_type_ = option[ basic::options::OptionKeys::stepwise::min_type ]();
	min_tolerance_ = option[ basic::options::OptionKeys::stepwise::min_tolerance ]();
	vary_rna_bond_geometry_ = option[ basic::options::OptionKeys::rna::vary_geometry ]();
	vary_polar_hydrogen_geometry_ = option[ basic::options::OptionKeys::stepwise::vary_polar_hydrogen_geometry ]();
	use_packer_instead_of_rotamer_trials_ = option[ basic::options::OptionKeys::stepwise::protein::use_packer_instead_of_rotamer_trials ]();
	output_minimized_pose_list_ = option[ basic::options::OptionKeys::stepwise::output_minimized_pose_list ]();
	output_cluster_size_ = option[ basic::options::OptionKeys::stepwise::output_cluster_size ]();

}


} //options
} //stepwise
} //protocols
