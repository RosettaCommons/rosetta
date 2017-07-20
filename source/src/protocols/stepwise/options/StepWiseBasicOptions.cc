// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/options/StepWiseBasicOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/options/StepWiseBasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>

#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.options.StepWiseBasicOptions" );

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
	atr_rep_screen_for_docking_ = false;
	rmsd_screen_ = 0.0;
	output_minimized_pose_list_ = false;
	output_cluster_size_ = false;
	min_type_ = "lbfgs_armijo_nonmonotone"; // used to be dfpmin
	min_tolerance_ = 0.000025;
	vary_rna_bond_geometry_ = false;
	vary_polar_hydrogen_geometry_ = false;
	disallow_pack_polar_hydrogens_ = false;
	use_packer_instead_of_rotamer_trials_ = false;
	mapfile_activated_ = false;
	lores_ = false;
	verbose_sampler_ = false;
	minimize_waters_ = false;
	hydrate_magnesiums_ = false;
	test_all_mg_hydration_frames_ = false;
	minimizer_mode_ = modeler::TRADITIONAL_MINIMIZER;
	n_cycles_ = 100;
	thermal_sampler_temperature_ = 0.5;
	thermal_sampler_output_min_pose_ = false;
	sample_pH_ = false;
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
	initialize_from_options_collection( option );
}

void
StepWiseBasicOptions::initialize_from_options_collection( utility::options::OptionCollection const & the_options ) {
	using namespace basic::options;

	sampler_silent_file_ = the_options[ OptionKeys::stepwise::sampler_silent_file ]();
	num_random_samples_ = the_options[ OptionKeys::stepwise::num_random_samples ]();
	max_tries_multiplier_for_ccd_ = the_options[ OptionKeys::stepwise::max_tries_multiplier_for_ccd ]();
	sampler_num_pose_kept_ = the_options[ OptionKeys::stepwise::rna::sampler_num_pose_kept ]();
	if ( the_options[ OptionKeys::cluster::radius ].user() ) cluster_rmsd_ = the_options[ OptionKeys::cluster::radius ]();
	atr_rep_screen_ = the_options[ OptionKeys::stepwise::atr_rep_screen ]();
	atr_rep_screen_for_docking_ = the_options[ OptionKeys::stepwise::atr_rep_screen_for_docking ]();
	rmsd_screen_ = the_options[ OptionKeys::stepwise::rmsd_screen ]();
	VDW_rep_screen_info_ = the_options[ OptionKeys::stepwise::rna::VDW_rep_screen_info ]();
	if ( the_options[ OptionKeys::stepwise::num_pose_minimize ].user() ) num_pose_minimize_ = the_options[ OptionKeys::stepwise::num_pose_minimize ]();
	min_type_ = the_options[ OptionKeys::stepwise::min_type ]();
	min_tolerance_ = the_options[ OptionKeys::stepwise::min_tolerance ]();
	vary_rna_bond_geometry_ = the_options[ OptionKeys::rna::vary_geometry ]();
	vary_polar_hydrogen_geometry_ = the_options[ OptionKeys::stepwise::polar_hydrogens::vary_polar_hydrogen_geometry ]();
	disallow_pack_polar_hydrogens_ = the_options[ OptionKeys::stepwise::polar_hydrogens::disallow_pack_polar_hydrogens ]();
	use_packer_instead_of_rotamer_trials_ = the_options[ OptionKeys::stepwise::protein::use_packer_instead_of_rotamer_trials ]();
	output_minimized_pose_list_ = the_options[ OptionKeys::stepwise::output_minimized_pose_list ]();
	output_cluster_size_ = the_options[ OptionKeys::stepwise::output_cluster_size ]();
	mapfile_activated_ = the_options[ OptionKeys::edensity::mapfile ].user();

	lores_ = the_options[ OptionKeys::stepwise::lores ]();
	verbose_sampler_ = the_options[ OptionKeys::stepwise::verbose_sampler ]();
	minimize_waters_ = the_options[ OptionKeys::stepwise::minimize_waters ]();
	hydrate_magnesiums_ = the_options[ OptionKeys::magnesium::hydrate ]();
	test_all_mg_hydration_frames_ = the_options[ magnesium::all_hydration_frames ]();

	if ( the_options[ OptionKeys::stepwise::minimizer_mode ]() == "THERMAL_SAMPLER" ) {
		minimizer_mode_ = modeler::THERMAL_SAMPLER;
	} else {
		minimizer_mode_ = modeler::TRADITIONAL_MINIMIZER;
	}
	n_cycles_ = the_options[ OptionKeys::recces::n_cycle ]();
	auto temps = the_options[ OptionKeys::recces::temps ]();
	thermal_sampler_temperature_ = ( temps.size() > 0 )  ? temps[ 1 ] : 0.5;
	thermal_sampler_output_min_pose_ = the_options[ OptionKeys::recces::output_min_pose ]();
	sample_pH_ = the_options[ OptionKeys::pH::pH_mode ]();
}

void
StepWiseBasicOptions::list_options_read( utility::options::OptionKeyList & opt ) {

	using namespace basic::options;

	opt + OptionKeys::stepwise::sampler_silent_file
		+ OptionKeys::stepwise::num_random_samples
		+ OptionKeys::stepwise::max_tries_multiplier_for_ccd
		+ OptionKeys::stepwise::rna::sampler_num_pose_kept
		+ OptionKeys::cluster::radius
		+ OptionKeys::stepwise::atr_rep_screen
		+ OptionKeys::stepwise::atr_rep_screen_for_docking
		+ OptionKeys::stepwise::rmsd_screen
		+ OptionKeys::stepwise::rna::VDW_rep_screen_info
		+ OptionKeys::stepwise::num_pose_minimize
		+ OptionKeys::stepwise::min_type
		+ OptionKeys::stepwise::min_tolerance
		+ OptionKeys::rna::vary_geometry
		+ OptionKeys::stepwise::polar_hydrogens::vary_polar_hydrogen_geometry
		+ OptionKeys::stepwise::polar_hydrogens::disallow_pack_polar_hydrogens
		+ OptionKeys::stepwise::protein::use_packer_instead_of_rotamer_trials
		+ OptionKeys::stepwise::output_minimized_pose_list
		+ OptionKeys::stepwise::output_cluster_size
		+ OptionKeys::edensity::mapfile
		+ OptionKeys::stepwise::lores
		+ OptionKeys::stepwise::verbose_sampler
		+ OptionKeys::stepwise::minimize_waters
		+ OptionKeys::magnesium::hydrate
		+ magnesium::all_hydration_frames
		+ OptionKeys::stepwise::minimizer_mode
		+ OptionKeys::recces::n_cycle
		+ OptionKeys::recces::temps
		+ OptionKeys::recces::output_min_pose
		+ OptionKeys::pH::pH_mode;
}



} //options
} //stepwise
} //protocols
