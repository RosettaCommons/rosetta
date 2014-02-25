// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::stepwise::sampling::rna;

static basic::Tracer TR( "protocols.stepwise.monte_carlo.rna.StepWiseRNA_MonteCarloOptions" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	//Constructor
	StepWiseRNA_MonteCarloOptions::StepWiseRNA_MonteCarloOptions():
		verbose_scores_( false ),
		force_centroid_interaction_( true ),
		use_phenix_geo_( true ),
		skip_deletions_( false ),
		erraser_( true ),
		allow_internal_hinge_moves_( true ),
		allow_internal_local_moves_( false ),
		num_random_samples_( 20 ),
		cycles_( 500 ),
		add_delete_frequency_( 0.5 ),
		minimize_single_res_frequency_( 0.0 ),
		minimizer_allow_variable_bond_geometry_( true ),
		minimizer_vary_bond_geometry_frequency_( 0.1 ),
		switch_focus_frequency_( 0.5 ),
		just_min_after_mutation_frequency_( 0.5 ),
		temperature_( 1.0 ),
		allow_skip_bulge_( false ),
		allow_from_scratch_( false ),
		allow_split_off_( true ),
		virtual_sugar_keep_base_fixed_( false ),
		constraint_x0_( 0.0 ),
		constraint_tol_( 0.0 ),
		make_movie_( false ),
		sampler_perform_phosphate_pack_( true ),
		rebuild_bulge_mode_( false )
	{}

	//Destructor
	StepWiseRNA_MonteCarloOptions::~StepWiseRNA_MonteCarloOptions()
	{}


	///////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_MonteCarloOptions::initialize_from_command_line() {
		set_verbose_scores( option[ OptionKeys::stepwise::monte_carlo::verbose_scores ]() );
		set_skip_deletions( option[ OptionKeys::stepwise::monte_carlo::skip_deletions ]() );
		set_num_random_samples( option[ OptionKeys::stepwise::rna::num_random_samples ]() );
		set_erraser( option[ OptionKeys::stepwise::rna::erraser ]() );
		set_cycles( option[ OptionKeys::stepwise::monte_carlo::cycles ]() );
		set_add_delete_frequency( option[ OptionKeys::stepwise::monte_carlo::add_delete_frequency ]() );
		set_minimize_single_res_frequency( option[ OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency ]() );
		set_switch_focus_frequency( option[ OptionKeys::stepwise::monte_carlo::switch_focus_frequency ]() );
		set_sample_res( option[ OptionKeys::stepwise::rna::sample_res ]() );
		set_just_min_after_mutation_frequency( option[ OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency ]() );
		set_allow_internal_hinge_moves( option[ OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves ]() );
		set_allow_internal_local_moves( option[ OptionKeys::stepwise::monte_carlo::allow_internal_local_moves ]() );
		set_allow_skip_bulge( option[ OptionKeys::stepwise::monte_carlo::allow_skip_bulge ]() );
		set_allow_from_scratch( option[ OptionKeys::stepwise::monte_carlo::allow_from_scratch ]() );
		set_allow_split_off( option[ OptionKeys::stepwise::monte_carlo::allow_split_off ]() );
		set_temperature( option[ OptionKeys::stepwise::monte_carlo::temperature ]() );
		set_extra_minimize_res( option[ OptionKeys::stepwise::monte_carlo::extra_min_res ]() );
		set_syn_chi_res_list( option[ OptionKeys::stepwise::rna::force_syn_chi_res_list]() );
		set_terminal_res( option[ OptionKeys::stepwise::rna::terminal_res]() );
		set_bulge_res( option[ basic::options::OptionKeys::stepwise::rna::bulge_res ]() );
		set_virtual_sugar_keep_base_fixed( option[ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed ]() );
		force_centroid_interaction_ = true; // note default is different from stepwise enumeration
		if ( option[ OptionKeys::stepwise::rna::force_centroid_interaction ].user() ) set_force_centroid_interaction( option[ OptionKeys::stepwise::rna::force_centroid_interaction ]() );
		set_minimizer_allow_variable_bond_geometry( option[ OptionKeys::stepwise::monte_carlo::allow_variable_bond_geometry ]() );
		set_use_phenix_geo(  option[ OptionKeys::rna::corrected_geo ] );
		set_constraint_x0(	option[ OptionKeys::stepwise::monte_carlo::constraint_x0 ]() );
		set_constraint_tol(	option[ OptionKeys::stepwise::monte_carlo::constraint_tol ]() );
		set_make_movie(	option[ OptionKeys::stepwise::monte_carlo::make_movie ]() );
		sampler_perform_phosphate_pack_ = option[ OptionKeys::stepwise::rna::sampler_perform_phosphate_pack ]();
		rebuild_bulge_mode_ = option[ OptionKeys::stepwise::rna::rebuild_bulge_mode ]();
	}


	//////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_ModelerOptionsOP
	StepWiseRNA_MonteCarloOptions::setup_modeler_options() const{

		StepWiseRNA_ModelerOptionsOP modeler_options = new StepWiseRNA_ModelerOptions;
		modeler_options->set_choose_random( true );
		modeler_options->set_force_centroid_interaction( force_centroid_interaction() );
		modeler_options->set_use_phenix_geo( use_phenix_geo() );
		modeler_options->set_kic_sampling_if_relevant( erraser() );
		modeler_options->set_num_random_samples( num_random_samples() );
		modeler_options->set_num_pose_minimize( 1 );
		modeler_options->set_minimizer_allow_variable_bond_geometry( minimizer_allow_variable_bond_geometry() );
		modeler_options->set_minimizer_vary_bond_geometry_frequency( minimizer_vary_bond_geometry_frequency() );
		modeler_options->set_virtual_sugar_keep_base_fixed( virtual_sugar_keep_base_fixed() );
		modeler_options->set_sampler_perform_phosphate_pack( sampler_perform_phosphate_pack() );
		if ( rebuild_bulge_mode_ )	modeler_options->set_force_centroid_interaction( false );

		return modeler_options;
	}


} //rna
} //monte_carlo
} //stepwise
} //protocols
