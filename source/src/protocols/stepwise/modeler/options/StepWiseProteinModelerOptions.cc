// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.cc
/// @brief
/// @detailed
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

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseProteinModelerOptions" );

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
	// If you add a variable, initialize it here, and include in operator = definition below!
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
		use_packer_instead_of_rotamer_trials_ = false;
		min_type_ = "dfpmin_armijo_nonmonotone"; // used to be dfpmin
		min_tolerance_ = 0.000025;
		expand_loop_takeoff_ = false;
		skip_coord_constraints_ = false;
	}

	/// @brief clone the options
	StepWiseProteinModelerOptionsOP
	StepWiseProteinModelerOptions::clone() const
	{
		return new StepWiseProteinModelerOptions( *this );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseProteinModelerOptions &
	StepWiseProteinModelerOptions::operator = ( StepWiseProteinModelerOptions const & src )
	{
		global_optimize_ = src.global_optimize_;
		mapfile_activated_ = src.mapfile_activated_;
		sample_beta_ = src.sample_beta_;
		move_jumps_between_chains_ = src.move_jumps_between_chains_;
		disable_sampling_of_loop_takeoff_ = src.disable_sampling_of_loop_takeoff_;
		cart_min_ = src.cart_min_;
		n_sample_ = src.n_sample_;
		filter_native_big_bins_ = src.filter_native_big_bins_;
		allow_virtual_side_chains_ = src.allow_virtual_side_chains_;
		prepack_ = src.prepack_;
		centroid_output_ = src.centroid_output_;
		centroid_screen_ = src.centroid_screen_;
		centroid_score_diff_cut_ = src.centroid_score_diff_cut_;
		centroid_weights_ = src.centroid_weights_;
		nstruct_centroid_ = src.nstruct_centroid_;
		ghost_loops_ = src.ghost_loops_;
		ccd_close_ = src.ccd_close_;
		cluster_by_all_atom_rmsd_ = src.cluster_by_all_atom_rmsd_;
		pack_weights_ = src.pack_weights_;
		use_packer_instead_of_rotamer_trials_ = src.use_packer_instead_of_rotamer_trials_;
		min_type_ = src.min_type_;
		min_tolerance_ = src.min_tolerance_;
		expand_loop_takeoff_ = src.expand_loop_takeoff_;
		skip_coord_constraints_ = src.skip_coord_constraints_;
		frag_files_ = src.frag_files_;
		bridge_res_ = src.bridge_res_;
		return *this;
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
		use_packer_instead_of_rotamer_trials_ = option[ basic::options::OptionKeys::stepwise::protein::use_packer_instead_of_rotamer_trials ]();
		min_type_ = option[ basic::options::OptionKeys::stepwise::protein::min_type ]();
		min_tolerance_ = option[ basic::options::OptionKeys::stepwise::protein::min_tolerance ]();
		expand_loop_takeoff_ = option[ basic::options::OptionKeys::stepwise::protein::expand_loop_takeoff ]() ;
		skip_coord_constraints_ = option[ OptionKeys::stepwise::protein::skip_coord_constraints ]();
		frag_files_ = option[ OptionKeys::in::file::frag_files ]();
		bridge_res_ = option[ OptionKeys::stepwise::protein::bridge_res ]();

		//		pack_weights_ = "stepwise/protein/pack_no_hb_env_dep.wts";
		if ( option[ basic::options::OptionKeys::score::pack_weights ].user() )	pack_weights_ = option[ score::pack_weights ]();

	}

} //options
} //modeler
} //stepwise
} //protocols
