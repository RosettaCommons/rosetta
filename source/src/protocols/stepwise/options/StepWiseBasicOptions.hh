// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/options/StepWiseBasicOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_options_StepWiseBasicOptions_HH
#define INCLUDED_protocols_stepwise_options_StepWiseBasicOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/options/StepWiseBasicOptions.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseMinimizer.hh> // for MinimizerMode enum
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/options/OptionCollection.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace stepwise {
namespace options {

class StepWiseBasicOptions: public virtual basic::resource_manager::ResourceOptions {

public:

	//constructor
	StepWiseBasicOptions();

	StepWiseBasicOptions( StepWiseBasicOptions const & src );

	//destructor
	~StepWiseBasicOptions();

public:

	StepWiseBasicOptionsOP clone() const;

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "StepWiseBasicOptions";}

	std::string const & silent_file() const { return silent_file_; }
	void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

	std::string const & sampler_silent_file() const { return sampler_silent_file_; }
	void set_sampler_silent_file( std::string const & setting ){ sampler_silent_file_ = setting; }

	core::Size sampler_num_pose_kept() const { return sampler_num_pose_kept_; }
	void set_sampler_num_pose_kept( core::Size const & setting ){ sampler_num_pose_kept_ = setting; }

	core::Real cluster_rmsd() const { return cluster_rmsd_; }
	void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

	core::Size num_pose_minimize() const { return num_pose_minimize_; }
	void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

	core::Size num_random_samples() const { return num_random_samples_; };
	void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }

	core::Size max_tries_multiplier_for_ccd() const { return max_tries_multiplier_for_ccd_; };
	void set_max_tries_multiplier_for_ccd( core::Size const & setting ){ max_tries_multiplier_for_ccd_ = setting; }

	void set_atr_rep_screen( bool const & setting ){ atr_rep_screen_ = setting; }
	bool atr_rep_screen() const{ return atr_rep_screen_; }

	void set_atr_rep_screen_for_docking( bool const & setting ){ atr_rep_screen_for_docking_ = setting; }
	bool atr_rep_screen_for_docking() const{ return atr_rep_screen_for_docking_; }

	utility::vector1< std::string > const & VDW_rep_screen_info() const { return VDW_rep_screen_info_; }
	void set_VDW_rep_screen_info( utility::vector1< std::string > const & setting ){ VDW_rep_screen_info_ = setting; }

	void set_rmsd_screen( core::Real const & setting ){ rmsd_screen_ = setting; }
	core::Real rmsd_screen() const{ return rmsd_screen_; }

	bool output_minimized_pose_list() const { return output_minimized_pose_list_; }
	void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

	bool output_cluster_size() const { return output_cluster_size_; }
	void set_output_cluster_size( bool const & setting ){ output_cluster_size_ = setting; }

	void set_min_type( std::string const & setting ){ min_type_ = setting; }
	std::string min_type() const{ return min_type_; }

	void set_min_tolerance( core::Real const & setting ){ min_tolerance_ = setting; }
	core::Real min_tolerance() const{ return min_tolerance_; }

	bool vary_rna_bond_geometry() const { return vary_rna_bond_geometry_; }
	void set_vary_rna_bond_geometry( bool const & setting ){ vary_rna_bond_geometry_ = setting; }

	bool vary_polar_hydrogen_geometry() const { return vary_polar_hydrogen_geometry_; }
	void set_vary_polar_hydrogen_geometry( bool const & setting ){ vary_polar_hydrogen_geometry_ = setting; }

	void set_disallow_pack_polar_hydrogens( bool const & setting ){ disallow_pack_polar_hydrogens_ = setting; }
	bool disallow_pack_polar_hydrogens() const { return disallow_pack_polar_hydrogens_; }

	void set_use_packer_instead_of_rotamer_trials( bool const & setting ){ use_packer_instead_of_rotamer_trials_ = setting; }
	bool use_packer_instead_of_rotamer_trials() const{ return use_packer_instead_of_rotamer_trials_; }

	bool mapfile_activated() const { return mapfile_activated_; }
	void set_mapfile_activated( bool const & setting ){ mapfile_activated_ = setting; }

	void set_lores( bool const & setting ){ lores_ = setting; }
	bool lores() const{ return lores_; }

	void set_verbose_sampler( bool const & setting ){ verbose_sampler_ = setting; }
	bool verbose_sampler() const{ return verbose_sampler_; }

	void set_minimize_waters( bool const & setting ){ minimize_waters_ = setting; }
	bool minimize_waters() const{ return minimize_waters_; }

	void set_hydrate_magnesiums( bool const & setting ){ hydrate_magnesiums_ = setting; }
	bool hydrate_magnesiums() const{ return hydrate_magnesiums_; }

	void set_test_all_mg_hydration_frames( bool const & setting ){ test_all_mg_hydration_frames_ = setting; }
	bool test_all_mg_hydration_frames() const{ return test_all_mg_hydration_frames_; }

	modeler::MinimizerMode const & minimizer_mode() const { return minimizer_mode_; }
	void set_minimizer_mode( modeler::MinimizerMode const & setting ){ minimizer_mode_ = setting; }

	Size n_cycles() const { return n_cycles_; }
	void set_n_cycles( Size const setting ) { n_cycles_ = setting; }

	core::Real thermal_sampler_temperature() const { return thermal_sampler_temperature_; }

	void set_thermal_sampler_output_min_pose( bool const setting ) { thermal_sampler_output_min_pose_ = setting; }
	bool thermal_sampler_output_min_pose() const{ return thermal_sampler_output_min_pose_; }

	void set_sample_pH( bool const setting ){ sample_pH_ = setting; }
	bool sample_pH() const{ return sample_pH_; }

	static void
	list_options_read( utility::options::OptionKeyList & opt );

protected:

	void
	initialize_from_command_line();
	void
	initialize_from_options_collection( utility::options::OptionCollection const & options );

	void
	initialize_variables();

private:

	std::string silent_file_;
	core::Size sampler_num_pose_kept_;
	core::Real cluster_rmsd_;
	core::Size num_pose_minimize_;
	core::Size num_random_samples_;
	core::Size max_tries_multiplier_for_ccd_;
	bool atr_rep_screen_;
	bool atr_rep_screen_for_docking_;
	core::Real rmsd_screen_;
	utility::vector1< std::string > VDW_rep_screen_info_;

	bool output_minimized_pose_list_;
	bool output_cluster_size_;
	std::string min_type_;
	core::Real min_tolerance_;
	bool vary_rna_bond_geometry_;
	bool vary_polar_hydrogen_geometry_;
	bool disallow_pack_polar_hydrogens_;
	bool use_packer_instead_of_rotamer_trials_;
	bool lores_;
	bool verbose_sampler_;
	std::string sampler_silent_file_;
	bool minimize_waters_;
	bool hydrate_magnesiums_;
	bool test_all_mg_hydration_frames_;
	bool mapfile_activated_;

	// thermal_sampler
	modeler::MinimizerMode minimizer_mode_;
	Size n_cycles_;
	core::Real thermal_sampler_temperature_;
	bool thermal_sampler_output_min_pose_;
	bool sample_pH_;

};

} //options
} //stepwise
} //protocols

#endif
