// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/options/StepWiseBasicOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_options_StepWiseBasicOptions_HH
#define INCLUDED_protocols_stepwise_options_StepWiseBasicOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/options/StepWiseBasicOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

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

	core::Size const & sampler_num_pose_kept() const { return sampler_num_pose_kept_; }
	void set_sampler_num_pose_kept( core::Size const & setting ){ sampler_num_pose_kept_ = setting; }

	core::Real const & cluster_rmsd() const { return cluster_rmsd_; }
	void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

	core::Size const & num_pose_minimize() const { return num_pose_minimize_; }
	void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

	core::Size const & num_random_samples() const { return num_random_samples_; };
	void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }

	core::Size const & max_tries_multiplier_for_ccd() const { return max_tries_multiplier_for_ccd_; };
	void set_max_tries_multiplier_for_ccd( core::Size const & setting ){ max_tries_multiplier_for_ccd_ = setting; }

	void set_atr_rep_screen( bool const & setting ){ atr_rep_screen_ = setting; }
	bool atr_rep_screen() const{ return atr_rep_screen_; }

	void set_atr_rep_screen_for_docking( bool const & setting ){ atr_rep_screen_for_docking_ = setting; }
	bool atr_rep_screen_for_docking() const{ return atr_rep_screen_for_docking_; }

	utility::vector1< std::string > const & VDW_rep_screen_info() const { return VDW_rep_screen_info_; }
	void set_VDW_rep_screen_info( utility::vector1< std::string > const & setting ){ VDW_rep_screen_info_ = setting; }

	void set_rmsd_screen( core::Real const & setting ){ rmsd_screen_ = setting; }
	core::Real rmsd_screen() const{ return rmsd_screen_; }

	bool const & output_minimized_pose_list() const { return output_minimized_pose_list_; }
	void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

	bool const & output_cluster_size() const { return output_cluster_size_; }
	void set_output_cluster_size( bool const & setting ){ output_cluster_size_ = setting; }

	void set_min_type( std::string const & setting ){ min_type_ = setting; }
	std::string min_type() const{ return min_type_; }

	void set_min_tolerance( core::Real const & setting ){ min_tolerance_ = setting; }
	core::Real min_tolerance() const{ return min_tolerance_; }

	bool const & vary_rna_bond_geometry() const { return vary_rna_bond_geometry_; }
	void set_vary_rna_bond_geometry( bool const & setting ){ vary_rna_bond_geometry_ = setting; }

	bool const & vary_polar_hydrogen_geometry() const { return vary_polar_hydrogen_geometry_; }
	void set_vary_polar_hydrogen_geometry( bool const & setting ){ vary_polar_hydrogen_geometry_ = setting; }

	void set_use_packer_instead_of_rotamer_trials( bool const & setting ){ use_packer_instead_of_rotamer_trials_ = setting; }
	bool use_packer_instead_of_rotamer_trials() const{ return use_packer_instead_of_rotamer_trials_; }

	void set_lores( bool const & setting ){ lores_ = setting; }
	bool lores() const{ return lores_; }

	void set_verbose_sampler( bool const & setting ){ verbose_sampler_ = setting; }
	bool verbose_sampler() const{ return verbose_sampler_; }

protected:

	void
	initialize_from_command_line();

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
	bool use_packer_instead_of_rotamer_trials_;
	bool lores_;
	bool verbose_sampler_;
	std::string sampler_silent_file_;

};

} //options
} //stepwise
} //protocols

#endif
