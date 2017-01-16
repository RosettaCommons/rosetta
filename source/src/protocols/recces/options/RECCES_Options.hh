// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/options/RECCES_Options.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_options_RECCES_Options_HH
#define INCLUDED_protocols_recces_options_RECCES_Options_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/recces/options/RECCES_Options.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace recces {
namespace options {

class RECCES_Options: public virtual basic::resource_manager::ResourceOptions {

public:

	//constructor
	RECCES_Options();

	RECCES_Options( RECCES_Options const & src );

	//destructor
	~RECCES_Options();

public:

	RECCES_OptionsOP clone() const;

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "RECCES_Options";}

	void
	initialize_from_command_line();


	void set_n_cycle( core::Size const & setting ){ n_cycle_ = setting; }
	core::Size n_cycle() const { return n_cycle_; }

	void set_n_dump( core::Size const & setting ){ n_dump_ = setting; }
	core::Size n_dump() const { return n_dump_; }

	void set_a_form_range( core::Real const & setting ){ a_form_range_ = setting; }
	core::Real a_form_range() const { return a_form_range_; }

	void set_base_pair_rmsd_cutoff( core::Real const & setting ){ base_pair_rmsd_cutoff_ = setting; }
	core::Real base_pair_rmsd_cutoff() const { return base_pair_rmsd_cutoff_; }

	void set_base_pair_rotation_mag( core::Real const & setting ){ base_pair_rotation_mag_ = setting; }
	core::Real base_pair_rotation_mag() const { return base_pair_rotation_mag_; }

	void set_base_pair_translation_mag( core::Real const & setting ){ base_pair_translation_mag_ = setting; }
	core::Real base_pair_translation_mag() const { return base_pair_translation_mag_; }

	void set_legacy_turner_mode( bool const & setting ){ legacy_turner_mode_ = setting; }
	bool legacy_turner_mode() const { return legacy_turner_mode_; }

	void set_dump_pdb( bool const & setting ){ dump_pdb_ = setting; }
	bool dump_pdb() const { return dump_pdb_; }

	void set_save_scores( bool const & setting ){ save_scores_ = setting; }
	bool save_scores() const { return save_scores_; }

	void set_block_stack( bool const & setting ){ block_stack_ = setting; }
	bool block_stack() const { return block_stack_; }

	void set_dump_freq( core::Size const & setting ){ dump_freq_ = setting; }
	core::Size dump_freq() const { return dump_freq_; }

	void set_dump_silent( bool const & setting ){ dump_silent_ = setting; }
	bool dump_silent() const { return dump_silent_; }

	void set_out_torsions( bool const & setting ){ out_torsions_ = setting; }
	bool out_torsions() const { return out_torsions_; }

	void set_setup_base_pair_constraints( bool const & setting ){ setup_base_pair_constraints_ = setting; }
	bool setup_base_pair_constraints() const { return setup_base_pair_constraints_; }

	void set_seq1( std::string const & setting ){ seq1_ = setting; }
	std::string seq1() const { return seq1_; }

	void set_seq2( std::string const & setting ){ seq2_ = setting; }
	std::string seq2() const { return seq2_; }

	void set_infile( std::string const & setting ){ infile_ = setting; }
	std::string infile() const { return infile_; }

	void set_xyz_file( std::string const & setting ){ xyz_file_ = setting; }
	std::string xyz_file() const { return xyz_file_; }

	void set_out_prefix( std::string const & setting ){ out_prefix_ = setting; }
	std::string out_prefix() const { return out_prefix_; }

	void set_temperatures( utility::vector1< core::Real > const & setting ){ temperatures_ = setting; }
	utility::vector1< core::Real > temperatures() const { return temperatures_; }

	void set_st_weights( utility::vector1< core::Real > const & setting ){ st_weights_ = setting; }
	utility::vector1< core::Real > st_weights() const { return st_weights_; }

	void set_sample_residues( utility::vector1< core::Size > const & setting ){ sample_residues_ = setting; }
	utility::vector1< core::Size > sample_residues() const { return sample_residues_; }

	void set_free_residues( utility::vector1< core::Size > const & setting ){ free_residues_ = setting; }
	utility::vector1< core::Size > free_residues() const { return free_residues_; }

	void set_loop_residues( utility::vector1< core::Size > const & setting ){ loop_residues_ = setting; }
	utility::vector1< core::Size > loop_residues() const { return loop_residues_; }

	void set_angle_range_bb( core::Real const & setting ){ angle_range_bb_ = setting; }
	core::Real angle_range_bb() const { return angle_range_bb_; }

	void set_angle_range_free_bb( core::Real const & setting ){ angle_range_free_bb_ = setting; }
	core::Real angle_range_free_bb() const { return angle_range_free_bb_; }

	void set_angle_range_chi( core::Real const & setting ){ angle_range_chi_ = setting; }
	core::Real angle_range_chi() const { return angle_range_chi_; }

	void set_angle_range_free_chi( core::Real const & setting ){ angle_range_free_chi_ = setting; }
	core::Real angle_range_free_chi() const { return angle_range_free_chi_; }

	void set_chi_stdev( core::Real const & setting ){ chi_stdev_ = setting; }
	core::Real chi_stdev() const { return chi_stdev_; }

	void set_bb_stdev( core::Real const & setting ){ bb_stdev_ = setting; }
	core::Real bb_stdev() const { return bb_stdev_; }

	void set_standard_bb_stdev( core::Real const & setting ){ standard_bb_stdev_ = setting; }
	core::Real standard_bb_stdev() const { return standard_bb_stdev_; }

	void set_silent_file( std::string const & setting ){ silent_file_ = setting; }
	std::string silent_file() const { return silent_file_; }

	void set_histogram_min( core::Real const & setting ){ histogram_min_ = setting; }
	core::Real histogram_min() const { return histogram_min_; }

	void set_histogram_max( core::Real const & setting ){ histogram_max_ = setting; }
	core::Real histogram_max() const { return histogram_max_; }

	void set_histogram_spacing( core::Real const & setting ){ histogram_spacing_ = setting; }
	core::Real histogram_spacing() const { return histogram_spacing_; }

	void set_thermal_sampler_mode( bool const & setting ){ thermal_sampler_mode_ = setting; }
	bool thermal_sampler_mode() const { return thermal_sampler_mode_; }

	void set_blank_score_terms( bool const & setting ){ blank_score_terms_ = setting; }
	bool blank_score_terms() const { return blank_score_terms_; }

	void set_skip_last_accept( bool const & setting ){ skip_last_accept_ = setting; }
	bool skip_last_accept() const { return skip_last_accept_; }

	void set_suppress_sampler_display( bool const & setting ){ suppress_sampler_display_ = setting; }
	bool suppress_sampler_display() const { return suppress_sampler_display_; }

	void set_prefix_each_output_pdb( bool const & setting ){ prefix_each_output_pdb_ = setting; }
	bool prefix_each_output_pdb() const { return prefix_each_output_pdb_; }

	void set_show_more_pose_scores( bool const & setting ){ show_more_pose_scores_ = setting; }
	bool show_more_pose_scores() const { return show_more_pose_scores_; }

	void set_output_simple_text_files( bool const & setting ){ output_simple_text_files_ = setting; }
	bool output_simple_text_files() const { return output_simple_text_files_; }

	void set_accept_no_op_moves( bool const & setting ){ accept_no_op_moves_ = setting; }
	bool accept_no_op_moves() const { return accept_no_op_moves_; }

	void set_sample_jump( bool const & setting ){ sample_jump_ = setting; }
	bool sample_jump() const { return sample_jump_; }

private:

	core::Size n_cycle_;
	core::Size n_dump_;

	core::Real a_form_range_;
	core::Real base_pair_rmsd_cutoff_;
	core::Real base_pair_rotation_mag_;
	core::Real base_pair_translation_mag_;

	bool legacy_turner_mode_;
	bool dump_pdb_;
	bool save_scores_;
	bool block_stack_;

	core::Size dump_freq_;
	bool dump_silent_;
	bool out_torsions_;
	bool setup_base_pair_constraints_;

	std::string seq1_;
	std::string seq2_;
	std::string infile_;
	std::string out_prefix_;
	std::string xyz_file_;

	utility::vector1< core::Real > temperatures_;
	utility::vector1< core::Real > st_weights_;

	core::Real histogram_min_;
	core::Real histogram_max_;
	core::Real histogram_spacing_;

	// for thermal_sampler
	core::Real angle_range_bb_;
	core::Real angle_range_free_bb_;
	core::Real angle_range_chi_;
	core::Real angle_range_free_chi_;
	core::Real chi_stdev_;
	core::Real bb_stdev_;
	core::Real standard_bb_stdev_;
	std::string silent_file_;
	// eventually hold these in FullModelInfo?
	utility::vector1< core::Size > sample_residues_;
	utility::vector1< core::Size > free_residues_;
	utility::vector1< core::Size > loop_residues_;

	bool thermal_sampler_mode_;

	// legacy match to thermal_sampler app
	bool blank_score_terms_;
	bool skip_last_accept_;
	bool suppress_sampler_display_;
	bool prefix_each_output_pdb_;
	bool show_more_pose_scores_;
	bool output_simple_text_files_;

	// legacy match to rb_recces app
	bool accept_no_op_moves_;
	bool sample_jump_;

};

} //options
} //recces
} //protocols

#endif
