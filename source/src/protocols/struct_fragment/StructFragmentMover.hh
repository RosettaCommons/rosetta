// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   StructFragmentMover.hh
/// @brief  Creating fragments from supplied structures.
/// @author andreas scheck (graziano.bud@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_struct_fragment_StructFragmentMover_hh
#define INCLUDED_protocols_struct_fragment_StructFragmentMover_hh


// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/struct_fragment/StructFragmentMover.fwd.hh>

// Package Headers
#include <core/fragment/FragSet.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// Project Headers
#include <basic/datacache/DataMap.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace struct_fragment {

class StructFragmentMover : public  protocols::moves::Mover {
public:
	/// @brief Empty Constructor
	StructFragmentMover();

	/// @brief Destructor
	~StructFragmentMover();

	/// @brief Apply Mover
	void apply( core::pose::Pose & pose );

	// Functions necessary for RosettaScripts
	/// @brief Create clone of this mover
	inline protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP( new StructFragmentMover( *this ) ); };

	/// @brief Create a fresh instance of this mover
	inline protocols::moves::MoverOP fresh_instance() const{ return protocols::moves::MoverOP( new StructFragmentMover() ); };

	/// @brief Parse RosettaScripts options for this mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
);

static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// getters

	/// @brief Get mover name
	inline std::string get_name() const { return mover_name(); };

	/// @brief Return mover name
	inline static std::string mover_name() { return "StructFragmentMover"; };

	/// @brief Get the number of fragment candidates per position
	core::Size get_n_candidates() const { return structPicker_->n_candidates_; };

	/// @brief Get the number of fragments per position
	core::Size get_n_frags() const { return structPicker_->n_frags_; };

	/// @brief Get the confidence for phi and psi angles
	core::Real get_loop_angle_conf() const { return loop_angle_conf_; };

	/// @brief Get the size of small fragments
	core::Size get_small_frag_size() const { return small_frag_size_; };

	/// @brief Get the size of large fragments
	core::Size get_large_frag_size() const { return large_frag_size_; };

	/// @brief Get the path to the small fragment file
	std::string get_small_frag_file() const { return small_frag_file_; };

	/// @brief Get the parth to the large fragment file
	std::string get_large_frage_file() const { return large_frag_file_; };

	/// @brief Get the prefix for the fragment file output
	std::string get_prefix() const { return structPicker_->prefix_; };

	/// @brief Check if fragment files should be outputted
	bool get_output_frag_files() const { return output_frag_files_; };

	/// @brief Check if small fragments are being stolen
	bool get_steal_small_frags() const { return steal_small_frags_; };

	/// @brief Check if large fragments are being stolen
	bool get_steal_large_frags() const { return steal_large_frags_; };

	/// @brief Get the path to the fragment weight file
	std::string get_frag_weight_file() const { return frag_weight_file_; };

	/// @brief Get the path to a sequence profile
	std::string get_sequence_profile() const { return sequence_profile_; };

	/// @brief Get the path to the vall file
	std::string get_vall_file() const { return vall_file_; };


	// setters

	/// @brief Set the number of fragment candidates per position
	void set_n_candidates( core::Size n_candidates ) { structPicker_->n_candidates_ = n_candidates; };

	/// @brief Set the number of fragments per position
	void set_n_frags( core::Size n_frags ) { structPicker_->n_frags_ = n_frags; };

	/// @brief Set the confidence for phi and psi angles
	void set_loop_angle_conf( core::Real loop_angle_conf ) { loop_angle_conf_ = loop_angle_conf; };

	/// @brief Set the size of small fragments
	void set_small_frag_size( core::Size small_frag_size ) { small_frag_size_ = small_frag_size; };

	/// @brief Set the size of large fragments
	void set_large_frag_size( core::Size large_frag_size ) { large_frag_size_ = large_frag_size; };

	/// @brief Set the path to the small fragment file
	void set_small_frag_file( std::string small_frag_file ) { small_frag_file_ = small_frag_file; };

	/// @brief Set the parth to the large fragment file
	void set_large_frag_file( std::string large_frag_file ) { large_frag_file_ = large_frag_file; };

	/// @brief Set the prefix for the fragment file output
	void set_prefix( std::string prefix ) { structPicker_->prefix_ = prefix; };

	/// @brief Set whether fragment files should be outputted
	void set_output_frag_files( bool output_frag_files ) { output_frag_files_ = output_frag_files; };

	/// @brief Set whether small fragments are being stolen
	void set_steal_small_frags( bool steal_small_frags ) { steal_small_frags_ = steal_small_frags; };

	/// @brief Set whether large fragments are being stolen
	void set_steal_large_frags( bool steal_large_frags ) { steal_large_frags_ = steal_large_frags; };

	/// @brief Set the path to the fragment weight file
	void set_frag_weight_file( std::string frag_weight_file ) { frag_weight_file_ = frag_weight_file; };

	/// @brief Set the path to a sequence profile
	void set_sequence_profile( std::string sequence_profile ) { sequence_profile_ = sequence_profile; };

	/// @brief Set the path to the vall file
	void set_vall_file( std::string vall_file ) { vall_file_ = vall_file; };



private:
	protocols::frag_picker::FragmentPickerOP structPicker_;
	core::Real loop_angle_conf_;			// confidence of phi and psi angles
	core::Size small_frag_size_;			// size of small fragments
	core::Size large_frag_size_;			// size of large fragments
	std::string small_frag_file_;			// path to small fragment file
	std::string large_frag_file_;			// path to large fragment file
	bool output_frag_files_;
	bool steal_small_frags_;
	bool steal_large_frags_;
	bool changed_frags_;
	std::string frag_weight_file_;		// path to fragment weight file
	std::string sequence_profile_;		// path to sequence profile
	std::string vall_file_;						// path to vall file
	core::fragment::FragSetOP smallF_;
	core::fragment::FragSetOP largeF_;


	// default values

	/// @brief Default value for the confidence of phi and psi angles
	inline static core::Real default_value_for_loop_angle_conf() { return 0.8; };

	/// @brief Default value for the small fragment size
	inline static core::Size default_value_for_small_frag_size() { return 3; };

	/// @brief Default value for the large fragment size
	inline static core::Size default_value_for_large_frag_size() { return 9; };

	/// @brief Default value for outputting fragment files
	inline static bool default_value_for_output_frag_files() { return false; };

	/// @brief Default value for stealing small fragments
	inline static bool default_value_for_steal_small_frags() { return false; };

	/// @brief Default value for stealing large fragments
	inline static bool default_value_for_steal_large_frags() { return false; };

	/// @brief Default value for empty file
	inline static std::string default_value_for_empty_file() { return std::string(); };

	/// @brief Default value for fragment candidates per position
	inline static core::Size default_value_for_n_candidates() { return 1000; };

	/// @brief Default value for the number of fragments per position
	inline static core::Size default_value_for_n_frags() { return 200; };

	/// @brief Default prefix for fragment output files
	inline static std::string default_value_for_prefix() { return "structFrag"; };


	/// @brief evaluates which process to run
	core::Size evaluate_job();
	/// @brief Collect fragments for the specified input structure
	utility::vector1< core::fragment::FragSetOP > get_fragments( core::pose::Pose const & pose );

	/// @brief Create fragment picker for the specified input structure
	protocols::frag_picker::FragmentPickerOP make_fragment_picker( core::pose::Pose pose, std::string vall_file_name );

	core::fragment::FragSetOP steal_fragments_by_size( core::pose::Pose const & pose, core::fragment::FragSetOP fset, core::Size size );

	utility::vector1< core::fragment::FragSetOP > steal_fragments( core::pose::Pose const & pose, utility::vector1< core::fragment::FragSetOP > fsets );

	void collector_to_picker( protocols::frag_picker::FragmentPickerOP pickIt, core::Size size, core::Size seqlength );

	core::fragment::FragSetOP get_fragset ( protocols::frag_picker::FragmentPickerOP pickIt, core::Size position, core::Size size );

	core::fragment::FragSetOP shift_fragset( core::fragment::FragSetOP fset, core::id::SequenceMapping const & seqmap );

	/// @brief Read existing fragment file
	core::fragment::FragSetOP read_frag_file(std::string frag_file );

	/// @brief Read small and large fragment file and combine the fragments
	utility::vector1< core::fragment::FragSetOP > read_frag_files(std::string small_frag_file, std::string large_frag_file );

}; // end class StructFragmentMover

} // end namespace struct_fragment
} // end namespace protocols

#endif
