// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/SewAnythingAddMover.hh
/// @brief perform SEWING-like chimerization-based addition of supersecondary structural segments into previously virtual residues
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_movers_SewAnythingAddMover_HH
#define INCLUDED_protocols_pose_sewing_movers_SewAnythingAddMover_HH

// Unit headers
#include <protocols/pose_sewing/movers/SewAnythingAddMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh>
#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Basic/Utility headers
#include <numeric/xyzVector.hh>
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility>

namespace protocols {
namespace pose_sewing {
namespace movers {

///@brief Settings for SewAnything.  This is stored for each SS.

///@brief perform SEWING-like chimerization-based addition of supersecondary structural segments into previously virtual residues
class SewAnythingAddMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SewAnythingAddMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SewAnythingAddMover( SewAnythingAddMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SewAnythingAddMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	std::pair<core::Size,core::Size>
	find_sewable_region(core::pose::Pose& pose);

	char
	find_mod_terminus(core::pose::Pose& pose,std::pair<core::Size,core::Size>);

	bool
	align(
		core::pose::Pose const & stationary_pose,
		core::Size stationary_residue,
		core::pose::Pose & mobile_pose,
		core::Size mobile_residue,
		core::Size window_width,
		core::Real tolerance,
		bool check_alignment = true) const;

	bool
	clashes(
		core::pose::Pose const & nterm_pose,
		core::Size first_nterm_res,
		core::Size last_nterm_res,
		core::pose::Pose const & cterm_pose,
		core::Size first_cterm_res,
		core::Size last_cterm_res,
		core::Real clash_radius,
		bool check_all_backbone) const;


	data_storage::TerminalDSSPSortedPoseVectorOP
	get_pose_vector() const;

	core::Real
	get_clash_radius() const;

	//std::string
	//get_pdb_file_name() const;

	utility::vector1< std::string >
	get_segment_file_paths() const;

	//bool
	//get_read_segments_from_segment_file() const;

	bool
	get_trim_terminal_loops() const;

	std::string
	get_permissible_segment_ends() const;

	void
	set_permissible_segment_ends( std::string );

	std::string
	get_permissible_termini() const;

	core::Size
	get_max_attempts() const;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//SewAnythingAddMover & operator=( SewAnythingAddMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	std::pair< std::string, std::string>
	segment_vector_dssp(char mod_terminus, char mod_terminus_dssp, std::string const & permissible_segment_ends) const;

public:

	void
	set_segment_file_paths(utility::vector1<std::string> segment_file_paths);

private: // methods

private: // data

	utility::vector1< std::string > segment_file_paths_;

	core::Real clash_radius_ = 4.0;
	data_storage::TerminalDSSPSortedPoseVectorOP pose_vector_;
	std::string permissible_segment_ends_;
	core::Size max_attempts_ = 1000;
	core::Size max_filter_attempts_ = 1000;
	bool trim_terminal_loops_ = true;
	bool check_all_backbone_atoms_ = false;
	std::string permissible_termini_ = "NC";

	bool fail_on_no_match_ = false;

	//Settings for each SS if set.  "A" is for ALL and indicates default.
	// Could add Enums for SS chars, but that seems like too much for now.
	utility::vector1< protocols::filters::FilterCOP > all_seg_filters_;
	utility::vector1< protocols::filters::FilterCOP > all_post_filters_;
	std::map< std::string, utility::vector1< protocols::filters::FilterCOP >> seg_filters_;

	core::simple_metrics::RealMetricCOP sort_metric_ = nullptr;

	bool enable_clash_check_v2_ = true;
	core::pose::PoseCOP ht_ref_pose_;
	bool even_sampling_ = false;
	bool positive_scores_are_better_ = false;

	bool use_absolute_sizes_ = true;
	bool use_relative_sizes_ = false;
	bool allow_DSSP_insertion_ = true;

	core::Size hashable_element_max_size_ = 1000;
	core::Size hashable_element_min_size_ = 0;

	core::Size hashable_element_relative_max_size_ = 1000;
	core::Size hashable_element_relative_min_size_ = 1000;

	core::Size window_width_ = 0;
	core::Real alignment_max_distance_ = 10.0;
	core::Size max_cter_len_ = 100;
	core::Size max_nter_len_ = 100;
};

std::ostream &
operator<<( std::ostream & os, SewAnythingAddMover const & mover );

} //protocols
} //pose_sewing
} //movers

#endif //protocols_pose_sewing_movers_SewAnythingAddMover_HH
