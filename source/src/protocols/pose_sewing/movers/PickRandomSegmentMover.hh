// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/PickRandomSegmentMover.hh
/// @brief replaces pose with a random starting segment from a segment file
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_movers_PickRandomSegmentMover_HH
#define INCLUDED_protocols_pose_sewing_movers_PickRandomSegmentMover_HH

// Unit headers
#include <protocols/pose_sewing/movers/PickRandomSegmentMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <numeric/xyzVector.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh>
#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace pose_sewing {
namespace movers {

struct SeedSettings {

	core::Size min_terminal_length_ = 0;
	core::Size max_terminal_length_ = 0;
};

///@brief replaces pose with a random starting segment from a segment file
class PickRandomSegmentMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PickRandomSegmentMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PickRandomSegmentMover( PickRandomSegmentMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PickRandomSegmentMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	data_storage::TerminalDSSPSortedPoseVectorOP
	get_pose_vector() const;

	utility::vector1< std::string >
	get_segment_file_paths() const;

	void
	set_segment_file_paths(utility::vector1< std::string > paths);

public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//PickRandomSegmentMover & operator=( PickRandomSegmentMover const & src );

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

public:

	///@brief Get all SS-specific settings
	std::map< char, SeedSettings >
	get_ss_settings() const;

	///@brief Create and return the SS-tunable settings.
	SeedSettings
	create_ss_defaults() const;

	///@brief Get the defaults of seed settings.
	SeedSettings
	get_ss_defaults() const;

	///@brief Set the configurable options as defaults.
	void
	set_ss_defaults();

private: // methods

private: // data

	utility::vector1< std::string > segment_file_paths_;
	data_storage::TerminalDSSPSortedPoseVectorOP pose_vector_;

	bool match_terminal_sheet_lengths_ = false;
	core::Size terminal_sheet_length_max_delta_=1;

	//Settings for each SS if set.  "A" is for ALL and indicates default.
	std::map < char, SeedSettings > ss_settings_;

	utility::vector1< protocols::filters::FilterCOP > filters_;
	std::map< std::string, utility::vector1< protocols::filters::FilterCOP >> seg_filters_;

	bool even_sampling_ = true;

};

std::ostream &
operator<<( std::ostream & os, PickRandomSegmentMover const & mover );

} //movers
} //pose_sewing
} //protocols

#endif //protocols_pose_sewing_movers_PickRandomSegmentMover_HH
