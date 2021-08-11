// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMover.hh
/// @brief reports information about disulfides
/// @author Frank Teets (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_movers_OmnibusDisulfideAnalysisLabelerMover_HH
#define INCLUDED_protocols_pose_sewing_movers_OmnibusDisulfideAnalysisLabelerMover_HH

// Unit headers
#include <protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <map>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace pose_sewing {
namespace movers {

///@brief reports information about disulfides
class OmnibusDisulfideAnalysisLabelerMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	OmnibusDisulfideAnalysisLabelerMover();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~OmnibusDisulfideAnalysisLabelerMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//OmnibusDisulfideAnalysisLabelerMover & operator=( OmnibusDisulfideAnalysisLabelerMover const & src );

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

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns Frank Teets as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private: // methods

private: // data
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Size min_helix_length_ = 6;
	core::Real cutoff_ = 2;
	core::Real crit_dist_ = 4;
	core::Real crit_score_ = -0.1;
	core::Real qualifying_score_ = -1.0;
	bool use_motifs_ = false;

};

std::ostream &
operator<<( std::ostream & os, OmnibusDisulfideAnalysisLabelerMover const & mover );

} //movers
} //pose_sewing
} //protocols

#endif //protocols_pose_sewing_movers_OmnibusDisulfideAnalysisLabelerMover_HH
