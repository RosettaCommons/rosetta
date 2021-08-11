// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/BlockwiseAnalysisMover.hh
/// @brief analyzes pairs of interacting elements
/// @author Frank Teets (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_movers_BlockwiseAnalysisMover_HH
#define INCLUDED_protocols_pose_sewing_movers_BlockwiseAnalysisMover_HH

// Unit headers
#include <protocols/pose_sewing/movers/BlockwiseAnalysisMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/selection.hh>
#include <core/types.hh>

namespace protocols {
namespace pose_sewing {
namespace movers {

///@brief analyzes pairs of interacting elements
class BlockwiseAnalysisMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	BlockwiseAnalysisMover();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~BlockwiseAnalysisMover() override;


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

	//BlockwiseAnalysisMover & operator=( BlockwiseAnalysisMover const & src );

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

	void
	provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;
private: // methods

private: // data
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Real crit_dist_ = 4.0;

};

std::ostream &
operator<<( std::ostream & os, BlockwiseAnalysisMover const & mover );

} //movers
} //pose_sewing
} //protocols

#endif //protocols_pose_sewing_movers_BlockwiseAnalysisMover_HH
