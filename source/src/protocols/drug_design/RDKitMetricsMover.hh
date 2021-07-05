// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/RDKitMetricsMover.hh
/// @brief Add metrics from RDKit to the pose extra data.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_RDKitMetricsMover_hh
#define INCLUDED_protocols_drug_design_RDKitMetricsMover_hh

// Unit header
#include <protocols/drug_design/RDKitMetricsMover.fwd.hh>

#include <protocols/moves/Mover.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Project Headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// External headers

// C/C++ headers
#include <string>

namespace protocols {
namespace drug_design {

class RDKitMetricsMover : public protocols::moves::Mover {
public:
	/// @brief default constructor
	RDKitMetricsMover();

	/// @brief destructor
	~RDKitMetricsMover() override;

	/// @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	void chains(utility::vector1< std::string > const & setting ) {
		chains_ = setting;
	}

	utility::vector1< std::string > const & chains() const { return chains_; }

	/// @brief apply RDKitMetricsMover
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	/// @brief parse xml file
	void
	parse_my_tag( TagCOP tag, basic::datacache::DataMap & data ) override;

	void
	add_scores_for_chain( Pose & pose, std::string const & chain ) const;

	void
	add_scores_for_residue( Pose & pose, core::Size resid ) const;

	std::string
	get_tag( Pose & pose, core::Size resid ) const;

	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // Data

	utility::vector1< std::string > chains_;
};

} // namespace drug_design
} // namespace protocols

#endif
