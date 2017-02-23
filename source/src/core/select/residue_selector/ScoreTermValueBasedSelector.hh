// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ScoreTermValueBasedSelector.hh
/// @brief  Selects residues based on per residue score of the given score_type.
/// Residues with scores equal to or below the given threshold are selected.
/// @author Gerard Daniel (gerardda@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_ScoreTermValueBasedSelector_hh
#define INCLUDED_core_select_residue_selector_ScoreTermValueBasedSelector_hh

// Unit headers
#include <core/select/residue_selector/ScoreTermValueBasedSelector.fwd.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/tag/XMLSchemaGeneration.hh>


namespace core {
namespace select {
namespace residue_selector {

class ScoreTermValueBasedSelector : public core::select::residue_selector::ResidueSelector {

public:

	/// @brief Constructor.
	ScoreTermValueBasedSelector();

	//ScoreTermValueBasedSelector( core::pose::Pose const & pose );

	/// @breif Destructor
	virtual ~ScoreTermValueBasedSelector();

	/// @brief Clone operator.
	/// @details Creates a copy of the object and return a pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	/// @brief Return a ResidueSubset indicating a selection of Residues from the input Pose.
	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief Initialize any data members of this instance from an input tag and a DataMap object
	virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datacache );

	virtual std::string get_name() const;

	static std::string class_name();

	/// @brief Define the structure of the XML file for this ResidueSelector
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & );

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::scoring::ScoreType score_type_;
	ResidueSelectorCOP input_residues_selector_;
	std::string residue_nums_string_;
	core::Real lower_threshold_;
	core::Real upper_threshold_;

};

}
}
}


#endif /* INCLUDED_core_select_residue_selector_ScoreTermValueBasedSelector_hh */
