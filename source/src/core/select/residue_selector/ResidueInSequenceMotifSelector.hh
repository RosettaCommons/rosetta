// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueInSequenceMotifSelector.hh
/// @brief  Select residues by motif search
/// @author raemisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueInSequenceMotifSelector_HH
#define INCLUDED_core_select_residue_selector_ResidueInSequenceMotifSelector_HH

// Unit headers
#include <core/select/residue_selector/ResidueInSequenceMotifSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>
#include <regex>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief Select residues by motif search
/// @detail This selector takes a regular expression and a position and then finds all matches of the regular expression
/// and selects the n-th residue in each match, where 'n' is the specified position.
class ResidueInSequenceMotifSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	ResidueInSequenceMotifSelector();
	/// @brief Constructor for derived classes.
	ResidueInSequenceMotifSelector(std::string regex ,core::Size position_in_motif) : regex_(regex), position_in_motif_(position_in_motif) {};

public:

	/// @brief Destructor.
	~ResidueInSequenceMotifSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:

	// Regular expression describing the sequence motif
	std::string regex_;
	// Position of residue withing the motif to be selected
	core::Size position_in_motif_;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_ResidueInSequenceMotifSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_ResidueInSequenceMotifSelector_HH
