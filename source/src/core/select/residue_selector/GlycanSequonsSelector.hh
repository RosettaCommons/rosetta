// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/GlycanSequonsSelector.hh
/// @brief  Find glycosylation sequons in pose
/// @author raemisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_select_residue_selector_GlycanSequonsSelector_HH
#define INCLUDED_core_select_residue_selector_GlycanSequonsSelector_HH

// Unit headers
#include <core/select/residue_selector/GlycanSequonsSelector.fwd.hh>
#include <core/select/residue_selector/ResidueInSequenceMotifSelector.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Find glycosylation sequons in pose
/// @detail Glycosylation sites are typically recognized by enzymes by means of a recognition sequence, or 'motif'.
/// This ResidueSelector selects all residues that can be glycosylated using one of a few sequence motifs that are
/// pre-defined and can be selected.
/// For selecting residues in other motifs, use the ResidueInSequenceMotifSelctor, which takes regular expressions as input.
class GlycanSequonsSelector : public core::select::residue_selector::ResidueInSequenceMotifSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	GlycanSequonsSelector();

public:

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual ResidueSelectorOP clone() const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	/// @brief Get the mover class name.
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	// Known glycosylation motifs that can be chosen
	bool NxST_ = true;
	bool NxC_ = false;
	bool WxxW_ = false;
	bool WSTxC_ = false;

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
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_GlycanSequonsSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_GlycanSequonsSelector_HH
