// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/NotResidueSelector.hh
/// @brief  The NotResidueSelector negates the logic of its loaded ResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

#ifndef INCLUDED_core_select_residue_selector_NotResidueSelector_HH
#define INCLUDED_core_select_residue_selector_NotResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/NotResidueSelector.fwd.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The NotResidueSelector negates the input of one loaded ResidueSelector, i.e., it is a logical NOT -
/// it selects all unselected residues and deselects the selected ones.  True becomes false, false becomes true.
/// The ResidueSelector to be negated can be pulled in through RosettaScipt using the selector option, subtags for
/// ResidueSelectors known to the ResidueSelectorFactory or programmatically using set_residue_selector.
/// Note that since most ResidueSelectors clear the input ResidueSubset, NOT can be thought of as simply selecting
/// the opposite of the passed in selector.
class NotResidueSelector : public ResidueSelector {
public:
	// derived from base class
	NotResidueSelector();

	/// @brief Copy constructor
	///
	NotResidueSelector( NotResidueSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	NotResidueSelector( ResidueSelectorCOP selector );
	~NotResidueSelector() override;

	ResidueSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief sets a ResidueSelector
	*/
	void set_residue_selector(ResidueSelectorCOP selector);

	std::string
	debug_string() const override;

public: //Functions needed for the citation manager

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private: // data members
	ResidueSelectorCOP selector_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_NotResidueSelector )
#endif // SERIALIZATION


#endif
