// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/AndResidueSelector.hh
/// @brief  The AndResidueSelector combines logic from multiple ResidueSelectors
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_AndResidueSelector_HH
#define INCLUDED_core_select_residue_selector_AndResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/AndResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The AndResidueSelector combines the output of multiple ResidueSelectors using AND
/// logic, i.e., only residues selected by ALL contained ResidueSelectors will be selected.
/// ResidueSelecters can be pulled in from a DataMap, from subtags (for ResidueSelectors
/// known to the ResidueSelectorFactory) or programmatically through %add_residue_selector.
class AndResidueSelector : public ResidueSelector {
public:
	// derived from base class
	AndResidueSelector();

	/// @brief Copy constructor
	///
	AndResidueSelector( AndResidueSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	AndResidueSelector( ResidueSelectorCOP selector1);
	AndResidueSelector( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2 );

	virtual ~AndResidueSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief adds a ResidueSelector
	*/
	void add_residue_selector(ResidueSelectorCOP selector);

	Size num_selectors() const;

	/**
	* @brief applies newSubset to existingSubset and thereby modifies the latter
	*/
	void apply_and_to_subset(ResidueSubset const & newSubset, ResidueSubset & existingSubset) const;

private: // data members

	std::list< ResidueSelectorCOP > selectors_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_AndResidueSelector )
#endif // SERIALIZATION


#endif
