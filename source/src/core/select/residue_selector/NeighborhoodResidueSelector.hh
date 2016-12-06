// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/NeighborhoodResidueSelector.hh
/// @brief  The NeighborhoodResidueSelector selects residues in a given proximity of set focus residues.
///  Clears the Passed ResidueSubset.
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - 10 A neighbor graph, simplification, ResidueSubset as focus, clean up, etc.


#ifndef INCLUDED_core_select_residue_selector_NeighborhoodResidueSelector_HH
#define INCLUDED_core_select_residue_selector_NeighborhoodResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/NeighborhoodResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The NeighborhoodResidueSelector selects residues neighboring a defined set of residues
/// (the focus). The focus residue set can be obtained from another ResidueSelector, from a
/// set of positions.  Focus is included in subset by default.  Use include_focus_in_subset to change this!
///
/// Note:  Uses the 10 A neighbor graph by default.  If neighbor distance is great than this, we use slow-style double for loop.
///
class NeighborhoodResidueSelector : public ResidueSelector {
public:
	// derived from base class
	NeighborhoodResidueSelector();

	/// @brief Copy constructor
	///
	NeighborhoodResidueSelector( NeighborhoodResidueSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	NeighborhoodResidueSelector( utility::vector1< bool > const & focus, Real distance, bool include_focus_in_subset = true);

	NeighborhoodResidueSelector( ResidueSelectorCOP selector, Real distance, bool include_focus_in_subset = true);

	virtual ~NeighborhoodResidueSelector();

public:

	/// @brief Set the focus, which is the residues for which we will be getting neighbors of.
	void
	set_focus( std::string const & focus_str );

	/// @brief Set the focus, which is the residues for which we will be getting neighbors of.
	void
	set_focus( utility::vector1< bool > const & focus);

	/// @brief Set a Residue Selector for the focus
	void
	set_focus_selector( ResidueSelectorCOP rs );

	/// @brief Set the distance we will be measuring to get neighbors
	void
	set_distance( Real distance );

	///@brief Setting to include the fucus in the resulting subset or not. Default is TRUE
	void
	set_include_focus_in_subset( bool include_focus);

public:


	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string
	class_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );




private:

	void
	get_focus( core::pose::Pose const &, ResidueSubset &, utility::vector1< bool > &) const;

	void
	set_defaults();

	void
	clear_focus();

private:

	// data in focus and focus_string will be stitched together.
	// focus_str_ can only be parsed when pose is available (PDB mappings)
	// think of either-or behavior also between set and string
	utility::vector1< bool > focus_;
	std::string focus_str_;
	Real distance_;

	// focus residues may be selected directly be another ResidueSelector
	ResidueSelectorCOP focus_selector_;

	bool include_focus_in_subset_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_NeighborhoodResidueSelector )
#endif // SERIALIZATION


#endif
