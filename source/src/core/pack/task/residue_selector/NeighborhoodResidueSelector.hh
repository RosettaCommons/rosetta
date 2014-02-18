// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/NeighborhoodResidueSelector.hh
/// @brief  The NeighborhoodResidueSelector selects residues in a given proximity of set focus residues
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

#ifndef INCLUDED_core_pack_task_residue_selector_NeighborhoodResidueSelector_HH
#define INCLUDED_core_pack_task_residue_selector_NeighborhoodResidueSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/NeighborhoodResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

/// @brief The NeighborhoodResidueSelector selects residues neighboring a defined set of residues
/// (the focus). The focus residue set can be obtained from another ResidueSelector, from a
/// std::set of residue positions or from a string specifying residue positions.
class NeighborhoodResidueSelector : public ResidueSelector {
public:
	// derived from base class
	NeighborhoodResidueSelector();
	NeighborhoodResidueSelector( std::set<core::Size> const & focus, Real distance );
	virtual ~NeighborhoodResidueSelector();

	virtual void apply( core::pose::Pose const & pose, ResidueSubset & subset ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();


	//unit-specific
	/**
	* @brief adds a ResidueSelector
	*/
	void set_focus( std::set<Size> const &focus );
	void set_focus( std::string const & focus_str );
	void set_focus_selector( ResidueSelectorCOP rs );
	void set_distance( Real distance );

private: //functions
	void get_focus( core::pose::Pose const &, ResidueSubset &, std::set< Size > &) const;

private: // data members
	// data in focus and focus_string will be stitched together.
	// focus_str_ can only be parsed when pose is available (PDB mappings)
	// think of either-or behavior also between set and string
	std::set< Size > focus_;
	std::string focus_str_;
	Real distance_;

	// focus residues may be selected directly be another ResidueSelector
	ResidueSelectorCOP focus_selector_;

	// has any focus been set
	bool focus_set_;

	// is the focus selector the last source of focus that has been set
	bool use_focus_selector_;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
