// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/OrResidueSelector.hh
/// @brief  The OrResidueSelector combines logic from multiple ResidueSelectors
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_OrResidueSelector_HH
#define INCLUDED_core_pack_task_residue_selector_OrResidueSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/OrResidueSelector.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreator.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

/// @brief The OrResidueSelector combines the output of multiple ResidueSelectors using OR
/// logic, i.e., residues selected by ANY of the contained ResidueSelectors will be selected.
/// ResidueSelecters can be pulled in from a DataMap, from subtags (for ResidueSelectors
/// known to the ResidueSelectorFactory) or programmatically through %add_residue_selector.
class OrResidueSelector : public ResidueSelector {
public:
	// derived from base class
	OrResidueSelector();
	OrResidueSelector( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2 );
	virtual ~OrResidueSelector();

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
	void add_residue_selector(ResidueSelectorCOP selector);

	Size num_selectors() const;

	/**
	* @brief applies newSubset to existingSubset and thereby modifies the latter
	*/
	void apply_or_to_subset(ResidueSubset const & newSubset, ResidueSubset & existingSubset) const;

private: // data members

	std::list< ResidueSelectorCOP > selectors_;

};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
