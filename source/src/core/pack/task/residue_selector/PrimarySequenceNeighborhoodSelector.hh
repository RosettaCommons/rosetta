// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/PrimarySequenceNeighborhoodSelector.hh
/// @brief  The PrimarySequenceNeighborhoodSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_core_pack_task_residue_selector_PrimarySequenceNeighborhoodSelector_HH
#define INCLUDED_core_pack_task_residue_selector_PrimarySequenceNeighborhoodSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/PrimarySequenceNeighborhoodSelector.fwd.hh>

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

/// @brief The PrimarySequenceNeighborhoodSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are located near the given selected residues in primary sequence space
class PrimarySequenceNeighborhoodSelector : public ResidueSelector {
public:
	PrimarySequenceNeighborhoodSelector();

	PrimarySequenceNeighborhoodSelector(
		core::Size const lower_residues,
		core::Size const upper_residues,
		core::pack::task::residue_selector::ResidueSelectorCOP const selector );

	virtual ~PrimarySequenceNeighborhoodSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	virtual
	std::string
	get_name() const;

	static std::string class_name();

	//unit-specific
	void set_lower_residues( core::Size const nres );
	void set_upper_residues( core::Size const nres );
	void set_selector( core::pack::task::residue_selector::ResidueSelectorCOP const selector );

private: // data members
	core::Size lower_residues_;
	core::Size upper_residues_;
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
