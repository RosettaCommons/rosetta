// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResiduePDBInfoHasLabelSelector.hh
/// @brief  The ResiduePDBInfoHasLabelSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_core_pack_task_residue_selector_ResiduePDBInfoHasLabelSelector_HH
#define INCLUDED_core_pack_task_residue_selector_ResiduePDBInfoHasLabelSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/ResiduePDBInfoHasLabelSelector.fwd.hh>

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

/// @brief The ResiduePDBInfoHasLabelSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions with the given label in pose PDB info for that residue.
class ResiduePDBInfoHasLabelSelector : public ResidueSelector {
public:
	// derived from base class
	ResiduePDBInfoHasLabelSelector();
	ResiduePDBInfoHasLabelSelector( std::string const & label_str );
	virtual ~ResiduePDBInfoHasLabelSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();

	void set_label( std::string const & label_str );

private: // data members
	std::string label_;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
