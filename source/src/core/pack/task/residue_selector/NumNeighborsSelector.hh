// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/NumNeighborsSelector.hh
/// @brief  The NumNeighborsSelector identifies residues that have at least X neighbors within a distance Y
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_NumNeighborsSelector_HH
#define INCLUDED_core_pack_task_residue_selector_NumNeighborsSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/NumNeighborsSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class NumNeighborsSelector : public ResidueSelector {
public:
	// derived from base class
	NumNeighborsSelector();
	NumNeighborsSelector( Size threshold, core::Real distance_cutoff );
	//// Undefined, commenting out to fix PyRosetta build  NumNeighborsSelector( bool count_water, Size threshold, core::Real distance_cutoff );
	virtual ~NumNeighborsSelector();

	virtual void apply( core::pose::Pose const & pose, ResidueSubset & subset ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();

	bool count_water() const;
	Size threshold() const;
	Real distance_cutoff() const;
	void count_water( bool setting );
	void threshold( Size setting );
	void distance_cutoff( Size setting );

private:
	bool count_water_;
	Size threshold_;
	Real distance_cutoff_;

};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
