// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/SecondaryStructureSelector.hh
/// @brief  The SecondaryStructureSelector selects residues based on their secondary structure
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_pack_task_residue_selector_SecondaryStructureSelector_HH
#define INCLUDED_core_pack_task_residue_selector_SecondaryStructureSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/SecondaryStructureSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <set>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class SecondaryStructureSelector : public ResidueSelector {
public:
	SecondaryStructureSelector();
	SecondaryStructureSelector( std::string const & selected );

	virtual ~SecondaryStructureSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & datamap );

	virtual std::string get_name() const;

	static std::string class_name();

public:
	// mutators

	/// @brief sets number of residues around the given SS elements to select (default=0)
	void set_overlap( core::Size const overlapval );
	/// @brief sets ss characters to select -- must be set prior to apply()
	void set_selected_ss( std::string const & selected );
	/// @brief if true, one-residue terminal "loops" will be included (default=false)
	void set_include_terminal_loops( bool const inc_term );

private:
	void add_overlap(
			ResidueSubset & matching_ss,
			pose::Pose const & pose,
			std::string const & ss ) const;

	bool check_ss( std::string const & ss ) const;

private:
	core::Size overlap_;
	bool include_terminal_loops_;
	std::set< char > selected_ss_;
};

typedef std::pair< Size, Size > Interval;
typedef utility::vector1< Interval > IntervalVec;
IntervalVec subset_to_intervals( ResidueSubset const & subset );

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
