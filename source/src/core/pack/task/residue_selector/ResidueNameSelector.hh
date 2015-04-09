// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueNameSelector.hh
/// @brief  The ResidueNameSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueNameSelector_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueNameSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/ResidueNameSelector.fwd.hh>

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

/// @brief The ResidueNameSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which match the given residue index. The index is read as comma-separated
/// list of either Rosetta indices (e.g. 10) or PDB numbers (e.g. 10A, residue 10 of chain A). Detection
/// and mapping from PDB to Rosetta residue numbers is done internally.
class ResidueNameSelector : public ResidueSelector {
public:
	// derived from base class
	ResidueNameSelector();
	ResidueNameSelector( std::string const & res_name_str );
	virtual ~ResidueNameSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();


	//unit-specific
	/**
	* @brief sets the comma-separated string of residue names to be selected
	*/
	void set_residue_names( std::string const & res_name_str );

private: // data members
	std::string res_name_str_;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
