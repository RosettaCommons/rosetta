// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ByMembraneDepthSelector.hh
/// @brief  Select residues according to their depth in the membrane bilayer relative to the center
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_ByMembraneDepthSelector_HH
#define INCLUDED_core_select_residue_selector_ByMembraneDepthSelector_HH

// Unit headers
#include <core/select/residue_selector/ByMembraneDepthSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Select residues according to their depth in the membrane bilayer relative to the center
class ByMembraneDepthSelector : public core::select::residue_selector::ResidueSelector {

public: // constructors & methods

	// derived from base class
	ByMembraneDepthSelector();
	virtual ~ByMembraneDepthSelector();

	virtual core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;
	
	virtual ResidueSelectorOP clone() const;

	static std::string class_name();
	static void provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd );

	// Access & Mdify data
	core::Real minimum_depth() const;
	void minimum_depth( core::Real min_depth );
	
	core::Real maximum_depth() const;
	void maximum_depth( core::Real max_depth );
	
private: // methods

	/// @brief Determine whether the given residue falls within the user specified depth range
	bool residue_within_depth( core::pose::Pose const & pose, core::Size rsdnum ) const;
	
	/// @brief Check that the user specified bounds are within a reasonable range
	void check_valid_bounds() const;

private:

	// Mimumum and maximum residue depths
	core::Real min_depth_;
	core::Real max_depth_;

};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_ByMembraneDepthSelector_hh
