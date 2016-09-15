// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ModifyVariantTypeMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_simple_moves_ModifyVariantTypeMover_HH
#define INCLUDED_protocols_simple_moves_ModifyVariantTypeMover_HH

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.fwd.hh>


// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief Adds variant types to selected residues
class ModifyVariantTypeMover : public protocols::moves::Mover
{
public:
	// default constructor (nmoves=1)
	ModifyVariantTypeMover();

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	/// @brief Set the ResidueSelector used by this mover.
	///
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in);

	/// @brief Get the ResidueSelector used by this mover.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

private:
	/// @brief List of types to ADD.
	///
	utility::vector1<std::string> add_target_types_;

	/// @brief List of types to REMOVE.
	///
	utility::vector1<std::string> remove_target_types_;

	/// @brief ResidueSelector specifying residues to which this should be applied.
	///
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;
};

} // moves
} // protocols


#endif
