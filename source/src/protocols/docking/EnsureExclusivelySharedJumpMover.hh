// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/EnsureExclusivelySharedJumpMover.hh
/// @brief Set up foldtree such that there exists a jump that builds all selected residues and does not build any unselected residues
/// @author Jack Maguire

#ifndef INCLUDED_protocols_docking_EnsureExclusivelySharedJumpMover_hh
#define INCLUDED_protocols_docking_EnsureExclusivelySharedJumpMover_hh

// Unit Headers
#include <protocols/docking/EnsureExclusivelySharedJumpMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace docking {

class EnsureExclusivelySharedJumpMover : public moves::Mover {
public:
	/// @brief Default constructor
	EnsureExclusivelySharedJumpMover() = default;

	EnsureExclusivelySharedJumpMover( EnsureExclusivelySharedJumpMover const & ) = default;
	EnsureExclusivelySharedJumpMover & operator=( EnsureExclusivelySharedJumpMover const & rhs ) = default;

	EnsureExclusivelySharedJumpMover( core::select::residue_selector::ResidueSelectorCOP selector );

	~EnsureExclusivelySharedJumpMover() override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void apply( core::pose::Pose & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

public:
	///@brief ExclusivelySharedJumpSelector selects the jump that builds ALL and ONLY the residues selected by this selector (getter)
	core::select::residue_selector::ResidueSelectorCOP
	selector() const {
		return selector_;
	}

	///@brief ExclusivelySharedJumpSelector selects the jump that builds ALL and ONLY the residues selected by this selector (setter)
	void
	set_selector( core::select::residue_selector::ResidueSelectorCOP setting ){
		selector_ = setting;
	}

private:
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
};

}
}
#endif
