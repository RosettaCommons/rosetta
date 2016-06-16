// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/cyclic_peptide/FlipChiralityMover.hh
/// @brief Mirrors a pose.
/// @author Parisa Hosseiznadeh (parisah@uw.edu) and Vikram Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_FlipChiralityMover_hh
#define INCLUDED_protocols_cyclic_peptide_FlipChiralityMover_hh

// Unit headers

// Package headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#include <protocols/moves/Mover.hh>
#include <protocols/cyclic_peptide/FlipChiralityMover.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

class FlipChiralityMover : public moves::Mover {

public:
	FlipChiralityMover();
	virtual ~FlipChiralityMover();
	FlipChiralityMover( FlipChiralityMover const &src );

	virtual void apply( core::pose::Pose & );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

	numeric::xyzVector <core::Real> const & get_normal();
	numeric::xyzVector <core::Real> const & get_center(core::select::residue_selector::ResidueSubset, core::pose::Pose const &);
	virtual void set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );
private:
	//Private functions:
	numeric::xyzVector<core::Real> normal_;
	numeric::xyzVector<core::Real> center_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
	bool normal_assigned_;
	bool center_assigned_;

	void set_center (core::Real, core::Real, core::Real);
	void set_normal (core::Real,core::Real, core::Real);
	numeric::xyzVector <core::Real> const calculate_reflect (core::pose::Pose const &, core::Size, core::Size);
	numeric::xyzVector <core::Real> center_mass (core::select::residue_selector::ResidueSubset, core::pose::Pose const &);

};

} // moves
} // protocols

#endif
