// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PSSM2BfactorMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_PSSM2BfactorMover_hh
#define INCLUDED_protocols_simple_moves_PSSM2BfactorMover_hh

#include <protocols/simple_moves/PSSM2BfactorMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class PSSM2BfactorMover : public moves::Mover {
public:
	PSSM2BfactorMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	core::Real offset() const{ return offset_; }
	void offset( core::Real const o ){ offset_ = o; }

	core::Real scaling_factor() const { return scaling_factor_; }
	void scaling_factor( core::Real const s ){ scaling_factor_ = s; }

private:
	core::Real offset_, scaling_factor_; // dflt +8, x3; scaling factor to convert PSSM scores to temperature factor
};


} // simple_moves
} // protocols

#endif
