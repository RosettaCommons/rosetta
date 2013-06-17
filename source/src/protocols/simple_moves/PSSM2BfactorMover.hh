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
	PSSM2BfactorMover(core::Real const min_in , core::Real const max_in);

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


	core::Real min_value() const { return min_value_; }
	void min_value( core::Real const s ){ min_value_ = s; }
	
	core::Real max_value() const { return max_value_; }
	void max_value( core::Real const s ){ max_value_ = s; }

private:
	core::Real min_value_,max_value_; // dflt -1, 5
};


} // simple_moves
} // protocols

#endif
