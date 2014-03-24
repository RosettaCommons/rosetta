// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Assigns a ResidueTypeConstraint to a pose. 
///@author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_ResidueTypeConstraintMover_hh
#define INCLUDED_protocols_simple_moves_ResidueTypeConstraintMover_hh

#include <protocols/simple_moves/ResidueTypeConstraintMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ResidueTypeConstraint.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace simple_moves {

class ResidueTypeConstraintMover : public protocols::moves::Mover {

public:
	typedef core::scoring::constraints::ResidueTypeConstraint ResidueTypeConstraint;
	typedef core::scoring::constraints::ResidueTypeConstraintOP ResidueTypeConstraintOP;
	typedef core::scoring::constraints::ResidueTypeConstraintCOP ResidueTypeConstraintCOP;

public:
	ResidueTypeConstraintMover();
	virtual ~ResidueTypeConstraintMover();
	ResidueTypeConstraintMover( std::string const & );

	void constraint_file( std::string const & );

	void constraint_set( ResidueTypeConstraintCOP );
	ResidueTypeConstraintOP constraint_set();
	ResidueTypeConstraintCOP constraint_set() const;

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, Pose const & );

private:
	std::string AA_name3_;
	core::Real	favor_bonus_;
};

} // moves
} // protocols

#endif
