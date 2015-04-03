// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/moves/IteratedConvergenceMover.hh
/// @brief  Mover class to repeatedly apply a submover until filter convergence is reached
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_moves_IteratedConvergenceMover_hh
#define INCLUDED_protocols_moves_IteratedConvergenceMover_hh

// Unit headers
#include <protocols/moves/IteratedConvergenceMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers

#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @brief A mover that repeatedly applies a sub-mover (up to a given maximum)
/// until the given filter returns values within a certain delta for a given
/// number of cycles.
class IteratedConvergenceMover : public Mover {
public:

public:
	/// @brief default constructor
	IteratedConvergenceMover();

	IteratedConvergenceMover( MoverOP submover, filters::FilterCOP filter, core::Real delta=0.1, core::Size cycles=1, core::Size maxcycles=1000 );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~IteratedConvergenceMover();

	// copy constructor
	IteratedConvergenceMover( IteratedConvergenceMover const & other );

	// methods
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );

	/// @brief required in the context of the parser/scripting scheme
	virtual MoverOP fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual MoverOP clone() const;

	// setters
	void submover( MoverOP mover );
	void filter( filters::FilterCOP filter );
	void delta( core::Real delta );
	void cycles( core::Size cycles );
	void maxcycles( core::Size maxcycles );

	// accessors
	MoverOP submover() const;
	filters::FilterCOP filter() const;
	core::Real delta() const;
	core::Size cycles() const;
	core::Size maxcycles() const;

private:
	MoverOP submover_;
	filters::FilterCOP filter_;
	core::Real delta_;
	core::Size cycles_;
	core::Size maxcycles_;

};


} // moves
} // protocols

#endif
