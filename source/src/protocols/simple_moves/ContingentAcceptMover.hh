// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ContingentAcceptMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_ContingentAcceptMover_hh
#define INCLUDED_protocols_simple_moves_ContingentAcceptMover_hh

#include <protocols/simple_moves/ContingentAcceptMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class ContingentAcceptMover : public moves::Mover {
public:
	ContingentAcceptMover();
	// Undefined, commenting out to fix PyRosetta build  ContingentAcceptMover(core::Real const min_in , core::Real const max_in);

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

		protocols::filters::FilterOP filter() const;
		void filter( protocols::filters::FilterOP f );

    protocols::moves::MoverOP mover() const;
    void mover( protocols::moves::MoverOP m );

   // protocols::filters::mover_FilterOP filter() const { return filter_; };
	//	void filter( protocols::filters::FilterOP f ){ filter_ = f; }

  //  protocols::moves::MoverOP mover() const { return mover_; };
	//	void mover( protocols::moves::MoverOP m ){ mover_ = m; }

	core::Real delta() const { return delta_; }
	void delta( core::Real const s ){ delta_ = s; }

private:
    protocols::filters::FilterOP filter_;
	  moves::MoverOP mover_;
    core::Real delta_;
};


} // simple_moves
} // protocols

#endif
