// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AddChainMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_AddChainMover_hh
#define INCLUDED_protocols_simple_moves_AddChainMover_hh

#include <protocols/simple_moves/AddChainMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class AddChainMover : public moves::Mover {
public:
	AddChainMover();

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

	void fname( std::string const f ){ fname_ = f; }
	std::string fname() const{ return fname_; }

	void new_chain( bool const n ){ new_chain_ = n; }
	bool new_chain() const{ return new_chain_; }
    
    void swap_chain_number( core::Size const wc ){ swap_chain_number_ = wc; }
	core::Size swap_chain_number() const { return swap_chain_number_; }
    
	void scorefxn( core::scoring::ScoreFunctionOP s );
	core::scoring::ScoreFunctionOP scorefxn() const;

	bool random_access() const{ return random_access_; }
	void random_access( bool const b ){ random_access_ = b; }
    
    void add_new_chain( core::pose::Pose & pose ) const; // Adds new chain to pose
    void swap_chain( core::pose::Pose & pose ) const; // Adds new chain to pose
    
private:
	std::string fname_; //dflt ""; pdb names to load (can accept a comma-separated list)
	bool new_chain_; //dflt true; add as a new chain?
	bool random_access_; //dflt false; if true randomly choose one file name from a list and work with that throughout the run.
  core::Size swap_chain_number_; //dflt 2; swap chain with specified chain number
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12; used to score the new pose
};


} // simple_moves
} // protocols

#endif
