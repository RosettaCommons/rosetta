// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   SetupNCSMover.hh
/// @brief  Sets up NCS restraints
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_moves_symmetry_SetupNCSMover_hh
#define INCLUDED_protocols_moves_symmetry_SetupNCSMover_hh

#include <protocols/moves/symmetry/SetupNCSMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace moves {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////

class SetupNCSMover : public Mover {
public:
	// default constructor
	SetupNCSMover();

	// initialize using a single NCS pair
	SetupNCSMover( std::string src, std::string tgt );

	// one->many
	SetupNCSMover( std::string src, utility::vector1<std::string> tgt );

	// many->many
	SetupNCSMover( utility::vector1<std::string> src, utility::vector1<std::string> tgt );

	~SetupNCSMover();

	// add an ncs group
	void add_group( std::string src, std::string tgt );


	void set_defaults();

	// clear all groups
	void clear();

	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new SetupNCSMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose );

	virtual void parse_my_tag( 
			utility::tag::TagPtr const tag,
			moves::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

	virtual std::string get_name() const;

	// getters/setters
	void set_bb( bool bb_in ) { bb_ = bb_in; }
	bool bb( ) { return bb_; }
	void set_chi( bool chi_in ) { chi_ = chi_in; }
	bool chi( ) { return chi_; }

	void set_limit( core::Real limit_in) { limit_ = limit_in; }
	core::Real limit( ) { return limit_; }
	void set_weight( core::Real wt_in) { wt_ = wt_in; }
	core::Real weight( ) { return wt_; }

private:
	utility::vector1< std::string > src_;
	utility::vector1< std::string > tgt_;

	bool bb_, chi_;
	core::Real limit_, wt_;
};

} // symmetry
} // moves
} // rosetta
#endif
