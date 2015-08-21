// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ScoreMover.hh
/// @brief
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_simple_moves_ScoreMover_hh
#define INCLUDED_protocols_simple_moves_ScoreMover_hh

// Unit headers
#include <protocols/simple_moves/ScoreMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class ScoreMover : public moves::Mover {

public:
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	/// @brief
	///  empty constructor fills values with the values
	///  read in from the commandline
	ScoreMover();

	~ScoreMover();

	/// @brief constructor
	/// @details creates the ScoreMover with the names passed in rather than
	///  taken from the commandline
	///  patch is not necessary
	ScoreMover( std::string const &, std::string const & patch = "" );

	/// @brief constructor
	/// @details creates the ScoreMover with the scorefunction itself passed in
	ScoreMover( ScoreFunctionOP );

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & );
	void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP cache );

	static void register_options();

	/// @brief add an rms to the score_map
	/// TODO possibly find a better way to do this?
	/// for now, there are too many different rmsd calculation
	/// functions to be able to do the actual calculation in the protocols::moves::Mover
	void insert_rms( core::Real rms ) { score_map_["rms"] = rms; }

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	virtual void test_move( Pose & pose )
	{
		apply(pose);
	}

	void set_verbose( bool value ) { verbose_ = value; }

	ScoreFunctionOP score_function() const;

	void set_score_file( std::string scorefile ) { scorefile_ = scorefile; }

private:
	ScoreFunctionOP score_function_;
	std::map< std::string, core::Real > score_map_;
	bool verbose_;
	std::string scorefile_;

};

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_ScoreMover_HH
