// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TrialMover
/// @brief performs a move and accepts it according to Monte Carlo accept/reject criterion.
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_moves_TrialMover_hh
#define INCLUDED_protocols_moves_TrialMover_hh

// Unit headers
#include <protocols/moves/TrialMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverStatistics.hh>

#include <core/scoring/ScoreType.hh>

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>
#include <string>

#include <utility/vector1.hh>


// Utility Headers


namespace protocols {
namespace moves {

///////////////////////////////////////////////////////////////////////////////
/// @brief A TrialMover applies a Mover and then accepts or rejects the move
/// according to a MonteCarlo object.
///
/// @detailed:
/// Each derived class should define its own apply() statement
/// the apply (mc) which requires a monte carlo object and only keeps
/// the move if the monte carlo test allows it
///
/// @authors Monica Berrondo
///
/// @last_modified August 16 2007
///////////////////////////////////////////////////////////////////////////////

// silly little enum to keep track of how we want to
// keep track of statistics from this Mover. It's
// important with the current design that these are
// kept in a hierarchically increasing order, so that
// the StatsType of i-1 keeps all of the same
// statistics as StatsType.
enum StatsType {
	all_stats = 1, // scores, accepts and rejects
	accept_reject, // only accept and reject counts
	no_stats	     // no stats. Always keep last.
};

	
/// @brief the MCResetMover applies a monte carlo reset
	
class MonteCarloUtil : public Mover {
public:
	MonteCarloUtil();
	MonteCarloUtil(protocols::moves::MonteCarloOP mc);
	virtual ~MonteCarloUtil();
	virtual void apply(Pose & pose);
	
	virtual void parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /* data */,
		protocols::filters::Filters_map const & /* filters */,
		protocols::moves::Movers_map const & /* movers */,
		core::pose::Pose const & /* pose */
	);
	
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;
	
	
private:
	protocols::moves::MonteCarloOP mc_;
	std::string mode_;
	
};
	
/// @brief A TrialMover applies a Mover and then accepts or rejects the move
/// according to a MonteCarlo object.
///
/// Common Methods:
///     TrialMover.apply
///     TrialMover.set_mc
class TrialMover : public Mover {
public:

	typedef core::Real Real;

public:

	TrialMover();

	// constructor with arguments -- BAD, MOVE TO .CC
	/// @brief Constructs a TrialMover
	/// trialmover = TrialMover( mover_in , mc_in )
	///
	/// Mover        mover_in   /object defining what move to make
	/// MonteCarlo   mc_in      /object defining acceptance
	TrialMover( MoverOP mover_in, MonteCarloOP mc_in );

	/// @brief Copy constructor
	TrialMover(TrialMover const & object_to_copy);

	~TrialMover();


	// set the weights for the score_type for ramping -- BAD, MOVE TO .CC
	virtual void initialize_weights(
		Real const start_weight,
		Real const end_weight,
		core::scoring::ScoreType score_type,
		int const ramp_cycles
	);

	/// @brief Performs a single trial, apply the Mover and accept or
	/// reject the move based on the MonteCarlo acceptance criteria
	///
	/// example(s):
	///     trialmover.apply(pose)
	/// See Also:
	///     Pose
	///     MonteCarlo
	///     MonteCarlo.boltzmann
	///     RepeatMover
	///     SequenceMover
	///     TrialMover
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	Real acceptance_rate() const;

	/// @brief Returns the number of moves accepted by this TrialMover
	int num_accepts() const;

	/// @brief Sets the MonteCarlo object to  <mc_in>
	///
	/// example(s):
	///     trialmover.set_mc(mc)
	/// See Also:
	///     MonteCarlo
	///     MonteCarlo.boltzmann
	///     TrialMover
	void set_mc( MonteCarloOP mc_in );

	/// @brief Returns the underlying Mover
	MoverOP mover() const {
		return mover_;
	}

	/// @brief Returns the underlying MonteCarlo object
	MonteCarlo const& mc () const {
		return *mc_;
	}

	StatsType keep_stats_type() const {
		return stats_type_;
	}

	void keep_stats_type( StatsType setting ) {
		stats_type_ = setting;
	}

	// @brief Sets the input Pose also for the contained Mover
	virtual void set_input_pose( PoseCOP pose );

	// @brief Sets the native Pose also for the contained Mover
	virtual void set_native_pose( PoseCOP pose );

	/// @brief Requires that a MonteCarlo object has already been
	/// loaded into the basic::datacache::DataMap in a prior MONTECARLOS objects
	/// section.
	virtual void parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const &
	);
	friend std::ostream &operator<< (std::ostream &os, TrialMover const &mover);

protected:

	MoverOP mover_;
	MonteCarloOP mc_;
	MoverStatistics stats_;

private:
	Real start_weight_;
	Real original_weight;
	Real ramp_weight;
	Real delta;
	StatsType stats_type_;
}; // TrialMover base class

} // moves
} // protocols


#endif //INCLUDED_protocols_moves_TrialMover_HH
