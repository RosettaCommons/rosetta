// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RampingMover.cc
/// @brief Mover class for ramping a score function over the course of
/// apply evaluations.
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_moves_RampingMover_hh
#define INCLUDED_protocols_moves_RampingMover_hh

// Unit Headers
#include <protocols/moves/RampingMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

namespace protocols {
namespace moves {

class RampingFunc : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real Real;
public:
	virtual
	~RampingFunc();

	/// @brief Func must be defined (and finite) over the range from 0 and 1.
	virtual
	Real
	func( Real ) const = 0;

};

class LinearFunc : public RampingFunc
{
public:
	virtual ~LinearFunc();

	virtual
	Real
	func( Real x ) const
	{
		return x;
	}
};

class FastLinearFunc : public RampingFunc
{
public:
	FastLinearFunc(
		Real xval_start_ramp,
		Real xval_end_ramp
	);

	virtual ~FastLinearFunc();
	virtual
	Real
	func( Real x ) const;

private:
	Real xval_start_ramp_;
	Real xval_end_ramp_;

};


/// @brief Ramps rapidly from the starting value to the final value.
/// Not 1 at x=1.  Doesn't really finish at (1,1).
/// func(x) = 1 - exp( -1 * x * inv_xval_at_0p5 * 0.6931 );
class GeometricFunc : public RampingFunc
{
public:
	GeometricFunc( Real xval_at_0p5 );
	GeometricFunc();

	virtual ~GeometricFunc();

	virtual
	Real
	func( Real x ) const;

private:
	Real inv_xval_at_0p5_;

};

/// @brief Ramps slowly from the starting value to the final value
/// Non-zero for x = 0.  Doesn't really start at (0,0).
/// func(x) = exp( -1 * ( 1 - x ) / ( 1 - xval_at_0p5 ) * 0.6931 );
class InvGeometricFunc : public RampingFunc
{
public:
	InvGeometricFunc( Real xval_at_0p5 );
	InvGeometricFunc();

	virtual ~InvGeometricFunc();

	virtual
	Real
	func( Real x ) const;

private:
	Real inv_one_minus_xval_at_0p5_;

};


class RampingMover : public Mover {
public:
	// default constructor -- what good is a default constructor if you don't
	// have mutators for the majority of the parameters?
	RampingMover();

	~RampingMover();

	/// Takes a scorefunciton OP because the score function it modifies is shared
	/// between many movers.
	RampingMover(
		MoverOP mover_in,
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::scoring::ScoreType score_type_in,
		int outer_cycles_in,
		int inner_cycles_in,
		MonteCarloOP  mc_in,
		bool geometric_in = false
	);

	RampingMover(
		MoverOP mover_in,
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::scoring::EnergyMap start_weights,
		core::scoring::EnergyMap end_weights,
		int outer_cycles_in,
		int inner_cycles_in,
		MonteCarloOP  mc_in
	);

	virtual MoverOP clone() const;

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tags,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose);

	void sfxn( core::scoring::ScoreFunctionOP );
	void start_weight( core::Real start_weight_in );
	void end_weight  ( core::Real end_weight_in   );

	void set_func_for_weight( core::scoring::ScoreType scoretype, RampingFuncOP func );

private:

	// private functions used with weight ramping
	void update_weights( int round );
	void set_weights( core::scoring::EnergyMap const & emap );

	RampingFuncOP instantiate_rampfunc(
		std::string const & func_name,
		utility::tag::TagCOP tag_ptr ) const;

private:

	MoverOP mover_;
	core::scoring::ScoreFunctionOP scorefxn_;

	bool ramp_one_weight_; /// Are we ramping one weight or several weights?
	core::scoring::ScoreType score_type_; // The single weight that will be ramped

	core::scoring::EnergyMap start_weights_;
	core::scoring::EnergyMap end_weights_;
	core::scoring::EnergyMap intermediate_weights_;

	// Once it's set, its immutable.
	utility::vector1< RampingFuncCOP > ramping_funcs_for_weights_;

	int outer_cycles_, inner_cycles_;

	MonteCarloOP mc_;
}; // end class RampingMover

} // moves
} // protocols


#endif //INCLUDED_protocols_moves_RampingMover_HH
