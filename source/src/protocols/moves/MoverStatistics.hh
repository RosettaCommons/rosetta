// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Keeps track of the acceptance_rate of a Mover
/// @author Monica Berrondo August 16 2007


#ifndef INCLUDED_protocols_moves_MoverStatistics_hh
#define INCLUDED_protocols_moves_MoverStatistics_hh


// Package headers
#include <core/types.hh>

//#include <basic/basic.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <vector>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace moves {


///////////////////////////////////////////////////////////////////////////////
// @brief MoverStatistics keeps track of the acceptance_rate for a mover


class MoverStatistics : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MoverStatistics();

	// default constructor
	MoverStatistics():
		accepted_(0),
		rejected_(0),
		thermal_accepts_(0),
		downhill_accepts_(0)

	{}

	void accepted( bool result ) const
	{
		if ( result ) accepted_++;
		else rejected_++;
	}

	void add_score( core::Real score_in ) { score_.push_back(score_in); }

	void clear_score() { score_.clear(); }

	void print( MonteCarloOP mc, std::string const & type );

	core::Real acceptance_rate() const { return (core::Real)accepted_/(accepted_+rejected_+1e-100); }

	int num_accepted() const { return accepted_; }

	void clear()
	{
		accepted_ = 0;
		rejected_ = 0;
		thermal_accepts_ = 0;
		downhill_accepts_ = 0;

	}

private:
	mutable int accepted_; ///< number of accepted moves
	mutable int rejected_; ///< number of rejected moves
	std::vector < core::Real > score_;
	Size thermal_accepts_;
	Size downhill_accepts_;
}; // MoverStatistics class
} // moves
} // rosetta


#endif
