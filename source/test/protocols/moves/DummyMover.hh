// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/DummyMover.cxxtest.hh
/// @brief  dummy mover object intendet to help test other Movers classes.
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_DummyMover_HH
#define INCLUDED_protocols_moves_DummyMover_HH


#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>

#include <vector>


/// @details This is special mover class for testing other movers.
///  contain information regardin number of call/sequence of call's etc.
///
class DummyMover : public protocols::moves::Mover
{
public:
	DummyMover(int id=0) : id_(id), call_count_(0) {};

	virtual void apply( core::pose::Pose& ) { call_records_.push_back(id_); call_count_++; };
	virtual std::string get_name() const { return "DummyMover"; }
	int call_count() { return call_count_; };

	static void reset() ///< Reset call record to []
	{ call_records_.resize(0); };

	static std::vector<int> & call_records() { return call_records_; };

private:
	int id_; ///< id of this mover, not nessery uniq
	int call_count_; ///< number of times 'apply' get called.

	static std::vector<int> call_records_;
};

std::vector<int> DummyMover::call_records_;

#endif // INCLUDED_test_protocols_moves_DummyMover_HH
