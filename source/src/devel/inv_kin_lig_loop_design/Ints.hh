// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Ints.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_INTS_HH
#define DEVEL_INVKINLIGLOOPDESIGN_INTS_HH

#include <string>
#include <vector>

namespace devel {
namespace inv_kin_lig_loop_design {

using namespace std;

class Ints {
private:
	vector<int> vInts;
	vector<pair<int,int> > vRanges;

	mutable bool b_is_init;
	mutable vector<int> vInit;

public:
	Ints();
	Ints(const string& s);

	void clear();

	void add(const int i);
	void add(const pair<int,int>& range);

	int getRandomInt() const;

	int operator()() const { return getRandomInt(); }

	const string toString() const;
	void fromString(const string& s);

}; // class Ints

istream& operator>>(istream& in, Ints& ints);
ostream& operator<<(ostream& out, const Ints& ints);

} // namespace LoopDesign
} // namespace devel

#endif // DEVEL_LOOPDESIGN_INTS_HH
