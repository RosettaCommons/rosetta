// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Ints.cc
///
/// @brief
/// @author
#include <devel/inv_kin_lig_loop_design/Ints.hh>
#include <numeric/random/random.hh>

#include <iostream>
#include <sstream>
#include <string>

#include <utility/vector1.hh>


#define FORVC(Iter,Type,Vec)  for( vector<Type>::const_iterator Iter  = (Vec).begin(); Iter != (Vec).end(); ++Iter)

namespace devel {

namespace inv_kin_lig_loop_design {

Ints::Ints() {
	b_is_init = false;
} // Int::Ints

Ints::Ints(const string& s) {
	b_is_init = false;
	fromString(s);
} // Int::Ints

void Ints::clear() {
	vInts.clear();
	b_is_init = false;
} // Ints::clear

void Ints::add(const int i) {
	vInts.push_back(i);
	b_is_init = false;
} // Ints::add

void Ints::add(const pair<int,int>& range) {
	if ( range.first <= range.second ) {
		vRanges.push_back(range);
	} else {
		vRanges.push_back( make_pair(range.second,range.first) );
	}
	b_is_init = false;
} // Ints::add

int Ints::getRandomInt() const {

	vInit.clear();

	if ( !b_is_init ) {
		vInit.insert( vInit.end(), vInts.begin(), vInts.end() );
		typedef pair<int,int> pii_t;
		FORVC ( i,pii_t,vRanges ) {
			for ( int k = i->first; k <= i->second; ++k ) {
				vInit.push_back(k);
			} // k
		} // i
		b_is_init = true;
	}

	if ( vInit.size() == 0 ) {
		return 0;
	} else {
		return vInit[ numeric::random::random_range(0,vInit.size()-1) ];
	}

} // Int::getRandomInt

void Ints::fromString(const string& s) {
	clear();

	enum State { STATE_0 , STATE_1 };

	State state = STATE_0;

	pair<string,string> p;

	for ( size_t k = 0; k < s.size(); ++k ) {

		switch( state ) {
		case STATE_0 :

			switch( s[k] ) {
			case ',' :
				add( atoi( p.first.c_str() ) );
				p.first.clear();
				break;
			case '-' :
				state = STATE_1;
				break;
			default :
				p.first += s[k];
				break;
			} // switch
			break;

		case STATE_1 :

			switch( s[k] ) {
			case ',' :
				add( make_pair( atoi(p.first.c_str()), atoi(p.second.c_str()) ) );
				p.first.clear();
				p.second.clear();
				state = STATE_0;
				break;
			default :
				p.second += s[k];
				break;
			} // switch

			break;
		default :
			assert( false );
		} // state

	} // k

	switch( state ) {
	case STATE_0 :
		add( atoi(p.first.c_str()));
		p.first.clear();
		break;
	case STATE_1 :
		add( make_pair( atoi(p.first.c_str()), atoi(p.second.c_str()) ) );
		p.first.clear();
		p.second.clear();
		break;
	default :
		assert( false );
	} // state

} // Ints::fromString

//#define FORVC(Iter,Type,Vector) for( vector<Type>::const_iterator Iter = (Vector).begin(); Iter != (Vector).end(); ++Iter )

const string Ints::toString() const {

	ostringstream out;
	FORVC ( k,int,vInts ) {
		if ( k != vInts.begin() ) {
			out << ",";
		}
		out << *k;
	} // k

	if ( vInts.size() != 0 && vRanges.size() != 0 ) {
		out << ",";
	}

	typedef pair<int,int> pii_t;
	FORVC ( k,pii_t,vRanges ) {
		if ( k != vRanges.begin() ) {
			out << ",";
		}
		out << k->first << "-" << k->second;
	} // k

	return out.str();

} // Ints::toString

ostream& operator<<(ostream& out, const Ints& ints) {
	out << ints.toString();
	return out;
} // operator<<

istream& operator>>(istream& in, Ints& ints) {
	string s;
	in >> s;
	ints.fromString(s);
	return in;
} // operator<<

} // namespace LoopDesign

} // namespace devel
