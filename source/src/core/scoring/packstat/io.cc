// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/io.cc
///
/// @brief
/// @author will sheffler

#include <core/scoring/packstat/io.hh>

#include <core/scoring/packstat/SimplePDB.hh>
#include <core/scoring/packstat/SimplePDB_Atom.hh>


#include <cstring>
#include <iostream>
#include <sstream>

#include <numeric/xyzVector.io.hh>


namespace core {
namespace scoring {
namespace packstat {


std::istream &
operator>> (
	std::istream & in,
	SimplePDB & pdb
) {
	SimplePDB_Atom atom;
	while ( in >> atom ) {
		// std::cerr << "reading line" << std::endl;
		if ( atom.num == -123456 ) {
			// std::cerr << "/*TER or*/ ENDMDL" << std::endl;
			break; // hit a /*TER or*/ ENDMDL, so stop
		}
		if ( atom.num != -12345 ) {
			pdb.atoms_.push_back(atom);
		}
	}
	pdb.remove_surface_waters();
	return in;
}

std::ostream &
operator<< (
	std::ostream & out,
	SimplePDB const & pdb
) {
	out << "SimplePDB:" << std::endl;
	for ( SPAtomCIter i = pdb.atoms_.begin(); i != pdb.atoms_.end(); ++i ) {
		out << "  " << *i << std::endl;
	}
	return out;
}

template<class T>
void read_stoopid(char * buf, size_t start, size_t end, T & field ) {
	if ( start == end ) {
		field = buf[start];
		return;
	}
	char smbuf[99];
	for ( int i = 0; i < 99; ++i ) smbuf[i] = 0;
	strncpy(smbuf,buf+start,end-start+1);
	std::istringstream iss(smbuf);
	iss >> field;
}

/*
reads an atom record like this, assuming chain is col 21 and all others present
11 13
ATOM     12  N   GLU A   2     -13.565  31.875  -5.182  1.00 51.33           N
ATOM    391  CG1AILE A  63       0.107  41.079  16.997  0.50 20.29           C

*/
std::istream &
operator>> (
	std::istream & in_raw,
	SimplePDB_Atom & atom
) {
	using namespace std;
	char buf[999];
	in_raw.getline(buf,999);
	atom.whole_line = buf;
	// std::cerr << "read line: " << buf << std::endl;
	//  istringstream in(buf);
	//  in >> atom.ATOM;
	atom.ATOM = atom.whole_line.substr( 0, 6 );
	if ( atom.ATOM != "ATOM  " && atom.ATOM != "HETATM" ) { // not an atom
		// std::cerr << "failed to read ATOM or HETATM! atom.ATOM is '" << atom.ATOM << "'" << std::endl;
		atom.num = -12345;
		if ( /*atom.ATOM == "TER" ||*/ atom.ATOM == "ENDMDL" ) {
			atom.num = -123456;
		}
		return in_raw;
	}
	read_stoopid<int   >(buf, 6,10,atom.num );
	read_stoopid<string>(buf,11,16,atom.type);
	read_stoopid<string>(buf,17,20,atom.res );
	read_stoopid<char  >(buf,21,21,atom.chain );
	read_stoopid<int   >(buf,22,26,atom.resnum );
	read_stoopid<char  >(buf,26,26,atom.icode);
	read_stoopid<PackstatReal >(buf,29,37,atom.x );
	read_stoopid<PackstatReal >(buf,38,45,atom.y );
	read_stoopid<PackstatReal >(buf,46,54,atom.z );
	read_stoopid<PackstatReal >(buf,55,59,atom.occ );
	read_stoopid<PackstatReal >(buf,60,65,atom.bfac );
	read_stoopid<string>(buf,66,99,atom.lastcol );
	return in_raw;
}


std::ostream &
operator<< (
	std::ostream & out,
	SimplePDB_Atom const & atom
) {
	out << atom.ATOM << " "
		<< atom.num << " "
		<< atom.type << " "
		<< atom.res << " "
		<< atom.resnum << " "
		<< atom.chain << " "
		<< atom.num << " "
		<< atom.x << " "
		<< atom.y << " "
		<< atom.z << " "
		<< atom.occ << " "
		<< atom.bfac << " "
		<< atom.lastcol;
	return out;
}

std::ostream &
operator<< (
	std::ostream & out,
	Sphere const & sphere
) {
	out << "Sphere( " << sphere.xyz << ", " << sphere.radius << " )";
	return out;
}


} // namespace packstat
} // namespace scoring
} // namespace core
