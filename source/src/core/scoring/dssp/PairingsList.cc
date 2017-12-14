// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @details
///
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <core/scoring/dssp/PairingsList.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1A.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

//// C++ headers
#include <core/types.hh>
#include <string>

#include <utility/vector1.hh>


static basic::Tracer TR( "core.scoring.dssp" );

namespace core {
namespace scoring {
namespace dssp {

using core::Real;
using namespace basic;
using namespace ObjexxFCL;


Pairing::Pairing( ObjexxFCL::FArray1A_int data) {
	pos1_        = data(1);
	pos2_        = data(2);
	orientation_ = data(3);
	pleating_    = data(4);
}

Pairing
Pairing::reverse() {
	Size tmp = pos2_;
	pos2_ = pos1_;
	pos1_ = tmp;

	if ( orientation_ == PARALLEL ) { // flip pleating
		if ( pleating_ == 1 ) {
			pleating_ = 2;
		} else if ( pleating_ == 2 ) {
			pleating_ = 1;
		} else {
			std::cout << "unrecognized pleating:" << format::SS( pleating_ ) << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
	return *this;
}

Pairing
Pairing::generate_reversed() const {
	Pairing p(*this);
	p.reverse();
	return p;
}
void read_pairing_list( std::string pairing_file, PairingsList& pairings)
{
	utility::io::izstream pairing_stream( pairing_file );
	if ( !pairing_stream ) {
		TR.Fatal << "can't open pairings file!!!" << pairing_file << std::endl;
		pairing_stream.close();
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}
	read_pairing_list( pairing_stream, pairings );
	pairing_stream.close();
}

void read_pairing_list( std::istream& pairing_stream, PairingsList& pairings) {
	std::string line;
	Size a, b, c, d;
	while ( getline( pairing_stream, line ) ) {
		std::istringstream line_stream( line );
		// a=i, b=j, c=orientation(1 or 2), d=pleating(1 or 2)
		std::string o, pleat;
		line_stream >> a >> b >> o >> pleat;
		if ( line_stream.fail() || o.size() != 1 ) {
			TR.Error << "unable to parse " << line << std::endl;
			continue;
		}

		if ( o == "A" || o == "1" ) {
			c = 1;
		} else if ( o == "P" || o == "2" ) {
			c = 2;
		} else if ( o == "X" ) {
			c = 0;
		} else {
			std::cout << "bad orientation: " << o << std::endl;
			continue;
		}

		if ( pleat == "O" || pleat== "1" ) {
			d = 1;
		} else if ( pleat == "I" || pleat == "2" ) {
			d = 2;
		} else if ( pleat == "X" ) {
			d = 0;
		} else {
			d = 3;  // avoids [-Wuninitialized] warning ~Labonte
			TR.Warning << "bad pleating: " << pleat << std::endl;
		}

		if ( ( a < 1 || b < 1 ) || ( a == b ) || ( c != 1 && c != 2 && c != 0 ) ||
				( d != 1 && d != 2 && d != 0 ) ) {
			std::cout << "bad pairing:" <<
				format::SS( a ) << format::SS( b ) << format::SS( c ) << format::SS( d ) << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		Pairing p( a, b, c, d);
		if ( a > b ) p.reverse();
		pairings.push_back( p );
	}
} // read_pairings

std::ostream& operator<< ( std::ostream& out, Pairing const& p) {
	out << format::RJ(5, p.Pos1() ) << format::RJ(5, p.Pos2() ) << " "
		<< ( p.Orientation() ? ( p.is_parallel() ? "P" : "A") : "X" ) << " "
		<< ( p.Pleating() ? ( p.is_inwards() ? "I" : "O" ) : "X" );
	return out;
}

std::ostream& operator<< ( std::ostream& out, PairingsList const& p) {
	for ( auto const & it : p ) {
		out << it << "\n";
	}
	return out;
}

bool has_orientation_and_pleating( PairingsList const& pairings ) {
	for ( auto const & pairing : pairings ) {
		if ( pairing.Orientation() == 0 || pairing.Pleating() == 0 ) return false;
	}
	return true;
}

} // dssp
} // scoring
} // core
