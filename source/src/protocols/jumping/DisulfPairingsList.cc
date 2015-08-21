// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @details
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/jumping/DisulfPairingsList.hh>

// Package Headers

// Project Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1A.hh>


// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

// #include <basic/options/option.hh>
// #include <basic/options/keys/OptionKeys.hh>

// numeric headers
// #include <numeric/random/random.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.jumping" );

using core::Real;
using namespace core;
using namespace basic;
using namespace ObjexxFCL;
//using namespace basic::options;

namespace protocols {
namespace jumping {

DisulfPairing::DisulfPairing( ObjexxFCL::FArray1A_int data) {
	pos1        = data(1);
	pos2        = data(2);
	seq_sep     = data(3);
	ss_type     = data(4);
}

void read_pairing_list( std::string disulf_pairing_file, DisulfPairingsList& disulf_pairings)
{
	utility::io::izstream disulf_pairing_stream( disulf_pairing_file );
	if ( !disulf_pairing_stream ) {
		tr.Fatal << "can't open disulf_pairings file!!!" << disulf_pairing_file << std::endl;
		disulf_pairing_stream.close();
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}
	read_disulf_pairing_list( disulf_pairing_stream, disulf_pairings );
	disulf_pairing_stream.close();
}

void read_disulf_pairing_list( std::istream& disulf_pairing_stream, DisulfPairingsList& disulf_pairings) {
	std::string line;
	Size a,b,c,d;
	while ( getline( disulf_pairing_stream, line ) ) {
		std::istringstream line_stream( line );
		// a=i, b=j, c=orientation(1 or 2), d=pleating(1 or 2)
		//std::string o, pleat;
		line_stream >> a >> b >> c >> d;
		if ( line_stream.fail() ) {
			std::cout << "[ERROR] unable to parse " << line << std::endl;
			continue;
		}

		DisulfPairing p( a, b, c, d);
		//if ( a > b ) p.reverse();
		disulf_pairings.push_back( p );

	}
} // read_disulf_pairings

std::ostream& operator<< ( std::ostream& out, DisulfPairing const& p) {
	out << format::RJ(5, p.pos1 ) << format::RJ(5, p.pos2 ) << " "
		<< format::RJ(5, p.seq_sep) << format::RJ(5, p.ss_type );
	return out;
}

std::ostream& operator<< ( std::ostream& out, DisulfPairingsList const& p) {
	for ( DisulfPairingsList::const_iterator it= p.begin(),
			eit = p.end(); it!=eit; ++it ) {
		out << (*it) << "\n";
	}
	return out;
}

} // jumping
} // protocols
