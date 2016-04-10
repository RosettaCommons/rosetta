// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/jumping/PairingTemplate
/// @brief header file for ClassicAbinitio protocol
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_jumping_DisulfPairingsList_hh
#define INCLUDED_protocols_jumping_DisulfPairingsList_hh

// Unit Headers
#include <protocols/jumping/DisulfPairingsList.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility headers
// we need utility/vector1.hh because of declarations in fwd.hh
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.fwd.hh>

//// C++ headers
#include <string>

#ifdef PYROSETTA
#include <ObjexxFCL/FArray1A.hh>
#endif

namespace protocols {
namespace jumping {

class DisulfPairing {
public:
	core::Size pos1;
	core::Size pos2;
	core::Size seq_sep;
	core::Size ss_type;

	DisulfPairing() :
		pos1( 0 ),
		pos2( 0 ),
		seq_sep( 0 ),
		ss_type( 0 )
	{}

	DisulfPairing( core::Size pos1_in, core::Size pos2_in ) :
		pos1( pos1_in),
		pos2( pos2_in),
		seq_sep( 0 ),
		ss_type( 0 )
	{}

	//c'stor to translate from old-style version of pairing
	DisulfPairing( ObjexxFCL::FArray1A_int );

	DisulfPairing( core::Size pos1_in, core::Size pos2_in, core::Size ori_in, core::Size pleat_in ) :
		pos1( pos1_in),
		pos2( pos2_in),
		seq_sep( ori_in),
		ss_type( pleat_in )
	{}

	bool operator ==( DisulfPairing const& p ) const {
		return ( (p.pos1 == pos1)
			&& ( p.pos2 == pos2 )
			&& ( p.seq_sep == seq_sep )
			&& ( p.ss_type == ss_type )
		);
	};

	// bool operator < ( DisulfPairing const& p ) const {
	//  return p.pos1 != pos1  ?  pos1 < p.pos1  :
	//    (  p.pos2 != pos2  ? pos2 < p.pos2 :
	//     ( p.seq_sep != seq_sep ? seq_sep < p.seq_sep : ss_type < p.ss_type ) );
	// };
};

/// @brief list of pairings
//typedef utility::vector1<DisulfPairing> DisulfDisulfPairingsList;

extern std::ostream& operator<< ( std::ostream& out, DisulfPairing const& );
extern std::ostream& operator<< ( std::ostream& out, DisulfPairingsList const& p);

/// @brief add pairings in pairing_file to list "pairings"
//extern void read_disulf_pairing_list( std::string disulf_pairing_file, DisulfPairingsList& disulf_pairings);
extern void read_disulf_pairing_list( std::istream &is, DisulfPairingsList& disulf_pairings);


} //protocols
} //jumping

#endif
