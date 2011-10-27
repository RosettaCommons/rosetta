// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @detailed
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author olange: ported from original bblum-rosetta++ version $

#ifndef INCLUDED_core_scoring_dssp_PairingsList_hh
#define INCLUDED_core_scoring_dssp_PairingsList_hh

#include <core/scoring/dssp/PairingsList.fwd.hh>
#include <core/types.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>

#include <utility/vector1.fwd.hh>


//// C++ headers
#include <cmath>
// AUTO-REMOVED #include <cstdlib> //required by GCC 4.3.2
#include <string>

namespace core {
namespace scoring {
namespace dssp {

class Pairing {
public:
  core::Size pos1;
  core::Size pos2;
  core::Size orientation;
  core::Size pleating;

	Pairing() :
		pos1( 0 ),
		pos2( 0 ),
		orientation( 0 ),
		pleating( 0 )
	{}

	Pairing( core::Size pos1_in, core::Size pos2_in ) :
		pos1( pos1_in),
		pos2( pos2_in),
		orientation( 0 ),
		pleating( 0 )
	{}

	//c'stor to translate from old-style version of pairing
	//Pairing( ObjexxFCL::FArray1A_int );

	Pairing( core::Size pos1_in, core::Size pos2_in, core::Size ori_in, core::Size pleat_in ) :
		pos1( pos1_in),
		pos2( pos2_in),
		orientation( ori_in),
		pleating( pleat_in )
	{}

	///@brief constant values that define orientation
	static core::Size const ANTI = 1;
	static core::Size const PARALLEL = 2;
	static core::Size const OUTWARDS = 1;
	static core::Size const INWARDS = 2;

	///@brief reverses the Pairing
	//Pairing reverse();

	///@brief returns a new reversed pairing
	//Pairing generate_reversed() const;

	///
	bool is_parallel() const {
		runtime_assert( orientation );
		return orientation == PARALLEL;
	}

	///
	bool is_anti() const {
		runtime_assert( orientation );
		return orientation == ANTI;
	}

	bool is_inwards() const {
		runtime_assert( pleating );
		return pleating == INWARDS;
	}

	bool is_outwards() const {
		runtime_assert( pleating );
		return pleating == OUTWARDS;
	}

	bool operator ==( Pairing const& p ) const {
		return ( (p.pos1 == pos1)
			&& ( p.pos2 == pos2 )
			&& ( p.orientation == orientation )
			&& ( p.pleating == pleating )
		);
	};

	core::Size get_register()  {
		return is_anti() ? pos1 + pos2 : std::abs( (int) pos1 - (int) pos2 );
	}


	bool operator < ( Pairing const& p ) const {
		return p.pos1 != pos1  ?  pos1 < p.pos1  :
			 (  p.pos2 != pos2  ? pos2 < p.pos2 :
				 ( p.orientation != orientation ? orientation < p.orientation : pleating < p.pleating ) );
	};
};

///@brief list of pairings
//typedef utility::vector1<Pairing> PairingsList;

//extern std::ostream& operator<< ( std::ostream& out, Pairing const& );
extern std::ostream& operator<< ( std::ostream& out, PairingsList const& p);

///@brief add pairings in pairing_file to list "pairings"
//extern void read_pairing_list( std::string pairing_file, PairingsList& pairings);
//extern void read_pairing_list( std::istream &is, PairingsList& pairings);

//extern bool has_orientation_and_pleating( PairingsList const& );

} //dssp
} //scoring
} //core

#endif

