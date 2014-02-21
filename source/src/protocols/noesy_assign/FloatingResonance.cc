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
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/FloatingResonance.hh>

// Package Headers
#include <protocols/noesy_assign/PeakCalibrator.hh>
#include <protocols/noesy_assign/FoldResonance.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

//// C++ headers
#include <string>
#include <utility/vector1.hh>



static basic::Tracer tr("protocols.noesy_assign.resonances");

using core::Real;
using namespace core;
using namespace basic;

namespace protocols {
namespace noesy_assign {

FloatingResonance::FloatingResonance() {}

FloatingResonance::FloatingResonance( Resonance const& resin, FloatList const& partner, ResonanceList* reslist ) :
	Resonance( resin ),
	partner_ids_( partner ),
	res_list_( reslist )
{}

FloatingResonance::~FloatingResonance() {}

core::Real FloatingResonance::pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const {
	Real min_match = 100000;
	for ( FloatList::const_iterator it = partner_ids_.begin(); it!=partner_ids_.end(); ++it ) {
		Real m = (*res_list_)[ *it ]._pmatch( peakfreq, error, folder );
		if ( m < min_match ) min_match = m;
	}
	return min_match;
}

void FloatingResonance::write_to_stream( std::ostream& os ) const {
	Parent::write_to_stream( os );
	_write_partner_ids( os );
}

void FloatingResonance::_write_partner_ids( std::ostream& os ) const {
	if ( partner_ids_.size() > 1 ) {
		os << " [";
		long pos( os.tellp() );
		for ( FloatList::const_iterator it = partner_ids_.begin(); it != partner_ids_.end(); ++it ) {
			os << *it;
			pos = os.tellp();
			os << ",";
		}
		os.seekp( pos );
		os << "]";
	}
}

void FloatingResonance::write_to_stream( std::ostream& os, core::chemical::AA aa ) const {
	Parent::write_to_stream( os, aa );
	_write_partner_ids( os );
}

} //NoesyAssign
} //devel
