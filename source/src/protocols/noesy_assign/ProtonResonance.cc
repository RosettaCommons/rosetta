// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ProtonResonance.hh>

#include <core/id/NamedAtomID.hh>

// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/exit.hh>


static basic::Tracer tr( "protocols.noesy_assign.resonances" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

ProtonResonance::ProtonResonance() {}

ProtonResonance::ProtonResonance(
	core::Size label,
	core::Real freq,
	core::Real error,
	core::id::NamedAtomID const& id,
	core::chemical::AA aa,
	core::Real intensity )
: Resonance( label, freq, error, id, aa, intensity )
{}

ProtonResonance::~ProtonResonance() = default;


/// @brief match the proton and corresponding label atom at same time
bool ProtonResonance::match2D(
	core::Real proton_freq,
	core::Real proton_error,
	FoldResonance const& proton_folder,
	core::Real label_freq,
	core::Real label_error,
	FoldResonance const& label_folder,
	ResonancePairs& matches
) const {
	// pseudo4D peaks with two heavy-atoms and one proton are matched via the heavy-atom
	if ( proton_error >= 99 ) return false;

	// if the proton doesn't match ...
	if ( !match( proton_freq, proton_error, proton_folder ) ) return false;

	// if we have no label -- we are not matching
	if ( !has_connected_resonances() ) return false;

	// does label match ?
	if ( first_connected_resonance().match( label_freq, label_error, label_folder ) ) {
		matches.push_back( std::make_pair( label(), first_connected_resonance().label() ) );
		return true;
	}
	return false;
}

} //NoesyAssign
} //devel
