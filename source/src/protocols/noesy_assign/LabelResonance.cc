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
/// @details
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/LabelResonance.hh>
#include <protocols/noesy_assign/ProtonResonance.hh>

#include <core/id/NamedAtomID.hh>

// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
//// C++ headers
#include <string>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.resonances" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

LabelResonance::LabelResonance() {}

LabelResonance::LabelResonance(
	core::Size label,
	core::Real freq,
	core::Real error,
	core::id::NamedAtomID const& id,
	core::chemical::AA aa,
	core::Real intensity
) : Resonance( label, freq, error, id, aa, intensity )
{}

LabelResonance::~LabelResonance() {}

/// @brief match the proton and corresponding label atom at same time
bool LabelResonance::match2D(
	core::Real /*proton_freq*/, //proton_frequenices are ignored
	core::Real proton_error,
	FoldResonance const& /*proton_folder*/,
	core::Real label_freq,
	core::Real label_error,
	FoldResonance const& label_folder,
	ResonancePairs& matches
) const {
	// only pseudo4D peaks are matched via the heavy-atom
	if ( proton_error < 99 ) return false;

	// if the label doesn't match ...
	if ( !match( label_freq, label_error, label_folder ) ) return false;

	// if we have no protons ... we cannot match
	if ( !has_connected_resonances() ) return false;

	// does the label match?
	for ( ResonanceAPs::const_iterator pit = connected_resonances().begin(); pit != connected_resonances().end(); ++pit ) {
		ResonanceOP r( *pit );
		Resonance const& proton( *r );
		matches.push_back( std::make_pair( proton.label(), label() ) );
	}
	return true;
}

} //NoesyAssign
} //devel
