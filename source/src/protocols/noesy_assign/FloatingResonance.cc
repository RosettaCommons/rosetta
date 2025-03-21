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
#include <protocols/noesy_assign/FloatingResonance.hh>

// Package Headers
#include <protocols/noesy_assign/FoldResonance.fwd.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>

// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers
#include <string>

#include <protocols/noesy_assign/ResonanceList.hh> // AUTO IWYU For ResonanceList


static basic::Tracer tr( "protocols.noesy_assign.resonances" );

using core::Real;
using namespace core;
using namespace basic;

namespace protocols {
namespace noesy_assign {

FloatingResonance::FloatingResonance() = default;

FloatingResonance::FloatingResonance( Resonance const& resin, FloatList const& partner, ResonanceList* reslist ) :
	Resonance( resin ),
	partner_ids_( partner ),
	res_list_( reslist )
{
	is_representative_resonance_ = true;
	for ( core::Size it : partner ) {
		is_representative_resonance_ = is_representative_resonance_ && !( it < resin.label() );
	}
}

FloatingResonance::~FloatingResonance() = default;


core::Size FloatingResonance::float_label( core::Size ifloat ) const {
	for ( auto it = partner_ids_.begin(); it!=partner_ids_.end(); ++it, --ifloat ) {
		if ( ifloat==1 ) return *it;
	}
	tr.Fatal << "Cannot find float_label for value " << ifloat << std::endl;
	utility_exit_with_message("Cannot obtain float label.");
	return 99999999;
}

void FloatingResonance::write_to_stream( std::ostream& os ) const {
	Parent::write_to_stream( os );
	_write_partner_ids( os );
}

void FloatingResonance::_write_partner_ids( std::ostream& os ) const {
	if ( partner_ids_.size() > 1 ) {
		os << " [";
		long pos( os.tellp() );
		for ( core::Size partner_id : partner_ids_ ) {
			os << partner_id;
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

/// @brief match the proton and corresponding label atom at same time
bool FloatingResonance::match2D(
	core::Real proton_freq,
	core::Real proton_error,
	FoldResonance const& proton_folder,
	core::Real label_freq,
	core::Real label_error,
	FoldResonance const& label_folder,
	ResonancePairs &matches
) const {
	// pseudo4D peaks with two heavy-atoms and one proton are matched via the heavy-atom
	if ( !is_representative_resonance_ ) return false;
	if ( proton_error >= 99 && is_proton() ) return false;
	if ( proton_error < 99 && !is_proton() ) return false;

	for ( core::Size partner_id : partner_ids_ ) {
		if ( is_proton() ) { // proton-matching
			Resonance const& proton( (*res_list_)[ partner_id ] );
			//   tr.Trace << "try proton float-match with reso: " <<  proton.atom() << std::endl;
			bool proton_match = proton._pmatch( proton_freq, proton_error, proton_folder ) <= 1.0;
			//   tr.Trace << "after proton match " <<  proton.atom() << std::endl;
			if ( !proton_match ) continue;
			if ( !proton.has_connected_resonances() ) continue;
			if ( proton.first_connected_resonance()._pmatch( label_freq, label_error, label_folder ) <= 1.0 ) {
				//now store the representative of the float-group as match
				matches.push_back( std::make_pair( label(), first_connected_resonance().label() ) );
				return true; //both, proton and corresponding label matched.
			}
		} else { //label-matching
			Resonance const& float_label( (*res_list_)[ partner_id ] );
			//   tr.Trace << "try label float-match with reso: " <<  label.atom() << std::endl;
			//bool has_any_match = false;
			bool label_match = float_label._pmatch( label_freq, label_error, label_folder ) <= 1.0;
			if ( label_match ) {
				ResonanceAPs const& protons( connected_resonances() );
				for ( auto const & proton : protons ) {
					ResonanceOP r( proton );
					matches.push_back( std::make_pair( r->label(), label() ) );
				}
				return true;
			}
		} // end label-matching
	}
	return false; //nothing matched
}

core::Real FloatingResonance::pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const {
	Real min_match = 100000;
	if ( !is_representative_resonance_ ) return min_match;
	for ( core::Size partner_id : partner_ids_ ) {
		Real m = (*res_list_)[ partner_id ]._pmatch( peakfreq, error, folder );
		if ( m < min_match ) min_match = m;
	}
	return min_match;
}

} //NoesyAssign
} //devel
