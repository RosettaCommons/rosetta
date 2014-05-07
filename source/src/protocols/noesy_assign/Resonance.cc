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
#include <protocols/noesy_assign/Resonance.hh>

// Package Headers
#include <protocols/noesy_assign/PeakCalibrator.hh>
#include <protocols/noesy_assign/FoldResonance.hh>

// AUTO-REMOVED #include <protocols/noesy_assign/Exceptions.hh>
#include <core/id/NamedAtomID.hh>
// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>

// #include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <deque>

#include <utility/vector1.hh>



static basic::Tracer tr("protocols.noesy_assign.resonances");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

Resonance::Resonance() {}

Resonance::Resonance(  core::Size label, core::Real freq, core::Real error, core::id::NamedAtomID const& id, core::chemical::AA aa, core::Real intensity ) :
  label_ ( label ),
  freq_( freq ),
  error_( error ),
  atom_( id ),
	aa_( aa ),
	intensity_( intensity )
{
	is_proton_ = ( id.atom()[ 0 ]=='Q' || ( id.atom().find("H") != std::string::npos && id.atom()[ 0 ] != 'C' /*not CH2 on TRP*/ ) );
	calibration_atom_type_ = PeakCalibrator::atom_type( id, aa_ );
}

Resonance::~Resonance() {}

core::Real Resonance::pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const {
	return _pmatch( peakfreq, error, folder );
}

core::Real Resonance::_pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const {
	return std::abs( folder( freq() ) - peakfreq ) / std::max( error, tolerance() );
}

void Resonance::_write_to_stream( std::ostream& os ) const {
  os << ObjexxFCL::format::RJ( 10, label_ ) << " ";
  os << ObjexxFCL::format::F( 10, 3, freq_ ) << " " << ObjexxFCL::format::F( 10, 3, error_ ) << " ";
  os << ObjexxFCL::format::RJ( 5, atom_.atom() ) << " " << ObjexxFCL::format::RJ( 8, atom_.rsd() );
}

void Resonance::write_to_stream( std::ostream& os ) const {
	_write_to_stream( os );
	//	os << " "<< intensity();
}

void Resonance::write_to_stream( std::ostream& os, core::chemical::AA aa ) const {
	_write_to_stream( os );
	os << " " << name_from_aa( aa ) << " " << oneletter_code_from_aa( aa );
	//	os << " "<< intensity();
}

void Resonance::combine( std::deque< ResonanceOP >& last_resonances, bool drain=false  ) {
	//is this a HB1, HG1 or a 1HG2 type of name ?
  bool single_char ( name().size() == 3 );
	std::string combine_name;
  Size str_cmp_length( 1 );

  std::string pseudo_letter( "Q" ); //default, single methyl group, proton
  if ( name().substr(0,1)=="C" ) pseudo_letter = "C"; //not default, we don't have a proton
	if ( name().substr(0,1)=="Q" ) pseudo_letter = "QQ"; //not default, this is a double-methyl...
	//if the two protons stuffed into the ambiguity queue are, .e.g, QG1 + QG2 --> QQG

	core::Real intensity_sum( intensity() );

//do we have enough protons for a QQX type of combined atom. (2 methyl groups... )
  bool double_methyl( last_resonances.size() == (6 - (drain ? 1 : 0))  );

// construct combined atom name, QX, QQX
	if ( single_char ) { 	//things like HB1+HB2 or QG1+QG2
    combine_name=pseudo_letter+name().substr(1,1);
  } else if ( double_methyl ) { //-->QQX
    combine_name=pseudo_letter+pseudo_letter+name().substr(1,1);
  } else { // HG11+HG12+HG13 --> QG1
    combine_name=pseudo_letter+name().substr(1,2);
    str_cmp_length=2;
  }

	//now figure out, how many atoms can be combined...
	Size limit( drain ? 0 : 1 ); //should we go to the very end, or leave last atom behind...
  while( last_resonances.size() > limit ) {//check others
    if ( name().substr( 1, str_cmp_length ) == last_resonances.front()->name().substr( 1, str_cmp_length ) ) {
			intensity_sum+=last_resonances.front()->intensity();
			last_resonances.pop_front();
		} else { //could not combine...
			combine_name = name();
			break;
		}
  } //now replace front of deque with the combined atom
	atom_ = core::id::NamedAtomID( combine_name, resid() );
	intensity_ = intensity_sum;
}

void Resonance::add_connected_resonance( ResonanceAP ptr ) {
	assert( ptr );
	connected_resonance_ids_.push_back( ptr->label() );
	connected_resonance_ptrs_.push_back( ptr );
}

void Resonance::clear_connected_resonances() {
	connected_resonance_ids_.clear();
	connected_resonance_ptrs_.clear();
}

Resonance const& Resonance::first_connected_resonance() const {
	runtime_assert( connected_resonance_ptrs_.size() );
	return *connected_resonance_ptrs_.front();
}

Resonance::ResonanceAPs const& Resonance::connected_resonances() const {
	return connected_resonance_ptrs_;
}


} //NoesyAssign
} //devel
