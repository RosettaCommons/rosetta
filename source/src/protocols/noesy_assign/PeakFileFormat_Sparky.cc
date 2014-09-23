// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/PeakFileFormat_Sparky.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakInfo.hh>
#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/noesy_assign/Exceptions.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// AUTO-REMOVED #include <utility/string_util.hh>
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

// AUTO-REMOVED #include <basic/options/option_macros.hh>

#include <protocols/noesy_assign/PeakAssignment.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.noesy_assign.io" );

namespace protocols {
namespace noesy_assign {

using namespace core;

//PeakFileFormat_Sparky::PeakFileFormat_Sparky( ResonanceListOP const& res )  :
//	PeakFileFormat( res )
//{}

void PeakFileFormat_Sparky::write_peak( std::ostream& os, Size ct, CrossPeak const& cp ) const {
	std::ostringstream line_end;
	line_end << " #d " << ObjexxFCL::format::F( 6, 3, cp.distance_bound() );
	if ( cp.eliminated( false /*recompute*/, true /*do_not_compute*/) ) line_end << " #eliminated: " << cp.elimination_reason();


  // cp.write_to_stream( os );
  write_assignments( os, cp, "" );
  write_resonances( os, cp );
	os << "# Peak " << ObjexxFCL::format::RJ( 6, ct ) << " ";
	//  write_strength( os, cp );

	os << line_end.str();


}

void PeakFileFormat_Sparky::set_format_from_peak( CrossPeak const& cp ) {
  info1_ = CrossPeakInfoOP( new CrossPeakInfo( cp.info( 1 ) ) );
  info2_ = CrossPeakInfoOP( new CrossPeakInfo( cp.info( 2 ) ) );
  col2proton_.clear();
  col2islabel_.clear();

	//dimension 2 - label
  if ( info2_->has_label() ) {
    col2proton_.push_back( 2 );
    col2islabel_.push_back( true );
  }

  //dimension 2
  col2proton_.push_back( 2 );
  col2islabel_.push_back( false );

	//dimension 1 - label
  if ( info1_->has_label() ) {
    col2proton_.push_back( 1 );
    col2islabel_.push_back( true );
  }

  //dimension 1
  col2proton_.push_back( 1 );
  col2islabel_.push_back( false );

}

void PeakFileFormat_Sparky::write_header( std::ostream& ) {}

// void PeakFileFormat_Sparky::write_resonances( std::ostream& os, CrossPeak const& cp ) const {
//   //  cp.write_to_stream( os );
//   Size const ncol( col2proton_.size() );
//   runtime_assert( col2islabel_.size() == ncol );

//   for ( Size icol=1; icol<=ncol; ++icol ) {
//     Real val;
//     Size iproton( col2proton_[ icol ] );
//     bool is_label( col2islabel_[ icol ] );
//     if ( !is_label ) val = cp.proton( iproton ).freq();
//     else {
//       runtime_assert( cp.has_label( iproton ) );
//       val = cp.label( iproton ).freq();
//     }
//     os << ObjexxFCL::format::F( 8, 3, val ) << " ";
//   }
// }


// void PeakFileFormat_Sparky::write_strength( std::ostream& os, CrossPeak const& cp ) const {
//   os << ObjexxFCL::format::E( 10, 3, cp.volume() ) << " " << ObjexxFCL::format::E( 10, 3, 0.0 ) << " ";
// }

void PeakFileFormat_Sparky::write_assignment_indent( std::ostream& os, CrossPeak const&) const {
	os << std::endl;// << "                                                        ";
	//	if ( cp.has_label( 1 ) && cp.has_label( 2 ) ) os << "         ";
}

void PeakFileFormat_Sparky::write_nil_assignment( std::ostream& os ) const {
	os << ObjexxFCL::format::RJ( 25, "?-?-?" ) << " ";
}

void PeakFileFormat_Sparky::write_assignment( std::ostream& os, PeakAssignment const& pa ) const {
	Size resid( 0 );
	std::ostringstream buf;
	for ( Size icol=1; icol<=ncol(); ++icol ) {
		Size val;
		Size iproton( col2proton_[ icol ] );
		bool is_label( col2islabel_[ icol ] );
		if ( write_atom_names() ) {
			core::id::NamedAtomID atom;
			if ( !is_label ) atom = pa.atom( iproton );
			else atom = pa.label_atom( iproton );
			if ( resid != atom.rsd() ) {
				resid = atom.rsd();
				core::chemical::AA aa( pa.resonances().aa_from_resid( atom.rsd() ));
				buf << oneletter_code_from_aa(aa) << resid;
			}
			buf << atom.atom();
		} else {
			if ( !is_label ) val = pa.resonance_id( iproton );
			else {
				val = pa.label_resonance_id( iproton );
			}
			buf << ObjexxFCL::format::RJ( 6, val );
		}
		if ( icol < ncol() ) buf << "-";
	}
	os << ObjexxFCL::format::RJ( 25, buf.str() ) << " ";
}


}
}
