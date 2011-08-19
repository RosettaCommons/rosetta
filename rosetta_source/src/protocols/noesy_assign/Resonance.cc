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
#include <protocols/noesy_assign/Exceptions.hh>
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
#include <cstdlib>
#include <string>
#include <deque>


static basic::Tracer tr("protocols.noesy_assign.resonances");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

Resonance::Resonance() {}

Resonance::Resonance(  core::Size label, core::Real freq, core::Real error, core::id::NamedAtomID id) :
  label_ ( label ),
  freq_( freq ),
  error_( error ),
  atom_( id )
{
	calibration_atom_type_ = PeakCalibrator::atom_type( id );
}

Resonance::~Resonance() {}

void Resonance::write_to_stream( std::ostream& os ) const {
  os << ObjexxFCL::fmt::RJ( 10, label_ ) << " ";
  os << ObjexxFCL::fmt::F( 10, 3, freq_ ) << " " << ObjexxFCL::fmt::F( 10, 3, error_ ) << " ";
  os << ObjexxFCL::fmt::RJ( 5, atom_.atom() ) << " " << ObjexxFCL::fmt::RJ( 8, atom_.rsd() );
}

} //NoesyAssign
} //devel
