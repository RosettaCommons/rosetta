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
#include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>

// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/string_util.hh>
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


static thread_local basic::Tracer tr( "protocols.noesy_assign.assignment" );

PeakAssignment::PeakAssignment( CrossPeakAP cp, core::Size assign_spin1, core::Size assign_spin2 )
  : crosspeak_( cp ),
    spin_assign_index1_( assign_spin1 ),
    spin_assign_index2_( assign_spin2 )
{
  runtime_assert( cp );
  update_resonances_from_peak();
}

PeakAssignment::update_resonances_from_peak() {
  resonance1_ = crosspeak_->proton( 1 ).assignment( spin_id( 1 ) );
  resonance2_ = crosspeak_->proton( 2 ).assignment( spin_id( 2 ) );
}

