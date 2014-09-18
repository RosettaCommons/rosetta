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
#include <protocols/noesy_assign/NoesyModule.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>
#include <protocols/noesy_assign/DistanceScoreMover.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.hh>

// for switching residue type set to centroid

#include <core/chemical/ChemicalManager.fwd.hh>

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers
#include <math.h> //for isnan

static thread_local basic::Tracer tr( "protocols.noesy_assign.NoesyModule" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

NoesyModule::NoesyModule( std::string const& fasta_sequence ) :
  crosspeaks_( NULL ),
  resonances_( new ResonanceList( fasta_sequence ) )
{
  utility::io::izstream input_file(basic::options::option[ basic::options::OptionKeys::r ]() );
  utility::io::ozstream output_file( basic::options::option[ basic::options::OptionKeys::ro ]() );
  if ( input_file.good() ) {
    resonances->read_from_stream( input_file );
    resonances->write_to_stream( output_file );
  } else {
    tr.Error << "cannot read " << basic::options::option[ basic::options::OptionKeys::r ]() << std::endl;
  }

}

} // namespace noesy_assign
} // namespace protocols
