// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NoesyModule.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_CyanaModule_hh
#define INCLUDED_protocols_noesy_assign_CyanaModule_hh


// Unit Header
#include <protocols/noesy_assign/NoesyModule.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/PeakFileFormat.fwd.hh>
#include <protocols/noesy_assign/PeakAssignmentList.hh>
//#include <devel/NoesyAssign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
// #include <cstdlib>
// #include <string>
// #include <list>
// #include <map>

namespace protocols {
namespace noesy_assign {

class NoesyModule : public utility::pointer::ReferenceCount {
public:
  NoesyModule( std::string const& fasta_sequence );
  void assign( core::Size cycle );
  void generate_constraint_files( core::pose::Pose const& pose, std::string const& cst_fa_file, std::string const& cst_centroid_file ) const;

private:
  CrossPeakListOP crosspeaks_;
  ResonanceListOP resonances_;
};

}
}

#endif
