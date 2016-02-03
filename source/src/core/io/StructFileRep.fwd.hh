// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileRep.cc
/// @brief  The representation of a structure file
/// @author Andy Watkins

#ifndef INCLUDED_core_io_StructFileRep_fwd_hh
#define INCLUDED_core_io_StructFileRep_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace io {

class AtomInformation;
class ResidueInformation;

struct ModifiedResidueInformation;
struct LinkInformation;
struct SSBondInformation;
struct CisPeptideInformation;

class StructFileRep;

typedef utility::pointer::shared_ptr< StructFileRep > StructFileRepOP;
typedef utility::pointer::shared_ptr< StructFileRep const > StructFileRepCOP;


} // namespace io
} // namespace core


#endif // INCLUDED_core_io_StructFileRep_fwd_hh
