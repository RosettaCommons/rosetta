// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileRep.fwd.hh
/// @brief  Class/structure declarations for StructFileRep
/// @author Andy Watkins


#ifndef INCLUDED_core_io_StructFileRep_FWD_HH
#define INCLUDED_core_io_StructFileRep_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace io {

/// @brief  A structure for storing information from .pdb MODRES records.
struct ModifiedResidueInformation;

/// @brief  A structure for storing information from .pdb LINK records.
struct LinkInformation;

/// @brief  A structure for storing information from .pdb SSBOND records.
struct SSBondInformation;

/// @brief  A structure for storing information from .pdb CISPEP records.
struct CisPeptideInformation;

/// @brief  An Intermediate representation of data for ease of reading or writing a structural file.
class StructFileRep;

typedef utility::pointer::shared_ptr< StructFileRep > StructFileRepOP;
typedef utility::pointer::shared_ptr< StructFileRep const > StructFileRepCOP;

} // namespace io
} // namespace core

#endif // INCLUDED_core_io_StructFileRep_FWD_HH
