// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Class to write kinemage-formatted output for ResidueType
/// @file   core/chemical/ResidueTypeKinWriter.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_ResidueKinWriter_fwd_hh
#define INCLUDED_core_chemical_ResidueKinWriter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {

class ResidueTypeKinWriter;

typedef utility::pointer::owning_ptr< ResidueTypeKinWriter > ResidueTypeKinWriterOP;
typedef utility::pointer::owning_ptr< ResidueTypeKinWriter const > ResidueTypeKinWriterCOP;

} // chemical
} // core


#endif
