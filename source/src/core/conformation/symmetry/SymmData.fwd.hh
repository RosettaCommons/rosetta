// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmData.fwd.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymmData_fwd_hh
#define INCLUDED_core_conformation_symmetry_SymmData_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace conformation {
namespace symmetry {

class SymmData;

typedef  utility::pointer::shared_ptr< SymmData >  SymmDataOP;
typedef  utility::pointer::shared_ptr< SymmData const >  SymmDataCOP;

} // namespace symmetry
} // namespace conformation
} // namespace core


#endif // INCLUDED_core_conformation_symmetry_SymmData_FWD_HH
