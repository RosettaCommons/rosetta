// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
