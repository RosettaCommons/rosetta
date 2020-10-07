// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sap/SapDatabase.fwd.hh
/// @brief Hdf5 version of rosetta database
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_sap_SapDatabase_fwd_hh
#define INCLUDED_core_sap_SapDatabase_fwd_hh

// Unit Headers

// Package headers

// Project headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

// Forward
class SapDatabase;

// Types
typedef utility::pointer::shared_ptr< SapDatabase > SapDatabaseOP;
typedef utility::pointer::shared_ptr< SapDatabase const > SapDatabaseCOP;


} //sap
} //guidance_scoreterms
} //pack
} //core
#endif
