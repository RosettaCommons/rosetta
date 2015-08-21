// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpSampleData.fwd.hh
/// @brief definition of the JumpSampleData class
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_JumpSampleData_fwd_hh
#define INCLUDED_protocols_abinitio_abscript_JumpSampleData_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace abinitio {
namespace abscript {

class JumpSampleData;
typedef utility::pointer::shared_ptr< JumpSampleData > JumpSampleDataOP;
typedef utility::pointer::shared_ptr< JumpSampleData const > JumpSampleDataCOP;

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_moves_JumpSampleData_fwd_HH
