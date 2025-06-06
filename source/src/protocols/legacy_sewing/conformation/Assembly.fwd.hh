// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/conformation/Assembly.fwd.hh
///
/// @brief
/// @author Tim Jacobs



#ifndef INCLUDED_protocols_legacy_sewing_conformation_Assembly_FWD_HH
#define INCLUDED_protocols_legacy_sewing_conformation_Assembly_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace legacy_sewing  {

class Assembly;
typedef utility::pointer::shared_ptr< Assembly > AssemblyOP;
typedef utility::pointer::shared_ptr< Assembly const > AssemblyCOP;

} //legacy_sewing namespace
} //protocols namespace

#endif


