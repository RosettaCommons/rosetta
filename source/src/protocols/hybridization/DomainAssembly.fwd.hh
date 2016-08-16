// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TemlateAssembly.fwd.hh
/// @brief
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_TemlateAssembly_fwd_hh
#define INCLUDED_protocols_hybridization_TemlateAssembly_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
//namespace comparative_modeling {
namespace hybridization {
			
class TemlateAssembly;
typedef utility::pointer::shared_ptr< TemlateAssembly > TemlateAssemblyOP;
typedef utility::pointer::shared_ptr< TemlateAssembly const > TemlateAssemblyCOP;


} // hybridize 
//} // comparative_modeling 
} // protocols

#endif
