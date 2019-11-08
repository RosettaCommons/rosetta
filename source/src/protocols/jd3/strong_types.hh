// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/strong_types.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_jd3_standard_StrongTypes_HH
#define INCLUDED_protocols_jd3_standard_StrongTypes_HH

//utility headers
#include <utility/strong_aliasing.hh>

namespace protocols {
namespace jd3 {

using JobDAGNodeID = utility::StrongSize< struct JobDAGNodeID_ >;

//preliminary job node (used in standard/)
using PrelimJobNodeID = utility::StrongSize< struct PrelimJobNodeID_ >;

//The ID for a job in the context of the entire JD3 run
using GlobalJobID = utility::StrongSize< struct GlobalJobID_ >;

//The ID for a job in the context of its individual stage
using LocalJobID = utility::StrongSize< struct LocalJobID_ >;

using NStructIndex = utility::StrongSize< struct NStructIndex_ >;

using ResultIndex = utility::StrongSize< struct ResultIndex_ >;

#ifdef SERIALIZATION
SERIALIZE_STRONG_SIZE_HH( JobDAGNodeID );
SERIALIZE_STRONG_SIZE_HH( PrelimJobNodeID );
SERIALIZE_STRONG_SIZE_HH( GlobalJobID );
SERIALIZE_STRONG_SIZE_HH( LocalJobID );
SERIALIZE_STRONG_SIZE_HH( NStructIndex );
SERIALIZE_STRONG_SIZE_HH( ResultIndex );
#endif

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_standard_StrongTypes_HH
