// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody2/Ab_Info.fwd.hh
/// @brief  Ab_Info class forward declarations header
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_Ab_Info_fwd_hh
#define INCLUDED_protocols_antibody2_Ab_Info_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


// C++ Headers
namespace protocols{
namespace antibody2{

// Forward
class Ab_Info;

typedef utility::pointer::owning_ptr< Ab_Info > Ab_InfoOP;
typedef utility::pointer::owning_ptr< Ab_Info const > Ab_InfoCOP;



} //namespace antibody2
} //namespace protocols

#endif //INCLUDED_protocols_Ab_Info_FWD_HH
