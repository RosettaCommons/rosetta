// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
 //////////////////////////////////////////////
 /// @file GridSearchIterator.fwd.hh
 ///
 /// @authorsv Christophe Schmitz , Kala Bharath Pilla
 /// ,
 /// @last_modified Mar 2014
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcsTs3_GridSearchIterator_fwd_hh
#define INCLUDED_protocols_scoring_methods_pcsTs3_GridSearchIterator_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs3{

class GridSearchIterator_Ts3;

typedef utility::pointer::owning_ptr< GridSearchIterator_Ts3 > GridSearchIterator_Ts3OP;
typedef utility::pointer::owning_ptr< GridSearchIterator_Ts3 const > GridSearchIterator_Ts3COP;

}//namespace pcsTs3
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
