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
 /// @authorsv Christophe Schmitz //kalabharath
 /// //kalabharath
 /// @last_modified Aug 2011
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcsTs8_GridSearchIterator_fwd_hh
#define INCLUDED_protocols_scoring_methods_pcsTs8_GridSearchIterator_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs8{

class GridSearchIterator_Ts8;

typedef utility::pointer::owning_ptr< GridSearchIterator_Ts8 > GridSearchIterator_Ts8OP;
typedef utility::pointer::owning_ptr< GridSearchIterator_Ts8 const > GridSearchIterator_Ts8COP;

}//namespace pcsTs8
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
