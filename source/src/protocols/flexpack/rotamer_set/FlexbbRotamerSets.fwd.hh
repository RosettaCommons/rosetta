// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/
/// @brief
/// @author Florian Richter (floric@u.washington.edu), sep 08

#ifndef INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSets_fwd_hh
#define INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSets_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

//#include <utility/pointer/all.fwd.hh>

namespace protocols {
namespace flexpack {
namespace rotamer_set {

class FlexbbRotamerSets;

typedef utility::pointer::shared_ptr< FlexbbRotamerSets > FlexbbRotamerSetsOP;
typedef utility::pointer::shared_ptr< FlexbbRotamerSets const > FlexbbRotamerSetsCOP;
typedef utility::pointer::weak_ptr< FlexbbRotamerSets > FlexbbRotamerSetsAP;
typedef utility::pointer::weak_ptr< FlexbbRotamerSets const > FlexbbRotamerSetsCAP;


} //namespace rotamer_set
} //namespace flexpack
} //namespace protocols

#endif
