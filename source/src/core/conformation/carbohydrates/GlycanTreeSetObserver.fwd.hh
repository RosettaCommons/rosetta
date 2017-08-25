// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/carbohydrates/GlycanTreeSetObserver.fwd.hh
/// @brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_fwd_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace core {
namespace conformation {
namespace carbohydrates {

class GlycanTreeSetObserver;

typedef utility::pointer::shared_ptr< GlycanTreeSetObserver > GlycanTreeSetObserverOP;
typedef utility::pointer::shared_ptr< GlycanTreeSetObserver const > GlycanTreeSetObserverCOP;



} //core
} //conformation
} //carbohydrates


#endif //INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_fwd_hh





