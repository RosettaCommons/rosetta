// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_memb_etable_MembEtable_fwd_hh
#define INCLUDED_core_scoring_memb_etable_MembEtable_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace etable {

class MembEtable;

typedef utility::pointer::shared_ptr< MembEtable > MembEtableOP;
	//typedef utility::pointer::owning_ptr< Etable const > EtableCOP;
typedef utility::pointer::weak_ptr< MembEtable const > MembEtableCAP;

} // etable
} // scoring
} // core

#endif
