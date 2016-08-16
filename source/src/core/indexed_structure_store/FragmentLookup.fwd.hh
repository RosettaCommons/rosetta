// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#ifndef INCLUDED_core_indexed_structure_store_FragmentLookup_fwd_hh
#define INCLUDED_core_indexed_structure_store_FragmentLookup_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core
{
namespace indexed_structure_store
{

class FragmentLookupResult;

class FragmentLookup;

typedef utility::pointer::shared_ptr<FragmentLookup> FragmentLookupOP;
typedef utility::pointer::shared_ptr<FragmentLookup const> FragmentLookupCOP;
typedef utility::pointer::weak_ptr<FragmentLookup> FragmentLookupAP;
typedef utility::pointer::weak_ptr<FragmentLookup const> FragmentLookupCAP;

}
}

#endif
