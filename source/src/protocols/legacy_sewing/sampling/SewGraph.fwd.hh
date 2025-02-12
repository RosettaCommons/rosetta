// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SewGraph.fwd.hh
///
/// @brief
/// @author Tim Jacobs



#ifndef INCLUDED_SewGraph_FWD_HH
#define INCLUDED_SewGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace legacy_sewing  {

class ModelNode;
class HashEdge;

class SewGraph;
typedef utility::pointer::shared_ptr< SewGraph > SewGraphOP;
typedef utility::pointer::shared_ptr< SewGraph const > SewGraphCOP;

} //legacy_sewing namespace
} //protocols namespace

#endif


