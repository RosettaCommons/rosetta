// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextGraph.cc
/// @brief  Context graph class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/scoring/ContextGraph.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

ContextGraph::~ContextGraph() = default;


ContextGraph::ContextGraph()
:
	parent()
{}

ContextGraph::ContextGraph(Size num_nodes)
:
	parent( num_nodes )
{}

ContextGraph::ContextGraph( ContextGraph const & )
:
	parent()
{}

ContextGraph &
ContextGraph::operator = ( ContextGraph const & source ) {
	return static_cast< ContextGraph & > ( parent::operator = ( source ) );
}


Size ContextGraph::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}


} // scoring
} // core

