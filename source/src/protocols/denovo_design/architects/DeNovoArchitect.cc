// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/DeNovoArchitect.cc
/// @brief Designs topologies
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.denovo_design.architects.DeNovoArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

std::string const
	DeNovoArchitect::DATA_MAP_NAME = "DeNovoArchitects";

DeNovoArchitect::DeNovoArchitect( std::string const & id ):
	StructureArchitect( id )
{}

DeNovoArchitect::~DeNovoArchitect() = default;

components::StructureDataOP
DeNovoArchitect::apply( core::pose::Pose const & pose ) const
{
	core::Real random = numeric::random::rg().uniform();
	return design( pose, random );
}

} //protocols
} //denovo_design
} //architects
