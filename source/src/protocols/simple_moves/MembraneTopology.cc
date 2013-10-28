// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

// Unit headers
#include <protocols/simple_moves/MembraneTopology.hh>
#include <protocols/simple_moves/MembraneTopologyCreator.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Project Headers
#include <core/pose/Pose.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <basic/Tracer.hh>

// Utility Headers

// Unit Headers

// C++ headers

namespace protocols {
namespace simple_moves {

using namespace std;

using core::pose::Pose;

static basic::Tracer TR( "protocols.simple_moves.MembraneTopology" );

std::string
MembraneTopologyCreator::keyname() const
{
	return MembraneTopologyCreator::mover_name();
}

protocols::moves::MoverOP
MembraneTopologyCreator::create_mover() const {
	return new MembraneTopology;
}

std::string
MembraneTopologyCreator::mover_name()
{
	return "MembraneTopology";
}

MembraneTopology::~MembraneTopology() {}

///@brief default ctor
MembraneTopology::MembraneTopology() :
	parent(),
	span_file_("")
{}

void MembraneTopology::apply( Pose & pose ) {
	core::scoring::MembraneTopologyOP membrane_topology = new core::scoring::MembraneTopology;
	membrane_topology->initialize( span_file_ );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, membrane_topology );
	TR<<"Setting pose's membrane topology according to span file "<<span_file()<<std::endl;
}

std::string
MembraneTopology::get_name() const {
	return MembraneTopologyCreator::mover_name();
}

void MembraneTopology::parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & )
{

	span_file( tag->getOption< std::string >( "span_file" ) );
	TR<<"Span file defined as "<<span_file()<<std::endl;
}


} // moves
} // protocols
