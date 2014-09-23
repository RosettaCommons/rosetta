// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandDockingLoaders.cc
/// @brief  Implementation of the InterfaceBuilderLoader and MoveMapBuilderLoader classes
/// @author Gordon Lemmon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) -- moved here from DockDesignParser.cc

// Unit Headers
#include <protocols/ligand_docking/LigandDockingLoaders.hh>
#include <protocols/ligand_docking/LigandDockingLoaderCreators.hh>

// Project Headers
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>


#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer TR( "protocols.ligand_docking.LigandDockingLoaders" );

InterfaceBuilderLoader::InterfaceBuilderLoader() {}
InterfaceBuilderLoader::~InterfaceBuilderLoader() {}

void InterfaceBuilderLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	BOOST_FOREACH(TagCOP interface_builder_tag, tag->getTags()){
		std::string const name( interface_builder_tag->getName() );

		if ( data.has("interface_builders", name)) {
			TR << "WARNING WARNING movemap_builder of name \"" << name
				<< ") already exists. Skipping\n" << interface_builder_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::InterfaceBuilderOP interface_builder( new protocols::ligand_docking::InterfaceBuilder() );
		interface_builder->parse_my_tag( interface_builder_tag, data );
		data.add( "interface_builders" , name, interface_builder);
	}
	TR.flush();
}

jd2::parser::DataLoaderOP
InterfaceBuilderLoaderCreator::create_loader() const { return jd2::parser::DataLoaderOP( new InterfaceBuilderLoader ); }

std::string
InterfaceBuilderLoaderCreator::keyname() const { return "INTERFACE_BUILDERS"; }

MoveMapBuilderLoader::MoveMapBuilderLoader() {}
MoveMapBuilderLoader::~MoveMapBuilderLoader() {}

void MoveMapBuilderLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	BOOST_FOREACH(TagCOP movemap_builder_tag, tag->getTags()){
		std::string const name( movemap_builder_tag->getName() );

		if ( data.has("movemap_builders", name)) {
			TR << "WARNING WARNING movemap_builder of name \"" << name
				<< ") already exists. Skipping\n" << movemap_builder_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::MoveMapBuilderOP movemap_builder( new protocols::ligand_docking::MoveMapBuilder() );
		movemap_builder->parse_my_tag( movemap_builder_tag, data );
		data.add( "movemap_builders" , name, movemap_builder);
	}
	TR.flush();
}

jd2::parser::DataLoaderOP
MoveMapBuilderLoaderCreator::create_loader() const { return jd2::parser::DataLoaderOP( new MoveMapBuilderLoader ); }

std::string
MoveMapBuilderLoaderCreator::keyname() const { return "MOVEMAP_BUILDERS"; }

LigandAreaLoader::LigandAreaLoader() {}
LigandAreaLoader::~LigandAreaLoader() {}

void LigandAreaLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	BOOST_FOREACH(TagCOP ligand_area_tag, tag->getTags()){
		std::string const name( ligand_area_tag->getName() );

		if ( data.has("ligand_areas", name)) {
			TR << "WARNING WARNING ligand_area of name \"" << name
				<< ") already exists. Skipping\n" << ligand_area_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::LigandAreaOP ligand_area( new protocols::ligand_docking::LigandArea() );
		ligand_area->parse_my_tag( ligand_area_tag );
		data.add( "ligand_areas" , name, ligand_area);
	}
	TR.flush();
}

jd2::parser::DataLoaderOP
LigandAreaLoaderCreator::create_loader() const { return jd2::parser::DataLoaderOP( new LigandAreaLoader ); }

std::string
LigandAreaLoaderCreator::keyname() const { return "LIGAND_AREAS"; }



} //namespace jd2
} //namespace protocols
