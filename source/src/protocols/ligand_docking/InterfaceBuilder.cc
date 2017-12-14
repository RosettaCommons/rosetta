// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/InterfaceBuilder.hh>

// Package headers
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <protocols/ligand_docking/ligand_options/interface_distance_functions.hh>
#include <protocols/ligand_docking/LigandDockingLoaders.hh>
#include <protocols/ligand_docking/LigandArea.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>

namespace protocols {
namespace ligand_docking {

static basic::Tracer interface_builder_tracer( "protocols.ligand_docking.ligand_options.InterfaceBuilder", basic::t_debug );

InterfaceBuilder::InterfaceBuilder():
	ReferenceCount(),
	ligand_areas_(),
	extension_window_(0)
{}

InterfaceBuilder::InterfaceBuilder(utility::vector1<LigandAreaOP> ligand_areas, core::Size extension_window):
	ReferenceCount(),
	extension_window_(extension_window)
{
	for ( LigandAreaOP ligand_area : ligand_areas ) {
		ligand_areas_[ligand_area->chain_] = ligand_area;
	}
}

InterfaceBuilder::InterfaceBuilder(InterfaceBuilder const & that):
	ReferenceCount(),
	ligand_areas_(that.ligand_areas_),
	extension_window_(that.extension_window_)
{}

InterfaceBuilder::~InterfaceBuilder() = default;

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
InterfaceBuilder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
){
	if ( tag->hasOption("extension_window") ) {
		extension_window_= tag->getOption<core::Size>("extension_window");
	}
	if ( ! tag->hasOption("ligand_areas") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "InterfaceBuilders need at least one ligand area to build around");

	std::string ligand_areas_string= tag->getOption<std::string>("ligand_areas");
	utility::vector1<std::string> ligand_area_strings= utility::string_split(ligand_areas_string, ',');

	for ( std::string const & ligand_area_string : ligand_area_strings ) {
		LigandAreaOP ligand_area = datamap.get_ptr< protocols::ligand_docking::LigandArea >( "ligand_areas", ligand_area_string);
		ligand_areas_[ ligand_area->chain_ ] = ligand_area;
	}
}

ligand_options::Interface InterfaceBuilder::build(core::pose::Pose const & pose) const{
	ligand_options::Interface interface( pose.size(), ligand_options::InterfaceInfo() ); // init all positions to false
	find_interface_residues(interface, pose);
	if ( extension_window_ > 0 ) enforce_minimum_length(interface, pose);
	interface_builder_tracer.Debug<< "built interface: "<< interface<< std::endl;
	return interface;
}

void InterfaceBuilder::find_interface_residues(
	ligand_options::Interface & interface,
	core::pose::Pose const & pose
)const{
	auto ligand_area= ligand_areas_.begin();
	for ( ; ligand_area != ligand_areas_.end(); ++ligand_area ) {
		char const & chain= ligand_area->first;
		utility::vector1<core::Size> chain_ids= core::pose::get_chain_ids_from_chain(chain, pose);
		for ( core::Size const chain_id : chain_ids ) {
			core::Size ligand_residue_id= pose.conformation().chain_begin(chain_id);
			core::Size const & end= pose.conformation().chain_end(chain_id);
			for ( ; ligand_residue_id <= end; ++ligand_residue_id ) {
				interface[ligand_residue_id].type=  ligand_options::InterfaceInfo::is_interface;
				interface[ligand_residue_id].chain_id= chain_id;
				interface[ligand_residue_id].chain= chain;
				find_protein_residues(interface, ligand_residue_id, pose);
			}
		}
	}
	interface_builder_tracer.Debug << "interface: " << interface << std::endl;
}

// First call "find_ligand_residues"
void InterfaceBuilder::find_protein_residues(
	ligand_options::Interface & interface,
	core::Size ligand_residue_id,
	core::pose::Pose const & pose
)const{
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( pose.residue(i).is_protein() &&
				interface[i].type != ligand_options::InterfaceInfo::is_interface ) {
			set_interface_residue(interface, i, ligand_residue_id, pose);
		}
	}
}

void InterfaceBuilder::set_interface_residue(
	ligand_options::Interface & interface,
	core::Size const potential_interface_residue_id,
	core::Size const ligand_interface_residue_id,
	core::pose::Pose const & pose
)const{
	core::conformation::Residue const & potential_interface_residue = pose.residue(potential_interface_residue_id);
	core::conformation::Residue const & ligand_residue = pose.residue(ligand_interface_residue_id);

	core::Size const & chain_id= interface[ligand_interface_residue_id].chain_id;
	core::Size const & chain= interface[ligand_interface_residue_id].chain;

	if ( is_interface_residue(potential_interface_residue, ligand_residue, chain) ) {
		interface[potential_interface_residue_id].type= ligand_options::InterfaceInfo::is_interface;
		interface[potential_interface_residue_id].chain_id= chain_id;
		interface[potential_interface_residue_id].chain= chain;
	}
}

bool InterfaceBuilder::is_interface_residue(
	core::conformation::Residue const & potential_interface_residue,
	core::conformation::Residue const & ligand_interface_residue,
	char const chain
)const {
	core::Size const potential_interface_neighbor_atom_id= potential_interface_residue.nbr_atom();
	core::Vector const potential_interface_vector= potential_interface_residue.xyz(potential_interface_neighbor_atom_id);

	auto found= ligand_areas_.find(chain);
	debug_assert(found != ligand_areas_.end());
	LigandAreaOP const ligand_area = found->second;

	double cutoff = ligand_area->add_nbr_radius_ ?
		potential_interface_residue.nbr_radius() + ligand_area->cutoff_ :
		ligand_area->cutoff_;

	if ( ligand_area->all_atom_mode_ ) {
		return ligand_options::check_all_ligand_atoms(ligand_interface_residue, potential_interface_vector, cutoff);
	} else {
		return ligand_options::check_neighbor_ligand_atom(ligand_interface_residue, potential_interface_vector, cutoff);
	}
}

void InterfaceBuilder::enforce_minimum_length(
	ligand_options::Interface & interface,
	core::pose::Pose const & pose
) const{
	for ( core::Size chain_id=1; chain_id <= pose.conformation().num_chains(); ++chain_id ) {
		core::Size start= pose.conformation().chain_begin(chain_id);
		core::Size end= pose.conformation().chain_end(chain_id);
		if ( pose.residue(start).is_polymer() ) { // because ligand residues don't know about terminus
			interface.enforce_minimum_length(start, end, extension_window_);
		}
	}
}

LigandAreas
InterfaceBuilder::get_ligand_areas() const{
	return ligand_areas_;
}

std::string InterfaceBuilder::element_name() { return "InterfaceBuilder"; }
void InterfaceBuilder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using Attr = XMLSchemaAttribute;
	AttributeList attributes;
	attributes + Attr( "extension_window", xsct_non_negative_integer, "XRW TO DO" )
		+ Attr::required_attribute( "ligand_areas", xs_string, "The comma-separated list of ligand areas, which are read out of the datamap" )
		+ required_name_attribute();

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( element_name() )
		.complex_type_naming_func( & InterfaceBuilderLoader::interface_builder_ct_namer )
		.description( "Each interface builder is composed from one or more ligand areas, which must have been previously defined" )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}

} //namespace ligand_docking
} //namespace protocols
