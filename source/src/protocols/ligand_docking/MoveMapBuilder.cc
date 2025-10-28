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
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@gmail.com) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <core/pose/chains_util.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/LigandDockingLoaders.hh>

//Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <utility>
#include <basic/Tracer.hh>
#include <core/types.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <utility/vector1.hh>

#include <basic/datacache/DataMap.hh> // AUTO IWYU For DataMap


//STL headers

namespace protocols {
namespace ligand_docking {

static basic::Tracer MoveMapBuilder_tracer( "protocols.ligand_docking.ligand_options.MoveMapBuilder", basic::t_debug );

void
set_jumps(
	core::pose::Pose const & pose, core::kinematics::MoveMapOP movemap,
	LigandAreas ligand_areas
){
	for ( LigandAreas::value_type const & ligand_area_pair : ligand_areas ) {
		std::string const & chain= ligand_area_pair.first;
		utility::vector1<core::Size> jump_ids= core::pose::get_jump_ids_from_chain(chain, pose);
		for ( core::Size const jump_id : jump_ids ) {
			movemap->set_jump(jump_id, true);
		}
	}
}

MoveMapBuilder::MoveMapBuilder():
	VirtualBase(),
	sc_interface_builder_(/* NULL */),
	bb_interface_builder_(/* NULL */),
	minimize_water_(false)
{}

MoveMapBuilder::MoveMapBuilder(MoveMapBuilder const & that):
	VirtualBase(),
	sc_interface_builder_(that.sc_interface_builder_),
	bb_interface_builder_(that.sc_interface_builder_),
	minimize_water_(that.minimize_water_)
{}

MoveMapBuilder::MoveMapBuilder(
	InterfaceBuilderOP sc,
	InterfaceBuilderOP bb,
	bool minimize_water
):
	VirtualBase(),
	sc_interface_builder_(std::move(sc)),
	bb_interface_builder_(std::move(bb)),
	minimize_water_(minimize_water)
{}


MoveMapBuilder::~MoveMapBuilder() = default;

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
MoveMapBuilder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
){
	if ( tag->hasOption("sc_interface") ) {
		std::string sc_interface_name= tag->getOption<std::string>("sc_interface");
		sc_interface_builder_= datamap.get_ptr<protocols::ligand_docking::InterfaceBuilder>( "interface_builders", sc_interface_name);
	}
	if ( tag->hasOption("bb_interface") ) {
		std::string bb_interface_name= tag->getOption<std::string>("bb_interface");
		bb_interface_builder_= datamap.get_ptr<protocols::ligand_docking::InterfaceBuilder>( "interface_builders", bb_interface_name);
	}

	if ( tag->hasOption("minimize_water") ) {
		if ( tag->getOption< bool >( "minimize_water" ) ) {
			minimize_water_= true;
		}
	}
}

core::kinematics::MoveMapOP
MoveMapBuilder::build(core::pose::Pose const & pose) const{
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );

	LigandAreas const & ligand_areas1 = sc_interface_builder_->get_ligand_areas();
	LigandAreas const & ligand_areas2 = sc_interface_builder_->get_ligand_areas();

	set_jumps(pose, movemap, ligand_areas1);
	set_jumps(pose, movemap, ligand_areas2);

	if ( sc_interface_builder_ ) set_all_chi(pose, movemap);
	if ( bb_interface_builder_ ) set_all_bb(pose, movemap);

	if ( minimize_water_ ) {
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( ! pose.residue(i).has_property("WATER") ) continue;
			core::kinematics::Edge const & e = pose.fold_tree().get_residue_edge(i);
			if ( ! e.is_jump() ) continue;
			movemap->set_jump( e.label(), true );
			MoveMapBuilder_tracer << "Minimize water jump " << e.label() << " to residue " << i << " " << pose.residue_type(i).name3() << std::endl;
		}
	}
	return movemap;
}

InterfaceBuilderOP
MoveMapBuilder::get_sc_interface_builder()const{
	debug_assert( sc_interface_builder_); // does the pointer point
	return sc_interface_builder_;
}

InterfaceBuilderOP
MoveMapBuilder::get_bb_interface_builder()const{
	debug_assert( bb_interface_builder_); // does the pointer point
	return bb_interface_builder_;
}

void
MoveMapBuilder::set_all_chi(
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP movemap
)const{
	ligand_options::Interface side_chain_interface= sc_interface_builder_->build(pose);
	MoveMapBuilder_tracer.Debug<< "moveMap interface: "<< side_chain_interface << std::endl;
	for ( core::Size i=1; i <= side_chain_interface.size(); ++i ) {
		if ( side_chain_interface[i].type != ligand_options::InterfaceInfo::non_interface ) { // Allow residue to minimize
			movemap->set_chi(i, true);
		}
	}
	// remove ligands with a minimize_ligand value of zero
	for ( core::Size chain_id = 1; chain_id <= pose.conformation().num_chains(); ++chain_id ) {
		std::string const chain= core::pose::get_chain_from_chain_id(chain_id, pose);
		LigandAreas const & ligand_areas= sc_interface_builder_->get_ligand_areas();
		auto ligand_area= ligand_areas.find(chain);
		if ( ligand_area == ligand_areas.end() ) continue;
		if ( ligand_area->second->minimize_ligand_ > 0 ) continue;
		// else if we aren't minimizing the ligand set_chi for this ligand to false
		core::Size begin = pose.conformation().chain_begin( chain_id);
		core::Size const end = pose.conformation().chain_end( chain_id );
		for ( ; begin <= end; ++begin ) movemap->set_chi(begin, false);
	}
}

/// @details You MUST call set_all_chi first
void
MoveMapBuilder::set_all_bb(
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP movemap
)const{
	ligand_options::Interface bb_interface = bb_interface_builder_->build(pose);

	for ( core::Size i=1; i <= bb_interface.size(); ++i ) {
		// bb allowed only if sc is allowed.
		if ( bb_interface[i].type != ligand_options::InterfaceInfo::non_interface
				&& pose.residue(i).is_protein()
				&& movemap->get_chi(i)
				) { // Allow residue to minimize
			movemap->set_bb(i, true);
		}
	}
}

bool
MoveMapBuilder::minimize_backbone(){
	return (bool) bb_interface_builder_;
}

std::string MoveMapBuilder::element_name() { return "MoveMapBuilder"; }
void MoveMapBuilder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "sc_interface", xs_string, "The name of the previously declared"
		" interface builder to use when deciding what sidechains to allow to move when building a MoveMap" )
		+ XMLSchemaAttribute( "bb_interface", xs_string, "The name of the previously declared interface"
		" builder to use when deciding what backbone torsions to allow to move when building a MoveMap" )
		+ XMLSchemaAttribute( "minimize_water", xsct_rosetta_bool, "Should water residues be allowed to move?" )
		+ required_name_attribute();

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( element_name() )
		.complex_type_naming_func( & MoveMapBuilderLoader::movemap_builder_ct_namer )
		.description( "Movemap builders can build on previousy-defined InterfaceBuilders to define the"
		" flexible sidechain and backbone dihedrals by creating a MoveMap." )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}

} //namespace ligand_docking
} //namespace protocols
