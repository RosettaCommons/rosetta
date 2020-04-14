// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ConcatenatePosesMover.cc
/// @brief links supplied Poses by their termini
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/simple_moves/ConcatenatePosesMover.hh>
#include <protocols/simple_moves/ConcatenatePosesMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <core/conformation/Atom.hh>
#include <numeric/random/random.hh>
#include <numeric/HomogeneousTransform.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <algorithm>
#include <regex>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.ConcatenatePosesMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ConcatenatePosesMover::ConcatenatePosesMover():
	protocols::moves::Mover( ConcatenatePosesMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ConcatenatePosesMover::ConcatenatePosesMover( ConcatenatePosesMover const & src ):
	protocols::moves::Mover( src )
{
	component_file_ = src.get_component_file();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ConcatenatePosesMover::~ConcatenatePosesMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ConcatenatePosesMover::apply( core::pose::Pose& pose){
	this->concatenate_poses(pose);
}

void
ConcatenatePosesMover::concatenate_poses(core::pose::Pose& pose){

	core::pose::PoseOP current_pose;
	core::pose::PoseOP last_pose;
	utility::vector1<core::pose::PoseOP > poses;
	utility::io::izstream pdb_source_file( component_file_ );
	utility::vector1<std::pair<core::Size,core::Size>> fixed_spans;
	std::string line;
	std::string fold_tree = "FOLD_TREE ";
	utility::vector1<std::string> tokens;
	std::string new_sequence;
	core::Size current_chain;
	std::string output_prefix;
	numeric::HomogeneousTransform< core::Real > immobile_ht;
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	utility::vector1<std::tuple<std::string,core::Size,core::pose::PoseOP>> new_chain_components; // new component name, chain, pointer
	utility::vector1<std::tuple<std::string,core::Size,core::pose::PoseOP>> binder_pairs; // new component name, chain, binder pointer
	utility::vector1<std::pair<core::Size,core::Size>> binder_root_residues; // residue in new component, binder number

	bool good_file = false;
	while ( getline( pdb_source_file, line) ) {
		tokens = utility::string_split(line,' ');
		if ( tokens[1]=="DOMAIN" ) {
			good_file=true;
			std::pair<core::Size,core::Size> new_span;
			current_pose = core::import_pose::pose_from_file(tokens[2]);
			current_chain = std::stoi(tokens[3]);
			current_pose = current_pose->split_by_chain(current_chain);
			new_chain_components.push_back(std::make_tuple(tokens[2],current_chain,current_pose));
			new_span.first = new_sequence.length()+1;
			new_sequence = new_sequence + std::get<2>(new_chain_components.back())->sequence();
			new_span.second = new_sequence.length();
			fixed_spans.push_back(new_span);
		} else if ( tokens[1]=="LINKER" ) {
			good_file=true;
			current_pose = core::pose::PoseOP(new core::pose::Pose());
			core::pose::make_pose_from_sequence(*current_pose, tokens[2], *res_type_set,false);
			for ( core::Size i = 1; i <= current_pose->size(); i++ ) {
				current_pose->set_phi(i,-180);
				current_pose->set_psi(i,180);
				current_pose->set_omega(i,180);
			}
			new_chain_components.push_back(std::make_tuple("LINKER",1,current_pose));
			new_sequence = new_sequence + std::get<2>(new_chain_components.back())->sequence();
		} else if ( tokens[1]=="BINDER" ) {
			current_pose = core::import_pose::pose_from_file(tokens[4]);
			current_chain = std::stoi(tokens[5]);
			current_pose = current_pose->split_by_chain(current_chain);
			binder_pairs.push_back(std::make_tuple(tokens[2],std::stoi(tokens[3]),current_pose));
		} else if ( tokens[1]!="IGNORE" ) {
			TR << "Bad token: " << tokens[1] << std::endl;
			return;
		}
	}
	if ( !good_file ) {
		TR << "nonexistent or malformed component file: " << component_file_ <<std::endl;
		return;
	}



	for ( core::Size current_pose_number = 1; current_pose_number <= new_chain_components.size(); current_pose_number++ ) {
		current_pose = std::get<2>(new_chain_components[current_pose_number]);
		for ( core::Size current_residue_number = 1; current_residue_number<=current_pose->size(); current_residue_number++ ) {
			if ( current_pose->residue(current_residue_number).is_upper_terminus() ) {
				core::pose::remove_variant_type_from_pose_residue(*current_pose, core::chemical::UPPER_TERMINUS_VARIANT, current_residue_number);
			}
			if ( current_pose->residue(current_residue_number).is_lower_terminus() ) {
				core::pose::remove_variant_type_from_pose_residue(*current_pose, core::chemical::LOWER_TERMINUS_VARIANT, current_residue_number);
			}
		}
	}

	core::pose::PoseOP working_pose = core::pose::PoseOP(new core::pose::Pose(*std::get<2>(new_chain_components[1])));
	for ( core::Size current_binder = 1; current_binder <= binder_pairs.size(); current_binder++ ) {
		if ( std::get<0>(new_chain_components[1])==std::get<0>(binder_pairs[current_binder]) && std::get<1>(new_chain_components[1])==std::get<1>(binder_pairs[current_binder]) ) {
			binder_root_residues.push_back(std::make_pair(working_pose->size()-1,current_binder));
		}
	}

	for ( core::Size current_pose_number = 2; current_pose_number <= new_chain_components.size(); current_pose_number++ ) {
		last_pose = std::get<2>(new_chain_components[current_pose_number-1]);
		current_pose = std::get<2>(new_chain_components[current_pose_number]);
		last_pose->append_residue_by_bond(current_pose->residue(1),true,0,0,0,false,false);

		core::conformation::Residue stationary_basis_residue = last_pose->residue(last_pose->size());
		utility::vector1< core::conformation::Atom > & stationary_basis_atoms = stationary_basis_residue.atoms();
		core::conformation::Residue mobile_basis_residue = current_pose->residue(1);
		utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue.atoms();


		numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
		numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;
		current_pose->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());
		working_pose->append_pose_by_jump(*current_pose,working_pose->size(),"CA","CA");
		for ( core::Size current_binder = 1; current_binder <= binder_pairs.size(); current_binder++ ) {
			if ( std::get<0>(new_chain_components[current_pose_number])==std::get<0>(binder_pairs[current_binder]) && std::get<1>(new_chain_components[current_pose_number])==std::get<1>(binder_pairs[current_binder]) ) {
				current_pose = std::get<2>(binder_pairs[current_binder]);
				current_pose->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());
				binder_root_residues.push_back(std::make_pair(working_pose->size()-1,current_binder));
			}
		}
	}
	core::pose::remove_variant_type_from_pose_residue(*working_pose, core::chemical::UPPER_TERMINUS_VARIANT, 1);
	core::pose::remove_variant_type_from_pose_residue(*working_pose, core::chemical::LOWER_TERMINUS_VARIANT, working_pose->size());
	fold_tree += "EDGE 1 " + std::to_string(working_pose->size()) + " -1 ";


	core::pose::PoseOP output_pose = core::pose::PoseOP(new core::pose::Pose);

	core::pose::make_pose_from_sequence(*output_pose, working_pose->sequence(), *res_type_set,false);

	output_pose->copy_segment(output_pose->size(),*working_pose,1,1);

	for ( auto span : fixed_spans ) {
		TR << span.first << " " << span.second << std::endl;
	}

	for ( core::Size current_binder = 1; current_binder <= binder_root_residues.size(); current_binder++ ) {
		std::pair<core::Size,core::Size> current_binder_index = binder_root_residues[current_binder];
		current_pose = std::get<2>(binder_pairs[current_binder_index.second]);
		core::pose::append_pose_to_pose(*output_pose,*current_pose,true);
		TR << std::to_string(output_pose->size()-current_pose->size()+1)  << " " << output_pose->size() << std::endl;
		fold_tree += "EDGE " + std::to_string(current_binder_index.first) + " " + std::to_string(output_pose->size()-current_pose->size()+1) + " " + std::to_string(current_binder) + " ";
		fold_tree += "EDGE " + std::to_string(output_pose->size()-current_pose->size()+1) + " " + std::to_string(output_pose->size()) + " -1 ";
	}


	TR << fold_tree << std::endl;

	pose = *output_pose;
}
////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ConcatenatePosesMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ConcatenatePosesMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
) {
	if ( tag->hasOption( "component_file" ) ) {
		component_file_ = tag->getOption< std::string >( "component_file") ;
	}

}
void ConcatenatePosesMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "component_file", xs_string,  "path to the component file", "NONE" );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ConcatenatePosesMover::fresh_instance() const
{
	return utility::pointer::make_shared< ConcatenatePosesMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ConcatenatePosesMover::clone() const
{
	return utility::pointer::make_shared< ConcatenatePosesMover >( *this );
}

std::string ConcatenatePosesMover::get_name() const {
	return mover_name();
}

std::string ConcatenatePosesMover::mover_name() {
	return "ConcatenatePosesMover";
}

std::string ConcatenatePosesMover::get_component_file() const {
	return component_file_;
}
void ConcatenatePosesMover::set_component_file( std::string new_component_file ) {
	component_file_ = new_component_file;
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
ConcatenatePosesMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< ConcatenatePosesMover >();
}

std::string
ConcatenatePosesMoverCreator::keyname() const
{
	return ConcatenatePosesMover::mover_name();
}

void ConcatenatePosesMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConcatenatePosesMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ConcatenatePosesMover const & mover )
{
	mover.show(os);
	return os;
}

} //simple_moves
} //protocols
