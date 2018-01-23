// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/GlycanInfoMover.cc
/// @brief Simple class for outputting glycan information. Currently, it simply prints the information.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/analysis/GlycanInfoMover.hh>
#include <protocols/analysis/GlycanInfoMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.analysis.GlycanInfoMover" );

namespace protocols {
namespace analysis {


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycanInfoMover::GlycanInfoMover():
	protocols::moves::Mover( GlycanInfoMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanInfoMover::GlycanInfoMover( GlycanInfoMover const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanInfoMover::~GlycanInfoMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
GlycanInfoMover::apply( core::pose::Pose& pose){
	apply_const( pose );
}

void
GlycanInfoMover::apply_const(const core::pose::Pose &pose){
	using namespace core::pose::carbohydrates;

	core::Size protein_branches = 0;
	core::Size carbohydrate_residues = 0;
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		if ( pose.residue( resnum ).is_carbohydrate() ) {
			std::string attachment_points = get_attachment_point_string( pose, resnum);
			core::Size parent_res = pose.glycan_tree_set()->get_parent( resnum );
			bool bp = pose.residue( resnum ).is_branch_point();


			std::cout << "Carbohydrate: "<< resnum  <<" "<< pose.pdb_info()->pose2pdb(resnum) << " Parent: " << parent_res << " BP: "<<bp <<" "<< pose.pdb_info()->pose2pdb(resnum) << " " << " CON: " << utility::pad_right( attachment_points, 10) << " DIS: " << pose.glycan_tree_set()->get_distance_to_start( resnum )
				<< " ShortName: "<< pose.residue( resnum ).carbohydrate_info()->short_name() << std::endl;

			carbohydrate_residues += 1;

		} else if ( pose.residue( resnum ).is_branch_point() ) {
			std::cout << "Branch Point: " << pose.residue( resnum ).name3()<<" "<< resnum <<" " <<pose.pdb_info()->pose2pdb(resnum) << std::endl;
			protein_branches += 1;

		}
	}
	std::cout << "Glycan Residues: " << carbohydrate_residues <<std::endl;
	std::cout << "Protein BPs: " << protein_branches << std::endl;

	std::cout << "TREES" << std::endl;
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		if ( pose.glycan_tree_set()->has_tree(resnum) ) {
			core::Size length = pose.glycan_tree_set()->get_tree(resnum)->get_size();
			core::Size root = pose.glycan_tree_set()->get_tree(resnum)->get_root();
			if ( root != 0 ) {
				std::cout << root << " " <<pose.pdb_info()->pose2pdb(root) <<" "<<"Length: "<< length << std::endl;
			} else {
				std::cout << root <<" "<<"Length: "<< length << std::endl;
			}
		}
	}

}


std::string
GlycanInfoMover::get_attachment_point_string( core::pose::Pose const & pose, core::Size resnum){
	using utility::to_string;
	using namespace core::chemical::carbohydrates;

	CarbohydrateInfoCOP info = pose.residue(resnum).carbohydrate_info();
	std::string outstring = "";
	std::string attach = "_->";

	if ( info->mainchain_glycosidic_bond_acceptor() ) {
		outstring = attach + to_string(info->mainchain_glycosidic_bond_acceptor());
	}

	for ( uint i = 1; i <= info->n_branches(); ++i ) {
		outstring = outstring + "," +attach + to_string( info->branch_point( i ));
	}
	return outstring;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycanInfoMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
GlycanInfoMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}
void GlycanInfoMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A simple Mover that (currently) prints information about the glycan trees in the pose.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanInfoMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanInfoMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanInfoMover::clone() const
{
	return protocols::moves::MoverOP( new GlycanInfoMover( *this ) );
}

std::string GlycanInfoMover::get_name() const {
	return mover_name();
}

std::string GlycanInfoMover::mover_name() {
	return "GlycanInfoMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycanInfoMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GlycanInfoMover );
}

std::string
GlycanInfoMoverCreator::keyname() const
{
	return GlycanInfoMover::mover_name();
}

void GlycanInfoMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanInfoMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, GlycanInfoMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //analysis
