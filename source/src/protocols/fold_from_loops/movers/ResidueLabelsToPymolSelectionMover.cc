// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/RelaseConstraintFromResidueMoverCreator.hh
/// @brief  Prints a Pymol selection command for each label in a pose
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/ResidueLabelsToPymolSelectionMover.hh>
#include <protocols/fold_from_loops/movers/ResidueLabelsToPymolSelectionMoverCreator.hh>

// Protocol headers


// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.ResidueLabelsToPymolSelectionMover" );

namespace protocols {
namespace fold_from_loops {
namespace movers {

ResidueLabelsToPymolSelectionMover::ResidueLabelsToPymolSelectionMover():
	protocols::moves::Mover( ResidueLabelsToPymolSelectionMover::mover_name() )
{}

ResidueLabelsToPymolSelectionMover::~ResidueLabelsToPymolSelectionMover() = default;

void
ResidueLabelsToPymolSelectionMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	pdb_count( tag->getOption( "pdb_count" , default_pdb_count() ) );
}

protocols::moves::MoverOP
ResidueLabelsToPymolSelectionMover::clone() const
{
	return protocols::moves::MoverOP( new ResidueLabelsToPymolSelectionMover( *this ) );
}

protocols::moves::MoverOP
ResidueLabelsToPymolSelectionMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ResidueLabelsToPymolSelectionMover );
}

void
ResidueLabelsToPymolSelectionMover::apply( core::pose::Pose & pose )
{
	std::map< std::string, utility::vector1< std::string > > data;

	for ( core::Size i=1; i <= pose.size(); ++i ) {
		utility::vector1< std::string > labels = pose.pdb_info()->get_reslabels( i );
		for ( auto label : labels ) {
			if ( data.find( label ) == data.end() ) {
				data.insert( std::pair < std::string, utility::vector1< std::string > >( label, utility::vector1< std::string >()));
			}
			utility::vector1< std::string > resdata;
			utility::vector1< std::string > pdbdata = utility::split_whitespace( pose.pdb_info()->pose2pdb( i ) );
			resdata.push_back( "( c. " );
			resdata.push_back( pdbdata[2] );
			resdata.push_back( " and i. " );
			if ( pdb_count_ ) {
				resdata.push_back( pdbdata[1] );
			} else {
				resdata.push_back( std::to_string( i ) );
			}
			resdata.push_back( " )" );
			data[label].push_back( utility::join( resdata, " ") );
		}
	}
	for ( auto label : data ) {
		TR << "sele " << label.first << ", ";
		for ( core::Size i=1; i<= label.second.size(); ++i ) {
			if ( i != 1 ) {
				TR << " OR ";
			}
			TR << label.second[i];
		}
		TR << std::endl;
	}
}

std::string ResidueLabelsToPymolSelectionMover::get_name() const {
	return mover_name();
}

std::string ResidueLabelsToPymolSelectionMover::mover_name() {
	return "ResidueLabelsToPymolSelectionMover";
}

void ResidueLabelsToPymolSelectionMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "pdb_count", xsct_rosetta_bool, "Provide the PDB count instead of Rosetta's", std::to_string( default_pdb_count() ) );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Print selection commands for labeled residues", attlist );
}

std::string ResidueLabelsToPymolSelectionMoverCreator::keyname() const {
	return ResidueLabelsToPymolSelectionMover::mover_name();
}

protocols::moves::MoverOP
ResidueLabelsToPymolSelectionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ResidueLabelsToPymolSelectionMover );
}

void ResidueLabelsToPymolSelectionMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueLabelsToPymolSelectionMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
