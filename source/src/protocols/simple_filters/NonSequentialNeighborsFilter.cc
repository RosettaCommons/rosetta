// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NonSequentialNeighborsFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/NonSequentialNeighborsFilter.hh>
#include <protocols/simple_filters/NonSequentialNeighborsFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <core/conformation/Conformation.hh>
#include <numeric/xyzVector.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.NonSequentialNeighborsFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP NonSequentialNeighborsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NonSequentialNeighborsFilter ); }

// XRW TEMP std::string
// XRW TEMP NonSequentialNeighborsFilterCreator::keyname() const { return "NonSequentialNeighbors"; }

//default ctor
NonSequentialNeighborsFilter::NonSequentialNeighborsFilter() :
	protocols::filters::Filter( "NonSequentialNeighbors" ),
	distance_threshold_( 8.0 ),
	neighbor_cutoff_( 10 ),
	bound_( false ),
	resnum_( 0 ),
	jump_( 1 )
{
}

NonSequentialNeighborsFilter::~NonSequentialNeighborsFilter() = default;

void
NonSequentialNeighborsFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	distance_threshold( tag->getOption< core::Real >( "distance_threshold", 8.0 ) );
	neighbor_cutoff( tag->getOption< core::Size >( "neighbor_cutoff", 10 ));
	bound( tag->getOption< bool >( "bound", false ) );
	resnum( tag->getOption< core::Size >( "resnum", 0 ) );
	jump( tag->getOption< core::Size >( "jump", 1 ) );

	TR<<"jump: "<<jump()<<" distance_threshold: "<<distance_threshold()<<" neighbor_cutoff: "<<neighbor_cutoff()<<" bound: "<<bound()<<" resnum: "<<resnum()<<std::endl;
}

bool
NonSequentialNeighborsFilter::apply( core::pose::Pose const & pose ) const {
	compute( pose );
	return( true );
}

void
NonSequentialNeighborsFilter::report( std::ostream &, core::pose::Pose const & pose ) const {
	compute( pose );
}

core::Real
NonSequentialNeighborsFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Size
NonSequentialNeighborsFilter::residue_neighbors( core::pose::Pose const & pose, core::Size const resi ) const{
	core::Size count_neighbors = 0;
	core::Size const target_chain( pose.chain( resi ) );
	for ( core::Size resj = pose.conformation().chain_begin( target_chain ); resj <= pose.conformation().chain_end( target_chain ); ++resj ) {
		if ( resj >= resi - neighbor_cutoff() && resj <= resi + neighbor_cutoff() ) {
			continue;
		}
		core::Real const distance( pose.residue( resi ).xyz( pose.residue( resi ).nbr_atom() ).distance( pose.residue( resj ).xyz( pose.residue( resj ).nbr_atom() ) ) );
		if ( distance <= distance_threshold() ) {
			count_neighbors++;
		}
	}
	return( count_neighbors );
}

core::Real
NonSequentialNeighborsFilter::compute(
	core::pose::Pose const & pose
) const {
	core::pose::Pose copy_pose( pose );

	if ( !bound() ) {
		protocols::rigid::RigidBodyTransMover rbtm( copy_pose, jump() );
		rbtm.step_size( 10000.0 );
		rbtm.apply( copy_pose );
		TR.Debug<<"Unbound complex"<<std::endl;
	}
	if ( resnum() == 0 ) { // working on entire protein
		for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
			core::Size const count_neighbors( residue_neighbors( copy_pose, resi ) );
			TR.Debug<<"neighbors of residue "<<resi<<": "<<count_neighbors<<std::endl;
		}// for resi
		return 1.0; // dummy return
	}// fi resnum==0
	core::Size const count_neighbors( residue_neighbors( copy_pose, resnum() ) );
	TR.Debug<<"neighbors of residue "<<resnum()<<": "<<count_neighbors<<std::endl;
	return( count_neighbors );
}

std::string NonSequentialNeighborsFilter::name() const {
	return class_name();
}

std::string NonSequentialNeighborsFilter::class_name() {
	return "NonSequentialNeighbors";
}

void NonSequentialNeighborsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("distance_threshold", xsct_real, "distance from given residue from which to count the number of neighbors", "8.0")
		+ XMLSchemaAttribute::attribute_w_default("neighbor_cutoff", xsct_non_negative_integer, "number of neighboring residues that defines the cutoff", "10")
		+ XMLSchemaAttribute::attribute_w_default("bound", xsct_rosetta_bool, "is the protein bound or not?", "false")
		+ XMLSchemaAttribute::attribute_w_default("resnum", xsct_non_negative_integer, "The number indicating the residue of interest. If set to 0, the entire protein will be anlayzed.", "0")
		+ XMLSchemaAttribute::attribute_w_default("jump", xsct_non_negative_integer, "jump number", "1");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filters based upon the number of neighboring residues in a non-sequential fashion.", attlist );
}

std::string NonSequentialNeighborsFilterCreator::keyname() const {
	return NonSequentialNeighborsFilter::class_name();
}

protocols::filters::FilterOP
NonSequentialNeighborsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new NonSequentialNeighborsFilter );
}

void NonSequentialNeighborsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NonSequentialNeighborsFilter::provide_xml_schema( xsd );
}


}
}
