// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/AverageDegreeFilter.hh>
#include <protocols/protein_interface_design/filters/AverageDegreeFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.AverageDegreeFilter" );

/// @brief default ctor
AverageDegreeFilter::AverageDegreeFilter() :
	parent( "AverageDegree" ),
	task_factory_( /* NULL */ ),
	threshold_( 0 ),
	distance_threshold_( 10.0 )
{}

core::Real
AverageDegreeFilter::threshold() const{
	return threshold_;
}

void
AverageDegreeFilter::threshold( core::Real const t ){
	threshold_ = t;
}

core::pack::task::TaskFactoryOP
AverageDegreeFilter::task_factory() const
{
	return task_factory_;
}

void
AverageDegreeFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

bool
AverageDegreeFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const average_degree( compute( pose ) );
	return( average_degree >= threshold() );
}

core::Real
AverageDegreeFilter::distance_threshold() const{
	return distance_threshold_;
}

void
AverageDegreeFilter::distance_threshold( core::Real const d ){
	distance_threshold_ = d;
}

core::Real
AverageDegreeFilter::compute( core::pose::Pose const & pose ) const{
	runtime_assert( task_factory() != 0 );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	core::Size count_residues( 0 );
	core::Size count_neighbors( 0 );
	for ( core::Size resi=1; resi<=pose.size(); ++resi ) {
		if ( packer_task->being_packed( resi ) ) {
			core::Size resi_neighbors( 0 );
			++count_residues;
			/// which chain is resi on?
			core::Size chain( 1 );
			for ( ; chain <= pose.conformation().num_chains(); ++chain ) {
				if ( pose.conformation().chain_begin( chain ) <= resi && pose.conformation().chain_end( chain ) >= resi ) break;
			}
			core::Size const chain_begin( pose.conformation().chain_begin( chain ) );
			core::Size const chain_end( pose.conformation().chain_end( chain ) );
			core::conformation::Residue const res_target( pose.conformation().residue( resi ) );
			for ( core::Size j=chain_begin; j<=chain_end; ++j ) {
				core::conformation::Residue const resj( pose.residue( j ) );
				core::Real const distance( resj.xyz( resj.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
				if ( distance <= distance_threshold() ) {
					++count_neighbors;
					++resi_neighbors;
				}
			}
			TR<<"Connectivity of "<<res_target.name3()<<resi<<" is "<<resi_neighbors<<std::endl;
		}
	}
	return( (core::Real) count_neighbors / count_residues );
}

core::Real
AverageDegreeFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
AverageDegreeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"AverageDegreeFilter returns "<<compute( pose )<<std::endl;
}

void
AverageDegreeFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	TR << "AverageDegreeFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	distance_threshold( tag->getOption< core::Real >( "distance_threshold", 8.0 ) );
	TR<<"with options threshold: "<<threshold()<<" and distance_threshold "<<distance_threshold()<<std::endl;
}

protocols::filters::FilterOP
AverageDegreeFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new AverageDegreeFilter() );
}

AverageDegreeFilter::~AverageDegreeFilter(){}

protocols::filters::FilterOP
AverageDegreeFilter::clone() const{
	return protocols::filters::FilterOP( new AverageDegreeFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP AverageDegreeFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AverageDegreeFilter ); }

// XRW TEMP std::string
// XRW TEMP AverageDegreeFilterCreator::keyname() const { return "AverageDegree"; }

std::string AverageDegreeFilter::name() const {
	return class_name();
}

std::string AverageDegreeFilter::class_name() {
	return "AverageDegree";
}

void AverageDegreeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Lower limit below which the filter fails", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_threshold", xsct_real, "Count neighbors for residues closer than this distance", "8.0" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "What is the average degree connectivity of a subset of residues? Found to be useful for discriminating non-interacting designs from natural complexes.", attlist );
}

std::string AverageDegreeFilterCreator::keyname() const {
	return AverageDegreeFilter::class_name();
}

protocols::filters::FilterOP
AverageDegreeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new AverageDegreeFilter );
}

void AverageDegreeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AverageDegreeFilter::provide_xml_schema( xsd );
}


} // filters
} // protein_interface_design
} // protocols
