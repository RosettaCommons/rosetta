// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/RMSDMetric.cc
/// @brief A metric to calculate the RMSD between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/RMSDMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueRMSDMetric.hh>
#include <core/scoring/rms_util.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>
#include <core/id/AtomID.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.RMSDMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::pose;
using namespace core::simple_metrics;
using namespace core::scoring;

/////////////////////
/// Constructors  ///
/////////////////////




/// @brief Default constructor
RMSDMetric::RMSDMetric():
	core::simple_metrics::RealMetric()
{
	setup_name_mapping();
}

RMSDMetric::RMSDMetric( PoseCOP ref_pose):
	core::simple_metrics::RealMetric()
{
	set_comparison_pose( ref_pose );
	setup_name_mapping();
}

RMSDMetric::RMSDMetric( PoseCOP ref_pose, ResidueSelectorCOP selector ):
	core::simple_metrics::RealMetric()
{
	set_comparison_pose( ref_pose );
	set_residue_selector( selector );
	setup_name_mapping();
}

void
RMSDMetric::setup_name_mapping(){

	name_mapping_ = get_rmsd_type_name_map();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RMSDMetric::~RMSDMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
RMSDMetric::RMSDMetric( RMSDMetric const & src ):
	RealMetric( src ),
	rmsd_map_(src.rmsd_map_),
	rmsd_type_(src.rmsd_type_),
	override_atom_names_( src.override_atom_names_ ),
	robust_( src.robust_),
	name_mapping_( src.name_mapping_ )
{
	residue_selector_ = src.residue_selector_;
	residue_selector_ref_ = src.residue_selector_ref_;
	ref_pose_ = src.ref_pose_;
}


std::string
RMSDMetric::name() const {
	return name_static();
}

std::string
RMSDMetric::name_static() {
	return "RMSDMetric";

}
std::string
RMSDMetric::metric() const {
	return "rmsd";
}

void
RMSDMetric::set_rmsd_type( rmsd_atoms rmsd_type){
	rmsd_type_ = rmsd_type;
}

void
RMSDMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ = residue_selector;
}

void
RMSDMetric::set_residue_selector_reference(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ref_ = residue_selector;
}

void
RMSDMetric::set_residue_mapping(std::map<core::Size, core::Size> const & rmsd_map ){
	rmsd_map_ = rmsd_map;
}

void
RMSDMetric::set_comparison_pose( PoseCOP pose){
	ref_pose_ =  pose;
}

void
RMSDMetric::set_corresponding_atoms_robust(bool robust){
	robust_ = robust;
}

core::simple_metrics::SimpleMetricOP
RMSDMetric::clone() const {
	return SimpleMetricOP(new RMSDMetric( *this ) );

}

void
RMSDMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(parse_residue_selector( tag, datamap ));
	}

	if ( tag->hasOption("residue_selector_ref") ) {
		set_residue_selector_reference(parse_residue_selector( tag, datamap, "residue_selector_ref" ));
	}


	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption<bool>("use_native", false) && datamap.has_resource("native_pose") ) {
		ref_pose_ = saved_native_pose(datamap);
	} else {
		std::string msg = "A reference pose must be set. Please use the SavePoseMover (embed the RMSDMetric in RunSimpleMetrics ) or pass the native as in:file:native and set use_native to true.";
		utility_exit_with_message(msg);
	}

	if ( tag->hasOption("rmsd_type") ) {
		set_rmsd_type( name_mapping_[ tag->getOption<std::string>("rmsd_type")]);
	}

	set_corresponding_atoms_robust(tag->getOption< bool >("robust", robust_));
}



void
RMSDMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attributes_for_saved_reference_pose( attlist );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Calculate the RMSD for these residues for both reference and main pose." );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector_ref",
		"Selector for the reference pose (input native or passed reference pose. ).  Residues selected must be same number of residues selected for the main selector." );


	attlist + XMLSchemaAttribute::attribute_w_default(
		"robust", xsct_rosetta_bool, "Set whether we are robust to atom mismatches for selected residues."
		"  By default we only match atoms that are corresponding. (True).", "true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	utility::vector1< std::string > rmsd_type_names = get_rmsd_type_names();
	utility::tag::add_schema_restrictions_for_strings( xsd, "rmsd_types", rmsd_type_names);

	attlist + XMLSchemaAttribute("rmsd_type", "rmsd_types", "Type of calculation.  Current choices are: \n" + utility::to_string(rmsd_type_names) );

	std::string description = "\tThis is the RMSD between the input and the set comparison pose.\n"
		"  If native is set on the cmd-line, we will use that.\n"
		"  Default is to calculate all_heavy atoms - but this can be set\n."
		"\n"
		"  Make sure that reference pose and set pose are the same length ."
		"   We match all corresponding atoms for each residue to match."
		"   By default we do not fail and are robust - only matching what we can for each residue.";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

core::Real
RMSDMetric::calculate(const core::pose::Pose & pose) const {

	using namespace utility;
	if ( ! ref_pose_ ) {
		utility_exit_with_message( "Must pass in a reference pose for RMSDMetric.  See RS XSD or use the set_comparison_pose function");
	}

	per_residue_metrics::PerResidueRMSDMetric core_metric = per_residue_metrics::PerResidueRMSDMetric();
	core_metric.set_rmsd_type(rmsd_type_);
	core_metric.set_comparison_pose(ref_pose_);
	core_metric.set_residue_mapping(rmsd_map_);
	core_metric.set_corresponding_atoms_robust(robust_);
	core_metric.set_residue_selector_reference( residue_selector_ref_);
	core_metric.set_residue_selector( residue_selector_);

	std::map< id::AtomID, id::AtomID > atom_map = core_metric.create_atom_id_map( pose );
	return scoring::rms_at_corresponding_atoms_no_super( pose, *ref_pose_, atom_map);
}



void
RMSDMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RMSDMetric::provide_xml_schema( xsd );
}

std::string
RMSDMetricCreator::keyname() const {
	return RMSDMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
RMSDMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new RMSDMetric );

}

} //core
} //simple_metrics
} //metrics






