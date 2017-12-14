// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/FNatFilter.cc
/// @brief filtering on fraction of native contacts
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/FNatFilter.hh>
#include <protocols/protein_interface_design/filters/FNatFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <algorithm>
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

FNatFilter::FNatFilter() :
	protocols::filters::Filter( "FNat" ),
	threshold_( 5.0 ),
	reference_pose_( /* NULL */ )
{}

FNatFilter::FNatFilter(protocols::docking::DockJumps const movable_jumps,
	core::Real const threshold,
	core::pose::PoseOP reference_pose)
: protocols::filters::Filter( "FNat" ),
	threshold_(threshold),
	reference_pose_(std::move(reference_pose)),
	movable_jumps_(movable_jumps)
{}

FNatFilter::~FNatFilter() = default;

protocols::filters::FilterOP
FNatFilter::clone() const {
	return protocols::filters::FilterOP( new FNatFilter( *this ) );
}

static basic::Tracer TR( "protocols.protein_interface_design.filters.FNatFilter" );
core::Real
FNatFilter::compute( core::pose::Pose const & pose ) const
{
	return protocols::docking::calc_Fnat(pose, *reference_pose_, scorefxn_, movable_jumps_);
}

bool
FNatFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const f_nat( compute( pose ));
	TR << "f(nat): " << f_nat;
	if ( f_nat <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else TR<<" failing." << std::endl;
	return( false );
}

void
FNatFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const f_nat( compute( pose ));
	out<<"f(nat): " << f_nat <<'\n';
}

core::Real
FNatFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const f_nat( compute( pose ));
	return( (core::Real) f_nat );
}

void
FNatFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
{
	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	///against the ref pose
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map );
	} else {
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( reference_pose ) );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}
	}

	threshold_ = tag->getOption<core::Real>( "threshold", 5 );

	//TODO: support multiple jumps
	auto jump_num = tag->getOption<core::Size>( "jump", 1);

	//TODO: convert jump_num to movable_jumps_ (vector0?)
	movable_jumps_.push_back(jump_num);

	scorefxn_ = rosetta_scripts::parse_score_function( tag, data_map )->clone();

	TR<<"Built FNatFilter with threshold " << threshold_ << ", scorefxn " <<
		rosetta_scripts::get_score_function_name(tag) <<", jump "<< jump_num << std::endl;

}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP FNatFilterCreator::create_filter() const { return protocols::filters::FilterOP( new FNatFilter ); }

// XRW TEMP std::string
// XRW TEMP FNatFilterCreator::keyname() const { return "FNat"; }

std::string FNatFilter::name() const {
	return class_name();
}

std::string FNatFilter::class_name() {
	return "FNat";
}

void FNatFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_saved_reference_pose( attlist, "reference_name" );
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold for Fnat, a docking metric above which the filter fails",  "5" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_non_negative_integer, "Jump across which the computation is carried out, numbered sequentially from 1", "1" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filter for poor values of Fnat, a docking metric", attlist );
}

std::string FNatFilterCreator::keyname() const {
	return FNatFilter::class_name();
}

protocols::filters::FilterOP
FNatFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new FNatFilter );
}

void FNatFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FNatFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // devel


