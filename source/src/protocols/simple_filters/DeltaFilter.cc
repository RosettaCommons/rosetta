// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/simple_filters/DeltaFilter.hh>
#include <protocols/simple_filters/DeltaFilterCreator.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.DeltaFilter" );

/// @brief default ctor
DeltaFilter::DeltaFilter() :
	parent( "Delta" ),
	filter_( /* NULL */ ),
	baseline_( 0.0 ),
	range_( 0.0 ),
	lower_( false ),
	upper_( true ),
	unbound_( false ),
	relax_unbound_( false ),
	changing_baseline_( false ),
	jump_( 0 ),
	reference_pose_(/* NULL */),
	ref_baseline_( 1234567890.0 ), // unlikely "uninitialized" sentinel value
	scorefxn_( /* NULL */ )
{}

bool
DeltaFilter::unbound() const{
	return unbound_;
}

void
DeltaFilter::unbound( bool const u ){
	unbound_ = u;
}

bool
DeltaFilter::relax_unbound() const{
	return relax_unbound_;
}

void
DeltaFilter::relax_unbound( bool const u ){
	relax_unbound_ = u;
}

bool
DeltaFilter::changing_baseline() const{
	return changing_baseline_;
}

void
DeltaFilter::changing_baseline( bool const c ){
	changing_baseline_ = c;
}

core::Size
DeltaFilter::jump() const{
	return jump_;
}

void
DeltaFilter::jump( core::Size const j ){
	jump_ = j;
}

core::Real
DeltaFilter::baseline() const{
	if ( reference_pose_ ) {
		//Hack to avoid keep re-applying (potentially computationally expensive) filter to reference pose
		//This should probably be replaced by some pose observer magic instead
		if ( (ref_baseline_ == 1234567890.0) || (changing_baseline_) ) { // unlikely "uninitialized" sentinel value
			core::pose::Pose p( *reference_pose_ );
			if ( p.size() == 0 ) { // If reference pose wasn't properly initialized, fast fail with interpretable error message
				utility_exit_with_message("Reference pose used with DeltaFilter wasn't initialized properly!");
			}
			relax_mover()->apply( p );
			unbind( p );
			ref_baseline_ = filter_->report_sm( p );
			TR << "Reference pose baseline is " << ref_baseline_ << "." << std::endl;
		}
		return ref_baseline_;
	}

	return baseline_;
}

void
DeltaFilter::ref_baseline( core::Real const rb ){
	ref_baseline_ = rb;
}

void
DeltaFilter::baseline( core::Real const b ){
	baseline_ = b;
}

core::Real
DeltaFilter::range() const{
	return range_;
}

void
DeltaFilter::range( core::Real const r ){
	range_ = r;
}

bool
DeltaFilter::lower() const{
	return lower_;
}

void
DeltaFilter::lower( bool const l ){
	lower_ = l;
}

void
DeltaFilter::upper( bool const u ){
	upper_ = u;
}

bool
DeltaFilter::upper() const{
	return( upper_ );
}

void
DeltaFilter::filter( protocols::filters::FilterOP filter ){
	filter_ = filter;
}

protocols::filters::FilterOP
DeltaFilter::filter() const{
	return filter_;
}

bool
DeltaFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const delta( compute( pose ) );
	if ( upper() && lower() ) {
		return( delta <= range() && delta >= range() );
	}
	if ( ( upper() && delta <= range() ) || ( lower() && delta >= range() ) ) return true;
	return( false );
}

void
DeltaFilter::unbind( core::pose::Pose & pose ) const{
	if ( !unbound() ) return;
	protocols::rigid::RigidBodyTransMover rbtm( pose, jump() );
	rbtm.step_size( 10000.0 );
	rbtm.apply( pose );
	if ( relax_unbound() ) relax_mover()->apply( pose );
}

core::Real
DeltaFilter::compute( core::pose::Pose const & p ) const{
	core::pose::Pose pose( p );
	unbind( pose );
	core::Real const filter_val( filter()->report_sm( pose ) );
	TR<<"Filter "<<filter()->get_user_defined_name()<<" returns "<<filter_val<<". Baseline is "<<baseline()<<std::endl;
	return( filter_val - baseline() );
}

core::Real
DeltaFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::pose::Pose p( pose );
	unbind( p );
	return( compute( p ) );
}

void
DeltaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"DeltaFilter returns "<<compute( pose )<<std::endl;
}

void
DeltaFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	TR << "DeltaFilter"<<std::endl;
	range( tag->getOption< core::Real >( "range", 0.0 ) );
	lower( tag->getOption< bool >( "lower", false ) );
	upper( tag->getOption< bool >( "upper", true ) );
	runtime_assert( lower() || upper() );
	filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	relax_mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover", "null" ), movers ) );
	unbound( tag->getOption< bool >( "unbound", false ) );
	relax_unbound( tag->getOption< bool >( "relax_unbound", false ) );
	changing_baseline( tag->getOption< bool >( "changing_baseline", false ) );
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	if ( unbound() ) {
		jump( tag->getOption< core::Size >( "jump", 1 ) );
	}
	// need to score the pose before packing...
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data );
		TR << "baseline will be caculated once, when first needed..." << std::endl;
	} else if ( tag->hasOption("reference_pdb") ) {
		std::string reference_pdb_filename( tag->getOption< std::string >( "reference_pdb", "" ) );
		reference_pose_ = core::import_pose::pose_from_file( reference_pdb_filename , core::import_pose::PDB_file);
		TR << "baseline will be caculated once, when first needed..." << std::endl;
	} else if ( tag->hasOption( "relax_mover" ) ) {
		core::pose::Pose p( pose );
		(*scorefxn())(p);
		relax_mover()->apply( p );
		unbind( p );
		baseline( filter()->report_sm( p ) );
	}

	TR<<"with options baseline: ";
	if ( reference_pose_ ) TR << "(deferred)";
	else TR << baseline(); // only called if reference_pose_ is NULL
	TR <<" range: "<<range()<<" upper: "<<upper()<<" unbound: "<<unbound()<<" jump: "<<jump()<<" and lower: "<<lower()<<std::endl;
}

protocols::filters::FilterOP
DeltaFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new DeltaFilter() );
}

DeltaFilter::~DeltaFilter()= default;

protocols::filters::FilterOP
DeltaFilter::clone() const{
	return protocols::filters::FilterOP( new DeltaFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP DeltaFilterCreator::create_filter() const { return protocols::filters::FilterOP( new DeltaFilter ); }

// XRW TEMP std::string
// XRW TEMP DeltaFilterCreator::keyname() const { return "Delta"; }

protocols::moves::MoverOP
DeltaFilter::relax_mover() const{
	return relax_mover_;
}

void
DeltaFilter::relax_mover( protocols::moves::MoverOP const m ){
	relax_mover_ = m;
}

void
DeltaFilter::scorefxn( core::scoring::ScoreFunctionOP s ){ scorefxn_ = s; }

core::scoring::ScoreFunctionOP
DeltaFilter::scorefxn() const{ return scorefxn_; }

std::string DeltaFilter::name() const {
	return class_name();
}

std::string DeltaFilter::class_name() {
	return "Delta";
}

void DeltaFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("range", xsct_real, "how much above/below the baseline to allow?", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("lower", xsct_rosetta_bool, "the threshold is upper/lower? Use both if the threshold is an exact value.", "false")
		+ XMLSchemaAttribute::attribute_w_default("upper", xsct_rosetta_bool, "the threshold is upper/lower? Use both if the threshold is an exact value.", "true")
		+ XMLSchemaAttribute::required_attribute("filter", xs_string, "the name of a predefined filter for evaluation.")
		+ XMLSchemaAttribute::attribute_w_default("relax_mover", xs_string, "called at parse-time before setting the baseline.", "null")
		+ XMLSchemaAttribute::attribute_w_default("unbound", xsct_rosetta_bool, "translates the partners by 10000A before evaluating the baseline and the filters. Allows evaluation of the unbound pose.", "false")
		+ XMLSchemaAttribute::attribute_w_default("relax_unbound", xsct_rosetta_bool, "relax the unbound state w/ relax mover?", "false")
		+ XMLSchemaAttribute::attribute_w_default("changing_baseline", xsct_rosetta_bool, "reset baseline value to current value after every accept", "false");

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("jump", xsct_non_negative_integer, "if unbound is set, this can be used to set the jump along which to translate.", "1");

	protocols::rosetta_scripts::attributes_for_saved_reference_pose( attlist, "reference_name");

	attlist + XMLSchemaAttribute::attribute_w_default("reference_pdb", xs_string, "use reference pose from disk", "XRW TO DO");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the difference in a filter's value compared to the input structure", attlist );
}

std::string DeltaFilterCreator::keyname() const {
	return DeltaFilter::class_name();
}

protocols::filters::FilterOP
DeltaFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new DeltaFilter );
}

void DeltaFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeltaFilter::provide_xml_schema( xsd );
}


} // simple_filters
} // protocols
