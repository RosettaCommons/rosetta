// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/BindingStrainFilter.hh>
#include <protocols/protein_interface_design/filters/BindingStrainFilterCreator.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.BindingStrainFilter" );

/// @brief default ctor
BindingStrainFilter::BindingStrainFilter() :
	parent( "BindingStrain" ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	relax_mover_( /* NULL */ ),
	jump_( 1 ),
	threshold_( 3.0 )
{}

core::Size
BindingStrainFilter::jump() const{
	return jump_;
}

void
BindingStrainFilter::jump( core::Size const j ){
	jump_ = j;
}

core::pack::task::TaskFactoryOP
BindingStrainFilter::task_factory() const{
	return task_factory_;
}

void
BindingStrainFilter::task_factory( core::pack::task::TaskFactoryOP task ){
	task_factory_ = task;
}

core::scoring::ScoreFunctionOP
BindingStrainFilter::scorefxn() const{
	return scorefxn_;
}

void
BindingStrainFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::Real
BindingStrainFilter::threshold() const{
	return threshold_;
}

void
BindingStrainFilter::threshold( core::Real const t ){
	threshold_ = t;
}

bool
BindingStrainFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const strain( compute( pose ) );
	return( strain <= threshold() );
}

void
BindingStrainFilter::unbind( core::pose::Pose & pose ) const{
	protocols::rigid::RigidBodyTransMover rbtm( pose, jump() );
	rbtm.step_size( 10000.0 );
	rbtm.apply( pose );
}

core::Real
BindingStrainFilter::compute( core::pose::Pose const & p ) const{
	using namespace core::pack;
	using namespace core::pack::task;

	PackerTaskOP pack = task_factory()->create_task_and_apply_taskoperations( p );
	pack->restrict_to_repacking();
	pack->initialize_from_command_line();
	protocols::simple_filters::ScoreTypeFilter const stf( scorefxn(), core::scoring::total_score, 0.0 );
	core::pose::Pose pose( p );
	unbind( pose );
	core::Real const energy_before_pack( stf.compute( pose ));

	protocols::simple_moves::PackRotamersMoverOP prm;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		prm = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), pack ) );
	} else {
		prm = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( scorefxn(), pack ) );
	}
	prm->apply( pose );
	core::Real const energy_after_pack( stf.compute( pose ) );
	return( energy_before_pack - energy_after_pack );
}

core::Real
BindingStrainFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
BindingStrainFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"BindingStrainFilter returns "<<compute( pose )<<std::endl;
}

void
BindingStrainFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	using namespace protocols::rosetta_scripts;
	TR << "BindingStrainFilter"<<std::endl;
	task_factory( parse_task_operations( tag, data ) );
	scorefxn( parse_score_function( tag, data ) );
	jump( tag->getOption< core::Size >( "jump", 1 ) );
	runtime_assert( jump() > 0 );
	threshold( tag->getOption< core::Real >( "threshold", 3.0 ));
	relax_mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover", "null" ), movers ) );

	TR<<"with options jump: "<<jump()<<" and strain threshold: "<<threshold()<<std::endl;
}

protocols::filters::FilterOP
BindingStrainFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new BindingStrainFilter() );
}

BindingStrainFilter::~BindingStrainFilter(){}

protocols::filters::FilterOP
BindingStrainFilter::clone() const{
	return protocols::filters::FilterOP( new BindingStrainFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP BindingStrainFilterCreator::create_filter() const { return protocols::filters::FilterOP( new BindingStrainFilter ); }

// XRW TEMP std::string
// XRW TEMP BindingStrainFilterCreator::keyname() const { return "BindingStrain"; }

protocols::moves::MoverOP
BindingStrainFilter::relax_mover() const{
	return relax_mover_;
}

void
BindingStrainFilter::relax_mover( protocols::moves::MoverOP const m ){
	relax_mover_ = m;
}

std::string BindingStrainFilter::name() const {
	return class_name();
}

std::string BindingStrainFilter::class_name() {
	return "BindingStrain";
}

void BindingStrainFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "jump", xsct_non_negative_integer, "Jump across which to compute binding", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Strain must be less than this value to pass", "3.0" )
		+ XMLSchemaAttribute::attribute_w_default( "relax_mover", xs_string, "Relax mover employed to relieve the unbound state", "null" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the energetic strain in a bound monomer. Automatically respects symmetry", attlist );
}

std::string BindingStrainFilterCreator::keyname() const {
	return BindingStrainFilter::class_name();
}

protocols::filters::FilterOP
BindingStrainFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new BindingStrainFilter );
}

void BindingStrainFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BindingStrainFilter::provide_xml_schema( xsd );
}


} // filters
} // protein_interface_design
} // protocols
