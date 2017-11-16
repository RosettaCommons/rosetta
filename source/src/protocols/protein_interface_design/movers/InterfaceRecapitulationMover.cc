// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/DesignRepackMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/InterfaceRecapitulationMover.hh>
#include <protocols/protein_interface_design/movers/InterfaceRecapitulationMoverCreator.hh>

// Package headers
#include <protocols/protein_interface_design/ReportPSSMDifference.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/protein_interface_design/design_utils.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// Utility Headers

// Unit Headers

// C++ headers
#include <map>
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace core;

// option key includes

using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.InterfaceRecapitulationMover" );


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

// XRW TEMP std::string
// XRW TEMP InterfaceRecapitulationMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return InterfaceRecapitulationMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP InterfaceRecapitulationMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new InterfaceRecapitulationMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP InterfaceRecapitulationMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "InterfaceRecapitulation";
// XRW TEMP }

InterfaceRecapitulationMover::InterfaceRecapitulationMover() :
	Mover( InterfaceRecapitulationMover::mover_name() ),
	saved_pose_( /* NULL */ ),
	design_mover_( /* NULL */ ),
	design_mover2_( /* NULL */ ),
	pssm_( false )
{}

MoverOP
InterfaceRecapitulationMover::clone() const
{ return MoverOP( new InterfaceRecapitulationMover( *this ) ); }

MoverOP
InterfaceRecapitulationMover::fresh_instance() const
{ return protocols::moves::MoverOP( new InterfaceRecapitulationMover ); }

core::pose::PoseCOP
InterfaceRecapitulationMover::get_reference_pose() const
{
	return saved_pose_;
}

void
InterfaceRecapitulationMover::set_reference_pose( core::pose::PoseOP pose )
{
	saved_pose_ = pose;
}

InterfaceRecapitulationMover::~InterfaceRecapitulationMover() {}

void
InterfaceRecapitulationMover::apply( core::pose::Pose & pose ){
	runtime_assert( saved_pose_ != 0 );
	runtime_assert( design_mover_ || design_mover2_ );
	core::pack::task::PackerTaskCOP task;
	if ( design_mover_ ) {
		design_mover_->apply( pose );
		task = design_mover_->task();
	} else {
		design_mover2_->apply( pose );
		task = design_mover2_->task();
	}
	core::Size designable_positions( 0 );
	for ( Size i( 1 ); i<= pose.size(); ++i ) {
		if ( task->being_designed( i ) ) ++designable_positions;
	}

	ReportPSSMDifferences rsd;
	if ( pssm_ && !rsd.load_pssm_data( protocols::jd2::current_input_tag() ) ) {
		pssm_ = false;
	}

	if ( !pssm_ ) {
		ReportSequenceDifferences rsd( core::scoring::get_score_function() );
		rsd.calculate( *get_reference_pose(), pose );
		std::map< core::Size, std::string > const res_names1( rsd.res_name1() );
		core::Size const mutated( res_names1.size() );
		core::Real const rate( (core::Real) mutated / designable_positions );
		TR<<"Your design mover mutated "<<mutated<<" positions out of "<<designable_positions<<" designable positions. Sequence recovery is: "<<1-rate<<std::endl;
	} else {
		core::Real pssm = rsd.calculate( *get_reference_pose(), pose, task );
		TR << "PSSM-Score: " << pssm << " at " << designable_positions << " designable positions. Mean score is " << pssm / (core::Real)designable_positions << std::endl;
	}

}

// XRW TEMP std::string
// XRW TEMP InterfaceRecapitulationMover::get_name() const {
// XRW TEMP  return InterfaceRecapitulationMover::mover_name();
// XRW TEMP }

void
InterfaceRecapitulationMover::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const & movers, core::pose::Pose const & pose ){
	set_reference_pose( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
	std::string const mover_name( tag->getOption<std::string>( "mover_name" ) );
	std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( mover_name ));
	bool const mover_found( find_mover != movers.end() );
	if ( mover_found ) {
		design_mover_ = utility::pointer::dynamic_pointer_cast< simple_moves::DesignRepackMover > ( find_mover->second );
		if ( !design_mover_ ) {
			design_mover2_ = utility::pointer::dynamic_pointer_cast< protocols::simple_moves::PackRotamersMover > ( find_mover->second );
			if ( !design_mover2_ ) {
				throw utility::excn::EXCN_RosettaScriptsOption( "dynamic cast failed in tag in RecapitulateMover. Make sure that the mover is either PackRotamers or DesignRepackMover derived" );
			}
		}
		pssm_ = tag->getOption<bool>( "pssm", false );
	}
}

std::string InterfaceRecapitulationMover::get_name() const {
	return mover_name();
}

std::string InterfaceRecapitulationMover::mover_name() {
	return "InterfaceRecapitulation";
}

void InterfaceRecapitulationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "mover_name", xs_string, "Name of mover to be applied, either a DesignRepack or PackRotamers" )
		+ XMLSchemaAttribute::attribute_w_default( "pssm", xsct_rosetta_bool, "Should the pssm score be calculated against a possibly-provided pose reference?", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string InterfaceRecapitulationMoverCreator::keyname() const {
	return InterfaceRecapitulationMover::mover_name();
}

protocols::moves::MoverOP
InterfaceRecapitulationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new InterfaceRecapitulationMover );
}

void InterfaceRecapitulationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceRecapitulationMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
