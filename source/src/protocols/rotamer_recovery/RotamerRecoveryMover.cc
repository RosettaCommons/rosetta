// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecoveryMover.cc
/// @brief A wrapper that measures how similar the rotamers are between before and after running the child mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()

#include <protocols/rotamer_recovery/RotamerRecoveryMover.hh>
#include <string>


// Setup Mover
#include <protocols/rotamer_recovery/RotamerRecoveryMoverCreator.hh>
namespace protocols {
namespace rotamer_recovery {

// XRW TEMP std::string
// XRW TEMP RotamerRecoveryMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RotamerRecoveryMover::mover_name();
// XRW TEMP }

// XRW TEMP moves::MoverOP
// XRW TEMP RotamerRecoveryMoverCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new RotamerRecoveryMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RotamerRecoveryMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "RotamerRecoveryMover";
// XRW TEMP }

}
}


// Unit Headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRProtocolMover.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Option System Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>

// Numeric Headers

// C++ Headers
#include <algorithm>
#include <iostream>

//Auto Headers
#include <core/kinematics/Jump.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector0.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//using std::ios::app;
using std::endl;
using std::max;
using std::ostream;
using std::string;
using core::Real;
using core::Size;
using core::pack::pack_rotamers;
using core::pose::PoseOP;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::InitializeFromCommandline;
using core::pack::task::operation::RestrictToRepacking;
using core::scoring::get_score_function;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using basic::Tracer;
using utility::vector1;
using protocols::moves::MoverOP;
using protocols::rotamer_recovery::RRProtocolOP;
using protocols::rotamer_recovery::RRProtocolMover;
using protocols::rotamer_recovery::RRComparerOP;
using protocols::rotamer_recovery::RRReporterOP;
using protocols::rotamer_recovery::RotamerRecovery;
using protocols::rotamer_recovery::RotamerRecoveryOP;
using protocols::rotamer_recovery::RotamerRecoveryFactory;
using protocols::rosetta_scripts::parse_mover;
using protocols::rosetta_scripts::parse_score_function;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.rotamer_recovery.RotamerRecoveryMover");

RotamerRecoveryMover::RotamerRecoveryMover() :
	rotamer_recovery_( /* NULL */ ),
	scfxn_(/* NULL */),
	task_factory_(core::pack::task::TaskFactoryOP( new TaskFactory ))
{
	using core::pack::task::operation::TaskOperationCOP;
	task_factory_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	task_factory_->push_back( TaskOperationCOP( new RestrictToRepacking ) );
}

RotamerRecoveryMover::RotamerRecoveryMover(
	RotamerRecoveryOP rotamer_recovery,
	ScoreFunctionOP scfxn,
	TaskFactoryOP task_factory) :
	rotamer_recovery_(std::move(rotamer_recovery)),
	scfxn_(std::move( scfxn )),
	task_factory_(std::move( task_factory ))
{}

RotamerRecoveryMover::~RotamerRecoveryMover()= default;

RotamerRecoveryMover::RotamerRecoveryMover( RotamerRecoveryMover const & src):
	//utility::pointer::ReferenceCount(),
	Mover( src ),
	rotamer_recovery_( src.rotamer_recovery_ ),
	scfxn_( src.scfxn_ ),
	task_factory_( src.task_factory_ )
{}


void
RotamerRecoveryMover::register_options() const {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// Use RotamerRecovery to test new score functions eg. all the corrections
	option.add_relevant( corrections::correct );

	// Use full atom binary silent files for best io-performance
	option.add_relevant( in::file::fullatom );
	option.add_relevant( in::file::silent_struct_type );
	option.add_relevant( in::file::silent );

	// If using an outputter that writes to a database improve
	// io-performace by not writing out structures
	option.add_relevant( out::nooutput );

	rotamer_recovery_->register_options();

}

bool
RotamerRecoveryMover::reinitialize_for_each_job() const {
	return false;
}

bool
RotamerRecoveryMover::reinitialize_for_new_input() const {
	return false;
}

void
RotamerRecoveryMover::apply( Pose & pose
) {
	runtime_assert( rotamer_recovery_ != nullptr );
	ScoreFunctionOP scfxn( score_function());
	scfxn->setup_for_scoring(pose);
	PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ));
	rotamer_recovery_->run(pose, *scfxn, *packer_task);
}

// XRW TEMP string
// XRW TEMP RotamerRecoveryMover::get_name() const {
// XRW TEMP  return "RotamerRecoveryMover";
// XRW TEMP }

MoverOP
RotamerRecoveryMover::fresh_instance() const {
	return MoverOP( new RotamerRecoveryMover );
}


MoverOP
RotamerRecoveryMover::clone() const {
	return MoverOP( new RotamerRecoveryMover( *this ) );
}

void
RotamerRecoveryMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & /*filters*/,
	moves::Movers_map const & movers,
	Pose const & /*pose*/ )
{
	score_function( parse_score_function( tag, datamap ) );

	if ( rotamer_recovery_ ) {
		TR.Warning << "Attempting to redefine rotamer_recovery_ object from Parser Script" << endl;
		throw utility::excn::EXCN_RosettaScriptsOption("");
	}

	if ( tag->hasOption("protocol") && (tag->hasOption("mover") || tag->hasOption("mover_name")) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Please either the 'protocol' field or the 'mover' field but not both.");
	}

	RotamerRecoveryFactory * factory(RotamerRecoveryFactory::get_instance());

	RRProtocolOP protocol;
	if ( tag->hasOption("mover") || tag->hasOption("mover_name") ) {
		MoverOP mover = parse_mover(tag->hasOption("mover") ?
			tag->getOption<string>("mover") : tag->getOption<string>("mover_name"), movers);
		protocol = RRProtocolOP( new RRProtocolMover(mover) );
	} else {
		protocol = factory->get_rotamer_recovery_protocol(tag->getOption<string>("protocol", "RRProtocolMinPack"));
	}
	RRComparerOP comparer(
		factory->get_rotamer_recovery_comparer(
		tag->getOption<string>("comparer", "RRComparerAutomorphicRMSD")));

	RRReporterOP reporter(
		factory->get_rotamer_recovery_reporter(
		tag->getOption<string>("reporter", "RRReporterSimple")));

	rotamer_recovery_ = rotamer_recovery::RotamerRecoveryOP( new RotamerRecovery(protocol, comparer, reporter) );
}

ScoreFunctionOP
RotamerRecoveryMover::score_function(){
	if ( !scfxn_ ) {
		scfxn_ = get_score_function();
	}

	return scfxn_;
}

void
RotamerRecoveryMover::score_function(
	ScoreFunctionOP scorefunction
) {
	scfxn_ = scorefunction;
}

void
RotamerRecoveryMover::show() const{
	rotamer_recovery_->show();
}

void
RotamerRecoveryMover::show(ostream & out) const
{
	rotamer_recovery_->show( out );
}

std::string RotamerRecoveryMover::get_name() const {
	return mover_name();
}

std::string RotamerRecoveryMover::mover_name() {
	return "RotamerRecoveryMover";
}

void RotamerRecoveryMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		//Can have protocol OR mover/mover_name but not both
		+ XMLSchemaAttribute::attribute_w_default( "protocol", xs_string, "Name of protocol to use for rotamer recovery detection", "RRProtocolMinPack" )
		+ XMLSchemaAttribute( "mover", xs_string, "Mover to apply to pose. Interchangeable with mover_name." )
		+ XMLSchemaAttribute( "mover_name", xs_string, "Mover to apply to pose. Interchangeable with mover." )
		+ XMLSchemaAttribute::attribute_w_default( "comparer", xs_string, "Name of RRComparer to use when measuring rotamer recovery", "RRComparerAutomorphicRMSD" )
		+ XMLSchemaAttribute::attribute_w_default( "reporter", xs_string, "Name of RRReporter to use to report rotamer recovery", "RRReporterSimple" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A wrapper mover that measures how similar the rotamers are before and after running the child mover", attlist );
}

std::string RotamerRecoveryMoverCreator::keyname() const {
	return RotamerRecoveryMover::mover_name();
}

protocols::moves::MoverOP
RotamerRecoveryMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotamerRecoveryMover );
}

void RotamerRecoveryMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotamerRecoveryMover::provide_xml_schema( xsd );
}


} // namespace rotamer_recovery
} // namespace protocols
