// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/RotamerRecoveryMover.cc
/// @brief A wrapper that measures how similar the rotamers are between before and after running the child mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()

#include <protocols/moves/RotamerRecoveryMover.hh>
#include <string>


// Setup Mover
#include <protocols/moves/RotamerRecoveryMoverCreator.hh>
namespace protocols{
namespace moves{

std::string
RotamerRecoveryMoverCreator::keyname() const
{
	return RotamerRecoveryMoverCreator::mover_name();
}

MoverOP
RotamerRecoveryMoverCreator::create_mover() const {
	return new RotamerRecoveryMover;
}

std::string
RotamerRecoveryMoverCreator::mover_name()
{
	return "RotamerRecoveryMover";
}

}
}


// Unit Headers
// Auto-header: duplicate removed #include <protocols/moves/RotamerRecoveryMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/types.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>

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
#include <numeric/angle.functions.hh>

// C++ Headers
#include <algorithm>
#include <iostream>
#include <fstream>
// Auto-header: duplicate removed #include <string>

//using std::ios::app;
using std::endl;
using std::max;
using std::ofstream;
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
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::pack::dunbrack::RotVector;
using core::pack::dunbrack::rotamer_from_chi;
using core::scoring::TenANeighborGraph;
using basic::Tracer;
using utility::vector1;
using protocols::rotamer_recovery::RotamerRecovery;
using protocols::rotamer_recovery::RRComparerRotBins;
using protocols::rotamer_recovery::RRReporterSimple;




namespace protocols {
namespace moves {

static Tracer TR("protocol.moves.RotamerRecoveryMover");

RotamerRecoveryMover::RotamerRecoveryMover() :
	rotamer_recovery_( NULL ),
	scfxn_(NULL),
	task_factory_(new TaskFactory),
	output_fname_("rotamer_recovery.out")
{
	task_factory_->push_back( new InitializeFromCommandline );
	task_factory_->push_back( new RestrictToRepacking );
}

RotamerRecoveryMover::RotamerRecoveryMover(
	string const & reporter,
	string const & output_fname,
	string const & comparer,
	ScoreFunctionOP scfxn,
	TaskFactoryOP task_factory) :
	rotamer_recovery_( new RotamerRecovery( reporter, output_fname, comparer ) ),
	scfxn_( scfxn ),
	task_factory_( task_factory ),
	output_fname_( output_fname )
	{ }

RotamerRecoveryMover::~RotamerRecoveryMover(){}

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
RotamerRecoveryMover::apply( Pose & pose )
{

	assert( rotamer_recovery_ );

	ScoreFunctionOP scfxn( get_scorefunction());
	scfxn->setup_for_scoring(pose);


	PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ));


	rotamer_recovery_->rtmin_rotamer_recovery(
		pose,
		*scfxn,
		*packer_task);

}

string
RotamerRecoveryMover::get_name() const {
	return "RotamerRecoveryMover";
}

MoverOP
RotamerRecoveryMover::fresh_instance() const {
	return new RotamerRecoveryMover;
}


MoverOP
RotamerRecoveryMover::clone() const {
	return new RotamerRecoveryMover( *this );
}

void
RotamerRecoveryMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	DataMap & datamap,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{
	string const scorefxn_key( tag->getOption<std::string>("scorefxn", "score12" ));
	if ( datamap.has( "scorefxns", scorefxn_key ) ) {
		set_scorefunction( datamap.get< ScoreFunction * >("scorefxns", scorefxn_key) );
	} else {
		utility_exit_with_message("ScoreFunction " + scorefxn_key + " not found in OADataMap.");
	}


	if( rotamer_recovery_ ){
		TR << "WARNING: Attempting to redefine rotamer_recovery_ object from Parser Script" << endl;
	}

	string comparer;
	if( tag->hasOption("comparer") ) {
		comparer = tag->getOption<string>("comparer");
	}

	string reporter;
	if( tag->hasOption("reporter") ){
		reporter = tag->getOption<string>("reporter");
	}

	// store this for writing out at the end if requested
	if( tag->hasOption("output_fname") ){
		output_fname_ = tag->getOption<string>("output_fname");
	}

	rotamer_recovery_ = new RotamerRecovery( reporter, output_fname_, comparer );
}

ScoreFunctionOP
RotamerRecoveryMover::get_scorefunction(){
	if ( !scfxn_ )
		scfxn_ = getScoreFunction();

	return scfxn_;
}

void
RotamerRecoveryMover::set_scorefunction(
	ScoreFunctionOP scorefunction
) {
	scfxn_ = scorefunction;
}

void
RotamerRecoveryMover::show() {
	rotamer_recovery_->show();
}

void
RotamerRecoveryMover::show(
	ostream & out
) {
	rotamer_recovery_->show( out );
}


void
RotamerRecoveryMover::write_to_file(){
	write_to_file( output_fname_ );
}

void
RotamerRecoveryMover::write_to_file(
	string const & output_fname
) {
	ofstream fout;
	fout.open( output_fname.c_str(), std::ios::out );
	if( !fout.is_open() ){
		TR << "Unable to open output file '" << output_fname << "'." << endl;
		utility_exit();
	}

	rotamer_recovery_->show( fout );

	fout.close();

}



} // namespace moves
} // namespace protocols
