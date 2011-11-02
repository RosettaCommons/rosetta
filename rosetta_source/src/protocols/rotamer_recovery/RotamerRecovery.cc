// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecovery.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()


// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRReporterHuman.hh>
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
#include <core/pack/rtmin.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <ostream>
#include <string>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {

using std::string;
using std::endl;
using std::ostream;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pack::RTMin;
using core::pack::task::PackerTask;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using basic::Tracer;
using protocols::moves::Mover;
using utility::vector1;

static Tracer TR("protocol.rotamer_recovery.RotamerRecovery");



RotamerRecovery::RotamerRecovery() :
	reporter_( new RRReporterSimple),
	comparer_( new RRComparerRotBins),
	ignore_unrecognized_res_(false)
{}

RotamerRecovery::RotamerRecovery(
	RRReporterOP reporter,
	RRComparerOP comparer ) :
	reporter_( reporter ),
	comparer_( comparer ),
	ignore_unrecognized_res_(false)
{
	reporter_->set_comparer_info( comparer_->get_name(), comparer_->get_parameters() );
}

RotamerRecovery::RotamerRecovery(
	string const & reporter,
	string const & output_fname,
	string const & comparer) {

	if( comparer.compare("RRComparerRotBins") == 0 ){
		comparer_ = new RRComparerRotBins();
	} else if( comparer.compare("RRComparerAutomorphicRMSD") == 0 ){
		comparer_ = new RRComparerAutomorphicRMSD();
	} else {
		TR << "Unrecognized rotamer recovery comparer '" << comparer << "'." << endl;
		TR << "Using default comparer 'RRComparerRotBins'." << endl;
		comparer_ = new RRComparerRotBins();
	}

	if( reporter.compare("RRRepoterSimple") == 0 ){
		reporter_ = new RRReporterSimple();
	} else if ( reporter.compare( "RRReporterHuman" ) == 0 ) {
		reporter_ = new RRReporterHuman();
	}
	else if ( reporter.compare("RRReporterSQLite") == 0 ){
		reporter_ = new RRReporterSQLite( output_fname );
	}
	else {
		TR << "Unrecongized rotamer recovery reporter '" << reporter << "'." << endl;
		TR << "Using default reporter 'RRReporterSimple'." << endl;
		reporter_ = new RRReporterSimple();
	}

	reporter_->set_comparer_info( comparer_->get_name(), comparer_->get_parameters() );

}





RotamerRecovery::~RotamerRecovery() {}

RotamerRecovery::RotamerRecovery(
	RotamerRecovery const & src
) :
	ReferenceCount(),
	reporter_( src.reporter_ ),
	comparer_( src.comparer_ )
{}

void
RotamerRecovery::reset_recovery(
) {
	reporter_->reset_recovery();
	//comparer_->reset_recovery();
}

void
RotamerRecovery::register_options() const {
	using basic::options::option;
	using basic::options::OptionKeys::in::ignore_unrecognized_res;

	/// options from the main Option file that are relevant in this context ( and should appear in -help output )
	option.add_relevant( ignore_unrecognized_res );
}


bool
RotamerRecovery::measure_rotamer_recovery(
	Pose const & pose1,
	Pose const & pose2,
	Residue const & res1,
	Residue const & res2
) {

	Real score;
	bool recovered;
	if(comparer_->measure_rotamer_recovery(
			pose1, pose2, res1, res2, score, recovered)){
		reporter_->report_rotamer_recovery(
			pose1, pose2, res1, res2, score, recovered );
	}
	return recovered;
}

Real
RotamerRecovery::rtmin_rotamer_recovery(
	Pose const & pose,
	ScoreFunction const & score_function,
	PackerTask const & packer_task
) {

	using core::chemical::aa_unk;

	// Assume score_function.setup_for_scoring(pose) has already been called.

	PackerTaskOP one_res_task( packer_task.clone() );

	RTMin rtmin;

	// I don't know if rtmin looks at more than pack_residue(..)
	one_res_task->temporarily_fix_everything();

	// For each residue in the packer task,
	// rtmin residue -> and measure recovery
	for( Size ii = 1; ii <= pose.total_residue(); ++ii ){
		if ( !packer_task.pack_residue(ii) ) continue;
		if ( ignore_unrecognized_res_ && pose.residue(ii).aa() == aa_unk ){
			continue;
		}
		Pose working_pose = pose;  // deep copy
		one_res_task->temporarily_set_pack_residue( ii, true );
		rtmin.rtmin( working_pose, score_function, one_res_task );
		measure_rotamer_recovery(
			pose, working_pose,
			pose.residue(ii), working_pose.residue(ii) );
		one_res_task->temporarily_set_pack_residue( ii, false );
	}
	return reporter_->recovery_rate();
}


void
RotamerRecovery::compare_rotamers(
	Pose const & ref_pose,
	Mover & mover,
	vector1< Size > const & res_ids
) {

	Pose pose = ref_pose;
	mover.apply(pose);

	for( Size resi = 1; resi <= res_ids.size(); ++resi ){
		if (res_ids[resi]){
			measure_rotamer_recovery(
				ref_pose,
				pose,
				ref_pose.residue(resi),
				pose.residue(resi) );
		}
	}
}

void
RotamerRecovery::show( ostream & out ) const {
	reporter_->show( out );
}

void
RotamerRecovery::show(){
	reporter_->show( TR );
}

Real
RotamerRecovery::recovery_rate() const {
	return reporter_->recovery_rate();
}

void
RotamerRecovery::init_rotamer_recovery_with_options(
	RotamerRecovery & rotamer_recovery
) {

	using basic::options::option;
	using basic::options::OptionKeys::in::ignore_unrecognized_res;

	if( option[ ignore_unrecognized_res ].user() ){
		rotamer_recovery.set_ignore_unrecognized_res(
			option[ ignore_unrecognized_res ].value() );
	}
}

void
RotamerRecovery::init_with_options()
{
	init_rotamer_recovery_with_options(*this);
}

void
RotamerRecovery::set_ignore_unrecognized_res(
	bool const ignore_unrecognized_res
) {
	ignore_unrecognized_res_ = ignore_unrecognized_res;
}

bool
RotamerRecovery::get_ignore_unrecognized_res()
{
	return ignore_unrecognized_res_;
}

} // rotamer_recovery
} // protocols
