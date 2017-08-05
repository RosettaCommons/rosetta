// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecovery.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()


// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RRProtocolMinPack.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh>

// Project Headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility>
#include <utility/vector1.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <sstream>
#include <ostream>
#include <string>

//Auto Headers
#include <protocols/moves/Mover.fwd.hh>
namespace protocols {
namespace rotamer_recovery {

using std::string;
using std::stringstream;
using std::endl;
using std::ostream;
using core::Size;
using core::Real;
using core::pack::task::PackerTask;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using basic::Tracer;
using protocols::moves::Mover;
using utility::vector1;

static THREAD_LOCAL Tracer TR("protocol.rotamer_recovery.RotamerRecovery");

RotamerRecovery::RotamerRecovery() :
	protocol_( RRProtocolOP( new RRProtocolMinPack )),
	comparer_( RRComparerOP( new RRComparerAutomorphicRMSD )),
	reporter_( RRReporterOP( new RRReporterSimple )),
	ignore_unrecognized_res_(false)
{}

RotamerRecovery::RotamerRecovery(
	RRProtocolOP protocol,
	RRComparerOP comparer,
	RRReporterOP reporter) :
	protocol_(std::move( protocol )),
	comparer_(std::move( comparer )),
	reporter_(std::move( reporter )),
	ignore_unrecognized_res_(false)
{
	reporter_->set_protocol_info( protocol_->get_name(), protocol_->get_parameters() );
	reporter_->set_comparer_info( comparer_->get_name(), comparer_->get_parameters() );
}

RotamerRecovery::~RotamerRecovery() = default;

RotamerRecovery::RotamerRecovery(
	RotamerRecovery const & src
) :
	ReferenceCount(),
	protocol_( src.protocol_ ),
	comparer_( src.comparer_ ),
	reporter_( src.reporter_ )
{}

void
RotamerRecovery::reset_recovery(
) {
	reporter_->reset_recovery();
}

void
RotamerRecovery::register_options() const {
	using basic::options::option;
	using basic::options::OptionKeys::in::ignore_unrecognized_res;

	/// options from the main Option file that are relevant in this context ( and should appear in -help output )
	option.add_relevant( ignore_unrecognized_res );
}


Real
RotamerRecovery::run(
	Pose const & pose,
	ScoreFunction const &  score_function,
	PackerTask const & packer_task
) {
	protocol_->run(comparer_, reporter_, pose, score_function, packer_task);
	return reporter_->recovery_rate();
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

	if ( option[ ignore_unrecognized_res ].user() ) {
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
