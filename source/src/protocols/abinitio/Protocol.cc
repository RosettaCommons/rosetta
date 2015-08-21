// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @details
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/Protocol.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>


// Utility Headers
#include <basic/Tracer.hh>

#include <basic/options/option.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols" );
using namespace core;

void
protocols::abinitio::Protocol::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::debug_structures );
}

namespace protocols {
namespace abinitio {
// Project Header


protocols::abinitio::Protocol::Protocol() :
	checkpoints_("Protocol"),
	evaluator_( /* NULL */ ),
	return_centroid_( true )
{
	if ( basic::options::option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		silentout_file_name_ = basic::options::option[ basic::options::OptionKeys::out::file::silent ]();
	}
}

void
protocols::abinitio::Protocol::evaluate_pose(
	core::pose::Pose &pose,
	std::string tag,
	core::io::silent::SilentStruct &pss) const
{
	if ( evaluator_ ) {
		evaluator_->apply( pose, tag, pss );
	}
}

void
protocols::abinitio::Protocol::add_evaluation( evaluation::PoseEvaluatorOP ev ) {
	if ( !evaluator_ ) {
		evaluator_ = evaluation::MetaPoseEvaluatorOP( new evaluation::MetaPoseEvaluator );
	}
	evaluator_->add_evaluation( ev );
}

void
protocols::abinitio::Protocol::set_evaluation( evaluation::MetaPoseEvaluatorOP ev ) {
	evaluator_ = ev;
}

void
protocols::abinitio::Protocol::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		std::string silent_file="_"+prefix;
		if ( get_current_job() && get_current_job()->output_file_name() != "" ) {
			silent_file = get_current_job()->output_file_name()+silent_file;
		} else silent_file = silentout_file_name_+silent_file;

		SilentFileData sfd;
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//  ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		pss->fill_struct( pose, get_current_tag() );

		evaluate_pose( pose, get_current_tag(), *pss );

		sfd.write_silent_struct( *pss, silent_file, !option[ OptionKeys::abinitio::debug_structures ]()  /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}
void
protocols::abinitio::Protocol::apply( core::pose::Pose& ) {
	reset_status();
	// structure_store().clear(); // keep this to avoid memory leaks
	// if (! bInitialized_ ) {
	//   tr.Warning << "WARNING: ClassicAbinitio::apply is called without previous call to init(). Will self-initialize now... this overrides any"
	//         << " custom settings if they existed... ;" << std::endl; // could call init at
	//   init( pose );
	//  }
}

std::string
protocols::abinitio::Protocol::get_name() const {
	return "Protocol";
}


}
}
