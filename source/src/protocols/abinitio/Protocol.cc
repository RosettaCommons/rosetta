// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>


// Utility Headers
#include <basic/Tracer.hh>

#ifdef USE_TENSORFLOW
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.hh>
#else
#include <basic/tensorflow_manager/util.hh>
#endif

#include <basic/options/option.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols" );
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
		evaluator_ = utility::pointer::make_shared< evaluation::MetaPoseEvaluator >();
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

		SilentFileOptions opts;
		SilentFileData sfd(opts);
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//  ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out(opts);
		pss->fill_struct( pose, get_current_tag() );

		evaluate_pose( pose, get_current_tag(), *pss );

		sfd.write_silent_struct( *pss, silent_file, !option[ OptionKeys::abinitio::debug_structures ]()  /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}
void
protocols::abinitio::Protocol::apply( core::pose::Pose& ) {
	reset_status();
	if ( use_trRosetta_constraints_ ) {
		runtime_assert_string_msg( supports_trRosetta_constraints(), "Error in protocols::abinitio::Protocol::apply(): The " + get_name() + " protocol does not support trRosetta constraints at the present time." );
	}
	// structure_store().clear(); // keep this to avoid memory leaks
	// if (! bInitialized_ ) {
	//   tr.Warning << "ClassicAbinitio::apply is called without previous call to init(). Will self-initialize now... this overrides any"
	//         << " custom settings if they existed... ;" << std::endl; // could call init at
	//   init( pose );
	//  }
}

std::string
protocols::abinitio::Protocol::get_name() const {
	return "Protocol";
}

/// @brief Set whether we're using trRosetta constraints.
void
protocols::abinitio::Protocol::set_use_trRosetta_constraints(
	bool const setting
) {
	std::string const errmsg("Error in protocols::abinitio::Protocol::set_use_trRosetta_constraints():  ");
#ifndef USE_TENSORFLOW
	if ( setting == true ) {
		utility_exit_with_message( errmsg + "Cannot use trRosetta without compiling with Tensorflow support!  " + basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaConstraintGenerator" ) );
	}
#else // USE_TENSORFLOW *is* defined
	if ( setting ) {
		runtime_assert_string_msg( supports_trRosetta_constraints(), errmsg + "The " + get_name() + " protocol does not support trRosetta constraints at this time!" );
	}
	if ( use_trRosetta_constraints_ == false && setting == true ) {
		runtime_assert_string_msg( trRosetta_cst_generator_ == nullptr, errmsg + "The set_use_trRosetta_constraints function was called, but a trRosettaConstraintGenerator is already configured!" );
		trRosetta_cst_generator_ = utility::pointer::make_shared< protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGenerator >();
	} else if ( setting == false ) {
		trRosetta_cst_generator_ = nullptr;
	}
#endif // ifndef USE_TENSORFLOW
	use_trRosetta_constraints_ = setting;
}

/// @brief Report whether this protocol supports trRosetta constraints.  The base
/// class returns false; derived classes can override this and return true.
/*virtual*/
bool
protocols::abinitio::Protocol::supports_trRosetta_constraints() const {
	return false;
}

#ifdef USE_TENSORFLOW
/// @brief Provide an trRosetta constraint generator as input.  Owning pointer is stored
/// directly, not cloned!
/// @details Must set use_trRosetta_constraints to true before calling this!
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
protocols::abinitio::Protocol::set_trRosetta_cst_generator(
	protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGeneratorOP cst_gen_in
) {
	runtime_assert_string_msg( use_trRosetta_constraints_, "Error in protocols::abinitio::Protocol::set_trRosetta_cst_generator(): A trRosettaConstraintGenerator was supplied, but the " + get_name() + " protocol is not configured to use trRosetta constraints!" );
	runtime_assert( cst_gen_in != nullptr );
	trRosetta_cst_generator_ = cst_gen_in;
}
#endif


}
}
