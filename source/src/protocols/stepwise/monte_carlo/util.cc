// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarloUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.util" );

using ObjexxFCL::lead_zero_string_of;
using namespace core;
using namespace core::pose;
using namespace protocols::stepwise::modeler;

namespace protocols {
namespace stepwise {
namespace monte_carlo {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_to_silent_file( std::string const & out_tag,
	std::string const & silent_file,
	pose::Pose & pose,
	pose::PoseCOP native_pose,
	bool const superimpose_over_all_instantiated /* = false */,
	bool const do_rms_fill_calculation /* = false */ ){

	using namespace core::io::silent;
	SilentStructOP s = prepare_silent_struct( out_tag, pose, native_pose, superimpose_over_all_instantiated, do_rms_fill_calculation );

	// output silent file.
	static SilentFileData const silent_file_data;
	silent_file_data.write_silent_struct( *s, silent_file, false /*score_only*/ );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::io::silent::SilentStructOP
prepare_silent_struct( std::string const & out_tag,
	pose::Pose & pose,
	pose::PoseCOP native_pose,
	bool const superimpose_over_all_instantiated /* = false */,
	bool const do_rms_fill_calculation /* = false */,
	core::pose::PoseOP full_model_pose /* = 0*/ ){

	using namespace core::io::silent;
	using namespace protocols::stepwise::modeler::align;

	Real rms( 0.0 ), rms_fill( 0.0 );
	if ( native_pose != 0 ) {
		// if built from scratch, make sure to superimpose over everything.
		bool superimpose_over_all_instantiated_ = superimpose_over_all_instantiated || core::pose::full_model_info::check_all_residues_sampled( pose );
		rms = superimpose_with_stepwise_aligner( pose, *native_pose, superimpose_over_all_instantiated_ );

		if ( do_rms_fill_calculation ) {
			TR <<  "Generating filled-in model for rms_fill... " << std::endl;
			if ( full_model_pose == 0 ) full_model_pose = build_full_model( pose );
			rms_fill = superimpose_with_stepwise_aligner( *full_model_pose, *native_pose, superimpose_over_all_instantiated_ );
		}
	}

	SilentStructOP s( new BinarySilentStruct( pose, out_tag ) );
	s->add_string_value( "missing", ObjexxFCL::string_of( core::pose::full_model_info::get_number_missing_residues_and_connections( pose ) ) );

	bool const eval_base_pairs = basic::options::option[ basic::options::OptionKeys::rna::evaluate_base_pairs ]();

	if ( native_pose ) {
		if ( eval_base_pairs ) core::pose::rna::add_number_native_base_pairs( pose, *native_pose, *s );
		s->add_energy( "rms",      rms );
		if ( do_rms_fill_calculation ) s->add_energy( "rms_fill", rms_fill );
	} else {
		// Only use the garden-variety BP function if we can't look at the native pose as a point of reference.
		if ( eval_base_pairs  ) core::pose::rna::add_number_base_pairs( pose, *s );
	}

	return s;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_to_silent_file( std::string const & out_tag,
	std::string const & silent_file,
	pose::Pose const & pose ){
	pose::Pose pose_copy = pose;
	output_to_silent_file( out_tag, silent_file, pose_copy, 0 /*no native -- no superimpose*/ );
}

////////////////////////////////////////////////////////////////////////////////////////////
void
output_to_silent_file( std::string const & silent_file,
	utility::vector1< pose::PoseOP > & pose_list,
	pose::PoseCOP native_pose ) {
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		output_to_silent_file( tag_from_pose( *pose_list[n] ), silent_file,
			*pose_list[n], native_pose );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// given pose, build in other residues in a deterministic manner. No optimization.
// used to estimate 'rms_fill' in stepwise_monte_carlo.
void
build_full_model( pose::Pose const & start_pose, pose::Pose & full_model_pose ){
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	options::StepWiseMonteCarloOptionsOP options( new options::StepWiseMonteCarloOptions );
	options->set_skip_deletions( true );
	options->set_enumerate( true ); // prevent randomness.
	options->set_skip_bulge_frequency( 0.2 ); // to avoid getting stuck
	mover::StepWiseMasterMover master_mover( scorefxn, options );
	master_mover.build_full_model( start_pose, full_model_pose );
}

////////////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
build_full_model( pose::Pose const & start_pose ){
	pose::PoseOP full_model_pose( new pose::Pose );
	build_full_model( start_pose, *full_model_pose );
	return full_model_pose;
}

////////////////////////////////////////////////////////////////////////////////////////////
std::string
get_move_type_string( mover::StepWiseMove const & swa_move ) {
	std::string move_type_string = to_string( swa_move.move_type() );
	std::transform(move_type_string.begin(), move_type_string.end(), move_type_string.begin(), ::tolower); // this is why we love C
	return move_type_string;
}

//////////////////////////////////////////////////////////////////////////////////////
std::string
get_all_res_list( pose::Pose & pose ) {
	std::string out_string;
	out_string = make_tag_with_dashes( core::pose::full_model_info::get_res_list_from_full_model_info_const( pose ) );
	utility::vector1< pose::PoseOP > const & other_pose_list = core::pose::full_model_info::const_full_model_info( pose ).other_pose_list();
	if ( other_pose_list.size() == 0 ) return out_string;
	out_string += " [ other_pose: ";
	for ( Size n = 1; n <= other_pose_list.size(); n++ ) {
		out_string += get_all_res_list( *other_pose_list[n] );
		if ( n < other_pose_list.size() ) out_string += "; ";
	}
	out_string += "]";
	return out_string;
}

} //monte_carlo
} //stepwise
} //protocols
