// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarloUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.monte_carlo.util" );

using ObjexxFCL::lead_zero_string_of;
using namespace core;
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
		bool superimpose_over_all_instantiated_ = superimpose_over_all_instantiated || check_all_residues_sampled( pose );
		rms = superimpose_with_stepwise_aligner( pose, *native_pose, superimpose_over_all_instantiated_ );

		if ( do_rms_fill_calculation ){
			TR << TR.Blue << "Generating filled-in model for rms_fill... " << TR.Reset << std::endl;
			if ( full_model_pose == 0 ) full_model_pose = build_full_model( pose );
			rms_fill = superimpose_with_stepwise_aligner( *full_model_pose, *native_pose, superimpose_over_all_instantiated_ );
		}
	}

	SilentStructOP s( new BinarySilentStruct( pose, out_tag ) );
	s->add_string_value( "missing", ObjexxFCL::string_of( get_number_missing_residues_and_connections( pose ) ) );

	if ( native_pose != 0 ) {
		s->add_energy( "rms",      rms );
		if ( do_rms_fill_calculation ) s->add_energy( "rms_fill", rms_fill );
	}

	return s;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_to_silent_file( std::string const & out_tag,
											 std::string const & silent_file,
											 pose::Pose const & pose ){
	Pose pose_copy = pose;
	output_to_silent_file( out_tag, silent_file, pose_copy, 0 /*no native -- no superimpose*/ );
}

////////////////////////////////////////////////////////////////////////////////////////////
void
output_to_silent_file( std::string const & silent_file,
											 utility::vector1< pose::PoseOP > & pose_list,
											 pose::PoseCOP native_pose ) {
	for ( Size n = 1; n <= pose_list.size(); n++ ){
		output_to_silent_file( tag_from_pose( *pose_list[n] ), silent_file,
													 *pose_list[n], native_pose );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// given pose, build in other residues in a deterministic manner. No optimization.
// used to estimate 'rms_fill' in stepwise_monte_carlo.
void
build_full_model( pose::Pose const & start_pose, Pose & full_model_pose ){

	// pretty inelegant -- using stepwise monte carlo object, since it has
	// all the functionality we need. though we're not really doing monte carlo.
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	StepWiseMonteCarlo stepwise_monte_carlo( scorefxn );
	options::StepWiseMonteCarloOptionsOP options( new options::StepWiseMonteCarloOptions );
	options->set_skip_deletions( true );
	stepwise_monte_carlo.set_options( options );
	stepwise_monte_carlo.build_full_model( start_pose, full_model_pose );
}

////////////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
build_full_model( pose::Pose const & start_pose ){
	PoseOP full_model_pose( new Pose );
	build_full_model( start_pose, *full_model_pose );
	return full_model_pose;
}

////////////////////////////////////////////////////////////////////////////////////////////
void
filter_out_bulge_res(  utility::vector1< Size > & sample_res,
 											 utility::vector1< Size > const & bulge_res ) {
	utility::vector1< Size > sample_res_new;
	for ( Size n = 1; n <= sample_res.size(); n++ ) {
		if ( !bulge_res.has_value( sample_res[ n ] ) ) sample_res_new.push_back( sample_res[ n ] );
	}
	sample_res = sample_res_new;
}

} //monte_carlo
} //stepwise
} //protocols
