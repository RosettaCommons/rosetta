// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/WorkUnitBase.cc
/// @brief
/// @author Mike Tyka

#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

namespace protocols {
namespace relax {


static THREAD_LOCAL basic::Tracer TR( "WorkUnit_BatchRelax" );


///  WorkUnit_BatchRelax

WorkUnit_BatchRelax::WorkUnit_BatchRelax(): protocols::wum::WorkUnit_SilentStructStore() {}

WorkUnit_BatchRelax::~WorkUnit_BatchRelax() {}

protocols::wum::WorkUnitBaseOP
WorkUnit_BatchRelax::clone() const {
	return protocols::wum::WorkUnitBaseOP( new WorkUnit_BatchRelax( *this ) );
}


void
WorkUnit_BatchRelax::run(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// scorefxn and sequence file are controllable from commandline - unfortunately that means
	// there can only be one setting. The batch solution handles this better but is also much more
	// restricted in a way.
	if ( !scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn_ );
	}
	protocols::relax::FastRelax relax( scorefxn_,  option[ OptionKeys::relax::sequence_file ]() );
	if ( native_pose_ ) {
		relax.set_native_pose( native_pose_ );
	}

	TR << "Pre Processing: "  << std::endl;
	pre_process();

	// relax handles all the silent structures in situ
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		core::pose::Pose pose;
		core::io::silent::SilentStructOP ss = decoys().get_struct(0);
		ss->fill_pose( pose );
		using namespace core::scoring::constraints;
		ConstraintSetOP cstset = ConstraintIO::get_instance()->read_constraints( get_cst_fa_file_option(), ConstraintSetOP( new ConstraintSet ), pose  );
		relax.batch_apply( decoys().store(), cstset );
	} else {
		relax.batch_apply( decoys().store() );
	}

	TR << "Post Processing: "  << std::endl;

	post_process();

	TR << "NRELAXED: " << decoys().size() << std::endl;
}


void
WorkUnit_BatchRelax::set_native_pose( core::pose::PoseCOP native_pose){
	native_pose_ = native_pose;
}


// @brief Hook for post processing such as rescoring etc.
void
WorkUnit_BatchRelax::pre_process(){
	TR << "WorkUnit_BatchRelax::pre_process()" << std::endl;
}


// @brief Hook for post processing such as rescoring etc.
void
WorkUnit_BatchRelax::post_process(){
	TR << "WorkUnit_BatchRelax::post_process()" << std::endl;
}


// -----------------------------------------------------------------------
//
// WorkUnit_BatchRelax_and_PostRescore


WorkUnit_BatchRelax_and_PostRescore::WorkUnit_BatchRelax_and_PostRescore():
	WorkUnit_BatchRelax()
{
	set_defaults();
}


void
WorkUnit_BatchRelax_and_PostRescore::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	trim_proportion_ = option[ OptionKeys::wum::trim_proportion ]();

	if (  option[ OptionKeys::wum::extra_scorefxn_ref_structure ].user() ) {
		core::import_pose::pose_from_pdb( ref_pose_, option[ OptionKeys::wum::extra_scorefxn_ref_structure ]() );
		// We do the following to save memory, since we dont actually need anything but the backbone coordiantes from this pose.
		//core::util::switch_to_residue_type_set( ref_pose_, core::chemical::CENTROID);
		if ( !ref_pose_.is_fullatom() ) {
			core::util::switch_to_residue_type_set( ref_pose_, core::chemical::FA_STANDARD);
		}
	}
}


WorkUnit_BatchRelax_and_PostRescore::~WorkUnit_BatchRelax_and_PostRescore(){}

protocols::wum::WorkUnitBaseOP
WorkUnit_BatchRelax_and_PostRescore::clone() const
{
	return protocols::wum::WorkUnitBaseOP( new WorkUnit_BatchRelax_and_PostRescore( *this ) );
}


// Here the pre_process routing will take each structure in turn and apply a new score function including, potentially,
// the electron density score. THen it will remove the worst n percent of the structures
void
WorkUnit_BatchRelax_and_PostRescore::pre_process(){
	// if trimming is to occur rescore all the input decoys and trim on combined score
	if ( trim_proportion_ > 0.0 ) {
		rescore_all_decoys();
		trim();
	}
}

// Here the post_process routing will take each structure in turn and apply a new score function including, potentially,
// the electron density score.

void
WorkUnit_BatchRelax_and_PostRescore::post_process(){
	rescore_all_decoys();
}


void
WorkUnit_BatchRelax_and_PostRescore::rescore_all_decoys(){
	using basic::options::option;
	using namespace basic::options;
	using namespace core;
	//bool superimpose_to_ref = true;

	TR << "rescore_all_decoys" << std::endl;

	if ( !basic::options::option[ OptionKeys::wum::extra_scorefxn ].user() ) return;

	core::Size starttime = time(NULL);

	core::scoring::ScoreFunctionOP extra_scorefxn;
	std::string weight_set = option[ OptionKeys::wum::extra_scorefxn];
	TR << "weights_set " << weight_set << std::endl;
	extra_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( weight_set );

	core::scoring::ScoreFunctionOP combined_scorefxn;
	combined_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( weight_set );
	combined_scorefxn->merge( *scorefxn_ );

	for ( protocols::wum::SilentStructStore::iterator struc = decoys().begin(); struc != decoys().end(); ++ struc ) {
		core::pose::Pose pose;
		(*struc)->fill_pose( pose );

		if ( !pose.is_fullatom() ) {
			TR.Debug << "Switching struct to fullatom" << std::endl;
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);
		}


		if ( ref_pose_.total_residue() > 0 ) {
			protocols::simple_moves::SuperimposeMover sm;
			sm.set_reference_pose( ref_pose_ );
			sm.apply( pose );
		}


		if ( option[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::electron_density::SetupForDensityScoringMover setup_for_density;
			setup_for_density.apply( pose );
		}

		core::Real the_extra_score = (*extra_scorefxn)(pose);


		if ( option[ OptionKeys::wum::extra_scorefxn_relax]() > 0 ) {
			// combine the score functions into one


			core::Size pre_relax_time = time(NULL);
			relax::FastRelax final_relax( combined_scorefxn, option[ OptionKeys::wum::extra_scorefxn_relax]() );
			final_relax.apply( pose );

			// make sure the total score reflects the normal score, not the combined score. (that is saved seperately)
			core::Real normal_score = (*scorefxn_)( pose );
			(*struc)->fill_struct( pose );
			core::Size post_relax_time = time(NULL);

			TR << "extra_scorefxn_relax time: " << (post_relax_time - pre_relax_time ) << "  " << normal_score << std::endl;


		}


		(*struc)->add_energy( "extra_score", the_extra_score );
		(*struc)->add_energy( "combined_score", the_extra_score +  (*struc)->get_energy("score") );
		TR << "ExtraScore: " << the_extra_score << "  " << the_extra_score +  (*struc)->get_energy("score") << std::endl;
	}

	core::Size endtime = time(NULL);

	TR << "rescore_all_decoys.end " << (endtime - starttime) << std::endl;

}


void
WorkUnit_BatchRelax_and_PostRescore::trim(){
	decoys().sort_by("combined_score");
	decoys().limit( (core::Size)( core::Real(decoys().size() * (1.0 - trim_proportion_)) ));
}


}
}

