// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/RefineOneCDRLoopCentroid.cc
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/VariantType.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/antibody/RefineOneCDRLoopCentroid.hh>
#include <protocols/antibody/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>//JQX:: this header file took care of the "CUTPOINT_LOWER" options below

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.antibody.RefineOneCDRLoopCentroid" );

using namespace core;
namespace protocols {
namespace antibody {

RefineOneCDRLoopCentroid::RefineOneCDRLoopCentroid( AntibodyInfoCOP antibody_info,
        CDRNameEnum const & loop_name ) : Mover() {
	the_cdr_loop_ =  antibody_info->get_CDR_loop(loop_name);
	set_default();
}

RefineOneCDRLoopCentroid::RefineOneCDRLoopCentroid( AntibodyInfoCOP antibody_info,
        CDRNameEnum const & loop_name,
        core::scoring::ScoreFunctionCOP scorefxn) : Mover() {
	lowres_scorefxn_ = scorefxn->clone();
	the_cdr_loop_ =  antibody_info->get_CDR_loop(loop_name);
	set_default();
}


RefineOneCDRLoopCentroid::RefineOneCDRLoopCentroid( loops::Loop const & a_cdr_loop ) : Mover() {
	the_cdr_loop_ = a_cdr_loop;
	set_default();
}

RefineOneCDRLoopCentroid::RefineOneCDRLoopCentroid( loops::Loop const & a_cdr_loop,
        core::scoring::ScoreFunctionCOP scorefxn  ) : Mover() {
	lowres_scorefxn_ = scorefxn->clone();
	the_cdr_loop_ = a_cdr_loop;
	set_default();
}




void RefineOneCDRLoopCentroid::set_default() {

	benchmark_          = false;
	snug_fit_           = false;
	refine_input_loop_  = true;

	TR << "Finished Setting Defaults" << std::endl;

}



// default destructor
RefineOneCDRLoopCentroid::~RefineOneCDRLoopCentroid() {}





std::string RefineOneCDRLoopCentroid::get_name() const {
	return "RefineOneCDRLoopCentroid";
}

void RefineOneCDRLoopCentroid::set_score_function(core::scoring::ScoreFunctionCOP scorefxn) {
	lowres_scorefxn_ = scorefxn->clone();
}





void RefineOneCDRLoopCentroid::finalize_setup( core::pose::Pose const & /*pose*/ ) {
	TR<<"   start finalize_setup function ..."<<std::endl;
	if (!lowres_scorefxn_) {
		lowres_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
		lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
	}
	TR<<"   finish finalize_setup function !!!"<<std::endl;

}



void RefineOneCDRLoopCentroid::apply( pose::Pose & pose ) {
	TR<<"start apply function ..."<<std::endl;

	finalize_setup(pose);


	loop_centroid_relax(pose, the_cdr_loop_.start(), the_cdr_loop_.stop()   );


	TR<<"finish apply function !!!"<<std::endl;

	return;

}



///////////////////////////////////////////////////////////////////////////
/// @begin loop_centroid_relax
///
/// @brief actually relaxes the region specified
///
/// @detailed This is all done in low resolution. Intention was to give
///           camelid CDR H1 a larger perturbation.
///
/// @param[in] pose, loop begin position, loop end position
///
///
/// @authors Aroop 05/07/2010
///
/// @last_modified 05/07/2010
///////////////////////////////////////////////////////////////////////////
void RefineOneCDRLoopCentroid::loop_centroid_relax(
		pose::Pose & pose_in,
		Size const loop_begin,
		Size const loop_end )
{
	using namespace protocols;
	using namespace protocols::simple_moves;
	using namespace protocols::loops;
	using namespace protocols::moves;

	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	//setting MoveMap
	kinematics::MoveMapOP loop_map;
	loop_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	loop_map->clear();
	loop_map->set_chi( false );
	loop_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
	for( Size ii = loop_begin; ii <= loop_end; ii++ ) {
		allow_bb_move[ ii ] = true;
	}
	loop_map->set_bb( allow_bb_move );
	loop_map->set_jump( 1, false );


	Size loop_size = ( loop_end - loop_begin ) + 1;
	Size cutpoint = loop_begin + Size(loop_size/2);

	loops::Loop one_loop( loop_begin, loop_end,	cutpoint,	0, false );
	simple_one_loop_fold_tree( pose_in, one_loop );

	// set cutpoint variants for correct chainbreak scoring
	if( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
		if( !pose_in.residue( cutpoint ).has_variant_type(chemical::CUTPOINT_LOWER))
			core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, cutpoint );
		if( !pose_in.residue( cutpoint + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, cutpoint + 1 );
	}


	Real min_tolerance = 0.001;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
	bool nb_list = true;
	MinMoverOP loop_min_mover( new MinMover( loop_map, lowres_scorefxn_, min_type, min_tolerance, nb_list ) );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	Size inner_cycles( loop_size );
	Size outer_cycles( 1 );
	if( refine_input_loop_ ) {
		outer_cycles = 5;
	}
	if(  snug_fit_ ) {
		outer_cycles = 2;
	}
	if( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles = 1;
		outer_cycles = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
	BackboneMoverOP small_mover( new SmallMover( loop_map, high_move_temp, n_small_moves ) );
	BackboneMoverOP shear_mover( new ShearMover( loop_map, high_move_temp, n_small_moves ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );

	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CCDLoopClosureMoverOP ccd_moves( new CCDLoopClosureMover( one_loop, loop_map ) );
	RepeatMoverOP ccd_cycle( new RepeatMover(ccd_moves, n_small_moves) );

	SequenceMoverOP wiggle_cdr_centroid_loop_( new SequenceMover() );
	wiggle_cdr_centroid_loop_->add_mover( small_mover );
	wiggle_cdr_centroid_loop_->add_mover( shear_mover );
	wiggle_cdr_centroid_loop_->add_mover( ccd_cycle );


	loop_min_mover->apply( pose_in );

	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	MonteCarloOP mc;
	mc = MonteCarloOP( new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temperature ) );
	mc->reset( pose_in ); // monte carlo reset

	// outer cycle
	for(Size i = 1; i <= outer_cycles; i++) {
		mc->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc->set_temperature( temperature );
			wiggle_cdr_centroid_loop_->apply( pose_in );
			loop_min_mover->apply( pose_in );

			mc->boltzmann( pose_in );

		} // inner cycles
	} // outer cycles
	mc->recover_low( pose_in );


	// minimize
	if( !benchmark_ ) {
		loop_min_mover->apply( pose_in );
	}

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree


	return;
} // loop_centroid_relax





} // namespace antibody
} // namespace protocols



