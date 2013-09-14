// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/RNA_ResampleMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/monte_carlo/RNA_ResampleMover.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/TransientCutpointHandler.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <protocols/swa/monte_carlo/SWA_MonteCarloUtil.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>


static numeric::random::RandomGenerator RG(239145021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.swa.monte_carlo.RNA_ResampleMover" );

namespace protocols {
namespace swa {
namespace monte_carlo {

	//Constructor
	RNA_ResampleMover::RNA_ResampleMover(	core::scoring::ScoreFunctionOP scorefxn ):
		scorefxn_( scorefxn ),
		just_min_after_mutation_frequency_( 0.5 ),
		allow_internal_moves_( false ),
		num_random_samples_( 20 ),
		use_phenix_geo_( true ),
		erraser_( true ),
		minimize_single_res_( false )
	{}

	//Destructor
	RNA_ResampleMover::~RNA_ResampleMover()
	{}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ResampleMover::apply( pose::Pose & pose,
														std::string & move_type ){

		using namespace protocols::swa;
		using namespace protocols::swa::rna;
		using namespace protocols::swa::monte_carlo;
		using namespace core::pose::full_model_info;

		Size remodel_res( 0 );
		utility::vector1< Size > const moving_res = get_moving_res_from_full_model_info( pose );

		bool is_terminal( true );
		if ( allow_internal_moves_ ){

			utility::vector1< utility::vector1< Size > > const moving_chunks = get_moving_chunks_from_full_model_info( pose );
			if ( moving_chunks.size() == 0 ) return false;
			utility::vector1< Size > remodel_chunk =	RG.random_element( moving_chunks );


		} else {
			utility::vector1< SWA_Move > swa_moves;
			get_potential_resample_chunks( pose, swa_moves );
			if ( swa_moves.size() == 0 ) return false;
			SWA_Move const & swa_move = RG.random_element( swa_moves );

			utility::vector1< Size > const & remodel_chunk = swa_move.chunk();
			MovingResidueCase const & moving_residue_case = swa_move.moving_residue_case();
			runtime_assert( swa_move.add_or_delete_choice() == NO_ADD_OR_DELETE );

			TR << "Chose move: " << swa_move << std::endl;

			if ( remodel_chunk.size() == 1 ) { //a single residue hanging off the end
				remodel_res = remodel_chunk[ 1 ];
			} else { // define the suite of internal residue that will move.
				if ( moving_residue_case == CHAIN_TERMINUS_3PRIME ){
					remodel_res = remodel_chunk[ 1 ] - 1;
				} else{
					runtime_assert( moving_residue_case == CHAIN_TERMINUS_5PRIME );
					remodel_res = remodel_chunk[ remodel_chunk.size() ];
				}
			}
		}
		TR << "About to remodel residue " << remodel_res << std::endl;

		bool const did_mutation = mutate_res_if_allowed( pose, remodel_res ); // based on 'n' in full_model_info.full_sequence
		bool just_min_after_mutation_ = ( did_mutation && ( RG.uniform() < just_min_after_mutation_frequency_ ) );

		swa::rna::StepWiseRNA_Modeler stepwise_rna_modeler( remodel_res, scorefxn_ );
		stepwise_rna_modeler.set_use_phenix_geo( use_phenix_geo_ );
		stepwise_rna_modeler.set_force_centroid_interaction( true );
		stepwise_rna_modeler.set_choose_random( true );
		stepwise_rna_modeler.set_num_random_samples( num_random_samples_ );
		stepwise_rna_modeler.set_num_pose_minimize( 1 );

		if ( just_min_after_mutation_ ) stepwise_rna_modeler.set_skip_sampling( true );
		if ( ! minimize_single_res_ ) stepwise_rna_modeler.set_minimize_res( moving_res );

		if ( is_at_terminus( pose, remodel_res ) || !allow_internal_moves_ ){

			// need to update is_at_terminus to instead check/assert if *chunk* is at terminus.

			move_type = "swa";
			stepwise_rna_modeler.apply( pose );

		} else {
			runtime_assert( allow_internal_moves_ );
			move_type = "internal";

			stepwise_rna_modeler.set_kic_sampling( erraser_  );

			TransientCutpointHandler cutpoint_handler( remodel_res );
			if ( ! minimize_single_res_ ) cutpoint_handler.set_minimize_res( moving_res );
			cutpoint_handler.put_in_cutpoints( pose );

			stepwise_rna_modeler.apply( pose );

			cutpoint_handler.take_out_cutpoints( pose );
		}

		if ( did_mutation ) move_type += "-mut";
		if ( just_min_after_mutation_ ) move_type = "mut";

		return true;

	}


} //monte_carlo
} //swa
} //protocols
