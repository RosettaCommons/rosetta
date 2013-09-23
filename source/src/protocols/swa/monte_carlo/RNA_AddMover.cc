// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover
/// @brief Adds an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_TorsionMover.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>

#include <map>

using namespace core;
using core::Real;
using utility::make_tag_with_dashes;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose full_model_info object.
//////////////////////////////////////////////////////////////////////////

static numeric::random::RandomGenerator RG(2555512);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_add_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddMover::RNA_AddMover( scoring::ScoreFunctionOP scorefxn ):
		scorefxn_( scorefxn ),
		presample_added_residue_( true ),
		presample_by_swa_( false ),
		minimize_single_res_( false ),
		start_added_residue_in_aform_( false ),
		internal_cycles_( 50 ),
		rna_torsion_mover_( new RNA_TorsionMover ),
		sample_range_small_( 5.0 ),
		sample_range_large_( 40.0 ),
		kT_( 0.5 )
	{}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_AddMover::~RNA_AddMover()
  {}

	////////////////////////////////////////////////////////////////////
	std::string
	RNA_AddMover::get_name() const {
		return "RNA_AddMover";
	}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_AddMover::apply( core::pose::Pose &  )
	{
		std::cout << "not defined yet" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_AddMover::apply( core::pose::Pose & pose, Size const res_to_add_in_full_model_numbering, Size const res_to_build_off_in_full_model_numbering )
	{

		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::swa::rna;

		runtime_assert( pose.total_residue() > 1 );
		ResidueTypeSet const & rsd_set = pose.residue_type( 1 ).residue_type_set();

		Size suite_num( 0 ), nucleoside_num( 0 ); // will record which new dofs added.
		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		std::string const & full_sequence  = full_model_info.full_sequence();

		runtime_assert( res_list.has_value( res_to_build_off_in_full_model_numbering ) );
		runtime_assert( !res_list.has_value( res_to_add_in_full_model_numbering ) );

		Size const res_to_build_off = res_list.index( res_to_build_off_in_full_model_numbering );

		// need to encapsulate the following residue addition & domain addition functions...

		if ( res_to_add_in_full_model_numbering == (res_to_build_off_in_full_model_numbering + 1) ){
			// addition to strand ending (append)
			runtime_assert( res_to_build_off == pose.total_residue() || res_list[ res_to_build_off ] < res_list[ res_to_build_off+1 ] -1 );
			runtime_assert( res_list[ res_to_build_off ] < full_sequence.size() );

			Size const res_to_add = res_to_build_off + 1;
			Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( res_to_add_in_full_model_numbering );

			if ( other_pose_idx ){ // addition of a domain (a whole sister pose)

				Pose & other_pose = *(full_model_info.other_pose_list()[ other_pose_idx ]);
				Size const res_to_build_off_in_full_model_numbering = res_list[ res_to_build_off ];
				merge_in_other_pose( pose, other_pose, res_to_build_off_in_full_model_numbering );

				nonconst_full_model_info( pose ).remove_other_pose_at_idx( other_pose_idx );

				suite_num = get_res_list_from_full_model_info( pose ).index( res_to_build_off_in_full_model_numbering );
				nucleoside_num = 0; // don't sample sugar pucker or side chain -- that will screw up the (fixed) domain structure.

			} else { // single residue addition -- can this just be combine with above?

				char newrestype = full_sequence[ res_to_add_in_full_model_numbering - 1 ];
				choose_random_if_unspecified_nucleotide( newrestype );

				chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );
				chemical::ResidueTypeCOPs const & rsd_type_list( rsd_set.aa_map( my_aa ) );

				// iterate over rsd_types, pick one.
				chemical::ResidueType const & rsd_type = *rsd_type_list[1];
				core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

				remove_variant_type_from_pose_residue( pose, "UPPER_TERMINUS", res_to_build_off ); // got to be safe.

				pose.append_polymer_residue_after_seqpos( *new_rsd, res_to_build_off, true /*build ideal geometry*/ );

				reorder_full_model_info_after_append( pose, res_to_add );

				suite_num = res_to_add - 1;
				nucleoside_num = res_to_add;
			}
		} else {

			// addition to strand beginning (prepend)
			runtime_assert( res_to_add_in_full_model_numbering == (res_to_build_off_in_full_model_numbering - 1) );

			Size const res_to_add = res_to_build_off;
			Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( res_to_add_in_full_model_numbering );

			TR << "About to add onto " << res_to_build_off_in_full_model_numbering << " the following residue (in full model numbering) " << res_to_add_in_full_model_numbering << " which may be part of other pose " << other_pose_idx << std::endl;

			if ( other_pose_idx ){ // addition of a domain (a whole sister pose)

				Pose & other_pose = *(full_model_info.other_pose_list()[ other_pose_idx ]);
				merge_in_other_pose( pose, other_pose, res_to_build_off_in_full_model_numbering - 1 /*merge_res*/ );

				nonconst_full_model_info( pose ).remove_other_pose_at_idx( other_pose_idx );

				suite_num = get_res_list_from_full_model_info( pose ).index( res_to_add_in_full_model_numbering );
				nucleoside_num = 0; // don't sample sugar pucker or side chain -- that will screw up the (fixed) domain structure.

			} else {  // single residue addition -- can this just be combine with above?

				runtime_assert( res_list[ res_to_add ] > 1 );

				char newrestype = full_sequence[ res_to_add_in_full_model_numbering - 1 ];
				choose_random_if_unspecified_nucleotide( newrestype );

				TR << "I want to add: " << newrestype << " before " << res_to_build_off << std::endl;

				chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );

				chemical::ResidueTypeCOPs const & rsd_type_list( rsd_set.aa_map( my_aa ) );

				// iterate over rsd_types, pick one.
				chemical::ResidueType const & rsd_type = *rsd_type_list[1];

				core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

				remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_add ); // got to be safe.
				remove_variant_type_from_pose_residue( pose, "LOWER_TERMINUS", res_to_add ); // got to be safe.

				pose.prepend_polymer_residue_before_seqpos( *new_rsd, res_to_add, true /*build ideal geometry*/ );

				reorder_full_model_info_after_prepend( pose, res_to_add );

				// initialize with a random torsion... ( how about an A-form + perturbation ... or go to a 'reasonable' rotamer)
				suite_num = res_to_add;
				nucleoside_num = res_to_add;
			}

		}

		fix_up_residue_type_variants( pose );

		if ( start_added_residue_in_aform_ ){
			rna_torsion_mover_->apply_suite_torsion_Aform( pose, suite_num );
			if ( nucleoside_num > 0 ) rna_torsion_mover_->apply_nucleoside_torsion_Aform( pose, nucleoside_num );
		} else {
			rna_torsion_mover_->apply_random_suite_torsion( pose, suite_num );
			if ( nucleoside_num > 0 ) rna_torsion_mover_->apply_random_nucleoside_torsion( pose, nucleoside_num );
		}

		rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_large_);
		if ( nucleoside_num > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large_);

		///////////////////////////////////
		// Presampling added residue
		///////////////////////////////////
		if ( presample_added_residue_ ){
			if ( presample_by_swa_ ){
				Size swa_sample_res = nucleoside_num;
				if ( swa_sample_res == 0) swa_sample_res = suite_num; // nucleoside_num = 0 in domain addition.
				sample_by_swa( pose, swa_sample_res );
			} else {
				sample_by_monte_carlo_internal( pose, nucleoside_num, suite_num );
			}
		}

	}


	////////////////////////////////////////////////////////////////////
	void
	RNA_AddMover::sample_by_swa( pose::Pose & pose, Size const res_to_add ) const{

		using namespace core::pose::full_model_info;

		runtime_assert( stepwise_rna_modeler_ != 0 );

		TR.Debug << "presampling by swa " << res_to_add << std::endl;
		stepwise_rna_modeler_->set_moving_res_and_reset( res_to_add );
		if ( !minimize_single_res_ ) stepwise_rna_modeler_->set_minimize_res( get_moving_res_from_full_model_info( pose ) );
		stepwise_rna_modeler_->apply( pose );

	}


	////////////////////////////////////////////////////////////////////
	void
	RNA_AddMover::sample_by_monte_carlo_internal( pose::Pose &  pose, Size const nucleoside_num, Size const suite_num ) const {

		using namespace protocols::moves;

		TR.Debug << "presampling added residue! " << nucleoside_num << " over " << internal_cycles_ << " cycles " << std::endl;
		MonteCarloOP monte_carlo_internal = new MonteCarlo( pose, *scorefxn_, kT_ );

		std::string move_type( "" );
		for ( Size count_internal = 1; count_internal <= internal_cycles_; count_internal++ ){

			rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_large_);
			if ( nucleoside_num > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large_);
			monte_carlo_internal->boltzmann( pose, move_type );

			rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_small_);
			if ( nucleoside_num > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_small_);
			monte_carlo_internal->boltzmann( pose, move_type );
			//std::cout << "During presampling: " << (*scorefxn_)( pose );
		} // monte carlo cycles
	}

	///////////////////////////////////////////////////////////////////
	void
	RNA_AddMover::set_stepwise_rna_modeler( protocols::swa::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler ){
		stepwise_rna_modeler_ = stepwise_rna_modeler;
	}


}
}
}
