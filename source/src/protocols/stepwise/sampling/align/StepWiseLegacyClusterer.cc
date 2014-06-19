// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseLegacyClusterer
/// @brief Not particularly fancy, just filters a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/align/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampling/util.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <list>

#ifdef PYROSETTA
	#include <time.h>
#endif

//Auto Headers
#include <core/id/AtomID.hh>
using namespace core;
using core::Real;

static basic::Tracer TR( "protocols.stepwise.sampling.align.StepWiseLegacyClusterer" ) ;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//////////////////////////////////////////////////////////////////////////
// This Clusterer will be deprecated soon in favor of an
//  'on-the-fly' clusterer for stepwise monte carlo, based on what
//  Parin piloted for RNA -- might be a little dependent on order
//  of sampling, but should be very memory efficient, and fits well
//  with new sample-and-screen framework which unifies RNA and protein.
//
//         -- rhiju, 2014
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampling {
namespace align {


  StepWiseLegacyClusterer::StepWiseLegacyClusterer( utility::vector1< PoseOP > const & pose_list ):
		input_pose_list_( pose_list )
	{
		initialize_parameters_and_input();
	}

	// convenience constructor, used in StepWiseProteinModeler. Maybe should go into a util.
	StepWiseLegacyClusterer::StepWiseLegacyClusterer( utility::vector1< PoseOP > const & pose_list,
																				utility::vector1< Size > const & moving_res_list,
																				modeler_options::StepWiseModelerOptionsCOP options,
																				bool const force_align ):
		input_pose_list_( pose_list )
	{
		initialize_parameters_and_input();
		if ( options->sampler_num_pose_kept() > 0 ) set_max_decoys( options->sampler_num_pose_kept() );
		set_cluster_by_all_atom_rmsd( options->cluster_by_all_atom_rmsd() ); // false by default
		set_rename_tags( true );
		set_force_align( force_align );
		set_calc_rms_res( moving_res_list );
		if ( options->cluster_rmsd() > 0.0 ) set_cluster_radius( options->cluster_rmsd() );
	}


  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseLegacyClusterer::~StepWiseLegacyClusterer()
  {}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseLegacyClusterer::initialize_parameters_and_input(){
		max_decoys_ = 400;
		cluster_radius_ = 0.1;
		cluster_by_all_atom_rmsd_ = true;
		score_diff_cut_ = 1000000.0;
		auto_tune_ = false;
		rename_tags_ = false;
		force_align_ = false;

		runtime_assert(  input_pose_list_.size() > 0 );
		sort( input_pose_list_.begin(), input_pose_list_.end(), sort_pose_by_score );
		score_min_ = total_energy_from_pose( *input_pose_list_[1] );
		input_pose_counter_ = 0;
		hit_score_cutoff_ = false;

		initialize_auto_tune_cluster_rmsds();
		initialized_atom_id_map_for_rmsd_ = false;
	}


  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseLegacyClusterer::cluster()
	{
		using namespace core::scoring;
		using namespace core::import_pose::pose_stream;
		using namespace core::chemical;
		using namespace core::pose;

		clock_t const time_start( clock() );

		// basic initialization
		initialize_cluster_list();

		if ( auto_tune_ ) {
			cluster_with_auto_tune();
		} else {
			do_some_clustering();
		}

		TR.Debug << "Total time in StepWiseLegacyClusterer: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::initialize_corresponding_atom_id_map( core::pose::Pose const & pose ){
		using namespace core::scoring;

		// Only need to do this once!!!
		if( cluster_by_all_atom_rmsd_ ) {
			setup_matching_heavy_atoms( pose,	pose, corresponding_atom_id_map_ );
		} else {
			setup_matching_protein_backbone_heavy_atoms( pose, pose, corresponding_atom_id_map_ );
		}
		initialized_atom_id_map_for_rmsd_ = true;
	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::do_some_clustering() {

		using namespace core::pose;

		hit_score_cutoff_ = false;
		while( input_pose_counter_ < input_pose_list_.size() ){
			input_pose_counter_++;
			PoseOP pose = input_pose_list_[ input_pose_counter_ ];
			Real score = total_energy_from_pose( *pose );
			if ( score > score_min_ + score_diff_cut_ ) {
				hit_score_cutoff_ = true;
				break;
			}

			std::string tag = tag_from_pose( *pose );
			TR.Debug << "Checking: " << tag << " with score " << score << " against list of size " << output_pose_list_.size();

			// carve out subset of residues for rms calculation.
			PoseOP rms_pose = pose;
			if ( calc_rms_res_.size() > 0 ) {
				rms_pose = new Pose;
				pdbslice( *rms_pose, *pose, calc_rms_res_ );
			}
			Size const found_close_cluster = check_for_closeness( rms_pose );

			if ( found_close_cluster == 0 )  {
				PoseOP pose_save = pose->clone();
				rms_pose_list_.push_back( rms_pose );
				output_pose_list_.push_back( pose );
				num_pose_in_cluster_.push_back( 1 );
				TR.Debug << " ... added. " << std::endl;

				if ( output_pose_list_.size() >= max_decoys_ ) break;
			} else {
				num_pose_in_cluster_[ found_close_cluster ]++;
				TR.Debug << " ... not added. " << std::endl;
			}
		}

	  TR.Debug << "After clustering, number of decoys: " << output_pose_list_.size() << std::endl;
		return;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::initialize_cluster_list() {
		rms_pose_list_.clear();
		output_pose_list_.clear();
		num_pose_in_cluster_.clear();
		input_pose_counter_ = 0;
		hit_score_cutoff_ = false;
	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::cluster_with_auto_tune() {

		for ( Size n = 1; n <= cluster_rmsds_to_try_with_auto_tune_.size(); n++ ) {

			cluster_radius_ = cluster_rmsds_to_try_with_auto_tune_[ n ];

			//can current cluster center list be shrunk, given the new cluster_radius cutoff?
			recluster_current_pose_list();

			do_some_clustering();

			if ( hit_score_cutoff_ ) {
				TR.Debug << "Hit score cutoff: " << score_diff_cut_ << std::endl;
				break;
			}
			if ( input_pose_counter_ == input_pose_list_.size() ){
				TR.Debug << "Done with pose list. " << std::endl;
				break;
			}
		}

		TR.Debug << "Clustering radius after auto_tune: " << cluster_radius_ << std::endl;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::recluster_current_pose_list() {

		utility::vector1< core::pose::PoseOP > old_rms_pose_list = rms_pose_list_;
		utility::vector1< core::pose::PoseOP > old_output_pose_list = output_pose_list_;
		utility::vector1< core::Size > old_num_pose_in_cluster = num_pose_in_cluster_;

		rms_pose_list_.clear();
		output_pose_list_.clear();

		for ( Size i = 1; i <= old_rms_pose_list.size(); i++ ) {

			core::pose::PoseOP rms_pose = old_rms_pose_list[ i ];

			Size const found_close_cluster = check_for_closeness( rms_pose );
			if ( found_close_cluster == 0 )  {
				rms_pose_list_.push_back( rms_pose );
				output_pose_list_.push_back(  old_output_pose_list[ i ]  );
				num_pose_in_cluster_.push_back( old_num_pose_in_cluster[ i ] );
			} else {
				num_pose_in_cluster_[ found_close_cluster ] += old_num_pose_in_cluster[ i ];
			}

		}

		TR.Debug <<  "After reclustering with rmsd " << cluster_radius_ << ", number of clusters reduced from " <<
			old_rms_pose_list.size() << " to " << rms_pose_list_.size() << std::endl;

	}


	///////////////////////////////////////////////////////////////
	Size
	StepWiseLegacyClusterer::check_for_closeness( core::pose::PoseOP rms_pose )
	{
		using namespace core::scoring;

		if ( !initialized_atom_id_map_for_rmsd_ ) initialize_corresponding_atom_id_map( *rms_pose );

		// go through the list backwards, because poses may be grouped by similarity --
		// the newest pose is probably closer to poses at the end of the list.
		for ( Size n = rms_pose_list_.size(); n >= 1; n-- ) {

			Real rmsd( 0.0 );

			if ( calc_rms_res_.size() == 0 || force_align_ ) {
				rmsd = rms_at_corresponding_atoms( *(rms_pose_list_[ n ]), *rms_pose, corresponding_atom_id_map_ );
			} else {
				// assumes prealignment of poses!!!
				rmsd = rms_at_corresponding_atoms_no_super( *(rms_pose_list_[ n ]), *rms_pose,
																										corresponding_atom_id_map_ );
			}

			if ( rmsd < cluster_radius_ ) return n;

		}
		return 0;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyClusterer::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseLegacyClusterer::initialize_auto_tune_cluster_rmsds(){
		cluster_rmsds_to_try_with_auto_tune_.clear();
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.1 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.2 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.25 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.3 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.4 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.6 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 0.8 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 1.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 1.2 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 1.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 1.75 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 2.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 2.25 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 2.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 2.75 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 3.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 3.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 4.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 4.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 5.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 6.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 7.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 8.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 9.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 10 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 12.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 15.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 17.5 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 20.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 25.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 30.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 35.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 40.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 45.0 );
		cluster_rmsds_to_try_with_auto_tune_.push_back( 50.0 );
	}


} //align
} //sampling
} //stepwise
} //protocols
