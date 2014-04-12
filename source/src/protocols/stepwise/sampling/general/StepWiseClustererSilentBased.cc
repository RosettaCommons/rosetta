// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseClustererSilentBased
/// @brief Not particularly fancy, just filters a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/general/StepWiseClustererSilentBased.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.hh>
#include <protocols/stepwise/StepWiseUtil.hh>

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

static basic::Tracer TR( "protocols.stepwise.StepWiseClustererSilentBased" ) ;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//////////////////////////////////////////////////////////////////////////
//
// This is a silent-file-based clusterer originally developed for
//  stepwise assembly of proteins.
// I decided to shift back to
//  a clusterer that saves full poses, to simplify the data structures
//  that we were keeping track of -- see StepWiseClusterer.
//  That requires more memory, however.
// Both of these Clusterers will be deprecated soon in favor of an
//  'on-the-fly' clusterer for stepwise monte carlo, based on what
//  Parin piloted for RNA -- might be a little dependent on order
//  of sampling, but should be very memory efficient.
//
//         -- rhiju, 2014
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampling {
namespace general {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseClustererSilentBased::StepWiseClustererSilentBased( utility::vector1< std::string > const & silent_files_in )
  {
		initialize_parameters_and_input();
		input_->set_record_source( true );
		input_->filenames( silent_files_in ); //triggers read in of files, too.
  }

  StepWiseClustererSilentBased::StepWiseClustererSilentBased( std::string const & silent_file_in )
	{
		initialize_parameters_and_input();
		input_->set_record_source( true );

		utility::vector1< std::string > silent_files_;
		silent_files_.push_back( silent_file_in );
		input_->filenames( silent_files_ ); //triggers read in of files, too.
	}

  StepWiseClustererSilentBased::StepWiseClustererSilentBased( core::io::silent::SilentFileDataOP & sfd )
	{
		initialize_parameters_and_input();
		input_->set_silent_file_data( sfd ); // triggers reordering by energy and all that.
	}

	// convenience constructor, used in StepWiseProteinModeler. Maybe should go into a util.
	StepWiseClustererSilentBased::StepWiseClustererSilentBased( core::io::silent::SilentFileDataOP silent_file_data,
																				utility::vector1< Size > const & moving_res_list,
																				protein::StepWiseProteinModelerOptionsCOP options,
																				bool const force_align )
	{
		initialize_parameters_and_input();
		input_->set_silent_file_data( silent_file_data ); // triggers reordering by energy and all that.
		set_max_decoys( options->max_decoys() );
		set_cluster_by_all_atom_rmsd( options->cluster_by_all_atom_rmsd() ); // false by default
		set_rename_tags( true );
		set_force_align( force_align );
		set_calc_rms_res( moving_res_list );
		Real cluster_radius = options->rescore_only() ? 0.0 : options->cluster_radius();
		set_cluster_radius( cluster_radius );
	}


  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseClustererSilentBased::~StepWiseClustererSilentBased()
  {}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseClustererSilentBased::initialize_parameters_and_input(){
		max_decoys_ = 400;
		cluster_radius_ = 1.5;
		cluster_by_all_atom_rmsd_ = true;
		score_diff_cut_ = 1000000.0;
		auto_tune_ = false;
		rename_tags_ = false;
		force_align_ = false;

		score_min_ =  0.0 ;
		score_min_defined_ = false;

		input_  = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		input_->set_order_by_energy( true );

		initialize_auto_tune_cluster_rmsds();
		hit_score_cutoff_ = false;
		initialized_atom_id_map_for_rmsd_ = false;

		rsd_type_set_ = option[ in::file::residue_type_set ]();
	}


  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseClustererSilentBased::cluster()
	{
		using namespace core::scoring;
		using namespace core::import_pose::pose_stream;
		using namespace core::chemical;
		using namespace core::pose;

		clock_t const time_start( clock() );

		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( rsd_type_set_ );

		// basic initialization
		initialize_cluster_list();

		if ( auto_tune_ ) {
			cluster_with_auto_tune();
		} else {
			do_some_clustering();
		}

		TR.Debug << "Total time in StepWiseClustererSilentBased: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::initialize_corresponding_atom_id_map( core::pose::Pose const & pose ){
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
	StepWiseClustererSilentBased::do_some_clustering() {

		using namespace core::pose;

		hit_score_cutoff_ = false;

		//for loop modeling, little moving_element of pose used to calculate rms -- and saved.
		PoseOP pose_op( new Pose );
		Pose & pose = *pose_op;

		while ( input_->has_another_pose() ) {

			core::io::silent::SilentStructOP silent_struct( input_->next_struct() );
			silent_struct->fill_pose( pose );

			Real score( 0.0 );
			getPoseExtraScores( pose, "score", score );

			if ( !score_min_defined_ ){
				score_min_ = score;
				score_min_defined_ = true;
			}

			if ( score > score_min_ + score_diff_cut_ ) {
				hit_score_cutoff_ = true;
				break;
			}

			std::string tag( silent_struct->decoy_tag() );
			TR.Debug << "Checking: " << tag << " with score " << score << " against list of size " << pose_output_list_.size();

			// carve out subset of residues for rms calculation.
			if ( calc_rms_res_.size() > 0 )	pdbslice( pose, calc_rms_res_ );

			Size const found_close_cluster = check_for_closeness( pose_op );

			if ( found_close_cluster == 0 )  {
				PoseOP pose_save = pose.clone();
				tag_output_list_.push_back(  tag );
				pose_output_list_.push_back(  pose_save );
				silent_struct_output_list_.push_back(  silent_struct  );
				num_pose_in_cluster_.push_back( 1 );
				TR.Debug << " ... added. " << std::endl;
				if ( pose_output_list_.size() >= max_decoys_ ) break;
			} else{
				num_pose_in_cluster_[ found_close_cluster ]++;
				TR.Debug << " ... not added. " << std::endl;
			}
		}

	  TR.Debug << "After clustering, number of decoys: " << pose_output_list_.size() << std::endl;
		return;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::initialize_cluster_list() {

		pose_output_list_.clear();
		tag_output_list_.clear();
		silent_struct_output_list_.clear();
		num_pose_in_cluster_.clear();

		score_min_ =  0.0 ;
		score_min_defined_ = false;

		hit_score_cutoff_ = false;
	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::cluster_with_auto_tune() {

		for ( Size n = 1; n <= cluster_rmsds_to_try_with_auto_tune_.size(); n++ ) {

			cluster_radius_ = cluster_rmsds_to_try_with_auto_tune_[ n ];

			//can current cluster center list be shrunk, given the new cluster_radius cutoff?
			recluster_current_pose_list();

			do_some_clustering();

			if ( hit_score_cutoff_ ) {
				TR.Debug << "Hit score cutoff: " << score_diff_cut_ << std::endl;
				break;
			}
			if ( !input_->has_another_pose() ) {
				TR.Debug << "Done with pose list. " << std::endl;
				break;
			}
		}

		TR.Debug << "Clustering radius after auto_tune: " << cluster_radius_ << std::endl;

	}


	/////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::recluster_current_pose_list() {

		utility::vector1< core::pose::PoseOP > old_pose_output_list = pose_output_list_;
		utility::vector1< std::string > old_tag_output_list = tag_output_list_;
		utility::vector1< core::io::silent::SilentStructOP > old_silent_struct_output_list = silent_struct_output_list_;
		utility::vector1< core::Size > old_num_pose_in_cluster = num_pose_in_cluster_;

		pose_output_list_.clear();
		tag_output_list_.clear();
		silent_struct_output_list_.clear();

		for ( Size i = 1; i <= old_pose_output_list.size(); i++ ) {

			core::pose::PoseOP pose_op = old_pose_output_list[ i ];

			Size const found_close_cluster = check_for_closeness( pose_op );
			if ( found_close_cluster == 0 )  {
				tag_output_list_.push_back(  old_tag_output_list[ i ] );
				pose_output_list_.push_back(  old_pose_output_list[ i ] );
				silent_struct_output_list_.push_back(  old_silent_struct_output_list[ i ]  );
				num_pose_in_cluster_.push_back( old_num_pose_in_cluster[ i ] );
			} else {
				num_pose_in_cluster_[ found_close_cluster ] += old_num_pose_in_cluster[ i ];
			}

		}

		TR.Debug <<  "After reclustering with rmsd " << cluster_radius_ << ", number of clusters reduced from " <<
			old_pose_output_list.size() << " to " << pose_output_list_.size() << std::endl;

	}


	///////////////////////////////////////////////////////////////
	Size
	StepWiseClustererSilentBased::check_for_closeness( core::pose::PoseOP const & pose_op )
	{
		using namespace core::scoring;

		if ( !initialized_atom_id_map_for_rmsd_ ) initialize_corresponding_atom_id_map( *pose_op );

		// go through the list backwards, because poses may be grouped by similarity --
		// the newest pose is probably closer to poses at the end of the list.
		for ( Size n = pose_output_list_.size(); n >= 1; n-- ) {

			Real rmsd( 0.0 );


			if ( calc_rms_res_.size() == 0 || force_align_ ) {
				rmsd = rms_at_corresponding_atoms( *(pose_output_list_[ n ]), *pose_op, corresponding_atom_id_map_ );
			} else {
				// assumes prealignment of poses!!!
				rmsd = rms_at_corresponding_atoms_no_super( *(pose_output_list_[ n ]), *pose_op,
																										corresponding_atom_id_map_ );
			}

			if ( rmsd < cluster_radius_ )	{
				return n;
			}
		}
		return 0;
	}

	core::io::silent::SilentStructOP
	StepWiseClustererSilentBased::setup_silent_struct( Size const n ) {
		using namespace core::io::silent;
		SilentStructOP & s( silent_struct_output_list_[ n ] );
		s->add_string_value( "nclust", ObjexxFCL::format::I(8,num_pose_in_cluster_[ n ]) );
		if ( rename_tags_ ){
			s->add_comment( "PARENT_TAG", s->decoy_tag() );
			std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
			s->set_decoy_tag( tag );
		}
		return s;
	}


	/////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::output_silent_file( std::string const & silent_file ){
		using namespace core::io::silent;

		SilentFileData silent_file_data;
		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {
			SilentStructOP s = setup_silent_struct( n );
			silent_file_data.write_silent_struct( *s, silent_file, false /*write score only*/ );
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP
	StepWiseClustererSilentBased::silent_file_data(){
		using namespace core::io::silent;

		SilentFileDataOP silent_file_data = new SilentFileData;
		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {
			silent_file_data->add_structure( silent_struct_output_list_[ n ] );
		}
		return silent_file_data;
	}

	//////////////////////////////////////////////
	PoseList
	StepWiseClustererSilentBased::clustered_pose_list(){

		PoseList pose_list;

		for ( Size n = 1 ; n <= pose_output_list_.size(); n++ ) {
			pose_list[ tag_output_list_[n] ] = pose_output_list_[ n ];
		}

		return pose_list;
	}

  //////////////////////////////////////////////////////////////////////////
	// new -- this will be deprecated after we transition StepWiseClustererSilentBased
	// to act on poses (instead of on silentstructs)... may end up being a
	// separate object.
  //////////////////////////////////////////////////////////////////////////
	utility::vector1< PoseOP >
	StepWiseClustererSilentBased::get_pose_list() {
		using namespace core::io::silent;
		utility::vector1< PoseOP > pose_list;
		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {
			SilentStructOP s = setup_silent_struct( n );
			PoseOP pose = new Pose;
			s->fill_pose( *pose );
			tag_into_pose( *pose, s->decoy_tag() );
			pose_list.push_back( pose );
		}
		return pose_list;
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseClustererSilentBased::set_silent_file_data( core::io::silent::SilentFileDataOP & sfd ){
		input_->set_silent_file_data( sfd );
	}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseClustererSilentBased::initialize_auto_tune_cluster_rmsds(){
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


} //general
} //sampling
} //stepwise
} //protocols
