// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseClusterer
/// @brief Not particularly fancy, just filters a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/StepWiseClusterer.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableString.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <list>

using namespace core;
using core::Real;

static basic::Tracer TR( "protocols.swa.stepwise_clusterer" ) ;

namespace protocols {
namespace swa {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseClusterer::StepWiseClusterer( utility::vector1< std::string > const & silent_files_in ):
		max_decoys_( 400 ),
		cluster_radius_( 1.5 ),
		cluster_by_all_atom_rmsd_( true ),
		score_diff_cut_( 1000000.0 ),
		rename_tags_( false )
		//		scorefxn_( core::scoring::getScoreFunction() )
  {
		input_  = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		input_->set_order_by_energy( true );
		input_->set_record_source( true );
		input_->filenames( silent_files_in ); //triggers read in of files, too.
  }

  StepWiseClusterer::StepWiseClusterer( std::string const & silent_file_in ):
		max_decoys_( 400 ),
		cluster_radius_( 1.5 ),
		cluster_by_all_atom_rmsd_( true ),
		score_diff_cut_( 1000000.0 ),
		rename_tags_( false )
		//		scorefxn_( core::scoring::getScoreFunction() )
	{
		input_  = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		input_->set_order_by_energy( true );
		input_->set_record_source( true );

		utility::vector1< std::string > silent_files_;
		silent_files_.push_back( silent_file_in );
		input_->filenames( silent_files_ ); //triggers read in of files, too.
	}

  StepWiseClusterer::StepWiseClusterer( core::io::silent::SilentFileDataOP & sfd ):
		max_decoys_( 400 ),
		cluster_radius_( 1.5 ),
		cluster_by_all_atom_rmsd_( true ),
		score_diff_cut_( 1000000.0 ),
		rename_tags_( false )
		//		scorefxn_( core::scoring::getScoreFunction() )
	{
		input_  = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		input_->set_order_by_energy( true );
		//input_->set_record_source( true );
		input_->set_silent_file_data( sfd ); // triggers reordering by energy and all that.
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseClusterer::~StepWiseClusterer()
  {}


	///////////////////////////////////////////////////////////////
	bool
	StepWiseClusterer::check_for_closeness( core::pose::PoseOP const & pose_op )
	{
		using namespace core::scoring;

		// go through the list backwards, because poses may be grouped by similarity --
		// the newest pose is probably closer to poses at the end of the list.
		for ( Size n = pose_output_list_.size(); n >= 1; n-- ) {

			Real rmsd( 0.0 );

			//			std::cout << "checking agaist pose " << n << " out of " << pose_output_list_.size() << std::endl;
			//			std::cout << " pose_output_list_[ n ]->size " << pose_output_list_[ n ]->total_residue() << std::endl;
			//			std::cout << " pose_op->size " << pose_op->total_residue() << std::endl;

			if (cluster_by_all_atom_rmsd_ ) {
				rmsd = all_atom_rmsd( *(pose_output_list_[ n ]), *pose_op );
			} else {
				rmsd = rmsd_with_super( *(pose_output_list_[ n ]), *pose_op, core::scoring::is_protein_backbone_including_O );
			}

			if ( rmsd < cluster_radius_ )	return false;
		}
		return true;
	}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseClusterer::cluster()
	{
		using namespace core::scoring;
		using namespace core::import_pose::pose_stream;
		using namespace core::chemical;
		using namespace core::pose;

		// setup residue types
		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		/////////////////////////////
		// Read in all the poses.

		Size count( 0 );
		pose_output_list_.clear();
		tag_output_list_.clear();

		Real score_min( 0.0 );
		bool score_min_defined( false );

		while ( input_->has_another_pose() ) {

			PoseOP pose_op( new Pose );
			core::io::silent::SilentStructOP silent_struct( input_->next_struct() );
			silent_struct->fill_pose( *pose_op, *rsd_set );

			core::Real score( 0.0 );
			getPoseExtraScores( *pose_op, "score", score );

			if ( !score_min_defined ){
				score_min = score;
				score_min_defined = true;
			}

			if ( score > score_min + score_diff_cut_ ) break;

			std::string tag( silent_struct->decoy_tag() );
			TR << "CHECKING " << tag << " with score " << score << " against list of size " << pose_output_list_.size() << std::endl;

			bool const OK = check_for_closeness( pose_op );
			if ( OK )  {
				count++;
				TR << "ADDING " << tag << std::endl;
				if ( count > max_decoys_ ) break;
				tag_output_list_.push_back(  tag );
				pose_output_list_.push_back(  pose_op );
				silent_struct_output_list_.push_back(  silent_struct  );
			}
		}

	  TR << "After clustering, number of decoys: " << pose_output_list_.size() << std::endl;

	}

	/////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseClusterer::output_silent_file( std::string const & silent_file ){

		using namespace core::io::silent;

		SilentFileData silent_file_data;

		for ( Size n = 1 ; n <= silent_struct_output_list_.size(); n++ ) {

			SilentStructOP & s( silent_struct_output_list_[ n ] );

			if ( rename_tags_ ){
				s->add_comment( "PARENT_TAG", s->decoy_tag() );
				std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
				s->set_decoy_tag( tag );
			}

			silent_file_data.write_silent_struct( *s, silent_file, false /*write score only*/ );

		}

	}

	//////////////////////////////////////////////
	PoseList
	StepWiseClusterer::clustered_pose_list(){

		PoseList pose_list;

		for ( Size n = 1 ; n <= pose_output_list_.size(); n++ ) {
			pose_list[ tag_output_list_[n] ] = pose_output_list_[ n ];
		}

		return pose_list;
	}


}
}
