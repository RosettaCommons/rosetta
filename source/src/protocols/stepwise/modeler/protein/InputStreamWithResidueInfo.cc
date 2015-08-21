// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinPoseSetup
/// @brief Sets up pose and job parameters for protein or RNA stepwise building.
/// @details
/// @author Rhiju Das
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

//////////////////////////////////
#include <core/chemical/ResidueType.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>
#include <core/import_pose/pose_stream/ExtendedPoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/stepwise/modeler/align/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// RNA stuff
#include <protocols/farna/util.hh>

#include <utility/exit.hh>
#include <string>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;
using core::Real;

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.protein.InputStreamWithResidueInfo" );

///////////////////////////////////////////////////////////////////////////////////
//
// This class was setup in 2009 for readin of protein silent/PDB files, including
// potential templates whose residues do not line up with target (for homology
// modeling).
//
// Later also set up ChunkLibrary (used in FragmentLibrary -- uses MiniPose for
//  memory savings) and CopyDofStepWiseSampler (used in stepwise sample-and-screen --
//  some useful cacheing). Probably should choose one of these, and get rid of the
//  rest -- need to define use cases & decide.
//
// Had to add 'advance_to_next_pose_segment', 'apply_current_pose_segment' operations
//  to enable wiring up to stepwise sample-and-screen code as InputStreamStepWiseSampler, so this
//  is really getting unwieldy.
//
//      -- rhiju, 2014
//
//////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

//////////////////////////////////////////////////////////////////////
// This is for file readin.
//////////////////////////////////////////////////////////////////////
import_pose::pose_stream::PoseInputStreamOP
setup_pose_input_stream(
	utility::options::StringVectorOption const & option_s1,
	utility::options::StringVectorOption const & option_silent1,
	utility::options::StringVectorOption const & option_tags1
){
	using namespace import_pose::pose_stream;

	PoseInputStreamOP input1;

	if ( option_s1.user() ) {
		// pdb input(s).
		input1 = PoseInputStreamOP( new PDBPoseInputStream( option_s1() ) );

	} else if ( option_silent1.size() > 0 ) {

		if ( option_tags1.user() > 0 ) {
			input1 = PoseInputStreamOP( new SilentFilePoseInputStream( option_silent1() ,
				option_tags1() ) );
		} else {
			input1 = PoseInputStreamOP( new SilentFilePoseInputStream( option_silent1() ) );
		}
	} else {
		// create a pose stream with a single blank pose...
		input1 = PoseInputStreamOP( new ExtendedPoseInputStream( "", 1 ) ); // hmm...
	}

	return input1;

}

////////////////////////////////////////////////////////////////////////////////
void
initialize_input_streams(   utility::vector1< protocols::stepwise::modeler::protein::InputStreamWithResidueInfoOP > & input_streams ){

	using namespace protocols::stepwise;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::stepwise;
	using namespace import_pose::pose_stream;

	// my options.
	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;


	if ( option[ s1 ].user()  || option[ silent1 ].user()  ) { // assume new style of input.
		utility::vector1< Size > slice_res_1 = blank_size_vector;
		if ( option[ slice_res1 ].user() ) slice_res_1 =  option[ slice_res1]();
		InputStreamWithResidueInfoOP stream1( new InputStreamWithResidueInfo(
			setup_pose_input_stream( option[ s1 ], option[ silent1 ], option[ tags1 ] ),
			option[ input_res1 ](),
			slice_res_1 ) );
		if ( option[ backbone_only1 ]() ) stream1->set_backbone_only( true );
		input_streams.push_back( stream1 );

		if ( option[ input_res2 ].user() ) {
			utility::vector1< Size > slice_res_2 = blank_size_vector;
			if ( option[ slice_res2 ].user() ) slice_res_2 =  option[ slice_res2]();
			InputStreamWithResidueInfoOP stream2( new InputStreamWithResidueInfo(
				setup_pose_input_stream( option[ s2 ], option[ silent2 ], option[ tags2 ] ),
				option[ input_res2 ](),
				slice_res_2 ) );
			if ( option[ backbone_only2 ]() ) stream2->set_backbone_only( true );
			input_streams.push_back( stream2 );
		}
	} else if ( option[ in::file::input_res ].user() ) {  //old style
		utility::vector1< std::string > silent_files_in, pdb_tags;

		// First read in any information on pdb read in from silent files.
		if ( option[ in::file::silent ].user() ) {
			silent_files_in = option[ in::file::silent ]();

			if ( option[ in::file::tags ].user() ) {
				pdb_tags = option[ in::file::tags ]();
			}
			if ( pdb_tags.size() < 1 ) return; //early finish!! no combo!
		}
		if ( option[ in::file::s ].user() ) {
			// Then any pdbs that need to be read in from disk.
			utility::vector1< std::string > const pdb_tags_from_disk( option[ in::file::s ]() );
			for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) pdb_tags.push_back( pdb_tags_from_disk[ n ] );
		}
		assert( pdb_tags.size() > 0 );
		initialize_input_streams_with_residue_info( input_streams,
			pdb_tags, silent_files_in,
			option[ in::file::input_res ](), option[ input_res2 ]() );
	}

}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
InputStreamWithResidueInfo::InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP pose_input_stream,
	utility::vector1< Size > const & input_res,
	utility::vector1< Size > const & slice_res ):
	pose_input_stream_( pose_input_stream ),
	input_res_( input_res ),
	slice_res_( slice_res ),
	backbone_only_( false )
{
	initialize_defaults();
}

//////////////////////////////////////////////////////////////////////////
InputStreamWithResidueInfo::InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP pose_input_stream,
	utility::vector1< Size > const & input_res ):
	pose_input_stream_( pose_input_stream ),
	input_res_( input_res ),
	backbone_only_( false )
{
	initialize_defaults();
}


//////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::initialize_defaults(){
	if ( slice_res_.size() == 0 ) {
		for ( Size i = 1; i <= input_res_.size(); i++ ) slice_res_.push_back( i );
	}
	for ( Size i = 1; i <= input_res_.size(); i++ ) full_to_sub_[ input_res_[i] ] = input_res_[i];
}


//////////////////////////////////////////////////////////////////////////
InputStreamWithResidueInfo::~InputStreamWithResidueInfo(){}
//////////////////////////////////////////////////////////////////////////
import_pose::pose_stream::PoseInputStreamOP &
InputStreamWithResidueInfo::pose_input_stream(){ return pose_input_stream_; }
//////////////////////////////////////////////////////////////////////////
utility::vector1< Size > const &
InputStreamWithResidueInfo::input_res(){ return input_res_; }
//////////////////////////////////////////////////////////////////////////
utility::vector1< Size > const &
InputStreamWithResidueInfo::slice_res(){ return slice_res_; }
//////////////////////////////////////////////////////////////////////////
std::map< Size, Size > &
InputStreamWithResidueInfo::full_to_sub(){ return full_to_sub_; }
//////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::set_slice_res(  utility::vector1< Size > const & slice_res ){ slice_res_ = slice_res; }
//////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::set_full_to_sub( std::map< Size, Size > const & full_to_sub ){ full_to_sub_ = full_to_sub; }
//////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::set_rsd_set( chemical::ResidueTypeSetCAP & rsd_set ){
	rsd_set_ = rsd_set;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::reset(){
	pose_input_stream_->reset();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
InputStreamWithResidueInfo::has_another_pose() const {
	return ( pose_input_stream_->has_another_pose() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::advance_to_next_pose_segment(){
	// Read in pose.
	core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ );
	runtime_assert( rsd_set != 0 );
	import_pose_ = pose::PoseOP( new core::pose::Pose );
	pose_input_stream_->fill_pose( *import_pose_, *rsd_set );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::copy_next_pose_segment( pose::Pose & pose ){
	advance_to_next_pose_segment();
	apply_current_pose_segment( pose, *import_pose_, false /*check_sequence_matches*/);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::copy_next_pose_segment( pose::Pose & pose,
	pose::Pose & import_pose,
	bool const check_sequence_matches,
	bool const align_pose_to_import_pose ){
	using namespace conformation;
	advance_to_next_pose_segment();
	import_pose = *import_pose_;
	apply_current_pose_segment( pose, import_pose, check_sequence_matches, align_pose_to_import_pose );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::apply_current_pose_segment( pose::Pose & pose ){
	apply_current_pose_segment( pose, *import_pose_, false /*check sequence matches*/ );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::apply_current_pose_segment( pose::Pose & pose,
	pose::Pose & import_pose,
	bool const check_sequence_matches,
	bool const align_pose_to_import_pose ){

	if ( check_sequence_matches && !backbone_only_ ) check_sequence( pose, import_pose );

	// Dec 2010 -- why do we need this? copy is based on atom names not indexes.
	//cleanup_pose( import_pose );

	std::map< Size, Size > res_map;  //This is map from sub numbering to input_res numbering..

	for ( Size n = 1; n <= input_res_.size(); n++ ) {
		res_map[ full_to_sub_[ input_res_[n] ] ] = slice_res_[ n ];
		//   std::cout << n << ' ' << input_res_[ n ] << ' ' << full_to_sub_[ input_res_[ n ] ] << " " << slice_res_[ n ] << std::endl;

		// silly, disulfides
		//    if ( import_pose.residue_type(n).is_protein() && import_pose.residue_type(n).has_variant_type( DISULFIDE ) ){
		//     //add_variant_type_to_pose_residue( pose, DISULFIDE, full_to_sub[ input_res[ n ] ] );
		//     change_cys_state( n, "CYD", pose.conformation() );
		//    }

	}

	//Does this work for the "overlap residue" case?? If there is a overlap residue, then order of input_res will manner...Parin Jan 2, 2010.
	//Is it possible to copy only backbone torsions?? Parin Jan 2, 2010.

	//std::cout << "IMPORT_POSE DURING COPY_NEXT_SEGMENT: " <<  import_pose.annotated_sequence( true ) << std::endl;
	//  std::cout << "FOLD TREE OF POSE: " << pose.fold_tree() << std::endl;

	// Dec 2010 -- why do we copy_dofs for virtual?
	core::pose::copydofs::copy_dofs_match_atom_names( pose, import_pose, res_map, backbone_only_, false /*ignore_sidechain*/, true /*ignore_virtual*/ );

	//can do an alignment too. Note that this does *not* happen by default. Need to set align_pose_to_import_pose = true.
	if ( align_pose_to_import_pose ) {
		id::AtomID_Map< id::AtomID > atom_ID_map = align::create_aligment_id_map_legacy( pose, import_pose, res_map );
		scoring::superimpose_pose( pose, import_pose, atom_ID_map );
	}

}


////////////////////////////////////////////////////////////////////////////////////////
Size
InputStreamWithResidueInfo::compute_size(){
	core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ );
	runtime_assert( rsd_set != 0 );
	utility::vector1< core::pose::PoseOP > import_pose_list_ = pose_input_stream_->get_all_poses( *rsd_set );
	Size const size = import_pose_list_.size();
	pose_input_stream_->reset();
	return size;
}

////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::cleanup_pose( pose::Pose & import_pose ) const {

	using namespace chemical;

	protocols::farna::make_phosphate_nomenclature_matches_mini( import_pose );

	// No virtual anything!
	utility::vector1< core::chemical::VariantType > remove_variants;
	remove_variants.push_back( core::chemical::VIRTUAL_PHOSPHATE );
	remove_variants.push_back( core::chemical::VIRTUAL_O2PRIME_HYDROGEN );
	remove_variants.push_back( core::chemical::CUTPOINT_LOWER );
	remove_variants.push_back( core::chemical::CUTPOINT_UPPER );
	// Also remove VIRTUAL_RESIDUE_VARIANT variant?

	for ( Size n = 1; n <= import_pose.total_residue(); n++  ) {
		for ( Size i = 1; i <= remove_variants.size(); i++ ) {
			if ( import_pose.residue_type( n ).has_variant_type( remove_variants[ i ] ) ) {
				pose::remove_variant_type_from_pose_residue( import_pose, remove_variants[ i ] , n );
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::check_sequence( pose::Pose const & pose, pose::Pose const & import_pose ){

	std::cout <<  pose.annotated_sequence( true ) << std::endl;
	std::cout << import_pose.annotated_sequence( true ) << std::endl;

	bool match( true );
	for ( Size n = 1; n <= slice_res_.size(); n++ ) {
		if ( ( slice_res_[ n ] > import_pose.total_residue() ) ||
				( full_to_sub_[ input_res_[ n ] ] > pose.total_residue() ) ||
				(  import_pose.sequence()[ slice_res_[n] - 1 ]  !=
				pose.sequence()[ full_to_sub_[ input_res_[n] ] - 1 ] ) ) {
			std::cout << " N         " << n << std::endl;
			std::cout << " SLICE_RES " << slice_res_[ n ] << std::endl;
			std::cout << " INPUT_RES " << input_res_[ n ] << std::endl;
			std::cout << " FULLTOSUB " << full_to_sub_[ input_res_[ n ] ] << std::endl;
			std::cout << " IMPORT  TOTRES "<< import_pose.total_residue()  << std::endl;
			std::cout << " DESIRED TOTRES "<< pose.total_residue()  << std::endl;
			match = false;
			break;
		}
	}
	if ( !match ) utility_exit_with_message( "mismatch in sequence between input pose and desired sequence, given input_res " );

}

//////////////////////////////////////////////////////////////////////////
void
initialize_input_streams_with_residue_info( utility::vector1< InputStreamWithResidueInfoOP > & input_streams_with_residue_info,
	utility::vector1< std::string > const & pdb_tags,
	utility::vector1< std::string > const & silent_files_in,
	utility::vector1< Size > const & input_res,
	utility::vector1< Size > const & input_res2
){

	using namespace import_pose::pose_stream;

	utility::vector1< utility::vector1< Size > > input_res_vectors;
	input_res_vectors.push_back( input_res );
	input_res_vectors.push_back( input_res2 );

	utility::vector1< Size  > blank_vector;

	// This was a silly convention (I shouldn't have used it.)
	// First set up silent file input
	for ( Size i = 1; i <= silent_files_in.size(); i++ ) {
		utility::vector1< std::string > silent_files1, pdb_tags1;
		silent_files1.push_back( silent_files_in[ i ] );
		pdb_tags1.push_back( pdb_tags[ i ] );
		InputStreamWithResidueInfoOP input_stream( new InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP( new SilentFilePoseInputStream( silent_files1, pdb_tags1 ) ), input_res_vectors[ i ] ) );
		input_streams_with_residue_info.push_back( input_stream );
	}
	// Then PDBs.
	for ( Size i = silent_files_in.size()+1; i <= pdb_tags.size(); i++ ) {
		std::string pose_name = pdb_tags[ i ];
		std::size_t found=pose_name.find(".pdb");
		if ( found==std::string::npos ) {
			pose_name.append(".pdb");
		}
		utility::vector1< std::string > pose_names;
		pose_names.push_back( pose_name );
		InputStreamWithResidueInfoOP input_stream( new InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP( new PDBPoseInputStream( pose_names ) ), input_res_vectors[ i ] ) );
		input_streams_with_residue_info.push_back( input_stream );
	}

}

//////////////////////////////////////////////////////////////////////////
void
InputStreamWithResidueInfo::set_backbone_only(  bool const setting ){
	backbone_only_ = setting;
}


} //protein
} //modeler
} //stepwise
} //protocols
