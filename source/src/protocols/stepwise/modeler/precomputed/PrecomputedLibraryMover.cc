// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

#include <core/kinematics/FoldTree.hh>

#include <numeric/random/random.hh>

#include <basic/database/open.hh>

#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.precomputed.PrecomputedLibraryMover" );
using namespace core::io::silent;
using namespace core::pose::full_model_info;

///////////////////////////////////////////////////////////////////////////////////////
//
// In stepwise modeling, sometimes we build dinucleotides, trinucleotides, etc.
//  and it doesn't make sense to resample those 'from scratch' every time. We
//  already know which conformations lead to reasonable models. This mover
//  allows readin & application of those precomputed models from disk.
//
// Caveats:
//
//   -- Having these precomputed low energy models almost certainly breaks detailed
//       balance. Could in principle be fixed by penalizing these models by a poor
//       weight, or lowering the probability of 'from_scratch' moves.
//
//   -- Currently, poses are stored in SilentFileData format. Might be more compact
//       to use MiniPose format.
//
//   -- Currently, the library of precomputed poses are indexed by sequences. Might
//       make better sense to index by the res_list that they correspond to.
//
//   -- Could pretty easily extend this mover to handle precomputed models that are not
//       totally 'free', but include residues built off a fixed domain. For example,
//       noncanonical base pairs built on a fixed helix.  Would allow 'stepwise assembly'
//       enumerative runs to be reused in 'stepwise monte carlo' runs. Not supported yet.
//
//                  -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace modeler {
namespace precomputed {

//Constructor
PrecomputedLibraryMover::PrecomputedLibraryMover()
{
	initialize_from_directory( basic::database::full_name( "sampling/rna/precomputed/" ) );
}

//Destructor
PrecomputedLibraryMover::~PrecomputedLibraryMover()
{}

void
PrecomputedLibraryMover::initialize_from_directory( std::string const & dir_name ){
	utility::vector1< std::string > filenames;
	utility::file::list_dir( dir_name, filenames );
	for ( Size n = 1; n <= filenames.size(); n++ ) {
		std::string const & filename = filenames[ n ];
		if ( filename.size() < 5 ) continue;
		if ( filename.substr(  filename.size()-4, 4 ) != ".out" ) continue;
		std::string const full_filename = dir_name + "/" + filename ;
		TR.Debug << TR.Magenta << "Reading in file: " << full_filename << TR.Reset << std::endl;

		SilentFileDataOP silent_file_data( new SilentFileData );
		silent_file_data->set_verbose( false );
		silent_file_data->read_file( full_filename );
		TR.Debug << TR.Magenta << "Number of models: " << silent_file_data->size() << TR.Reset << std::endl;
		std::string const sequence = silent_file_data->begin()->one_letter_sequence();
		library_map_[ sequence ] = silent_file_data;
	}
}

///////////////////////////////////////////////////////////////////
// for now, only can apply precomputed library moves to
//  poses that are totally 'free' -- could be generalized easily.

//
// Appropriate generalization is in SubMotifLibrary!
//
bool
PrecomputedLibraryMover::has_precomputed_move( core::pose::Pose const & pose ) const {
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	utility::vector1< Size > const & sample_res = const_full_model_info( pose ).sample_res();
	for ( Size i = 1; i <= res_list.size(); i++ ) {
		if ( !sample_res.has_value( res_list[ i ] ) ) return false;
	}
	for ( Size n = 1; n < pose.size(); n++ ) {
		if ( pose.fold_tree().is_cutpoint( n ) ) return false;
	}
	return ( library_map_.find( pose.sequence() ) != library_map_.end() );
}

///////////////////////////////////////////////////////////////////////////////////////
void
PrecomputedLibraryMover::apply( core::pose::Pose &  ){
	utility_exit_with_message( "call the const apply() function, not this one." );
}
///////////////////////////////////////////////////////////////////////////////////////
void
PrecomputedLibraryMover::apply( core::pose::Pose & pose ) const {
	runtime_assert( has_precomputed_move( pose ) );
	SilentFileDataOP library = library_map_.find( pose.sequence() )->second;
	SilentStructOP silent_struct = numeric::random::rg().random_element( library->structure_list() );

	Pose pose_scratch;

	silent_struct->fill_pose( pose_scratch, pose.residue( 1 ).residue_type_set() );
	pose.conformation() = pose_scratch.conformation();

	//simple fold_tree --> generalized in SubMotifLibrary.
	for ( Size n = 1; n < pose.size(); n++ ) runtime_assert( !pose.fold_tree().is_cutpoint( n ) );
}


} //precomputed
} //modeler
} //stepwise
} //protocols
