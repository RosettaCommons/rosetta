// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/rna/denovo/libraries/RNA_JumpLibrary.hh>
#include <protocols/rna/denovo/libraries/RNA_LibraryManager.hh>
#include <protocols/rna/denovo/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.monte_carlo.submotif.SubMotifLibrary" );

///////////////////////////////////////////////////////////////////////////////////
//
// Allow testing of U-turns, UA_handle's, K-turns, etc. that sit in a library in
//  the Rosetta database.
//
// Some overlap with PrecomputedLibraryMover -- may deprecate that one.
//
// This may also offer a viable alternative to Fragment and ChunkLibrary classes
//   in FARFAR (rna denovo).
//
// Note that this is was moved to stepwise::monte_carlo rather than stepwise::modeler
//   namespace as it is only used in StepWiseMonteCarlo. Conceptually SubMotifLibrary
//   offers a bunch of fixed poses for used in monte carlo.
//
//   -- rhiju, 2015
//
///////////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;
using namespace core::chemical::rna;
using namespace protocols::stepwise::monte_carlo::mover;
using core::kinematics::FoldTree;
using utility::vector1;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace submotif {

//Constructor
SubMotifLibrary::SubMotifLibrary(
	core::chemical::ResidueTypeSetCAP rsd_set,
	bool const include_submotifs_from_jump_library,
	bool const use_first_jump_for_submotif,
	utility::vector1< std::string > const & exclude_submotif_list ):
	rsd_set_( rsd_set ),
	include_submotifs_from_jump_library_( include_submotifs_from_jump_library ),
	use_first_jump_for_submotif_( use_first_jump_for_submotif ),
	exclude_submotif_list_( exclude_submotif_list )
{
	initialize();
}

//Destructor
SubMotifLibrary::~SubMotifLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////////////////
// pose is only used to figure out if we're doing an RNA only run, which can
// save us some time in instantiating the residue type set. Silly -- we need
// to refactor ResidueTypeSet.
void
SubMotifLibrary::initialize(){
	// could also create a protein library for, e.g., helices & beta strand pairings.
	initialize_from_directory( basic::database::full_name( "sampling/rna/submotif/" ) );
	if ( include_submotifs_from_jump_library_ ) initialize_from_jump_library();
}

//////////////////////////////////////////////////////////////////////////////////////////
void
SubMotifLibrary::initialize_from_directory( std::string const & directory ){
	utility::vector1< std::string > filenames;

	utility::io::izstream data( ( directory+"/submotifs.txt" ).c_str() );
	if ( data.good() ) { // add extra data
		std::string line;
		while ( getline( data, line ) ) {
			if ( line.size() && line[0] == '#' ) continue;
			filenames.push_back( line );
		}
	}

	for ( auto const & exclude_submotif : exclude_submotif_list_ ) {
		if ( !filenames.has_value( exclude_submotif ) ) {
			TR << TR.Red << exclude_submotif << " must be in " << filenames;
			utility_exit_with_message( "requested exclude submotif is not in list" );
		}
	}

	for ( Size n = 1; n <= filenames.size(); n++ ) {
		PoseTag const & tag = filenames[ n ];
		if ( exclude_submotif_list_.has_value( tag ) ) continue;
		if ( tag.substr(  tag.size()-4, 4 ) != ".pdb" ) continue;
		std::string const full_filename = directory + "/" + tag ;
		runtime_assert( utility::file::file_exists( full_filename ) );

		TR  << "Reading in submotif: " << full_filename << TR.Reset << std::endl;
		PoseOP pose = core::import_pose::get_pdb_and_cleanup( full_filename, rsd_set_ );
		save_pose_as_submotif( pose, tag );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
void
SubMotifLibrary::save_pose_as_submotif( pose::PoseOP pose, std::string const & tag )
{
	tag_into_pose( *pose, tag );
	submotif_poses_by_tag_[ tag ] = pose;
	SubMotifSequenceSet submotif_sequence_set = get_submotif_sequence_set( *pose );
	submotif_sequence_sets_.insert( submotif_sequence_set );
	utility::vector1< SequenceMapping > matches = get_matches_for_one_submotif_sequence_set( submotif_sequence_set, *pose, false /*use_full_model_info*/ );
	runtime_assert( matches.size() > 0 );
	for ( Size n = 1; n <= matches.size(); n++ ) {
		if ( n > 1 ) continue; // hmm. trying to avoid redundancy -- all threadings will be enumerated in target full_model in get_matches_for_one_submotif_sequence_set later.
		submotif_mappings_by_sequence_set_[ submotif_sequence_set ].push_back( std::make_pair( tag, matches[ n ] ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
void
SubMotifLibrary::initialize_from_jump_library()
{
	using namespace protocols::rna::denovo;
	using namespace protocols::rna::denovo::libraries;
	RNA_JumpLibraryCOP rna_jump_library( RNA_LibraryManager::get_instance()->rna_jump_library_cop() );
	for ( Size i = 0; i <= 3; i++ ) {
		for ( Size j = 0; j <= 3; j++ ) {
			for ( Size e1 = 1; e1 <= 3; e1++ ) {
				for ( Size e2 = 1; e2 <= 3; e2++ ) {
					for ( Size o = 1; o <= 2; o++ ) {
						BasePairType const base_pair_type( rna_nts[ i ], rna_nts[ j ], BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) );
						if ( !rna_jump_library->has_template( base_pair_type ) ) continue;
						TR << "filling out submotifs for: " << base_pair_type.tag() << std::endl;
						RNA_PairingTemplateList const & templates( rna_jump_library->rna_pairing_template_map( base_pair_type ) );
						Size const num_templates = use_first_jump_for_submotif_ ? 1 : templates.size();
						for ( Size n = 1; n <= num_templates; n++ ) {

							PoseOP pose( new Pose );
							std::string sequence;
							sequence += rna_nts[ i ];
							sequence += rna_nts[ j ];
							core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set_ );
							make_pose_from_sequence( *pose, sequence, *rsd_set_op, false /*auto termin*/  );

							FoldTree f( 2 );
							f.new_jump( 1, 2, 1 );
							fill_in_default_jump_atoms( f, *pose );
							pose->fold_tree( f );
							core::kinematics::Jump const & j( templates[ n ]->jump_forward() );
							pose->set_jump( 1, j );
							core::import_pose::cleanup( *pose, true /* force_cut_at_rna_chainbreak */ );

							std::string tag( base_pair_type.tag() + "_" + ObjexxFCL::lead_zero_string_of( n, 5 ) );
							//pose->dump_pdb( tag + ".pdb" );
							save_pose_as_submotif( pose, tag );

						} // n
					} // o
				} // e2
			} // e1
		} // j
	} // i
}

//////////////////////////////////////////////////////////////////////////////////////////
SubMotifSequenceSet
SubMotifLibrary::get_submotif_sequence_set( pose::Pose const & pose, bool sort_sequences /* = true */ ) const {

	utility::vector1< std::string > submotif_sequence_set;
	std::string chain_sequence( "" );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		chain_sequence += pose.sequence()[ n - 1 ];
		if ( pose.fold_tree().is_cutpoint( n ) ) {
			submotif_sequence_set.push_back( chain_sequence );
			chain_sequence = "";
		}
	}
	if ( chain_sequence.size() > 0 ) submotif_sequence_set.push_back( chain_sequence );

	// important -- need to match up submotifs if they have the same chain sequences.
	if ( sort_sequences ) std::sort( submotif_sequence_set.begin(), submotif_sequence_set.end() );

	return submotif_sequence_set;
}

//////////////////////////////////////////////////////////////////////////////////////////
// currently slow but a complete recursive search.
// Could be computed once, and then filtered by pose_domain_map...
//////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< StepWiseMove >
SubMotifLibrary::get_submotif_moves( pose::Pose const & pose ) const {

	utility::vector1< StepWiseMove > submotif_moves;

	for ( std::set< SubMotifSequenceSet >::const_iterator it = submotif_sequence_sets_.begin(),
			end = submotif_sequence_sets_.end(); it != end; ++it ) {

		SubMotifSequenceSet const & submotif_sequence_set( *it );
		utility::vector1< SequenceMapping > all_matches =
			get_matches_for_one_submotif_sequence_set( submotif_sequence_set, pose );

		for ( Size m = 1; m <= all_matches.size(); m++ ) {
			SequenceMapping const & mapping_to_pose = all_matches[ m ];

			// need to combine this mapping from submotif_sequence_set --> pose
			// with stored          mapping from submotif_sequence_set --> submotif_pose
			utility::vector1< std::pair< PoseTag, SequenceMapping > > const & submotif_tags_with_mapping = submotif_mappings_by_sequence_set_.find( submotif_sequence_set )->second;
			for ( Size q = 1; q <= submotif_tags_with_mapping.size(); q++ ) {

				std::string const & submotif_tag = submotif_tags_with_mapping[ q ].first;
				SequenceMapping const & mapping_to_submotif_pose = submotif_tags_with_mapping[ q ].second;
				utility::vector1< Size > moving_res;
				for ( Size k = 1; k <= mapping_to_submotif_pose.size(); k++ ) {
					moving_res.push_back( mapping_to_pose[ mapping_to_submotif_pose.index( k ) ] );
				}

				StepWiseMove submotif_move( moving_res, Attachments(), FROM_SCRATCH, submotif_tag );
				submotif_moves.push_back( submotif_move );
			}
		}
	}

	return submotif_moves;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< SequenceMapping >
SubMotifLibrary::get_matches_for_one_submotif_sequence_set( SubMotifSequenceSet const & submotif_sequence_set,
	pose::Pose const & pose,
	bool const use_full_model_info /* = true */ ) const {

	// a little precomputation...

	Size submotif_length = 0;
	vector1< Size > submotif_cutpoints;
	std::string submotif_full_sequence;
	for ( Size n = 1; n <= submotif_sequence_set.size(); n++ ) {
		std::string const & submotif_sequence = submotif_sequence_set[ n ];
		submotif_length += submotif_sequence.size();
		submotif_cutpoints.push_back( submotif_length );
		submotif_full_sequence += submotif_sequence;
	}

	std::string pose_full_sequence = pose.sequence();
	utility::vector1< Size > pose_cutpoints;
	for ( Size n = 1; n < pose.size(); n++ ) {
		if ( pose.fold_tree().is_cutpoint( n ) ) pose_cutpoints.push_back( n );
	}
	utility::vector1< Size > pose_domain_map( pose_full_sequence.size(), 0 );

	if ( use_full_model_info ) {
		runtime_assert( full_model_info_defined( pose ) );
		pose_full_sequence   = const_full_model_info( pose ).full_sequence();
		pose_cutpoints  = const_full_model_info( pose ).cutpoint_open_in_full_model();
		pose_domain_map = figure_out_pose_domain_map_const( pose );
		utility::vector1< Size > const & sample_res( const_full_model_info( pose ).sample_res() );
		utility::vector1< Size > const & bulge_res( const_full_model_info( pose ).rna_bulge_res() );
		for ( Size n = 1; n <= pose_domain_map.size(); n++ ) {
			if ( ( !sample_res.has_value( n ) || bulge_res.has_value( n ) )  &&
					pose_domain_map[ n ] == 0 ) pose_domain_map[ n ] = 999; // disallow match unless sample_res.
		}
	}

	// do via a recursive search.
	utility::vector1< utility::vector1< Size > > all_matches;
	SequenceMapping matching_residues;
	get_matches( all_matches,
		matching_residues,
		submotif_full_sequence,
		submotif_cutpoints,
		pose_full_sequence,
		pose_cutpoints,
		pose_domain_map );

	return all_matches;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
SubMotifLibrary::get_matches( utility::vector1< SequenceMapping > & all_matches /* stores matches */,
	SequenceMapping const & matching_residues /* working mapping */,
	std::string const & submotif_full_sequence,
	utility::vector1< Size > const & submotif_cutpoints,
	std::string const & pose_full_sequence,
	utility::vector1< Size > const & pose_cutpoints,
	utility::vector1< Size > const & pose_domain_map ) const
{

	if ( matching_residues.size() == submotif_full_sequence.size() ) { // done!
		all_matches.push_back( matching_residues );
		return;
	}

	// if not done, look for next match.
	utility::vector1< Size > possible_next_res;
	Size const i_submotif = matching_residues.size() + 1;
	Size const nfull = pose_full_sequence.size();
	if ( i_submotif == 1 || submotif_cutpoints.has_value( i_submotif - 1 ) ) {
		// starting a new chain -- only constraint is that it should *not* be the very next residue.
		for ( Size i_full = 1; i_full <= nfull; i_full++ ) {
			if ( i_submotif == 1 || i_full == 1 || pose_cutpoints.has_value( i_full - 1 ) ||
					( i_full != matching_residues[ i_submotif - 1 ] && i_full != matching_residues[ i_submotif - 1 ] + 1 ) ) {
				possible_next_res.push_back( i_full );
			}
		}
	} else {
		// must be contiguous
		runtime_assert( matching_residues.size() >= 1 );
		Size const i_full = matching_residues[ i_submotif - 1 ] + 1;
		if ( !pose_cutpoints.has_value( i_full - 1 ) ) {
			possible_next_res.push_back( i_full );
		}
	}

	////////////////////////////////////////////////////////////////////////////
	//  for full_model residue to match submotif residue, should
	//
	//    1. actually be in range
	//    2. match in sequence
	//    3. be sample-able
	//    4. not already sampled in full_model.
	//    5. not already matched to another residue in the submotif.
	//
	// info on what is already sampled in full_model is in pose_domain_map_ which
	///  is set separately.
	//
	////////////////////////////////////////////////////////////////////////////
	for ( Size n = 1; n <= possible_next_res.size(); n++ ) {
		Size const i_full = possible_next_res[ n ];

		if ( i_full <= pose_full_sequence.size() &&
				submotif_full_sequence[ i_submotif - 1 ] == pose_full_sequence[ i_full - 1 ] &&
				pose_domain_map[ i_full ] == 0 &&
				!matching_residues.has_value( i_full ) &&
				( !matching_residues.has_value( i_full + 1 ) ||
				pose_cutpoints.has_value( i_full ) ) ) {

			SequenceMapping matching_residues_new( matching_residues );
			matching_residues_new.push_back( i_full );
			get_matches( all_matches,
				matching_residues_new,
				submotif_full_sequence,
				submotif_cutpoints,
				pose_full_sequence,
				pose_cutpoints,
				pose_domain_map );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
void
SubMotifLibrary::output_tags() const {
	for ( std::map< PoseTag, core::pose::PoseCOP >::const_iterator it = submotif_poses_by_tag_.begin(),
			end = submotif_poses_by_tag_.end(); it != end; ++it ) {
		TR << TR.Green << it->first << " " << get_submotif_sequence_set( *(it->second), false /*sort_sequences*/ ) << TR.Reset << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
pose::PoseOP
SubMotifLibrary::create_new_submotif( SequenceMapping const & move_element,
	PoseTag const & submotif_tag,
	pose::Pose const & pose,
	bool const & seed /*= false*/ ) const {

	TR.Debug << TR.Magenta << "Creating: " << submotif_tag << TR.Reset << std::endl;

	if ( submotif_poses_by_tag_.find( submotif_tag ) == submotif_poses_by_tag_.end() ) {
		TR << std::endl << TR.Green << "Valid submotif tags " << TR.Reset << std::endl;
		output_tags();
		utility_exit_with_message( "Invalid SUBMOTIF tag: "+submotif_tag );
	}

	PoseOP new_pose = submotif_poses_by_tag_.find( submotif_tag )->second->clone();

	std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
	for ( Size n = 1; n <= move_element.size(); n++ ) {
		runtime_assert( full_sequence[ move_element[ n ]-1 ] == new_pose->sequence()[ n - 1 ] );
		// should also check that cutpoints match up.
	}

	FullModelInfoOP full_model_info_for_pose( new FullModelInfo( const_full_model_info( pose ).full_model_parameters() )  );
	full_model_info_for_pose->set_res_list( move_element );
	full_model_info_for_pose->add_submotif_info( move_element, submotif_tag, seed );
	full_model_info_for_pose->update_submotif_info_list();
	set_full_model_info( *new_pose, full_model_info_for_pose );

	core::pose::fix_up_residue_type_variants( *new_pose );

	return new_pose;
}

} //submotif
} //monte_carlo
} //stepwise
} //protocols
