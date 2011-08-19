 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
 // vi: set ts=2 noet:
 //  CVS information:
 //  $Revision: 1.1.2.1 $
 //  $Date: 2005/11/07 21:05:35 $
 //  $Author: rhiju $
 // (c) Copyright Rosetta Commons Member Institutions.
 // (c) This file is part of the Rosetta software suite and is made available under license.
 // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
 // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/rna/RNA_ChunkLibrary.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/conformation/Residue.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Numeric headers
#include <numeric/constants.hh>

// C++ headers
// AUTO-REMOVED #include <fstream>
#include <iostream>


static numeric::random::RandomGenerator RG(2380934);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.rna.rna_chunk_library" ) ;

namespace protocols{
namespace rna{

	using namespace core;
	using namespace ObjexxFCL;

	using core::Size;
	using core::Real;

	using core::pose::ResMap;

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
											ResMap const & res_map ) {
		mini_pose_list_ = mini_pose_list;
		res_map_ = res_map;
	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::PoseOP > const & pose_list,
											ResMap const & res_map ) {
		for ( Size n = 1; n <= pose_list.size(); n++ ) {
			mini_pose_list_.push_back( core::pose::MiniPoseOP( new core::pose::MiniPose( *(pose_list[n]) ) ) );
		}
		res_map_ = res_map;
	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::~ChunkSet()	{}


	///////////////////////////////////////////////////////////////////////
	void
	ChunkSet::insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index ) const{

		using namespace core::pose;

		core::pose::MiniPose const & scratch_pose ( *(mini_pose_list_[ chunk_pose_index ]) );

		//		TR << "SCRATCH_POSE " << scratch_pose.sequence() << ' ' << scratch_pose.fold_tree() << std::endl;
		copy_dofs( pose, scratch_pose, res_map_, false  /*copy_dofs_for_junction_residues*/ );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	RNA_ChunkLibrary::RNA_ChunkLibrary(){
		// currently nothing.
		chunk_coverage_ = 0.0;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	// constructor -- needs a list of silent files. Each silent file
	//  has solutions for a particular piece of the desired pose.
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & silent_files,
								std::string const & sequence_of_big_pose /* to figure out mapping to big pose */,
								std::map< Size, Size > const & connections_in_big_pose /* to figure out mapping to big pose*/ )
	{

		// allow_insert keeps track of where chunks are placed -- only allow
		// fragment insertions *outside* these regions.
		allow_insert_.dimension( sequence_of_big_pose.size(), true );
		covered_by_chunk_.dimension( sequence_of_big_pose.size(), false );

		for ( Size n = 1; n <= silent_files.size(); n++ ) {

			utility::vector1< pose::PoseOP > pose_list;
			process_silent_file( silent_files[n], pose_list );

			// There may be more than one part of the pose to which this sequence maps.
			utility::vector1< ResMap > res_maps;
			figure_out_possible_res_maps( res_maps, *(pose_list[1]), sequence_of_big_pose, connections_in_big_pose );

			for (Size k = 1; k <= res_maps.size(); k++ )  {
				check_res_map( res_maps[ k ], *(pose_list[1]), sequence_of_big_pose );

				ChunkSetOP chunk_set( new ChunkSet( pose_list, res_maps[ k ] ) );
				chunk_sets_.push_back( chunk_set );

				zero_out_allow_insert( res_maps[ k ], connections_in_big_pose );
			}

		}

		figure_out_chunk_coverage();

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::add_chunk_set(
												std::string const & silent_file,
												ResMap const & res_map,
												pose::Pose const & big_pose )
	{

		utility::vector1< pose::PoseOP > pose_list;

		process_silent_file( silent_file, pose_list );
		check_res_map( res_map, *(pose_list[1]), big_pose.sequence() );

		ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
		chunk_sets_.push_back( chunk_set );

	}

	////////////////////////////////////////////////////////////////////////////
	void RNA_ChunkLibrary::insert_chunk_into_pose(
					 pose::Pose & pose,
					 Size const & chunk_list_index,
					 Size const & chunk_pose_index ) const
	{
		chunk_sets_[ chunk_list_index ]->insert_chunk_into_pose( pose, chunk_pose_index );
	}


	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::random_chunk_insertion( core::pose::Pose & pose ) const{

		Size const chunk_set_index = static_cast <int> ( RG.uniform() * num_chunk_sets() ) + 1;

		ChunkSet const & chunk_set( *chunk_sets_[ chunk_set_index ] );
		Size const chunk_index = static_cast <int> ( RG.uniform() * chunk_set.num_chunks() ) + 1;

		chunk_set.insert_chunk_into_pose( pose, chunk_index );

		//		TR << "INSERTED CHUNK " << chunk_index << " FROM SET " << chunk_set_index << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::zero_out_allow_insert( ResMap const & res_map,
																			 std::map< Size, Size > const & connections_in_big_pose )
	{

		FArray1D< bool > connected( allow_insert_.size(), false );

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// THIS MASKS OUT STEMS AS STILL OK FOR FRAGMENT INSERTIONS ...
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for ( ResMap::const_iterator
						it=connections_in_big_pose.begin(), it_end = connections_in_big_pose.end(); it != it_end; ++it ) {
			Size const i = it->first; //Index in big pose.
			connected( i ) = true;
		}

		covered_by_chunk_ = false;

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
			Size const i = it->first; //Index in big pose.
			covered_by_chunk_( i ) = true;
			if ( !connected( i ) ) {
				allow_insert_( i ) = false;
			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::figure_out_chunk_coverage()
	{

		Size const tot_res( allow_insert_.size() );
		Size num_chunk_res( 0 );
		Size num_other_res( 0 );

		for ( Size n = 1; n <= tot_res; n++ ) {
			// Allow insert keeps track of where the chunk *aren't*, and
			// where other moves (fragments, jumps) can be carried out.
			if ( covered_by_chunk_(n) ){
				num_chunk_res++;
			} else {
				num_other_res++;
			}
		}
		chunk_coverage_ = Real( 3 * num_chunk_res ) / ( 3 * num_chunk_res +  tot_res );

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::get_component_sequences(
							 core::pose::Pose const & pose,
							 utility::vector1< std::string > & sequences,
							 utility::vector1< core::Size > & chain_id,
							 utility::vector1< core::Size > & sequence_start ) const{

		chain_id.clear();
		sequences.clear();

		std::string sequence = "";
		Size count( 1 );
		sequence_start.push_back( 1 );

		for ( Size i = 1; i <= pose.total_residue(); i++ ) {

			sequence += pose.residue(i).name1();
			chain_id.push_back( count );

			if ( pose.fold_tree().is_cutpoint( i ) ) {
				sequences.push_back( sequence );
				sequence = "";
				count++;
				if ( i < pose.total_residue() ) sequence_start.push_back( i+1 );
			}

		}


		for (Size n = 1; n <= sequences.size(); n++ ) {
			TR << "SEQUENCE " << n << " " << sequences[ n ] << std::endl;
		}

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::figure_out_possible_res_maps(
					utility::vector1< ResMap > & res_maps,
					pose::Pose const & scratch_pose,
					std::string const & sequence_of_big_pose,
					std::map< Size, Size > const & connections_in_big_pose ) const
	{

		// Note -- using zero-indexed vectors -- easier to do modulo, etc.
		utility::vector1< std::string > scratch_sequences;
		utility::vector1< core::Size >  chain_id;
		utility::vector1< core::Size >  scratch_sequence_start;

		get_component_sequences( scratch_pose, scratch_sequences, chain_id, scratch_sequence_start );

		// Go through each sequence and look for matches
		utility::vector1< utility::vector1< Size > > matches_to_each_scratch_sequence;
		get_sequence_matches( matches_to_each_scratch_sequence, scratch_sequences, sequence_of_big_pose );

		// now actually find some good res maps.
		find_res_maps( chain_id, scratch_sequence_start, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, res_maps );
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::find_res_maps(
								utility::vector1< Size > const & chain_id,
								utility::vector1< Size > const & scratch_sequence_start,
								utility::vector1< std::string > const & scratch_sequences,
								utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
								core::pose::Pose const & scratch_pose,
								std::map< Size, Size > const & connections_in_big_pose,
								utility::vector1< ResMap > & res_maps ) const
	{
		res_maps.clear();

		// Loop over matches for chain 1 -- if there are any chunks that match, they should be taggable by
		// their match to chain 1.

		Size num_chain( 1 );
		for (Size k = 1; k <= matches_to_each_scratch_sequence[ num_chain ].size(); k++ ) {
			ResMap res_map;
			fill_res_map( res_map, matches_to_each_scratch_sequence[ num_chain ][ k ],
										scratch_sequence_start[ num_chain ] /*should be 1*/,
										scratch_sequences[ num_chain ].size()  );

			// Dig deep into each connection from this chain.
			check_connections( num_chain, res_map,
												 chain_id, scratch_sequence_start, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, res_maps );
		}

		TR << "Number of matches found:  " << res_maps.size() << std::endl;
		if ( res_maps.size() == 0 )  utility_exit_with_message(  "Could not match silent file with sequence "+scratch_pose.sequence() );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::get_sequence_matches( 	utility::vector1< utility::vector1< Size > > & matches_to_each_scratch_sequence,
																						utility::vector1< std::string > const &	scratch_sequences,
																						std::string const & sequence_of_big_pose ) const
	{
		//		Size tot_matches = 1;
		for ( Size n = 1; n <= scratch_sequences.size(); n++ ) {
			utility::vector1< Size > matches;
			std::string const scratch_sequence( scratch_sequences[n] );
			Size const scratch_sequence_length = scratch_sequence.size();
			for ( Size i = 0; i <= sequence_of_big_pose.size() - scratch_sequence_length; i++ ) {
				bool does_it_match( true );
				for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
					if ( sequence_of_big_pose[ i + offset ] != scratch_sequence[ offset ] ) {
						does_it_match = false;
						break;
					}
				}
				if (does_it_match) {
					matches.push_back( i+1 ); // convert numbering to start with 1
					TR << "Found match to scratch_sequence " << n <<
						//						"  starting at scratch pose position " << scratch_sequence_start[ n ] << " " <<
						"   at big pose position: " << i+1 << std::endl;
				}
			}

			matches_to_each_scratch_sequence.push_back( matches );

			if ( matches.size() < 1 ) 	 utility_exit_with_message(  "Could not find match to sequence" );

		}
}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::check_connections( Size const & num_chain, ResMap & res_map,
																			 utility::vector1< Size > const & chain_id,
																			 utility::vector1< Size > const & scratch_sequence_start,
																			 utility::vector1< std::string > const & scratch_sequences,
																			 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																			 core::pose::Pose const & scratch_pose,
																			 std::map< Size, Size > const & connections_in_big_pose,
																			 utility::vector1< ResMap > & res_maps ) const
	{
		// We might already be done!
		if ( res_map.size() == scratch_pose.total_residue() ) { //everything assigned!
			res_maps.push_back( res_map );
			return;
		}

		// Look for connections coming off the current chain
		Size const n_jump( scratch_pose.num_jump() ) ;

		for ( Size n = 1; n <= n_jump; n++ ) {

			//			TR << "CHECK OUT JUMP  " << n << " for chain " << num_chain << std::endl;
			Size const res1 =  scratch_pose.fold_tree().upstream_jump_residue( n );
			Size const res2 =  scratch_pose.fold_tree().downstream_jump_residue( n );

			if ( chain_id[ res1 ] == num_chain ) {
				test_matches( res1, res2, res_map,
											chain_id, scratch_sequence_start,
											scratch_sequences, matches_to_each_scratch_sequence,
											scratch_pose, connections_in_big_pose, res_maps );
			} else if ( chain_id[ res2 ] == num_chain ) {
				test_matches( res2, res1, res_map,
											chain_id, scratch_sequence_start,
											scratch_sequences, matches_to_each_scratch_sequence,
											scratch_pose, connections_in_big_pose, res_maps );
			}

		}

	}

	//////////////////////////////////////////////
	void
	RNA_ChunkLibrary::test_matches( Size const & res1, Size const & res2, ResMap & res_map,
																	utility::vector1< Size > const & chain_id,
																	utility::vector1< Size > const & scratch_sequence_start,
																	utility::vector1< std::string > const & scratch_sequences,
																	utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																	core::pose::Pose const & scratch_pose,
																	std::map< Size, Size > const & connections_in_big_pose,
																	utility::vector1< ResMap > & res_maps ) const {

		//		TR << " TESTING : " << res1 << " " << res2 << "  from chains: " << chain_id[ res1 ] << " " << chain_id[ res2 ] <<
		//			". Length of resmap " << res_map.size() << " out of " << scratch_pose.total_residue() <<  std::endl;

		Size const next_chain = chain_id[ res2 ];
		// res2 may already be "taken care of" -- no need to check it out.
		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
			if ( it->second == res2 ) {
				//				TR << " ALREADY TESTED " << std::endl;
				return;
			}
		}

		// Cycle through potential matches -- are any good?
		for (Size k = 1; k <= matches_to_each_scratch_sequence[ next_chain ].size(); k++ ) {

			ResMap res_map_test( res_map );

			bool const res_map_ok = fill_res_map( res_map_test, matches_to_each_scratch_sequence[ next_chain ][ k ],
																						scratch_sequence_start[ next_chain ],scratch_sequences[next_chain].size()  );

			if (!res_map_ok) continue;

			// Are the res1 and res2 connected in the big pose?
			// Need to figure out where they map to in big pose.
			Size res1_map( 0 ), res2_map( 0 );
			for ( ResMap::const_iterator
							it=res_map_test.begin(), it_end = res_map_test.end(); it != it_end; ++it ) {
				if ( it->second == res1 ) res1_map = it->first;
				if ( it->second == res2 ) res2_map = it->first;
			}
			if (res1_map == 0 || res2_map == 0 ){
				TR << res1 << " " << res1_map << "        " << res2 << " " << res2_map << std::endl;
				utility_exit_with_message( "SHOULD NOT BE HERE! " );
			}


			bool connection_ok( false );
			for ( ResMap::const_iterator
							it=connections_in_big_pose.begin(), it_end = connections_in_big_pose.end(); it != it_end; ++it ) {

				Size const res1_in_big_pose = it->first;
				Size const res2_in_big_pose = it->second;

				if ( res1_in_big_pose == res1_map &&
						 res2_in_big_pose == res2_map ) {
					connection_ok = true; break;
				}

				if ( res2_in_big_pose == res1_map &&
						 res1_in_big_pose == res2_map ) {
					connection_ok = true; break;
				}
			}

			if (!connection_ok) {
				//				TR << "DENIED!!  " << std::endl;
				continue;
			} else {
				//				TR << "OK!!  " << std::endl;
				res_map = res_map_test;
				check_connections( next_chain, res_map,
													 chain_id, scratch_sequence_start,
													 scratch_sequences, matches_to_each_scratch_sequence,
													 scratch_pose, connections_in_big_pose, res_maps );
				break;
			}
		}


	}


	///////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::fill_res_map( ResMap & res_map, Size const & match_pos, Size const & scratch_start_pos, Size const & scratch_sequence_length ) const
	{
		bool one_to_one( true );
		for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
			Size const big_pose_pos = match_pos + offset;
			Size const scratch_pos = scratch_start_pos + offset;
			if ( res_map.find( big_pose_pos) == res_map.end() ) {
				res_map[ big_pose_pos ] = scratch_pos;
				TR << "MAPPING " << match_pos+offset << " --> " << scratch_start_pos+offset << std::endl;
			} else {
				// this is a problem -- expect res_maps to be one-to-one.
				one_to_one = false;
				break;
			}
		}
		return one_to_one;
	}



	///////////////////
	// DELETE FOLLOWING AFTER SHIT WORKS!
	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::check_res_map_recursively( ResMap const & res_map_old,
																							 utility::vector1< std::string > const & scratch_sequences,
																							 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																							 pose::Pose const & scratch_pose,
																							 std::map< Size, Size > const & connections_in_big_pose,
																							 utility::vector1< core::Size > const & chain_id,
																							 Size const & num_sequence, Size const & num_match, utility::vector1< ResMap > & res_maps ) const{
		ResMap res_map( res_map_old );

		Size const match_pos = matches_to_each_scratch_sequence[ num_sequence ][ num_match ];
		//		TR << "HEY MATCH_POS " << match_pos << std::endl;
		Size const scratch_sequence_length = scratch_sequences[num_sequence].size();
		Size i( res_map.size() );
		for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
			res_map[ match_pos + offset + 1 ] = i + 1; // Add back in 1 to match pose numbering.
			i++;
		}

		bool const jump_match = check_jump_match( scratch_pose, connections_in_big_pose, res_map, chain_id );

		//		for ( ResMap::const_iterator
		//						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
		//			TR << it->first << " mapped to " << it->second << std::endl;
		//		}

		//		TR << " CHECKING: " << num_sequence << " -- " << num_match << " : " << res_map.size() << " --> " << jump_match << std::endl;

		if ( !jump_match ) return;

		if ( num_sequence == matches_to_each_scratch_sequence.size() ) { //Success!
			res_maps.push_back( res_map );
		}

		//		TR << "ADVANCING to sequence " << num_sequence << std::endl;
		Size const num_sequence_next = num_sequence + 1;

		if ( num_sequence_next <= matches_to_each_scratch_sequence.size() ) {

			for ( Size k = 1; k <= matches_to_each_scratch_sequence[ num_sequence_next ].size(); k++ ) {
				//				TR << "ABOUT TO TRY  " << num_sequence_next << " " << k << std::endl;
				check_res_map_recursively( res_map, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, chain_id, num_sequence_next, k, res_maps );
			}

		}

		return;

	}


	///////////////
	///////////////
	///////////
	// DELETE THIS AFTER SHIT IS WORKING.
// 		// Go through each sequence and look for matches
// 		utility::vector1< std::vector< Size > > matches_to_each_scratch_sequence;
// 		Size tot_matches = 1;
// 		for ( Size n = 1; n <= scratch_sequences.size(); n++ ) {
// 			std::vector< Size > matches;
// 			std::string const scratch_sequence( scratch_sequences[n] );
// 			Size const scratch_sequence_length = scratch_sequence.size();
// 			for ( Size i = 0; i <= sequence_of_big_pose.size() - scratch_sequence_length; i++ ) {
// 				bool does_it_match( true );
// 				for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
// 					if ( sequence_of_big_pose[ i + offset ] != scratch_sequence[ offset ] ) {
// 						does_it_match = false;
// 						break;
// 					}
// 				}
// 				if (does_it_match) {
// 					matches.push_back( i );
// 					TR << "Found match to scratch_sequence " << n << "   at position: " << i+1 << std::endl;
// 				}
// 			}

// 			matches_to_each_scratch_sequence.push_back( matches );
// 			num_sequence_match.push_back( matches.size() );

// 			if ( matches.size() < 1 ) 	 utility_exit_with_message(  "Could not find match to sequence" );

// 			tot_matches *= num_sequence_match[ n ];

// 		}


// 		// Test out every possible combination of sequence match -- make sure that
// 		// jump-connected residues inside the mini "scratch pose" are connected in the
// 		// big pose.
// 		// Note that following is not optimal at all, but this code doesn't have to be --
// 		// just runs once.
// 		bool found_a_match( false );
// 		TR << "TOT_MATCHES " << tot_matches << std::endl;
// 		for ( Size k = 0; k < tot_matches; k++ ) {

// 			Size count = k;

// 			////////////////////////////////////////////////////////////////////////////////////////////
// 			////////////////////////////////////////////////////////////////////////////////////////////
// 			// This is a pretty ghetto decoding mechanism. The other option might be recursion.
// 			utility::vector1< Size > which_match;
// 			for ( Size n = 1; n <= scratch_sequences.size(); n++ ) {
// 				which_match.push_back(  count % num_sequence_match[ n ] );
// 				count = floor( count / num_sequence_match[ n ] );
// 			}

// 			// Fill out res_map
// 			ResMap res_map;
// 			 Size i( 0 );
// 			for ( Size n = 1; n <= scratch_sequences.size(); n++ ) {
// 				Size const match_pos = matches_to_each_scratch_sequence[ n ][ which_match[n] ];
// 				Size const scratch_sequence_length = scratch_sequences[n].size();
// 				for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
// 					res_map[ match_pos + offset + 1 ] = i + 1; // Add back in 1 to match pose numbering.
// 					i++;
// 				}
// 			}

// 			if ( check_jump_match( scratch_pose, connections_in_big_pose, res_map, chain_id ) ) {
// 				res_maps.push_back( res_map );
// 				found_a_match = true;
// 				for ( ResMap::const_iterator
// 								it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
// 					TR << it->first << " mapped to " << it->second << std::endl;
// 				}

// 			}
// 		}
// 		if ( !found_a_match )  utility_exit_with_message(  "Could not match silent file with sequence "+scratch_pose.sequence() );

// 	}


	/////////////
	///////////
	/////////////


	///////////////////////////////////////////////////////////////////////////
	// Any jump inside the scratch pose (which came from a user-inputted silent file with its own fold tree)
	// must correspond to a pairing defined in the new pose of interest. Could make this a little less restrictive,
	// just asking for any connection between the chain segments that are putatively matching.
	bool
	RNA_ChunkLibrary::check_jump_match(
									 pose::Pose const & scratch_pose,
									 std::map< Size, Size > const & connections_in_big_pose,
									 ResMap const & res_map,
									 utility::vector1< Size > const & chain_id ) const
	{

		Size const n_jump( scratch_pose.num_jump() ) ;

		for ( Size n = 1; n <= n_jump; n++ ) {

			Size const res1 =  scratch_pose.fold_tree().upstream_jump_residue( n );
			Size const res2 =  scratch_pose.fold_tree().downstream_jump_residue( n );

			// First check that res1 and res2 are in the resmap -- otherwise no point in looking for a matching
			// jump in the big pose.
			bool found_res1( false ), found_res2( false );
			for ( ResMap::const_iterator it = res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
				if ( res1 == it->second ) found_res1 = true;
				if ( res2 == it->second ) found_res2 = true;
			}
			if ( !found_res1 || !found_res2 ) {
				//consider it OK if the res_map isn't complete yet.
				continue;
			}

			bool connection_ok( false );

			for ( ResMap::const_iterator
							it=connections_in_big_pose.begin(), it_end = connections_in_big_pose.end(); it != it_end; ++it ) {

				Size const res1_in_big_pose = it->first;
				Size const res2_in_big_pose = it->second;

				//				TR << "Checking connection: " << res1_in_big_pose << " - " << res2_in_big_pose << std::endl;

				ResMap::const_iterator res1_map_id = res_map.find( res1_in_big_pose );
				ResMap::const_iterator res2_map_id = res_map.find( res2_in_big_pose );

				if ( res1_map_id != res_map.end() &&
						 res2_map_id != res_map.end() ) {

					Size const res1_in_scratch_pose = res1_map_id->second;
					Size const res2_in_scratch_pose = res2_map_id->second;
					//					TR << " Checking connection in scratch_pose: " << res1_in_scratch_pose << " - " << res2_in_scratch_pose << std::endl;

					if ( ( chain_id[ res1 ] == chain_id[ res1_in_scratch_pose ] &&
								 chain_id[ res2 ] == chain_id[ res2_in_scratch_pose ] ) ||
							 ( chain_id[ res1 ] == chain_id[ res2_in_scratch_pose ] &&
								 chain_id[ res2 ] == chain_id[ res1_in_scratch_pose ] ) ) {
						connection_ok = true;
						break;
					}

				}
			}

			if (!connection_ok) {
				//				TR << "FAIL! Could not find a match to " << res1 << " -- " << res2 << std::endl;
				return false;
			}

		}

		return true;
	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::check_res_map( ResMap const & res_map, pose::Pose const & scratch_pose, std::string const & sequence ) const{

		// CHECK SEQUENCE HERE!!! EXIT IF NO MATCH!!!!!!!

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

			// For now, just do bonded atoms...update later to do jumps too.
			Size const i = it->first; //Index in big pose.
			Size const i_scratch_pose = it->second; // Index in the little "chunk" or "scratch" pose

			if ( sequence[ i-1 ] != scratch_pose.residue( i_scratch_pose ).name1() ){
				utility_exit_with_message(  "Mismatched sequence!!" );
				return false;
			}

		}
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::process_silent_file( std::string const & silent_file,
																		 utility::vector1< pose::PoseOP > & pose_list ) const
	{
		using namespace  core::io::silent;

		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

		SilentFileData silent_file_data;
		silent_file_data.read_file( silent_file );

		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
							end = silent_file_data.end(); iter != end; ++iter ) {

			pose::PoseOP pose_op( new pose::Pose );

			std::string const tag = iter->decoy_tag();
			iter->fill_pose( *pose_op );

			pose_list.push_back( pose_op );
		}

		if ( pose_list.size() < 1)  {
			utility_exit_with_message(  "No structure found in silent file  " + silent_file );
		}

	}

	////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::initialize_random_chunks( pose::Pose & pose ) const{
		for ( Size n = 1; n <= num_chunk_sets(); n++ ) {

			ChunkSet const & chunk_set( *chunk_sets_[ n ] );
			Size const chunk_index = static_cast <int> ( RG.uniform() * chunk_set.num_chunks() ) + 1;
			//			TR << "NUM_CHUNKS " << chunk_index << " " << chunk_set.num_chunks() << std::endl;
			chunk_set.insert_chunk_into_pose( pose, chunk_index );

			//			pose.dump_pdb( "start_"+string_of(n)+".pdb" );
		}

	}

}
}

