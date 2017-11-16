// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/canonical_sampling/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/MPIHPool_ConvergenceCheck.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/pool_util.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.fwd.hh>


// MPI only headers (wrapping these headers in an ifdef block will prevent my #inclusion-removal script from removing them)
#ifdef USEMPI
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/canonical_sampling/mc_convergence_checks/HPool.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/toolbox/superimpose.hh>

#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/prof.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#endif

#include <core/io/silent/SilentStruct.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>


#include <utility/vector1.hh>


static basic::Tracer tr( "MPIHPool_ConvergenceCheck" );

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

#ifdef USEMPI

  int const FINISHED = 1;
  int const IN_PROGRESS = 0;
  int const MPI_OUTPUT_RANK = 0;
  int const OUTPUT_TAG = 1000;

  MPI_Comm   protocols::canonical_sampling::mc_convergence_checks::MPIHPool_RMSD::MPI_COMM_POOL;
  using namespace basic;


  //specified silent-file assumed to contain top-level structures only
  MPIHPool_RMSD::MPIHPool_RMSD( std::string const& silent_file, core::Size levels ):
    Pool_RMSD(),
    hlevel_( levels, silent_file ),
    pool_size_(0),
    num_structures_added_(0),
    npes_(-1),
    rank_(-1),
    new_decoys_out_("discovered_decoys.out"),
    nresidues_(-1),
    nlevels_( levels ),
    first_time_writing_( true ),
    level_radii_(),
    current_address_(),
    current_best_rmsds_(),
    best_address_(),
    current_address_str_("no_address"),
    buf_(),
    current_trajectory_state_(IN_PROGRESS)
  {
    initialize();
  }

  void
  MPIHPool_RMSD::receive_and_output_structures(
					       core::io::silent::SilentFileData& sfd,
					       core::Size num_structures_to_write
					       ) {
		core::io::silent::SilentFileOptions opts;
    while( num_structures_to_write > 0 ) {
      core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
      receive_silent_struct_any_source( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ );
      if( first_time_writing_ ) {
	first_time_writing_ = false;
	write_headers_to_hierarchy( ss );
      }
      write_decoys_to_hierarchy( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ );
      num_structures_to_write--;
    }
  }

  void
  MPIHPool_RMSD::resolve_address_and_assign_tag( Address& new_addr, core::Size& new_level_start, std::string& new_candidate_tag ) {
    new_level_start = 0;
    for( core::Size ii = 1; ii <= new_addr.size(); ii++ ) {
      if( new_addr[ ii ] == 0 ) {
	if( new_level_start == 0 ) new_level_start = ii;
	core::Size next_free_index = hlevel_.pool_size( new_addr, ii - 1 ) + 1;
	runtime_assert( next_free_index != 0 );
	new_addr[ ii ] = next_free_index;
      }
    }
    core::Size id_num = hlevel_.pool_size( new_addr, hlevel_.nlevels() ) + 1;
    if( new_level_start < new_addr.size() && new_level_start > 0 && id_num != 1 ) {
    	if ( tr.Error.visible() ) {
	   tr.Error << "new level start is " << new_level_start << " so new branches are created, but pool size is NOT zero! "  << hlevel_.pool_size( new_addr, hlevel_.nlevels() ) << " for address: ";
	for( core::Size ii = 1; ii <= new_addr.size(); ii++ ) {
	  tr.Error << new_addr[ ii ] << " ";
	}
	tr.Error << std::endl;
      }
    }
    if( tr.visible() ) tr.Debug << "about to check runtime-asserts: " << id_num << " " << hlevel_.first_zero_pos( new_addr ) << " " << hlevel_.nlevels() << std::endl;

    runtime_assert( ( id_num == 1 && hlevel_.first_zero_pos( new_addr ) <= hlevel_.nlevels() ) ||
		    ( id_num >= 1  && hlevel_.first_zero_pos( new_addr ) > hlevel_.nlevels()  ) ); //if new branches are created, make sure structure is first in pool
    if( tracer_visible_ ) {
      tr.Debug << "pool size is " << id_num << " for addr: ";
      for( core::Size ii = 1; ii <= hlevel_.nlevels(); ii++ ) {
	tr.Debug << new_addr[ ii ] << " | ";
      }
      for( core::Size ii = 1; ii <= hlevel_.nlevels(); ii++ ) {
	tr.Debug << buf_.candidate_address_[ ii ] << " ";
      }
      tr.Debug << std::endl;
    }

    assign_tag( new_addr, id_num, new_candidate_tag );
  }


  core::Size
  MPIHPool_RMSD::evaluate_and_add(
				  core::pose::Pose const& pose,
				  std::string& best_decoy,
				  core::Real& best_rmsd ) {

    using namespace basic;
    PROF_START( basic::MPIBARRIER_BEGIN );
    MPI_Barrier( MPI_COMM_POOL );
    PROF_STOP( basic::MPIBARRIER_BEGIN );

    PROF_START( basic::MPIH_EVAL_CHECK_PROGRESS );
    if( tracer_visible_ ) {
      tr.Debug << "now in evaluate and add, number of structures in top-level pool: " << hlevel_.top_level_pool_size() << " pool-rank: " << pool_rank_ << " pool-size: " << pool_npes_ << std::endl;
    }

    if( nresidues_ != pose.size() ) {
      //~buf_(); //call destructor?
      //if sizes are in-consistent, setup new arrays
      nresidues_ = pose.size();
      runtime_assert( pool_npes_ > 0 );
      buf_.setup( pool_npes_, nresidues_, nlevels_ );
    }

    //check progress. have any nodes finished?

    core::Size num_nodes_finished = any_node_finished();
    if( num_nodes_finished > 0 ) {
      if( pool_npes_ - num_nodes_finished == 0 ) {
      	if( tracer_visible_ ) {
	  tr.Info << "no more nodes running trajectories, returning " << std::endl;
	  tr.Info << "finishing trajectory" << std::endl;
	}
	return -1; //no more nstruct, so simply return
      }
      update_comm( pool_npes_ - num_nodes_finished );
      if( static_cast<int>(current_trajectory_state_) == IN_PROGRESS ) {
	MPI_Comm_rank( MPI_COMM_POOL, ( int* )( &pool_rank_ ) );
	MPI_Comm_size( MPI_COMM_POOL, ( int* )( &pool_npes_ ) );
      }
    }

    PROF_STOP( basic::MPIH_EVAL_CHECK_PROGRESS );
    core::Size best_index = -1;
    if( static_cast<int>(current_trajectory_state_) == IN_PROGRESS ) {
      //evaluate the structure
      best_index = hlevel_.evaluate( pose, best_decoy, current_best_rmsds_, best_address_ );
      if( tracer_visible_ ) {
	tr.Debug << "finished evaluating, best decoy has the tag: " << best_decoy << "\n";
	for( core::Size ii = 1; ii <= current_best_rmsds_.size() ; ii++ ) {
	  tr.Debug << "level=" << ii << " best-rmsd=" << current_best_rmsds_[ ii ] << " level-radii=" << level_radii_[ ii ] << " best-level-address=" << best_address_[ ii ] << std::endl;
	}
	tr.Debug << "done dumping out information about best-rmsd" << std::endl;
      }

      hlevel_.debug_print_size_per_level();
      PROF_START( basic::MPIH_EVAL_COMMUNICATE_NEW );
      //store the highest-resolution rmsd
      core::Size sendcounts = 0;
      //send coordinates and information about evaluation to other nodes
      protocols::toolbox::fill_CA_coords( pose, buf_.coords_ );
      if ( is_new_structure( best_address_, level_radii_,  current_best_rmsds_  ) ) {
	if( tr.visible() ) tr.Debug << "i've found a new structure!" << std::endl;
	runtime_assert(buf_.coords_.u1() > 0 && buf_.coords_.u2() > 0);
	runtime_assert(buf_.coords_.u1() > 0 && buf_.coords_.u2() > 0);
	prepare_send_new_coords( true ); //uses buf_.int_buf1_
	sendcounts = ( 3 * nresidues_ );
      } else {
	prepare_send_new_coords( false ); //uses buf_.int_buf1_
      }
      MPI_Allgather( buf_.int_buf1_, nlevels_, MPI_INT, buf_.neighbor_addresses_, nlevels_, MPI_INT, MPI_COMM_POOL );

      best_rmsd = current_best_rmsds_[ current_best_rmsds_.size() ]; //this is what is written out to traj files
      double candidate_best_rms = best_rmsd;
      //remove this, not needed?
      MPI_Allgather( &candidate_best_rms, 1, MPI_DOUBLE, buf_.candidate_best_rmsds_, 1, MPI_DOUBLE, MPI_COMM_POOL );
      if( tracer_visible_ ) {
	tr.Debug << "address of all-nbrs: " << (nlevels_ * pool_npes_ )
		 << " nlevel: " << nlevels_ << " npes " << pool_npes_ << std::endl;
	for( core::Size ii = 0; ii < ( nlevels_ * pool_npes_  ); ii++ ) {
	  tr.Debug << buf_.neighbor_addresses_[ ii ] << " ";
	}
	tr.Debug << std::endl;
      }
      if( tracer_visible_ ) tr.Debug << "scan output and setup to receive" << std::endl;
      scan_output_and_setup_to_receive(); //scans neighbor_addresses_ and determines nbrs
      if( tracer_visible_ ) tr.Debug << "done calling scan output and setup to receive" << std::endl;
      MPI_Allgatherv( buf_.coords_transfer_buffer_,
		      sendcounts,
		      MPI_DOUBLE,
		      buf_.coords_receiving_buffer_,
		      buf_.int_buf1_, //counts
		      buf_.memory_offset_,//displacements
		      MPI_DOUBLE,
		      MPI_COMM_POOL
		      ); //have to receive all structures

      PROF_STOP( basic::MPIH_EVAL_COMMUNICATE_NEW );

      //bool i_am_a_winning_rank = false;
      buf_.candidate_nbr_index_ = 0;
      utility::vector1< Address > prev_added_addresses;
      utility::vector1< core::Size > prev_added_start_indices;
      Address new_addr;
      num_structures_added_ = 0;

      if( buf_.num_new_neighbors_ > 0 ) {
        bool i_am_a_winning_rank = false;
	std::string new_candidate_tag;
	PROF_START( basic::MPIH_ADD_FIRST_STRUCTURE );
	bool has_new_structure = get_next_candidate(); //find next neighbor and copies into appropriate buffs
	runtime_assert( has_new_structure ); 	//first structure in list is automatically added
	++num_structures_added_;
	new_addr = buf_.candidate_address_;
	core::Size new_level_start = 0;
	Address unresolved = buf_.candidate_address_;
	resolve_address_and_assign_tag( new_addr, new_level_start, new_candidate_tag ); //added 7/28/10

	// save this information for sending to output-node later
	if( is_my_structure() ) {
	  i_am_a_winning_rank = true;
	  buf_.winning_address_ = new_addr;
	  buf_.winning_tag_ = new_candidate_tag;
	  buf_.new_level_begins_ = new_level_start;
	}

	hlevel_.add_new( buf_.candidate_coords_, new_candidate_tag, new_addr );
	if( tracer_visible_ ) {
	  tr.Debug << "adding this structure to hierarchy: ";
	  for( core::Size ii = 1; ii <= new_addr.size(); ii++ ) {
	    tr.Debug << new_addr[ ii ] << " ";
	  }
	  tr.Debug << " " << new_candidate_tag << std::endl;
	}

	prev_added_addresses.push_back( unresolved );
	//keep track of what addresses have been added so far for evaluation of new structures
	//prev_added_addresses.push_back( new_addr );
	std::list<PoolData>::iterator itr;
	hlevel_.level_find( new_addr, hlevel_.nlevels(), itr );
	prev_added_start_indices.push_back( (*itr).pool_->size() );

	PROF_STOP( basic::MPIH_ADD_FIRST_STRUCTURE );
	//after processing the first structure, go through the rest and evaluate against newly discovered
	while( get_next_candidate() ) {
	  PROF_START( basic::MPIH_EVAL_AGAINST_NBR );
	  //setup the universal address for evaluation
	  Address best_addr( hlevel_.nlevels(), 0);
	  best_addr[ 1 ] = 1;
	  //
	  //tr.Debug << "about to evaluate against previously added structures" << std::endl;
	  bool is_new_structure = true;
	  for( core::Size ii = 1; ii <= prev_added_addresses.size(); ii++ ) {
	    utility::vector1< core::Real > candidate_rmsd( hlevel_.nlevels(), 0.0 );
	    bool equal_addresses = true;
	    for( core::Size jj = 1; jj <= buf_.candidate_address_.size(); jj++ ) {
	      if( buf_.candidate_address_[ jj ] != (prev_added_addresses[ ii ])[ jj ] ) {
		equal_addresses = false;
		break;
	      }
	    }
	    if( equal_addresses ) {
	      Address test_addr = prev_added_addresses[ ii ];
	      hlevel_.evaluate( buf_.candidate_coords_,
				new_candidate_tag,
				candidate_rmsd,
				test_addr,
				false,
				false );
	      core::Size last_pos_nonzero = hlevel_.first_zero_pos( test_addr ) - 1;
	      if( candidate_rmsd[ last_pos_nonzero ] < level_radii_[ last_pos_nonzero ] ) {
		//tr.Debug << "lowest resolved level: " << last_pos_nonzero << " for evaluated address: ";
		if( tr.visible() ) {
		  for( core::Size ii = 1; ii <= candidate_rmsd.size(); ii++ ) {
		    tr.Debug << "level-addr=" << test_addr[ ii ] << " rms=" << candidate_rmsd[ ii ] << " radius=" << level_radii_[ ii ] << std::endl;
		  }
		  tr.Debug << "hence this structure is deemed a redundant structure, rejecting!" << std::endl;
		}
		is_new_structure = false;
		break;
	      }
	    } //if( equal_addresses )
	  }
	  //tr.Debug << "done evaluating against previously added structures, tag: " << new_candidate_tag << " and addr: "; print_address( best_addr );

	  if( is_new_structure  ) {
	    ++num_structures_added_;
	    utility::vector1< core::Real > best_test_rmsd( hlevel_.nlevels(), 0.0 );
	    Address new_addr = buf_.candidate_address_;
	    Address unresolved = new_addr;
	    hlevel_.evaluate( buf_.candidate_coords_,
			      new_candidate_tag,
			      best_test_rmsd,
			      new_addr,
			      false,
			      false ); //attempt to resolve any remaining levels which may be sub-clusters of  structures which were just added
	    core::Size new_level_start = 0;
	    if( hlevel_.address_exists_in_cache( new_addr ) ) {
	      new_candidate_tag = "";
	      resolve_address_and_assign_tag( new_addr, new_level_start, new_candidate_tag );
	      hlevel_.add_new( buf_.candidate_coords_, new_candidate_tag, new_addr );

	      if( is_my_structure() ) {
		i_am_a_winning_rank = true;
		buf_.winning_address_ = new_addr;
		buf_.winning_tag_ = new_candidate_tag;
		buf_.new_level_begins_ = new_level_start;
	      }
	      if( tracer_visible_ ) {
		tr.Debug << "adding structure to hierarchy ";
		for( core::Size ii = 1; ii <= new_addr.size(); ii++ ) {
		  tr.Debug << new_addr[ ii ] << " ";
		}
		tr.Debug << " " << new_candidate_tag << std::endl;
	      }

	      core::Size index = find_address( new_addr, prev_added_addresses );
	      if( index > prev_added_addresses.size() ) { //not seen before
		//save the address if it hasn't been seen before
		prev_added_addresses.push_back( unresolved );
		std::list<PoolData>::iterator itr;
		hlevel_.level_find( new_addr, hlevel_.nlevels(), itr );
		prev_added_start_indices.push_back( (*itr).pool_->size() );
	      }
	    }
	  }
	  PROF_STOP( basic::MPIH_EVAL_AGAINST_NBR );
	}
	send_receive_and_write_structures( i_am_a_winning_rank, pose );
	if ( num_structures_added_ > 0 ) {
	  //if new structures have been added after the first eval, update the potentially out-of-date assignment
	  PROF_START( basic::MPIH_UPDATE_EVAL );
	  //tr.Debug << "now doing second evaluation update" << std::endl;
	  //update previous evaluation, if need be.
	  std::string second_update_best_decoy;
	  utility::vector1< core::Real > second_update_rmsd( hlevel_.nlevels(), 0.0 );
	  for( core::Size ii =1 ; ii <= prev_added_addresses.size(); ii++ ) {
	    Address test_addr = prev_added_addresses[ ii ];
	    core::Size second_update_index = hlevel_.evaluate( buf_.coords_,//update original assignment
							       second_update_best_decoy,
							       second_update_rmsd,
							       test_addr,
							       false,
							       false );
	    if( hlevel_.first_zero_pos( test_addr ) == hlevel_.nlevels() + 1 && //only replace if fully resolved
		second_update_rmsd[ second_update_rmsd.size() ] < best_rmsd &&
		second_update_best_decoy.compare("") != 0 ) {
	      best_rmsd = second_update_rmsd[ second_update_rmsd.size() ];
	      best_decoy = second_update_best_decoy;
	      best_index = second_update_index;
	    }
	  }
	}
	PROF_STOP( basic::MPIH_UPDATE_EVAL );
      } else {
	send_receive_and_write_structures( false, pose );
	/**
	//if there are no new structures in my neighborhood, then report that I have no new structures
	PROF_START( basic::MPIH_PREPARE_WRITE_STRUCTURES );
	int num_to_print = 0;
	int num_structures_to_write;
	MPI_Reduce( &num_to_print, &num_structures_to_write, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_POOL );
	PROF_STOP( basic::MPIH_PREPARE_WRITE_STRUCTURES );
	if( pool_rank_ == MPI_OUTPUT_RANK ) {
	  PROF_START( basic::MPIH_WRITE_STRUCT );
	  core::io::silent::SilentFileData sfd;
	  sfd.strict_column_mode( true );
	  receive_and_output_structures( sfd, num_structures_to_write );
	  PROF_STOP( basic::MPIH_WRITE_STRUCT );
	}
	**/
      }
    }
    hlevel_.debug_print_size_per_level();
    PROF_START( basic::MPIBARRIER_END );
    MPI_Barrier( MPI_COMM_POOL );
    PROF_STOP( basic::MPIBARRIER_END );
    return best_index;
  }
  void
  MPIHPool_RMSD::send_receive_and_write_structures(bool i_am_a_winning_rank, core::pose::Pose const& pose) {

  //look in buf_.neighbor_addresses_ for addresses that match my query
    bool use_batch_write_out = false;
    if( use_batch_write_out ) {
      if( tr.visible() ) tr.Debug << "attempting to use batch write-out to dump decoys" << std::endl;
      int print = 0;
      if( i_am_a_winning_rank ) print = 1;
      int* have_structure_to_print = new int[ pool_npes_ ];
      MPI_Allgather( &print, 1, MPI_INT, have_structure_to_print, 1, MPI_INT, MPI_COMM_POOL );
      if( i_am_a_winning_rank ) {
	core::Size num_nbrs = 0;
	core::Size min_ranking_proc = pool_npes_;

	for( core::Size ii = 0; ii < pool_npes_; ii++ ) { //determine which group you belong to
	  if( tr.visible() ) tr.Debug << "looking at position " << ii << " of " << pool_npes_ << std::endl;
	  if( have_structure_to_print[ ii ] == 1 ) {
	    bool same_unresolved_addr = true;
	    for( core::Size jj = 0; jj < hlevel_.nlevels(); jj++ ) {
	      if( static_cast<int>(best_address_[ jj + 1 ]) != buf_.neighbor_addresses_[ (ii * hlevel_.nlevels()) + jj ] ) same_unresolved_addr = false;
	    }
	    if( same_unresolved_addr ) {  //only winning ranks with same unresolved addresses care
	      //DEBUG OUTPUT
	      if( tr.visible() ) {
		tr.Debug << "my current query address: ";
		for( core::Size kk = 1; kk <= best_address_.size(); kk++ ) {
		  tr.Debug << best_address_[ kk ] << " ";
		}
		tr.Debug << "  compared to nbr rank: " << ii << " ";
		for( core::Size kk = 0; kk < hlevel_.nlevels(); kk++ ) {
		  tr.Debug << buf_.neighbor_addresses_[ (ii * hlevel_.nlevels()) + kk ] << " ";
		}
		tr.Debug << "are the same and will be grouped together for writing to file" << std::endl;
	      }
	      //END DEBUG OUTPUT
	      if( ii < min_ranking_proc ) min_ranking_proc = ii;
	      num_nbrs++;
	    } else {
	      if( tr.visible() ) tr.Debug << "no matching addresses, will not group any processes together" << std::endl;
	      have_structure_to_print[ ii ] = 0;
	    }
	  }
	}
	//how many levels to print?
	if( pool_rank_ == min_ranking_proc ) {
	  if( tr.visible() ) tr.Debug << "I am the min ranking proc: " << pool_rank_ << " and i am printing out these rank neighbors: ";
	  for( core::Size ii = 0; ii < num_nbrs; ii++ ) {
	    if( have_structure_to_print[ ii ] == 1 ) {
	      if( tr.visible() ) {
		tr.Debug << " nbr rank: " << ii << " ";
		for( core::Size jj = 0; jj < hlevel_.nlevels(); jj++ ) {
		  tr.Debug << buf_.neighbor_addresses_[ ( ii * hlevel_.nlevels() ) + jj ] << " ";
		}
	      }
	    }
	  }
	  if( tr.visible() ) tr.Debug << std::endl;
		core::io::silent::SilentFileOptions opts;
	  core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	  core::io::silent::SilentFileData sfd( opts );
	  ss->fill_struct( pose, buf_.winning_tag_ );
	  if( first_time_writing_ ) {
	    first_time_writing_ = false;
	    write_headers_to_hierarchy( ss );
	  }
	  write_decoys_to_hierarchy( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ );
	  num_nbrs--;
	  if( num_nbrs > 0 ) {
	    receive_and_output_structures( sfd, num_nbrs );
	  }
	} else if( pool_rank_ != min_ranking_proc ) {
	  if( tr.visible() ) tr.Debug << "sending my structure to min-ranking-proc: " << min_ranking_proc << std::endl;
		core::io::silent::SilentFileOptions opts;
	  core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	  core::io::silent::SilentFileData sfd( opts );
	  ss->fill_struct( pose, buf_.winning_tag_ ); //ek debug 7/30/10 uncomment when everything is fixed
	  send_silent_struct_to_rank( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ , min_ranking_proc );
	}
      }
    } else {

      PROF_START( basic::MPIH_PREPARE_WRITE_STRUCTURES );
      int num_structures_to_write = 0;
      //report whether or not i have a structure to write out.
      if( i_am_a_winning_rank ) {
	int num_to_print = 1;
	MPI_Reduce( &num_to_print, &num_structures_to_write, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_POOL );
      } else {
	int num_to_print = 0;
	MPI_Reduce( &num_to_print, &num_structures_to_write, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_POOL );
      }
      PROF_STOP( basic::MPIH_PREPARE_WRITE_STRUCTURES );
      if( static_cast<int>(pool_rank_) == MPI_OUTPUT_RANK ) {
	//if I am the output rank, receive and write out structures
	if( tracer_visible_ ) { tr.Debug << "expecting to write out " << num_structures_to_write << " structures" << std::endl; }
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData sfd( opts );
	sfd.strict_column_mode( true );
	PROF_START( basic::MPIH_WRITE_STRUCT );

	if( i_am_a_winning_rank ) {
	  //write out my own structure first
	  core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	  ss->fill_struct( pose, buf_.winning_tag_ );
	  if( first_time_writing_ ) {
	    first_time_writing_ = false;
	    write_headers_to_hierarchy( ss );
	  }
	  write_decoys_to_hierarchy( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ );
	  num_structures_to_write--;
	}
	receive_and_output_structures( sfd, num_structures_to_write );
	PROF_STOP( basic::MPIH_WRITE_STRUCT );
      } else {
	//if I am not the output rank, then send my structures to the output rank
	PROF_START( basic::MPIH_WRITE_STRUCT );
	if( i_am_a_winning_rank ) {
		core::io::silent::SilentFileOptions opts;
	  core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	  core::io::silent::SilentFileData sfd( opts );

	  ss->fill_struct( pose, buf_.winning_tag_ ); //ek debug 7/30/10 uncomment when everything is fixed
	  send_silent_struct_to_rank( sfd, ss, buf_.winning_address_, buf_.new_level_begins_ );
	}
	PROF_STOP( basic::MPIH_WRITE_STRUCT );
      }
    }
  }

  void
  MPIHPool_RMSD::buf_to_address( Address & addr, int* addr_buf, core::Size index ) {
    for( core::Size ii = 1; ii <= addr.size(); ii++ ) {
      addr[ ii ] = addr_buf[ index + ii - 1 ];
    }
  }

  void
  MPIHPool_RMSD::address_to_buf( Address & addr, int* addr_buf, core::Size index ) {
    for( core::Size ii = 1; ii <= addr.size(); ii++ ) {
      addr_buf[ index + ii - 1 ] = addr[ ii ];
    }
  }

  bool
  MPIHPool_RMSD::is_my_structure() {
    return ( ( buf_.candidate_nbr_index_ - 1 ) == pool_rank_);
  }


  bool
  MPIHPool_RMSD::is_new_structure( Address & address,
				   utility::vector1< core::Real > &,
				   utility::vector1< core::Real > & rmsds ) {
    for( core::Size ii = 1; ii <= address.size(); ii++ ) {
      if( address[ ii ] == 0 ) {
	if( tr.visible() ) tr.Debug << "address at level " << ii << " is 0, so is a new structure " << std::endl;
	return true;
      }
    }
    if( rmsds[ rmsds.size() ] > hlevel_.level_radius( hlevel_.nlevels() ) ) {
      if( tr.visible() ) tr.Debug << "rms at last-level is: " << rmsds[ rmsds.size() ] << " which is greater than radius, " << hlevel_.level_radius( hlevel_.nlevels() ) << " so is a new structure" << std::endl;
      return true;
    }
    if( tr.visible() ) tr.Debug << " structure is not a new structure" << std::endl;
    return false;
  }

  core::Real
  MPIHPool_RMSD::resolved_level_best_rmsd( Address & addr, utility::vector1< core::Real > & rmsd ) {
    core::Real best_rms = 0;
    for( core::Size ii = 1; ii <= addr.size(); ii++ ) {
      if( addr[ ii ] != 0 ) {
	best_rms = rmsd[ ii ];
      }
    }
    return best_rms;
  }

  bool
  MPIHPool_RMSD::is_new_structure( Address & address,
				   utility::vector1< core::Real > & radii,
				   core::Real & rmsd ) {
    for( core::Size ii = 1; ii <= address.size(); ii++ ) {
      if( address[ ii ] == 0 ){
	//tr.Debug << "address at level " << ii << " is 0, so is a new structure " << std::endl;
	return true;
      }
    }
    if( rmsd > radii[ radii.size() ] ) {
      //tr.Debug << "rmsd at last level: " << rmsds[ rmsds.size() ] << " is greater than threshold: " << radii[ radii.size() ] << " so is a new structure " << std::endl;
      return true;
    }
    //tr.Debug << "structure is not a new structure" << std::endl;
    return false;
  }

  void
  MPIHPool_RMSD::print_address( Address & addr ) {
    if( tr.visible() ) {
      for( core::Size ii = 1; ii <= addr.size(); ii++ ) {
	tr.Debug << addr[ ii ] << " ";
      }
      tr.Debug << std::endl;
    }
  }

  core::Size
  MPIHPool_RMSD::find_address( Address & query, utility::vector1< Address > & addr_database ) {
    PROF_START( basic::HIERARCHY_FIND_ADDRESS );
    core::Size  index;
    for(  index = 1; index <= addr_database.size(); index++ ) {
      bool found_a_match = true;
      for( core::Size jj = 1; jj <= query.size(); jj++ ) {
	if( query[ jj ] != (( Address )addr_database[ index ])[ jj ] ) {
	  found_a_match = false;
	}
      }
      if( found_a_match ) {
	break;
      }
    }
    PROF_STOP( basic::HIERARCHY_FIND_ADDRESS );
    return index;
  }

  core::Size
  MPIHPool_RMSD::any_node_finished(){
    buf_.int_buf1_[ 0 ] = current_trajectory_state_;
    MPI_Allgather( buf_.int_buf1_, 1, MPI_INT, buf_.finished_, 1, MPI_INT, MPI_COMM_POOL );
    core::Size num_nodes_finished = 0;
    core::Size index_in_prog = 0;
    for( unsigned int ii = 0; ii < pool_npes_; ii++ ) {
      if( buf_.finished_[ ii ] == FINISHED ) {
	num_nodes_finished++;
      }else {
	buf_.int_buf1_[ index_in_prog++ ] = ii;
      }
    }
    if( tr.visible() ) {
      tr.Debug << "number of nodes finished this round: " << num_nodes_finished << std::endl;
    }
    return num_nodes_finished;
  }

  void
  MPIHPool_RMSD::update_comm( core::Size newsize ) {
    if( tr.visible() ) {
      tr.Debug << "some trajectories finished, creating new comm with " << newsize  << std::endl;
    }
    create_comm( buf_.int_buf1_, newsize );
    //~buf_();
    buf_.setup( newsize, nresidues_, nlevels_ );
    if( static_cast<int>(current_trajectory_state_) == IN_PROGRESS ) {
      MPI_Comm_rank( MPI_COMM_POOL, ( int* )( &pool_rank_ ) );
      MPI_Comm_size( MPI_COMM_POOL, ( int* )( &pool_npes_ ) );
      tr.Info << "remaining ranks has pool-size of " << pool_npes_ << " and rank: " << pool_rank_ << std::endl;
    }
  }


  void
  MPIHPool_RMSD::prepare_send_new_coords( bool send_coords ){
    if( send_coords ) {
      runtime_assert( nlevels_ == best_address_.size() );
      for( core::Size ii = 1; ii <= best_address_.size(); ii++ ) {
	//tr.Debug << "writing to " << (ii-1) << " of " << best_address_.size() - 1 << std::endl;
	buf_.int_buf1_[ ii - 1 ] = best_address_[ ii ];
      }
      buf_.farray_to_array( 0, buf_.coords_, buf_.coords_transfer_buffer_ );
      //only one structure at a time, index is always 0
    } else {
      for( core::Size ii = 0; ii < nlevels_; ii++ ){
	buf_.int_buf1_[ ii ] = -1;
      }
    }
  }


  bool
  MPIHPool_RMSD::get_next_candidate() {
    PROF_START( basic::HIERARCHY_GET_NEXT_CANDIDATE );
    //buf_.candidate_nbr_index_ starts from 1, equivalent to pool_rank_ + 1
    if( buf_.candidate_nbr_index_ < (buf_.is_a_neighbor_).size() ) { //not at the end
      core::Size itr;
      for( itr = buf_.candidate_nbr_index_ + 1; itr <= (buf_.is_a_neighbor_).size(); itr++ ) {
	if( buf_.is_a_neighbor_[ itr ] ) {
	  buf_.array_to_farray( buf_.memory_offset_[ itr - 1 ], buf_.candidate_coords_, buf_.coords_receiving_buffer_ );
	  for( core::Size ii = 1; ii <= nlevels_; ii++ ) {
	    buf_.candidate_address_[ ii ] = buf_.neighbor_addresses_[ ( ( itr - 1 ) * nlevels_ ) + ( ii - 1 ) ];
	  }
	  buf_.candidate_best_rmsd_ = buf_.candidate_best_rmsds_[ itr - 1 ];
	  buf_.candidate_nbr_index_ = itr;
	  if( tr.visible() ) {
	    tr.Debug << "next examining address: ";
	    for( core::Size ii = 1; ii <= buf_.candidate_address_.size(); ii++ ) {
	      tr.Debug << buf_.candidate_address_[ ii ] << " ";
	    }
	    tr.Debug << std::endl;
	  }
	  break;
	}
      }
      PROF_STOP(  basic::HIERARCHY_GET_NEXT_CANDIDATE );
      if( itr > (buf_.is_a_neighbor_).size()) { return false; } // no neighbors
      return true;
    } else {
      return false;
    }
  }

  void
  MPIHPool_RMSD::receive_silent_struct_any_source( core::io::silent::SilentFileData& recv_ss, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level_begins ) {
    using namespace core::io::silent;
    using namespace basic;
    PROF_START( basic::HIERARCHY_RECV_COORDS );
    recv_ss.clear();
    std::istringstream os;
    int string_size = 0;
    MPI_Status stat;
    std::string received_string;
    int* output_info = new int[ 2 + ss_addr.size() ];
    //maybe receive address string_size

    //tr.Debug << "receiving this many ints: " << (2 + ss_addr.size() ) << std::endl;
    MPI_Recv( output_info, ( 2 + ss_addr.size() ), MPI_INT, MPI_ANY_SOURCE, OUTPUT_TAG, MPI_COMM_POOL, &stat );
    string_size = output_info[ 0 ];
    char *cbuf = new char[ string_size + 1 ];
    int sending_rank = stat.MPI_SOURCE;
    for( core::Size ii =1; ii <= ss_addr.size(); ii++ ) {
      ss_addr[ ii ] = output_info[ ii ];
    }
    new_level_begins = output_info[ ss_addr.size() + 1 ];
    //tr.Debug << "receiving this many chars: " << string_size << std::endl;
    MPI_Recv( cbuf, string_size, MPI_CHAR, sending_rank, OUTPUT_TAG, MPI_COMM_POOL, &stat );
    received_string.assign( cbuf, string_size );
    os.str( received_string );
    utility::vector1< std::string > tags;
    recv_ss.read_stream( os, tags, false );
    ss = *(recv_ss.begin());
    if( tr.visible() ) {
      tr.Debug << "just received structure from rank:  " << sending_rank
	       << " with tag " << ss->decoy_tag() << " new_level begins at: " << new_level_begins << std::endl;
    }
    delete[] cbuf;
    PROF_STOP( basic::HIERARCHY_RECV_COORDS );
  }

  void
  MPIHPool_RMSD::send_silent_struct_to_rank( core::io::silent::SilentFileData& send_ss, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level, core::Size receiving_rank ) {
    using namespace core::io::silent;
    using namespace basic;
    PROF_START( basic::HIERARCHY_SEND_COORDS );
    send_ss.clear();
    std::ostringstream os;
    send_ss._write_silent_struct( *ss, os );
    int* output_info = new int[ 2 + ss_addr.size() ];
    //std::string decoy_tag = ss->decoy_tag();
    int string_size =  (os.str()).size();
    output_info[ 0 ] = string_size;
    for( core::Size ii = 1; ii <= ss_addr.size(); ii++ ) {
      output_info[ ii ] = ss_addr[ ii ];
    }
    //runtime_assert( new_level > 0 );
    output_info[ ss_addr.size() + 1 ] = new_level;
    if ( tr.visible() ) {
      tr.Debug << "sending decoy with tag: " << ss->decoy_tag() << " new-level is " << new_level << std::endl;
    }
    MPI_Send(output_info, ( 2 + ss_addr.size() ), MPI_INT, receiving_rank, OUTPUT_TAG, MPI_COMM_POOL );
    //tr.Debug << "sending this many chars: " << string_size << std::endl;
    MPI_Send(const_cast<char*> (os.str().data()), string_size, MPI_CHAR, receiving_rank, OUTPUT_TAG, MPI_COMM_POOL);
    PROF_STOP( basic::HIERARCHY_SEND_COORDS );
  }

  void
  MPIHPool_RMSD::send_silent_struct_to_rank( core::io::silent::SilentFileData& send_ss, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level ) {
    send_silent_struct_to_rank( send_ss, ss, ss_addr, new_level, MPI_OUTPUT_RANK );
  }


  void
  MPIHPool_RMSD::scan_output_and_setup_to_receive(){ //determines the number of new neighbors
    buf_.num_new_neighbors_ = 0;
    PROF_START( basic::HIERARCHY_SETUP_TO_RECV );
    buf_.is_a_neighbor_.resize( pool_npes_, false );

    runtime_assert(buf_.is_a_neighbor_.size() == pool_npes_ );
    //convert to proper address
    Address nbr_address( nlevels_, 0 );

    for( unsigned int i = 0; i < pool_npes_; i++  ) {
      if( buf_.neighbor_addresses_[ (i * nlevels_ ) ] != -1 ) { //quick indication that no pose discovered this round
	for( core::Size ii = 1; ii <= nbr_address.size(); ii++ ) {
	  nbr_address[ ii ] = buf_.neighbor_addresses_[ ( i * nlevels_ ) + ( ii - 1 ) ];
	}
	if( hlevel_.address_exists_in_cache( nbr_address ) ) {
	  buf_.is_a_neighbor_[ i + 1 ] = true;
	  buf_.num_new_neighbors_++;
	} else {
	  buf_.is_a_neighbor_[ i + 1 ] = false;
	}
      } else {
	buf_.is_a_neighbor_[ i + 1 ] = false;
      }
    }
    if( tr.visible() ) {
      tr.Debug << "my current query address: ";
      for( core::Size ii = 1; ii <= best_address_.size(); ii++ ) {
	tr.Debug << best_address_[ ii ] << " ";
      }
      tr.Debug << std::endl;

      hlevel_.debug_print_size_per_level();
      tr.Debug << "these are the addresses I will be examining: ";
      for( core::Size ii = 1; ii <= buf_.is_a_neighbor_.size(); ii++ ) {
	if( buf_.is_a_neighbor_[ ii ] == 1 ) {
	  for( core::Size jj = 1; jj <= nbr_address.size(); jj++ ) {
	    tr.Debug << buf_.neighbor_addresses_[ ( (ii-1) * nlevels_ ) + ( jj - 1 ) ] << " ";
	  }
	  tr.Debug << ", ";
	}
      }
      tr.Debug << std::endl;
    }

    core::Size index = 0;
    core::Size receive_counts = 0;
    //setup coords to receive stuff
    for( core::Size ii = 0; ii < (pool_npes_ * nlevels_ ); ii+=nlevels_ ) {
      if ( buf_.neighbor_addresses_[ ii ] > 0 ) { //is it a valid address?
	buf_.int_buf1_[ index ] = ( nresidues_ * 3 ); //counts
	buf_.memory_offset_[ index ] = receive_counts; //displacement
	receive_counts += ( nresidues_ * 3 );
      }else{
	buf_.int_buf1_[ index ] = 0;
	buf_.memory_offset_[ index ] = 0;
      }
      index++;
    }
    PROF_STOP( basic::HIERARCHY_SETUP_TO_RECV );
  }

  void
  MPIHPool_RMSD::set_discovered_out( std::string const& newout){
    new_decoys_out_ = newout;
  }

  std::string const&
  MPIHPool_RMSD::get_discovered_out(){
    return new_decoys_out_;
  }

  void
  MPIHPool_RMSD::set_transition_threshold( core::Real threshold ){
    //no need for this function in this class??
    runtime_assert( threshold >= 0 );
  }

  void
  MPIHPool_RMSD::set_nresidues( core::Size nres ){
    nresidues_ = nres;
  }

  core::Size
  MPIHPool_RMSD::get_nresidues(){
    return nresidues_;
  }

  void
  MPIHPool_RMSD::write_headers_to_hierarchy( core::io::silent::SilentStructOP& ss ) {
    hlevel_.write_headers_to_hierarchy( ss );
  }

  void
  MPIHPool_RMSD::write_decoys_to_hierarchy( core::io::silent::SilentFileData& sfd, core::io::silent::SilentStructOP& ss, Address& ss_addr, core::Size new_level_begins ) {
    using namespace basic;
    PROF_START( basic::WRITE_DECOYS_TO_HIERARCHY );
    Address tmp_addr = ss_addr;
    if( new_level_begins == 0 ) new_level_begins = tmp_addr.size() + 1;
    if( tr.visible() ) {
      tr.Debug << "writing decoy to hierarchy: " << ss->decoy_tag() << " ";
      for( core::Size ii =1; ii <= ss_addr.size(); ii++ ) {
	tr.Debug << ss_addr[ ii ] << " ";
      }
      tr.Debug << " new_level: " << new_level_begins << std::endl;
    }

    for( core::Size ii = new_level_begins; ii <= tmp_addr.size(); ii++ ) {
      tmp_addr[ ii ] = 0;
    }
    core::Size index = new_level_begins;
    do {
      std::string file_in_hierarchy = hlevel_.lib_full_path( tmp_addr );
      utility::file::FileName file( file_in_hierarchy );
      if( !utility::file::file_exists( file_in_hierarchy ) ){
	utility::file::create_directory_recursive( file.path() );
	utility::file::create_blank_file( file.name() );
	std::ofstream os;
	os.open( (file.name()).c_str() );
	ss->print_header( os ); //temporary fix
	os.close();
      }
      if( tr.visible() ) {
	tr.Debug << index << " writing decoy " << ss->decoy_tag() << " to file: " << file_in_hierarchy << std::endl;
      }
      std::ofstream os;
      os.open( (file_in_hierarchy).c_str(), std::ios::app );
      sfd._write_silent_struct( *ss, os, false );
      os.close();
      if( index <= ss_addr.size() ) {
	tmp_addr[ index ] = ss_addr[ index ];
      }
      //index++;
    } while( index++ <= ss_addr.size() );
    sfd.write_silent_struct( *ss, new_decoys_out_, false );
    //tr.Debug << "done dumping structure into hierarchy " << std::endl;
    PROF_STOP( basic::WRITE_DECOYS_TO_HIERARCHY );
  }

  void
  MPIHPool_RMSD::max_cache_size( core::Size max_cache ) {
    hlevel_.max_cache_size( max_cache );
  }

  void
  MPIHPool_RMSD::finalize(){
    PROF_START( basic::FINALIZE );
    MPI_Barrier( MPI_COMM_POOL );
    current_trajectory_state_ = FINISHED;
    tr.Info << "rank: " << pool_rank_ << " calling finalized on trajectory " << std::endl;
    core::Size num_nodes_finished = any_node_finished();
    if( tr.visible() ) {
      tr.Debug << "dumping hierarchy (debug) " << std::endl;
      hlevel_.debug_print_hierarchy();
    }
    if( num_nodes_finished > 0 ) {
      tr.Info << "num nodes finished this round: " << num_nodes_finished << " of " << pool_npes_ << std::endl;
      if( ( pool_npes_ - num_nodes_finished ) == 0 ) {
	return;
      }
      update_comm( pool_npes_ - num_nodes_finished );
    }
    PROF_STOP( basic::FINALIZE );
  }


  void
  MPIHPool_RMSD::initialize(){
    using namespace core;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    using namespace basic;

    tr.Info << "initializing hierarchical MPI class" << std::endl;
    //time setup
    clock_t start_setup = clock();
    rank_ = 0;
    npes_ = 0;

    PROF_START( basic::CHECK_COMM_SIZE );
    MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank_ ) );
    MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &npes_ ) );
    PROF_STOP( basic::CHECK_COMM_SIZE );

    PROF_START( basic::INITIALIZE );
    pool_rank_ = rank_;
    pool_npes_ = npes_;
    //tr.Debug << "just checked rank: it's " << rank_ << " out of " << npes_ << std::endl;

    //initialize radii
    level_radii_ = option[ cluster::K_radius ]();
    runtime_assert( level_radii_.size() == nlevels_ );

    protocols::jd2::MPIFileBufJobDistributor* jd2 =
      dynamic_cast< protocols::jd2::MPIFileBufJobDistributor * > (protocols::jd2::JobDistributor::get_instance());
    core::Size min_client_rank, new_size;
    if( jd2 ) {
      min_client_rank = jd2->min_client_rank();
      new_size = npes_ - min_client_rank;
    } else {
      utility_exit_with_message("cannot use MPIHPool_RMSD without using the MPIFileBufJobDistributor! try again!");
    }

    pool_size_ = Pool_RMSD::size();
    nresidues_ =  Pool_RMSD::natom();
    current_address_.resize( nlevels_, 0 );

    PROF_STOP( basic::INITIALIZE );
    //create new MPI_COMM_WORLD based on sub-set of nodes
    tr.Info << "initializing comm" << std::endl;
    PROF_START( basic::MPICOMMCREATION );
    int index = 0;
    if( tr.visible() ) tr.Debug << "now trying to set up new communicator" << std::endl;
    buf_.setup( pool_npes_, nresidues_, nlevels_ );
    for(unsigned int ii = min_client_rank; ii < npes_; ii++){
      (buf_.int_buf1_)[ index++ ] = ii;
      if( tr.visible() ) tr.Debug << "including " << ii << " in the new pool communicator " << std::endl;
    }

    //initialize all num_slave dependent buffers for MPI transfers
    MPI_Group pool_group, all;
    int returnval;

    returnval = MPI_Comm_group( MPI_COMM_WORLD, &all);
    if ( returnval != MPI_SUCCESS ) {
      utility_exit_with_message("failed in creating a new communicator!");
    }
    if( tr.visible() ) tr.Debug << "set up MPI_COMM_WORLD group" << std::endl;

    returnval = MPI_Group_incl( all, (new_size), buf_.int_buf1_, &pool_group );
    if ( returnval != MPI_SUCCESS ) {
      utility_exit_with_message("failed in creating a new communicator!");
    }
    if( tr.visible() )tr.Debug << "created the pool group" << std::endl;

    returnval = MPI_Comm_create( MPI_COMM_WORLD, pool_group, &MPI_COMM_POOL );
    if ( returnval != MPI_SUCCESS ) {
      utility_exit_with_message("failed in creating a new communicator!");
    }
    if( tr.visible() ) tr.Debug << "created the MPI_COMM_POOL communicator " << std::endl;

    if( rank_ >= min_client_rank ) {
      MPI_Comm_rank( MPI_COMM_POOL, ( int* )( &pool_rank_ ) );
      MPI_Comm_size( MPI_COMM_POOL, ( int* )( &pool_npes_ ) );
      if( tr.visible() ) tr.Debug << "new ranks from MPI_COMM_POOL: " << pool_rank_ << " " << pool_npes_ << std::endl;
    }

    PROF_STOP( basic::MPICOMMCREATION );
    if( tr.visible() ) tr.Debug << "finished initializing, setting up new comm, and ready to go!" << std::endl;
    PROF_START( basic::INITIALIZE_HIERARCHY );
    tracer_visible_ = tr.visible();

    Address universal_address( nlevels_, 0 );
    universal_address[ 1 ] = 1;

    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
    tr.Info << "checking for hierarchy..." << std::endl;
    if( tr.visible() ) tr.Debug << "does pool exist? " << hlevel_.lib_full_path( universal_address ) << " : " << hlevel_.pool_exists( universal_address ) << std::endl;
    if( !hlevel_.pool_exists( universal_address ) ) { //checks for directory assumed to contain hierarchy
      if( rank_ >= min_client_rank ) {
	MPI_Barrier( MPI_COMM_POOL );
	tr.Info << "hierarchical pool doesn't exist. creating hierarchy " << std::endl;
	//create a hierarchical pool to go along with this silent-file
	if( option[ mc::read_structures_into_pool ].user() ) {
	  if( !utility::file::file_exists( option[ mc::read_structures_into_pool ]() ) ) {
	    utility_exit_with_message("you specified a file for option mc::read_structures_into_pool that does not exist! try again!");
	  } else {
			core::io::silent::SilentFileOptions opts;
	    core::io::silent::SilentFileData sfd( opts );
	    sfd.strict_column_mode( true );
	    core::pose::Pose pose;
	    bool successfully_read_file = false;
	    core::Size max_attempts = 100;
	    core::Size attempt = 0;
	    while( !successfully_read_file && attempt < max_attempts ) {
	      successfully_read_file = sfd.read_file( hlevel_.filename() );
	    }
	    if ( successfully_read_file ) {
	      utility::vector1< std::string > tags = sfd.tags();
	      for( core::Size ii = 1; ii <= tags.size(); ii++ ) {
		//core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_in();
		core::io::silent::SilentStructOP ss = sfd[ tags[ ii ] ];
		std::string best_decoy;
		utility::vector1< core::Real > best_rms;
		Address best_addr;
		core::Size new_level_starts_at = 0;
		hlevel_.evaluate( (*ss), best_decoy, best_rms, best_addr );
		for( core::Size jj = 1; jj <= best_addr.size(); jj++ ) {
		  if( best_addr[ jj ] == 0 ) {
		    if( new_level_starts_at == 0 ) new_level_starts_at = jj;
		    best_addr[ jj ] = hlevel_.pool_size( best_addr, jj - 1 ) + 1;
		  }
		}
		if( rank_ == min_client_rank ) {
		  hlevel_.add_new( (*ss), tags[ ii ], best_addr, true, 2 );
		} else {
		  hlevel_.add_new( (*ss), tags[ ii ], best_addr, false, 2 );
		}
	      }
	    } else {
	      utility_exit_with_message("cannot read silent file. thus cannot create hierarchy! filename: " + hlevel_.filename());
	    }
	  }
	}
      }
    } else {
      if( rank_ >= min_client_rank ) {
	tr.Info << "pool already exists!! loading file as top-level" << std::endl;
	Pool_RMSD_OP top_level_pool( new Pool_RMSD( hlevel_.lib_full_path( universal_address ) ) );
	if( tr.visible() ) tr.Debug << "finished reading pool from file: " << hlevel_.lib_full_path( universal_address ) << std::endl;
	hlevel_.fill_top_level( top_level_pool );
      }
      if( rank_ >= min_client_rank ) {
	MPI_Barrier( MPI_COMM_POOL );
      }
    }
    tr.Info << "this rank: " << pool_rank_ << " starting out with " << hlevel_.top_level_pool_size() << " structures in the top-most-level of hierarchy" << std::endl;
    if( tr.visible() ) {
      hlevel_.debug_print_size_per_level();
      tr.Debug << "END PRINTING SIZE PER LEVEL FOR INITIATION" << std::endl;
    }
    clock_t finish_setup = clock();
    double time_to_setup = ( double(finish_setup) - double(start_setup) ) / CLOCKS_PER_SEC;
    tr << "time to setup " << pool_npes_ << " nodes: " << time_to_setup << std::endl;
    PROF_STOP( basic::INITIALIZE_HIERARCHY );
  }

  bool
  MPIHPool_RMSD::is_in_neighborhood( Address & q_address, Address & ref_address ) {

    runtime_assert( q_address.size() == ref_address.size() );

    for( core::Size ii = 1; ii <= ref_address.size(); ii++ ) {
      if( q_address[ ii ] != 0 && ref_address[ ii ] != 0 && //either is un-resolved
	  q_address[ ii ] != ref_address[ ii ] ) {
	return false;
      }
    }
    return true;
  }

  void
  MPIHPool_RMSD::address_to_string( Address & address_buf, core::Size, std::string & address_tag ) {
    address_tag = "";
    std::ostringstream q;
    q.width(5);
    q.fill('0');
    for( core::Size ii = 1; ii <= address_buf.size(); ii++) {
      q << address_buf[ ii ];
      address_tag += q.str() + ".";
    }
    //tr.Debug << "resulting address from address_to_string is: " << address_tag << std::endl;
  }

  void
  MPIHPool_RMSD::string_to_address( Address & address_buf, core::Size index, std::string & address_tag ){
    //erase prefix, should be something like "c." or "new."
    core::Size pos = 0, newpos;
    std::string subtag="";
    core::Size prefix_pos = address_tag.find('.',0);
    address_tag.erase(0,prefix_pos);
    core::Size counted_levels = 0;
    while( counted_levels < nlevels_ ) {
      newpos = address_tag.find( '.', pos+1 );
      if( newpos == address_tag.length() ) {
	break; // no more "."
      }

      subtag = address_tag.substr( pos + 1, (newpos-pos - 1) );

      pos = newpos;
      core::Size first_nonzero_pos = subtag.find_first_not_of('0',0);
      if (first_nonzero_pos > 0 ) {
	subtag.erase(0,first_nonzero_pos);
      }
      address_buf[ index + counted_levels ] = atoi(subtag.c_str());
      counted_levels++;
    }
    while( counted_levels < nlevels_ ){
      address_buf[ index + counted_levels ] = 0;
      counted_levels++;
    }

  }

  void
  MPIHPool_RMSD::assign_tag( Address& address_tag, core::Size id_num, std::string & newtag ){
    std::ostringstream q;
    for( core::Size ii = 1; ii <= address_tag.size(); ii++ ) {
      if( address_tag[ ii ] == 0 ) {
	core::Size prev_level = ii - 1;
	q << hlevel_.pool_size( address_tag, (prev_level) ) + 1 << ".";
	if( tr.visible() ) tr.Debug << " at level " << ii << " assigning " << ( hlevel_.pool_size( address_tag, ii - 1 ) + 1 ) << std::endl;
      } else {
	q << address_tag[ ii ] << ".";
      }
    }
    q << id_num;
    newtag = "new." + q.str();
    //tr.Debug << " id_num: " << id_num << " newtag: " << newtag << " address: ";
    //for( core::Size ii = 1; ii <= address_tag.size(); ii++ ) {
      //tr.Debug << address_tag[ ii ] << " ";
    //}
    //tr.Debug << std::endl;

  }

  void
  MPIHPool_RMSD::assign_tag( std::string const& address_tag, core::Size assigned_id_num, std::string& newtag ){
    std::ostringstream q;
    q << assigned_id_num;
    newtag = "new." + address_tag + "." + q.str();
    //tr.Debug << "from address: " << address_tag << " id: " << assigned_id_num << " produced: " << newtag << std::endl;
  }


  void  MPIHPool_RMSD::create_comm( int* ranks_to_include, int new_size ){
    int returnval;
    MPI_Group new_pool_group, old_pool_group;
    MPI_Comm dup_pool_comm;
    PROF_START( basic::MPICOMMCREATION );
    MPI_Comm_dup( MPI_COMM_POOL, &dup_pool_comm );
    returnval = MPI_Comm_group( dup_pool_comm, &old_pool_group );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved!" );

    returnval = MPI_Group_incl( old_pool_group, new_size, ranks_to_include, &new_pool_group );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved!" );

    returnval = MPI_Comm_create( dup_pool_comm, new_pool_group, &MPI_COMM_POOL );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved!" );
    PROF_STOP( basic::MPICOMMCREATION );

  }


#endif //useMPI?
}
}
}
