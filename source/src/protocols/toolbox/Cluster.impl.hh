// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file toolbox/Cluster.impl.hh
/// @brief template implementation for clustering of silentstructs that are provided by iterators
/// @author Oliver Lange

#ifndef INCLUDED_protocols_toolbox_Cluster_impl_hh
#define INCLUDED_protocols_toolbox_Cluster_impl_hh


// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <deque>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/toolbox/Cluster.hh>
// AUTO-REMOVED #include <protocols/toolbox/DecoySetEvaluation.impl.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/toolbox/DecoySetEvaluation.hh>

namespace protocols {
namespace toolbox {

static thread_local basic::Tracer _impl_tr( "protocols.toolbox.cluster" );

template< typename SilentStructIterator, typename StructureContainer >
void cluster_silent_structs(
			    core::Size n_decoys,
			    SilentStructIterator input_decoys_begin,
			    SilentStructIterator input_decoys_end, //it->SilentStructOP
			    StructureContainer& new_structs, //provides a "push_back" method
			    ClusterOptions opts
) {
  //read CA coords into DecoySetEvaluation
  DecoySetEvaluation CA_set;
  CA_set.push_back_CA_xyz_from_silent_file( n_decoys, input_decoys_begin, input_decoys_end, true /*store energies*/ );

  //cluster the decoys according to CA_rmsd
  cluster_silent_structs( CA_set, input_decoys_begin, input_decoys_end, new_structs, opts );
}

template< typename SilentStructIterator, typename StructureContainer >
void cluster_silent_structs( DecoySetEvaluation const& CA_set,
  SilentStructIterator input_decoys_begin,
  SilentStructIterator input_decoys_end,
  StructureContainer& new_structs,
  ClusterOptions opts
) {

  core::Size const n_decoys( CA_set.n_decoys() );

  // initialize cluster object
  toolbox::ClusterPhilStyle cluster( n_decoys , opts.cluster_radius );

  _impl_tr.Info << "compute distance matrix" << std::endl;
  // compute distance matrix
  CA_set.compute_distance_matrix( cluster.distance_matrix() );

  // don't limit the maximum cluster size
  cluster.set_n_max_cluster( 1000000 );

  _impl_tr.Info << "compute clusters" << std::endl;
  cluster.do_clustering();

  // for the next steps we need energies --- could be made optional
  _impl_tr.Info << "sort clusters by energy" << std::endl;
  runtime_assert( CA_set.all_energies().size() == n_decoys );
  cluster.sort_each_group_by_energy( CA_set.all_energies(), opts.keep_center );
  if ( opts.limit_cluster_size > 0 ) {
    cluster.limit_groupsize( opts.limit_cluster_size );
  }

  // now go thru the clusters and write out c.XXX.NNN tags
  toolbox::ClusterBase::ClusterList const & clusterlist=cluster.clusterlist();

  _impl_tr.Info << " clustering: " << clusterlist.size() << " clusters found. ";
  if ( opts.limit_cluster_size ) {
    _impl_tr.Info << " cluster size limited to "
	    << opts.limit_cluster_size << "\n";
  }
  _impl_tr.Info << std::endl;
  // one entry per decoys corresponding to the running number in the input data which is also used as index in the clusters
  // the "kept_tags" entry will remain "" for removed structures (limit_cluster_size).
  // this is slightly memory intensive... the alternative is that for each decoy in input_decoys.begin...end we have
  // to find the corresponding entry in the clusterlist.
  utility::vector1< std::string > kept_tags( n_decoys, "" );
  utility::vector1< std::string > kept_orig_tags( n_decoys, "" );

  _impl_tr.Info << " generate new tags... " << std::endl;
  for ( core::Size i=1; i<=clusterlist.size(); i++ ) {
    for ( core::Size j=1; j<=clusterlist[i].size(); j++ ) {
      using namespace ObjexxFCL;
      //form tags of type: c.<cluster_number>.<decoy_in_group_nr>
      kept_tags[ clusterlist[i][j-1] ] = "c."+lead_zero_string_of( i-1,4 )+"."+lead_zero_string_of(  j-1,3 );
    }
  }

  _impl_tr.Info << "copy remaining structures to output... " << std::endl;
  utility::vector1< std::string >::const_iterator tags_it = kept_tags.begin();
  core::Size ct( 1 );
  for ( SilentStructIterator it=input_decoys_begin; it!=input_decoys_end; ++it, ++tags_it, ++ct )	{

    //if tag is "" this decoy has been filtered out by limit_cluster_size
    if ( tags_it->size() == 0 ) {
      _impl_tr.Info << "removed decoy " << it->decoy_tag() << " with score " <<  it->get_energy( "score" ) << "\n";
      continue;
    }

    //this decoy will be transfered into output
    _impl_tr.Info << "keep decoy "<< it->decoy_tag() << " with score " << it->get_energy( "score" ) << " as " << *tags_it << "\n";
    core::io::silent::SilentStructOP new_decoy = it->clone();

    //rename decoy to "c.XXX.NNN" ?
    if ( opts.assign_new_cluster_tag ) {
      new_decoy->add_string_value( "decoy_tag", it->decoy_tag() );
      new_decoy->set_decoy_tag( *tags_it );
    }

    new_structs.push_back( new_decoy );

    //store these for the "print_summary" output --- definitly memory that could be saved
    kept_orig_tags[ ct ] = it->decoy_tag();
  }

  _impl_tr.Info << std::endl;
  cluster.print_summary( kept_orig_tags, CA_set.all_energies() );
}


}
}


#endif
