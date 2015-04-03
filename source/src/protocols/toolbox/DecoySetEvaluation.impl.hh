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

#ifndef INCLUDED_protocols_toolbox_DecoySetEvaluation_impl_hh
#define INCLUDED_protocols_toolbox_DecoySetEvaluation_impl_hh

#include <protocols/toolbox/DecoySetEvaluation.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {


template< typename SilentStructIterator >
void DecoySetEvaluation::push_back_CA_xyz_from_silent_file( Size n_decoys_in, SilentStructIterator begin, SilentStructIterator end, bool store_energies ) {
static thread_local basic::Tracer _impl_tr( "protocols.toolbox.DecoySetEvaluation" );

  Size n_new_decoys( n_decoys_in );

  if ( begin == end ) return;
  core::pose::Pose pose;
  begin->fill_pose( pose);
  Size pos = 1;
  for ( ; pos <= pose.total_residue(); ++pos ) {
    if ( !pose.residue( pos ).is_protein() ) break;
  }
  --pos;
  if ( pos < pose.total_residue() ) {
    _impl_tr.Warning << "Found no CA in sequnce at position " << pos+1 << " Will stop pool before this position... " << std::endl;
    set_n_atom( pos );
  }

  //if you read more than one sfd, don't switch the behaviour with energy storage around
  runtime_assert( store_energies == store_energies_ || n_decoys() == 0 );
  store_energies_ = store_energies;

  if ( n_decoys_max() <= n_decoys() + n_new_decoys ) {
    reserve( n_decoys() + n_new_decoys );
  }

  for ( SilentStructIterator it=begin; it!=end && n_new_decoys>0; ++it, --n_new_decoys )	{
    push_back_CA_xyz( it->get_CA_xyz(), n_atoms() > 0 ? n_atoms() : it->nres() );
    if ( store_energies_ ) {
      all_energies_.push_back( it->get_energy( "score" ) );
    }
  }
}

}
}

#endif
