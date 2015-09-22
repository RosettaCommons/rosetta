// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .hh file for enzdes sequence recovery cache
/// @brief
/// @author sinibjelic@gmail.com

//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.hh>

//package headers

//project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <set>
#include <string>

#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static THREAD_LOCAL basic::Tracer TR( "protocols.enzdes.EnzdesSeqRecoveryCache" );

EnzdesSeqRecoveryCache::EnzdesSeqRecoveryCache() {
	sequence_.clear();
	designable_residues_.clear();
}

//copy constructor
EnzdesSeqRecoveryCache::EnzdesSeqRecoveryCache( EnzdesSeqRecoveryCache const & other ) :
	ReferenceCount( other ) {
	sequence_ = other.sequence_;
	designable_residues_ = other.designable_residues_;
}

EnzdesSeqRecoveryCache::~EnzdesSeqRecoveryCache() {}

void
EnzdesSeqRecoveryCache::set_sequence(
	core::pose::Pose & native_pose
) {
	for ( core::Size jj=1; jj <= native_pose.total_residue(); ++jj ) {
		sequence_[jj] = native_pose.residue( jj ).name1();
	}
}

std::map< core::Size, char >
EnzdesSeqRecoveryCache::get_sequence() {
	return sequence_;
}

void
EnzdesSeqRecoveryCache::set_designable_residues(
	std::set< core::Size > des_res
) {
	std::set< core::Size >::const_iterator it;
	for ( it = des_res.begin(); it != des_res.end(); ++it ) {
		//sequnce_ keeps indirectly track of what residues are wt
		//as insertions and deletions are never added and always deleted
		if ( sequence_.find(*it) != sequence_.end() ) designable_residues_.insert( *it );
	}
}

std::set< core::Size >
EnzdesSeqRecoveryCache::get_designable_residues() {
	return  designable_residues_;
}

core::Real
EnzdesSeqRecoveryCache::sequence_recovery(
	core::pose::Pose const & designed_pose
) const
{
	core::Size n_residues_recovered(0);
	core::Size n_residues_total(0);

	//check if the container is full/empty
	if ( !designable_residues_.empty() ) {
		std::set< core::Size >::const_iterator it;
		for ( it = designable_residues_.begin(); it != designable_residues_.end(); ++it ) {
			if ( sequence_.find(*it)->second  == designed_pose.residue(*it).name1() )  {
				++n_residues_recovered;
			}
		}
		n_residues_total = designable_residues_.size();
		return ( static_cast< core::Real > ( n_residues_recovered) / n_residues_total );
	}

	//No residues have been set to designable or there is something seriously wrong
	//No sequence change => 1.0
	return 1.0;
}

void
EnzdesSeqRecoveryCache::remap_residues(
	core::id::SequenceMapping const & smap
){
	std::map< core::Size, char > remap_sequence;
	std::set< core::Size >remap_designable_residues;

	//smap( old res number ) = new res number
	for ( core::Size it=1; it <= smap.mapping().size(); ++it ) {
		//remap sequence_
		if ( smap[it] != 0 && sequence_.find( it ) != sequence_.end() ) {
			remap_sequence[ smap[it] ] = sequence_.find( it ) -> second;
		}
		//remap designable_residues_
		if (  smap[it] != 0 && designable_residues_.find( it ) != designable_residues_.end() ) {
			remap_designable_residues.insert( smap[it] );
		}
	}//smap for loop

	//modify wt sequence accordingly
	sequence_ = remap_sequence;

	//assign new modified designable_residues_ from remap_designable_residues.
	//designable_residues_ are cleared automatically by = operator
	designable_residues_ = remap_designable_residues;
}


} //match_enzdes_util
} //toolbox
} //protocols
