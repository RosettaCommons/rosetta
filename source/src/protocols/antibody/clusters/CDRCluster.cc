// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRCluster.cc
/// @brief Simple class to hold and access CDRCluster information of a region in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/CDRCluster.hh>
#include <core/pose/PDBInfo.hh>
#include <math.h>
#include <numeric/NumericTraits.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace antibody {
namespace clusters {
using namespace protocols::antibody;


CDRCluster::CDRCluster(core::pose::Pose const & pose, CDRNameEnum const cdr, core::Size const cdr_length,
	CDRClusterEnum const cluster, core::Size const start, core::Real const distance,
	bool cis_trans_match /*true*/):
	utility::pointer::ReferenceCount(),
	cdr_(cdr),
	cluster_(cluster),
	distance_(distance),
	start_(start),
	length_(cdr_length),
	cis_trans_match_(cis_trans_match)
{
	end_ = start+cdr_length-1;
	normalized_distance_ = distance/(cdr_length*2);
	set_pdb_numbering(pose, start, end_);
}

CDRCluster::CDRCluster(const CDRCluster& src):
	utility::pointer::ReferenceCount(src),
	cdr_(src.cdr_),
	cluster_(src.cluster_),
	distance_(src.distance_),
	normalized_distance_(src.normalized_distance_),
	pdb_start_(src.pdb_start_),
	pdb_end_(src.pdb_end_),
	pdb_start_insertion_code_(src.pdb_start_insertion_code_),
	pdb_end_insertion_code_(src.pdb_end_insertion_code_),
	start_(src.start_),
	end_(src.end_),
	length_(src.length_),
	chain_(src.chain_),
	cis_trans_match_(src.cis_trans_match_)

{

}

CDRClusterOP
CDRCluster::clone() const {
	return CDRClusterOP( new CDRCluster(*this) );
}

CDRCluster::~CDRCluster(){}

core::Real
CDRCluster::normalized_distance_in_degrees() const {

	core::Real result = acos(1 - (normalized_distance_/2)) *(180/numeric::NumericTraits<core::Real>::pi());
	return result;
}

void
CDRCluster::set_pdb_numbering(core::pose::Pose const & pose, core::Size start, core::Size end){
	pdb_start_ = pose.pdb_info()->number(start);
	pdb_end_ = pose.pdb_info()->number(end);

	pdb_start_insertion_code_ = pose.pdb_info()->icode(start);
	pdb_end_insertion_code_ = pose.pdb_info()->icode(end);
	chain_ = pose.pdb_info()->chain(start);

}

} //design
} //antibody
} //protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::antibody::clusters::CDRCluster::CDRCluster() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::antibody::clusters::CDRCluster::save( Archive & arc ) const {
	arc( CEREAL_NVP( cdr_ ) ); // enum protocols::antibody::CDRNameEnum
	arc( CEREAL_NVP( cluster_ ) ); // enum protocols::antibody::clusters::CDRClusterEnum
	arc( CEREAL_NVP( distance_ ) ); // core::Real
	arc( CEREAL_NVP( normalized_distance_ ) ); // core::Real
	arc( CEREAL_NVP( pdb_start_ ) ); // core::Size
	arc( CEREAL_NVP( pdb_end_ ) ); // core::Size
	arc( CEREAL_NVP( pdb_start_insertion_code_ ) ); // char
	arc( CEREAL_NVP( pdb_end_insertion_code_ ) ); // char
	arc( CEREAL_NVP( start_ ) ); // core::Size
	arc( CEREAL_NVP( end_ ) ); // core::Size
	arc( CEREAL_NVP( length_ ) ); // core::Size
	arc( CEREAL_NVP( chain_ ) ); // char
	arc( CEREAL_NVP( cis_trans_match_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::antibody::clusters::CDRCluster::load( Archive & arc ) {
	arc( cdr_ ); // enum protocols::antibody::CDRNameEnum
	arc( cluster_ ); // enum protocols::antibody::clusters::CDRClusterEnum
	arc( distance_ ); // core::Real
	arc( normalized_distance_ ); // core::Real
	arc( pdb_start_ ); // core::Size
	arc( pdb_end_ ); // core::Size
	arc( pdb_start_insertion_code_ ); // char
	arc( pdb_end_insertion_code_ ); // char
	arc( start_ ); // core::Size
	arc( end_ ); // core::Size
	arc( length_ ); // core::Size
	arc( chain_ ); // char
	arc( cis_trans_match_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::antibody::clusters::CDRCluster );
CEREAL_REGISTER_TYPE( protocols::antibody::clusters::CDRCluster )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_antibody_clusters_CDRCluster )
#endif // SERIALIZATION
