// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/MutateFrameworkForCluster.hh
/// @brief Mutates Framework regions after insertion of a particular cluster
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/design/NativeAntibodySeq.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/design/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/Tracer.hh>

#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <utility/vector1.srlz.hh>

#endif // SERIALIZATION

static basic::Tracer TR("protocols.antibody.NativeAntibodySeq");

namespace protocols {
namespace antibody {
namespace design {
using namespace core::chemical;
using namespace basic::datacache;
using basic::datacache::DataCache_CacheableData;

NativeAntibodySeq::NativeAntibodySeq(const core::pose::Pose &pose,
	protocols::antibody::AntibodyInfo const & ab_info) :
	CacheableData()
{
	set_sequence(pose, ab_info);
}

NativeAntibodySeq::NativeAntibodySeq(NativeAntibodySeq const &src):
	CacheableData(),
	seq_(src.seq_),
	cdr_seq_(src.cdr_seq_)
{

}

NativeAntibodySeq::~NativeAntibodySeq() = default;

basic::datacache::CacheableDataOP
NativeAntibodySeq::clone() const {
	return utility::pointer::make_shared< NativeAntibodySeq >(*this);
}

void
NativeAntibodySeq::set_sequence(const core::pose::Pose &pose, AntibodyInfo const & ab_info ) {
	seq_.clear();
	cdr_seq_.clear();
	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		AA res = pose.aa(i);
		//PDBNumbering info;
		//info.resnum = pose.pdb_info()->number(i);
		//info.chain = pose.pdb_info()->chain(i);
		//info.icode = pose.pdb_info()->icode(i);

		if ( ab_info.get_region_of_residue(pose, i, false /* count CDR4 as framework */) != cdr_region ) {
			seq_[ pose.pdb_info()->pose2pdb( i )] = res;

		}
	}

	//Setup CDR Regions
	for ( auto const & cdr : ab_info.get_all_cdrs_present( true /* include CDR4 */) ) {
		set_from_cdr( pose, ab_info, cdr);
	}
}

void
NativeAntibodySeq::set_to_pose(core::pose::Pose &pose) {
	// Do NOT attempt to put 'this' into an owning pointer. We don't know where it's allocaed (e.g. on the stack).
	// Make a copy of this object and use that, instead.
	pose.data().set(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ, clone());
}

void
NativeAntibodySeq::set_from_cdr(const core::pose::Pose &pose, AntibodyInfo const & ab_info, CDRNameEnum cdr) {

	utility::vector1<AA> s;
	core::Size cdr_start_resnum = ab_info.get_CDR_start( cdr, pose);
	for ( core::Size cdr_pose_index = cdr_start_resnum; cdr_pose_index <= cdr_start_resnum + ab_info.get_CDR_length( cdr, pose ) - 1; ++cdr_pose_index ) {
		AA res = pose.aa( cdr_pose_index );
		s.push_back( res );
	}

	cdr_seq_[ cdr ] = s;


}

std::map< std::string, core::chemical::AA> const &
NativeAntibodySeq::get_full_sequence() const {
	return seq_;
}

std::map< CDRNameEnum, utility::vector1< core::chemical::AA>> const &
NativeAntibodySeq::get_cdr_sequence() const {
	return cdr_seq_;
}

std::string
NativeAntibodySeq::get_native_sequence_matching_current_length(const core::pose::Pose & pose, AntibodyInfo const & ab_info) const {

	//Use any sequence that is set here.  If not, use the current one.
	utility::vector1<AA> local_seq;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		std::string info = pose.pdb_info()->pose2pdb( i );
		if ( seq_.find( info )!= seq_.end() ) {
			AA res = seq_.find( info )->second;
			local_seq.push_back( res );
		} else {
			AA res = pose.aa( i );
			local_seq.push_back( res );
		}
	}

	//Go through each of the CDRs.  Update the final sequence with what is stored here.
	for ( auto const & cdr : ab_info.get_all_cdrs_present() ) {
		core::Size cdr_start_resnum = ab_info.get_CDR_start( cdr, pose );

		if ( cdr_seq_.find( cdr ) != cdr_seq_.end() ) continue;

		utility::vector1< AA > local_cdr_seq = cdr_seq_.find( cdr )->second;
		core::Size local_index = 1;
		for ( core::Size cdr_pose_index = cdr_start_resnum; cdr_pose_index <= cdr_start_resnum + ab_info.get_CDR_length( cdr, pose ) - 1; ++cdr_pose_index ) {
			local_seq[ cdr_pose_index ] = local_cdr_seq[ local_index ];
			local_index +=1;
		}
	}

	//Convert vector to string and return
	std::string return_seq = "";
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		return_seq += core::chemical::oneletter_code_from_aa( local_seq[ i ]);
	}
	return return_seq;
}

#ifdef    SERIALIZATION
NativeAntibodySeq::NativeAntibodySeq():
	basic::datacache::CacheableData(){

}
template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION



} //design
} //antibody
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::antibody::design::NativeAntibodySeq::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( seq_ ) );
	arc( CEREAL_NVP( cdr_seq_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::antibody::design::NativeAntibodySeq::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( seq_ );
	arc( cdr_seq_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::antibody::design::NativeAntibodySeq );
CEREAL_REGISTER_TYPE( protocols::antibody::design::NativeAntibodySeq )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_antibody_design_NativeAntibodySeq )
#endif // SERIALIZATION


