// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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


static THREAD_LOCAL basic::Tracer TR("protocols.antibody.NativeAntibodySeq");

namespace protocols {
namespace antibody {
namespace design {
using namespace core::chemical;
using namespace basic::datacache;
using basic::datacache::DataCache_CacheableData;

NativeAntibodySeq::NativeAntibodySeq(const core::pose::Pose &pose,
	protocols::antibody::AntibodyInfoCOP ab_info) :
	CacheableData(),
	ab_info_(ab_info)
{
	ab_info_ = ab_info->clone();
	set_sequence(pose);
}

NativeAntibodySeq::NativeAntibodySeq(NativeAntibodySeq const &src):
	CacheableData(),
	ab_info_(src.ab_info_),
	seq_(src.seq_),
	cdr_seq_(src.cdr_seq_)
{}

NativeAntibodySeq::~NativeAntibodySeq() {}

basic::datacache::CacheableDataOP
NativeAntibodySeq::clone() const {
	return CacheableDataOP( new NativeAntibodySeq(*this) );
}

void
NativeAntibodySeq::set_sequence(const core::pose::Pose &pose) {
	seq_.clear();
	cdr_seq_.clear();
	TR << "Setting full pose sequence" << std::endl;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		AA res = pose.aa(i);
		//PDBNumbering info;
		//info.resnum = pose.pdb_info()->number(i);
		//info.chain = pose.pdb_info()->chain(i);
		//info.icode = pose.pdb_info()->icode(i);

		if ( ab_info_->get_region_of_residue(pose, i) != cdr_region ) {
			seq_[ pose.pdb_info()->pose2pdb( i )] = res;

		}
	}

	//Setup CDR Regions
	for ( core::Size i = 1; i <= core::Size(ab_info_->get_total_num_CDRs()); ++i ) {
		CDRNameEnum cdr = static_cast< CDRNameEnum >( i );
		set_from_cdr( pose, cdr );
	}
}

void
NativeAntibodySeq::set_to_pose(core::pose::Pose &pose) {
	// Do NOT attempt to put 'this' into an owning pointer. We don't know where it's allocaed (e.g. on the stack).
	// Make a copy of this object and use that, instead.
	pose.data().set(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ, clone());
}

void
NativeAntibodySeq::set_from_cdr(const core::pose::Pose &pose, CDRNameEnum cdr) {

	utility::vector1<AA> s;
	core::Size cdr_start_resnum = ab_info_->get_CDR_start( cdr, pose);
	for ( core::Size cdr_pose_index = cdr_start_resnum; cdr_pose_index <= cdr_start_resnum + ab_info_->get_CDR_length( cdr, pose ) - 1; ++cdr_pose_index ) {
		AA res = pose.aa( cdr_pose_index );
		s.push_back( res );
	}

	cdr_seq_[ cdr ] = s;


}

std::string
NativeAntibodySeq::get_sequence(const core::pose::Pose & pose) const {

	//Use any sequence that is set here.  If not, use the current one.
	utility::vector1<AA> local_seq;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
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
	for ( core::Size i = 1; i <= core::Size( ab_info_->get_total_num_CDRs()); ++i ) {
		CDRNameEnum cdr = static_cast< CDRNameEnum >( i );
		core::Size cdr_start_resnum = ab_info_->get_CDR_start( cdr, pose );

		if ( cdr_seq_.find( cdr ) != cdr_seq_.end() ) continue;

		utility::vector1< AA > local_cdr_seq = cdr_seq_.find( cdr )->second;
		core::Size local_index = 1;
		for ( core::Size cdr_pose_index = cdr_start_resnum; cdr_pose_index <= cdr_start_resnum + ab_info_->get_CDR_length( cdr, pose ) - 1; ++cdr_pose_index ) {
			local_seq[ cdr_pose_index ] = local_cdr_seq[ local_index ];
			local_index +=1;
		}
	}

	//Convert vector to string and return
	std::string return_seq = "";
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		return_seq += core::chemical::oneletter_code_from_aa( local_seq[ i ]);
	}
	return return_seq;
}



}
}
}


