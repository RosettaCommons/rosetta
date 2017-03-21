// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/matdes/MotifHitsRotamersOperation.hh
///
/// @brief
/// @author Jorge Fallas (jaf18@uw.edu)

#include <protocols/matdes/MotifHitsRotamersOperation.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <utility/io/ozstream.hh>
#include <numeric/xyzTransform.hh>

namespace protocols {
namespace matdes {


static basic::Tracer TR("protocols.matdes.MotifHitsRotamersOperation");

//Constructor

MotifHitsRotamersOperation::MotifHitsRotamersOperation(
	core::scoring::motif::MotifHits const &  motif_hits
):
	RotamerSetOperation(),
	motif_hits_(motif_hits),
	total_count_(0)
{
}

core::pack::rotamer_set::RotamerSetOperationOP
MotifHitsRotamersOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new MotifHitsRotamersOperation( *this ));
}


/// @brief Helper function, combines existing's metadata with conformer's conformation.
core::conformation::ResidueOP
dup_residue(
	core::conformation::Residue const & existing,
	core::conformation::Residue const & conformer
)
{
	// This is bad:  fields like seqpos, chain, etc. don't match existing residue!
	//conformation::ResidueOP newrsd = rotamers_[i]->clone();

	// Could start by cloning either one, but I think people are more likely to introduce
	// new metadata than new conformational data, so I'll let clone() copy the metadata.
	//conformation::ResidueOP newrsd = existing.clone();
	//newrsd->atoms() = conformer.atoms();
	//newrsd->chi() = conformer.chi();
	//newrsd->mainchain_torsions() = conformer.mainchain_torsions();
	//newrsd->actcoord() = conformer.actcoord();

	// The above is also bad:  existing may not be the same residue type as conformer!
	core::conformation::ResidueOP newrsd = conformer.clone();
	newrsd->chain( existing.chain() );
	newrsd->seqpos( existing.seqpos() );
	newrsd->copy_residue_connections_from( existing ); // this is probably not good enough if residue types diverge more than protonation state...

	return newrsd;
}


void
MotifHitsRotamersOperation::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & /*sfxn*/,
	core::pack::task::PackerTask const & ptask,
	utility::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
)
{
	// Size const seqnum = (Size) rotamer_set.resid();
	// debug_assert( seqnum <= ptask.total_residue() );
	// core::pack::task::ResidueLevelTask const & rtask = ptask.residue_task(seqnum);
	// for(Size i = 1; i <= poses_.size(); ++i) {
	//  core::pose::Pose const & ubr_pose = *(poses_[i]);
	//  if(seqnum > ubr_pose.total_residue()) continue;
	//  core::chemical::ResidueType const & restype = ubr_pose.residue_type(seqnum);
	//  bool type_is_allowed = false;
	//  for(core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter j = rtask.allowed_residue_types_begin(),
	//   j_end = rtask.allowed_residue_types_end(); j != j_end; ++j)
	//  {
	//   if( restype.name() == (**j).name() ) {
	//    type_is_allowed = true;
	//    break;
	//   }
	//  }
	//  if( type_is_allowed ) {
	//   TR.Debug << "Adding 'unbound' rotamer at position " << seqnum << std::endl;
	//   conformation::ResidueOP newrsd = dup_residue( pose.residue(seqnum), ubr_pose.residue(seqnum) );
	//   newrsd->place( pose.residue(seqnum), pose.conformation() );
	//   rotamer_set.add_rotamer( *newrsd );
	//  } else {
	//   TR.Debug << "Residue names do not match. Skipping 'unbound' rotamer at position " << seqnum << std::endl;
	//  }
	// }
	//utility::io::ozstream out("test_residues.pdb", std::ios::app) ;
	//core::Size atoms = 1  ;
	using namespace core::pose::symmetry ;
	using namespace core::scoring::motif  ;
	Size const seqnum = (Size) rotamer_set.resid();
	debug_assert( seqnum <= ptask.total_residue() );
	core::conformation::symmetry::SymmetryInfoCOP syminfo = 0 ;
	if ( is_symmetric(pose) ) syminfo = symmetry_info(pose) ;
	int count = 0 ;
	std::set<core::Size> sym_pos ;
	for ( core::Size i = 1 ; i <= pose.size() ; ++i ) {
		if ( i == seqnum ) sym_pos.insert(i) ;
		if ( syminfo && syminfo->bb_follows(i) == seqnum ) sym_pos.insert(i) ;
	}
	//std::cout <<  " set of matching residues"  ;
	//for (core::Size const i : sym_pos) std::cout << " " << i ;
	//std::cout << std::endl  ;
	for ( MotifHit const & hit : motif_hits_ ) {
		core::conformation::ResidueOP res = NULL ;
		core::Size motif_seqnum = 0   ;
		bool is_bb_frame=false  ;
		if ( sym_pos.find( (core::Size) hit.residue1) != sym_pos.end() )  {
			res = dup_residue(pose.residue(seqnum), hit.mpose().residue(1)) ;
			is_bb_frame = ( hit.motif.type1() == core::scoring::motif::RM_BB) ;
			motif_seqnum = hit.residue1 ;
		}
		if ( sym_pos.find( (core::Size) hit.residue2) != sym_pos.end() ) {
			res = dup_residue(pose.residue(seqnum), hit.mpose().residue(2)) ;
			is_bb_frame = ( hit.motif.type2() == core::scoring::motif::RM_BB) ;
			motif_seqnum = hit.residue2  ;
		}
		if ( res && is_bb_frame ) {
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				numeric::xyzTransform<core::Real> from(
					pose.residue(motif_seqnum).xyz("N"),
					pose.residue(motif_seqnum).xyz("CA"),
					pose.residue(motif_seqnum).xyz("C")) ;
				numeric::xyzTransform<core::Real> to(
					pose.residue(seqnum).xyz("N"),
					pose.residue(seqnum).xyz("CA"),
					pose.residue(seqnum).xyz("C")) ;
				numeric::xyzTransform<core::Real> X(to*~from);
				for ( core::Size ia = 1; ia <= res->natoms(); ++ia ) {
					if ( res->atom_is_backbone(ia) ) {
						res->set_xyz(ia, pose.residue(seqnum).xyz(ia)) ;
					} else res->set_xyz(ia, X*res->xyz(ia)) ;
				}
			}
			// out << "MODEL" << std::endl ;
			// core::io::pdb::dump_pdb_residue(*res,atoms,out)    ;
			// out << "ENDMDL" << std::endl ;
			res->place( pose.residue(seqnum), pose.conformation() );
			rotamer_set.add_rotamer( *res );
			++count  ;
			++total_count_ ;
		}
	}
	if ( count ) {
		TR << "added " << count << " motif rotamers " << "at position "  <<  seqnum  <<std::endl ;
		TR << "added " << total_count_ << " motif rotamers " << "so far "  << std::endl ;
	}
	//out.close()   ;
}


} // namespace matdes
} // namespace protocols
