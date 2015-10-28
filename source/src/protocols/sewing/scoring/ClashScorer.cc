// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ClashScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/ClashScorer.hh>

//Core headers
#include <core/types.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/conformation/Residue.hh>

//Utility headers
#include <basic/Tracer.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.ClashScorer");

ClashScorer::ClashScorer(){}

core::Real
ClashScorer::score(
	AssemblyCOP assembly
) {

	core::Real score = 0.0;

	core::scoring::AtomVDW const & atom_vdw( core::scoring::ScoringManager::get_instance()->get_AtomVDW( core::chemical::CENTROID ));

	ModelConstIterator<SewSegment> it1=assembly->assembly_begin();
	ModelConstIterator<SewSegment> end=assembly->assembly_end();

	//Foreach atom in the assembly (only backbone atoms)
	for ( ; it1 != end; ++it1 ) {

		//Get the xyz coordinates and the vanderwaals size for the current atom
		numeric::xyzVector<core::Real> i_xyz ( it1.atom()->coords_ );
		Size const i_type( it1.atom()->atomno_ );
		utility::vector1< core::Real > const & i_atom_vdw( atom_vdw( i_type ) );

		ModelConstIterator<SewSegment> it2 = it1;
		++it2;
		for ( ; it2 != end; ++it2 ) {

			//Don't score clashes within the same model, they are either natives or
			//were generated during chimerization of segments. Either way, we'll give
			//refinment a chance to fix these.
			if ( it1.segment()->model_id_ == it2.segment()->model_id_ ) {
				continue;
			}

			//If these are two adjacent residues that don't share a model id due to chimerization
			//then it isn't a clash
			//   if(it1.segment()->model_id_ == CHIMERA_SEGMENT &&
			//    std::find(it1.segment()->parent_segments_.begin(), it1.segment()->parent_segments_.end(), std::make_pair(it2.segment()->model_id_, it2.segment()->segment_id_) ) != it1.segment()->parent_segments_.end()) {
			//    continue;
			//   }
			//   if(it2.segment()->model_id_ == CHIMERA_SEGMENT &&
			//    std::find(it2.segment()->parent_segments_.begin(), it2.segment()->parent_segments_.end(), std::make_pair(it1.segment()->model_id_, it1.segment()->segment_id_) ) != it2.segment()->parent_segments_.end()) {
			//    continue;
			//   }
			if ( (it1.segment()->model_id_ <= 0 || it2.segment()->model_id_ <= 0)
					&& it1.residue()->resnum_ == it2.residue()->resnum_ - 1 ) {
				continue;
			}

			//Get the xyz coordinates and the vanderwaals size for the current atom
			numeric::xyzVector<core::Real> j_xyz ( it2.atom()->coords_ );
			Size const j_type( it2.atom()->atomno_ );

			//This bump distance seems weird, double check this...
			core::Real const bump_dsq = i_atom_vdw[ j_type ] ;
			core::Real const clash = bump_dsq - i_xyz.distance_squared( j_xyz );
			if ( clash > 0.0 ) {
				score += 1;
				if ( TR.Debug.visible() ) {
					int model_id_1 = it1.segment()->model_id_;
					core::Size resnum_1 = it1.residue()->resnum_;
					int model_id_2 = it2.segment()->model_id_;
					core::Size resnum_2 = it2.residue()->resnum_;

					core::Size num_1 = assembly->pose_num(model_id_1, resnum_1);
					core::Size num_2 = assembly->pose_num(model_id_2, resnum_2);

					TR.Debug << "Assembly clash "
						<< model_id_1 << ":" << resnum_1 << " - "
						<< model_id_2 << ":" << resnum_2 << std::endl
						<< num_1 << std::endl
						<< num_2 << std::endl;
				}
			} //clash
		} //it2
	} //it1

	core::pose::PoseCOP partner_pose = assembly->get_partner();
	if ( partner_pose ) {
		ModelConstIterator<SewSegment> assembly_it=assembly->assembly_begin();
		ModelConstIterator<SewSegment> assembly_it_end=assembly->assembly_end();

		//Foreach atom in the assembly (only backbone atoms)
		for ( ; assembly_it != assembly_it_end; ++assembly_it ) {

			//Get the xyz coordinates and the vanderwaals size for the current atom
			numeric::xyzVector<core::Real> i_xyz ( assembly_it.atom()->coords_ );
			Size const i_type( assembly_it.atom()->atomno_ );
			utility::vector1< core::Real > const & i_atom_vdw( atom_vdw( i_type ) );

			for ( core::Size i=1; i<=partner_pose->total_residue(); ++i ) {
				core::conformation::Residue const & partner_res( partner_pose->residue( i ) );
				core::chemical::AtomIndices const & bb_atoms = partner_pose->residue(i).all_bb_atoms();

				//Foreach bb-atom in the current residue
				for ( core::Size j=1; j <= bb_atoms.size(); ++j ) {
					numeric::xyzVector< core::Real > const & j_xyz =
						partner_pose->residue(i).atom(bb_atoms[j]).xyz();
					core::Real const bump_dsq( i_atom_vdw[ partner_res.atom_type_index(j) ] );
					core::Real const clash( bump_dsq - i_xyz.distance_squared( j_xyz ) );
					if ( clash > 0.0 ) {
						score += 1;
						if ( TR.Debug.visible() ) {
							int model_id = assembly_it.segment()->model_id_;
							core::Size resnum = assembly_it.residue()->resnum_;

							TR.Debug << "Assembly clash with partner residue"
								<< model_id << ":" << resnum << " - "
								<< i << ":" << j << std::endl;
						}
					}
				}
			}
		}
	}

	return score;
}


} //scoring namespace
} //sewing namespace
} //protocols namespace
