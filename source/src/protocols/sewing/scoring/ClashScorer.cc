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


#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.ClashScorer");

ClashScorer::ClashScorer(){}



core::Real
ClashScorer::score(
	AssemblyCOP assembly
) {

	using namespace basic::options;
	// using namespace basic::options::OptionKeys;

	core::Real score = 0.0;

	core::scoring::AtomVDW const & atom_vdw( core::scoring::ScoringManager::get_instance()->get_AtomVDW( core::chemical::CENTROID ));
	//TR << AtomVDW
	ModelConstIterator<SewSegment> it1=assembly->assembly_begin();
	ModelConstIterator<SewSegment> end=assembly->assembly_end();

	//Foreach atom in the assembly (only backbone atoms)
	for ( ; it1 != end; ++it1 ) {

		//Get the xyz coordinates and the vanderwaals size for the current atom
		numeric::xyzVector<core::Real> i_xyz ( it1.atom()->coords_ );
		Size const i_type( it1.atom()->atomno_ );


		utility::vector1< core::Real > const & i_atom_vdw( atom_vdw( i_type ) );

		//  if(TR.Debug.visible()){
		//   TR << "it1.segment()->dssp_: " << it1.segment()->dssp_ << std::endl;
		//  }
		//  if(TR.Debug.visible()){
		////   TR.Debug << "it1.atom()->atomno_ : " << it1.atom()->atomno_ << " " ;
		////   if (it1.atom()->atomno_ == 1){ TR.Debug << "element: N" << std::endl; }
		////   else if (it1.atom()->atomno_ == 2){ TR.Debug << "element: CA" << std::endl; }
		////   else if (it1.atom()->atomno_ == 3){ TR.Debug << "element: C" << std::endl; }
		////   else if (it1.atom()->atomno_ == 4){ TR.Debug << "element: O" << std::endl; }
		//
		//  // TR.Debug << "i_type( it1.atom()->atomno_ ): " << i_type << std::endl;
		//  // TR.Debug << "i_atom_vdw[ i_type ]: " << i_atom_vdw[ i_type ] << std::endl;
		////   TR.Debug << "sqrt of i_atom_vdw[ i_type ]: " << sqrt(i_atom_vdw[ i_type ]) << std::endl;
		////   TR.Debug << std::endl;
		//  }


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
			} // Hi Tim! How come this case is even possible? (Doonam)



			//Get the xyz coordinates and the vanderwaals size for the current atom
			numeric::xyzVector<core::Real> j_xyz ( it2.atom()->coords_ );
			Size const j_type( it2.atom()->atomno_ );

			//This bump distance seems weird, double check this...

			// original bump_dsq
			//core::Real const bump_dsq = i_atom_vdw[ j_type ] ;

			// offset to bump_dsq
			core::Real offset_bump_dsq_ = option[OptionKeys::sewing::offset_bump_dsq];

			//   if(TR.Debug.visible()){
			//    //TR.Debug << "offset_bump_dsq_: " << offset_bump_dsq_ << std::endl;
			//   }

			// bump_dsq => (bump_distance)^2
			//core::Real const bump_dsq = i_atom_vdw[ j_type ] + offset_bump_dsq_;
			core::Real const bump_dsq = i_atom_vdw[ j_type ];

			//utility::vector1< core::Real > const & j_atom_vdw( atom_vdw( j_type ) );

			/*
			if(TR.Debug.visible()){
			TR.Debug << "it1.atom()->atomno_ : " << it1.atom()->atomno_ << " " ;
			if (it1.atom()->atomno_ == 1){ TR.Debug << "element: N" << std::endl; }
			else if (it1.atom()->atomno_ == 2){ TR.Debug << "element: CA" << std::endl; }
			else if (it1.atom()->atomno_ == 3){ TR.Debug << "element: C" << std::endl; }
			else if (it1.atom()->atomno_ == 4){ TR.Debug << "element: O" << std::endl; }

			TR.Debug << "it2.atom()->atomno_ : " << it2.atom()->atomno_ << " " ;
			// TR.Debug << "j_type( it2.atom()->atomno_ ): " << j_type << std::endl;
			if (it2.atom()->atomno_ == 1){ TR.Debug << "element: N" << std::endl; }
			else if (it2.atom()->atomno_ == 2){ TR.Debug << "element: CA" << std::endl; }
			else if (it2.atom()->atomno_ == 3){ TR.Debug << "element: C" << std::endl; }
			else if (it2.atom()->atomno_ == 4){ TR.Debug << "element: O" << std::endl; }

			if ( ((it1.atom()->atomno_ == 1) || (it1.atom()->atomno_ == 4) ) && ( (it2.atom()->atomno_ == 1) || (it2.atom()->atomno_ == 4) ) ){
			TR.Debug << "N, O" << std::endl;
			}

			// TR.Debug << "bump_dsq (=i_atom_vdw[ j_type ]): " << bump_dsq << std::endl;
			//TR.Debug << "sqrt of bump_dsq: " << sqrt(bump_dsq) << std::endl;
			}


			if(TR.Debug.visible()){
			//TR.Debug << "i_xyz.distance_squared( j_xyz ): " << i_xyz.distance_squared( j_xyz ) << std::endl;
			//TR.Debug << "sqrt of i_xyz.distance_squared( j_xyz ): " << sqrt(i_xyz.distance_squared( j_xyz )) << std::endl;
			//TR.Debug << "clash = bump_dsq - i_xyz.distance_squared( j_xyz ): " << clash << std::endl;
			//TR.Debug << "sqrt of clash: " << sqrt(clash) << std::endl;
			TR.Debug << std::endl;
			}
			*/

			core::Real const clash = bump_dsq - i_xyz.distance_squared( j_xyz ) + offset_bump_dsq_;

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

				if ( TR.Debug.visible() ) {
					//TR << "it1.segment()->dssp_: " << it1.segment()->dssp_ << std::endl;
					//TR << "it2.segment()->dssp_: " << it2.segment()->dssp_ << std::endl;
					if ( it1.segment()->dssp_ == 'E' && it2.segment()->dssp_ == 'E' ) {
						TR.Debug << "sqrt of i_xyz.distance_squared( j_xyz ): " << sqrt(i_xyz.distance_squared( j_xyz )) << std::endl;
						TR.Debug << "sqrt of bump_dsq: " << sqrt(bump_dsq) << std::endl;
						//      TR.Debug << "clash = bump_dsq - i_xyz.distance_squared( j_xyz ): " << clash << std::endl;
						TR.Debug << "it1.atom()->atomno_ : " << it1.atom()->atomno_ << " " ;
						TR.Debug << "it2.atom()->atomno_ : " << it2.atom()->atomno_ << " " << std::endl;
					}
				} //if ( TR.Debug.visible() )

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
	// if(TR.Debug.visible()){
	//  utility_exit_with_message("quit now!"); // temporary for devel
	//  //exit(1);
	// }
	return score;
}


} //scoring namespace
} //sewing namespace
} //protocols namespace
