// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#include <protocols/comparative_modeling/hybridize/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>
#include <core/types.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/SOGFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/PDBInfo.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <list>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {
		
using namespace core;
using namespace kinematics;
using namespace core::scoring::constraints;

void setup_centroid_constraints( 
		core::pose::Pose &pose,
		utility::vector1 < core::pose::PoseCOP > templates,
		utility::vector1 < core::Real > template_weights,
		std::string cen_cst_file ) {
	if (cen_cst_file == "AUTO") {
		// automatic constraints
		generate_centroid_constraints( pose, templates, template_weights );
	} else if (!cen_cst_file.empty() && cen_cst_file != "NONE") {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( cen_cst_file, new ConstraintSet, pose );
		pose.constraint_set( constraint_set );  //reset constraints
	}
}

void setup_fullatom_constraints(
		core::pose::Pose &pose,
		utility::vector1 < core::pose::PoseCOP > templates,
		utility::vector1 < core::Real > template_weights,
		std::string cen_cst_file,
		std::string fa_cst_file  ) {

	// use fa if specified, otherwise centroid
	if (fa_cst_file == "AUTO") {
		// automatic fa constraints
		generate_fullatom_constraints( pose, templates, template_weights );
	}	else if (!fa_cst_file.empty() && fa_cst_file != "NONE") {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( fa_cst_file, new ConstraintSet, pose );
		pose.constraint_set( constraint_set );  //reset constraints
	} else if (cen_cst_file == "AUTO") {
		// automatic constraints
		generate_centroid_constraints( pose, templates, template_weights );
	} else if (!cen_cst_file.empty() && cen_cst_file != "NONE") {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( cen_cst_file, new ConstraintSet, pose );
		pose.constraint_set( constraint_set );  //reset constraints
	}
}

void generate_centroid_constraints( 
		core::pose::Pose &pose,
		utility::vector1 < core::pose::PoseCOP > templates,
		utility::vector1 < core::Real > template_weights )
{

	core::Size MINSEQSEP = 8;
	core::Real MAXDIST = 12.0;
	core::Size GAPBUFFER = 3;
	core::Real COORDDEV = 1.0;

	pose.remove_constraints();

	// number of residues
	core::Size nres_tgt = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}
	if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;

	//
	utility::vector1< utility::vector1< core::Real > > tgt_dists(nres_tgt);
	utility::vector1< utility::vector1< core::Real > > tgt_weights(nres_tgt);
	for (int i=1; i<=(int)templates.size(); ++i) {
		utility::vector1< bool > passed_gapcheck(nres_tgt,false);
		for (int j=1; j<(int)templates[i]->total_residue(); ++j ) {
			bool includeme=true;
			for (int k=1; k<=(int)GAPBUFFER && includeme; ++k) {
				if ( j-k < 1 || j+k > (int)templates[i]->total_residue() ) includeme=false;
				else if (templates[i]->pdb_info()->number(j+k) - templates[i]->pdb_info()->number(j) != k ) includeme=false;
				else if (templates[i]->pdb_info()->number(j-k) - templates[i]->pdb_info()->number(j) != -k ) includeme=false;
			}
			passed_gapcheck[j] = includeme;
		}

		for (core::Size j=1; j<templates[i]->total_residue(); ++j ) {
			if (!passed_gapcheck[j]) continue;
			for (core::Size k=j+1; k<templates[i]->total_residue(); ++k ) {
				if (!passed_gapcheck[k]) continue;
				if (templates[i]->pdb_info()->number(k) - templates[i]->pdb_info()->number(j) < (int)MINSEQSEP) continue;
				core::Real dist = templates[i]->residue(j).xyz(2).distance( templates[i]->residue(k).xyz(2) );
				if ( dist <= MAXDIST ) {
					pose.add_constraint(
						new AtomPairConstraint( core::id::AtomID(2,templates[i]->pdb_info()->number(j)),
						                        core::id::AtomID(2,templates[i]->pdb_info()->number(k)), 
							new ScalarWeightedFunc( template_weights[i], new SOGFunc( dist, COORDDEV )  )
						)
					);
				}
			}
		}
	}
}

void generate_fullatom_constraints(
		core::pose::Pose &pose,
		utility::vector1 < core::pose::PoseCOP > templates,
		utility::vector1 < core::Real > template_weights ) {
	//fpd ... for now just use centroid variant
	generate_centroid_constraints( pose, templates, template_weights);
}



bool discontinued_upper(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);
	
	if (seqpos == 1) return true;
	if (!pose.residue_type(seqpos).is_polymer()) return true;
	if (!pose.residue_type(seqpos-1).is_polymer()) return true;
	if (pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos-1).is_protein()) {
		if ( pose.residue(seqpos).xyz("N").distance(pose.residue(seqpos-1).xyz("C")) > N_C_cutoff ) {
			return true;
		}
	}
	return false;
}

bool discontinued_lower(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);
	
	if (seqpos == pose.total_residue()) return true;
	if (!pose.residue_type(seqpos).is_polymer()) return true;
	if (!pose.residue_type(seqpos+1).is_polymer()) return true;
	if ( pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos+1).is_protein()) {
		if ( pose.residue(seqpos).xyz("C").distance(pose.residue(seqpos+1).xyz("N")) > N_C_cutoff ) {
			return true;
		}
	}
	return false;
}

std::list < Size >
downstream_residues_from_jump(core::pose::Pose const & pose, Size const jump_number) {
	std::list < Size > residue_list;
	Size downstream_res = pose.fold_tree().jump_edge(jump_number).stop();
	utility::vector1< Edge > edges = pose.fold_tree().get_outgoing_edges(downstream_res);

	// for jumps to singletons
	residue_list.push_back(downstream_res);

	for (Size i_edge = 1; i_edge <= edges.size(); ++i_edge) {
		if ( !edges[i_edge].is_polymer() ) continue;
		Size start = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].start() : edges[i_edge].stop();
		Size stop  = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].stop()  : edges[i_edge].start();
		for ( Size ires = start; ires <= stop; ++ires ) {
			residue_list.push_back(ires);
		}
	}
	residue_list.sort();
	residue_list.unique();
	return residue_list;
}
	
void
partial_align(
		core::pose::Pose & pose,
		core::pose::Pose const & ref_pose,
		id::AtomID_Map< id::AtomID > const & atom_map,
		std::list <Size> const & residue_list,
		bool iterate_convergence,
		utility::vector1<core::Real> distance_thresholds,
		core::Real min_coverage )
{
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT;
	numeric::xyzVector< core::Real > postT;

	// default
	if (distance_thresholds.size() == 0) {
		distance_thresholds.push_back(6);
		distance_thresholds.push_back(4);
		distance_thresholds.push_back(3);
		distance_thresholds.push_back(2);
		distance_thresholds.push_back(1.5);
		distance_thresholds.push_back(1);
	}

	get_superposition_transformation( pose, ref_pose, atom_map, R, preT, postT );
	apply_transformation( pose, residue_list, R, preT, postT );

	if (iterate_convergence) {
		core::id::AtomID_Map< core::id::AtomID > updated_atom_map(atom_map);
		core::Real coverage = 1.0;
		core::Size natoms_aln = atom_map_valid_size(pose, updated_atom_map);

		//std::cout << "coverage: " << coverage  << " " << natoms_aln << std::endl;

		for (int i=1; i<=(int)distance_thresholds.size() && coverage>=min_coverage; ++i) {
			core::Real my_d_sq = distance_thresholds[i]*distance_thresholds[i];
			bool converged = false;
			while (!converged) {
				core::id::AtomID_Map< core::id::AtomID > updated_atom_map_last_round = updated_atom_map;
				updated_atom_map = update_atom_map(pose, ref_pose, atom_map, my_d_sq);
				coverage = ((core::Real)(atom_map_valid_size(pose, updated_atom_map)))/natoms_aln;
				//std::cout << "coverage: " << coverage  << " " << natoms_aln << std::endl;
				if (updated_atom_map == updated_atom_map_last_round || coverage<min_coverage) {
					converged = true;
				} else {
					get_superposition_transformation( pose, ref_pose, updated_atom_map, R, preT, postT );
					apply_transformation( pose, residue_list, R, preT, postT );
				}
			}
		}
	}
}

core::Size atom_map_valid_size(
							   core::pose::Pose const & pose,
							   core::id::AtomID_Map< core::id::AtomID > const & atom_map
							   )
{
	core::Size n_valid = 0;
	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			++n_valid;
		}
	}
	return n_valid;
}

core::id::AtomID_Map< core::id::AtomID >
update_atom_map(
		core::pose::Pose & pose,
		core::pose::Pose const & ref_pose,
		id::AtomID_Map< id::AtomID > const & atom_map,
		core::Real distance_squared_threshold )
{
	core::id::AtomID_Map< core::id::AtomID > updated_atom_map;

	core::pose::initialize_atomid_map( updated_atom_map, pose, core::id::BOGUS_ATOM_ID );

	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			
			if (pose.xyz(id::AtomID( iatom,ires )).distance_squared( ref_pose.xyz(aid) ) < distance_squared_threshold)
				updated_atom_map[ id::AtomID( iatom,ires ) ] = aid;
		}
	}
	return updated_atom_map;
}

Size
natom_aligned(
		  core::pose::Pose & pose,
		  core::pose::Pose const & ref_pose,
		  id::AtomID_Map< id::AtomID > const & atom_map,
		  core::Real distance_squared_threshold
		  )
{
	Size n_align=0;
	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			
			if (pose.xyz(id::AtomID( iatom,ires )).distance_squared( ref_pose.xyz(aid) ) < distance_squared_threshold) {
				++n_align;
			}
		}
	}
	return n_align;
	
}
	
// atom_map: from mod_pose to ref_pose
void
get_superposition_transformation(
								 pose::Pose const & mod_pose,
								 pose::Pose const & ref_pose,
								 core::id::AtomID_Map< core::id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT )
{
	using namespace core;
	using namespace core::id;
	// count number of atoms for the array
	Size total_mapped_atoms(0);
	for ( Size ires=1; ires<= mod_pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			
			++total_mapped_atoms;
		}
	}
	
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	if (total_mapped_atoms <= 2) {
		R.xx() = R.yy() = R.zz() = 1;
		R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
		return;
	}
	
	ObjexxFCL::FArray2D< core::Real > final_coords( 3, total_mapped_atoms );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, total_mapped_atoms );
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	Size atomno(0);
	for ( Size ires=1; ires<= mod_pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			++atomno;
			
			numeric::xyzVector< core::Real > x_i = mod_pose.residue(ires).atom(iatom).xyz();
			preT += x_i;
			numeric::xyzVector< core::Real > y_i = ref_pose.xyz( aid );
			postT += y_i;
			
			for (int j=0; j<3; ++j) {
				init_coords(j+1,atomno) = x_i[j];
				final_coords(j+1,atomno) = y_i[j];
			}
		}
	}
	
	preT /= (float) total_mapped_atoms;
	postT /= (float) total_mapped_atoms;
	for (int i=1; i<=(int)total_mapped_atoms; ++i) {
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i) -= preT[j];
			final_coords(j+1,i) -= postT[j];
		}
	}
	
	// get optimal superposition
	// rotate >init< to >final<
	ObjexxFCL::FArray1D< numeric::Real > ww( total_mapped_atoms, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;
	
	numeric::model_quality::findUU( init_coords, final_coords, ww, total_mapped_atoms, uu, ctx );
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
}

void
apply_transformation(
	pose::Pose & mod_pose,
	std::list <Size> const & residue_list,
	numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
) {
	using namespace ObjexxFCL;
	// translate xx2 by COM and fill in the new ref_pose coordinates
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;
	
	for (std::list<Size>::const_iterator it = residue_list.begin();
		 it != residue_list.end();
		 ++it) {
		Size ires = *it;
		for ( Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
			ids.push_back(core::id::AtomID(iatom,ires));
			positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
		}
	}
	mod_pose.batch_set_xyz(ids,positions);
}

    
} // hybridize 
} // comparative_modeling 
} // protocols
