// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Align a random jump to template
/// @details
/// @author Yifan Song

#include <protocols/hybridization/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>

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
#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMover.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/USOGFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/PDBInfo.hh>

// dynamic fragpick
#include <protocols/moves/DsspMover.hh>
#include <core/fragment/picking_old/vall/util.hh>


// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// evaluation
#include <core/scoring/rms_util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <protocols/comparative_modeling/coord_util.hh>

#include <list>
#include <numeric/xyzVector.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR( "protocols.hybridization.util" );

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

using namespace core;
using namespace kinematics;
using namespace core::scoring::constraints;

core::Size
get_num_residues_nonvirt( core::pose::Pose const & pose ) {
	core::Size nres_tgt = pose.size();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}
	while ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;
	return nres_tgt;
}

core::Size
get_num_residues_prot( core::pose::Pose const & pose ) {
	core::Size nres_tgt = pose.size();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}
	while ( !pose.residue(nres_tgt).is_protein() ) nres_tgt--;
	return nres_tgt;
}

void setup_centroid_constraints(
	core::pose::Pose &pose,
	utility::vector1< core::pose::PoseCOP > const & templates,
	utility::vector1< core::Real > const & template_weights,
	std::string const & cen_cst_file,
	std::set< core::Size > const & ignore_res_for_AUTO) {
	if ( cen_cst_file == "AUTO" ) {
		// automatic constraints
		generate_centroid_constraints( pose, templates, template_weights, ignore_res_for_AUTO );
	} else if ( !cen_cst_file.empty() && cen_cst_file != "NONE" ) {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( cen_cst_file, ConstraintSetOP( new ConstraintSet ), pose );
		pose.constraint_set( constraint_set );  //reset constraints
	}
}


void setup_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1< core::pose::PoseCOP > const & templates,
	utility::vector1< core::Real > const & template_weights,
	std::string const & cen_cst_file,
	std::string const & fa_cst_file  ) {

	// use fa if specified, otherwise centroid
	if ( fa_cst_file == "AUTO" ) {
		// automatic fa constraints
		generate_fullatom_constraints( pose, templates, template_weights );
	} else if ( fa_cst_file == "SELF" ) {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover add_constraints;
		add_constraints.apply(pose);
	} else if ( !fa_cst_file.empty() && fa_cst_file != "NONE" ) {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( fa_cst_file, ConstraintSetOP( new ConstraintSet ), pose );
		pose.constraint_set( constraint_set );  //reset constraints
	} else if ( cen_cst_file == "AUTO" ) {
		// automatic constraints
		generate_centroid_constraints( pose, templates, template_weights );
	} else if ( !cen_cst_file.empty() && cen_cst_file != "NONE" ) {
		ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( cen_cst_file, ConstraintSetOP( new ConstraintSet ), pose );
		pose.constraint_set( constraint_set );  //reset constraints
	}
}

void setup_constraints(
	core::pose::Pose &pose,
	std::string & cst_in) {
	std::istringstream cst_ss(cst_in);
	ConstraintSetOP constraint_set = ConstraintIO::get_instance()->read_constraints_new( cst_ss, ConstraintSetOP( new ConstraintSet ), pose );
	pose.constraint_set( constraint_set );  //reset constraints
}


void generate_centroid_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > /*template_weights*/,
	std::set< core::Size > ignore_res )
{
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Size MINSEQSEP = 8;
	core::Real MAXDIST = 12.0;
	core::Size GAPBUFFER = 3;
	core::Real COORDDEV = 1.0;

	pose.remove_constraints();

	// number of residues
	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	utility::vector1< utility::vector1< core::Real > > tgt_dists(nres_tgt);
	utility::vector1< utility::vector1< core::Real > > tgt_weights(nres_tgt);
	for ( int i=1; i<=(int)templates.size(); ++i ) {
		utility::vector1< bool > passed_gapcheck(nres_tgt,false);
		for ( int j=1; j<(int)templates[i]->size(); ++j ) {
			bool includeme=true;
			for ( int k=1; k<=(int)GAPBUFFER && includeme; ++k ) {
				if ( j-k < 1 || j+k > (int)templates[i]->size() ) includeme=false;
				else if ( templates[i]->pdb_info()->number(j+k) - templates[i]->pdb_info()->number(j) != k ) includeme=false;
				else if ( templates[i]->pdb_info()->number(j-k) - templates[i]->pdb_info()->number(j) != -k ) includeme=false;
			}
			passed_gapcheck[j] = includeme;
		}

		for ( core::Size j=1; j<templates[i]->size(); ++j ) {
			if ( !templates[i]->residue_type(j).is_protein() ) continue;
			if ( !passed_gapcheck[j] ) continue;
			for ( core::Size k=j+1; k<templates[i]->size(); ++k ) {
				if ( !templates[i]->residue_type(k).is_protein() ) continue;
				if ( !passed_gapcheck[k] ) continue;
				if ( templates[i]->pdb_info()->number(k) - templates[i]->pdb_info()->number(j) < (int)MINSEQSEP ) continue;
				if ( ignore_res.count(templates[i]->pdb_info()->number(j)) ||
						ignore_res.count(templates[i]->pdb_info()->number(k)) ) continue;
				core::Real dist = templates[i]->residue(j).xyz(2).distance( templates[i]->residue(k).xyz(2) );
				if ( dist <= MAXDIST ) {
					core::Size resid_j = templates[i]->pdb_info()->number(j);
					core::Size resid_k = templates[i]->pdb_info()->number(k);

					if ( symm_info && !symm_info->bb_is_independent( resid_j ) ) {
						resid_j = symm_info->bb_follows( resid_j );
					}
					if ( symm_info && !symm_info->bb_is_independent( resid_k ) ) {
						resid_k = symm_info->bb_follows( resid_k );
					}

					using namespace core::scoring::func;
					FuncOP fx( new ScalarWeightedFunc( 1.0, FuncOP( new USOGFunc( dist, COORDDEV ) ) ) );
					pose.add_constraint(
						scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint( core::id::AtomID(2,resid_j), core::id::AtomID(2,resid_k), fx ) ) )
					);
				}
			}
		}
	}
}


void setup_user_coordinate_constraints(
	core::pose::Pose &pose,
	utility::vector1<Size> reses )
{
	core::Real COORDDEV = 0.39894;

	for ( core::Size i=1; i<=reses.size(); ++i ) {
		core::scoring::func::FuncOP fx( new core::scoring::func::USOGFunc( 0, COORDDEV ) );
		pose.add_constraint(
			scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint(
			core::id::AtomID(2,reses[i]),
			core::id::AtomID(2,pose.size()),
			pose.residue(reses[i]).atom(2).xyz(),
			fx
			) ) ) );
	}
}


void generate_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights ) {
	//fpd ... for now just use centroid variant
	generate_centroid_constraints( pose, templates, template_weights);
}

void add_strand_pairs_cst(core::pose::Pose & pose, utility::vector1< std::pair< core::Size, core::Size > > const strand_pairs) {
	core::Real MAXDIST = 12.0;
	core::Real COORDDEV = 1.0;
	for ( core::Size i=1; i<=strand_pairs.size(); ++i ) {
		std::pair< core::Size, core::Size > strand_pair = strand_pairs[i];
		core::Real dist = pose.residue(strand_pair.first).xyz(2).distance( pose.residue(strand_pair.second).xyz(2) );
		if ( dist <= MAXDIST ) {
			using namespace core::scoring::func;
			FuncOP fx( new ScalarWeightedFunc( 4.0, FuncOP( new USOGFunc( dist, COORDDEV ) ) ) ); // try to lock it down with a high weight
			pose.add_constraint(
				scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint( core::id::AtomID(2,strand_pair.first), core::id::AtomID(2,strand_pair.second), fx ) ) )
			);
		}
	}
}

void add_non_protein_cst(core::pose::Pose & pose, core::pose::Pose & tmpl, core::Real const self_cst_weight, core::Real const het_prot_cst_weight) {
	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Real MAXDIST = 8.0;
	core::Real COORDDEV = 3.0;

	//fpd -- generate constraints between the residues in tmpl with the hetatms in pose
	if ( het_prot_cst_weight > 1e-7 ) {
		// constraints protein<->substrate
		for ( Size ires=1; ires<=tmpl.size(); ++ires ) {
			if ( !tmpl.residue(ires).is_protein() ) continue;

			core::Size iatom;
			if ( tmpl.residue(ires).aa() == core::chemical::aa_gly ) {
				iatom = tmpl.residue_type(ires).atom_index(" CA ");
			} else {
				iatom = tmpl.residue_type(ires).atom_index(" CB ");
			}

			core::Size ires_tgt = (tmpl.pdb_info()) ? tmpl.pdb_info()->number(ires) : ires;
			if ( symm_info && !symm_info->bb_is_independent(ires_tgt) ) {
				ires_tgt = symm_info->bb_follows( ires_tgt );
			}

			for ( Size jres=1; jres<=tmpl.size(); ++jres ) {
				if ( tmpl.residue(jres).is_protein() || tmpl.residue(jres).aa() == core::chemical::aa_vrt ) continue;

				core::Size jres_tgt = (tmpl.pdb_info()) ? tmpl.pdb_info()->number(ires) : ires;
				if ( symm_info && !symm_info->bb_is_independent(jres_tgt) ) {
					jres_tgt = symm_info->bb_follows( jres_tgt );
				}

				for ( Size jatom=1; jatom<=tmpl.residue(jres).nheavyatoms(); ++jatom ) {
					core::Real dist = tmpl.residue(ires).xyz(iatom).distance( tmpl.residue(jres).xyz(jatom) );

					if ( dist <= MAXDIST ) {
						using namespace core::scoring::func;
						FuncOP fx( new ScalarWeightedFunc( het_prot_cst_weight, FuncOP( new USOGFunc( dist, COORDDEV ) ) ) );
						pose.add_constraint(
							scoring::constraints::ConstraintCOP(
							new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires_tgt), core::id::AtomID(jatom,jres_tgt), fx ) )
						);
					}
				}
			}
		}
	}

	MAXDIST = 8.0;
	COORDDEV = 1.0;

	if ( self_cst_weight > 1e-7 ) {
		// constraints within substrate
		for ( Size ires=1; ires<=tmpl.size(); ++ires ) {
			if ( tmpl.residue(ires).is_protein() || tmpl.residue(ires).aa() == core::chemical::aa_vrt ) continue;

			core::Size ires_tgt = (tmpl.pdb_info()) ? tmpl.pdb_info()->number(ires) : ires;
			if ( symm_info && !symm_info->bb_is_independent(ires_tgt) ) {
				ires_tgt = symm_info->bb_follows( ires_tgt );
			}

			for ( Size iatom=1; iatom<=tmpl.residue(ires).nheavyatoms(); ++iatom ) {
				for ( Size jres=1; jres<=tmpl.size(); ++jres ) {
					if ( tmpl.residue(jres).is_protein() || tmpl.residue(jres).aa() == core::chemical::aa_vrt ) continue;

					core::Size jres_tgt = (tmpl.pdb_info()) ? tmpl.pdb_info()->number(jres) : jres;
					if ( symm_info && !symm_info->bb_is_independent(jres_tgt) ) {
						jres_tgt = symm_info->bb_follows( jres_tgt );
					}

					for ( Size jatom=1; jatom<=tmpl.residue(jres).nheavyatoms(); ++jatom ) {
						if ( ires == jres && iatom <= jatom ) continue;

						core::Real dist1 = tmpl.residue(ires).xyz(iatom).distance( tmpl.residue(jres).xyz(jatom) );

						if ( dist1 <= MAXDIST ) {
							using namespace core::scoring::func;
							FuncOP fx( new ScalarWeightedFunc( self_cst_weight, FuncOP( new USOGFunc( dist1, COORDDEV ) ) ) );
							pose.add_constraint(
								scoring::constraints::ConstraintCOP(
								new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires_tgt), core::id::AtomID(jatom,jres_tgt), fx ) )
							);
						}
					}
				}
			}
		}
	}
}

bool discontinued_upper(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);

	if ( seqpos == 1 ) return true;
	if ( !pose.residue_type(seqpos).is_polymer() ) return true;
	if ( !pose.residue_type(seqpos-1).is_polymer() ) return true;
	if ( pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos-1).is_protein() ) {
		if ( pose.residue(seqpos).xyz("N").distance(pose.residue(seqpos-1).xyz("C")) > N_C_cutoff ) {
			return true;
		}
	}
	return false;
}

bool discontinued_lower(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);

	if ( seqpos == pose.size() ) return true;
	if ( !pose.residue_type(seqpos).is_polymer() ) return true;
	if ( !pose.residue_type(seqpos+1).is_polymer() ) return true;
	if ( pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos+1).is_protein() ) {
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

	TR.Debug << "Downstream residue " << downstream_res << std::endl;

	for ( Size i_edge = 1; i_edge <= edges.size(); ++i_edge ) {
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

// iterative partial alignment
void
partial_align(
	core::pose::Pose & pose,
	core::pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map,
	bool iterate_convergence,
	utility::vector1<core::Real> distance_thresholds,
	core::Real min_coverage )
{
	std::list <core::Size> residue_list;
	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		residue_list.push_back(ires);
	}

	partial_align( pose, ref_pose, atom_map, residue_list, iterate_convergence, distance_thresholds, min_coverage );
}

// iterative partial alignment over a subset of residues
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
	if ( distance_thresholds.size() == 0 ) {
		distance_thresholds.push_back(6);
		distance_thresholds.push_back(4);
		distance_thresholds.push_back(3);
		distance_thresholds.push_back(2);
		distance_thresholds.push_back(1.5);
		distance_thresholds.push_back(1);
	}

	get_superposition_transformation( pose, ref_pose, atom_map, R, preT, postT );
	apply_transformation( pose, residue_list, R, preT, postT );

	if ( iterate_convergence ) {
		core::id::AtomID_Map< core::id::AtomID > updated_atom_map(atom_map);
		core::Real coverage = 1.0;
		core::Size natoms_aln = atom_map_valid_size(pose, updated_atom_map);

		for ( int i=1; i<=(int)distance_thresholds.size() && coverage>=min_coverage; ++i ) {
			core::Real my_d_sq = distance_thresholds[i]*distance_thresholds[i];
			bool converged = false;
			while ( !converged ) {
				core::id::AtomID_Map< core::id::AtomID > updated_atom_map_last_round = updated_atom_map;
				updated_atom_map = update_atom_map(pose, ref_pose, atom_map, my_d_sq);
				coverage = ((core::Real)(atom_map_valid_size(pose, updated_atom_map)))/natoms_aln;
				//std::cout << "coverage: " << coverage  << " " << natoms_aln << std::endl;
				if ( updated_atom_map == updated_atom_map_last_round || coverage<min_coverage ) {
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
	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;
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

	core::pose::initialize_atomid_map( updated_atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() );

	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;

			if ( pose.xyz(id::AtomID( iatom,ires )).distance_squared( ref_pose.xyz(aid) ) < distance_squared_threshold ) {
				updated_atom_map[ id::AtomID( iatom,ires ) ] = aid;
			}
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
	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= pose.residue(ires).natoms(); ++iatom ) {
			core::id::AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;

			if ( pose.xyz(id::AtomID( iatom,ires )).distance_squared( ref_pose.xyz(aid) ) < distance_squared_threshold ) {
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
	for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;

			++total_mapped_atoms;
		}
	}

	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	if ( total_mapped_atoms <= 2 ) {
		R.xx() = R.yy() = R.zz() = 1;
		R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
		return;
	}

	ObjexxFCL::FArray2D< core::Real > final_coords( 3, total_mapped_atoms );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, total_mapped_atoms );
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	Size atomno(0);
	for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if ( !aid.valid() ) continue;
			++atomno;

			numeric::xyzVector< core::Real > x_i = mod_pose.residue(ires).atom(iatom).xyz();
			preT += x_i;
			numeric::xyzVector< core::Real > y_i = ref_pose.xyz( aid );
			postT += y_i;

			for ( int j=0; j<3; ++j ) {
				init_coords(j+1,atomno) = x_i[j];
				final_coords(j+1,atomno) = y_i[j];
			}
		}
	}

	preT /= (float) total_mapped_atoms;
	postT /= (float) total_mapped_atoms;
	for ( int i=1; i<=(int)total_mapped_atoms; ++i ) {
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

	for ( core::Size ires : residue_list ) {
		for ( Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
			ids.push_back(core::id::AtomID(iatom,ires));
			positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
		}
	}
	mod_pose.batch_set_xyz(ids,positions);
}

core::fragment::FragSetOP
create_fragment_set( core::pose::Pose const & pose, core::Size len, core::Size nfrag ) {
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( len ) );
	core::scoring::dssp::Dssp dssp( pose );

	// number of residues
	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	// sequence
	std::string tgt_seq = pose.sequence();
	std::string tgt_ss = dssp.get_dssp_secstruct();

	// pick from vall based on template SS + target sequence
	for ( core::Size j=1; j<=nres_tgt-len+1; ++j ) {
		bool protein_only = true;
		for ( core::Size k=j; k<j+len; ++k ) {   // it's alright if the last residue is a cutpoint
			protein_only &= pose.residue(k).is_protein(  );
		}

		bool crosses_cut = false;
		for ( core::Size k=j; k<j+len-1; ++k ) {   // it's alright if the last residue is a cutpoint
			crosses_cut |= pose.fold_tree().is_cutpoint( k );
		}

		if ( !crosses_cut && protein_only ) {
			core::fragment::FrameOP frame( new core::fragment::Frame( j, len ) );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(
				tgt_ss.substr( j-1, len ), tgt_seq.substr( j-1, len ), nfrag, true, core::fragment::IndependentBBTorsionSRFD()
				)
			);
			fragset->add( frame );
		}
	}
	return fragset;
}

core::fragment::FragSetOP
create_fragment_set_no_ssbias( core::pose::Pose const & pose, core::Size len, core::Size nfrag, char force_ss ) {
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( len ) );

	// number of residues
	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	// sequence
	std::string tgt_seq = pose.sequence();
	std::string tgt_ss (len, force_ss);

	// pick from vall based on template SS + target sequence
	for ( core::Size j=1; j<=nres_tgt-len+1; ++j ) {
		bool protein_only = true;
		for ( core::Size k=j; k<j+len; ++k ) {   // it's alright if the last residue is a cutpoint
			protein_only &= pose.residue(k).is_protein(  );
		}

		bool crosses_cut = false;
		for ( core::Size k=j; k<j+len-1; ++k ) {   // it's alright if the last residue is a cutpoint
			crosses_cut |= pose.fold_tree().is_cutpoint( k );
		}

		if ( !crosses_cut && protein_only ) {
			core::fragment::FrameOP frame( new core::fragment::Frame( j, len ) );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(
				tgt_ss, tgt_seq.substr( j-1, len ), nfrag, true, core::fragment::IndependentBBTorsionSRFD()
				)
			);
			fragset->add( frame );
		}
	}
	return fragset;
}


core::fragment::FragSetOP
create_fragment_set_no_ssbias( core::pose::Pose const & pose, std::set<core::Size> user_pos, core::Size len, core::Size nfrag, char force_ss ) {
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( len ) );

	// number of residues
	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	// sequence
	std::string tgt_seq = pose.sequence();
	std::string tgt_ss (len, force_ss);

	bool anyLegalpos = false;

	// pick from vall based on template SS + target sequence
	for ( core::Size j=1; j<=nres_tgt-len+1; ++j ) {
		bool legalpos = true;
		for ( core::Size k=j; k<j+len-1; ++k ) {
			legalpos &= !pose.fold_tree().is_cutpoint( k );
		}

		for ( core::Size k=j; k<j+len; ++k ) {
			legalpos &= pose.residue(k).is_protein(  );
		}

		// allow 1 position feathering
		for ( core::Size k=j+1; k<j+len-1; ++k ) {
			legalpos &= (user_pos.find( k ) != user_pos.end());
		}

		if ( legalpos ) {
			anyLegalpos = true;
			core::fragment::FrameOP frame( new core::fragment::Frame( j, len ) );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(
				tgt_ss, tgt_seq.substr( j-1, len ), nfrag, true, core::fragment::IndependentBBTorsionSRFD()
				)
			);
			fragset->add( frame );
		}
	}

	if ( !anyLegalpos ) {
		//utility_exit_with_message( "No legal positions allowed in create_fragment_set_no_ssbias!" );
		TR << "WARNING! No legal positions allowed in create_fragment_set_no_ssbias!" << std::endl;
	}

	return fragset;
}

core::fragment::FragSetOP
create_fragment_set_no_ssbias( std::string tgt_seq, core::Size len, core::Size nfrag, char force_ss ) {
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( len ));

	// number of residues
	core::Size nres_tgt = tgt_seq.length();

	// sequence
	std::string tgt_ss (len, force_ss);

	// pick from vall based on template SS + target sequence
	for ( core::Size j=1; j<=nres_tgt-len+1; ++j ) {
		bool crosses_cut = false;
		for ( core::Size k=j; k<j+len-1; ++k ) {   // it's alright if the last residue is a cutpoint
			crosses_cut |= (tgt_seq[k]=='X' || tgt_seq[k]=='Z');
		}

		if ( !crosses_cut ) {
			core::fragment::FrameOP frame( new core::fragment::Frame( j, len ) );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(
				tgt_ss, tgt_seq.substr( j-1, len ), nfrag, true, core::fragment::IndependentBBTorsionSRFD()
				)
			);
			fragset->add( frame );
		}
	}
	return fragset;
}

protocols::loops::Loops renumber_with_pdb_info(
	protocols::loops::Loops & template_chunk,
	core::pose::PoseCOP template_pose
) {
	protocols::loops::Loops renumbered_template_chunks(template_chunk);
	for ( core::Size ichunk = 1; ichunk<=template_chunk.num_loop(); ++ichunk ) {
		Size seqpos_start_templ = template_chunk[ichunk].start();
		Size seqpos_start_target = template_pose->pdb_info()->number(seqpos_start_templ);
		renumbered_template_chunks[ichunk].set_start( seqpos_start_target );

		Size seqpos_stop_templ = template_chunk[ichunk].stop();
		Size seqpos_stop_target = template_pose->pdb_info()->number(seqpos_stop_templ);
		renumbered_template_chunks[ichunk].set_stop( seqpos_stop_target );
	}
	return renumbered_template_chunks;
}


core::Real get_gdtmm( core::pose::Pose & native, core::pose::Pose & pose, core::sequence::SequenceAlignmentOP & aln ) {
	runtime_assert( native.size() > 0 && pose.size() > 0 );
	if ( !aln ) {
		core::sequence::SequenceOP model_seq( new core::sequence::Sequence( pose.sequence(),  "model",  1 ) );
		core::sequence::SequenceOP native_seq( new core::sequence::Sequence( native.sequence(), "native", 1 ) );
		aln = core::sequence::SequenceAlignmentOP( new core::sequence::SequenceAlignment );
		*aln = align_naive(model_seq,native_seq);
	}

	int n_atoms;
	ObjexxFCL::FArray2D< core::Real > p1a, p2a;
	protocols::comparative_modeling::gather_coords( pose, native, *aln, n_atoms, p1a, p2a );

	core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	return core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
}


} // hybridize
//} // comparative_modeling
} // protocols
