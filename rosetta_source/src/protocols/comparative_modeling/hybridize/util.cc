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
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>
#include <core/types.hh>

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
	for (int i=1; i<=templates.size(); ++i) {
		utility::vector1< bool > passed_gapcheck(nres_tgt,false);
		for (int j=1; j<templates[i]->total_residue(); ++j ) {
			bool includeme=true;
			for (int k=1; k<=GAPBUFFER && includeme; ++k) {
				if ( j-k < 1 || j+k > templates[i]->total_residue() ) includeme=false;
				else if (templates[i]->pdb_info()->number(j+k) - templates[i]->pdb_info()->number(j) != k ) includeme=false;
				else if (templates[i]->pdb_info()->number(j-k) - templates[i]->pdb_info()->number(j) != -k ) includeme=false;
			}
			passed_gapcheck[j] = includeme;
		}

		for (core::Size j=1; j<templates[i]->total_residue(); ++j ) {
			if (!passed_gapcheck[j]) continue;
			for (core::Size k=j+1; k<templates[i]->total_residue(); ++k ) {
				if (!passed_gapcheck[k]) continue;
				if (templates[i]->pdb_info()->number(k) - templates[i]->pdb_info()->number(j) < MINSEQSEP) continue;
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
    
} // hybridize 
} // comparative_modeling 
} // protocols
