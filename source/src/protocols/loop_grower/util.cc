// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief utility functions for the LoopGrowing protocol
/// @details
/// @author Brandon Frenz

#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

// Symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <protocols/rigid/RB_geometry.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace loop_grower {

static THREAD_LOCAL basic::Tracer TR( "protocols.loop_grower.util" );

using namespace protocols;
using namespace core;

void
transform_res_to_subunit( core::pose::Pose &pose, core::conformation::ResidueOP xformed, core::Size symmcopy) {
				core::conformation::symmetry::SymmetricConformation & SymmConf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
				core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
        core::Size nres_asu = symm_info->num_independent_residues();

        for (Size i=1; i<=xformed->atoms().size(); ++i) {
                numeric::xyzVector< core::Real > X = xformed->xyz(i);
                numeric::xyzVector< core::Real > sX = SymmConf.apply_transformation( X, 1, symmcopy*nres_asu);
                xformed->set_xyz(i, sX);
        }
}

void
transform_to_closest_symmunit(core::pose::Pose & cen_pose, core::pose::Pose & fa_pose, Size lower_pose){
        core::conformation::symmetry::SymmetricConformation & SymmConf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( cen_pose.conformation()) );
				core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
        Size nsubunits = symm_info->subunits();
				//add terminus checks!!!!
				if( lower_pose == symm_info->num_independent_residues() ) return;
				if( lower_pose == 1 ) return;
				Size rescount = lower_pose;
				while (!cen_pose.fold_tree().is_cutpoint(rescount)) rescount++;
				core::conformation::ResidueOP xformed = cen_pose.conformation().residue_op(rescount+1)->clone();
				std::string mobileatom = "C";
				std::string staticatom = "N";
				Real shortest_dist = 9999;
				Size bestsubunit = 1;
				for(Size i=1; i<=nsubunits; i++){
					transform_res_to_subunit(cen_pose, xformed, i);
					Real distance = (cen_pose.residue(rescount).atom(staticatom).xyz()-xformed->atom(mobileatom).xyz()).length();
					if(distance < shortest_dist){
						shortest_dist = distance;
						bestsubunit = i;
					}
				}
				//replace with best symmetric unit
				bool ischainbreak = false;
				for(Size i=rescount+1; i<=cen_pose.total_residue(); i++){
						if(cen_pose.fold_tree().is_cutpoint(i)){
								ischainbreak = true;
						}
						xformed = cen_pose.conformation().residue_op(i)->clone();
						transform_res_to_subunit(cen_pose, xformed, bestsubunit);
						cen_pose.replace_residue(i,*xformed,false);
						xformed = fa_pose.conformation().residue_op(i)->clone();
						transform_res_to_subunit(fa_pose, xformed, bestsubunit);
						fa_pose.replace_residue(i,*xformed,false);
						if( ischainbreak ) break;
				}
}


}
}
