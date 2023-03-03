// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/ga_dock/util.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/util.hh>

#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
//#include <utility/graph/Graph.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

using namespace core;

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.util" );

core::Size
count_neighbors_on_coord( core::pose::Pose const &pose,
	core::Vector const &xyz1,
	std::string const atomname,
	core::Real const dcut )
{
	core::Size neighbor_counts( 0 );
	core::Real const dcut2( dcut*dcut );

	for ( core::Size resno = 1; resno <= pose.size(); ++resno ) {
		if ( pose.residue(resno).is_virtual_residue() ) continue;

		core::Vector xyz2;
		if ( atomname == "nbr" ) {
			xyz2 = pose.residue(resno).nbr_atom_xyz();
		} else {
			xyz2 = pose.residue(resno).has(atomname)? pose.residue(resno).xyz(atomname)
				: pose.residue(resno).nbr_atom_xyz();
		}
		core::Real d2 = xyz1.distance_squared(xyz2);
		if ( d2 < dcut2 ) neighbor_counts++;
	}

	return neighbor_counts;
}

utility::vector1< core::Size >
count_neighbors( core::pose::Pose const &pose,
	std::string const atomname,
	core::Real const dcut )
{
	utility::vector1< core::Size > neighbor_counts( pose.size(), 0 );

	core::scoring::TwelveANeighborGraph const & graph = pose.energies().twelveA_neighbor_graph();
	core::Real const dcut2( dcut*dcut );

	for ( core::Size resno = 1; resno <= pose.size(); ++resno ) {
		if ( pose.residue(resno).is_virtual_residue() ) continue;
		core::Vector xyz1;
		if ( atomname == "nbr" ) {
			xyz1 = pose.residue(resno).nbr_atom_xyz();
		} else {
			xyz1 = pose.residue(resno).has(atomname)? pose.residue(resno).xyz(atomname) : pose.residue(resno).nbr_atom_xyz();
		}

		for ( utility::graph::Graph::EdgeListConstIter
				ir  = graph.get_node( resno )->const_edge_list_begin(),
				ire = graph.get_node( resno )->const_edge_list_end();
				ir != ire; ++ir ) {
			core::Size const j( (*ir)->get_other_ind( resno ) );

			// supporting residue level only yet...
			core::Vector xyz2;
			if ( atomname == "nbr" ) {
				xyz2 = pose.residue(j).nbr_atom_xyz();
			} else {
				xyz2 = pose.residue(j).has(atomname)? pose.residue(j).xyz(atomname) : pose.residue(j).nbr_atom_xyz();
			}
			core::Real d2 = xyz1.distance_squared(xyz2);
			if ( d2 < dcut2 ) neighbor_counts[resno]++;
		}
	}

	return neighbor_counts;
}

utility::vector1< core::Size >
get_atomic_contacting_sidechains(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const &ligids,
	core::Real const atomic_distance_cut
) {
	utility::vector1< core::Size > contact_scs;
	core::Real const D2_COARSE( 225.0 );
	core::Real const D2_FINE( atomic_distance_cut*atomic_distance_cut );

	utility::vector1< core::Size > flexscs;

	//core::Vector const &ligcom = pose.residue(ligid_).nbr_atom_xyz();
	utility::vector1< core::Vector > ligCOMs;
	for ( auto ligid : ligids ) {
		ligCOMs.push_back ( pose.residue(ligid).nbr_atom_xyz() );
	}
	utility::vector1< bool > is_close( pose.size(), false );
	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_gly ||
				pose.residue(ires).aa() == core::chemical::aa_ala ||
				pose.residue(ires).aa() == core::chemical::aa_pro ) continue;

		core::Vector const & resnbr = pose.residue(ires).nbr_atom_xyz();

		core::Real d2min = D2_COARSE+1;
		for ( auto com : ligCOMs ) {
			d2min = std::min( d2min, com.distance_squared( resnbr ) );
		}
		if ( d2min < D2_COARSE ) is_close[ires] = true;
	}

	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		if ( !is_close[ires] ) continue;
		if ( std::find( ligids.begin(), ligids.end(), ires ) != ligids.end() ) continue;

		bool i_is_contacting = false;
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).nheavyatoms(); ++iatm ) {
			if ( pose.residue(ires).atom_is_backbone(iatm) ) continue;
			core::Vector const &xyz_i = pose.residue(ires).xyz(iatm);

			for ( auto ligid : ligids ) {
				for ( core::Size jatm = 1; jatm <= pose.residue(ligid).nheavyatoms(); ++jatm ) {
					core::Vector const &xyz_j = pose.residue(ligid).xyz(jatm);

					core::Real d2 = xyz_i.distance_squared( xyz_j );
					if ( d2 < D2_FINE ) {
						i_is_contacting = true;
						break;
					}
					if ( i_is_contacting ) break;
				}
				if ( i_is_contacting ) break;
			}
			if ( i_is_contacting ) break;
		} // iatm

		if ( i_is_contacting ) {
			contact_scs.push_back( ires );
		}
	} // ires

	return contact_scs;
}

void
constraint_relax( core::pose::Pose &pose,
	utility::vector1< core::Size > const &ligids,
	utility::vector1< core::Size > const &movable_scs
) {
	core::pose::Pose pose0(pose); //copy input pose for rmsd calculation
	core::scoring::ScoreFunctionOP scfxn = core::scoring::get_score_function();
	core::Real w_cart = (*scfxn)[core::scoring::cart_bonded];
	if ( w_cart<1.0e-5 ) {
		scfxn->set_weight(core::scoring::cart_bonded, 0.5);
		scfxn->set_weight(core::scoring::pro_close, 0.0);
	}
	core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
	mm->set_bb(false); mm->set_chi(false); mm->set_jump(false);
	for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
		mm->set_chi(movable_scs[i], true);
	}
	for ( core::Size i=1; i<=ligids.size(); ++i ) {
		mm->set_chi(ligids[i], true);
		mm->set_bb(ligids[i], true);
	}

	pose.remove_constraints();
	core::id::AtomID anchorid( 1, pose.fold_tree().root() );
	for ( core::Size ires = 1; ires <= movable_scs.size(); ++ires ) {
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).natoms(); ++iatm ) {
			core::id::AtomID atomid( iatm, ires );
			core::Vector const &xyz = pose.xyz( atomid );
			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
			pose.add_constraint( core::scoring::constraints::ConstraintCOP
				( core::scoring::constraints::ConstraintOP
				( new core::scoring::constraints::CoordinateConstraint( atomid, anchorid, xyz, fx ) )));
		}
	}

	for ( core::Size ires = 1; ires <= ligids.size(); ++ires ) {
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).natoms(); ++iatm ) {
			core::id::AtomID atomid( iatm, ires );
			core::Vector const &xyz = pose.xyz( atomid );
			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
			pose.add_constraint( core::scoring::constraints::ConstraintCOP
				( core::scoring::constraints::ConstraintOP
				( new core::scoring::constraints::CoordinateConstraint( atomid, anchorid, xyz, fx ) )));
		}
	}

	core::optimization::CartesianMinimizer minimizer;
	core::optimization::MinimizerOptions options("lbfgs_armijo", 0.0001, true, false);
	options.max_iter(50);
	minimizer.run(pose, *mm, *scfxn, options);

	core::Real rms_scs = core::scoring::all_atom_rmsd_nosuper( pose0, pose,
		movable_scs, movable_scs );
	core::Real rms_lig = core::scoring::all_atom_rmsd_nosuper( pose0, pose,
		ligids, ligids );
	core::pose::setPoseExtraScore( pose, "cst_rms_scs", rms_scs );
	core::pose::setPoseExtraScore( pose, "cst_rms_lig", rms_lig );
	if ( TR.Debug.visible() ) {
		TR << "After minimization, sidechain rms: " << rms_scs << ", " << "ligand rms: " << rms_lig << std::endl;
	}

}

void
make_ligand_only_pose(
	core::pose::PoseOP pose_new,
	core::pose::PoseCOP pose,
	utility::vector1< core::Size > const& lig_resnos
) {

	TR << "pose.fold_tree() " << pose->fold_tree() << std::endl;
	core::kinematics::FoldTree f(lig_resnos.size());
	core::pose::create_subpose(*pose, lig_resnos, f, *pose_new);
	//core::pose::pdbslice(*pose_new, *pose, lig_resnos);

	core::Real dH, TdS, dG;
	std::string ligandname;
	core::Real rms, complexscore, ligscore, recscore, ranking_prerelax;
	core::pose::getPoseExtraScore( *pose, "dH", dH);
	core::pose::getPoseExtraScore( *pose, "-TdS", TdS);
	core::pose::getPoseExtraScore( *pose, "dG", dG);
	core::pose::getPoseExtraScore( *pose, "lig_rms", rms);
	core::pose::getPoseExtraScore( *pose, "ligscore", ligscore);
	core::pose::getPoseExtraScore( *pose, "recscore", recscore);
	core::pose::getPoseExtraScore( *pose, "complexscore", complexscore );
	core::pose::getPoseExtraScore( *pose, "ranking_prerelax", ranking_prerelax);
	core::pose::getPoseExtraScore( *pose, "ligandname", ligandname);

	core::pose::setPoseExtraScore( *pose_new, "-TdS", TdS );
	core::pose::setPoseExtraScore( *pose_new, "dG", dG );
	core::pose::setPoseExtraScore( *pose_new, "dH", dH );
	core::pose::setPoseExtraScore( *pose_new, "lig_rms", rms);
	core::pose::setPoseExtraScore( *pose_new, "ligscore", ligscore );
	core::pose::setPoseExtraScore( *pose_new, "recscore", recscore );
	core::pose::setPoseExtraScore( *pose_new, "complexscore", complexscore );
	core::pose::setPoseExtraScore( *pose_new, "ranking_prerelax", ranking_prerelax );
	core::pose::setPoseExtraScore( *pose_new, "ligandname", ligandname);

}

} // ga_dock
} // ligand_docking
} // protocols


