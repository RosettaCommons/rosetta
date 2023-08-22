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
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>

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
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/conversions.hh>
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
	utility::vector1< core::Size > const &movable_scs,
	core::Real maxiter
) {
	core::Size jumpid = get_ligand_jumpid( pose, ligids );
	core::pose::Pose pose0(pose); //copy input pose for rmsd calculation
	core::scoring::ScoreFunctionOP scfxn = core::scoring::get_score_function();
	core::Real w_cst = (*scfxn)[core::scoring::coordinate_constraint];
	if ( TR.Debug.visible() ) TR.Debug << "weight coordinate_constraint: " << w_cst << std::endl;
	core::Real w_cart = (*scfxn)[core::scoring::cart_bonded];
	if ( w_cart<1.0e-5 ) {
		scfxn->set_weight(core::scoring::cart_bonded, 0.5);
		scfxn->set_weight(core::scoring::pro_close, 0.0);
	}
	core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
	mm->set_bb(false); mm->set_chi(false); mm->set_jump(false);
	mm->set_jump(jumpid, true);
	for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
		mm->set_chi(movable_scs[i], true);
	}
	for ( core::Size i=1; i<=ligids.size(); ++i ) {
		mm->set_chi(ligids[i], true);
		mm->set_bb(ligids[i], true);
	}

	pose.remove_constraints();
	core::id::AtomID anchorid( 1, pose.fold_tree().root() );
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	for ( core::Size ires : movable_scs ) {
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).natoms(); ++iatm ) {
			core::id::AtomID atomid( iatm, ires );
			core::Vector const &xyz = pose.xyz( atomid );
			pose.add_constraint( core::scoring::constraints::ConstraintCOP
				( core::scoring::constraints::ConstraintOP
				( new core::scoring::constraints::CoordinateConstraint( atomid, anchorid, xyz, fx ) )));
		}
	}

	for ( core::Size ires : ligids ) {
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).natoms(); ++iatm ) {
			core::id::AtomID atomid( iatm, ires );
			core::Vector const &xyz = pose.xyz( atomid );
			pose.add_constraint( core::scoring::constraints::ConstraintCOP
				( core::scoring::constraints::ConstraintOP
				( new core::scoring::constraints::CoordinateConstraint( atomid, anchorid, xyz, fx ) )));
		}
	}

	core::optimization::CartesianMinimizer minimizer;
	core::optimization::MinimizerOptions options("lbfgs_armijo", 0.0001, true, false);
	options.max_iter(maxiter);
	minimizer.run(pose, *mm, *scfxn, options);
	pose.remove_constraints();

	core::Real rms_scs = core::scoring::all_atom_rmsd_nosuper( pose0, pose,
		movable_scs, movable_scs );
	core::Real rms_lig = core::scoring::all_atom_rmsd_nosuper( pose0, pose,
		ligids, ligids );
	core::pose::setPoseExtraScore( pose, "cst_rms_scs", rms_scs );
	core::pose::setPoseExtraScore( pose, "cst_rms_lig", rms_lig );
	if ( TR.Debug.visible() ) {
		TR.Debug << "After minimization, sidechain rms: " << rms_scs << ", " << "ligand rms: " << rms_lig << std::endl;
	}
}

void
make_ligand_only_pose(
	core::pose::PoseOP pose_new,
	core::pose::PoseCOP pose,
	utility::vector1< core::Size > const& lig_resnos
) {

	core::kinematics::FoldTree f(lig_resnos.size());
	core::pose::create_subpose(*pose, lig_resnos, f, *pose_new);
	//core::pose::pdbslice(*pose_new, *pose, lig_resnos);

	core::Real dH, TdS, dG, buns, n_hbonds_total, n_hbonds_max1;
	std::string ligandname;
	core::Real rms, complexscore, ligscore, recscore, ranking_prerelax;
	bool retval;
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

	retval = core::pose::getPoseExtraScore( *pose, "buns", buns);
	if ( retval ) core::pose::setPoseExtraScore( *pose_new, "buns", buns );

	retval = core::pose::getPoseExtraScore( *pose, "n_hbonds_total", n_hbonds_total);
	if ( retval ) core::pose::setPoseExtraScore( *pose_new, "n_hbonds_total", n_hbonds_total );

	retval = core::pose::getPoseExtraScore( *pose, "n_hbonds_max1", n_hbonds_max1);
	if ( retval ) core::pose::setPoseExtraScore( *pose_new, "n_hbonds_max1", n_hbonds_max1 );

}

// perturb ligand rigid body
void
perturb_ligand_rb(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const& ligids,
	core::Real trans_step,
	core::Real rot_step
) {
	core::Vector Raxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	core::Real angle = rot_step * numeric::NumericTraits<core::Real>::deg2rad() * numeric::random::rg().gaussian();
	numeric::xyzMatrix< core::Real > R = numeric::rotation_matrix( Raxis, angle );

	core::Vector Tcom( 0.0,0.0,0.0 );
	core::Size natm = 0;
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			Tcom += pose.xyz(core::id::AtomID(iatm_src,ligid));
		}
		natm += pose.residue_type(ligid).natoms();
	}
	Tcom /= natm;

	core::Real len = trans_step * numeric::random::rg().uniform();
	core::Vector T( len*numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	for ( auto ligid : ligids ) {
		for ( core::Size iatm_src=1; iatm_src<=pose.residue_type(ligid).natoms(); ++iatm_src ) {
			pose.set_xyz(
				core::id::AtomID(iatm_src,ligid),
				R*( pose.xyz(core::id::AtomID(iatm_src,ligid)) - Tcom) + T + Tcom
			);
		}
	}

}

// perturb ligand internal torsions
void
perturb_ligand_torsions(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const& ligids,
	utility::vector1< core::Size > const& freeze_chi,
	core::Real chi_step
) {
	if ( ligids.size() >1 && freeze_chi.size() >0 ) {
		utility_exit_with_message("perturb_ligand_torsions: freeze_chi not implemented for multiple ligands");
	}
	// internal torsions
	for ( auto ligid : ligids ) {
		core::Size nchi = pose.residue_type(ligid).nchi();
		for ( core::Size ichi_src=1; ichi_src<=nchi; ++ichi_src ) {
			if ( freeze_chi.size() > 0 && freeze_chi[ichi_src] ) continue;
			core::chemical::AtomIndices const & chiatoms = pose.residue_type(ligid).chi_atoms( ichi_src );
			core::chemical::BondName bondtype( pose.residue_type(ligid).bond_type( chiatoms[2], chiatoms[3] ) );

			core::Real angle_i = pose.chi( ichi_src, ligid );
			if ( bondtype == core::chemical::DoubleBond ) {
				// keep amide bonds at 0 or 180
				if ( numeric::random::rg().uniform() < (chi_step / 180.0) ) {
					angle_i = std::fmod( angle_i + 180.0, 360.0);
				}
			} else {
				angle_i = std::fmod( angle_i + chi_step * numeric::random::rg().gaussian(), 360.0);
			}
			pose.set_chi( ichi_src, ligid, angle_i );
		}
		if ( pose.residue_type(ligid).is_protein() ) {
			core::Real phi_i = std::fmod( pose.phi(ligid) + chi_step * numeric::random::rg().gaussian(), 360.0);
			core::Real psi_i = std::fmod( pose.psi(ligid) + chi_step * numeric::random::rg().gaussian(), 360.0);
			core::Real omega_i = pose.omega(ligid);
			if ( numeric::random::rg().uniform() < (chi_step / 180.0) ) {
				omega_i = std::fmod( omega_i + 0.1 * chi_step * numeric::random::rg().gaussian(), 360.0);
			}
			pose.set_phi( ligid, phi_i );
			pose.set_psi( ligid, psi_i );
			pose.set_omega( ligid, omega_i );
		}
	}
}

core::Size
get_ligand_jumpid(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const& ligids
){
	core::Size jumpid(0);
	utility::vector1< core::kinematics::Edge > jumps = pose.fold_tree().get_jump_edges();
	// skip if ligand-only case
	if ( jumps.size() >= 1 ) {
		for ( auto j : jumps ) {
			core::Size up = j.start(), down = j.stop();
			if ( std::find( ligids.begin(), ligids.end(), up) == ligids.end()
					&& std::find( ligids.begin(), ligids.end(), down) != ligids.end()
					) {
				runtime_assert(jumpid==0);
				jumpid = j.label();
			}
		}
		runtime_assert(jumpid!=0);
	}
	return jumpid;
}

void
get_ligand_resids(core::pose::Pose const& pose,
	utility::vector1 < core::Size >& lig_resids)
{
	lig_resids.clear();
	core::Size lastres = pose.total_residue();
	while ( pose.residue(lastres).is_virtual_residue() && lastres>1 ) lastres--;
	lig_resids.push_back( lastres );
	bool extending = true;
	while ( extending ) {
		if ( !pose.residue(lastres).is_polymer() ) {
			extending = false; // ligand not a polymer, we're done
		} else if ( pose.residue(lastres).has_lower_connect() ) {
			core::Size nextres = pose.residue(lastres).connected_residue_at_resconn( pose.residue(lastres).type().lower_connect_id() );
			if ( std::find (lig_resids.begin(), lig_resids.end(), nextres) == lig_resids.end() ) {
				lastres = nextres;
				lig_resids.push_back( lastres );
			} else {
				extending = false;  // finished a cyclic peptide, we're done
			}
		} else {
			extending = false; // hit n-term
		}
	}
}

bool
is_hb_satisfied(core::scoring::ScoreFunctionOP sf, core::scoring::hbonds::HBondDatabaseCOP hb_database,
	core::scoring::hbonds::HBondOptions const & hbopt, hbAcc const& acc,
	hbDon const& don, core::Real const& maxHbdis2,
	core::Real const& hb_energy_cutoff, std::string const& metric)
{
	numeric::xyzVector<core::Real> const& D = don.D;
	numeric::xyzVector<core::Real> const& H = don.H;
	core::scoring::hbonds::HBEvalTuple hbt;
	hbt.don_type( don.dontype );
	numeric::xyzVector<core::Real> const& A = acc.A;
	numeric::xyzVector<core::Real> const& B = acc.B;
	numeric::xyzVector<core::Real> const& B_0 = acc.B_0;
	hbt.acc_type( acc.acctype );

	core::Real dist2 = (A-D).length_squared();
	if ( metric == "simple" ) {
		if ( dist2 > maxHbdis2 ) return false;
		return true;
	} else {
		if ( dist2 > core::scoring::hbonds::MAX_R2 ) return false;
		core::Real hb_energy = get_hbond_score_weighted(sf, hb_database, hbopt, hbt, D, H, A, B, B_0);
		if ( hb_energy < hb_energy_cutoff ) return true;
	}
	return false;
}

void
compute_nhbonds(core::pose::Pose const& pose,
	utility::vector1<core::Size> const& ligids,
	utility::vector1<core::Size> const& resids,
	core::Size & nhbonds_total, core::Size & nhbonds_max1,
	bool const& include_bb, std::string const& hb_metric)
{

	utility::vector1<hbDon> lig_hbDons, rec_hbDons;
	utility::vector1<hbAcc> lig_hbAccs, rec_hbAccs;

	for ( auto ligid : ligids ) {
		core::conformation::Residue const& ires = pose.residue(ligid);
		// ligand acceptors
		for ( auto anum=ires.accpt_pos().begin(), anume=ires.accpt_pos().end(); anum!=anume; ++anum ) {
			hbAcc acc;
			core::Size bnum = ires.atom_base( *anum );
			core::Size b0num = ires.abase2( *anum );
			acc.A =  ires.xyz( *anum );
			acc.B = ires.xyz( bnum );
			acc.B_0 = ires.xyz( b0num );
			acc.acctype = core::scoring::hbonds::get_hb_acc_chem_type( *anum, ires );
			lig_hbAccs.push_back(acc);
		}
		// ligand donors
		for ( auto hnum=ires.Hpos_polar().begin(), hnume=ires.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			hbDon don;
			core::Size dnum = ires.atom_base( *hnum );
			don.H =  ires.xyz( *hnum );
			don.D = ires.xyz( dnum );
			don.dontype = core::scoring::hbonds::get_hb_don_chem_type( dnum, ires );
			lig_hbDons.push_back(don);
		}
	}

	for ( auto resid : resids ) {
		core::conformation::Residue const& ires = pose.residue(resid);
		// receptor acceptors
		for ( auto anum=ires.accpt_pos().begin(), anume=ires.accpt_pos().end(); anum!=anume; ++anum ) {
			if ( !include_bb && !(*anum>ires.last_backbone_atom() && *anum<=ires.nheavyatoms()) ) continue;
			hbAcc acc;
			core::Size bnum = ires.atom_base( *anum );
			core::Size b0num = ires.abase2( *anum );
			acc.A =  ires.xyz( *anum );
			acc.B = ires.xyz( bnum );
			acc.B_0 = ires.xyz( b0num );
			acc.acctype = core::scoring::hbonds::get_hb_acc_chem_type( *anum, ires );
			rec_hbAccs.push_back(acc);
		}
		// receptor donors
		for ( auto hnum=ires.Hpos_polar().begin(), hnume=ires.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			if ( !include_bb && !(*hnum>ires.first_sidechain_hydrogen()) ) continue;
			hbDon don;
			core::Size dnum = ires.atom_base( *hnum );
			don.H =  ires.xyz( *hnum );
			don.D = ires.xyz( dnum );
			don.dontype = core::scoring::hbonds::get_hb_don_chem_type( dnum, ires );
			rec_hbDons.push_back(don);
		}
	}
	nhbonds_total = 0;
	core::Real maxHbdis = 4.2;
	core::Real maxHbdis2 = maxHbdis*maxHbdis;
	core::Real hb_energy_cutoff = -0.3;
	utility::vector1< bool > rec_hbAcc_satisfied( rec_hbAccs.size(), false );
	utility::vector1< bool > rec_hbDon_satisfied( rec_hbDons.size(), false );

	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
	core::scoring::hbonds::HBondOptions const & hbopt = sf->energy_method_options().hbond_options();
	core::scoring::hbonds::HBondDatabaseCOP hb_database = core::scoring::hbonds::HBondDatabase::get_database( hbopt.params_database_tag() );
	//compute the number of unsatisfied hbonds of receptor donors
	bool is_satisfied(false);
	for ( core::Size i=1; i<=rec_hbDons.size(); ++i ) {
		for ( auto lig_hbAcc : lig_hbAccs ) {
			is_satisfied = is_hb_satisfied(sf, hb_database, hbopt, lig_hbAcc, rec_hbDons[i], maxHbdis2, hb_energy_cutoff, hb_metric );
			if ( is_satisfied ) {
				nhbonds_total++;
				if ( rec_hbDon_satisfied[i] ) continue;
				rec_hbDon_satisfied[i] = true;
				nhbonds_max1++;
			}
		}
	}
	//compute the number of unsatisfied hbonds of receptor acceptors
	for ( core::Size i=1; i<=rec_hbAccs.size(); ++i ) {
		for ( auto lig_hbDon : lig_hbDons ) {
			is_satisfied = is_hb_satisfied(sf, hb_database, hbopt, rec_hbAccs[i], lig_hbDon, maxHbdis2, hb_energy_cutoff, hb_metric );
			if ( is_satisfied ) {
				nhbonds_total++;
				if ( rec_hbAcc_satisfied[i] ) continue;
				rec_hbAcc_satisfied[i] = true;
				nhbonds_max1++;
			}
		}
	}
}

} // ga_dock
} // ligand_docking
} // protocols


