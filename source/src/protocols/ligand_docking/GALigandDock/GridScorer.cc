// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GridScorer.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/GriddedAtomTreeMultifunc.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/conformation/residue_datacache.hh>
#include <core/pack/optimizeH.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/constants.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <protocols/relax/FastRelax.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>

#include <basic/options/util.hh> // HACK
#include <basic/options/option.hh> // HACK
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

//////////////////////////////
namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

//using basic::T;
using basic::Tracer;
static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.GridScorer" );

static core::Real OUT_OF_BOUND_E = 100;

// fast integer powers
inline
core::Real ipow(core::Real base, core::Size exp) {
	core::Real result = 1;
	while ( exp ) {
		if ( exp & 1 ) {
			result *= base;
		}
		exp >>= 1;
		base *= base;
	}
	return result;
}

// nonlinear xforms
inline core::Real xform_rep(core::Real rep) {
	return ( std::pow(rep+1, -1.0/4.0) );
	//return rep;
}
inline core::Real ixform_rep(core::Real xfrep) {
	return (1.0/ipow(xfrep, 4)-1);
	//return xfrep;
}
inline core::Real dxform_rep(core::Real xfrep) {
	return (-4.0/ipow(xfrep, 5));
	//return 1.0;
}

// nonlinear xforms
inline core::Real xform_atr(core::Real atr) {
	return atr;
}
inline core::Real ixform_atr(core::Real xfatr) {
	return xfatr;
}
inline core::Real dxform_atr(core::Real /*xfatr*/) {
	return 1.0;
}

// helper function for getting hbond energies
inline
core::Real
get_hbond_score_weighted (
	core::scoring::ScoreFunctionOP sf,
	core::scoring::hbonds::HBondDatabaseCOP database,
	core::scoring::hbonds::HBondOptions const & hbopt,
	core::scoring::hbonds::HBEvalTuple const &hbt,
	numeric::xyzVector<core::Real> const &D,
	numeric::xyzVector<core::Real> const &H,
	numeric::xyzVector<core::Real> const &A,
	numeric::xyzVector<core::Real> const &B,
	numeric::xyzVector<core::Real> const &B0
) {
	core::Real energy=0;
	core::scoring::hbonds::hb_energy_deriv(
		*database, hbopt, hbt,
		D,H,A,B,B0,
		energy, false, core::scoring::hbonds::DUMMY_DERIVS
	);
	if ( energy<core::scoring::hbonds::MAX_HB_ENERGY ) {
		switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
		case core::scoring::hbonds::hbw_SR_BB_SC :
		case core::scoring::hbonds::hbw_LR_BB_SC :
			energy *= sf->get_weight( core::scoring::hbond_bb_sc );
			break;
		case core::scoring::hbonds::hbw_SC :
		default :
			energy *= sf->get_weight( core::scoring::hbond_sc );
		}
		return energy;
	}
	return 0.0;
}


GridScorer::GridScorer( core::scoring::ScoreFunctionOP sfxn ) :
	etable_( *core::scoring::ScoringManager::get_instance()->etable( sfxn->energy_method_options() ).lock() )
{
	// set up scorefunctions
	sfxn_ = sfxn->clone();
	sfxn_soft_ = sfxn->clone();
	sfxn_clash_ = core::scoring::ScoreFunctionFactory::create_score_function("empty.wts");
	sfxn_1b_ = core::scoring::ScoreFunctionFactory::create_score_function("empty.wts");
	sfxn_1b_soft_ = core::scoring::ScoreFunctionFactory::create_score_function("empty.wts");
	sfxn_1b_clash_ = core::scoring::ScoreFunctionFactory::create_score_function("empty.wts");


	// copy options!
	core::scoring::methods::EnergyMethodOptions e_opts = sfxn_->energy_method_options();
	sfxn_clash_->set_energy_method_options(e_opts);
	sfxn_1b_->set_energy_method_options(e_opts);
	sfxn_1b_soft_->set_energy_method_options(e_opts);
	sfxn_1b_clash_->set_energy_method_options(e_opts);

	// set up scorefunction-derived parameters
	LKBe_ = core::scoring::lkball::LK_BallEnergyOP( new core::scoring::lkball::LK_BallEnergy(e_opts) );
	coulomb_ = core::scoring::etable::coulomb::CoulombOP ( new core::scoring::etable::coulomb::Coulomb(e_opts) );
	coulomb_->initialize();
	core::scoring::hbonds::HBondOptions const & hbopt = e_opts.hbond_options();
	hb_database_ = core::scoring::hbonds::HBondDatabase::get_database( hbopt.params_database_tag() );

	core::Real soft_rep_weight( 0.001 ), soft_intrarep_weight( 0.01 ), soft_gen_bonded_weight( 0.1 ); //run9 params

	// modify sf soft
	if ( (*sfxn_)[ core::scoring::fa_rep ] > 0.0 ) {
		sfxn_soft_->set_weight( core::scoring::fa_rep, soft_rep_weight );
	}
	if ( (*sfxn_)[ core::scoring::fa_intra_rep ] > 0.0 ) {
		sfxn_soft_->set_weight( core::scoring::fa_intra_rep, soft_intrarep_weight );
	}
	if ( (*sfxn_)[ core::scoring::fa_intra_rep_xover4 ] > 0.0 ) {
		sfxn_soft_->set_weight( core::scoring::fa_intra_rep_xover4, soft_intrarep_weight );
	}
	if ( (*sfxn_)[ core::scoring::gen_bonded ] > soft_gen_bonded_weight ) {
		sfxn_soft_->set_weight( core::scoring::gen_bonded, soft_gen_bonded_weight );
	}
	// post-run9g: also increase fa_elec & hbond strength
	core::Real w_fa_elec = (*sfxn_soft_)[ core::scoring::fa_elec ];
	core::Real w_hbond_sc = (*sfxn_soft_)[ core::scoring::hbond_sc ];
	core::Real w_hbond_bb_sc = (*sfxn_soft_)[ core::scoring::hbond_bb_sc ];
	sfxn_soft_->set_weight( core::scoring::fa_elec, w_fa_elec*2.0 );
	sfxn_soft_->set_weight( core::scoring::hbond_sc, w_hbond_sc*5.0 );
	sfxn_soft_->set_weight( core::scoring::hbond_bb_sc, w_hbond_bb_sc*5.0 );

	// modify sf clash
	if ( (*sfxn_)[ core::scoring::fa_rep ] > 0.0 ) {
		sfxn_clash_->set_weight( core::scoring::fa_rep, 0.01 );
	}
	if ( (*sfxn_)[ core::scoring::fa_intra_rep ] > 0.0 ) {
		sfxn_clash_->set_weight( core::scoring::fa_intra_rep, 0.001 );
	}

	// modify sfxn 1body
	sfxn_1b_->set_weight( core::scoring::fa_dun, sfxn_->get_weight( core::scoring::fa_dun ) );
	sfxn_1b_->set_weight( core::scoring::fa_dun_dev, sfxn_->get_weight( core::scoring::fa_dun_dev ) );
	sfxn_1b_->set_weight( core::scoring::fa_dun_rot, sfxn_->get_weight( core::scoring::fa_dun_rot ) );
	sfxn_1b_->set_weight( core::scoring::fa_dun_semi, sfxn_->get_weight( core::scoring::fa_dun_semi ) );
	sfxn_1b_->set_weight( core::scoring::fa_intra_rep, sfxn_->get_weight( core::scoring::fa_intra_rep ) );
	sfxn_1b_->set_weight( core::scoring::fa_intra_atr_xover4, sfxn_->get_weight( core::scoring::fa_intra_atr_xover4 ) );
	sfxn_1b_->set_weight( core::scoring::fa_intra_rep_xover4, sfxn_->get_weight( core::scoring::fa_intra_rep_xover4 ) );
	sfxn_1b_->set_weight( core::scoring::fa_intra_sol_xover4, sfxn_->get_weight( core::scoring::fa_intra_sol_xover4 ) );
	sfxn_1b_->set_weight( core::scoring::fa_intra_elec, sfxn_->get_weight( core::scoring::fa_intra_elec ) );
	sfxn_1b_->set_weight( core::scoring::gen_bonded, sfxn_->get_weight( core::scoring::gen_bonded ) );

	// modify sfxn 1body soft
	sfxn_1b_soft_->set_weight( core::scoring::fa_dun, sfxn_soft_->get_weight( core::scoring::fa_dun ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_dun_dev, sfxn_soft_->get_weight( core::scoring::fa_dun_dev ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_dun_rot, sfxn_soft_->get_weight( core::scoring::fa_dun_rot ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_dun_semi, sfxn_soft_->get_weight( core::scoring::fa_dun_semi ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_intra_rep, sfxn_soft_->get_weight( core::scoring::fa_intra_rep ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_intra_atr_xover4, sfxn_soft_->get_weight( core::scoring::fa_intra_atr_xover4 ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_intra_rep_xover4, sfxn_soft_->get_weight( core::scoring::fa_intra_rep_xover4 ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_intra_sol_xover4, sfxn_soft_->get_weight( core::scoring::fa_intra_sol_xover4 ) );
	sfxn_1b_soft_->set_weight( core::scoring::fa_intra_elec, sfxn_soft_->get_weight( core::scoring::fa_intra_elec ) );
	sfxn_1b_soft_->set_weight( core::scoring::gen_bonded, sfxn_soft_->get_weight( core::scoring::gen_bonded ) );

	// modify sfxn 1body clash
	sfxn_1b_clash_->set_weight( core::scoring::fa_intra_rep, sfxn_soft_->get_weight( core::scoring::fa_intra_rep ) );
	sfxn_1b_clash_->set_weight( core::scoring::fa_intra_rep_xover4, sfxn_soft_->get_weight( core::scoring::fa_intra_rep_xover4 ) );

	bbox_padding_ = 5.0;
	voxel_spacing_ = 0.5;
	hash_grid_ = 6.0;
	hash_subgrid_ = 1;

	maxdis_ = 0; // initialized in build_grid

	debug_ = exact_ = false;

	//LK_fade_ = 0.01; // ? do a softmin with lk terms
	LK_fade_ = 1.0; // ? do a softmin with lk terms

	dims_ = numeric::xyzVector< core::Size >(0,0,0);

	maxiter_minimize_ = 100; // maxiter for minimization

	w_rep_ = 1.0;
	smoothing_ = 0.0;
	packer_cycles_ = 0;

	pack_time_ = min_time_ = std::chrono::duration<double>{};
}

GridScorer::~GridScorer(){}

void
GridScorer::prepare_grid( core::pose::Pose const &pose, core::Size const lig_resid ) {
	// define bounding box
	core::conformation::Residue const &resLig = pose.residue( lig_resid );
	numeric::xyzVector< core::Real > com(0,0,0);
	for ( core::Size i=1; i<=resLig.natoms(); ++i ) com += resLig.xyz(i);
	com /= resLig.natoms();
	core::Real maxRad = 0.0;
	for ( core::Size i=1; i<=resLig.natoms(); ++i ) maxRad = std::max( maxRad, com.distance(resLig.xyz(i)) );

	if ( maxRad < bbox_padding_*1.5 ) { maxRad = bbox_padding_*1.5; } // modification for very small ligands

	origin_ = com - (maxRad + bbox_padding_);
	dims_[0] = (core::Size)std::ceil( 2*(maxRad + bbox_padding_) / voxel_spacing_ );
	dims_[2] = dims_[1] = dims_[0];

	ref_pose_ = core::pose::PoseOP( new core::pose::Pose(pose) ); // reference pose (used in scoring)

	TR << "  ... built grid with " << maxRad << " + " << bbox_padding_ << " A radius; "
		<< voxel_spacing_ << " A grid spacing" << std::endl;
}

void
GridScorer::calculate_grid(
	core::pose::Pose const &pose,
	core::Size const lig_resid,
	utility::vector1< core::Size > const &movingSCs
) {
	auto start = std::chrono::steady_clock::now();

	// initialize the grid if it has not been already
	if ( dims_[2]*dims_[1]*dims_[0] == 0 ) prepare_grid( pose, lig_resid );

	core::pose::PoseOP trim_pose( new core::pose::Pose(pose) ); // pose with moving sidechains trimmed to GLY (used in grid calcs)

	core::chemical::AtomTypeSetCOP atom_set
		= core::chemical::ChemicalManager::get_instance()->atom_type_set(core::chemical::FA_STANDARD);

	// build grid
	(*sfxn_)(*trim_pose); // make sure all pose-derived info initialized

	// compress all five lkball weights into two tables
	// note this does not allow dynamic reweighing of lkball subterms
	core::Real weightLKb=sfxn_->get_weight( core::scoring::lk_ball );
	core::Real weightLKbi=sfxn_->get_weight( core::scoring::lk_ball_iso );
	core::Real weightLKbr=sfxn_->get_weight( core::scoring::lk_ball_bridge );
	core::Real weightLKbru=sfxn_->get_weight( core::scoring::lk_ball_bridge_uncpl );
	bool useLKB = ( weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 );
	bool useLKBr = ( weightLKbr != 0 || weightLKbru != 0 );

	if ( sfxn_->get_weight( core::scoring::lk_ball_wtd ) != 0.0 ) {
		TR << "ERROR! The grid-based scheme does not work with lk_ball_wtd.  Please use beta_nov16 or later." << std::endl;
		utility_exit();
	}

	core::scoring::lkball::LK_BallEnergy lkballE( sfxn_->energy_method_options() );
	maxdis_ = std::max( etable_.max_dis(), coulomb_->max_dis( ) );
	if ( useLKB ) maxdis_ = std::max( maxdis_, lkballE.lkb_max_dis( ) );

	core::Real subhash_buffer = sqrt(3.0)/2.0 * (hash_subgrid_ - 1) * voxel_spacing_;

	TR << "Building " << dims_[0] << " x " << dims_[1] << " x " << dims_[2] << " grid; "
		<< "origin = ( " << origin_[0] << "," << origin_[1] << "," << origin_[2] << " )" << std::endl;
	TR << "   maxdis = " << maxdis_ << std::endl;
	TR << "   subhash_buffer = " << subhash_buffer << std::endl;
	TR << "   hash_grid = " << hash_grid_ << std::endl;

	// build 3d hashmaps
	//    one for all heavyatoms for etable gridding, only used in this function
	//    two for hb donors/acceptors used in scoring
	GridHash3D<atmInfo> gridhash;
	gridhash.set_resolution( hash_grid_ );
	hbdonors_.set_resolution( hash_grid_ );
	hbacceptors_.set_resolution( hash_grid_ );

	for ( core::Size b_res=1; b_res<=pose.total_residue(); ++b_res ) {
		if ( b_res == lig_resid ) continue;

		bool is_moving = (std::find (movingSCs.begin(),movingSCs.end(), b_res) != movingSCs.end());

		core::conformation::Residue const &resB = trim_pose->residue( b_res );

		// 0.1 : hbonds
		for ( auto anum=resB.accpt_pos().begin(),anume=resB.accpt_pos().end(); anum!=anume; ++anum ) {
			// add bb hbonds from moving residues
			if ( is_moving && !(*anum>resB.last_backbone_atom() && *anum<=resB.nheavyatoms()) ) continue;

			hbAcc acc_b;
			core::Size bnum = resB.atom_base( *anum );
			core::Size b0num = resB.abase2( *anum  );
			acc_b.A = resB.xyz( *anum );
			acc_b.B = resB.xyz( bnum );
			acc_b.B0 = resB.xyz( b0num );
			acc_b.acctype = core::scoring::hbonds::get_hb_acc_chem_type( *anum, resB );
			hbacceptors_.add_point( acc_b );
		}

		for ( auto hnum=resB.Hpos_polar().begin(),hnume=resB.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			// add bb hbonds from moving residues
			if ( is_moving && !(*hnum>resB.first_sidechain_hydrogen()) ) continue;

			hbDon don_b;
			core::Size dnum = resB.atom_base( *hnum );
			don_b.H = resB.xyz( *hnum );
			don_b.D = resB.xyz( dnum );
			don_b.dontype = core::scoring::hbonds::get_hb_don_chem_type( dnum, resB );
			hbdonors_.add_point( don_b ); //? or H?
		}

		if ( is_moving ) continue;

		// 0.2 : grid interactions
		utility::vector1< core::Size > rsdBn_waters;
		utility::vector1< core::scoring::lkball::WaterCoords > rsdBwaters;

		if ( useLKB ) {
			core::scoring::lkball::LKB_ResidueInfo const &lkbinfo =
				static_cast< core::scoring::lkball::LKB_ResidueInfo const & > (
				resB.data_ptr()->get( core::conformation::residue_datacache::LK_BALL_INFO ));
			rsdBn_waters = lkbinfo.n_attached_waters();
			rsdBwaters = lkbinfo.waters();
		}

		for ( core::Size b_atm=1; b_atm<=resB.natoms(); ++b_atm ) {
			atmInfo info_b;
			info_b.atom = resB.atom(b_atm);
			info_b.atomicCharge = resB.atomic_charge(b_atm);

			if ( b_atm <= rsdBwaters.size() && rsdBwaters[b_atm].size() > 0 ) {
				info_b.atomWaters = rsdBwaters[b_atm];
				info_b.nAtomWaters = rsdBn_waters[b_atm];
			} else {
				info_b.nAtomWaters = 0;
			}
			gridhash.add_point( info_b );
		}

	}

	// find all unique atoms & initialize grids
	// find don/acc _ligand_ heavyatoms
	core::Size nSCs = movingSCs.size();
	utility::vector1< utility::vector1<bool> > hbdon_acc( nSCs+1 );
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig_resid : movingSCs[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		hbdon_acc[i+1].resize( res_i.natoms(), false );
		for ( auto anum=res_i.accpt_pos().begin(),anume=res_i.accpt_pos().end(); anum!=anume; ++anum ) {
			hbdon_acc[i+1][*anum] = true;
		}
		for ( auto hnum=res_i.Hpos_polar().begin(),hnume=res_i.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			core::Size datm = res_i.atom_base( *hnum );
			hbdon_acc[i+1][datm] = true;
		}
	}

	std::map< int, core::conformation::Atom > uniq_atoms; // KEEP AS INT, IT IS NOT A SIZE
	core::Size nLK=0;
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig_resid : movingSCs[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );

		for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
			bool j_is_backbone = (i!=0) && (
				j<=res_i.last_backbone_atom() || (j>res_i.nheavyatoms() && j<res_i.first_sidechain_hydrogen()) );

			if ( j_is_backbone ) continue;

			int atmtype_ij = res_i.atom(j).type();
			if ( uniq_atoms.find(atmtype_ij) == uniq_atoms.end() ) {
				uniq_atoms[atmtype_ij] = res_i.atom(j);

				raw_faatr_[ atmtype_ij ] = ObjexxFCL::FArray3D< float >();
				raw_faatr_[ atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);

				raw_farep_[ atmtype_ij ] = ObjexxFCL::FArray3D< float >();
				raw_farep_[ atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);

				raw_fasol_[ atmtype_ij ] = ObjexxFCL::FArray3D< float >();
				raw_fasol_[ atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);

				if ( !useLKB ) continue;

				//fd we need this even for nonpolars (since they will desolvate polars through lkb)
				raw_lkball_[ atmtype_ij ] = ObjexxFCL::FArray3D< float >();
				raw_lkball_[ atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);

				if ( hbdon_acc[i+1][j] ) {
					raw_lkball_[ -atmtype_ij ] = ObjexxFCL::FArray3D< float >();
					raw_lkball_[ -atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);
					if ( useLKBr ) {
						raw_lkbridge_[ atmtype_ij ] = ObjexxFCL::FArray3D< float >();
						raw_lkbridge_[ atmtype_ij ].dimension(dims_[0],dims_[1],dims_[2], 0.0);
					}
					nLK++;
				}
			}
		}
	}
	raw_faelec_.dimension(dims_[0],dims_[1],dims_[2], 0.0);

	// report uniq atom index
	TR.Debug << "atom type: " << std::endl;
	for ( auto it=uniq_atoms.begin(); it != uniq_atoms.end(); ++it ) {
		core::conformation::Atom a_i = it->second;
		TR.Debug << it->first << ": " << (*atom_set)[a_i.type()].name() << std::endl;
	}

	TR << "Calculating gridded energies for " << uniq_atoms.size() << " unique atom types" << std::endl;
	TR << "  and lkb virtual parameters for " << nLK << " unique polar atom types" << std::endl;
	utility::vector1<atmInfo> neighborlist;

	// lkball setup stuff
	core::Real heavyWatDist = 2.65;
	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].user() ) {
		heavyWatDist = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ]()[1];
	}
	core::Real watWatDist = heavyWatDist * 4*sqrt(2)/3; // assume perfect tetrahedral geom
	core::Real ramp_width_A2=lkballE.get_ramp_width_A2();
	core::Real overlap_width_A2=lkballE.get_overlap_width_A2();
	core::Real overlap_gap_A2=lkballE.get_overlap_gap2();
	core::Real d2_low1,d2_low2;

	for ( core::Size z=1; z<=dims_[2]; z+=hash_subgrid_ ) {
		for ( core::Size y=1; y<=dims_[1]; y+=hash_subgrid_ ) {
			for ( core::Size x=1; x<=dims_[0]; x+=hash_subgrid_ ) {
				numeric::xyzVector< core::Real > subhash_mdpt (
					voxel_spacing_*(x-1+(hash_subgrid_-1)/2.0),
					voxel_spacing_*(y-1+(hash_subgrid_-1)/2.0),
					voxel_spacing_*(z-1+(hash_subgrid_-1)/2.0)
				);
				gridhash.get_neighbors( subhash_mdpt+origin_, neighborlist, maxdis_ + subhash_buffer );

				for ( core::Size xo=0; xo<hash_subgrid_; xo++ ) {
					if ( x+xo > dims_[0] ) continue;
					for ( core::Size yo=0; yo<hash_subgrid_; yo++ ) {
						if ( y+yo > dims_[1] ) continue;
						for ( core::Size zo=0; zo<hash_subgrid_; zo++ ) {
							if ( z+zo > dims_[2] ) continue;

							numeric::xyzVector< core::Real > cart_xyz (
								voxel_spacing_*(x+xo-1), voxel_spacing_*(y+yo-1), voxel_spacing_*(z+zo-1) );
							cart_xyz = cart_xyz+origin_;

							//////
							// 1) etable
							for ( auto it=uniq_atoms.begin(); it != uniq_atoms.end(); ++it ) {
								core::conformation::Atom a_i = it->second;
								a_i.xyz( cart_xyz );

								bool calc_virt_sites;
								calc_virt_sites = (raw_lkball_.find( -it->first ) != raw_lkball_.end());

								core::scoring::lkball::WaterCoords water_xyz;
								water_xyz[1] = cart_xyz;

								core::Real faatrSum=0, farepSum=0, fasolSum=0;
								core::Real lkballHeavySum=0, lkballWatSum=0, lkbridgeSum=0;

								for ( core::Size i_neigh=1; i_neigh<=neighborlist.size(); ++i_neigh ) {
									atmInfo const & neighborAtom = neighborlist[i_neigh];

									// a) fasol / farep / faatr
									core::Real faatr,farep,fasol1,fasol2;
									fast_eval_etable_split_fasol( etable_, a_i, neighborAtom.atom, faatr, farep, fasol1, fasol2 );

									faatrSum += faatr;
									farepSum += farep;
									fasolSum += fasol1+fasol2;

									if ( !useLKB ) continue;

									// b) lk* part 1: occlusion _by_ a heavyatom at xyz
									// atom 1 = "grid", atom2 = "background"
									core::Size nWaters_i = neighborAtom.nAtomWaters;
									core::scoring::lkball::WaterCoords const &waters_i = neighborAtom.atomWaters;
									d2_low1 = lkballE.get_d2_low( a_i.type() );
									if ( ( nWaters_i > 0) ) {
										core::Real fasol2_lkball =
											fasol2 * fast_get_lk_fractional_contribution( cart_xyz, nWaters_i, waters_i, ramp_width_A2, d2_low1, 0.0 );
										lkballHeavySum += weightLKb * ( fasol2_lkball );
										lkballHeavySum += weightLKbi * ( fasol2 );
									}

									// c) lk* part 2: occlusion of (or bridging to) a virt at xyz
									if ( calc_virt_sites ) {
										lkballHeavySum += weightLKbi * ( fasol1 );

										// now the estimate ...
										// we know lk_frac exactly but need to estimate fasol
										// we will assume occusion is directly over the VRT H2O
										core::conformation::Atom b_i = neighborAtom.atom;
										b_i.xyz( numeric::xyzVector<core::Real> ( cart_xyz[0]+heavyWatDist, cart_xyz[1], cart_xyz[2] ) );
										etable_.analytic_lk_energy( a_i, b_i, fasol1, fasol2 );

										d2_low2 = lkballE.get_d2_low( neighborAtom.atom.type() );
										core::Real fasol1_lkball =
											fasol1 * fast_get_lk_fractional_contribution( neighborAtom.atom.xyz(), 1, water_xyz, ramp_width_A2, d2_low2, 0.0 );
										lkballWatSum += weightLKb * ( fasol1_lkball );

										// finally lkbridge
										if ( useLKBr &&  ( nWaters_i > 0) ) {
											b_i.xyz( numeric::xyzVector<core::Real> ( cart_xyz[0]+watWatDist, cart_xyz[1], cart_xyz[2] ) );
											etable_.analytic_lk_energy( a_i, b_i, fasol1, fasol2 );
											core::Real lkbrfrac = fast_get_lkbr_fractional_contribution( 1, nWaters_i, water_xyz, waters_i,
												overlap_gap_A2, overlap_width_A2, 0.0 );
											lkbridgeSum += (weightLKbr*(fasol1+fasol2) + weightLKbru) * ( lkbrfrac );
										}
									}
								} // i_neigh

								raw_faatr_[ it->first ](x+xo,y+yo,z+zo) = float(faatrSum);
								raw_farep_[ it->first ](x+xo,y+yo,z+zo) = float(farepSum);
								raw_fasol_[ it->first ](x+xo,y+yo,z+zo) = float(fasolSum);
								if ( lkballHeavySum != 0 ) {
									raw_lkball_[ it->first ](x+xo,y+yo,z+zo) = float(lkballHeavySum);
								}
								if ( lkballWatSum != 0 ) {
									raw_lkball_[ -it->first ](x+xo,y+yo,z+zo) = float(lkballWatSum);
								}
								if ( lkbridgeSum != 0 ) {
									raw_lkbridge_[ it->first ](x+xo,y+yo,z+zo) = float(lkbridgeSum);
								}
							}

							//////
							// 2) fa_elec
							core::Real faelecSum=0;
							for ( core::Size i_neigh=1; i_neigh<=neighborlist.size(); ++i_neigh ) {
								faelecSum += coulomb_->eval_atom_atom_fa_elecE(cart_xyz, 1.0, neighborlist[i_neigh].atom.xyz(), neighborlist[i_neigh].atomicCharge );
							}
							raw_faelec_(x+xo,y+yo,z+zo) = float(faelecSum);
						}
					}
				}
			}
		}
	}

	// transform data so spline interp is in a different space
	for ( core::Size z=1; z<=dims_[2]; z+=1 ) {
		for ( core::Size y=1; y<=dims_[1]; y+=1 ) {
			for ( core::Size x=1; x<=dims_[0]; x+=1 ) {
				//int iatm(0);
				for ( auto it=uniq_atoms.begin(); it != uniq_atoms.end(); ++it ) {
					float farep = raw_farep_[ it->first ](x,y,z);
					raw_farep_[ it->first ](x,y,z) = xform_rep(farep);

					float faatr = raw_faatr_[ it->first ](x,y,z);
					raw_faatr_[ it->first ](x,y,z) = xform_atr(faatr);
				}
			}
		}
	}

	for ( auto it=uniq_atoms.begin(); it != uniq_atoms.end(); ++it ) {
		int atmtype = (int)it->first;

		do_convolution_and_compute_coeffs( raw_faatr_[ atmtype ], coeffs_faatr_[ atmtype ], 0.0 );
		do_convolution_and_compute_coeffs( raw_farep_[ atmtype ], coeffs_farep_[ atmtype ], 0.0 );
		do_convolution_and_compute_coeffs( raw_fasol_[ atmtype ], coeffs_fasol_[ atmtype ], 0.0 );

		if ( raw_lkball_.find( atmtype ) != raw_lkball_.end() ) {
			do_convolution_and_compute_coeffs( raw_lkball_[ atmtype ], coeffs_lkball_[ atmtype ], 0.0 );
		}
		if ( raw_lkball_.find( -atmtype ) != raw_lkball_.end() ) {
			do_convolution_and_compute_coeffs( raw_lkball_[ -atmtype ], coeffs_lkball_[ -atmtype ], 0.0 );
		}
		if ( raw_lkbridge_.find( atmtype ) != raw_lkbridge_.end() ) {
			do_convolution_and_compute_coeffs( raw_lkbridge_[ atmtype ], coeffs_lkbridge_[ atmtype ], 0.0 );
		}
	}
	do_convolution_and_compute_coeffs( raw_faelec_, coeffs_faelec_, 0.0 );

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
	TR << "Grid calculation took " << (diff).count() << " seconds." << std::endl;
}

// check to see all rotamers lie within grid
bool
GridScorer::is_residue_in_grid(
	core::conformation::Residue const &res,
	core::Real angle_buffer,
	core::Real edge_buffer
) {
	using namespace ObjexxFCL::format;

	// check 1: is the residue CB not too close to the edge?
	numeric::xyzVector< core::Real > midpoint( (dims_[0]-1.0)/2.0, (dims_[1]-1.0)/2.0, (dims_[2]-1.0)/2.0 );
	midpoint = origin_ + midpoint * voxel_spacing_;
	core::Real maxradius = 0.5*(dims_[0]-1.0)*voxel_spacing_ - edge_buffer;

	numeric::xyzVector< core::Real > capos = res.xyz( "CA" );
	numeric::xyzVector< core::Real > cbpos = res.xyz( res.nbr_atom() );
	core::Real dist = (cbpos - midpoint).length();
	if ( dist > maxradius ) return false; // CB too far from center

	// check 2: is the residue pointing "inward"?
	//    at closer than inner radius, be most permissive
	//    at max radius, be most strict
	if ( angle_buffer>0 ) {
		core::Real innerradius = maxradius - (bbox_padding_+edge_buffer);
		core::Real fracdist = std::min( 1.0, std::max( 0.0, (dist-innerradius)/(bbox_padding_+edge_buffer) ) );
		core::Real minangle=0, maxangle=angle_buffer;
		core::Real angle_tgt = minangle+(maxangle-minangle)*fracdist;
		core::Real angle = numeric::angle_degrees( capos, cbpos, midpoint );
		if ( angle < angle_tgt ) {
			TR << "residue " << res.seqpos() << " fail angle check (" << angle << "<" <<  angle_tgt << ")" << std::endl;
			return false;
		}
	}

	// check 3: are all atoms in the box?
	for ( core::Size i=res.first_sidechain_atom(); i<=res.natoms(); ++i ) {
		if ( i>res.nheavyatoms() && i<res.first_sidechain_hydrogen() ) continue;
		dist = (res.xyz( i ) - midpoint).length();
		if ( dist > maxradius ) {
			TR << "residue " << res.seqpos() << " fail placement check" << std::endl;
			return false;
		}
	}

	return true;
}

bool
GridScorer::is_residue_in_grid(
	core::conformation::Residue const &res,
	core::Real edge_buffer,
	numeric::xyzVector< core::Real > const &eigval,
	numeric::xyzMatrix< core::Real > const &eigvec
) {
	using namespace ObjexxFCL::format;

	// check 1: does the residue neighbor-sphere overlap ellipsoid defined by input?
	numeric::xyzVector< core::Real > midpoint( (dims_[0]-1.0)/2.0, (dims_[1]-1.0)/2.0, (dims_[2]-1.0)/2.0 );
	midpoint = origin_ + midpoint * voxel_spacing_;
	core::Real maxradius = 0.5*(dims_[0]-1.0)*voxel_spacing_ - edge_buffer;
	numeric::xyzVector< core::Real > offset = res.xyz( res.nbr_atom() ) - midpoint;

	// scale for each eigendimension
	numeric::xyzVector< core::Real > eigvalN (
		1.0,
		std::sqrt( eigval[0] / eigval[1] ),
		std::sqrt( eigval[0] / eigval[2] ) );

	numeric::xyzVector< core::Real > offsetEigspace = eigvec.transposed()*offset;

	// proper intersection of nbr_radius sphere and eigenspace ellipse is difficult (and we don't need it)
	// approximate overlap by moving nbr_radius to the middle
	core::Real len = offsetEigspace.length();
	offsetEigspace *= (1 - (res.nbr_radius()/len));

	core::Real dEigenspace = numeric::xyzVector< core::Real > (
		offsetEigspace[0]*eigvalN[0],
		offsetEigspace[1]*eigvalN[1],
		offsetEigspace[2]*eigvalN[2]).length();

	if ( dEigenspace > maxradius ) {
		//TR << "residue " << res.seqpos() << " fail intersect " << dEigenspace << " "
		// << offsetEigspace.length() << " " << offsetEigspace[0]<<" "<<offsetEigspace[1]<<" "<<offsetEigspace[2]
		// << " " << eigvalN[0]<<" "<<eigvalN[1]<<" "<<eigvalN[2]
		// << std::endl;
		return false;
	}

	// check 2: are all atoms in the box?
	for ( core::Size i=res.first_sidechain_atom(); i<=res.natoms(); ++i ) {
		if ( i>res.nheavyatoms() && i<res.first_sidechain_hydrogen() ) continue;
		core::Real dist = (res.xyz( i ) - midpoint).length();
		if ( dist > maxradius ) {
			// TR << "residue " << res.seqpos() << " fails placement check (atom " << i << " d=" << dist << ")" << std::endl;
			return false;
		}
	}

	return true;
}


core::Real
GridScorer::move_to_boundary( numeric::xyzVector< core::Real > &idxX ) {
	core::Real penalty = 0.0;
	for ( core::Size dim=0; dim<3; ++dim ) {
		if ( idxX[dim]<0 ) {
			core::Real offset = idxX[dim];
			penalty += (OUT_OF_BOUND_E*voxel_spacing_)*offset*offset;
			idxX[dim] = 0;
		} else if ( idxX[dim]>dims_[dim]-1 ) {
			core::Real offset = (idxX[dim]-dims_[dim]+1);
			penalty += (OUT_OF_BOUND_E*voxel_spacing_)*offset*offset;
			idxX[dim] = dims_[dim]-1;
		}
	}
	return penalty;
}

core::Real
GridScorer::move_to_boundary( numeric::xyzVector< core::Real > &idxX, numeric::xyzVector< core::Real > &dpen_dx ) {
	core::Real penalty = 0.0;
	dpen_dx = numeric::xyzVector< core::Real >(0.0,0.0,0.0);
	for ( core::Size dim=0; dim<3; ++dim ) {
		if ( idxX[dim]<0 ) {
			core::Real offset = idxX[dim];
			penalty += (OUT_OF_BOUND_E*voxel_spacing_)*offset*offset;
			dpen_dx[dim] += 2*(OUT_OF_BOUND_E*voxel_spacing_)*offset;
			idxX[dim] = 0;
		} else if ( idxX[dim]>dims_[dim]-1 ) {
			core::Real offset = (idxX[dim]-dims_[dim]+1);
			penalty += (OUT_OF_BOUND_E*voxel_spacing_)*offset*offset;
			dpen_dx[dim] += 2*(OUT_OF_BOUND_E*voxel_spacing_)*offset;
			idxX[dim] = dims_[dim]-1;
		}
	}
	return penalty;
}

bool
GridScorer::set_smoothing( core::Real setting ) {
	if ( setting == smoothing_ ) return false;
	smoothing_=setting;

	TR << "Recalculating grids with smoothing of " << setting << std::endl;

	// farep
	for ( auto iter = raw_farep_.begin(); iter != raw_farep_.end(); ++iter ) {
		// invert (since high xformed farep = low farep)
		do_convolution_and_compute_coeffs( raw_farep_[ iter->first ] , coeffs_farep_[ iter->first ], setting, true );
	}

	// fd: the following unfortunately seem to hurt sampling
	// faatr
	//for (auto iter = raw_faatr_.begin(); iter != raw_faatr_.end(); ++iter) {
	// do_convolution_and_compute_coeffs( raw_faatr_[ iter->first ] , coeffs_faatr_[ iter->first ], setting );
	//}
	// fasol
	//for (auto iter = raw_fasol_.begin(); iter != raw_fasol_.end(); ++iter) {
	// do_convolution_and_compute_coeffs( raw_fasol_[ iter->first ] , coeffs_fasol_[ iter->first ], setting );
	//}
	// lkball
	//for (auto iter = raw_lkball_.begin(); iter != raw_lkball_.end(); ++iter) {
	// do_convolution_and_compute_coeffs( raw_lkball_[ iter->first ] , coeffs_lkball_[ iter->first ], setting );
	//}
	// lkbridge
	//for (auto iter = raw_lkbridge_.begin(); iter != raw_lkbridge_.end(); ++iter) {
	// do_convolution_and_compute_coeffs( raw_lkbridge_[ iter->first ] , coeffs_lkbridge_[ iter->first ], setting );
	//}
	//faelec
	//do_convolution_and_compute_coeffs( raw_faelec_ , coeffs_faelec_, setting );

	return true;
}

void
GridScorer::do_convolution_and_compute_coeffs(
	ObjexxFCL::FArray3D< float > const &rawdata,
	ObjexxFCL::FArray3D< float > &smoothed_coeffs,
	core::Real smoothing,
	bool inverted
) {
	ObjexxFCL::FArray3D< double > tempdata, tempcoeffs;

	int N = rawdata.u1()*rawdata.u2()*rawdata.u3();

	if ( smoothing > 0 ) {
		// v2: min over window
		int window_width = std::ceil( smoothing / voxel_spacing_ );

		ObjexxFCL::FArray3D< bool > voxel_mask( 2*window_width+1, 2*window_width+1, 2*window_width+1 , false );
		for ( int dk=-window_width; dk<=window_width; dk++ ) {
			for ( int dj=-window_width; dj<=window_width; dj++ ) {
				for ( int di=-window_width; di<=window_width; di++ ) {
					numeric::xyzVector< core::Real > X( di*voxel_spacing_, dj*voxel_spacing_, dk*voxel_spacing_ );
					voxel_mask( di+window_width+1, dj+window_width+1, dk+window_width+1 ) = (X.length_squared() <= smoothing*smoothing);
				}
			}
		}

		tempdata.dimension(rawdata.u1(),rawdata.u2(),rawdata.u3());
		for ( int i=0; i<N; ++i ) tempdata[i] = rawdata[i];

		for ( int k=1; k<=rawdata.u3(); k++ ) {
			for ( int j=1; j<=rawdata.u2(); j++ ) {
				for ( int i=1; i<=rawdata.u1(); i++ ) {
					for ( int dk=-window_width; dk<=window_width; dk++ ) {
						if ( k+dk <= 0 || k+dk > rawdata.u3() ) continue;
						for ( int dj=-window_width; dj<=window_width; dj++ ) {
							if ( j+dj <= 0 || j+dj > rawdata.u2() ) continue;
							for ( int di=-window_width; di<=window_width; di++ ) {
								if ( i+di <= 0 || i+di > rawdata.u1() ) continue;

								if ( voxel_mask( di+window_width+1, dj+window_width+1, dk+window_width+1 ) ) {
									if ( inverted ) {
										tempdata(i,j,k) = std::max( tempdata(i,j,k), (double)(rawdata(i+di,j+dj,k+dk)) );
									} else {
										tempdata(i,j,k) = std::min( tempdata(i,j,k), (double)(rawdata(i+di,j+dj,k+dk)) );
									}
								}
							}
						}
					}
				}
			}
		}
	} else {
		tempdata.dimension(rawdata.u1(),rawdata.u2(),rawdata.u3());
		for ( int i=0; i<N; ++i ) tempdata[i] = rawdata[i];
	}

	core::scoring::electron_density::spline_coeffs( tempdata, tempcoeffs, true);

	smoothed_coeffs.dimension(tempcoeffs.u1(),tempcoeffs.u2(),tempcoeffs.u3());
	for ( int i=0; i<N; ++i ) smoothed_coeffs[i] = float(tempcoeffs[i]);
}


// score
core::Real
GridScorer::score( LigandConformer const &lig, bool soft ) {
	lig.to_pose( ref_pose_ );
	return score( *ref_pose_, lig, soft );
}

// score
core::Real
GridScorer::score( core::pose::Pose &pose, LigandConformer const &lig, bool soft ) {
	core::Real score_grid = 0.0;

	// exact
	if ( exact_ ) {
		core::scoring::ScoreFunctionOP sf = soft? sfxn_soft_ : sfxn_;
		core::Real score_exact = (*sf)( *ref_pose_ );
		return score_exact;
	}

	core::Real weightLKb = sfxn_->get_weight( core::scoring::lk_ball );
	core::Real weightLKbi = sfxn_->get_weight( core::scoring::lk_ball_iso );
	core::Real weightLKbr = sfxn_->get_weight( core::scoring::lk_ball_bridge );
	core::Real weightLKbru = sfxn_->get_weight( core::scoring::lk_ball_bridge_uncpl );
	bool useLKB = ( weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 );

	core::Size nSCs = lig.moving_scs().size();
	utility::vector1< core::scoring::lkball::LKB_ResidueInfoOP > alllkbrinfo(nSCs+1);

	// 1-body energies
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		if ( useLKB ) {
			alllkbrinfo[i+1] = core::scoring::lkball::LKB_ResidueInfoOP( new core::scoring::lkball::LKB_ResidueInfo (res_i) );
		}
		ReweightableRepEnergy score_i = get_1b_energy( res_i, alllkbrinfo[i+1], soft );
		score_grid += score_i.score( w_rep_ );
	}

	// 2-body energies
	//    [brute force calculations for interaction graph ... should be reasonable as long as too many sidechains are not allowed to move]
	if ( nSCs > 0 ) {
		for ( core::Size i=0; i<=nSCs; ++i ) {
			core::Size resid_i = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
			core::conformation::Residue const &res_i = pose.residue( resid_i );

			for ( core::Size j=i+1; j<=nSCs; ++j ) {
				core::Size resid_j = lig.moving_scs()[j];
				core::conformation::Residue const &res_j = pose.residue( resid_j );

				ReweightableRepEnergy score_ij = get_2b_energy( res_i, alllkbrinfo[i+1], res_j, alllkbrinfo[j+1], soft );
				score_grid += score_ij.score( w_rep_ );
			}
		}
	}

	// debugging output
	if ( debug_ ) {
		using namespace ObjexxFCL::format;
		core::scoring::ScoreFunctionOP sf = soft? sfxn_soft_ : sfxn_;
		core::Real score_exact = (*sf)(pose);
		TR << "scores (grid/exact)" << " " << F( 8, 3, score_grid) << " / " << F( 8, 3, score_exact ) << std::endl;
	}

	return score_grid;
}


// one-body energies
ReweightableRepEnergy
GridScorer::get_1b_energy(
	core::conformation::Residue const &res_i,
	core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo,
	bool soft
) {
	core::scoring::ScoreFunctionOP sf = soft? sfxn_soft_ : sfxn_;
	core::scoring::ScoreFunctionOP sf1b = soft? sfxn_1b_soft_ : sfxn_1b_;

	ReweightableRepEnergy score_grid;

	core::Real weightatr = sf->get_weight( core::scoring::fa_atr );
	core::Real weightrep = sf->get_weight( core::scoring::fa_rep );
	core::Real weightsol = sf->get_weight( core::scoring::fa_sol );
	core::Real weightelec = sf->get_weight( core::scoring::fa_elec );
	core::Real weightLKb = sf->get_weight( core::scoring::lk_ball );
	core::Real weightLKbi = sf->get_weight( core::scoring::lk_ball_iso );
	core::Real weightLKbr = sf->get_weight( core::scoring::lk_ball_bridge );
	core::Real weightLKbru = sf->get_weight( core::scoring::lk_ball_bridge_uncpl );
	bool useLKB = ( weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 );
	bool useLKBr = ( weightLKbr != 0 || weightLKbru != 0 );

	// 1 etable grid
	for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
		bool j_is_backbone =
			( j<=res_i.last_backbone_atom() || (j>res_i.nheavyatoms() && j<res_i.first_sidechain_hydrogen()) );
		if ( j_is_backbone ) continue;

		int atmtype_ij = res_i.atom(j).type();
		numeric::xyzVector< core::Real > idxX = (res_i.xyz(j) - origin_) / voxel_spacing_;

		core::Real penalty = move_to_boundary(idxX);
		core::Real atr( 0.0 ), rep( 0.0 ), sol( 0.0 ), elec( 0.0 ), lkb( 0.0 );

		atr = core::scoring::electron_density::interp_spline(coeffs_faatr_[ atmtype_ij ], idxX, true) ;
		rep = core::scoring::electron_density::interp_spline(coeffs_farep_[ atmtype_ij ], idxX, true) ;
		sol = core::scoring::electron_density::interp_spline(coeffs_fasol_[ atmtype_ij ], idxX, true);
		elec = core::scoring::electron_density::interp_spline(coeffs_faelec_, idxX, true);
		lkb = 0.0;
		if ( useLKB ) {
			lkb = core::scoring::electron_density::interp_spline(coeffs_lkball_[ atmtype_ij ], idxX, true); // already weighted
		}

		// xform
		atr = ixform_atr( atr );
		rep = ixform_rep( rep );

		// penalty added as a separate scoreterm with weight 1
		score_grid.fa_rep_ += weightrep*rep;
		//score_grid.energy_ += penalty + weightatr*atr + weightsol*sol + weightelec*res_i.atomic_charge(j)*elec + lkb;
		score_grid.penalty_wtd_ += penalty;
		score_grid.fa_atr_wtd_ += weightatr*atr;
		score_grid.fa_sol_wtd_ += weightsol*sol;
		score_grid.fa_elec_wtd_ += weightelec*res_i.atomic_charge(j)*elec;
		score_grid.lk_ball_wtd_ += lkb;

		if ( j > lkbrinfo->n_attached_waters().size() || lkbrinfo->n_attached_waters()[j] == 0 ) continue;

		core::Real lk0=0, lkbr0=0, lk=0, lkbr=0;
		for ( core::Size k=1; k<=lkbrinfo->n_attached_waters()[j]; ++k ) {
			numeric::xyzVector< core::Real > idxXj = (lkbrinfo->waters()[j][k] - origin_) / voxel_spacing_;
			move_to_boundary(idxXj);

			core::Real lk_ijk = core::scoring::electron_density::interp_spline(coeffs_lkball_[ -atmtype_ij ], idxXj, true);
			lk0 += std::exp( lk_ijk / LK_fade_ ); // softmax

			if ( useLKBr ) {
				core::Real lkbr_ijk = core::scoring::electron_density::interp_spline(coeffs_lkbridge_[ atmtype_ij ], idxXj, true);
				// sign reversed
				lkbr0 += std::exp( -lkbr_ijk / LK_fade_ ); // softmin
			}
		}

		if ( lk0 > 0 ) lk = LK_fade_ * log( lk0 ); // softmax
		if ( lkbr0 > 0 ) lkbr = -LK_fade_ * log( lkbr0 ); // softmin

		//score_grid.energy_ += lk + lkbr;
		score_grid.lk_ball_wtd_ += lk + lkbr;
	}

	// 2 hbond
	utility::vector1<hbDon> hbdon_neighborlist;
	utility::vector1<hbAcc> hbacc_neighborlist;
	numeric::xyzVector<core::Real> D,H,A,B,B0;
	core::scoring::hbonds::HBondOptions const & hbopt = sf->energy_method_options().hbond_options();

	core::Real maxHbdis = 4.2;
	for ( auto anum=res_i.accpt_pos().begin(),anume=res_i.accpt_pos().end(); anum!=anume; ++anum ) {
		if ( !(*anum>res_i.last_backbone_atom() && *anum<=res_i.nheavyatoms()) ) continue;

		hbdonors_.get_neighbors( res_i.xyz( *anum ), hbdon_neighborlist, maxHbdis );
		core::Size bnum = res_i.atom_base( *anum );
		core::Size b0num = res_i.abase2( *anum );
		A = res_i.xyz( *anum );
		B = res_i.xyz( bnum );
		B0 = res_i.xyz( b0num );

		core::scoring::hbonds::HBEvalTuple hbt;
		hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_i ));

		for ( core::Size i=1; i<=hbdon_neighborlist.size(); ++i ) {
			D = hbdon_neighborlist[i].D;
			H = hbdon_neighborlist[i].H;
			hbt.don_type( hbdon_neighborlist[i].dontype );
			if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;
			//score_grid.energy_ += get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
			score_grid.hbond_wtd_ += get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
		}
	}
	for ( auto hnum=res_i.Hpos_polar().begin(),hnume=res_i.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
		if ( !(*hnum>res_i.first_sidechain_hydrogen()) ) continue;

		core::Size dnum = res_i.atom_base( *hnum );
		hbacceptors_.get_neighbors( res_i.xyz( dnum ), hbacc_neighborlist, maxHbdis );
		H = res_i.xyz( *hnum );
		D = res_i.xyz( dnum );

		core::scoring::hbonds::HBEvalTuple hbt;
		hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_i ));

		for ( core::Size i=1; i<=hbacc_neighborlist.size(); ++i ) {
			A = hbacc_neighborlist[i].A;
			B = hbacc_neighborlist[i].B;
			B0 = hbacc_neighborlist[i].B0;
			hbt.acc_type( hbacc_neighborlist[i].acctype );

			if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;
			//score_grid.energy_ += get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
			score_grid.hbond_wtd_ += get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
		}
	}

	// 3 internal energy for ligands / moving SCs
	core::pose::Pose pose;
	//pose.append_residue_by_bond( res_i );

	core::scoring::EnergyMap emap;
	sf1b->eval_ci_1b( res_i, pose, emap );
	sf1b->eval_cd_1b( res_i, pose, emap );
	sf1b->eval_ci_intrares_energy( res_i, pose, emap );
	sf1b->eval_cd_intrares_energy( res_i, pose, emap );
	core::Real intraE = emap.dot( sf1b->weights() );
	//score_grid.energy_ += intraE;
	score_grid.oneb_wtd_ += intraE;

	return score_grid;
}


// two-body energies
ReweightableRepEnergy
GridScorer::get_2b_energy(
	core::conformation::Residue const &res_i,
	core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo_i,
	core::conformation::Residue const &res_j,
	core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo_j,
	bool soft
) {
	core::scoring::ScoreFunctionOP sf = soft? sfxn_soft_ : sfxn_;

	ReweightableRepEnergy score_grid;

	core::Real weightatr = sf->get_weight( core::scoring::fa_atr );
	core::Real weightrep = sf->get_weight( core::scoring::fa_rep );
	core::Real weightsol = sf->get_weight( core::scoring::fa_sol );
	core::Real weightelec = sf->get_weight( core::scoring::fa_elec );
	core::Real weightLKb = sf->get_weight( core::scoring::lk_ball );
	core::Real weightLKbi = sf->get_weight( core::scoring::lk_ball_iso );
	core::Real weightLKbr = sf->get_weight( core::scoring::lk_ball_bridge );
	core::Real weightLKbru = sf->get_weight( core::scoring::lk_ball_bridge_uncpl );

	utility::vector1<hbDon> hbdon_neighborlist;
	utility::vector1<hbAcc> hbacc_neighborlist;
	numeric::xyzVector<core::Real> D,H,A,B,B0;
	core::scoring::hbonds::HBondOptions const & hbopt = sf->energy_method_options().hbond_options();

	//core::Real maxHbdis = 4.2; unused variable

	// etable/lk/elec
	for ( core::Size ii=1; ii<=res_i.natoms(); ++ii ) {
		bool ii_is_backbone =
			( ii<=res_i.last_backbone_atom() || (ii>res_i.nheavyatoms() && ii<res_i.first_sidechain_hydrogen()) );

		for ( core::Size jj=1; jj<=res_j.natoms(); ++jj ) {
			bool jj_is_backbone =
				( jj<=res_j.last_backbone_atom() || (jj>res_j.nheavyatoms() && jj<res_j.first_sidechain_hydrogen()) );

			// new version: only skip BB:BB interactions
			if ( ii_is_backbone && jj_is_backbone ) continue;

			core::Real faatr, farep, fasol1, fasol2;
			fast_eval_etable_split_fasol( etable_, res_i.atom(ii), res_j.atom(jj), faatr, farep, fasol1, fasol2 );
			core::Real faelec = coulomb_->eval_atom_atom_fa_elecE(res_i.xyz(ii), res_i.atomic_charge(ii), res_j.xyz(jj), res_j.atomic_charge(jj) );

			score_grid.fa_rep_ += weightrep*farep;
			//score_grid.energy_ += weightatr*faatr + weightsol*(fasol1+fasol2) + weightelec*faelec;
			score_grid.fa_atr_wtd_ += weightatr*faatr;
			score_grid.fa_sol_wtd_ += weightsol*(fasol1+fasol2);
			score_grid.fa_elec_wtd_ += weightelec*faelec;

			if ( res_i.atom_is_hydrogen( ii ) || res_j.atom_is_hydrogen( jj ) ) continue;

			// ii polar, jj occlusion
			if ( ii<=lkbrinfo_i->n_attached_waters().size() && lkbrinfo_i->n_attached_waters()[ii]>0 ) {
				core::Real fasol1_lkball = fasol1 * LKBe_->get_lk_fractional_contribution(
					res_j.xyz(jj), res_j.atom(jj).type(), lkbrinfo_i->n_attached_waters()[ii], lkbrinfo_i->waters()[ii] );
				//score_grid.energy_ += weightLKb*fasol1_lkball + weightLKbi*fasol1;
				score_grid.lk_ball_wtd_ += weightLKb*fasol1_lkball + weightLKbi*fasol1;
			}

			// jj polar, ii collusion
			if ( jj<=lkbrinfo_j->n_attached_waters().size() && lkbrinfo_j->n_attached_waters()[jj]>0 ) {
				core::Real fasol2_lkball = fasol2 * LKBe_->get_lk_fractional_contribution(
					res_i.xyz(ii), res_i.atom(ii).type(), lkbrinfo_j->n_attached_waters()[jj], lkbrinfo_j->waters()[jj] );
				//score_grid.energy_ += weightLKb*fasol2_lkball + weightLKbi*fasol2;
				score_grid.lk_ball_wtd_ += weightLKb*fasol2_lkball + weightLKbi*fasol2;
			}

			// ii+jj polar, bridge term
			if ( ii<=lkbrinfo_i->n_attached_waters().size() && lkbrinfo_i->n_attached_waters()[ii]>0 &&
					jj<=lkbrinfo_j->n_attached_waters().size() && lkbrinfo_j->n_attached_waters()[jj]>0 ) {
				core::Real fasol1_lkbridge = LKBe_->get_lkbr_fractional_contribution(
					res_i.xyz(ii), res_j.xyz(jj),
					lkbrinfo_i->n_attached_waters()[ii], lkbrinfo_j->n_attached_waters()[jj],
					lkbrinfo_i->waters()[ii], lkbrinfo_j->waters()[jj] );
				//score_grid.energy_ += weightLKbr * (fasol1+fasol2) * fasol1_lkbridge + weightLKbru * fasol1_lkbridge;
				score_grid.lk_ball_wtd_ += weightLKbr * (fasol1+fasol2) * fasol1_lkbridge + weightLKbru * fasol1_lkbridge;
			}
		}
	}

	// hbond
	// donor j->acceptor i
	for ( auto anum=res_i.accpt_pos().begin(),anume=res_i.accpt_pos().end(); anum!=anume; ++anum ) {
		// keep as is : bb hbonds are in the background grid still
		if ( !(*anum>res_i.last_backbone_atom() && *anum<=res_i.nheavyatoms()) ) continue;
		core::Size bnum = res_i.atom_base( *anum );
		core::Size b0num = res_i.abase2( *anum );
		A = res_i.xyz( *anum );
		B = res_i.xyz( bnum );
		B0 = res_i.xyz( b0num );

		core::scoring::hbonds::HBEvalTuple hbt;
		hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_i ));

		for ( auto hnum=res_j.Hpos_polar().begin(),hnume=res_j.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			if ( !(*hnum>res_j.first_sidechain_hydrogen()) ) continue;
			core::Size dnum = res_j.atom_base( *hnum );
			H = res_j.xyz( *hnum );
			D = res_j.xyz( dnum );
			hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_j ) );
			if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;
			core::Real hbe = get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
			//score_grid.energy_ += hbe;
			score_grid.hbond_wtd_ += hbe;
		}
	}

	// donor i->acceptor j
	for ( auto anum=res_j.accpt_pos().begin(),anume=res_j.accpt_pos().end(); anum!=anume; ++anum ) {
		// keep as is : bb hbonds are in the background grid still
		if ( !(*anum>res_j.last_backbone_atom() && *anum<=res_j.nheavyatoms()) ) continue;
		core::Size bnum = res_j.atom_base( *anum );
		core::Size b0num = res_j.abase2( *anum );
		A = res_j.xyz( *anum );
		B = res_j.xyz( bnum );
		B0 = res_j.xyz( b0num );

		core::scoring::hbonds::HBEvalTuple hbt;
		hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_j ));

		for ( auto hnum=res_i.Hpos_polar().begin(),hnume=res_i.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			if ( !(*hnum>res_i.first_sidechain_hydrogen()) ) continue;
			core::Size dnum = res_i.atom_base( *hnum );
			H = res_i.xyz( *hnum );
			D = res_i.xyz( dnum );
			hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_i ) );
			if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;
			core::Real hbe = get_hbond_score_weighted( sf, hb_database_, hbopt, hbt, D,H,A,B,B0 );
			//score_grid.energy_ += hbe;
			score_grid.hbond_wtd_ += hbe;
		}
	}

	return score_grid;
}


core::Real
GridScorer::clash_score( LigandConformer const &lig ) {
	lig.to_pose( ref_pose_ );

	// exact
	if ( exact_ ) {
		core::Real score_exact = (*sfxn_clash_)( *ref_pose_ );
		return score_exact;
	}

	core::conformation::Residue const &resLig = ref_pose_->residue( lig.ligand_id() );

	// 1 grid scores
	core::Real score_grid = 0.0;
	core::Real weightrep = sfxn_1b_clash_->get_weight( core::scoring::fa_rep );
	for ( core::Size i=1; i<=resLig.natoms(); ++i ) {
		int atmtype_i = resLig.atom(i).type();
		numeric::xyzVector< core::Real > idxX = (resLig.xyz(i) - origin_) / voxel_spacing_;

		core::Real penalty = move_to_boundary(idxX);

		core::Real rep = core::scoring::electron_density::interp_spline(coeffs_farep_[ atmtype_i ], idxX, true);

		rep = ixform_rep( rep ) + penalty;
		score_grid += weightrep * rep;
	}

	// 2 internal ligand energy
	core::scoring::EnergyMap emap;
	sfxn_1b_clash_->eval_ci_1b( ref_pose_->residue(lig.ligand_id()), *ref_pose_, emap );
	sfxn_1b_clash_->eval_cd_1b( ref_pose_->residue(lig.ligand_id()), *ref_pose_, emap );
	sfxn_1b_clash_->eval_ci_intrares_energy( ref_pose_->residue(lig.ligand_id()), *ref_pose_, emap );
	sfxn_1b_clash_->eval_cd_intrares_energy( ref_pose_->residue(lig.ligand_id()), *ref_pose_, emap );
	// sum in
	score_grid += emap.dot( sfxn_1b_clash_->weights() );

	return score_grid;
}

void
GridScorer::derivatives(
	core::pose::Pose &pose,
	LigandConformer const &lig,
	core::optimization::MinimizerMap &min_map
) {
	core::Real weightatr = sfxn_->get_weight( core::scoring::fa_atr );
	core::Real weightrep = sfxn_->get_weight( core::scoring::fa_rep );
	core::Real weightsol = sfxn_->get_weight( core::scoring::fa_sol );
	core::Real weightelec = sfxn_->get_weight( core::scoring::fa_elec );
	core::Real weightLKb = sfxn_->get_weight( core::scoring::lk_ball );
	core::Real weightLKbi = sfxn_->get_weight( core::scoring::lk_ball_iso );
	core::Real weightLKbr = sfxn_->get_weight( core::scoring::lk_ball_bridge );
	core::Real weightLKbru = sfxn_->get_weight( core::scoring::lk_ball_bridge_uncpl );
	bool useLKB = ( weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 );
	bool useLKBr = ( weightLKbr != 0 || weightLKbru != 0 );

	// 0 build ligand virtual sites
	core::Size nSCs = lig.moving_scs().size();
	utility::vector1< utility::vector1< core::Size > > allNwaters(nSCs+1);
	utility::vector1< utility::vector1< core::scoring::lkball::WaterCoords > > allwaters(nSCs+1);
	utility::vector1< core::scoring::lkball::LKB_ResidueInfoOP > alllkbrinfo(nSCs+1);
	if ( useLKB ) {
		for ( core::Size i=0; i<=nSCs; ++i ) {
			core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
			core::conformation::Residue const &res_i = pose.residue( resid );
			alllkbrinfo[i+1] = core::scoring::lkball::LKB_ResidueInfoOP( new core::scoring::lkball::LKB_ResidueInfo (res_i) );
			alllkbrinfo[i+1]->build_waters( res_i, true );
			allNwaters[i+1] = alllkbrinfo[i+1]->n_attached_waters();
			allwaters[i+1] = alllkbrinfo[i+1]->waters();
		}
	}

	// 1 grid scores
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
			min_map.atom_derivatives( resid )[j].f2() = 0;
		}
	}

	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );

		for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
			if ( i>0 && !( (j>res_i.last_backbone_atom() && j<=res_i.nheavyatoms()) || j>=res_i.first_sidechain_hydrogen() ) ) continue;

			int atmtype_ij = res_i.atom(j).type();
			numeric::xyzVector< core::Real > idxX = (res_i.xyz(j) - origin_) / voxel_spacing_;

			numeric::xyzVector< core::Real> dpenalty;
			core::Real penalty = move_to_boundary(idxX, dpenalty);
			//core::Real penalty = (i==0)? move_to_boundary(idxX, dpenalty) : 0.0 ; // no penalty for sidechain

			numeric::xyzVector< core::Real> datr(0,0,0),drep(0,0,0),dsol(0,0,0),delec(0,0,0),dlk(0,0,0);

			core::Real datrscale = dxform_atr( core::scoring::electron_density::interp_spline(coeffs_faatr_[ atmtype_ij ], idxX, true) );
			datr = (1/voxel_spacing_) * ( datrscale * core::scoring::electron_density::interp_dspline(coeffs_faatr_[ atmtype_ij ], idxX, true) );

			core::Real drepscale = w_rep_ * dxform_rep( core::scoring::electron_density::interp_spline(coeffs_farep_[ atmtype_ij ], idxX, true) );
			drep = (1/voxel_spacing_) * ( drepscale * core::scoring::electron_density::interp_dspline(coeffs_farep_[ atmtype_ij ], idxX, true) );

			dsol = (1/voxel_spacing_) * ( core::scoring::electron_density::interp_dspline(coeffs_fasol_[ atmtype_ij ], idxX, true) );
			delec = (1/voxel_spacing_) * res_i.atomic_charge(j) * ( core::scoring::electron_density::interp_dspline(coeffs_faelec_, idxX, true) );

			if ( useLKB ) {
				dlk = (1/voxel_spacing_) * (core::scoring::electron_density::interp_dspline(coeffs_lkball_[ atmtype_ij ], idxX, true));
			}

			if ( dpenalty[0] != 0 ) {
				datr[0] = drep[0] = dsol[0] = delec[0] = dlk[0] = 0.0;
			}
			if ( dpenalty[1] != 0 ) {
				datr[1] = drep[1] = dsol[1] = delec[1] = dlk[1] = 0.0;
			}
			if ( dpenalty[2] != 0 ) {
				datr[2] = drep[2] = dsol[2] = delec[2] = dlk[2] = 0.0;
			}

			min_map.atom_derivatives( resid )[j].f2() += dpenalty + weightatr*datr + weightrep*drep + weightsol*dsol + weightelec*delec + dlk;

			// VRT derivs
			if ( j > allNwaters[i+1].size() || allNwaters[i+1][j] == 0 ) continue;

			utility::vector1< numeric::xyzVector< core::Real > > nums_lk( allNwaters[i+1][j], numeric::xyzVector< core::Real >(0.,0.,0.) );
			utility::vector1< numeric::xyzVector< core::Real > > nums_lkbr( allNwaters[i+1][j], numeric::xyzVector< core::Real >(0.,0.,0.) );
			core::Real denom_lk=0.0, denom_lkbr=0.0;

			for ( core::Size k=1; k<=allNwaters[i+1][j]; ++k ) {
				numeric::xyzVector< core::Real > idxXj = (allwaters[i+1][j][k] - origin_) / voxel_spacing_;
				numeric::xyzVector<core::Real> dpenalty;
				move_to_boundary(idxXj,dpenalty);

				core::Real score_lk_ij( 0 );
				score_lk_ij = core::scoring::electron_density::interp_spline(coeffs_lkball_[ -atmtype_ij ], idxXj, true);
				score_lk_ij = std::exp( score_lk_ij / LK_fade_ );
				nums_lk[k] = score_lk_ij * (core::scoring::electron_density::interp_dspline(coeffs_lkball_[ -atmtype_ij ], idxXj, true));
				denom_lk += score_lk_ij;

				if ( dpenalty[0] != 0 ) nums_lk[k][0] = 0.0;
				if ( dpenalty[1] != 0 ) nums_lk[k][1] = 0.0;
				if ( dpenalty[2] != 0 ) nums_lk[k][2] = 0.0;

				if ( useLKBr ) {
					core::Real score_lkbr_ij( 0 );
					score_lkbr_ij = core::scoring::electron_density::interp_spline(coeffs_lkbridge_[ atmtype_ij ], idxXj, true);
					score_lkbr_ij = std::exp( -score_lkbr_ij / LK_fade_ );
					if ( penalty==0 ) {
						nums_lkbr[k] = score_lkbr_ij * (core::scoring::electron_density::interp_dspline(coeffs_lkbridge_[ atmtype_ij ], idxXj, true) );
					}
					denom_lkbr += score_lkbr_ij;
				}
			}

			core::scoring::lkball::WaterBuilders const & res_i_wb( alllkbrinfo[i+1]->get_water_builder( res_i , j ) );
			for ( core::Size k=1; k<=allNwaters[i+1][j]; ++k ) {
				numeric::xyzVector< core::Real > dall_dj(0.,0.,0.);

				if ( denom_lk!=0 ) dall_dj += nums_lk[k]/denom_lk;
				if ( denom_lkbr!=0 ) dall_dj += nums_lkbr[k]/denom_lkbr;
				if ( dall_dj.length_squared() < 1e-9 ) continue;

				core::Size atom1 = res_i_wb[k].atom1();
				numeric::xyzMatrix< core::Real >const & dwater_datom1 = alllkbrinfo[i+1]->atom1_derivs()[j][k];
				numeric::xyzVector< core::Real > dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
				numeric::xyzVector< core::Real > dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
				numeric::xyzVector< core::Real > dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
				numeric::xyzVector< core::Real > dRdatom1( dall_dj.dot( dwater_datom1x ), dall_dj.dot( dwater_datom1y ), dall_dj.dot( dwater_datom1z ) );
				min_map.atom_derivatives( resid )[atom1].f2() +=  (1/voxel_spacing_) * dRdatom1;

				core::Size atom2 = res_i_wb[k].atom2();
				numeric::xyzMatrix< core::Real >const & dwater_datom2 = alllkbrinfo[i+1]->atom2_derivs()[j][k];
				numeric::xyzVector< core::Real > dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
				numeric::xyzVector< core::Real > dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
				numeric::xyzVector< core::Real > dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
				numeric::xyzVector< core::Real > dRdatom2( dall_dj.dot( dwater_datom2x ), dall_dj.dot( dwater_datom2y ), dall_dj.dot( dwater_datom2z ) );
				min_map.atom_derivatives( resid )[atom2].f2() +=  (1/voxel_spacing_) * dRdatom2;

				core::Size atom3 = res_i_wb[k].atom3();
				numeric::xyzMatrix< core::Real >const & dwater_datom3 = alllkbrinfo[i+1]->atom3_derivs()[j][k];
				numeric::xyzVector< core::Real > dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
				numeric::xyzVector< core::Real > dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
				numeric::xyzVector< core::Real > dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
				numeric::xyzVector< core::Real > dRdatom3( dall_dj.dot( dwater_datom3x ), dall_dj.dot( dwater_datom3y ), dall_dj.dot( dwater_datom3z ) );
				min_map.atom_derivatives( resid )[atom3].f2() +=  (1/voxel_spacing_) * dRdatom3;
			}
		}
	}

	// 2 hbond derivs
	utility::vector1<hbDon> hbdon_neighborlist;
	utility::vector1<hbAcc> hbacc_neighborlist;
	numeric::xyzVector<core::Real> D,H,A,B,B0;

	core::scoring::methods::EnergyMethodOptions const &e_opts = sfxn_->energy_method_options();
	core::scoring::hbonds::HBondOptions const & hbopt = e_opts.hbond_options();
	core::scoring::hbonds::HBondDatabaseCOP database = core::scoring::hbonds::HBondDatabase::get_database( hbopt.params_database_tag() );

	core::Real maxHbdis = 4.2;

	core::scoring::hbonds::HBondDerivs dhb_dxs;

	core::Real sfwt=0;
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		for ( auto anum=res_i.accpt_pos().begin(),anume=res_i.accpt_pos().end(); anum!=anume; ++anum ) {
			if ( i>0 && !(*anum>res_i.last_backbone_atom() && *anum<=res_i.nheavyatoms()) ) continue;

			hbdonors_.get_neighbors( res_i.xyz( *anum ), hbdon_neighborlist, maxHbdis );
			core::Size bnum = res_i.atom_base( *anum );
			core::Size b0num = res_i.abase2( *anum );
			A = res_i.xyz( *anum );
			B = res_i.xyz( bnum );
			B0 = res_i.xyz( b0num );

			core::scoring::hbonds::HBEvalTuple hbt;
			hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_i ));

			for ( core::Size i=1; i<=hbdon_neighborlist.size(); ++i ) {
				D = hbdon_neighborlist[i].D;
				H = hbdon_neighborlist[i].H;
				hbt.don_type( hbdon_neighborlist[i].dontype );

				if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;

				core::Real energy=0;
				core::scoring::hbonds::hb_energy_deriv(
					*database, hbopt, hbt,
					D,H,A,B,B0,
					energy, true, dhb_dxs
				);
				if ( energy<core::scoring::hbonds::MAX_HB_ENERGY ) {
					switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
					case core::scoring::hbonds::hbw_SR_BB_SC :
					case core::scoring::hbonds::hbw_LR_BB_SC :
						sfwt = sfxn_->get_weight( core::scoring::hbond_bb_sc );
						break;
					case core::scoring::hbonds::hbw_SC :
					default :
						sfwt = sfxn_->get_weight( core::scoring::hbond_sc );
					}
					min_map.atom_derivatives( resid )[*anum].f2() += sfwt * dhb_dxs.acc_deriv.f2();
					min_map.atom_derivatives( resid )[bnum].f2()  += sfwt * dhb_dxs.abase_deriv.f2();
					min_map.atom_derivatives( resid )[b0num].f2() += sfwt * dhb_dxs.abase2_deriv.f2();
				}
			}
		}
		for ( auto hnum=res_i.Hpos_polar().begin(),hnume=res_i.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			if ( i>0 && !(*hnum>res_i.first_sidechain_hydrogen()) ) continue;

			core::Size dnum = res_i.atom_base( *hnum );
			hbacceptors_.get_neighbors( res_i.xyz( dnum ), hbacc_neighborlist, maxHbdis );
			H = res_i.xyz( *hnum );
			D = res_i.xyz( dnum );

			core::scoring::hbonds::HBEvalTuple hbt;
			hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_i ));

			for ( core::Size i=1; i<=hbacc_neighborlist.size(); ++i ) {
				A = hbacc_neighborlist[i].A;
				B = hbacc_neighborlist[i].B;
				B0 = hbacc_neighborlist[i].B0;
				hbt.acc_type( hbacc_neighborlist[i].acctype );

				if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;

				core::Real energy=0;
				core::scoring::hbonds::hb_energy_deriv(
					*database, hbopt, hbt,
					D,H,A,B,B0,
					energy, true, dhb_dxs
				);
				if ( energy<core::scoring::hbonds::MAX_HB_ENERGY ) {
					switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
					case core::scoring::hbonds::hbw_SR_BB_SC :
					case core::scoring::hbonds::hbw_LR_BB_SC :
						sfwt = sfxn_->get_weight( core::scoring::hbond_bb_sc );
						break;
					case core::scoring::hbonds::hbw_SC :
					default :
						sfwt = sfxn_->get_weight( core::scoring::hbond_sc );
					}
					min_map.atom_derivatives( resid )[*hnum].f2() += sfwt * dhb_dxs.h_deriv.f2();
					min_map.atom_derivatives( resid )[dnum].f2()  += sfwt * dhb_dxs.don_deriv.f2();
				}
			}
		}
	}

	// 3 1b derivs
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		core::scoring::MinimizationGraph g( 1 );
		core::scoring::EnergyMap fixed_energies;
		core::scoring::MinimizationNode & minnode = *( g.get_minimization_node( 1 ) );
		sfxn_1b_->setup_for_minimizing_for_node(
			minnode, pose.residue(resid), pose.residue_data( resid ),
			min_map, pose, false, fixed_energies ); //?
		core::scoring::eval_atom_derivatives_for_minnode(
			minnode, res_i, pose, sfxn_1b_->weights(), min_map.atom_derivatives( resid ) );
	}

	// 4 : 2b derivs for ligand : sidechains
	//      brute force calculations for interaction graph ... should be reasonable as long as too many sidechains are not allowed to move...
	if ( nSCs > 0 ) {
		for ( core::Size i=0; i<=nSCs; ++i ) {
			core::Size resid_i = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
			core::conformation::Residue const &res_i = pose.residue( resid_i );

			for ( core::Size j=i+1; j<=nSCs; ++j ) {
				core::Size resid_j = lig.moving_scs()[j];
				core::conformation::Residue const &res_j = pose.residue( resid_j );

				// TO DO: interaction check here

				// 4a etable
				for ( core::Size ii=1; ii<=res_i.natoms(); ++ii ) {
					bool ii_is_backbone = (i!=0) && (
						ii<=res_i.last_backbone_atom() || (ii>res_i.nheavyatoms() && ii<res_i.first_sidechain_hydrogen()) );

					for ( core::Size jj=1; jj<=res_j.natoms(); ++jj ) {
						bool jj_is_backbone = (j!=0) && (
							jj<=res_j.last_backbone_atom() || (jj>res_j.nheavyatoms() && jj<res_j.first_sidechain_hydrogen()) );

						// new version: only skip BB:BB interactions
						if ( ii_is_backbone && jj_is_backbone ) continue;

						core::Real d_faatr, d_farep, d_fasol, d_faelec, invD;
						numeric::xyzVector<core::Real> x_ij = (res_i.xyz(ii)-res_j.xyz(jj));
						etable_.analytic_etable_derivatives(res_i.atom(ii), res_j.atom(jj),  d_faatr, d_farep, d_fasol, invD );
						d_farep *= w_rep_;

						core::Real dis2 = x_ij.length_squared();
						d_faelec = coulomb_->eval_dfa_elecE_dr_over_r( dis2, res_i.atomic_charge(ii), res_j.atomic_charge(jj) );

						core::Real dE_dR_over_r = ( weightatr*d_faatr + weightrep*d_farep + weightsol*d_fasol) * invD + weightelec*d_faelec; //?

						min_map.atom_derivatives( resid_i )[ii].f2() += dE_dR_over_r*x_ij;
						min_map.atom_derivatives( resid_j )[jj].f2() -= dE_dR_over_r*x_ij;

						if ( res_i.atom_is_hydrogen( ii ) || res_j.atom_is_hydrogen( jj ) ) continue;

						bool ii_has_waters = (ii<=allNwaters[i+1].size() && allNwaters[i+1][ii]>0);
						bool jj_has_waters = (jj<=allNwaters[j+1].size() && allNwaters[j+1][jj]>0);

						if ( ii_has_waters || jj_has_waters ) {
							LKBe_->sum_deriv_contributions_for_heavyatom_pair(
								dis2,
								ii,res_i,*(alllkbrinfo[i+1]),
								jj,res_j,*(alllkbrinfo[j+1]),
								pose, sfxn_->weights(),
								1.0,
								min_map.atom_derivatives( resid_i ),
								min_map.atom_derivatives( resid_j )
							);
						}
					}
				}
			}
		}

		for ( core::Size i=0; i<=nSCs; ++i ) {
			core::Size resid_i = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
			core::conformation::Residue const &res_i = pose.residue( resid_i );

			for ( core::Size j=i+1; j<=nSCs; ++j ) {
				core::Size resid_j = lig.moving_scs()[j];
				core::conformation::Residue const &res_j = pose.residue( resid_j );

				// 4b hbond
				// donor j->acceptor i
				for ( auto anum=res_i.accpt_pos().begin(),anume=res_i.accpt_pos().end(); anum!=anume; ++anum ) {
					if ( i>0 && !(*anum>res_i.last_backbone_atom() && *anum<=res_i.nheavyatoms()) ) continue;
					core::Size bnum = res_i.atom_base( *anum );
					core::Size b0num = res_i.abase2( *anum );
					A = res_i.xyz( *anum );
					B = res_i.xyz( bnum );
					B0 = res_i.xyz( b0num );

					core::scoring::hbonds::HBEvalTuple hbt;
					hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_i ));

					for ( auto hnum=res_j.Hpos_polar().begin(),hnume=res_j.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
						if ( j>0 && !(*hnum>res_j.first_sidechain_hydrogen()) ) continue;
						core::Size dnum = res_j.atom_base( *hnum );
						H = res_j.xyz( *hnum );
						D = res_j.xyz( dnum );
						hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_j ) );
						if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;

						core::Real energy=0;
						core::scoring::hbonds::hb_energy_deriv(
							*database, hbopt, hbt,
							D,H,A,B,B0,
							energy, true, dhb_dxs
						);
						if ( energy<core::scoring::hbonds::MAX_HB_ENERGY ) {
							switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
							case core::scoring::hbonds::hbw_SR_BB_SC :
							case core::scoring::hbonds::hbw_LR_BB_SC :
								sfwt = sfxn_->get_weight( core::scoring::hbond_bb_sc );
								break;
							case core::scoring::hbonds::hbw_SC :
							default :
								sfwt = sfxn_->get_weight( core::scoring::hbond_sc );
							}
							min_map.atom_derivatives( resid_i )[*anum].f2() += sfwt * dhb_dxs.acc_deriv.f2();
							min_map.atom_derivatives( resid_i )[bnum].f2()  += sfwt * dhb_dxs.abase_deriv.f2();
							min_map.atom_derivatives( resid_i )[b0num].f2() += sfwt * dhb_dxs.abase2_deriv.f2();
							min_map.atom_derivatives( resid_j )[*hnum].f2() += sfwt * dhb_dxs.h_deriv.f2();
							min_map.atom_derivatives( resid_j )[dnum].f2()  += sfwt * dhb_dxs.don_deriv.f2();
						}
					}
				}

				// donor i->acceptor j
				for ( auto anum=res_j.accpt_pos().begin(),anume=res_j.accpt_pos().end(); anum!=anume; ++anum ) {
					if ( i>0 && !(*anum>res_j.last_backbone_atom() && *anum<=res_j.nheavyatoms()) ) continue;
					core::Size bnum = res_j.atom_base( *anum );
					core::Size b0num = res_j.abase2( *anum );
					A = res_j.xyz( *anum );
					B = res_j.xyz( bnum );
					B0 = res_j.xyz( b0num );

					core::scoring::hbonds::HBEvalTuple hbt;
					hbt.acc_type( core::scoring::hbonds::get_hb_acc_chem_type( *anum, res_j ));

					for ( auto hnum=res_i.Hpos_polar().begin(),hnume=res_i.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
						if ( i>0 && !(*hnum>res_i.first_sidechain_hydrogen()) ) continue;
						core::Size dnum = res_i.atom_base( *hnum );
						H = res_i.xyz( *hnum );
						D = res_i.xyz( dnum );
						hbt.don_type( core::scoring::hbonds::get_hb_don_chem_type( dnum, res_i ) );
						if ( (A-H).length_squared() > core::scoring::hbonds::MAX_R2 ) continue;

						core::Real energy=0;
						core::scoring::hbonds::hb_energy_deriv(
							*database, hbopt, hbt,
							D,H,A,B,B0,
							energy, true, dhb_dxs
						);
						if ( energy<core::scoring::hbonds::MAX_HB_ENERGY ) {
							switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
							case core::scoring::hbonds::hbw_SR_BB_SC :
							case core::scoring::hbonds::hbw_LR_BB_SC :
								sfwt = sfxn_->get_weight( core::scoring::hbond_bb_sc );
								break;
							case core::scoring::hbonds::hbw_SC :
							default :
								sfwt = sfxn_->get_weight( core::scoring::hbond_sc );
							}
							min_map.atom_derivatives( resid_j )[*anum].f2() += sfwt * dhb_dxs.acc_deriv.f2();
							min_map.atom_derivatives( resid_j )[bnum].f2()  += sfwt * dhb_dxs.abase_deriv.f2();
							min_map.atom_derivatives( resid_j )[b0num].f2() += sfwt * dhb_dxs.abase2_deriv.f2();
							min_map.atom_derivatives( resid_i )[*hnum].f2() += sfwt * dhb_dxs.h_deriv.f2();
							min_map.atom_derivatives( resid_i )[dnum].f2()  += sfwt * dhb_dxs.don_deriv.f2();
						}
					}
				}
			}
		}
	}

	// now convert all to f1/f2s
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );

		for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
			if ( i>0 && !( (j>res_i.last_backbone_atom() && j<=res_i.nheavyatoms()) || j>=res_i.first_sidechain_hydrogen() ) ) continue;
			numeric::xyzVector< core::Real > deriv_ij = min_map.atom_derivatives( resid )[j].f2();
			min_map.atom_derivatives( resid )[j].f1() = res_i.xyz(j).cross( res_i.xyz(j)-deriv_ij );
		}
	}
}

core::Real
GridScorer::dof_derivative(
	core::pose::Pose &pose,
	core::optimization::MinimizerMap &min_map,
	core::id::DOF_ID const & dof_id,
	core::id::TorsionID const & torsion_id
) {
	core::Size const resid = torsion_id.valid() ? torsion_id.rsd() : dof_id.atom_id().rsd();
	core::conformation::Residue const & rsd = pose.residue( resid );

	core::scoring::EnergyMap fixed_energies;
	core::scoring::MinimizationGraph g( 1 );
	core::scoring::MinimizationNode & minnode = *( g.get_minimization_node( 1 ) );
	sfxn_1b_->setup_for_minimizing_for_node(
		minnode, pose.residue(resid), pose.residue_data( resid ),
		min_map, pose, false, fixed_energies );
	core::scoring::eval_atom_derivatives_for_minnode(
		minnode, rsd, pose, sfxn_1b_->weights(), min_map.atom_derivatives( resid ) );
	return core::scoring::eval_dof_deriv_for_minnode( minnode, rsd, pose, dof_id, torsion_id, *sfxn_1b_, sfxn_1b_->weights() );
}

void
GridScorer::debug_deriv(
	core::pose::Pose &pose,
	LigandConformer const &lig,
	core::optimization::MinimizerMap &min_map
) {
	// first get analytic deriv
	derivatives( pose, lig, min_map );

	// DEBUG DERIVS
	core::Size nSCs = lig.moving_scs().size();

	core::pose::Pose posecopy = pose;
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? lig.ligand_id() : lig.moving_scs()[i]);
		core::conformation::Residue const &res_i = pose.residue( resid );
		for ( core::Size j=1; j<=res_i.natoms(); ++j ) {
			bool j_is_backbone = (i!=0) && (
				j<=res_i.last_backbone_atom() || (j>res_i.nheavyatoms() && j<res_i.first_sidechain_hydrogen()) );

			if ( j_is_backbone ) continue;

			core::Real scpx,scmx,scpy,scmy,scpz,scmz;
			numeric::xyzVector<core::Real> x_ij = posecopy.xyz( core::id::AtomID(j,resid) );

			x_ij[0] += 0.0001;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scpx = score(posecopy,lig);
			x_ij[0] -= 0.0002;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scmx = score(posecopy,lig);
			x_ij[0] += 0.0001;

			x_ij[1] += 0.0001;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scpy = score(posecopy,lig);
			x_ij[1] -= 0.0002;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scmy = score(posecopy,lig);
			x_ij[1] += 0.0001;

			x_ij[2] += 0.0001;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scpz = score(posecopy,lig);
			x_ij[2] -= 0.0002;
			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );
			scmz = score(posecopy,lig);
			x_ij[2] += 0.0001;

			posecopy.set_xyz( core::id::AtomID(j,resid), x_ij );

			core::Vector f2( min_map.atom_derivatives( resid )[j].f2() );
			core::Real nfx( (scpx-scmx)/0.0002 ), nfy( (scpy-scmy)/0.0002 ), nfz((scpz-scmz)/0.0002);
			printf("DERIV %3d %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
				int(resid),int(j), nfx,nfy,nfz,f2[0],f2[1],f2[2]);

			//numeric::xyzVector< core::Real > deriv_ij = numeric::xyzVector<core::Real>(nfx,nfy,nfz);
			//min_map.atom_derivatives( resid )[j].f1() = res_i.xyz(j).cross( res_i.xyz(j)-deriv_ij );
			//min_map.atom_derivatives( resid )[j].f2() = deriv_ij;
		}
	}
}

// run an optimization cycle (ramped repack/min) for a single ligand conformer
core::Real
GridScorer::optimize(
	LigandConformer &lig,
	utility::vector1<core::Real> ramp_schedule, // YES I KNOW THIS IS PASSED BY VALUE.  IT IS INTENTIONAL
	utility::vector1< PlaceableRotamers > &placeable_rotdb,  // by non-const ref since we update "in place"
	RotamerPairEnergies &rot_energies                        // by non-const ref since we update "in place"
) {
	// maybe make the configurable?
	core::optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 1e-3, true, false, false );
	minopt.max_iter( maxiter_minimize_ );

	core::Real finalscore = 0.0;
	if ( exact_ ) {
		//fpd relax will idealize sidechains, which are no longer accurately represented by the gene!
		TR.Error << "Exact optimization unsupported!" << std::endl;
		utility_exit();
	} else {
		// ensure reasonable ramp schedule (move to parse time?)
		if ( ramp_schedule.size() == 0 ) {
			ramp_schedule.push_back(1.0);
		}
		if ( ramp_schedule[ ramp_schedule.size() ] != 1.0 ) {
			TR.Warning << "Warning! ramp_schedule does not end with full weight (1.0)." << std::endl;
		}

		// temporarily turn off debug info
		bool debug_orig = debug_;
		debug_ = false;

		// ramping is w.r.t. original rep weight
		core::Real wbase = get_w_rep();

		for ( core::Size iramp = 1; iramp <= ramp_schedule.size(); ++iramp ) {
			//set_w_rep( wbase*ramp_schedule[iramp] );
			// FD: commented out => always run for maxiter_minimize_
			//if (iramp < ramp_schedule.size()) {
			// minopt.max_iter( (core::Size) (maxiter_minimize_*0.25) ); // ? not sure about the 0.25 here
			//} else {
			// minopt.max_iter( maxiter_minimize_ );
			//}

			auto timer1 = std::chrono::system_clock::now();
			set_w_rep( wbase*ramp_schedule[iramp] );
			packer_loop( lig, placeable_rotdb, rot_energies );

			auto timer2 = std::chrono::system_clock::now();
			set_w_rep( wbase );
			minimizer_loop( lig, minopt );

			auto timer3 = std::chrono::system_clock::now();

			pack_time_ += (timer2 - timer1);
			min_time_ += (timer3 - timer2);
		} // everything stored dofs

		// restore original rep weight / debug state
		//set_w_rep(wbase);
		debug_ = debug_orig;

		finalscore = score(lig);
	}

	return finalscore;
}

//
core::Real
GridScorer::minimizer_loop(
	LigandConformer &lig,
	core::optimization::MinimizerOptions const &minopt
) {
	core::Real finalscore=0.0;
	core::Size nSCs = lig.moving_scs().size();

	core::pose::PoseOP minipose( new core::pose::Pose );
	LigandConformer minilig;
	lig.to_minipose( minipose, minilig );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	core::Size lig_jump = minipose->fold_tree().get_jump_that_builds_residue( minilig.ligand_id() );
	mm->set_jump( lig_jump, true );
	for ( core::Size i=0; i<=nSCs; ++i ) {
		core::Size resid = (i==0 ? minilig.ligand_id() : minilig.moving_scs()[i]);
		mm->set_chi( resid, true );
	}

	// fd : minimize on the grid
	//      use standard minimization machinery but not scoring machinery
	core::optimization::MinimizerMap min_map;
	min_map.setup( *minipose, *mm );

	core::optimization::Multivec dofs( min_map.nangles() );
	min_map.copy_dofs_from_pose( *minipose, dofs );

	GriddedAtomTreeMultifunc f_i( minilig, *minipose, *this, min_map );
	core::optimization::Minimizer min_i( f_i, minopt );

	min_i.run( dofs );

	min_map.reset_jump_rb_deltas( *minipose, dofs );
	min_map.copy_dofs_to_pose( *minipose, dofs );

	lig.update_conf_from_minipose( minipose);
	finalscore = score(lig);

	return finalscore;
}

//
core::Real
GridScorer::packer_loop(
	LigandConformer &lig,
	utility::vector1< PlaceableRotamers > &placeable_rotdb, // by non-const ref since we update "in place"
	RotamerPairEnergies &rotamer_energies                   // by non-const ref since we update "in place"
) {
	lig.to_pose( ref_pose_ );
	core::Real start_score = score( *ref_pose_, lig );

	core::Size npos = lig.moving_scs().size();
	if ( npos == 0 || placeable_rotdb.size() == 0 ) return start_score;

	// [0] precompute current rotamers with all other rotamers
	for ( core::Size ii=1; ii<=npos; ++ii ) {
		core::Size nrot_i = placeable_rotdb[ii].size();
		core::Size resid_i = lig.moving_scs()[ii];
		placeable_rotdb[ii][nrot_i].restype = lig.get_protein_restype(ii);
		placeable_rotdb[ii][nrot_i].chis = lig.get_protein_chis(ii);
		placeable_rotdb[ii][nrot_i].lkbrinfo =
			core::scoring::lkball::LKB_ResidueInfoOP(new core::scoring::lkball::LKB_ResidueInfo (ref_pose_->residue( resid_i )) );
		rotamer_energies.energy1b( ii, nrot_i ) = this->get_1b_energy(
			ref_pose_->residue( resid_i ), placeable_rotdb[ii][nrot_i].lkbrinfo );
	}

	for ( core::Size ii=1; ii<=npos; ++ii ) {
		core::Size nrot_i = placeable_rotdb[ii].size();
		core::Size resid_i = lig.moving_scs()[ii];

		// fd this could be more efficient...
		if ( ref_pose_->residue_type(resid_i).name() != placeable_rotdb[ii][nrot_i].restype->name() ) {  // better way?
			core::conformation::Residue newres( placeable_rotdb[ii][nrot_i].restype, false);
			ref_pose_->replace_residue( resid_i, newres, true );
		}
		for ( core::Size ichi = 1; ichi <= ref_pose_->residue_type( resid_i ).nchi(); ++ichi ) {
			ref_pose_->set_chi( ichi, resid_i, placeable_rotdb[ii][nrot_i].chis[ichi] );
		}

		for ( core::Size jj=1; jj<=npos; ++jj ) {
			if ( ii==jj ) continue;
			core::Size nrot_j = placeable_rotdb[jj].size();
			core::Size resid_j = lig.moving_scs()[jj];
			for ( core::Size jrot=1; jrot<=nrot_j; ++jrot ) {

				if ( ref_pose_->residue_type(resid_j).name() != placeable_rotdb[jj][jrot].restype->name() ) {
					core::conformation::Residue newres( placeable_rotdb[jj][jrot].restype, false);
					ref_pose_->replace_residue( resid_j, newres, true );
				}
				for ( core::Size jchi = 1; jchi <= ref_pose_->residue_type( resid_j ).nchi(); ++jchi ) {
					ref_pose_->set_chi( jchi, resid_j, placeable_rotdb[jj][jrot].chis[jchi] );
				}

				rotamer_energies.energy2b( ii, nrot_i, jj, jrot ) = this->get_2b_energy(
					ref_pose_->residue( resid_i ), placeable_rotdb[ii][nrot_i].lkbrinfo,
					ref_pose_->residue( resid_j ), placeable_rotdb[jj][jrot].lkbrinfo
				);
			}
		}
	}

	// [1] precompute ligand with all other rotamers
	core::scoring::lkball::LKB_ResidueInfoOP lkbr_lig (
		new core::scoring::lkball::LKB_ResidueInfo(ref_pose_->residue(lig.ligand_id())) );
	for ( core::Size ii=1; ii<=npos; ++ii ) {
		core::Size nrot_i = placeable_rotdb[ii].size();
		core::Size resid_i = lig.moving_scs()[ii];
		for ( core::Size irot=1; irot<=nrot_i; ++irot ) { // include placeholder rotamer
			// fd this could be more efficient...
			if ( ref_pose_->residue_type(resid_i).name() != placeable_rotdb[ii][irot].restype->name() ) {
				core::conformation::Residue newres(*placeable_rotdb[ii][irot].restype, false);
				ref_pose_->replace_residue( resid_i, newres, true );
			}
			for ( core::Size ichi = 1; ichi <= ref_pose_->residue_type( resid_i ).nchi(); ++ichi ) {
				ref_pose_->set_chi( ichi, resid_i, placeable_rotdb[ii][irot].chis[ichi] );
			}
			rotamer_energies.energyBG( ii, irot ) = this->get_2b_energy(
				ref_pose_->residue( resid_i ), placeable_rotdb[ii][irot].lkbrinfo,
				ref_pose_->residue( lig.ligand_id() ), lkbr_lig
			);
		}
	}

	// [2] setup
	// # cycles and temperature schedule
	core::Size nmacro = 25, nstep = packer_cycles_*npos;
	core::Real temp = 10, temp0 = 1;
	core::Real tmult = std::pow(temp0/temp, 1.0/(nmacro-1));

	// initial assignment & MC parameters
	utility::vector1< core::Real > rotamer_assignment( npos ), bestRot_assignment;
	for ( core::Size ii=1; ii<=npos; ++ii ) {
		core::Size nrot_i = placeable_rotdb[ii].size();
		rotamer_assignment[ii] = nrot_i;
	}
	bestRot_assignment = rotamer_assignment;

	core::Real curr_score=0.0, best_score, test_score;
	for ( core::Size ii=1; ii<=npos; ++ii ) {
		curr_score += rotamer_energies.energyBG( ii, rotamer_assignment[ii] ).score( w_rep_ );
		curr_score += rotamer_energies.energy1b( ii, rotamer_assignment[ii] ).score( w_rep_ );
		for ( core::Size jj=ii+1; jj<=npos; ++jj ) {
			curr_score += rotamer_energies.energy2b( ii, rotamer_assignment[ii], jj, rotamer_assignment[jj] ).score( w_rep_ );
		}
	}
	best_score = curr_score;

	// [3] main MC loop
	for ( core::Size imacro=1; imacro<=nmacro; ++imacro ) {
		for ( core::Size istep=1; istep<=nstep; ++istep ) {
			core::Size i_res = numeric::random::rg().random_range( 1, npos );
			core::Size i_nrot( placeable_rotdb[i_res].size() );

			core::Size i_rot(1), old_i_rot = rotamer_assignment[i_res];

			// weighted
			//core::Real irot_prob( numeric::random::rg().uniform() );
			//while (i_rot < i_nrot && placeable_rotdb[i_res][i_rot].prob_accum < irot_prob) { ++i_rot; }
			// uniform
			i_rot = numeric::random::rg().random_range( 1, i_nrot );

			// possible speedup by only computing delta energies?
			rotamer_assignment[i_res] = i_rot;
			test_score = 0.0;
			for ( core::Size ii=1; ii<=npos; ++ii ) {
				test_score += rotamer_energies.energyBG( ii, rotamer_assignment[ii] ).score( w_rep_ );
				test_score += rotamer_energies.energy1b( ii, rotamer_assignment[ii] ).score( w_rep_ );
				for ( core::Size jj=ii+1; jj<=npos; ++jj ) {
					test_score += rotamer_energies.energy2b( ii, rotamer_assignment[ii], jj, rotamer_assignment[jj] ).score( w_rep_ );
				}
			}

			core::Real scoreDel = std::max( -20.0, std::min( 20.0, (curr_score-test_score)/temp) );
			core::Real paccept = std::exp( scoreDel );
			if ( paccept > 1 || numeric::random::rg().uniform() < paccept ) {
				curr_score = test_score;
			} else {
				rotamer_assignment[i_res] = old_i_rot; // revert
			}

			if ( curr_score < best_score ) {
				best_score = curr_score;
				bestRot_assignment = rotamer_assignment;
			}
		}

		// recover low
		curr_score = best_score;
		rotamer_assignment = bestRot_assignment;

		temp *= tmult;
	}

	// [4] postprocessing
	for ( core::Size ii=1; ii<=npos; ++ii ) {
		lig.set_protein_restype( ii, placeable_rotdb[ii][ rotamer_assignment[ii] ].restype );
		lig.set_protein_chis( ii, placeable_rotdb[ii][ rotamer_assignment[ii] ].chis );
	}
	best_score = score(lig);

	//TR << "packer_loop: [w_rep=" << w_rep_ << "] " << curr_score << "/" << best_score << std::endl;

	return best_score;
}


} // ga_dock
} // ligand_docking
} // protocols
