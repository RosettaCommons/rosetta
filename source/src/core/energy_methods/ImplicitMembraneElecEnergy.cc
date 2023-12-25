// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/ImplicitMembraneElecEnergy.cc
/// @brief Energy method for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)
///@author rsamant2(rituparna@utexas.edu) after Jun2021

//This was made following FA_ElecEnergyAroAro.cc by ralford.
//RS has added some functions from FA_elecEnergy.cc

// Unit headers
#include <core/energy_methods/ImplicitMembraneElecEnergy.hh>
#include <core/energy_methods/ImplicitMembraneElecEnergyCreator.hh>

#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>

// Project Headers
#include <core/energy_methods/ImplicitMembraneCoulomb.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/Atom.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/NeighborList.tmpl.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MinimizerMapBase.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/vector1.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>

static basic::Tracer TR( "core.energy_methods.ImplicitMembraneElecEnergy" );

namespace core {
namespace energy_methods {
// namespace elec {

using namespace core::scoring;
using namespace core::scoring::elec;

/// @details This must return a fresh instance of the ImplicitMembraneElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ImplicitMembraneElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< ImplicitMembraneElecEnergy >( options );
}

ScoreTypes
ImplicitMembraneElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_imm_elec );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
ImplicitMembraneElecEnergy::ImplicitMembraneElecEnergy(
	methods::EnergyMethodOptions const & options
):
	parent( options ),
	membrane_coulomb_( utility::pointer::make_shared< core::energy_methods::ImplicitMembraneCoulomb >() )
{
	// this is a EnergyMethod Function. membrane_coulomb_ calls membrane_coulomb_initialize()
	set_score_types( utility::pointer::make_shared< ImplicitMembraneElecEnergyCreator >() );
}

////////////////////////////////////////////////////////////////////////////
ImplicitMembraneElecEnergy::ImplicitMembraneElecEnergy( ImplicitMembraneElecEnergy const & src ):
	parent( src ),
	membrane_coulomb_( utility::pointer::make_shared< core::energy_methods::ImplicitMembraneCoulomb >() )
{
	set_score_types( utility::pointer::make_shared< ImplicitMembraneElecEnergyCreator >() );
}

// destructor (important for properly forward-declaring smart-pointer members)
ImplicitMembraneElecEnergy::~ImplicitMembraneElecEnergy() {}

////////////////////////////////////////////////////////////////////////////
methods::EnergyMethodOP
ImplicitMembraneElecEnergy::clone() const
{
	//return methods::EnergyMethodOP(utility::pointer::make_shared< ImplicitMembraneElecEnergy >( *this ));
	return utility::pointer::make_shared< ImplicitMembraneElecEnergy >( *this );
}

////////////////////////////////////////////////////////////////////////////
void
ImplicitMembraneElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

////////////////////////////////////////////////////////////////////////////
void
ImplicitMembraneElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn) const
{
	TR.Debug << "setup_for_scoring" << std::endl;
	TR.Debug << "====use_nblist()====" << pose.energies().use_nblist() << std::endl;
	TR.Debug << "teh weight of this term: " << scfxn.weights()[ core::scoring::fa_imm_elec ];
	pose.update_residue_neighbors();
	//**this energy is calculated only when the nb_list is changed**

}


/////////////////////////////////////////////////////////////////////////////
// scoring
////////////////////////////////////////////////////////////////////////////
void
ImplicitMembraneElecEnergy::residue_pair_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	/*TR << "residue_pair_energy" << std::endl;*/
	if ( pose.energies().use_nblist() ) return;
	using namespace etable::count_pair;

	Real score(0.0);

	core::conformation::Conformation const & conf( pose.conformation() );
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );

	// temporary
	bool count_pair_full = false;
	bool count_pair_hybrid = false;

	Real attached_h_max_dis2 = hydrogen_interaction_cutoff2();

	if ( rsd1.name3() == "MEM" || !rsd1.is_protein() || rsd2.name3() == "MEM" || !rsd2.is_protein() ) return;
	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;
	//defines_score_for_residue_pair is an function inheritted from FaElecEnergy.cc.

	// NULL if no info
	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topologies
		// also assuming crossover of 4, should be closest (?) to classic rosetta

		bool const is_rsd1_ncaa_polymer = rsd1.is_polymer() && !rsd1.type().is_base_type();
		bool const is_rsd2_ncaa_polymer = rsd2.is_polymer() && !rsd2.type().is_base_type();
		bool const is_cp3full = count_pair_full ||
			(count_pair_hybrid && ( rsd1.is_ligand() || rsd2.is_ligand() || is_rsd1_ncaa_polymer || is_rsd2_ncaa_polymer )); // if any of rsd1/rsd2 is ncaa_polymer

		CountPairFunctionOP cpfxn = is_cp3full ?
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3FULL ) :
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		Real d2;
		for ( Size ii = 1; ii <= rsd1.nheavyatoms(); ++ii ) {
			Size ii_rep = get_countpair_representative_atom( rsd1.type(), ii );
			for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
				Size jj_rep = get_countpair_representative_atom( rsd2.type(), jj );

				// compute env hydration for each position
				// Real i_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(ii) ) );
				// Real j_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(jj) ) );
				Real i_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), ii );
				Real j_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), jj );

				Real weight=1.0;
				Size path_dist( 0 );
				if ( cpfxn->count( ii_rep, jj_rep, weight, path_dist ) ) {
					score += score_atom_pair( rsd1, rsd2, ii, jj, i_hyd, j_hyd, emap, weight, d2 );
				} else {
					d2 = rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) );
				}

				if ( d2 > attached_h_max_dis2 ) continue;

				Size ii_hatbegin( rsd1.attached_H_begin()[ ii ] ), ii_hatend( rsd1.attached_H_end()[ ii ] );
				Size jj_hatbegin( rsd2.attached_H_begin()[ jj ] ), jj_hatend( rsd2.attached_H_end()[ jj ] );
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {

					Size kk_rep = get_countpair_representative_atom( rsd1.type(), kk );
					//returns kk
					// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(kk) ) );
					Real k_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), kk );
					weight = 1.0;
					path_dist = 0;
					if ( cpfxn->count( kk_rep, jj_rep, weight, path_dist ) ) {
						score += score_atom_pair( rsd1, rsd2, kk, jj, k_hyd, j_hyd, emap, weight, d2 );
					}
				}
				for ( Size kk = jj_hatbegin; kk <= jj_hatend; ++kk ) {

					Size kk_rep = get_countpair_representative_atom( rsd2.type(), kk );
					// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(kk) ) );
					Real k_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), kk );

					weight = 1.0;
					path_dist = 0;
					if ( cpfxn->count( ii_rep, kk_rep, weight, path_dist ) ) {
						score += score_atom_pair( rsd1, rsd2, ii, kk, i_hyd, k_hyd, emap, weight, d2 );
					}
				}
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {

					Size kk_rep = get_countpair_representative_atom( rsd1.type(), kk );
					// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(kk) ) );
					Real k_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), kk );
					for ( Size ll = jj_hatbegin; ll <= jj_hatend; ++ll ) {
						Size ll_rep = get_countpair_representative_atom( rsd2.type(), ll );
						// Real l_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(ll) ) );
						Real l_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), ll );

						weight = 1.0;
						path_dist = 0;
						if ( cpfxn->count( kk_rep, ll_rep, weight, path_dist ) ) {
							score += score_atom_pair( rsd1, rsd2, kk, ll, k_hyd, l_hyd, emap, weight, d2 );

						}
					}
				}
			}
		}
	} else {
		Real d2;

		for ( Size ii = 1; ii <= rsd1.nheavyatoms(); ++ii ) {
			for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {

				// Real i_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(ii) ) );
				Real i_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), ii );
				// Real j_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(jj) ) );
				Real j_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), jj );

				Real weight=1.0;
				score += score_atom_pair( rsd1, rsd2, ii, jj, i_hyd, j_hyd, emap, weight, d2 );

				if ( d2 > attached_h_max_dis2 ) continue;

				Size ii_hatbegin( rsd1.attached_H_begin()[ ii ] ), ii_hatend( rsd1.attached_H_end()[ ii ] );
				Size jj_hatbegin( rsd2.attached_H_begin()[ jj ] ), jj_hatend( rsd2.attached_H_end()[ jj ] );
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(kk) ) );
					Real k_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), kk );
					score += score_atom_pair( rsd1, rsd2, kk, jj, k_hyd, j_hyd, emap, weight, d2 );
				}
				for ( Size kk = jj_hatbegin; kk <= jj_hatend; ++kk ) {
					// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(kk) ) );
					Real k_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), kk );
					score += score_atom_pair( rsd1, rsd2, ii, kk, i_hyd, k_hyd, emap, weight, d2 );
				}
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					for ( Size ll = jj_hatbegin; ll <= jj_hatend; ++ll ) {
						// Real k_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd1.xyz(kk) ) );
						Real k_hyd = mp_geometry->f_transition( conf, rsd1.seqpos(), kk);
						// Real l_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( rsd2.xyz(ll) ) );
						Real l_hyd = mp_geometry->f_transition( conf, rsd2.seqpos(), ll);

						score += score_atom_pair( rsd1, rsd2, kk, ll, k_hyd, l_hyd, emap, weight, d2 );
					}
				}
			}
		}
	}
	// TR << "fa_imm_elec terms:: ";
	// TR << rsd1.name3() << ' ' << rsd2.name3() << ' ' << score << std::endl;
	emap[ core::scoring::fa_imm_elec ] += score;
}

//eval_inrares_energy in fa_elec is to calculate fa_intra_elec
////////////////////////////////////////////////////////////////////////////
Real
ImplicitMembraneElecEnergy::score_atom_pair(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::Size const at1,
	core::Size const at2,
	core::Real const at1_hyd,
	core::Real const at2_hyd,
	EnergyMap &,
	core::Real const cpweight,
	core::Real & d2
) const {

	core::Real energy;
	energy = cpweight * membrane_coulomb_->eval_atom_atom_fa_elecE(
		rsd1.xyz(at1), rsd1.atomic_charge(at1), at1_hyd,
		rsd2.xyz(at2), rsd2.atomic_charge(at2), at2_hyd, d2);

	return energy;

}

////////////////////////////////////////////////////////////////////////////
void
ImplicitMembraneElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,// domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {

	using namespace etable::count_pair;
	if ( ! pose.energies().use_nblist_auto_update() ) return;
	core::conformation::Conformation const & conf( pose.conformation() );
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );


	// what is my charge?
	Size const i( atom_id.rsd() );
	Size const ii( atom_id.atomno() );
	conformation::Residue const & irsd( pose.residue( i ) );
	Real const ii_charge( irsd.atomic_charge( ii ) );
	// Real const ii_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( pose.residue( atom_id.rsd() ).xyz( atom_id.atomno() )));
	Real const ii_hyd( mp_geometry->f_transition( conf, pose.residue(atom_id.rsd()).seqpos(), atom_id.atomno() ) );
	core::Vector dii_hyd_df1 ( mp_geometry->f_transition_f1( conf, pose.residue(atom_id.rsd()).seqpos(), atom_id.atomno() ) );
	core::Vector dii_hyd_df2 ( mp_geometry->f_transition_f2( conf, pose.residue(atom_id.rsd()).seqpos(), atom_id.atomno() ) );

	if ( ii_charge == 0.0 ) return;

	Vector const & ii_xyz( irsd.xyz(ii) );
	bool const ii_isbb( irsd.atom_is_backbone( ii ) );

	debug_assert( pose.energies().use_nblist() );
	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::ELEC_NBLIST ) );
	AtomNeighbors const & nbrs( nblist.atom_neighbors(i,ii) );

	weight_triple wtrip;
	wtrip.wbb_bb_ = weights[ fa_elec ] + weights[ fa_elec_bb_bb ];
	wtrip.wbb_sc_ = weights[ fa_elec ] + weights[ fa_elec_bb_sc ];
	wtrip.wsc_bb_ = weights[ fa_elec ] + weights[ fa_elec_bb_sc ];
	wtrip.wsc_sc_ = weights[ fa_elec ] + weights[ fa_elec_sc_sc ];
	wtrip.w_intra_ = weights[ fa_intra_elec ];

	for ( auto const & nbr : nbrs ) {
		Size const j( nbr.rsd() );
		Size const jj( nbr.atomno() );
		conformation::Residue const & jrsd( pose.residue( j ) );

		Real const jj_charge( jrsd.atomic_charge(jj) );
		// Real const jj_hyd( pose.conformation().membrane_info()->implicit_lipids()->f_hydration( pose.residue( nbr.rsd() ).xyz( nbr.atomno() )));
		Real const jj_hyd( mp_geometry->f_transition( conf, pose.residue( nbr.rsd() ).seqpos(), nbr.atomno()) );
		core::Vector djj_hyd_df1( mp_geometry->f_transition_f1( conf, pose.residue( nbr.rsd() ).seqpos(), nbr.atomno()) );
		core::Vector djj_hyd_df2( mp_geometry->f_transition_f2( conf, pose.residue( nbr.rsd() ).seqpos(), nbr.atomno()) );

		if ( jj_charge == 0.0 ) continue; /// should prune out such atoms when constructing the neighborlist!
		Vector const & jj_xyz( jrsd.xyz( jj ) );
		Vector f2 = ( ii_xyz - jj_xyz );
		Real const dis2( f2.length_squared() );
		Real const dE_dr_over_r = nbr.weight() * membrane_coulomb_->eval_dfa_elecE_dr_over_r( dis2, ii_charge, jj_charge, ii_hyd, jj_hyd );
		if ( dE_dr_over_r == 0.0 ) continue;

		Real sfxn_weight = elec_weight( i == j, ii_isbb, jrsd.atom_is_backbone( jj ), wtrip );
		Vector f1 = ii_xyz.cross( jj_xyz );
		f1 *= dE_dr_over_r * sfxn_weight;
		f2 *= dE_dr_over_r * sfxn_weight;
		f1 += nbr.weight() * sfxn_weight * membrane_coulomb_->eval_dfa_elecE_df( dis2, ii_charge, jj_charge, ii_hyd, jj_hyd, dii_hyd_df1, djj_hyd_df1 );
		f2 += nbr.weight() * sfxn_weight * membrane_coulomb_->eval_dfa_elecE_df( dis2, ii_charge, jj_charge, ii_hyd, jj_hyd, dii_hyd_df2, djj_hyd_df2 );

		F1 += f1;
		F2 += f2;
	}
}

void
ImplicitMembraneElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{}

core::Size
ImplicitMembraneElecEnergy::version() const
{
	return 1; // Initial versioning
}

// } // elec
} // energy_methods
} // core

