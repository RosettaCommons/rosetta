// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:	 core/scoring/facts/FACTSPotential.cc
// @brief:	The definitions of 3 classes of the FACTS algorithm resides here (see devel/khorvash/FACTSPotential.hh
// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSResidue.hh>
#include <core/scoring/facts/FACTSPose.hh>
#include <core/scoring/facts/FACTSPotential.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <math.h>
#include <stdio.h>
#include <cassert>
#include <utility/assert.hh>

static basic::Tracer TR("core.scoring.FACTSPotential");

using namespace std;

namespace core {
namespace scoring {


// fast math -- from https://code.google.com/p/fastapprox/downloads/detail?name=fastapprox-0.3.2.tar.gz
inline float
fastpow2 (float p) {
	float offset = (p < 0.0) ? 1.0f : 0.0f;
	float clipp = (p < -126.0) ? -126.0f : p;
	int w = (int)(clipp);
	float z = clipp - (float)(w) + offset;
	union { uint32_t i; float f; } v = { (uint32_t) ( (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) ) };

	return v.f;
}

inline float
fastlog2 (float x) {
	union { float f; uint32_t i; } vx = { x };
	union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	float y = vx.i;
	y *= 1.1920928955078125e-7f;

	return y - 124.22551499f
					 - 1.498030302f * mx.f
					 - 1.72587999f / (0.3520887068f + mx.f);
}

inline float fastexp (float p) { return fastpow2 (1.442695040f * p); }
inline float fastpow (float x, float p) { return fastpow2 (p * fastlog2 (x)); }

/**************************************************************************************************/
///
///		@brief: The FACTSPotential class provides all the functions, constants, and parameters
///				common to all atoms required to calculate the free energy of solvation of a
///				(macro)molecule embedded in a continuum solvent using FACTS method
///
/**************************************************************************************************/
FACTSPotential::FACTSPotential (): MultiplicitiveFactor_(332.07156)
{
	set_default();
}

void
FACTSPotential::set_default()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Kappa_ = option[ score::facts_kappa ]();
	GBpair_cut_ = option[ score::facts_GBpair_cut ]();
	Tau_ = (1.0/option[ score::facts_die ]() ) - (1.0/20.0);
	inv_die_ = 1.0/option[ score::facts_die ]();

	do_apprx_ = false;
	eq_type_ = option[ score::facts_eq_type ]();
	if( eq_type_.compare("apprx") == 0 ) do_apprx_ = true;

	intrascale_by_level_ = option[ score::facts_intrascale_by_level ]();

	saltbridge_correction_ = option[ score::facts_saltbridge_correction ]();
	utility::vector1<Real> dshift = option[ score::facts_dshift ]();
	dshift2_bb_ = dshift[1]*dshift[1];
	dshift2_bs_ = dshift[2]*dshift[2];
	dshift2_sc_ = dshift[3]*dshift[3];
	dshift2_saltbridge_ = dshift[4]*dshift[4];

	assert( option[score::facts_adjbb_elec_scale ]().size() == 5 );
	assert( option[score::facts_adjbb_solv_scale ]().size() == 5 );
	assert( option[score::facts_adjbs_elec_scale ]().size() == 5 );
	assert( option[score::facts_adjbs_solv_scale ]().size() == 5 );

	adjbb_elec_scale_.resize( 5 );
	adjbb_elec_scale_ = option[ score::facts_adjbb_elec_scale ]();
	adjbb_solv_scale_.resize( 5 );
	adjbb_solv_scale_ = option[ score::facts_adjbb_solv_scale ]();

	adjbs_elec_scale_.resize( 5 );
	adjbs_elec_scale_ = option[ score::facts_adjbs_elec_scale ]();
	adjbs_solv_scale_.resize( 5 );
	adjbs_solv_scale_ = option[ score::facts_adjbs_solv_scale ]();

	if( do_apprx_ ) TR << "FACTS Approximation turned on." << endl;

	TR << "FACTS ASP set using: " << option[ score::facts_asp_patch ]()<< endl;
}


//This function "pre-calculates" all the atomic contents (Born radius, sasa,...) and energy values
void FACTSPotential::setup_for_scoring(pose::Pose & pose, bool const & packing) const
{
	Size res1;

	PROF_START( basic::FACTS_GET_ALL_BORN_RADII );
	Size const nres( pose.total_residue() );
	FACTSPoseInfoOP facts_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* >
				( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
		// commenting out below would always set whole residues as enumeration_shell
		//facts_info->update_enumeration_shell( pose, true );
	} else {
		facts_info = new FACTSPoseInfo();
	}

	facts_info->initialize( pose, FACTSrsdtypemap_ );

	// Set enumeration shell / store current coordinate
	TR.Debug << "Enumeration shell: ";
	Size nshell( 0 );
	for ( res1 = 1; res1 <= nres; ++res1 ) {
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		if( packing ){	// Nothing to take care of; just do full enumeration
			facts1.set_enumeration_shell( true );

		} else { // Scoring, minimizing; refresh energy cache for enumeration shell
			if( facts1.enumeration_shell() ){
				facts1.refresh_energy_cache( nres );
				TR.Debug << " " << res1;
				nshell++;
			}
			facts1.store_xyz( pose.residue(res1) );
		}
	}
	TR.Debug << std::endl;
	TR.Debug << "nres_shell/nres_total " << nshell << "/" << nres << std::endl;

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// 1. First get Born radius, Solvation energy, SASA for all atoms
	for ( res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		// Iter res1 over enumeration_shell only
		if( facts1.natoms() == 0 ) continue;

		if( facts1.enumeration_shell() )
			res_res_burial_for_scoring( rsd1, facts1, rsd1, facts1 );

		// changed again from upperedge to edge to take care of enumeration_shell stuffs
		for ( graph::Graph::EdgeListConstIter
			iru	= energy_graph.get_node( res1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( res1 )->const_edge_list_end();
		iru != irue; ++iru ) {
			Size const res2( (*iru)->get_other_ind( res1 ) );
			FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );

			// pass if both are not in enumeration shell
			if ( !facts1.enumeration_shell() && !facts2.enumeration_shell() ) continue;

			// pass double counting cases
			if ( facts1.enumeration_shell() && facts2.enumeration_shell() && res2 <= res1 ) continue;

			res_res_burial_for_scoring( rsd1, facts1, pose.residue( res2 ), facts2 );
		}
	}

	// 2. Refresh Born radii / SASA
	for ( res1 = 1; res1 <= nres; ++res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		if( facts1.enumeration_shell() ){
			FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
			get_self_terms( factstype1, facts1, packing );
		}
	}

	// 3. Then Pre-calculate all GB-pair-related parts and store them in FACTSINFO
	Size nenum( 0 );
	Size nenum_full( 0 );
	for ( res1 = 1; res1 <= nres; ++ res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		Residue const & rsd1( pose.residue( res1 ) );

		if( facts1.enumeration_shell() ){
			if( eq_type_.compare("apprx") == 0 ){
				calculate_GBpair_fast( rsd1, rsd1, facts1, facts1 );
			} else if( eq_type_.compare("v1trunk") == 0 ){
				calculate_GBpair_v1trunk( rsd1, rsd1, facts1, facts1 );
			} else {
				calculate_GBpair_exact( rsd1, rsd1, facts1, facts1 );
			}
		}

		// changed again from upperedge to edge to take care of enumeration_shell stuffs
		for ( graph::Graph::EdgeListConstIter
						iru	= energy_graph.get_node( res1 )->const_edge_list_begin(),
						irue = energy_graph.get_node( res1 )->const_edge_list_end();
					iru != irue; ++iru ) {

			Size const res2( (*iru)->get_other_ind( res1 ) );
			Residue const & rsd2( pose.residue( res2 ) );
			FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );

			nenum_full ++;
			// option 1(slow). Enumerate all the shell-nonshell pairs
			if ( !facts1.enumeration_shell() && !facts2.enumeration_shell() ) continue;

			// option 2(fast). Enumerate only shell-shell pairs: this seems like having error on derivatives...
			//if ( !facts1.enumeration_shell() || !facts2.enumeration_shell() ) continue;

			// pass double counting cases
			if ( facts1.enumeration_shell() && facts2.enumeration_shell() && res2 <= res1 ) continue;

			nenum ++;
			if( eq_type_.compare("apprx") == 0 ){
				calculate_GBpair_fast( rsd1, rsd2, facts1, facts2 );
			} else if( eq_type_.compare("v1trunk") == 0 ){
				calculate_GBpair_v1trunk( rsd1, rsd2, facts1, facts2 );
			} else {
				calculate_GBpair_exact( rsd1, rsd2, facts1, facts2 );
			}
		}
	}
	TR.Debug << "nrespair for GBpair: " << nenum << "/" << nenum_full << std::endl;

	// 3. Finally store everything into pose
	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );
}// END setup_for_scoring


// This function evaluates Ai (volume of each atom) and components of Bi (symmetry of each atom),
// which are converted to Born Radius & SASA in the next step
void FACTSPotential::res_res_burial(
						conformation::Residue const & rsd1,
						FACTSResidueInfo & facts1,
						conformation::Residue const & rsd2,
						FACTSResidueInfo const & facts2
						) const
{
	bool const same_res( rsd1.seqpos() == rsd2.seqpos() );
	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();
	Vector xyz2, xyz1, dyxz;

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

	// START: This part calculates the volume and symmetry of the solute around atom atm
	// and maintains it in Ai_ and nmtr_ & dnmtr_
	for ( Size atm1 = 1; atm1 <= natoms1; ++atm1 ) {
		xyz1 = rsd1.xyz(atm1);

		if ( factstype1->not_using(atm1) ) continue;

		for ( Size atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			if (same_res && (atm1 == atm2 )) continue;

			xyz2 = (rsd2.xyz( atm2 ));
			Vector const dxyz( xyz1 - xyz2 );

			if ( factstype2->not_using(atm2) || factstype2->volume(atm2) < 1e-3) continue;

			Real dis2 = xyz1.distance_squared( xyz2 );

			// take special care for distance to prevent from exploding
			if( dis2 != dis2 ) continue;
			if( dis2 < 0.01 ) dis2 = 0.01;

			Real dis = std::sqrt(dis2);

			Real CutOff_sqr = 64.0; //factstype1->COradius2(atm1);

			//this is a redundant check as there is a stricter check on the next if statement ...
			if( dis2 >= CutOff_sqr ) continue;

			// Equation 5 on page 704 of FACTS paper
			Real theta_sqrt = 1.0 - (dis2 / CutOff_sqr);
			Real thetaij = theta_sqrt*theta_sqrt;

			// The term within the sigma of equation 3 on page 704 of FACTS'
			Real Vi = factstype2->volume(atm2) * thetaij;

			// 1. Ai in equation 3 of page 704 of FACTS paper
			facts1.Ai_[atm1] += Vi;

			// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
			facts1.nmtr_[atm1]	+= (Vi/dis2)*dxyz;
			facts1.dnmtr_[atm1] +=	Vi/dis;
		}
	}
}//end res_res_burial

// This function has same logic with res_res_burial, but modified for efficient scoring
// and for derivative evaluation
void FACTSPotential::res_res_burial_for_scoring(
						conformation::Residue const & rsd1,
						FACTSResidueInfo & facts1,
						conformation::Residue const & rsd2,
						FACTSResidueInfo & facts2
						) const
{
	bool const same_res( rsd1.seqpos() == rsd2.seqpos() );

	Real thetaij, theta_sqrt;
	Real MAX_SELFDCUT2 = min(100.0, GBPair_cut()*GBPair_cut() );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
	Size natoms1 = factstype1->natoms();
	Size natoms2 = factstype2->natoms();

	runtime_assert( natoms1 == rsd1.natoms() );
	runtime_assert( natoms2 == rsd2.natoms() );

	// START: This part calculates the volume and symmetry of the solute around atom atm
	// and maintains it in Ai_ and nmtr_ & dnmtr_
	for ( Size atm1 = 1; atm1 <= natoms1; ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		if ( factstype1->not_using(atm1) ) continue;

		Real CutOff_sqr1 = min( factstype1->COradius2(atm1), MAX_SELFDCUT2 );
		Real const &V1 = factstype1->volume(atm1);

		// Iterate only for upper diagonal :)
		for ( Size atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			if ( same_res && (atm1 >= atm2) ) continue;
			if ( factstype2->not_using(atm2) ) continue;

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Vector const dxyz( xyz1 - xyz2 );
			Real dis2 = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];

			// take special care for distance to prevent from exploding
			if( dis2 != dis2 || dis2 > MAX_SELFDCUT2 ){
				continue;
			} else if( dis2 < 0.01 ) {
				dis2 = 0.01;
			}

			Real dis = std::sqrt(dis2);
			Real CutOff_sqr2 = min( factstype2->COradius2(atm2), MAX_SELFDCUT2 );

			// Consideration for the first atom
			if( facts1.enumeration_shell() && dis2 <= CutOff_sqr1 ){
	
				// Equation 5 on page 704 of FACTS paper
				theta_sqrt = 1.0 - (dis2 / CutOff_sqr1);
				thetaij = theta_sqrt*theta_sqrt;
				// The term within the sigma of equation 3 on page 704 of FACTS'
				Real Vi = factstype2->volume(atm2) * thetaij;
				
				// 1. Ai in equation 3 of page 704 of FACTS paper
				facts1.Ai_[atm1] += Vi;
				
				// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
				facts1.nmtr_[atm1] += (Vi/dis2)*dxyz;
				facts1.dnmtr_[atm1] += Vi/dis;
			}

			// Consideration for the second atom
			if( facts2.enumeration_shell() &&	dis2 <= CutOff_sqr2 ){
				
				// Equation 5 on page 704 of FACTS paper
				theta_sqrt = 1.0 - (dis2 / CutOff_sqr2);
				thetaij = theta_sqrt*theta_sqrt;
				
				// The term within the sigma of equation 3 on page 704 of FACTS'
				Real const Vi = V1*thetaij;
				
				// 1. Ai in equation 3 of page 704 of FACTS paper
				facts2.Ai_[atm2] += Vi;
				
				// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
				facts2.nmtr_[atm2]	-= (Vi/dis2)*dxyz;
				facts2.dnmtr_[atm2] +=	Vi/dis; // the dnmtr of equation 4 on page 704 of FACTS paper
			}
		}
	}
}//end res_res_burial_for_scoring

// Converts Ai & Bi (actually its nmtr & dnmtr) into Born Radius & SASA
void FACTSPotential::get_self_terms(
						FACTSRsdTypeInfoCOP factstype1,
						FACTSResidueInfo & facts1,
						bool const packing
						) const
{
	for( Size atm1 = 1; atm1<=factstype1->natoms(); atm1++ ){
		if ( factstype1->not_using(atm1) ) continue;

		// Confirm B denominator isn't any strange value
		if( !(facts1.dnmtr(atm1) >= 1.0 ) ) facts1.dnmtr_[atm1] = 1.0;

		// Bi needs to be calculated after all vectors are collected in res_res_burial
		// Bi in equation 4 of page 704 of FACTS paper
		facts1.Bi_[atm1] = (facts1.nmtr(atm1)).norm() / facts1.dnmtr(atm1);

		// 1. Born Radius & Self Polar energy!
		facts1.Ci_[atm1] = facts1.Ai(atm1) + (factstype1->b1( atm1 ) * facts1.Bi(atm1)) +
			(factstype1->b2(atm1)	* facts1.Ai(atm1) * facts1.Bi(atm1) );

		Real expterm = exp( -factstype1->a2(atm1) * (facts1.Ci(atm1) - factstype1->a3(atm1)) );
		Real tmp = factstype1->a1(atm1)/(1.0 + expterm);

		// Equation 7 of page 704 of FACTS paper
		facts1.esolvE_[atm1] = factstype1->a0(atm1) + tmp;

		facts1.BR_[atm1] = -0.5*Tau()*MultiplicitiveFactor()/facts1.esolvE_[atm1];
		// Correction for charged hydrogen
		if( factstype1->is_chargedH(atm1) ) facts1.BR_[atm1] *= saltbridge_correction_;

		// Save partial derivative for nonpolar term here
		if( !packing )
			facts1.dG_dCi_[atm1] = factstype1->a2(atm1)*expterm*tmp/( 1.0 + expterm );

		// 2. SASA!
		// Equation 10 of page 706 of FACTS paper
		facts1.Di_[atm1] = facts1.Ai(atm1) + (factstype1->d1( atm1 ) * facts1.Bi(atm1)) +
			(factstype1->d2(atm1)	* facts1.Ai(atm1) * facts1.Bi(atm1) );

		// Equation 11 of page 706 of FACTS paper
		expterm = exp( -factstype1->c2(atm1) * (facts1.Di(atm1) - factstype1->c3(atm1)) );
		tmp = factstype1->c1(atm1)/(1.0 + expterm);

		facts1.sasa_[atm1] = factstype1->c0(atm1) + tmp;

		// Save partial derivative for nonpolar term here
		if( !packing )
			facts1.dSA_dDi_[atm1] = factstype1->c2(atm1)*expterm*tmp/( 1.0 + expterm );

		// Assert before proceed
		runtime_assert( facts1.BR(atm1) < 1.0e10 );
		//if( facts1.BR(atm1) > 1.0e10 )
		//	std::cout << atm1 << " " << tmp << " " << factstype1->a0(atm1) << " " << facts1.BR(atm1) << std::endl;
		runtime_assert( facts1.BR(atm1) > 1.0e-1 );
		runtime_assert( std::abs( facts1.esolvE(atm1) ) < 1.0e4 );
		runtime_assert( facts1.sasa(atm1) < 1.0e10 );
		runtime_assert( facts1.sasa(atm1) >= 0.0 );
	}
}

// Collect the relationship information to correct the scales
void FACTSPotential::atompair_scale( FACTSRsdTypeInfoCOP factstype1, 
		FACTSRsdTypeInfoCOP , //factstype2,
		scoring::etable::count_pair::CountPairFunctionCOP cpfxn14,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Size const atm1,
		Size const atm2,
		Real &scale_solv,
		Real &scale_elec,
		bool &self_pair,
		bool const same_res,
		bool const adjacent ) const
{
	// Default definitions
	scale_solv = 1.0;
	scale_elec = 1.0;
	self_pair = false;

	// Code becomes messy just because Rosetta treats intra and inter in different way...
	// but the rule is not complicated actually :)
	Size path_dist( 0 );
	Real cpweight4( 0.0 );

	bool const is_atm1_CO =
		( rsd1.atom_is_backbone(atm1) &&
			( rsd1.atom_name(atm1).compare(" C  ") == 0 || rsd1.atom_name(atm1).compare(" O  ") == 0 ));
	bool const is_atm2_NH =
		( rsd2.atom_is_backbone(atm2) &&
			( rsd2.atom_name(atm2).compare(" N  ") == 0 || rsd2.atom_name(atm2).compare(" H  ") == 0 ));

	if( same_res ){
		scale_solv = factstype1->intra_solv_scale(atm1,atm2);
		scale_elec = factstype1->intra_elec_scale(atm1,atm2);

		if( scale_solv >= 0.0 ) self_pair = true;

	} else if ( adjacent ){
		// i) consider as intrashell for i-1 C=O & i+1 N-H dipoles
		// since they are exactly determined by i-th phi/psi (unless omega is cis)

		// Turn off elec for angle pairs
		if( path_dist <= 2 ){
			self_pair = true;
			scale_elec = 0.0;
			scale_solv = 1.0;

			// Rule from 1-4
			// BB-BB
		} else if (rsd1.atom_is_backbone(atm1) && rsd2.atom_is_backbone(atm2) ) {
			self_pair = true;

			if( path_dist <= 5 ){
				scale_solv = adjbb_solv_scale( path_dist - 2 );
				scale_elec = adjbb_elec_scale( path_dist - 2 );

			} else { // Decoupled case ( spanned by more than phi+psi )
				scale_solv = adjbb_solv_scale( 5 );
				scale_elec = adjbb_elec_scale( 5 );
			}
	
		// BB-SC
		} else if (( !rsd1.atom_is_backbone(atm1) && rsd2.atom_is_backbone(atm2)) ||
				( rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2)) ){
			// i) Coupled atoms
			if( is_atm1_CO || is_atm2_NH ){
				self_pair = true;

				// option1. Use path_dist b/w sidechain to its C-alpha
				if( intrascale_by_level_ ){
					if( is_atm1_CO && rsd2.has( "CA" ) ){ // Sidechain for rsd2
						cpfxn14->count( atm1, atm2, cpweight4, path_dist );

					} else if( is_atm2_NH && rsd1.has( "CA" ) ){ // Sidechain for rsd1
						Size const i_CA1( rsd1.atom_index( "CA" ) );
						cpfxn14->count( i_CA1, atm2, cpweight4, path_dist );
					}

					path_dist ++; // To convert indexing start from Gamma
				}

				if( path_dist <= 5 ){
					scale_solv = adjbs_solv_scale( path_dist - 2 );
					scale_elec = adjbs_elec_scale( path_dist - 2 );
				} else {
					scale_solv = adjbs_solv_scale( 4 );
					scale_elec = adjbs_elec_scale( 4 );
				}

			// ii) Otherwise apply scale given as last parameter (typically 1.0)
			} else {
				scale_solv = adjbs_solv_scale( 5 );
				scale_elec = adjbs_elec_scale( 5 );
			}
		}

	} else { // 2-residue separation
		if(( rsd1.atom_name(atm1).compare(" C  ") == 0 || rsd1.atom_name(atm1).compare(" O  ") == 0 ) &&
			 ( rsd2.atom_name(atm2).compare(" N  ") == 0 || rsd2.atom_name(atm2).compare(" H  ") == 0 ) ){
			self_pair = true;
			scale_solv = adjbb_solv_scale( 4 );
			scale_elec = adjbb_elec_scale( 4 );
		}
	}
}

// Calculate polar interaction between rsd1 & rsd2 using Born Radius information
void FACTSPotential::calculate_GBpair_exact(
							conformation::Residue const & rsd1,
							conformation::Residue const & rsd2,
							FACTSResidueInfo & facts1,
							FACTSResidueInfo & facts2
							) const
{
	using namespace core::scoring::etable::count_pair;

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool const adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real Esolv( 0.0 ), Eelec( 0.0 ), Esolv_self( 0.0 ), Esolv_pair( 0.0 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
	CountPairFunctionOP cpfxn13 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	// Special dshift for salt-bridging pairs
	bool const is_res1_negative( rsd1.aa() == core::chemical::aa_glu || rsd1.aa() == core::chemical::aa_asp );
	bool const is_res2_negative( rsd2.aa() == core::chemical::aa_glu || rsd2.aa() == core::chemical::aa_asp );
	bool const is_res1_positive( rsd1.aa() == core::chemical::aa_arg || rsd1.aa() == core::chemical::aa_lys );
	bool const is_res2_positive( rsd2.aa() == core::chemical::aa_arg || rsd2.aa() == core::chemical::aa_lys );

	bool is_saltbridge_pair( false );
	if( ( is_res1_negative && is_res2_positive ) || (is_res1_positive && is_res2_negative ) ) is_saltbridge_pair = true;

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		if( factstype1->not_using(atm1) || !factstype1->charged(atm1) ) continue;
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			if( factstype2->not_using(atm2) || !factstype2->charged(atm2) ) continue;
			Real const &q2 = factstype2->q( atm2 );

			// Get scaling factor for atm1&atm2
			Real scale_solv( 1.0 ), scale_elec( 1.0 );
			bool self_pair( false );

			if( same_res || adjacent ||
					rsd2.seqpos() - rsd1.seqpos() == 2 ){
				atompair_scale( factstype1, factstype2, cpfxn14, 
												rsd1, rsd2, atm1, atm2, 
												scale_solv, scale_elec, self_pair, same_res, adjacent );
			}

			/// 0. Start evaluation: Common for elec and solv
			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );
			if ( !(dis2 < cut_off_square) ) continue;

			Real const dis = std::sqrt(dis2);
			Vector const dxyz = xyz1 - xyz2;

			// Shift function (required for truncation at cut_off)
			Real const sf1 = 1.0 - dis2/cut_off_square;
			Real const sf2 = sf1*sf1;

			// 1. GB part
			Real const BRi = facts1.BR(atm1);
			Real const BRj = facts2.BR(atm2);
			Real const BRij = BRi*BRj;

			// Distance shift for GB formula - only works on > 1-5 pairs
			Real dshift2( 0.0 );
			if( is_saltbridge_pair && !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_saltbridge_;
			} else if ( rsd1.atom_is_backbone(atm1) && rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_bb_;
			} else if ( !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_sc_;
			} else {
				dshift2 = dshift2_bs_;
			}
			if( self_pair )	dshift2 = 0.0;

			bool fully_exposed( false );

			// Set boundary for extreme case: converge to Coulomb term
			// In all case, self pair should use the exact formula
			Real const arg = MultiplicitiveFactor()*q1*q2;
			Real fsolv( 0.0 );
			Real tmp1 = dis2/Kappa();
			Real tmp2 = exp(-tmp1/BRij);
			//Real tmp2 = (Real)fastexp((float)(-tmp1/BRij));
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;

			if( !self_pair && tmp3 < dis2 ){
				fully_exposed = true;
				tmp1 = 0.0; tmp2 = 0.0; tmp3 = dis2;
				if( dis > 1.0 ){
					fsolv = 1.0/dis;
				} else {
					fsolv = -dis+2.0;
				}
			} else {
				fsolv = 1.0/std::sqrt(tmp3);
			}

			fsolv *= scale_solv*arg*Tau();

			// 2. Coulomb part
			// Take care of short-range singularity
			Real felec( 0.0 );
			if( scale_elec > 0.0 ){
				if( dis > 1.0 ){
					felec = 1.0/dis;
				} else {
					felec = -dis+2.0;
				}
				felec *= scale_elec*arg*inv_die();
			}

			// 3. Derivative stuffs
			Real g1, g2;
			Real dsolv_drij;
			if( fully_exposed ){
				g1 = 0.0; g2 = 0.0;
				if( dis > 1.0 ){
					dsolv_drij = -fsolv/dis2;
				} else {
					dsolv_drij = -scale_solv*arg*Tau()/dis;
				}
			} else {
				g1 = 0.5*fsolv/tmp3;
				g2 = 2.0 - 2.0*tmp2/Kappa(); // for original GB formula and Way2
				dsolv_drij = -g1*g2;
			}

			Real delec_drij( 0.0 );
			if( scale_elec > 0.0 ){
				if( dis > 1.0 ){
					delec_drij = -felec/dis2;
				} else {
					delec_drij = -scale_elec*arg*inv_die()/dis;
				}
			}

			Real dsf2_drij = 4.0*sf1/cut_off_square;
			Real dsolvsf_drij = - dsolv_drij*sf2 + fsolv*dsf2_drij;
			Real delecsf_drij = delec_drij*sf2 - felec*dsf2_drij;

			// Make sure that fpair isn't any weird value
			if( fsolv != fsolv || std::abs(fsolv) > 1.0e6 || felec != felec || std::abs(felec) > 1.0e6 ){
				TR << "Bad pair interaction score(fsolv/felec)! " << fsolv << " " << felec << std::endl;
				continue;
			}

			// 4. Finally sum up to total energy

			// To avoid double counting for intrares
			Real sintra( 1.0 );
			if( same_res ) sintra = 0.5;

			Eelec += sintra*felec*sf2;
			Esolv -= sintra*fsolv*sf2;
			if( self_pair ){
				Esolv_self -= sintra*fsolv*sf2;
			} else {
				Esolv_pair -= sintra*fsolv*sf2;
			}

			// Derivatives
			facts1.elecF2_[atm1] += sintra*delecsf_drij*dxyz;
			facts2.elecF2_[atm2] -= sintra*delecsf_drij*dxyz;

			facts1.solvF2d_[atm1] += sintra*dsolvsf_drij*dxyz;
			facts2.solvF2d_[atm2] -= sintra*dsolvsf_drij*dxyz;

			if( !fully_exposed ){
				facts1.dsolv_dBR_[atm1] -= sintra*sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dsolv_dBR_[atm2] -= sintra*sf2*g1*tmp2*(BRi + tmp1/BRj);
			}
		}//atm2
	}//atm1

	// Finally store into facts residue info
	facts1.E_solv_[rsd2.seqpos()] = Esolv;
	facts1.E_solv_self_[rsd2.seqpos()] = Esolv_self;
	facts1.E_solv_pair_[rsd2.seqpos()] = Esolv_pair;
	facts1.E_elec_[rsd2.seqpos()] = Eelec;
} //END FACTSPotential::calculate_GBpair

// Calculate polar interaction between rsd1 & rsd2 using Born Radius information
void FACTSPotential::calculate_GBpair_fast(
						 conformation::Residue const & rsd1,
						 conformation::Residue const & rsd2,
						 FACTSResidueInfo & facts1,
						 FACTSResidueInfo & facts2
						 ) const
{
	using namespace core::scoring::etable::count_pair;

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool const adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real Esolv( 0.0 ), Eelec( 0.0 ), Esolv_self( 0.0 ), Esolv_pair( 0.0 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
	CountPairFunctionOP cpfxn13 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	// Special dshift for salt-bridging pairs
	bool const is_res1_negative( rsd1.aa() == core::chemical::aa_glu || rsd1.aa() == core::chemical::aa_asp );
	bool const is_res2_negative( rsd2.aa() == core::chemical::aa_glu || rsd2.aa() == core::chemical::aa_asp );
	bool const is_res1_positive( rsd1.aa() == core::chemical::aa_arg || rsd1.aa() == core::chemical::aa_lys );
	bool const is_res2_positive( rsd2.aa() == core::chemical::aa_arg || rsd2.aa() == core::chemical::aa_lys );

	bool is_saltbridge_pair( false );
	if( ( is_res1_negative && is_res2_positive ) || (is_res1_positive && is_res2_negative ) ) is_saltbridge_pair = true;

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		if( factstype1->not_using(atm1) || !factstype1->charged(atm1) ) continue;
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			if( factstype2->not_using(atm2) || !factstype2->charged(atm2) ) continue;
			Real const &q2 = factstype2->q( atm2 );

			// Get scaling factor for atm1&atm2
			Real scale_solv( 1.0 ), scale_elec( 1.0 );
			bool self_pair( false );
			if( same_res || adjacent )
				atompair_scale( factstype1, factstype2, cpfxn14, 
												rsd1, rsd2, atm1, atm2, 
												scale_solv, scale_elec, self_pair, same_res, adjacent );

			/// 0. Start evaluation: Common for elec and solv
			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );
			if ( !(dis2 < cut_off_square) ) continue;

			Real const dis = std::sqrt(dis2);
			Vector const dxyz = xyz1 - xyz2;

			// Shift function (required for truncation at cut_off)
			Real const sf1 = 1.0 - dis2/cut_off_square;
			Real const sf2 = sf1*sf1;

			// 1. GB part
			Real const BRi = facts1.BR(atm1);
			Real const BRj = facts2.BR(atm2);
			Real const BRij = BRi*BRj;

			// Distance shift for GB formula - only works on > 1-5 pairs
			// Distance shift for GB formula - only works on > 1-5 pairs
			Real dshift2( 0.0 );
			if( is_saltbridge_pair && !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_saltbridge_;
			} else if ( rsd1.atom_is_backbone(atm1) && rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_bb_;
			} else if ( !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_sc_;
			} else {
				dshift2 = dshift2_bs_;
			}
			if( self_pair )	dshift2 = 0.0;

			Real tmp1, tmp2, tmp3;
			bool fully_exposed( false );

			tmp1 = dis2/Kappa();
			tmp2 = BRij - tmp1;
			tmp3 = dis2 + tmp2 - dshift2;

			// Set boundary for extreme case: converge to Coulomb term
			// In all case, self pair should use the exact formula
			if( !self_pair && tmp3 < dis2 ){
				fully_exposed = true;
				tmp1 = 0.0; tmp2 = 0.0; tmp3 = dis2;
			}

			Real const arg = MultiplicitiveFactor()*q1*q2;
			Real const fsolv = scale_solv*arg*Tau()/std::sqrt(tmp3);

			// 2. Coulomb part
			// Take care of short-range singularity
			Real felec( 0.0 );
			if( scale_elec > 0.0 ){
				if( dis > 1.0 ){
					felec = 1.0/dis;
				} else {
					felec = -dis+2.0;
				}
				felec *= scale_elec*arg*inv_die();
			}

			// 3. Derivative stuffs
			Real g1 = 0.5*fsolv/tmp3;
			Real g2 = 2.0*(1.0 - 1.0/Kappa()); // for original GB formula and Way2
			Real dsolv_drij = g1*g2;
			Real delec_drij( 0.0 );
			if( scale_elec > 0.0 ){
				if( dis > 1.0 ){
					delec_drij = -felec/dis2;
				} else {
					delec_drij = -scale_elec*arg*inv_die()/dis;
				}
			}

			Real dsf2_drij = 4.0*sf1/cut_off_square;
			Real dsolvsf_drij = dsolv_drij*sf2 + fsolv*dsf2_drij;
			Real delecsf_drij = delec_drij*sf2 - felec*dsf2_drij;

			// Make sure that fpair isn't any weird value
			if( fsolv != fsolv || std::abs(fsolv) > 1.0e6 || felec != felec || std::abs(felec) > 1.0e6 ){
				TR << "Bad pair interaction score(fsolv/felec)! " << fsolv << " " << felec << std::endl;
				continue;
			}

			// 4. Finally sum up to total energy

			// To avoid double counting for intrares
			Real sintra( 1.0 );
			if( same_res ) sintra = 0.5;

			Eelec += sintra*felec*sf2;
			Esolv -= sintra*fsolv*sf2;
			if( self_pair ){
				Esolv_self -= sintra*fsolv*sf2;
			} else {
				Esolv_pair -= sintra*fsolv*sf2;
			}

			// Derivatives
			facts1.elecF2_[atm1] += sintra*delecsf_drij*dxyz;
			facts2.elecF2_[atm2] -= sintra*delecsf_drij*dxyz;

			facts1.solvF2d_[atm1] += sintra*dsolvsf_drij*dxyz;
			facts2.solvF2d_[atm2] -= sintra*dsolvsf_drij*dxyz;

			if( !fully_exposed ){
				facts1.dsolv_dBR_[atm1] -= sintra*sf2*g1*BRj;
				facts2.dsolv_dBR_[atm2] -= sintra*sf2*g1*BRi;
			}
		}//atm2
	}//atm1

	// Finally store into facts residue info
	facts1.E_solv_[rsd2.seqpos()] = Esolv;
	facts1.E_solv_self_[rsd2.seqpos()] = Esolv_self;
	facts1.E_solv_pair_[rsd2.seqpos()] = Esolv_pair;
	facts1.E_elec_[rsd2.seqpos()] = Eelec;
} //END FACTSPotential::calculate_GBpair_fast

// Calculate derivatives for both polar & nonpolar interactions
// Note:
// "res_res_burial, get_self_terms, calculate_GBpair" should precede this function,
// otherwise derivative will be inconsistent to current structure
//
// Then derivative for distance dependent part will be took into account during setup_for_scoring
// This is only for Context-dependent part: Born Radius & SASA dependent terms
//

// Calculate polar interaction between rsd1 & rsd2 using Born Radius information
void FACTSPotential::calculate_GBpair_v1trunk(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2,
			FACTSResidueInfo & facts1,
			FACTSResidueInfo & facts2
			) const
{
	using namespace core::scoring::etable::count_pair;

	Real const min_dis(1.5);
	Real const elec_sh_exponent(2.0);
	Real const saltbridge_correction(1.2);

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real Esolv( 0.0 ), Eelec( 0.0 ), Esolv_self( 0.0 ), Esolv_pair( 0.0 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

	CountPairFunctionOP cpfxn13 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);
		if( factstype1->not_using(atm1) || std::fabs( q1 ) < 1.0e-6 ) continue;
		Size path_dist( 0 );

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );
			if( factstype2->not_using(atm2) || std::fabs( q2 ) < 1.0e-6 ) continue;

			// 1. Selfpair definition: up to 1-3 (connected by angle) whatever respair relation is
			Real dumm = 1.0;
			bool self_pair =
				(adjacent && !(cpfxn13->count( atm1, atm2, dumm, path_dist ))) ||
				(same_res && rsd1.path_distance(atm1,atm2) <= 2);

			// Adjacent respair 1-4, 1-5: special care to avoid overlap with Rama term
			// 1-4 weight is 0.0, 1-5 weight is 0.2, and else 1.0
			// This is just to be consistent with the "hacked way" in hack_elec
			Real cpweight = 1.0;
			if( adjacent ){
				bool is_cp = cpfxn14->count( atm1, atm2, cpweight, path_dist );
				if( !self_pair && !is_cp ) cpweight = 0.0; // 1-4
			}

			// Collect the relationship information to correct the scales
			Real scale_solv( 1.0 );
			Real scale_elec( 1.0 );
			if( self_pair ){
				scale_solv = 1.0; scale_elec = 0.0;
			} else if ( same_res ){
				scale_solv = 0.4; scale_elec = 0.0;
			} else if ( adjacent ){
				scale_solv = scale_elec = cpweight;
			}

			// Start evaluation
			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );
			if ( !(dis2 < cut_off_square) ) continue;

			Real dis = std::sqrt(dis2);
			Vector dxyz;
			Real BRi = facts1.BR(atm1);
			Real BRj = facts2.BR(atm2);
			Real dshift2( 0.0 );

			// Fading short-distance solvation effect, to be in harmonied with hack_elec
			// Consider self-energy (and pseudo self-energy b/w 1-2, 1-3) to be free from hacking.
			if( self_pair ){
				dshift2 = 0.0;
				dxyz = xyz1 - xyz2;
			} else {
				if ( dis < min_dis ){
					dis = min_dis;
					dis2 = dis*dis;
					dxyz[0] = dxyz[1] = dxyz[2] = 0.0;
				} else {
					dxyz = xyz1 - xyz2;
				}
				
				dshift2 = std::min( 2.25, dis2 );
				if( factstype1->is_chargedH(atm1) ) BRi *= saltbridge_correction;
				if( factstype2->is_chargedH(atm2) ) BRj *= saltbridge_correction;
			}
			
			// Main run
			Real BRij = BRi*BRj;
			Real tmp1 = dis2/Kappa();
			Real tmp2 = (Real)fastexp((float)(-tmp1/BRij));
			//Real tmp2 = exp((float)(-tmp1/BRij));
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;
			if ( !(tmp3 > 1.0e-3) ) continue;

			Real const arg = MultiplicitiveFactor()*q1*q2;
			Real const fsolv = scale_solv*arg*Tau()/std::sqrt(tmp3);
			Real const felec = scale_elec*arg*inv_die()/dis;

			// Shift function (required for truncation at cut_off)
			Real sf1 = 1.0 - dis2/cut_off_square;
			Real sf2 = sf1*sf1;
			//Real sf_elec = (Real)fastpow( (float)sf1, (float)elec_sh_exponent );
			Real sf_elec = sf1*sf1;

			// Derivative stuffs
			Real g1 = 0.5*fsolv/tmp3;
			Real g2 = 2.0 - 2.0*tmp2/Kappa();
			Real dsolv_drij = g1*g2;
			Real delec_drij = -felec/dis2;

			Real dsf2_drij = 4.0*sf1/cut_off_square;
			Real dsf_elec_drij = -2.0*elec_sh_exponent*sf_elec/(sf1*cut_off_square);

			Real dsolvsf_drij = dsolv_drij*sf2 + fsolv*dsf2_drij;
			Real delecsf_drij = delec_drij*sf_elec + felec*dsf_elec_drij;

			// Make sure that fpair isn't any weird value
			if( fsolv != fsolv || std::abs(fsolv) > 1.0e6 || felec != felec || std::abs(felec) > 1.0e6 ){
	TR << "Bad pair interaction score(fsolv/felec)! " << fsolv << " " << felec << std::endl;
	continue;
			}

			// Store here and reuse at res_pair scoring & derivative call
			if ( same_res ){ // No elec part
				Esolv -= 0.5*fsolv*sf2;
				if( self_pair ){
					Esolv_self -= 0.5*fsolv*sf2;
				} else {
					Esolv_pair -= 0.5*fsolv*sf2;
				}
				
				facts1.solvF2d_[atm1] += dsolvsf_drij*dxyz;
				facts1.dsolv_dBR_[atm1] -= 0.5*sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dsolv_dBR_[atm2] -= 0.5*sf2*g1*tmp2*(BRi + tmp1/BRj);
			} else {
				Esolv -= fsolv*sf2;
				Eelec += felec*sf_elec;
				if( self_pair ){
					Esolv_self -= fsolv*sf2;
				} else {
					Esolv_pair -= fsolv*sf2;
				}
				
				facts1.elecF2_[atm1] += delecsf_drij*dxyz;
				facts2.elecF2_[atm2] -= delecsf_drij*dxyz;
				
				facts1.solvF2d_[atm1] += dsolvsf_drij*dxyz;
				facts2.solvF2d_[atm2] -= dsolvsf_drij*dxyz;
				facts1.dsolv_dBR_[atm1] -= sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dsolv_dBR_[atm2] -= sf2*g1*tmp2*(BRi + tmp1/BRj);
			}
		}//atm2
	}//atm1

	// Finally store into facts residue info
	facts1.E_solv_[rsd2.seqpos()] = Esolv;
	facts1.E_solv_self_[rsd2.seqpos()] = Esolv_self;
	facts1.E_solv_pair_[rsd2.seqpos()] = Esolv_pair;
	facts1.E_elec_[rsd2.seqpos()] = Eelec;
} //END FACTSPotential::calculate_GBpair_v1trunk

void FACTSPotential::setup_for_derivatives( pose::Pose & pose ) const
{
	FACTSPoseInfoOP facts_info;
	Vector virtualcrd;
	Vector cross_v;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	} else {
		facts_info = new FACTSPoseInfo();
	}

	// Check whether information is changed from the given pose - if changed, call setup_for_scoring again
	if( facts_info->is_changed( pose ) ){
		TR.Debug << "Pose changed since last scoring, call setup_for_scoring..." << std::endl;
		setup_for_scoring( pose, false );
		facts_info = static_cast< FACTSPoseInfo* >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	}

	// This is to make sure that context-dependent derivative arrays are initialized at least once
	bool full_update( facts_info->context_derivative_empty() );

	// First update dBi/dX
	for ( Size res1 = 1; res1 <= facts_info->size(); ++res1){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();

		for ( Size atm1 = 1; atm1 <= facts1.natoms(); ++atm1){

			// Warning: Below can be dangerous for SASA...
			if ( factstype1->not_using(atm1) ) continue;

			facts1.dB_dBdnmtr_[atm1] = -facts1.Bi(atm1)/facts1.dnmtr(atm1);
			facts1.dB_dBnmtr_[atm1]	= 1.0/(facts1.nmtr(atm1).norm() * facts1.dnmtr(atm1));
			facts1.dBR_dG_[atm1] = facts1.BR(atm1)/facts1.esolvE(atm1);
		}
	}

	// Iter
	Real MAX_SELFDCUT2 = min(100.0, GBPair_cut()*GBPair_cut() );
	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( Size res1 = 1; res1 <= facts_info->size(); ++res1){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
		core::conformation::Residue const &rsd1 = pose.residue(res1);

		// Pass if res1 is not in enumeration_shell and not full_update
		if( !(full_update || facts1.enumeration_shell() ) ) continue;

		for ( Size atm1 = 1; atm1 <= facts1.natoms(); ++atm1){
			Vector const &crd1 = rsd1.xyz(atm1);
			Real const CutOff_sqr1( min( factstype1->COradius2(atm1), MAX_SELFDCUT2 ) );

			// intra-res
			for ( Size atm2 = atm1+1; atm2 <= facts1.natoms(); ++atm2){
	if( atm1 == atm2 || facts1.restypeinfo()->not_using(atm2) ) continue; 
	
	Vector const &crd2 = rsd1.xyz(atm2);
	Vector const dxyz( crd1 - crd2 );
	
	atom_atom_context_derivative( facts1, facts1, atm1, atm2, dxyz,
							full_update );
			}

			// inter-res
			for ( graph::Graph::EdgeListConstIter
				iru	= energy_graph.get_node( res1 )->const_upper_edge_list_begin(),
				irue = energy_graph.get_node( res1 )->const_upper_edge_list_end();
			iru != irue; ++iru ) {
	Size res2( (*iru)->get_other_ind( res1 ) );
	
	FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
	core::conformation::Residue const &rsd2 = pose.residue(res2);
	
	for ( Size atm2 = 1; atm2 <= facts2.natoms(); ++atm2){
		Vector const &crd2 = rsd2.xyz(atm2);

		Vector const dxyz( crd1 - crd2 );
		Real const dis2 ( dxyz.length_squared() );
		Real const CutOff_sqr2( factstype2->COradius2(atm2) );

		if( dis2 > MAX_SELFDCUT2 || (dis2 > CutOff_sqr1 && dis2 > CutOff_sqr2 )|| facts2.restypeinfo()->not_using(atm2) ) continue;

		Vector dsolv_drij( 0.0 );
		Vector dSA_drij( 0.0 );
		atom_atom_context_derivative( facts1, facts2, atm1, atm2, dxyz,
					full_update );

	}
			} // end inter-res

		} // atm1
	} // res1

	// Save status before storing
	facts_info->context_derivative_empty_ = false;

	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );
}

// Calculate atom_atom context derivatives for "atm1" brought by "atm2"
// All the information is being stored into facts1 & facts2
void
FACTSPotential::atom_atom_context_derivative( FACTSResidueInfo & facts1,
								FACTSResidueInfo & facts2,
								Size const & atm1,
								Size const & atm2,
								Vector const & dxyz,
								bool const full_update
								) const
{
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

	Real const dis2 ( dxyz.length_squared() );
	Real const dis ( std::sqrt(dis2) );
	Real const i_dis2 ( 1.0/dis2 );
	Real const CutOff_sqr1( factstype1->COradius2(atm1) );
	Real const CutOff_sqr2( factstype2->COradius2(atm2) );

	Vector dsolv_drij( 0.0 );
	Vector dSA_drij( 0.0 );

	// 1. derivative for atm1 by atm2
	if( dis2 < CutOff_sqr1 ){
		Real const theta_sqrt = 1.0 - (dis2 / CutOff_sqr1);
		Real const dtheta_tmp1 = factstype2->volume(atm2) * theta_sqrt*theta_sqrt;
		Real const dtheta_tmp2 = 4.0*factstype2->volume(atm2) *theta_sqrt/CutOff_sqr1;

		Real const tmp1 = dtheta_tmp1*i_dis2; // = theta*Vj*xij/r ;
		Real const tmp2 = (tmp1 + dtheta_tmp2)/dis;

		// Derivative for Ai
		Real const dAi_drij = -dtheta_tmp2;

		// Derivative for Bi
		Real const &dB_dBnmtr = facts1.dB_dBnmtr(atm1);
		Real const &dB_dBdnmtr = facts1.dB_dBdnmtr(atm1);
		Real const &dBR_dG = facts1.dBR_dG(atm1);

		Real const arg12 = -2.0*dtheta_tmp1*i_dis2 - dtheta_tmp2;
		Vector const argv1 = i_dis2*arg12*dxyz;
		Real const dBn_drij = argv1[0]*facts1.nmtr(atm1)[0] + argv1[1]*facts1.nmtr(atm1)[1] + argv1[2]*facts1.nmtr(atm1)[2];
		Vector const dBn_drij2 = facts1.nmtr(atm1)*dtheta_tmp1*i_dis2;
		Vector const dBi_drij = (dB_dBnmtr*dBn_drij - dB_dBdnmtr*tmp2)*dxyz + dB_dBnmtr*dBn_drij2;

		// Derivatives for Polar Interaction
		Vector const dCi_drij_for_F2 = dAi_drij*dxyz + factstype1->b1(atm1)*dBi_drij +
			factstype1->b2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);

		dsolv_drij += facts1.dsolv_dBR(atm1)*dBR_dG*facts1.dG_dCi(atm1)*dCi_drij_for_F2;

		// 2-2. Derivatives for NonPolar Interaction
		Vector dDi_drij_for_F2 = dAi_drij*dxyz + factstype1->d1(atm1)*dBi_drij +
			factstype1->d2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);

		dSA_drij += factstype1->alpha(atm1)*facts1.dSA_dDi(atm1)*dDi_drij_for_F2;
	}

	// 2. derivative for atm2 by atm1
	if( dis2 < CutOff_sqr2 ){
		Real const theta_sqrt = 1.0 - (dis2 / CutOff_sqr2);
		Real const dtheta_tmp1 = factstype1->volume(atm1) * theta_sqrt*theta_sqrt;
		Real const dtheta_tmp2 = 4.0*factstype1->volume(atm1) *theta_sqrt/CutOff_sqr2;

		Real const tmp1 = dtheta_tmp1*i_dis2; // = theta*Vj*xij/r ;
		Real const tmp2 = (tmp1 + dtheta_tmp2)/dis;

		// Derivative for Ai
		Real const dAi_drij = -dtheta_tmp2;

		// Derivative for Bi
		Real const &dB_dBnmtr = facts2.dB_dBnmtr(atm2);
		Real const &dB_dBdnmtr = facts2.dB_dBdnmtr(atm2);
		Real const &dBR_dG = facts2.dBR_dG(atm2);

		Real const arg12 = -2.0*dtheta_tmp1*i_dis2 - dtheta_tmp2;
		Vector const argv1 = -i_dis2*arg12*dxyz;
		Real const dBn_drij = argv1[0]*facts2.nmtr(atm2)[0] + argv1[1]*facts2.nmtr(atm2)[1] + argv1[2]*facts2.nmtr(atm2)[2];
		Vector const dBn_drij2 = facts2.nmtr(atm2)*dtheta_tmp1*i_dis2;
		Vector const dBi_drij = -(dB_dBnmtr*dBn_drij - dB_dBdnmtr*tmp2)*dxyz + dB_dBnmtr*dBn_drij2;

		// Derivatives for Polar Interaction
		Vector const dCi_drij_for_F2 = -dAi_drij*dxyz + factstype2->b1(atm2)*dBi_drij +
			factstype2->b2(atm2)*(facts2.Ai(atm2)*dBi_drij - facts2.Bi(atm2)*dAi_drij*dxyz);

		dsolv_drij -= facts2.dsolv_dBR(atm2)*dBR_dG*facts2.dG_dCi(atm2)*dCi_drij_for_F2;

		// 2-2. Derivatives for NonPolar Interaction
		Vector dDi_drij_for_F2 = -dAi_drij*dxyz + factstype2->d1(atm2)*dBi_drij +
			factstype2->d2(atm2)*(facts2.Ai(atm2)*dBi_drij - facts2.Bi(atm2)*dAi_drij*dxyz);

		dSA_drij -= factstype2->alpha(atm2)*facts2.dSA_dDi(atm2)*dDi_drij_for_F2;
	}

	// Assert that derivative isn't any weird value
	Real const drv_dot1( dsolv_drij.dot( dsolv_drij ) );
	Real const drv_dot2( dSA_drij.dot( dSA_drij ) );
	if( drv_dot1 != drv_dot1 || drv_dot1 > 1.0e10 ){
		TR << "Bad solvation derivatives! " << drv_dot1 << std::endl;
		return;
	}

	if( drv_dot2 != drv_dot2 || drv_dot2 > 1.0e10 ){
		TR << "Bad nonpolar derivatives! " << drv_dot2 << std::endl;
		return;
	}

	if( full_update || facts1.enumeration_shell() ){
		facts1.solvF2BR_[atm1] += dsolv_drij;
		facts1.sasaF2_[atm1] += dSA_drij;
	}

	if( full_update || facts2.enumeration_shell() ){
		facts2.solvF2BR_[atm2] -= dsolv_drij;
		facts2.sasaF2_[atm2] -= dSA_drij;
	}

}

// Called at scoring step - for polar energy
// Just reuse scores calculated at setup_for_scoring
void FACTSPotential::evaluate_polar_energy(Residue const & rsd1,
						 FACTSResidueInfo const & facts1,
						 Residue const & rsd2,
						 Real & E_elec,
						 Real & E_solv_self,
						 Real & E_solv_pair
						 ) const {

	E_elec = facts1.E_elec( rsd2.seqpos() );
	E_solv_self = facts1.E_solv_self( rsd2.seqpos() );
	E_solv_pair = facts1.E_solv_pair( rsd2.seqpos() );

	TR.Debug << "Facts elec/solv_self/solv_pair:";
	TR.Debug << " " << std::setw(4) << rsd1.seqpos()	<< " " << std::setw(4) << rsd2.seqpos();
	TR.Debug << " " << std::setw(10) << E_elec << " " << std::setw(10) << E_solv_self;
	TR.Debug << " " << std::setw(10) << E_solv_pair;
	TR.Debug << " " << std::setw(4) << rsd1.name()	<< " " << std::setw(4) << rsd2.name();
	TR.Debug << std::endl;

}

// Called at scoring step - for nonpolar energy
// Just reuse scores calculated at setup_for_scoring
Real FACTSPotential::evaluate_nonpolar_energy(Residue const & rsd1,
								FACTSResidueInfo const & facts1,
								Residue const & rsd2
								) const {
	Real E_SA = 0.0;
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();

	if ( rsd1.seqpos() == rsd2.seqpos() ){
		for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++ atm1 ){
			E_SA += factstype1->alpha(atm1)*facts1.sasa(atm1);
		}

		TR.Debug << "Facts SA:";
		TR.Debug << " " << std::setw(4) << rsd1.seqpos();
		TR.Debug << " " << std::setw(10) << E_SA;
		TR.Debug << " " << rsd1.name();
		TR.Debug << std::endl;
	}

	return E_SA;
}

void FACTSPotential::eval_atom_polar_derivative(
						id::AtomID const & id,
						Real const weight_elec,
						Real const weight_solv,
						pose::Pose const & pose,
						kinematics::DomainMap const &, //domain_map,
						bool const, //exclude_DNA_DNA,
						Vector & F1,
						Vector & F2
						) const
{
	FACTSPoseInfo const & facts_info( static_cast< FACTSPoseInfo const & >
						( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));
	Size const atm1( id.atomno() );
	Size const res1( id.rsd() );

	// Pose stuff - this is not supported currently
	//int const i_map( domain_map( res2 ) );
	//bool const i_fixed( i_map != 0 );

	FACTSResidueInfo const & facts1( facts_info.residue_info( res1 ) );

	// Just copy from the scratch saved in facts1 pre-calculated at setup_for_derivates
	Vector tmpv =
		weight_elec * facts1.elecF2(atm1)
		+ weight_solv * (facts1.solvF2d(atm1) + facts1.solvF2BR(atm1));
	F2 += tmpv;

	Vector const &crd1 = pose.residue(res1).xyz(atm1);
	Vector virtualcrd = -tmpv + crd1;
	F1 += crd1.cross( virtualcrd );

}

// Called during minimization; Just call derivatives calculated at setup_for_derivative
void FACTSPotential::eval_atom_nonpolar_derivative(
							 id::AtomID const & id,
							 Real const weight,
							 pose::Pose const & pose,
							 kinematics::DomainMap const &,// domain_map,
							 bool const, //exclude_DNA_DNA,
							 Vector & F1,
							 Vector & F2
							 ) const
{
	FACTSPoseInfo const & facts_info( static_cast< FACTSPoseInfo const & >
						( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));
	Size const atm1( id.atomno() );
	Size const res1( id.rsd() );

	// Pose stuff - this is not supported currently
	//int const i_map( domain_map( res2 ) );
	//bool const i_fixed( i_map != 0 );

	FACTSResidueInfo const & facts1( facts_info.residue_info( res1 ) );

	// Just copy from the scratch saved in facts1 pre-calculated at setup_for_derivates
	Vector tmpv = weight * facts1.sasaF2(atm1);
	F2 += tmpv;

	Vector const &crd1 = pose.residue(res1).xyz(atm1);
	Vector virtualcrd = -tmpv + crd1;
	F1 += crd1.cross( virtualcrd );
}

/// @note when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
/// this ensures that we use all the information we have to compute the current set of radii
void FACTSPotential::setup_for_packing(
							 pose::Pose & pose,
							 utility::vector1< bool > const & repacking_residues ) const
{
	PROF_START( basic::FACTS_SETUP_FOR_PACKING );

	FACTSPoseInfoOP facts_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	} else {
		facts_info = new FACTSPoseInfo();
		setup_for_scoring( pose, true );
	}

	/// store info about which positions are moving
	facts_info->set_repack_list( repacking_residues );

	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_SETUP_FOR_PACKING );
}

void FACTSPotential::get_template_born_radii(pose::Pose const & pose, FACTSPoseInfo & facts_info) const{
	Size const nres( pose.total_residue() );
	runtime_assert( facts_info.size() == nres );

	for ( Size i=1; i<= nres; ++i ) {
		if ( facts_info.being_packed( i ) ) continue;

		Residue const & rsd1( pose.residue( i ) );
		FACTSResidueInfo & facts1( facts_info.residue_info( i ) );
		FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
		runtime_assert( rsd1.natoms()<1 || std::fabs(facts1.Ai(1)) < 1e-3 );

		for ( Size j=1; j<= nres; ++j ) {
			// we are not using placeholder in FACTS - this can be changed if we want to do "Design"
			//if ( facts_info.being_packed(j) ) {
			//	res_res_burial( rsd1, facts1, facts_info.placeholder_residue(j), facts_info.placeholder_info(j) );
			//} else {
			res_res_burial( rsd1, facts1, pose.residue(j), facts_info.residue_info(j) );
			//}
		}
		get_self_terms( factstype1, facts1, true );
	}

}
/// called eg after a rotamer substitution is accepted during rotamer trials
void FACTSPotential::update_residue_for_packing(pose::Pose & pose,Size const seqpos) const
{
	FACTSPoseInfo & facts_info( static_cast< FACTSPoseInfo & >
						( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) );
	FACTSResidueInfo & facts( facts_info.residue_info( seqpos ) );

	Residue const & rsd( pose.residue( seqpos ) );
	core::chemical::ResidueType const &rsdtype = rsd.type();
	FACTSRsdTypeMap::const_iterator it = FACTSrsdtypemap_.find( &rsdtype );

	if ( it == FACTSrsdtypemap_.end() ) {
		TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
		FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
		rsdtypeinfo->create_info( rsdtype );
		FACTSrsdtypemap_[ &rsdtype ] = rsdtypeinfo;
		it = FACTSrsdtypemap_.find( &rsdtype );
	}

	facts.initialize( rsd, it->second );
	get_single_rotamer_born_radii( rsd, pose, facts_info, facts );
}

void FACTSPotential::get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const {

	FACTSPoseInfo const & facts_info_pose( static_cast< FACTSPoseInfo const & >
					 (pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));

	// this will get cached in the rotamer set
	// this call should initialize the residue_info objects with the appropriate Residue info
	FACTSRotamerSetInfoOP facts_info_rotamers( new FACTSRotamerSetInfo( rotamer_set, FACTSrsdtypemap_ ) );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		get_single_rotamer_born_radii( *rotamer_set.rotamer(n),
																	 pose, facts_info_pose,
																	 facts_info_rotamers->residue_info( n ) );
	}

	rotamer_set.data().set( core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO, facts_info_rotamers );
}

void FACTSPotential::get_single_rotamer_born_radii(Residue const & rsd1,
							 pose::Pose const & pose,
							 FACTSPoseInfo const & facts_info,
							 FACTSResidueInfo	& facts1) 	const
{
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();

	assert( rsd1.natoms()<1 || std::fabs(facts1.Ai(1)) < 1e-3 );
	for (Size res2=1; res2<= pose.total_residue(); ++res2 ) {
		// we are not using placeholder in FACTS - this can be changed if we want to do "Design"
		//if ( facts_info.being_packed( res2 ) ) {
		//	res_res_burial( rsd1, facts1, facts_info.placeholder_residue( res2 ), facts_info.placeholder_info( res2));
		//} else {
		res_res_burial( rsd1, facts1, pose.residue( res2 ), facts_info.residue_info( res2 ) );
		//}
	}//end for loop 1
	get_self_terms( factstype1, facts1, true );
}

// Given precalculated born radius, called at packing - for polar energy
void FACTSPotential::evaluate_polar_otf_energy(Residue const & rsd1,
								 FACTSResidueInfo const & facts1,
								 Residue const & rsd2,
								 FACTSResidueInfo const & facts2,
								 Real & E_elec,
								 Real & E_solv_self,
								 Real & E_solv_pair
								 ) const {

	// Initialize
	E_elec = 0.0; E_solv_pair = 0.0; E_solv_self = 0.0;

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	using namespace core::scoring::etable::count_pair;

	CountPairFunctionOP cpfxn13 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

	// Special dshift for salt-bridging pairs
	bool const is_res1_negative( rsd1.aa() == core::chemical::aa_glu || rsd1.aa() == core::chemical::aa_asp );
	bool const is_res2_negative( rsd2.aa() == core::chemical::aa_glu || rsd2.aa() == core::chemical::aa_asp );
	bool const is_res1_positive( rsd1.aa() == core::chemical::aa_arg || rsd1.aa() == core::chemical::aa_lys );
	bool const is_res2_positive( rsd2.aa() == core::chemical::aa_arg || rsd2.aa() == core::chemical::aa_lys );

	bool is_saltbridge_pair( false );
	if( ( is_res1_negative && is_res2_positive ) || (is_res1_positive && is_res2_negative ) ) is_saltbridge_pair = true;

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);
		if( factstype1->not_using(atm1) || !factstype1->charged(atm1) ) continue;

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );

			if( factstype2->not_using(atm2) || !factstype2->charged(atm2) ) continue;

			// Get scaling factor for atm1&atm2
			Real cpweight4( 0.0 );
			Size path_dist( 0 );
			cpfxn14->count( atm1, atm2, cpweight4, path_dist );
			Real scale_solv( cpweight4 ); Real scale_elec( cpweight4 );
			bool self_pair( false );
			
			if( same_res || adjacent )
				atompair_scale( factstype1, factstype2, cpfxn14, 
												rsd1, rsd2, atm1, atm2, 
												scale_solv, scale_elec, self_pair, same_res, adjacent );

			// 0. Start evaluation: Common for elec and solv
			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );
			if ( !(dis2 < cut_off_square) ) continue;

			Real const dis = std::sqrt(dis2);
			Vector const dxyz = xyz1 - xyz2;

			// Shift function (required for truncation at cut_off)
			Real const sf1 = 1.0 - dis2/cut_off_square;
			Real const sf2 = sf1*sf1;

			// 1. GB part
			Real const BRi = facts1.BR(atm1);
			Real const BRj = facts2.BR(atm2);
			Real const BRij = BRi*BRj;

			// Distance shift for GB formula - only works on > 1-5 pairs
			// WARNING: high dshift will make wrong packing during relax - set upper bound as 1.0
			Real dshift2( 0.0 );
			if( is_saltbridge_pair && !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_saltbridge_;
			} else if ( rsd1.atom_is_backbone(atm1) && rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_bb_;
			} else if ( !rsd1.atom_is_backbone(atm1) && !rsd2.atom_is_backbone(atm2) ){
				dshift2 = dshift2_sc_;
			} else {
				dshift2 = dshift2_bs_;
			}
			if( self_pair )	dshift2 = 0.0;

			// Set boundary for extreme case: converge to Coulomb term
			// In all case, self pair should use the exact formula
			Real const arg = MultiplicitiveFactor()*q1*q2;
			Real fsolv( 0.0 );
			Real tmp1 = dis2/Kappa();
			//tmp2 = exp(-tmp1/BRij);
			Real tmp2 = (Real)fastexp((float)(-tmp1/BRij));
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;

			// Fast version
			//tmp2 = BRij - tmp1;
			//tmp3 = dis2 + tmp2 - dshift2;

			if( !self_pair && tmp3 < dis2 ){
				tmp1 = 0.0; tmp2 = 0.0; tmp3 = dis2;
				if( dis > 1.0 ){
					fsolv = 1.0/dis;
				} else {
					fsolv = -dis+2.0;
				}
			} else {
				fsolv = 1.0/std::sqrt(tmp3);
			}

			fsolv *= scale_solv*arg*Tau();

			// 2. Coulomb part
			// Take care of short-range singularity
			Real felec( 0.0 );
			if( scale_elec > 0.0 ){
				if( dis > 1.0 ){
					felec = 1.0/dis;
				} else {
					felec = -dis+2.0;
				}
				felec *= scale_elec*arg*inv_die();
			}

			// 4. Finally sum up to total energy

			// To avoid double counting for intrares
			Real sintra( 1.0 );
			if( same_res ) sintra = 0.5;

			// Don't use intraelec in packer!
			//if( !same_res )
			E_elec += sintra*sf2*felec;

			if( self_pair ){
				E_solv_self -= sintra*fsolv*sf2;
			} else {
				E_solv_pair -= sintra*fsolv*sf2;
			}
		}
	}
}

} // namespace scoring
} // namespace core
