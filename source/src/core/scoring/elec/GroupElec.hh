// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FA_ElecEnergy.hh
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Phil Bradley, modifed by James Gleixner


#ifndef INCLUDED_core_scoring_elec_GroupElec_hh
#define INCLUDED_core_scoring_elec_GroupElec_hh

// Package headers
//#include <core/scoring/elec/ElecAtom.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers

#include <cmath>

namespace core {
namespace scoring {
namespace elec {

struct ElecGroup
{
  utility::vector1< core::Size > atms;
  utility::vector1< core::Size > comatms;
  Size n_acceptor;
  Size n_donor;
  Real qeps;
};

typedef utility::vector1< ElecGroup > ResElecGroup;
	//typedef utility::pointer::shared_ptr< ResElecGroup const > ResElecGroupCOP;

class GroupElec : public utility::pointer::ReferenceCount
{
public:

	GroupElec( GroupElec const & src );
	GroupElec( methods::EnergyMethodOptions const & options );
	~GroupElec();

	void initialize( etable::coulomb::Coulomb const &coulomb );

	inline
  Real
  eval_respair_group_coulomb(
														 core::conformation::Residue const & rsd1,
														 core::conformation::Residue const & rsd2
														 //hbonds::HBondSet const & hbond_set
														 ) const;

	inline
  void
  eval_respair_group_derivatives(
 	  core::conformation::Residue const & rsd1,
	  core::conformation::Residue const & rsd2,
	  //hbonds::HBondSet const & hbond_set,
	  utility::vector1< DerivVectorPair > & r1_atom_derivs,
	  utility::vector1< DerivVectorPair > & r2_atom_derivs,
		core::Real const elec_weight,
		core::Real & Erespair
 	  ) const;

private:

  void
  build_groupinfo( std::string const group_file,
		   bool const extra = false );

  ResElecGroup const &
  get_group( core::chemical::ResidueType const &rsdtype ) const;

	bool
	fade_hbonding_group_score( ElecGroup const &grp1,
														 ElecGroup const &grp2,
														 Real & group_score,
														 Real & dw_dE ) const;

	Real
	get_grp_countpair( 
										utility::vector1< Size > const &grp1atms,
										utility::vector1< Size > const &grp2atms,
										etable::count_pair::CountPairFunctionCOP cpfxn,
										Size &path_dist
										 ) const;

	inline
	Real
	eval_atompair_deriv(
											core::conformation::Residue const & rsd1,
											core::conformation::Residue const & rsd2,
											core::Size const atm1,
											core::Size const atm2,
											Real &dE_dr
											) const;

	inline
	Real
	eval_atompair_score(
											conformation::Residue const & rsd1,
											conformation::Residue const & rsd2,
											core::Size const atm1,
											core::Size const atm2
											) const;

	inline
  Vector
  get_grpdis2( conformation::Residue const & rsd1,
	       conformation::Residue const & rsd2,
	       utility::vector1< Size > const &com1atms,
	       utility::vector1< Size > const &com2atms
						 //core::Vector &com1,
							 //core::Vector &com2
	       ) const ;

	inline
	Real
	eval_standard_coulomb( Real const &q1, Real const &q2,
												 Real const &dis2, 
												 bool const &eval_deriv,
												 Real &dE_dr
												 ) const;

	inline
	Real
	eval_grp_trunc( bool const &use_switch,
									Real const &grpdis2,
									bool const &eval_deriv,
									Real &dsw_dr
									) const;

	inline std::string fade_type() const { return fade_type_; }

protected:

  inline
  etable::coulomb::Coulomb const &
  coulomb() const {return *coulomb_; }

  etable::coulomb::CoulombCOP coulomb_;

private:
  // mutable std::map< std::string, ElecGroup > rsdgrps_;
  mutable std::map< std::string const , ResElecGroup > rsdgrps_;

	utility::vector1< Real > cpfxn_weight_;

	// for standard coulomb
	std::string fade_type_;
	core::Real fade_param1_;
	core::Real fade_param2_;
  bool fade_hbond_;
	std::string group_file_;
	bool grp_cpfxn_;
	bool use_subtract_, use_shift_;

	core::Real grp_maxdis2_, grp_swdis2_;
	//std::string grp_cpfxn_mode_;

}; // class GroupElec

inline
Real
GroupElec::eval_atompair_score(
					  core::conformation::Residue const & rsd1,
 					  core::conformation::Residue const & rsd2,
						core::Size const atm1,
						core::Size const atm2
						) const
{
	Real d2 = rsd1.xyz(atm1).distance_squared( rsd2.xyz(atm2) );
	if ( use_subtract_ && d2 > coulomb().max_dis2() ) return 0.0;

	Real atompair_score( 0.0 );
	Real const q1( rsd1.atomic_charge( atm1 ) );
	Real const q2( rsd2.atomic_charge( atm2 ) );
	Real dE_dr( 0.0 );

	if( use_subtract_ ){
		atompair_score = coulomb().eval_atom_atom_fa_elecE( rsd1.xyz( atm1 ), q1, 
																												rsd2.xyz( atm2 ), q2 );
	} else {
		atompair_score = eval_standard_coulomb( q1, q2, d2, false, dE_dr );
	}
	return atompair_score;
}

inline
Real
GroupElec::eval_atompair_deriv(
					  core::conformation::Residue const & rsd1,
 					  core::conformation::Residue const & rsd2,
						core::Size const atm1,
						core::Size const atm2,
						Real &dE_dr
						) const
{
	dE_dr = 0.0;
	Real d2 = rsd1.xyz(atm1).distance_squared( rsd2.xyz(atm2) );
	if ( use_subtract_ && d2 > coulomb().max_dis2() ) return 0.0;

	Real atompair_score( 0.0 );
	Real const q1( rsd1.atomic_charge( atm1 ) );
	Real const q2( rsd2.atomic_charge( atm2 ) );

	if( use_subtract_ ){
		atompair_score = coulomb().eval_atom_atom_fa_elecE( rsd1.xyz( atm1 ), q1, 
																												rsd2.xyz( atm2 ), q2 );
		dE_dr = coulomb().eval_dfa_elecE_dr_over_r( d2, q1, q2 );
	} else {
		atompair_score = eval_standard_coulomb( q1, q2, d2, true, dE_dr );
	}
	return atompair_score;
}

inline
Vector
GroupElec::get_grpdis2( 
											 core::conformation::Residue const & rsd1,
											 core::conformation::Residue const & rsd2,
											 utility::vector1< Size > const &com1atms,
											 utility::vector1< Size > const &com2atms
											 //core::Vector &com1,
											 //core::Vector &com2
												) const
{
  Vector com1( 0.0 );
  Vector com2( 0.0 );

  //core::Size n1( 0 ), n2( 0 );
  for( core::Size i = 1; i <= com1atms.size(); ++i ) com1 += rsd1.xyz( com1atms[i] ); 
  for( core::Size i = 1; i <= com2atms.size(); ++i ) com2 += rsd2.xyz( com2atms[i] ); 

  com1 /= (core::Real)( com1atms.size() );
  com2 /= (core::Real)( com2atms.size() );

  Vector dcom = com1 - com2;
  return dcom;
}

inline
Real
GroupElec::eval_grp_trunc( bool const &use_switch,
													 Real const &grpdis2,
													 bool const &eval_deriv,
													 Real &dsw_dr
													 ) const 
{
	dsw_dr = 0.0;
	//Real const &max_dis2( coulomb().max_dis2() );

	if( grpdis2 > grp_maxdis2_ ){ 
		return 0.0;
	} else if( grpdis2 < grp_swdis2_ ){ 
		return 1.0;
	} else {
		Real const &max_dis( coulomb().max_dis() );

		core::Real const grpdis( std::sqrt(grpdis2) );

		core::Real const dis_sw( fade_param1_ );
		core::Real const min_sw( max_dis - dis_sw );

		// switch function stuffs
		core::Real const sw_dis2 = (max_dis - dis_sw)*(max_dis - dis_sw);
		core::Real const R3on_off( grp_maxdis2_ - 3.0*sw_dis2 );
		core::Real Ron_off( 1.0/( grp_maxdis2_ - sw_dis2) );
		core::Real const Ron_off_3( Ron_off*Ron_off*Ron_off );
		core::Real dr1( grp_maxdis2_ - grpdis2);
		core::Real dr2(R3on_off + 2.0*grp_maxdis2_ );

		// shift function
		core::Size const nexp = (core::Size)( fade_param2_ );
		Real arg = 1.0;
		for( core::Size iexp = 1; iexp <= nexp; ++iexp ) 
			arg *= (grpdis-min_sw)/dis_sw * (grpdis-min_sw)/dis_sw;

		Real const sf1 = 1.0 - arg;
		Real const darg = 2.0*nexp*arg/((grpdis-min_sw)*grpdis); 

		// fade function
		core::Real sw( 1.0 );
		if( grpdis2 > sw_dis2 ){
			if( use_switch ){
				sw = dr1*dr1*dr2*Ron_off_3;
			} else {
				sw = sf1*sf1;
			}
		}

		if( eval_deriv ){
			// fade function
			if( grpdis2 > sw_dis2 ){
				if( use_switch ){
					dsw_dr = 4.0*dr1*Ron_off_3*(dr1-dr2);
				} else {
					dsw_dr = -2.0*darg*sf1;
				}
			}
		}
		return sw;
	}
}

// hpark 09/24/2014
// this was made to reproduce physical numbers by following molecular mechanics way;
// original implementation had deviation from the real Coulomb energy when constant-dielectric is turned on
// looks not too different when using distance-dependent dielectric though
inline
Real
GroupElec::eval_standard_coulomb( Real const &q1, Real const &q2,
																	Real const &dis2, 
																	bool const &eval_deriv,
																	Real &dE_dr
																	) const {
	dE_dr = 0.0;
	core::Real const d( std::sqrt( dis2 ) );

	Real const &die = coulomb().die();
	Real const C0( 322.0637 );

	// choose dielectric
	Real fdie, ddie_dr( 0.0 );
	//std::string const &diefunc( coulomb().dielectric_function() );
	if( dis2 < 1.0 ){
		fdie = die;
		ddie_dr = 0.0;
	} else if( !coulomb().no_dis_dep_die() ){  // ddd
		/*
		if( diefunc.compare("sigmoid") == 0 ){
			fdie = coulomb().die_sigmoid( d, eval_deriv, ddie_dr );
		} else if( diefunc.compare("warshel") == 0 ){
			fdie = coulomb().die_warshel( d, eval_deriv, ddie_dr );
		} else {
		*/
		fdie = die*d;
		ddie_dr = die;
		//}
	} else { // constant
		fdie = die;
		ddie_dr = 0.0;
	}

	core::Real e( 0.0 );
	if( dis2 < 1.0 ){
		e = -d+2.0;
	} else { //regular
		e = 1.0/d;
	}

	e *= q1*q2*C0;

	if( eval_deriv ){
		core::Real de_dr( 0.0 ); // on Coulomb part without fdie
		if( dis2 < 1.0 ){ //use linear softening at short-distance
			de_dr = -q1*q2*C0/fdie;
			//dE_dr = de_dr/(d+0.000001); // to prevent from being infinity
			dE_dr = de_dr/d; // to prevent from being infinity
		} else { // regular 
			de_dr = -e/d;
			dE_dr = ( de_dr/fdie - e*ddie_dr/(fdie*fdie) )/d;
		}
	}

	//std::cout << e << " " << d << " " << fdie << std::endl;

	return e/fdie;
}

inline
Real
GroupElec::eval_respair_group_coulomb( 
	core::conformation::Residue const & rsd1,
  core::conformation::Residue const & rsd2
	//hbonds::HBondSet const & hbond_set
) const
{
  using namespace etable::count_pair;

	Real score( 0.0 );
	//Real d2;
	ResElecGroup const &resgrp1 = get_group( rsd1.type() );
	ResElecGroup const &resgrp2 = get_group( rsd2.type() );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	bool const is_bonded = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );

	for ( Size ii = 1; ii <= resgrp1.size(); ++ii ) {
		ElecGroup const &grp1 = resgrp1[ii];
		for ( Size jj = 1; jj <= resgrp2.size(); ++jj ) {
			if( intrares && ii >= jj ) continue;

			ElecGroup const &grp2 = resgrp2[jj];

			core::Vector dcom;
			dcom = get_grpdis2( rsd1, rsd2, grp1.comatms, grp2.comatms
													);

			core::Real const grpdis2 = dcom.length_squared();

			if( use_shift_ && (grpdis2 > coulomb().max_dis2()) ) continue;

			Real grp_cpweight( 1.0 );
			Size path_dist_grp( 0 );
			if( is_bonded && grp_cpfxn_ )
				grp_cpweight = get_grp_countpair( grp1.atms, grp2.atms, cpfxn, 
																					path_dist_grp );

			//core::Real dE_dr( 0.0 );
			core::Real group_score( 0.0 );

			Real dsw_dr( 1.0 ), sw( 1.0 );
			//if( !use_subtract_ )
			if( use_shift_ )
				sw = eval_grp_trunc( false, grpdis2, false, dsw_dr );

			Size const n1( grp1.atms.size() );
			Size const n2( grp2.atms.size() );
			for ( Size kk = 1; kk <= n1; ++kk ){
				core::Size const &atm1( grp1.atms[kk] );

				for ( Size ll = 1; ll <= n2; ++ll ){
					core::Size const &atm2( grp2.atms[ll] );

					bool is_count( true );
					Real atom_cpweight( 1.0 );
					if( is_bonded && !grp_cpfxn_ ){
						Size path_dist( 0 );
						is_count = cpfxn->count( atm1, atm2, atom_cpweight, path_dist );
					}
					if( !is_count ) continue;
					Real atompair_score = eval_atompair_score( rsd1, rsd2, atm1, atm2 );

					group_score += atom_cpweight*atompair_score;
				}
			}


			group_score *= grp_cpweight;

			// converge to constant for hbonding groups
			bool is_hbond_pair = (grp1.n_donor*grp2.n_acceptor + grp2.n_donor*grp1.n_acceptor > 0) && (path_dist_grp > 5);

			Real dw_dE( 0.0 );
			if( is_hbond_pair && fade_hbond_ ) 
				fade_hbonding_group_score( grp1, grp2, group_score, dw_dE );

			group_score *= sw;

			score += group_score;
		} // grp2
	} // grp1

  return score;
}

inline
void
GroupElec::eval_respair_group_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	//hbonds::HBondSet const & hbond_set,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs,
	core::Real const elec_weight,
	core::Real & Erespair
	) const 
{
  using namespace etable::count_pair;

	Erespair = 0.0;
	ResElecGroup const &resgrp1 = get_group( rsd1.type() );
	ResElecGroup const &resgrp2 = get_group( rsd2.type() );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	bool const is_bonded = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );

	// dummy array
	//utility::vector1< Vector > f2r1( r1_atom_derivs.size(), Iv );
	//utility::vector1< Vector > f2r2( r2_atom_derivs.size(), Iv );

	//TR << "!!" << std::endl;

	for ( Size ii = 1; ii <= resgrp1.size(); ++ii ) {
		ElecGroup const &grp1 = resgrp1[ii];
		core::Size const &ncom1 = grp1.comatms.size();

		for ( Size jj = 1; jj <= resgrp2.size(); ++jj ) {
			if( intrares && (ii >= jj ) ) continue;

			ElecGroup const &grp2 = resgrp2[jj];
			core::Size const &ncom2 = grp2.comatms.size();

			core::Vector dcom;
			dcom = get_grpdis2( rsd1, rsd2, grp1.comatms, grp2.comatms
													);
			core::Real const grpdis2 = dcom.length_squared();

			if ( use_shift_ && (grpdis2 > coulomb().max_dis2()) ) continue;
			//if ( !use_subtract && (grpdis2 > coulomb().max_dis2()) ) continue;
			//if ( grpdis2 > coulomb().max_dis2() ) continue;

			Real grp_cpweight( 1.0 );
			Size path_dist_grp( 0 );

			if( is_bonded && grp_cpfxn_ )
				grp_cpweight = get_grp_countpair( grp1.atms, grp2.atms, cpfxn, 
																					path_dist_grp );

			//core::Real grpdis = std::sqrt( grpdis2 );

			core::Size const &n1 = grp1.atms.size();
			core::Size const &n2 = grp2.atms.size();

			// get group-truncation info first
			Real dsw_dr( 1.0 ), sw( 1.0 );
			//if( !use_subtract )
			if( use_shift_ )
				sw = eval_grp_trunc( false, grpdis2, true, dsw_dr );

			// 1. derivative on Coulomb part: dE*sw
			core::Real group_score( 0.0 );

			bool is_hbond_pair = (grp1.n_donor*grp2.n_acceptor + grp2.n_donor*grp1.n_acceptor > 0) && (path_dist_grp > 5);

			// temporary vector for hbonding pair; allocate only if hbonding
			utility::vector1< Vector > v1,v2; 
			Vector Iv( 0.0 );
			if( is_hbond_pair && fade_hbond_ ){
				v1 = utility::vector1< Vector >( n1, Iv );
				v2 = utility::vector1< Vector >( n2, Iv );
			}

			for ( Size kk = 1; kk <= n1; ++kk ){
				core::Size const &atm1( grp1.atms[kk] );
				Vector const & atom1xyz( rsd1.xyz( atm1 ) );
				//Real const q1( rsd1.atomic_charge( atm1 ) );

				for ( Size ll = 1; ll <= n2; ++ll ){
					core::Size const &atm2( grp2.atms[ll] );
					Vector const & atom2xyz( rsd2.xyz( atm2 ) );
					//Real const q2( rsd2.atomic_charge( atm2 ) );

					Vector f2 = ( atom1xyz - atom2xyz );
					//Real const &dis2( f2.length_squared() );
					//if ( use_subtract_ && dis2 > coulomb().max_dis2() ) continue;

					bool is_count( true );
					Real atom_cpweight( 1.0 );
					if( is_bonded && !grp_cpfxn_ ){
						Size path_dist( 0 );
						is_count = cpfxn->count( atm1, atm2, atom_cpweight, path_dist );
					}

					if( !is_count ) continue;

					atom_cpweight = 1.0;
					Real dE_dr( 0.0 );
					Real atompair_score = eval_atompair_deriv( rsd1, rsd2, atm1, atm2, dE_dr );

					Real sfxn_weight = atom_cpweight*elec_weight;

					group_score += atom_cpweight*atompair_score;

					f2 *= dE_dr*sw*sfxn_weight*grp_cpweight;
					Vector const f1 = atom1xyz.cross( -f2 );

					if( is_hbond_pair && fade_hbond_ ){
						v1[kk] += f2;
						v2[ll] -= f2;
					} else {
						//f2r1_[ atm1 ] += f2*grp_cpweight;
						//f2r2_[ atm2 ] -= f2*grp_cpweight;
						r1_atom_derivs[ atm1 ].f2() += f2;
						r2_atom_derivs[ atm2 ].f2() -= f2;
						r1_atom_derivs[ atm1 ].f1() += f1;
						r2_atom_derivs[ atm2 ].f1() -= f1;
					}
				}
			}

			group_score *= grp_cpweight;
			Erespair += group_score;

			// converge to constant for hbonding groups
			Real dw_dE( 1.0 );
			bool do_fade( false );

			if( is_hbond_pair && fade_hbond_ ){
				do_fade = fade_hbonding_group_score( grp1, grp2, group_score, dw_dE );
				if( dw_dE > 0.0 ){
					for ( Size kk = 1; kk <= n1; ++kk ){ v1[kk] *= dw_dE; }
					for ( Size kk = 1; kk <= n2; ++kk ){ v2[kk] *= dw_dE; }
				}

				if( !do_fade ){ // long enough
					// E = E
					for ( Size kk = 1; kk <= n1; ++kk ){
						core::Size const &atm1( grp1.atms[kk] );
						Vector const & atom1xyz( rsd1.xyz( atm1 ) );
						r1_atom_derivs[ atm1 ].f2() += v1[kk];
						r1_atom_derivs[ atm1 ].f1() += atom1xyz.cross( -v1[kk] );
					}
					for ( Size ll = 1; ll <= n2; ++ll ){
						core::Size const &atm2( grp2.atms[ll] );
						Vector const & atom2xyz( rsd2.xyz( atm2 ) );
						r2_atom_derivs[ atm2 ].f2() += v2[ll];
						r2_atom_derivs[ atm2 ].f1() += atom2xyz.cross( -v2[ll] );
					}
				}
			}

			// 2. derivative on truncation part: E*dsw
			// below will matter if applying different weight on bb/sc...

			//if( !use_subtract ){
			if( use_shift_ ){
				if( ncom1 > 0 ){
					core::Real const c_grp1_heavy = 1.0/((core::Real)(ncom1));

					for ( Size kk = 1; kk <= ncom1; ++kk ){
						core::Size const &atm1( grp1.comatms[kk] );
						Vector const & atom1xyz( rsd1.xyz( atm1 ) );
						Vector f2 = c_grp1_heavy*group_score*dsw_dr*dcom*elec_weight;
						//f2r1_[ atm1 ] += f2;
						r1_atom_derivs[ atm1 ].f2() += f2;
						r1_atom_derivs[ atm1 ].f1() += atom1xyz.cross( -f2 );
					}
				}

				if( ncom2 > 0 ){
					core::Real const c_grp2_heavy = 1.0/((core::Real)(ncom2));

					for ( Size kk = 1; kk <= ncom2; ++kk ){
						core::Size const &atm2( grp2.comatms[kk] );
						Vector const & atom2xyz( rsd2.xyz( atm2 ) );
						Vector f2 = -c_grp2_heavy*group_score*dsw_dr*dcom*elec_weight;
						//f2r2_[ atm2 ] += f2;
						r2_atom_derivs[ atm2 ].f2() += f2;
						r2_atom_derivs[ atm2 ].f1() += atom2xyz.cross( -f2 );
					}
				}
			}

		} // grp2
	} // grp2
  //TR << "end deriv" << std::endl;

	// finally get f1
	//for( Size ii = 1; ii <= rsd1.natoms(); ++ii ){
	//	Vector const & atom1xyz( rsd1.xyz( ii ) );
	//	Vector const & f2 =	f2r1_[ii];
	//		r1_atom_derivs[ ii ].f2() += f2;
	//	r1_atom_derivs[ ii ].f1() += atom1xyz.cross( -f2 );
	//}
	//for( Size ii = 1; ii <= rsd2.natoms(); ++ii ){
	//	Vector const & atom1xyz( rsd2.xyz( ii ) );
	//	Vector const & f2 =	f2r2_[ii];
	//	r2_atom_derivs[ ii ].f2() += f2;
	//		r2_atom_derivs[ ii ].f1() += atom1xyz.cross( -f2 );
	//}

	Erespair *= elec_weight;

}

} // namespace elec
} // namespace scoring
} // namespace core

#endif
