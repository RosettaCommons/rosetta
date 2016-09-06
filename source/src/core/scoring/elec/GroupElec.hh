// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

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

	Real
	eval_respair_group_coulomb(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
		//hbonds::HBondSet const & hbond_set
	) const;

	void
	eval_respair_group_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		//hbonds::HBondSet const & hbond_set,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs,
		core::Real const elec_weight,
		core::Real & Erespair
	) const;

private:

	void
	build_groupinfo( std::string const & group_file,
		bool const extra = false );

	ResElecGroup const &
	get_group( core::chemical::ResidueType const &rsdtype ) const;


	Vector
	get_grpdis2( conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< Size > const &com1atms,
		utility::vector1< Size > const &com2atms,
		core::Vector &com1,
		core::Vector &com2
	) const ;

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
	//std::string grp_cpfxn_mode_;

	/// @delete in its current form GroupElec is not assignable due to presense of std::map< const std::string, ...>
	/// but compiler tries to generate assigment operator anyway
	GroupElec & operator= ( const GroupElec & ) = delete;

}; // class GroupElec

inline
Real
GroupElec::eval_grp_trunc( bool const &use_switch,
	Real const &grpdis2,
	bool const &eval_deriv,
	Real &dsw_dr
) const
{
	dsw_dr = 0.0;
	Real const &max_dis2( coulomb().max_dis2() );
	Real const &max_dis( coulomb().max_dis() );

	if ( grpdis2 > max_dis2 ) return 0.0;

	core::Real const grpdis( std::sqrt(grpdis2) );

	core::Real const dis_sw( fade_param1_ );
	core::Real const min_sw( max_dis - dis_sw );

	// switch function stuffs
	core::Real const sw_dis2 = (max_dis - dis_sw)*(max_dis - dis_sw);
	core::Real const R3on_off( max_dis2 - 3.0*sw_dis2 );
	core::Real Ron_off( 1.0/(max_dis2 - sw_dis2) );
	core::Real const Ron_off_3( Ron_off*Ron_off*Ron_off );
	core::Real dr1(max_dis2 - grpdis2);
	core::Real dr2(R3on_off + 2.0*max_dis2);

	// shift function
	core::Size const nexp = (core::Size)( fade_param2_ );
	Real arg = 1.0;
	for ( core::Size iexp = 1; iexp <= nexp; ++iexp ) {
		arg *= (grpdis-min_sw)/dis_sw * (grpdis-min_sw)/dis_sw;
	}

	Real const sf1 = 1.0 - arg;
	Real const darg = 2.0*nexp*arg/((grpdis-min_sw)*grpdis);

	// fade function
	core::Real sw( 1.0 );
	if ( grpdis2 > sw_dis2 ) {
		if ( use_switch ) {
			sw = dr1*dr1*dr2*Ron_off_3;
		} else {
			sw = sf1*sf1;
		}
	}

	if ( eval_deriv ) {
		// fade function
		if ( grpdis2 > sw_dis2 ) {
			if ( use_switch ) {
				dsw_dr = 4.0*dr1*Ron_off_3*(dr1-dr2);
			} else {
				dsw_dr = -2.0*darg*sf1;
			}
		}
	}

	return sw;
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
	if ( dis2 < 1.0 ) {
		fdie = die;
		ddie_dr = 0.0;
	} else if ( !coulomb().no_dis_dep_die() ) {  // ddd
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
	if ( dis2 < 1.0 ) {
		e = -d+2.0;
	} else { //regular
		e = 1.0/d;
	}

	e *= q1*q2*C0;

	if ( eval_deriv ) {
		core::Real de_dr( 0.0 ); // on Coulomb part without fdie
		if ( dis2 < 1.0 ) { //use linear softening at short-distance
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

} // namespace elec
} // namespace scoring
} // namespace core

#endif
