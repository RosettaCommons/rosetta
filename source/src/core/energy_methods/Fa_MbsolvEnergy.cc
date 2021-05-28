// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/Fa_MbsolvEnergy.hh
/// @author Patrick Barth


// Unit headers
#include <core/energy_methods/Fa_MbsolvEnergy.hh>
#include <core/energy_methods/Fa_MbsolvEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh> //pba
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh> //pba
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <ObjexxFCL/formatted.o.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh> //pba


#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>


namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the Fa_MbsolvEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
Fa_MbsolvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< Fa_MbsolvEnergy >(
		*( core::scoring::ScoringManager::get_instance()->etable( options ).lock() ),
		*( core::scoring::ScoringManager::get_instance()->memb_etable( options.etable_type() ).lock() ),
		options.analytic_membetable_evaluation()
	);
}

core::scoring::ScoreTypes
Fa_MbsolvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( fa_mbsolv );
	return sts;
}

Fa_MbsolvEnergy::Fa_MbsolvEnergy( core::scoring::etable::Etable const & etable_in, core::scoring::etable::MembEtable const & memb_etable_in, bool const analytic_membetable_evaluation):
	parent( utility::pointer::make_shared< Fa_MbsolvEnergyCreator >() ),
	etable_(etable_in),
	memb_etable_(memb_etable_in),
	solv1_(memb_etable_in.solv1()),
	solv2_(memb_etable_in.solv2()),
	dsolv1_( memb_etable_in.dsolv1() ),
	dsolv2_( memb_etable_in.dsolv2() ),
	dsolv_( etable_in.dsolv() ),
	memb_solv1_(memb_etable_in.memb_solv1()),
	memb_solv2_(memb_etable_in.memb_solv2()),
	memb_dsolv1_( memb_etable_in.memb_dsolv1() ),
	memb_dsolv2_( memb_etable_in.memb_dsolv2() ),
	//added to calculate fampsolv analytically
	lk_dgfree_( memb_etable_in.lk_dgfree() ),
	memb_lk_dgfree_( memb_etable_in.memb_lk_dgfree() ),
	lj_radius_( memb_etable_in.lj_radius() ),
	lk_volume_( memb_etable_in.lk_volume() ),
	lk_lambda_( memb_etable_in.lk_lambda() ),
	//
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.get_bins_per_A2()),
	verbose_( false ),
	max_dis_(etable_in.max_dis()),
	max_normal_dis_( max_dis_ - 1.5 ),
	potential_( core::scoring::ScoringManager::get_instance()->get_Membrane_FAPotential() ),
	analytic_etable_evaluation_(analytic_membetable_evaluation)
{}

Distance
Fa_MbsolvEnergy::atomic_interaction_cutoff() const
{
	return etable_.max_dis();
}

/// clone
core::scoring::methods::EnergyMethodOP
Fa_MbsolvEnergy::clone() const
{
	return utility::pointer::make_shared< Fa_MbsolvEnergy >( *this );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
///
void
Fa_MbsolvEnergy::setup_for_scoring(
	pose::Pose & pose, core::scoring::ScoreFunction const &
) const
{
	potential_.compute_fa_projection( pose );
}


void
Fa_MbsolvEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	Real fa_mbsolv_score( 0.0 );

	get_residue_pair_energy( rsd1, rsd2, pose, fa_mbsolv_score);

	// hydrate/SPaDES protocol
	// track if we should ignore fa_sol
	bool no_fa_sol = false;

	// hydrate/SPaDES protocol
	// check to make sure water_hybrid_sf and ignore_fa_sol_at_positions is on
	if ( basic::options::option[ basic::options::OptionKeys::score::water_hybrid_sf ] ) {
		if ( basic::options::option[ basic::options::OptionKeys::hydrate::ignore_fa_sol_at_positions ].user() ) {

			// get the vector of positions to ignore
			utility::vector1< Size > const ignore_positions = basic::options::option[ basic::options::OptionKeys::hydrate::ignore_fa_sol_at_positions ]();

			// check if the current residues are supposed to be ignored
			for ( Size pos = 1; pos <= ignore_positions.size(); pos++ ) {
				Size ignore_fa_sol_pos = ignore_positions[ pos ];
				if ( rsd1.seqpos() == ignore_fa_sol_pos || rsd2.seqpos() == ignore_fa_sol_pos ) {
					no_fa_sol = true;
					break;
				}
			}

			// if current residues are to be ignored, set fa_sol to zero
			if ( no_fa_sol ) {
				fa_mbsolv_score = 0.0;
			}
		}

		// waters should also not have fa_mbsolv
		if ( rsd1.name() == "TP3" || rsd2.name() == "TP3" ) {
			fa_mbsolv_score = 0.0;
		}
	}

	emap[ core::scoring::fa_mbsolv ] += fa_mbsolv_score;
}


////////////////////////////////////////////////
void
Fa_MbsolvEnergy::get_residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & fa_mbsolv_score
) const
{

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	Real temp_score (0.0);

	using namespace core::scoring::etable::count_pair;
	CountPairFunctionOP cpfxn( nullptr );

	if ( same_res ) {
		cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
	} else {
		cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {
		Vector const heavy_atom_i( rsd1.xyz( i ) );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0; Size path_dist( 0 );
			if ( !cpfxn->count( i, j, cp_weight, path_dist ) ) continue;

			Vector const heavy_atom_j( rsd2.xyz( j ) );

			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Real dummy_deriv( 0.0 );

			//pbadebug WARNING
			bool debug(false);

			temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), d2, dummy_deriv,
				core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd1.seqpos(),i),
				core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd2.seqpos(),j),debug);

			if ( same_res ) temp_score *= 0.5;
			fa_mbsolv_score += temp_score;

			/*if ( verbose_ && std::abs( fa_mbsolv_score ) > 0.1 )
			std::cout << "fa_mbsolv_score: rsd1 " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " rsd2 " << rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,3,fa_mbsolv_score) << ' ' << std::
			endl;*/
		} // j
	} // i
}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
Fa_MbsolvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & scfxn
) const
{
	potential_.compute_fa_projection( pose );
	pose.update_residue_neighbors();
	fa_mbsolv_weight_ = scfxn.weights()[ core::scoring::fa_mbsolv ];
}


////////////////////////////////////////////////
Real
Fa_MbsolvEnergy::eval_lk(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const & d2,
	Real & deriv,
	Real const & f1,
	Real const & f2,
	bool & debug ) const
{

	Real temp_score( 0.0 );
	deriv = 0.0;
	//Make this an input option for efficiency
	bool const eval_deriv( true );

	if ( !analytic_etable_evaluation_ ) {

		if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) return 0.0;

		Real const d2_bin = d2 * get_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real frac = d2_bin - ( disbin - 1 );

		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin  ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		int const l1 = solv1_.index( disbin, atom2.type(), atom1.type() );
		int const l2 = l1 + 1;

		//pba Membrane specific solvation
		//pba solvation of atom1 based on its distance from the membrane center on the membrane normal

		Real e11 = f1 * solv1_[ l1 ] + (1 - f1) * memb_solv1_[ l1 ];
		Real e12 = f1 * solv1_[ l2 ] + (1 - f1) * memb_solv1_[ l2 ];

		//pba solvation of atom2 based on its distance from the membrane center on the membrane normal

		Real e21 = f2 * solv2_[ l1 ] + (1 - f2) * memb_solv2_[ l1 ];
		Real e22 = f2 * solv2_[ l2 ] + (1 - f2) * memb_solv2_[ l2 ];

		Real e1 = e11 + e21;
		Real e2 = e12 + e22;

		temp_score = e1 + frac * ( e2 - e1 ); //temp_score = weight * ( e1 + frac * ( e2 - e1 ) );

		if ( debug ) {
			std::cout << "f1 s1l1 mbs1l1 s1l2 mbs1l2 " << f1 << " " << solv1_[ l1 ] << " " << memb_solv1_[ l1 ] << " " <<
				solv1_[ l2 ] << " " << memb_solv1_[ l2 ] << std::endl;
			std::cout << "f2 s2l1 mbs2l1 s2l2 mbs2l2 " << f2 << " " << solv2_[ l1 ] << " " << memb_solv2_[ l1 ] << " " <<
				solv2_[ l2 ] << " " << memb_solv2_[ l2 ] << std::endl;
		}

		if ( eval_deriv ) {
			//   int const l1 = dsolv1_.index( disbin, atom2.type(), atom1.type() ),
			//    l2 = l1 + 1;
			//   Real e1 = dsolv1_[ l1 ];
			//   deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
			e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
			e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
			e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
			e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
			e1 = e11 + e21;
			e2 = e12 + e22;

			deriv = e1 + frac * ( e2 - e1 );
			deriv = deriv / std::sqrt( d2 );
		}

	} else {

		Real solve1(0.0), solve2(0.0), solvE1( 0.0 ), solvE2( 0.0 ), membsolvE1( 0.0 ), membsolvE2( 0.0 );
		solvationE( atom1, atom2, d2, solvE1, solvE2, membsolvE1, membsolvE2);

		solve1 = f1 * solvE1 + (1 - f1) * membsolvE1;
		solve2 = f2 * solvE2 + (1 - f2) * membsolvE2;
		temp_score = solve1 + solve2;
	}

	return temp_score;
}
////////////////////////////////////////////////////////////////////////////////////
Real
Fa_MbsolvEnergy::eval_dE_dR_over_r(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	core::scoring::EnergyMap const & /*weights*/,
	Vector & F1,
	Vector & F2,
	Real const & f1,
	Real const & f2
) const
{
	Real deriv;
	F1 = atom1.xyz().cross( atom2.xyz() );
	F2 = atom1.xyz() - atom2.xyz();
	Real d2,frac;

	d2 = atom1.xyz().distance_squared( atom2.xyz() );

	if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real(0.0) ) ) return 0.0;

	if ( !analytic_etable_evaluation_ ) {

		// bin by distance:
		Real const d2_bin = d2 * get_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		frac = d2_bin - ( disbin - 1 );

		int const l1 = dsolv1_.index( disbin, atom1.type(), atom2.type()),
			l2 = l1 + 1;

		Real e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
		Real e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
		Real e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
		Real e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
		Real e1 = e11 + e21;
		Real e2 = e12 + e22;

		//Real deriv = fa_mbsolv_weight_ * ( e1 + frac * ( e2 - e1 ) );
		deriv =( e1 + frac * ( e2 - e1 ) );

	} else {
		Real dsolve1( 0.0 ), dsolve2( 0.0 ), dmembsolve1( 0.0 ), dmembsolve2( 0.0 );
		dsolvationE( atom1, atom2, d2, dsolve1, dsolve2, dmembsolve1, dmembsolve2 );
		Real de1 = f1 * dsolve1 + (1 - f1) * dmembsolve1;
		Real de2 = f2 * dsolve2 + (1 - f2) * dmembsolve2;
		deriv = de1 + de2;
	}

	return deriv / std::sqrt( d2 );
}

//////////////////////////////////////////////////////////////////////////////////////
//solvation of atom i and j in water and chex
void
Fa_MbsolvEnergy::solvationE(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	core::Real dis2,
	core::Real &solvE1,
	core::Real &solvE2,
	core::Real &membsolvE1,
	core::Real &membsolvE2
) const {
	core::Real solv1, solv2;

	core::Real min_dis2 = etable_.min_dis2();
	core::Real max_dis2 = etable_.max_dis2();
	core::Real epsilon = etable_.epsilon();


	if ( dis2 > max_dis2 + epsilon ) {
		return;
	}
	if ( dis2 < min_dis2 ) {
		dis2 = min_dis2;
	}

	core::Real dis = std::sqrt(dis2); //distance between atom i and j


	//pulling sigma from MembEtable
	Real lk_min_dis2sigma = etable_.lk_min_dis2sigma();
	Real sigma = memb_etable_.lj_sigma(atom1.type(), atom2.type());

	Real thresh_dis = lk_min_dis2sigma * sigma;

	Real thresh_dis_min = thresh_dis/2.0;


	if ( dis > max_normal_dis_ && dis < max_dis_ ) {

		core::Real start_damping = max_normal_dis_;
		core::Real t = (dis - start_damping) / (max_dis_ - start_damping);
		core::Real pk1 = solv( atom1.type(), atom2.type(), start_damping );
		core::Real mk1 = solv_deriv( atom1, start_damping ) * pk1;

		core::Real pk2 = solv( atom2.type(), atom1.type(), start_damping );
		core::Real mk2 = solv_deriv( atom2, start_damping ) * pk2;

		core::Real a1 = mk1*(max_dis_ - start_damping) - ( -pk1 );
		core::Real b1 = -pk1;

		core::Real a2 = mk2*(max_dis_ - start_damping) - ( -pk2 );
		core::Real b2 = -pk2;

		solv1 = (1-t)*pk1 + t*(1-t)*((1-t)*a1 + t*b1);
		solv2 = (1-t)*pk2 + t*(1-t)*((1-t)*a2 + t*b2);

	} else if ( dis < thresh_dis && dis > thresh_dis_min ) {

		core::Real t = ( dis - thresh_dis_min ) / ( thresh_dis - thresh_dis_min );
		//pk11 is the value of solvation at thresh_dis_min for atom 1
		//pk12 is the value of solvation at thresh_dis for atom 1
		//mk12 is the value of the derivative at thresh_dis for atom 1
		//the value of the slope at thresh_dis_min is 0
		//pk21, pk22, mk22 are the same except it is for atom2

		core::Real pk11 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		core::Real pk12 = solv( atom1.type(), atom2.type(), thresh_dis );
		core::Real mk12 = solv_deriv( atom1, thresh_dis ) * pk12;

		core::Real pk21 = solv( atom2.type(), atom1.type(), thresh_dis_min );
		core::Real pk22 = solv( atom2.type(), atom1.type(), thresh_dis );
		core::Real mk22 = solv_deriv( atom2, thresh_dis ) * pk22;

		solv1 = ( 2*std::pow(t,3) - 3*std::pow(t,2) + 1 )*pk11 + (-2*std::pow(t,3) + 3*std::pow(t,2))*pk12 + ( std::pow(t,3) - std::pow(t,2) )*mk12;
		solv2 = ( 2*std::pow(t,3) - 3*std::pow(t,2) + 1 )*pk21 + (-2*std::pow(t,3) + 3*std::pow(t,2))*pk22 + ( std::pow(t,3) - std::pow(t,2) )*mk22;

	} else if ( dis <= thresh_dis_min ) {

		solv1 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		solv2 = solv( atom2.type(), atom1.type(), thresh_dis_min );

	} else if ( dis >= max_dis_ ) {
		solv1 = 0;
		solv2 = 0;
	} else {
		solv1 = solv( atom1.type(), atom2.type(), dis );
		solv2 = solv( atom2.type(), atom1.type(), dis );
	}

	solvE1 = lk_dgfree_[ atom1.type() ] * solv1;
	solvE2 = lk_dgfree_[ atom2.type() ] * solv2;
	membsolvE1 = memb_lk_dgfree_[ atom1.type() ] * solv1;
	membsolvE2 = memb_lk_dgfree_[ atom2.type() ] * solv2;

}

core::Real
Fa_MbsolvEnergy::solv(
	int atom1type,
	int atom2type,
	core::Real dis
) const {
	Real const k = { -0.089793561062583294 }; //inv_neg2_tms_pi_sqrt_pi
	return k * lk_volume_[ atom2type ] * solv_piece( atom1type, dis );
}

core::Real
Fa_MbsolvEnergy::solv_piece(
	int atom_type,
	core::Real d
) const {
	core::Real lambda = lk_lambda_[ atom_type ];
	core::Real denom = d * d * lambda;
	core::Real dis_rad = d - lj_radius_[ atom_type ];
	core::Real lk_inv_lambda2 = ( 1.0/lambda ) * ( 1.0/lambda );
	core::Real expo = ( dis_rad * dis_rad ) * lk_inv_lambda2;
	core::Real solv = ( 1 / denom ) * std::exp(-expo);
	return solv;
}

//solvation partial derivative wrt distance from atom i and j
void
Fa_MbsolvEnergy::dsolvationE(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	core::Real dis2,
	core::Real &dsolvE1,
	core::Real &dsolvE2,
	core::Real &dmembsolvE1,
	core::Real &dmembsolvE2
) const {

	core::Real min_dis2 = etable_.min_dis2();
	core::Real max_dis2 = etable_.max_dis2();
	core::Real epsilon = etable_.epsilon();


	if ( dis2 > max_dis2 + epsilon ) {
		return;
	}
	if ( dis2 < min_dis2 ) {
		dis2 = min_dis2;
	}

	core::Real solvE1, solvE2, membsolvE1, membsolvE2, deriv1, deriv2;
	core::Real dis = std::sqrt(dis2);

	//see comments above in SolvationE for notes on sigma and thresh_dis
	core::Real lk_min_dis2sigma = etable_.lk_min_dis2sigma();
	core::Real sigma = memb_etable_.lj_sigma(atom1.type(), atom2.type());

	core::Real thresh_dis = lk_min_dis2sigma * sigma;
	core::Real thresh_dis_min = thresh_dis/2.0;

	if ( dis > max_normal_dis_ && dis < max_dis_ ) {

		core::Real start_damping = max_normal_dis_;
		core::Real t = (dis - start_damping) / (max_dis_ - start_damping);

		core::Real pk1 = solv( atom1.type(), atom2.type(), start_damping );
		core::Real mk1 = solv_deriv( atom1, start_damping ) * pk1;

		core::Real pk2 = solv( atom2.type(), atom1.type(), start_damping );
		core::Real mk2 = solv_deriv( atom2, start_damping ) * pk2;

		core::Real a1 = mk1*(max_dis_ - start_damping) - ( -pk1 );
		core::Real b1 = -pk1;

		core::Real a2 = mk2*(max_dis_ - start_damping) - ( -pk2 );
		core::Real b2 = -pk2;

		deriv1 = (-pk1 + (1 - 2*t)*((1-t)*a1 + t*b1) + (t - t*t)*(-a1 + b1))*(1/(max_dis_ - start_damping));
		deriv2 = (-pk2 + (1 - 2*t)*((1-t)*a2 + t*b2) + (t - t*t)*(-a2 + b2))*(1/(max_dis_ - start_damping));

		dsolvE1 = lk_dgfree_[ atom1.type() ] * deriv1;
		dsolvE2 = lk_dgfree_[ atom2.type() ] * deriv2;
		dmembsolvE1 = memb_lk_dgfree_[ atom1.type() ] * deriv1;
		dmembsolvE2 = memb_lk_dgfree_[ atom2.type() ] * deriv2;

	} else if ( dis < thresh_dis && dis > thresh_dis_min ) {

		core::Real t = ( dis - thresh_dis_min ) / ( thresh_dis - thresh_dis_min );

		core::Real pk11 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		core::Real pk12 = solv( atom1.type(), atom2.type(), thresh_dis );
		core::Real mk12 = solv_deriv( atom1, thresh_dis ) * pk12;

		core::Real pk21 = solv( atom2.type(), atom1.type(), thresh_dis_min );
		core::Real pk22 = solv( atom2.type(), atom1.type(), thresh_dis );
		core::Real mk22 = solv_deriv( atom2, thresh_dis ) * pk22;

		deriv1 = (( 6*std::pow(t,2) - 6*t )*pk11 + ( -6*std::pow(t,2) + 6*t )*pk12 + ( 3*std::pow(t,2) - 2*t )*mk12 )*(1/(thresh_dis-thresh_dis_min));
		deriv2 = (( 6*std::pow(t,2) - 6*t )*pk21 + ( -6*std::pow(t,2) + 6*t )*pk22 + ( 3*std::pow(t,2) - 2*t )*mk22 )*(1/(thresh_dis-thresh_dis_min));

		dsolvE1 = lk_dgfree_[ atom1.type() ] * deriv1;
		dsolvE2 = lk_dgfree_[ atom2.type() ] * deriv2;
		dmembsolvE1 = memb_lk_dgfree_[ atom1.type() ] * deriv1;
		dmembsolvE2 = memb_lk_dgfree_[ atom2.type() ] * deriv2;

	} else if ( ( dis <= thresh_dis_min || dis >= max_dis_) ) {

		dsolvE1 = 0.0;
		dsolvE1 = 0.0;
		dmembsolvE1 = 0.0;
		dmembsolvE1 = 0.0;

	} else {
		solvationE( atom1, atom2, dis2, solvE1, solvE2, membsolvE1, membsolvE2 );

		deriv1 = solv_deriv( atom1, dis );
		deriv2 = solv_deriv( atom2, dis );

		dsolvE1 = solvE1 * deriv1;
		dsolvE2 = solvE2 * deriv2;
		dmembsolvE1 = membsolvE1 * deriv1;
		dmembsolvE2 = membsolvE2 * deriv2;

	}


}

core::Real
Fa_MbsolvEnergy::solv_deriv(
	conformation::Atom const & atom,
	core::Real dis
) const {
	core::Real deriv = -2.0 * (((dis - lj_radius_[ atom.type() ]) * ( 1/std::pow(lk_lambda_[ atom.type() ], 2) )) + (1/dis));
	return deriv;
}


//////////////////////////////////////////////////////////////////////////////////////
void
Fa_MbsolvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	core::scoring::ScoreFunction const &,// sfxn,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	if ( m > rsd1.nheavyatoms() ) return;

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	bool const pos1_fixed( domain_map( i ) != 0 );

	// cached energies object
	core::scoring::Energies const & energies( pose.energies() );

	// the neighbor/energy links
	core::scoring::EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( utility::graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );
		bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

		using namespace core::scoring::etable::count_pair;
		CountPairFunctionOP cpfxn( nullptr );

		if ( same_res ) {
			cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
		} else {
			cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
		}

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0; Size path_dist(0);

			if ( !cpfxn->count(m, n, cp_weight, path_dist ) ) continue;

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();
			//Vector const d_ij_norm = d_ij.normalized();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Vector f1( 0.0 ), f2( 0.0 );
			Real const dE_dR_over_r
				( eval_dE_dR_over_r( rsd1.atom(m), rsd2.atom(n), weights, f1, f2,
				core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd1.seqpos(),m),
				core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd2.seqpos(),n) ) );


			//now calculate f1 and f2 with respect to distance from membrane center
			//need to end up with F1 = f1d*dE/dd + f1z*dE/dz
			//                    F2 = f2d*dE/dd + f2z*dE/dz
			Real solvE1( 0.0 ), solvE2( 0.0 ), membsolvE1( 0.0 ), membsolvE2( 0.0 );
			solvationE( rsd1.atom(m), rsd2.atom(n), d2, solvE1, solvE2, membsolvE1, membsolvE2);

			//derivative with respect to distance from membrane center
			Real dE_dZ = Membrane_FAEmbed_from_pose( pose ).fa_proj_deriv(rsd1.seqpos(),m) * (solvE1 - membsolvE1 );

			Vector f1z( 0.0 ), f2z( 0.0 );

			//Vector const d_iz = Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m) - heavy_atom_i;
			Vector const d_iz = heavy_atom_i - Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m);
			Real const d_iz_norm = d_iz.length();
			if ( d_iz_norm != Real(0.0) ) {
				Real const invd = 1.0 / d_iz_norm;
				f2z = d_iz * invd;
				f1z = heavy_atom_i.cross(Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m));
				f1z *= invd;
			}

			if ( dE_dR_over_r == 0.0 && dE_dZ == 0.0 ) continue;

			if ( same_res ) {
				F1 += 0.5 * cp_weight * fa_mbsolv_weight_ * ( ( dE_dR_over_r * f1) + ( dE_dZ * f1z ) );
				F2 += 0.5 * cp_weight * fa_mbsolv_weight_ * ( ( dE_dR_over_r * f2) + ( dE_dZ * f2z ) );
			} else {
				F1 += cp_weight * fa_mbsolv_weight_ * ( ( dE_dR_over_r * f1) + ( dE_dZ * f1z ) );
				F2 += cp_weight * fa_mbsolv_weight_ * ( ( dE_dR_over_r * f2) + ( dE_dZ * f2z ) );
			}
		}
	}
}


////////////////////////////////////////////////
void
Fa_MbsolvEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}

////////////////////////////////////////////////
void
Fa_MbsolvEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{

	Real fa_mbsolv_score( 0.0 );

	get_residue_pair_energy( rsd, rsd, pose, fa_mbsolv_score);

	// hydrate/SPaDES protocol
	// track if we should ignore fa_sol
	bool no_fa_sol = false;

	// hydrate/SPaDES protocol
	// check to make sure water_hybrid_sf and ignore_fa_sol_at_positions is on
	if ( basic::options::option[ basic::options::OptionKeys::score::water_hybrid_sf ] ) {
		if ( basic::options::option[ basic::options::OptionKeys::hydrate::ignore_fa_sol_at_positions ].user() ) {

			// get the vector of positions to ignore
			utility::vector1< Size > const ignore_positions = basic::options::option[ basic::options::OptionKeys::hydrate::ignore_fa_sol_at_positions ]();

			// check if the current residues are supposed to be ignored
			for ( Size pos = 1; pos <= ignore_positions.size(); pos++ ) {
				Size ignore_fa_sol_pos = ignore_positions[ pos ];
				if ( rsd.seqpos() == ignore_fa_sol_pos ) {
					no_fa_sol = true;
					break;
				}
			}

			// if current residues are to be ignored, set fa_sol to zero
			if ( no_fa_sol ) {
				fa_mbsolv_score = 0.0;
			}
		}

		// waters should also not have fa_mbsolv
		if ( rsd.name() == "TP3" ) {
			fa_mbsolv_score = 0.0;
		}
	}

	emap[ core::scoring::fa_mbsolv ] += fa_mbsolv_score;

}

////////////////////////////////////////////////
void
Fa_MbsolvEnergy::finalize_total_energy(
	pose::Pose & /*pose*/,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & /*emap*/
) const
{
	if ( verbose_ ) std::cout << "DONE SCORING" << std::endl;
}


/// @details Pose must already contain a cenlist object or this method will fail.
core::scoring::Membrane_FAEmbed const &
Fa_MbsolvEnergy::Membrane_FAEmbed_from_pose( pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED;
	return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
}

core::scoring::MembraneTopology const &
Fa_MbsolvEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
core::Size
Fa_MbsolvEnergy::version() const
{
	return 1; // Initial versioning
}


}
}


