// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_hack.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/LK_hack.hh>
#include <core/scoring/methods/LK_hackCreator.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/constraints/AngleConstraint.hh>

#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the LK_hack class,
/// never an instance already in use
methods::EnergyMethodOP
LK_hackCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	etable::EtableCOP etable( ScoringManager::get_instance()->etable( options ) );
	return methods::EnergyMethodOP( new LK_hack( *etable ) );
}

ScoreTypes
LK_hackCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( lk_hack );
	return sts;
}


using namespace constraints;

class LK_SigmoidalFunc : public func::Func {
public:
	LK_SigmoidalFunc();

	func::FuncOP
	clone() const;

	virtual Real func( Real const x ) const;
	virtual Real dfunc( Real const x ) const;

	static Real const ANGLE_CUTOFF_HIGH;
	static Real const ANGLE_CUTOFF_LOW;

	static Real const cos_flipped_ANGLE_CUTOFF_HIGH;
	static Real const cos_flipped_ANGLE_CUTOFF_LOW;

};

using namespace numeric;
using namespace numeric::constants::d;

/// Ramp the score from 1 to 0 over the range from 100 degrees to 90 degrees
Real const LK_SigmoidalFunc::ANGLE_CUTOFF_HIGH( 100.0 * degrees_to_radians );
Real const LK_SigmoidalFunc::ANGLE_CUTOFF_LOW(   90.0 * degrees_to_radians );

Real const LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_HIGH( std::cos( pi - LK_SigmoidalFunc::ANGLE_CUTOFF_HIGH ));
Real const LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_LOW(  std::cos( pi - LK_SigmoidalFunc::ANGLE_CUTOFF_LOW  ));


LK_SigmoidalFunc::LK_SigmoidalFunc() {}

core::scoring::func::FuncOP LK_SigmoidalFunc::clone() const { return core::scoring::func::FuncOP( new LK_SigmoidalFunc ); }

/// @brief a Sigmoidal function that ramps from 1 to 0 over a certain range.
/// Thanks to Mike Tyka for having a sigmoidal function on the top of his head.
Real
LK_SigmoidalFunc::func( Real const x ) const
{
	//std::cout << "LK_SigmoidalFunc: x=" << x << " degx = " << radians_to_degrees * x  << " high " << ANGLE_CUTOFF_HIGH << " low " << ANGLE_CUTOFF_LOW << std::endl;
	/// x in radians, angle between base atom, atom, and desolvator atom
	if ( x >= ANGLE_CUTOFF_HIGH ) { return 1.0; }
	else if ( x <= ANGLE_CUTOFF_LOW ) { return 0.0; }
	else {
		static Real const range = ANGLE_CUTOFF_HIGH - ANGLE_CUTOFF_LOW;
		// x on a zero-to-one range = x'
		Real const xprime = ( ANGLE_CUTOFF_HIGH - x ) / range;
		//std::cout << "LK_SigmoidalFunc::func x= " << radians_to_degrees * x << " xprime= " << xprime << " " << ( 1 - xprime*xprime )*( 1 - xprime*xprime ) << std::endl;
		return ( 1 - xprime*xprime )*( 1 - xprime*xprime );
	}
}

Real
LK_SigmoidalFunc::dfunc( Real const x ) const
{
	/// x in radians, angle between base atom, atom, and desolvator atom
	if ( x >= ANGLE_CUTOFF_HIGH ) { return 0.0; }
	else if ( x <= ANGLE_CUTOFF_LOW ) { return 0.0; }
	else {
		// x on a zero-to-one range
		static Real const range = ANGLE_CUTOFF_HIGH - ANGLE_CUTOFF_LOW;
		Real const xprime = ( ANGLE_CUTOFF_HIGH - x ) / range;
		//Real const eps = 0.00001;
		//Real const numD = ( func( x + eps ) - func( x ) ) / eps ;
		//Real const anD  = ( 4.0 / (range) ) * ( 1 - xprime*xprime ) * ( xprime );
		//std::cout << "LK_SigmoidalFunc::dfunc x= " << radians_to_degrees * x << " xprime= " << xprime << " d: " <<  anD << " range: " << range << " invrange " << 1 / range << std::endl;
		//std::cout << "func( " << x << ") " << func( x ) << " func( " << x + eps << " )= " <<  func( x + eps ) << " numD : " << numD << " numD rat: "<< numD / anD << std::endl;
		return  ( 4.0 / (range) ) * ( 1 - xprime*xprime ) * ( xprime );
	}
}

LK_hack::LK_hack( etable::Etable const & etable_in ) :
	parent( EnergyMethodCreatorOP( new LK_hackCreator ) ),
	etable_(etable_in),
	solv1_(etable_in.solv1()),
	solv2_(etable_in.solv2()),
	dsolv1_( etable_in.dsolv1() ),
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.get_bins_per_A2())
{}

Distance
LK_hack::atomic_interaction_cutoff() const
{
	return etable_.max_dis();
}

/// clone
EnergyMethodOP
LK_hack::clone() const
{
	return EnergyMethodOP( new LK_hack( *this ) );
}

LK_hack::LK_hack( LK_hack const & src ):
	parent( src ),
	etable_(src.etable_),
	solv1_( src.solv1_ ),
	solv2_( src.solv2_ ),
	dsolv1_( src.dsolv1_ ),
	safe_max_dis2_( src.safe_max_dis2_ ),
	get_bins_per_A2_( src.get_bins_per_A2_   )
{
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
LK_hack::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;

	Real score(0.0);

	LK_SigmoidalFunc lksigmoidalfunc;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	/// allocation and deallocation of these vectors is somewhat slow.
	utility::vector1< Vector > res1_base_vectors( rsd1.nheavyatoms(),Vector( 0.0 ));
	utility::vector1< Vector > res2_base_vectors( rsd2.nheavyatoms(), Vector( 0.0));
	utility::vector1< Size >   res1_heavy_is_polar( rsd1.nheavyatoms(), false );
	utility::vector1< Size >   res2_heavy_is_polar( rsd2.nheavyatoms(), false );

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {
		if (rsd1.atom_type(i).is_acceptor() || rsd1.atom_type(i).is_donor()){
			res1_heavy_is_polar[ i ] = true;
			Size non_H_neib = 0;
			Vector  base_pseudo_atom(0);
			for (Size ii = 1; ii <=rsd1.bonded_neighbor(i).size(); ++ii){
				Size neighbor_id = rsd1.bonded_neighbor(i)[ii];
				if ( !  rsd1.atom_is_hydrogen(neighbor_id)){
					base_pseudo_atom += rsd1.xyz(neighbor_id);
					non_H_neib++;
				}
			}
			if ( rsd1.type().n_residue_connections_for_atom(i) > 0 ) {
				/// CONTEXT DEPENDENCY HERE -- e.g. if c_prev moves, rsd1 needs to be rescored.  Fortunately,
				/// if c_prev moves, the internal "psi" for rsd1 will be updated, and this residue will be rescored.
				for ( Size ii = 1; ii <= rsd1.type().residue_connections_for_atom(i).size(); ++ii ) {
					chemical::ResConnID const ii_conn = rsd1.connect_map( rsd1.type().residue_connections_for_atom(i)[ ii ] );
					Size const neighbor_res_id(  ii_conn.resid() );
					Size const nieghbor_atom_id( pose.residue( ii_conn.resid() ).residue_connection( ii_conn.connid() ).atomno() );
					if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
						base_pseudo_atom += pose.residue( neighbor_res_id ).xyz( nieghbor_atom_id );
						non_H_neib++;
					}
				}
			}
			base_pseudo_atom /= non_H_neib;
			res1_base_vectors[i] =  rsd1.xyz(i) - base_pseudo_atom;
			res1_base_vectors[i].normalize();
		}
	}
	//same for residue 2//
	for ( Size i = 1, i_end = rsd2.nheavyatoms(); i <= i_end; ++i ) {
		if (rsd2.atom_type(i).is_acceptor() || rsd2.atom_type(i).is_donor()){
			res2_heavy_is_polar[ i ] = true;
			Size non_H_neib = 0;
			Vector  base_pseudo_atom(0);
			for (Size ii = 1; ii <=rsd2.bonded_neighbor(i).size(); ++ii){
				Size neighbor_id = rsd2.bonded_neighbor(i)[ii];
				if ( ! rsd2.atom_is_hydrogen(neighbor_id) ){
					base_pseudo_atom += rsd2.xyz(neighbor_id);
					non_H_neib++;
				}
			}
			if ( rsd2.type().n_residue_connections_for_atom(i) > 0 ) {
				for ( Size ii = 1; ii <= rsd2.type().residue_connections_for_atom(i).size(); ++ii ) {
					chemical::ResConnID const ii_conn = rsd2.connect_map( rsd2.type().residue_connections_for_atom(i)[ ii ] );
					Size const neighbor_res_id( ii_conn.resid() );
					Size const nieghbor_atom_id( pose.residue( ii_conn.resid() ).residue_connection( ii_conn.connid() ).atomno() );
					if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
						base_pseudo_atom += pose.residue( neighbor_res_id ).xyz( nieghbor_atom_id );
						non_H_neib++;
					}
				}
			}

			base_pseudo_atom /= non_H_neib;
			res2_base_vectors[i] = rsd2.xyz(i) - base_pseudo_atom;
			res2_base_vectors[i].normalize();

		}
	}

	Real cp_weight=0.;
	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i){
		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {
				/// if ( weight < 0.1 ) continue; // apl --- I don't think this ever happens
				Real d2 =   rsd1.atom(i).xyz().distance_squared( rsd2.atom(j).xyz() );

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				Real const d2_bin = d2 * get_bins_per_A2_;
				int	disbin = static_cast< int >( d2_bin ) + 1;
				Real	frac = d2_bin - ( disbin - 1 );
				int const l1 = solv1_.index( disbin, rsd2.atom(j).type(), rsd1.atom(i).type() ); // atom i being desolvated by atom j
				int const l2 = l1 + 1;

				Vector i_to_j_vec( rsd2.xyz( j ) - rsd1.xyz( i ) );
				bool i_to_j_vec_normalized( false ); // don't normalize twice!

				Real i_to_j_angle_weight( 1.0 ), j_to_i_angle_weight( 1.0 );

				if ( res1_heavy_is_polar[i] ) {
					if (res1_base_vectors[i].dot( i_to_j_vec )  >=  0 ) {
						//std::cout << "Normalizing i_to_j" << std::endl;

						i_to_j_vec *= 1 / ( std::sqrt( d2 ) ); /// we already have d2, reuse it
						i_to_j_vec_normalized = true;
						Real dotprod = res1_base_vectors[i].dot( i_to_j_vec );
						//std::cout << "r1 base dot prod: " << dotprod << " angle " << radians_to_degrees * (pi - std::acos( dotprod )) << " with " << rsd2.seqpos() << " " << j << std::endl;
						if ( dotprod >= LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_HIGH ) {
							// noop, i_to_j weight is already 1.
						} else if ( dotprod <= LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_LOW ) {
							i_to_j_angle_weight = 0;
						} else {
							Real angle = pi - std::acos( dotprod );
							i_to_j_angle_weight = lksigmoidalfunc.func( angle );
							//std::cout << "Angle: " << radians_to_degrees * angle << " func " << i_to_j_angle_weight << std::endl;
						}
					} else { // Dot product is 0 or negative; do not count this interaction
						i_to_j_angle_weight = 0.0;
					}
				} // else, rsd1 is apolar, and its weight should be 1, which it already is

				//std::cout << "Desolvation of atom i: " << rsd1.atom_name( i ) << " on " << rsd1.name() << " with weight: " << i_to_j_angle_weight << std::endl;
				if ( i_to_j_angle_weight != 0 ) {
					score += i_to_j_angle_weight * cp_weight * ( (1.-frac)* solv1_[ l1 ] + frac * solv1_[ l2 ]);
				}

				if ( res2_heavy_is_polar[j] ) {
					Vector j_to_i_vec = -1 * i_to_j_vec;
					if (res2_base_vectors[j].dot( j_to_i_vec )  >  0 ) {
						if ( ! i_to_j_vec_normalized ) {
							//std::cout << "Normalizing j_to_i" << std::endl;
							j_to_i_vec *= 1 / ( std::sqrt( d2 ));
						}
						Real dotprod = res2_base_vectors[j].dot( j_to_i_vec );
						//std::cout << "r2 base dot prod: " << dotprod << " angle " << radians_to_degrees * (pi - std::acos( dotprod ))  << " with " << rsd1.seqpos() << " " << i << std::endl;

						if ( dotprod >= LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_HIGH ) {
							// noop, i_to_j weight is already 1.
						} else if ( dotprod <= LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_LOW ) {
							j_to_i_angle_weight = 0;
						} else {
							Real angle = pi - std::acos( dotprod );
							j_to_i_angle_weight = lksigmoidalfunc.func( angle );
							//std::cout << "Angle: " << radians_to_degrees * angle << " func " << j_to_i_angle_weight << std::endl;

						}
					} else { // Dot product is 0 or negative; do not count this interaction
						j_to_i_angle_weight = 0.0;
					}
				} // else, rsd2 is apolar, and its weight should be 1, which it already is

				//std::cout << "Desolvation of atom j: " << rsd2.atom_name( j ) << " on " << rsd2.name() << " with weight: " << j_to_i_angle_weight << std::endl;
				if ( j_to_i_angle_weight != 0 ) {
					score += j_to_i_angle_weight * cp_weight * ( (1.-frac)* solv2_[ l1 ] + frac * solv2_[ l2 ]);
				}
				//std::cout << "Finished pair " << i << " " << j << std::endl;
			}
		}
	}

	emap[ lk_hack ]+=score;

}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////

void
LK_hack::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const & sfxn
) const
{
	pose.update_residue_neighbors();

	lk_hack_weight_ = sfxn.weights()[ lk_hack ];
	allocate_appropriate_memory( pose );
	calculate_orientation_vectors_and_pseudo_base_atoms( pose );
	calculate_derivatives_for_atoms_and_pseudo_base_atoms( pose );
	distribute_pseudo_base_atom_derivatives( pose );
}

void
LK_hack::allocate_appropriate_memory( pose::Pose const & pose ) const
{
	Size const total_residue( pose.total_residue() );

	nneighbs_.resize( total_residue );
	orientation_vectors_.resize( total_residue );
	base_pseudo_atom_centers_.resize( total_residue );
	atom_f1_f2s_.resize( total_residue );
	base_pseudo_atom_f1_f2s_.resize( total_residue );

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		Size const ii_nheavy( pose.residue( ii ).nheavyatoms() );
		nneighbs_[ ii ].resize( ii_nheavy );
		orientation_vectors_[ ii ].resize( ii_nheavy );
		base_pseudo_atom_centers_[ ii ].resize( ii_nheavy );
		atom_f1_f2s_[ ii ].resize( ii_nheavy );
		base_pseudo_atom_f1_f2s_[ ii ].resize( ii_nheavy );
		for ( Size jj = 1; jj <= ii_nheavy; ++jj ) {
			nneighbs_[ ii ][ jj ] = 0;
			orientation_vectors_[ ii ][ jj ] = Vector( 0.0 );
			base_pseudo_atom_centers_[ ii ][ jj ] = Vector( 0.0 );
			atom_f1_f2s_[ ii ][ jj ].first = Vector(0.0);
			atom_f1_f2s_[ ii ][ jj ].second = Vector(0.0);
			base_pseudo_atom_f1_f2s_[ ii ][ jj ].first = Vector(0.0);
			base_pseudo_atom_f1_f2s_[ ii ][ jj ].second = Vector(0.0);
		}
	}
}

void
LK_hack::calculate_orientation_vectors_and_pseudo_base_atoms( pose::Pose const & pose ) const
{
	using namespace conformation;

	Size const total_residue( pose.total_residue() );

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		Size const ii_nheavy( pose.residue( ii ).nheavyatoms() );
		Residue const & ii_res( pose.residue( ii ));
		for ( Size jj = 1; jj <= ii_nheavy; ++jj ) {
			if (ii_res.atom_type(jj).is_acceptor() || ii_res.atom_type(jj).is_donor()) {
				for (Size kk = 1; kk <= ii_res.bonded_neighbor(jj).size(); ++kk){
					Size neighbor_id = ii_res.bonded_neighbor(jj)[kk];
					if ( ! ii_res.atom_is_hydrogen(neighbor_id)){
						base_pseudo_atom_centers_[ii][jj] += ii_res.xyz(neighbor_id);
						++nneighbs_[ii][jj];
					}
				}
				if ( ii_res.type().n_residue_connections_for_atom(jj) > 0 ) {
					for ( Size kk = 1; kk <= ii_res.type().residue_connections_for_atom(jj).size(); ++kk ) {
						chemical::ResConnID kk_conn = ii_res.connect_map( ii_res.type().residue_connections_for_atom(jj)[ kk ] );
						Size const neighbor_res_id( kk_conn.resid() );
						Size const nieghbor_atom_id( pose.residue( kk_conn.resid() ).residue_connection( kk_conn.connid() ).atomno() );
						if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
							base_pseudo_atom_centers_[ii][jj] += pose.residue( neighbor_res_id ).xyz( nieghbor_atom_id );
							++nneighbs_[ii][jj];
						}
					}
				}
				if ( nneighbs_[ii][jj] != 0 ) {
					base_pseudo_atom_centers_[ii][jj] /= nneighbs_[ii][jj];
					orientation_vectors_[ii][jj] = ii_res.xyz( jj ) - base_pseudo_atom_centers_[ii][jj];
				}
			} //end if
		}//end for jj
	} // end for ii
	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		for ( Size jj = 1, jje = orientation_vectors_[ii].size(); jj <= jje; ++jj ) {
			if ( nneighbs_[ii][jj] != 0 )  {
				orientation_vectors_[ ii ][ jj ].normalize_or_zero();
			}
		}
	}
}


void
LK_hack::calculate_derivatives_for_atoms_and_pseudo_base_atoms( pose::Pose const & pose ) const
{
	using namespace constraints;

	Size const total_residue( pose.total_residue() );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			Size jj = (*iru)->get_second_node_ind();
			calculate_derivatives_for_residue_pair( pose, ii, jj );
		}
	}
}

void
LK_hack::calculate_derivatives_for_residue_pair
(
	pose::Pose const & pose,
	Size lower_res_id,
	Size upper_res_id
) const
{
	using namespace etable::count_pair;

	core::scoring::func::FuncOP lkfunc( new LK_SigmoidalFunc );
	AngleConstraint lk_angle_cst( lkfunc ); //Using the stupid and dangerous version of the AngleConstraint ctor

	conformation::Residue const & lowerres( pose.residue( lower_res_id ) );
	conformation::Residue const & upperres( pose.residue( upper_res_id ) );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( lowerres, upperres, CP_CROSSOVER_4 );

	for ( Size ii = 1, iie = lowerres.nheavyatoms(), jje = upperres.nheavyatoms();
		ii <= iie; ++ii ) {
		for ( Size jj = 1; jj <= jje; ++jj ) {
			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( ! cpfxn->count( ii, jj, cp_weight, path_dist ) ) continue;

			Vector ii_to_jj = upperres.xyz( jj ) - lowerres.xyz( ii );
			DistanceSquared ii_jj_dis2 = ii_to_jj.norm_squared();

			if ( ii_jj_dis2 > safe_max_dis2_ || ii_jj_dis2 == 0.0 ) continue;
			Distance const ii_jj_one_over_d( 1 / (std::sqrt( ii_jj_dis2 )));
			ii_to_jj *= ii_jj_one_over_d; // normalize

			/// 1. Atom jj's desolvation of atom ii.
			if ( nneighbs_[ lower_res_id ][ ii ] != 0 ) { // ii polar
				Real const dot_prod = ii_to_jj.dot(orientation_vectors_[lower_res_id][ ii ]);
				//std::cout << "Derivative for lowerres " << lowerres.atom_name( ii ) << " on " << lowerres.name() << " dotprod: " << dot_prod << " with " << upper_res_id << " " << jj << std::endl;
				if ( dot_prod > 0 && dot_prod > LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_LOW ) {

					Vector f1_ii( 0.0 ), f2_ii( 0.0 ), f1_jj( 0.0 ), f2_jj( 0.0 );

					Real const dE_dR_over_r(
						eval_dE_dR_over_r(
						lowerres.atom(ii), upperres.atom(jj),
						ii_jj_one_over_d, ii_jj_dis2,
						f1_ii, f2_ii, f1_jj, f2_jj ) );

					if ( dot_prod > LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_HIGH ) {
						/// sigmoidal function evaluates to 1 and its derivative is 0
						atom_f1_f2s_[ lower_res_id ][ ii ].first  += dE_dR_over_r * cp_weight * f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += dE_dR_over_r * cp_weight * f2_ii;
						atom_f1_f2s_[ upper_res_id ][ jj ].first  += dE_dR_over_r * cp_weight * f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += dE_dR_over_r * cp_weight * f2_jj;
					} else { // ( dot_prod < cos_flipped_ANGLE_CUTOFF_HIGH && dot_prod > cos_flipped_ANGLE_CUTOFF_LOW ) {
						/// ramping range of sigmoidal function "w"
						/// Energy is w*lk; derivative is w'*lk + w*lk'
						/// lk' already computed

						// lk
						Real const d2_bin = ii_jj_dis2 * get_bins_per_A2_;
						int	disbin = static_cast< int >( d2_bin ) + 1;
						Real	frac = d2_bin - ( disbin - 1 );
						int const l1 = solv1_.index( disbin, upperres.atom(jj).type(), lowerres.atom(ii).type() );
						int const l2 = l1 + 1;
						Real const lk = lk_hack_weight_ * cp_weight * ( (1.-frac)* solv1_[ l1 ] + frac * solv1_[ l2 ]);

						// w'
						Vector w_f1_ii( 0.0 ), w_f2_ii( 0.0 ), w_f1_jj( 0.0 ), w_f2_jj( 0.0 ), w_f1_ii_pbase( 0.0 ), w_f2_ii_pbase( 0.0 );
						lk_angle_cst.p1_deriv(
							base_pseudo_atom_centers_[ lower_res_id ][ ii ], lowerres.xyz( ii ), upperres.xyz( jj ),
							w_f1_ii_pbase, w_f2_ii_pbase );

						lk_angle_cst.p2_deriv(
							base_pseudo_atom_centers_[ lower_res_id ][ ii ], lowerres.xyz( ii ), upperres.xyz( jj ),
							w_f1_ii, w_f2_ii );

						lk_angle_cst.p1_deriv(
							upperres.xyz( jj ), lowerres.xyz( ii ), base_pseudo_atom_centers_[ lower_res_id ][ ii ],
							w_f1_jj, w_f2_jj );

						/// w
						Real const w = lkfunc->func(
							angle_radians(
							base_pseudo_atom_centers_[ lower_res_id ][ ii ],
							lowerres.xyz( ii ),
							upperres.xyz( jj ) ) );

						//// Add everything up
						atom_f1_f2s_[ lower_res_id ][ ii ].first  += w * dE_dR_over_r * cp_weight * f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += w * dE_dR_over_r * cp_weight * f2_ii;
						atom_f1_f2s_[ upper_res_id ][ jj ].first  += w * dE_dR_over_r * cp_weight * f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += w * dE_dR_over_r * cp_weight * f2_jj;

						atom_f1_f2s_[ lower_res_id ][ ii ].first  += lk * w_f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += lk * w_f2_ii;
						atom_f1_f2s_[ upper_res_id ][ jj ].first  += lk * w_f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += lk * w_f2_jj;
						base_pseudo_atom_f1_f2s_[ lower_res_id ][ ii ].first  += lk * w_f1_ii_pbase;
						base_pseudo_atom_f1_f2s_[ lower_res_id ][ ii ].second += lk * w_f2_ii_pbase;

					}
				}
			} else {
				// ii apolar or water.
				Vector f1_ii, f2_ii, f1_jj, f2_jj;

				Real const dE_dR_over_r(
					eval_dE_dR_over_r(
					lowerres.atom(ii), upperres.atom(jj),
					ii_jj_one_over_d, ii_jj_dis2,
					f1_ii, f2_ii, f1_jj, f2_jj ) );

				if ( dE_dR_over_r != 0.0 ) {
					atom_f1_f2s_[ lower_res_id ][ ii ].first  += dE_dR_over_r * cp_weight * f1_ii;
					atom_f1_f2s_[ lower_res_id ][ ii ].second += dE_dR_over_r * cp_weight * f2_ii;
					atom_f1_f2s_[ upper_res_id ][ jj ].first  += dE_dR_over_r * cp_weight * f1_jj;
					atom_f1_f2s_[ upper_res_id ][ jj ].second += dE_dR_over_r * cp_weight * f2_jj;
				}
			}

			/// 2. Atom ii's desolvation of atom jj
			if ( nneighbs_[ upper_res_id ][ jj ] != 0 ) { // jj polar
				Vector jj_to_ii = -1*ii_to_jj;
				Real dot_prod = jj_to_ii.dot( orientation_vectors_[upper_res_id][ jj ] );
				//std::cout << "Derivative for upperres " << upperres.atom_name( jj ) << " on " << upperres.name() << " dotprod: " << dot_prod << " with  " << lower_res_id << " " << ii << std::endl;

				if ( dot_prod > 0  && dot_prod > LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_LOW) {

					Vector f1_jj( 0.0 ), f2_jj( 0.0 ), f1_ii( 0.0 ), f2_ii( 0.0 );

					Real const dE_dR_over_r(
						eval_dE_dR_over_r(
						upperres.atom(jj), lowerres.atom(ii),
						ii_jj_one_over_d, ii_jj_dis2,
						f1_jj, f2_jj, f1_ii, f2_ii ) );

					if ( dot_prod > LK_SigmoidalFunc::cos_flipped_ANGLE_CUTOFF_HIGH ) {
						/// sigmoidal function evaluates to 1 and its derivative is 0
						atom_f1_f2s_[ upper_res_id ][ jj ].first  += dE_dR_over_r * cp_weight * f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += dE_dR_over_r * cp_weight * f2_jj;
						atom_f1_f2s_[ lower_res_id ][ ii ].first  += dE_dR_over_r * cp_weight * f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += dE_dR_over_r * cp_weight * f2_ii;
					} else { // ( dot_prod < cos_flipped_ANGLE_CUTOFF_HIGH && dot_prod > cos_flipped_ANGLE_CUTOFF_LOW ) {
						/// ramping range of sigmoidal function "w"
						/// Energy is w*lk; chain rule says the derivative is w'*lk + w*lk'
						/// lk' already computed

						// lk
						Real const d2_bin = ii_jj_dis2 * get_bins_per_A2_;
						int	disbin = static_cast< int >( d2_bin ) + 1;
						Real	frac = d2_bin - ( disbin - 1 );
						int const l1 = solv1_.index( disbin, lowerres.atom(ii).type(), upperres.atom(jj).type() );
						int const l2 = l1 + 1;
						Real const lk = lk_hack_weight_ * cp_weight * ( (1.-frac)* solv1_[ l1 ] + frac * solv1_[ l2 ]);

						// w'
						Vector w_f1_jj(0.0), w_f2_jj(0.0), w_f1_ii(0.0), w_f2_ii(0.0), w_f1_jj_pbase(0.0), w_f2_jj_pbase(0.0);
						lk_angle_cst.p1_deriv(
							base_pseudo_atom_centers_[ upper_res_id ][ jj ], upperres.xyz( jj ), lowerres.xyz( ii ),
							w_f1_jj_pbase, w_f2_jj_pbase );

						lk_angle_cst.p2_deriv(
							base_pseudo_atom_centers_[ upper_res_id ][ jj ], upperres.xyz( jj ), lowerres.xyz( ii ),
							w_f1_jj, w_f2_jj );

						lk_angle_cst.p1_deriv(
							lowerres.xyz( ii ), upperres.xyz( jj ), base_pseudo_atom_centers_[ upper_res_id ][ jj ],
							w_f1_ii, w_f2_ii );

						/// w
						Real const w = lkfunc->func(
							angle_radians(
							base_pseudo_atom_centers_[ upper_res_id ][ jj ],
							upperres.xyz( jj ),
							lowerres.xyz( ii ) ) );

						//// Add everything up
						atom_f1_f2s_[ upper_res_id ][ jj ].first  += w * dE_dR_over_r * cp_weight * f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += w * dE_dR_over_r * cp_weight * f2_jj;
						atom_f1_f2s_[ lower_res_id ][ ii ].first  += w * dE_dR_over_r * cp_weight * f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += w * dE_dR_over_r * cp_weight * f2_ii;

						atom_f1_f2s_[ upper_res_id ][ jj ].first  += lk * w_f1_jj;
						atom_f1_f2s_[ upper_res_id ][ jj ].second += lk * w_f2_jj;
						atom_f1_f2s_[ lower_res_id ][ ii ].first  += lk * w_f1_ii;
						atom_f1_f2s_[ lower_res_id ][ ii ].second += lk * w_f2_ii;
						base_pseudo_atom_f1_f2s_[ upper_res_id ][ jj ].first  += lk * w_f1_jj_pbase;
						base_pseudo_atom_f1_f2s_[ upper_res_id ][ jj ].second += lk * w_f2_jj_pbase;

					}
				}
			} else {
				// jj apolar or water.
				Vector f1_jj, f2_jj, f1_ii, f2_ii;

				Real const dE_dR_over_r(
					eval_dE_dR_over_r(
					upperres.atom(jj), lowerres.atom(ii),
					ii_jj_one_over_d, ii_jj_dis2,
					f1_jj, f2_jj, f1_ii, f2_ii ) );

				if ( dE_dR_over_r != 0.0 ) {
					atom_f1_f2s_[ upper_res_id ][ jj ].first  += dE_dR_over_r * cp_weight * f1_jj;
					atom_f1_f2s_[ upper_res_id ][ jj ].second += dE_dR_over_r * cp_weight * f2_jj;
					atom_f1_f2s_[ lower_res_id ][ ii ].first  += dE_dR_over_r * cp_weight * f1_ii;
					atom_f1_f2s_[ lower_res_id ][ ii ].second += dE_dR_over_r * cp_weight * f2_ii;
				}
			}
		}
	}
}

/// @details Evaluate the weighted derivative of atom2's desolvatation of atom1
Real
LK_hack::eval_dE_dR_over_r(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Distance const one_over_d, // known
	DistanceSquared const d2,  // known
	Vector & f1_1,
	Vector & f2_1,
	Vector & f1_2,
	Vector & f2_2
) const
{
	f1_1 = atom1.xyz().cross( atom2.xyz() );
	f2_1 = atom1.xyz() - atom2.xyz();
	f1_2 = -1.0 * f1_1;
	f2_2 = -1.0 * f2_1;

	if ( ( d2 < safe_max_dis2_) && ( d2 != Real(0.0) ) ) {

		Real const d2_bin = d2 * get_bins_per_A2_;
		int	disbin = static_cast< int >( d2_bin ) + 1;
		Real	frac = d2_bin - ( disbin - 1 );

		Real deriv = 0.0;

		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin  ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		int const l1 = dsolv1_.index( disbin, atom2.type(), atom1.type() ),
			l2 = l1 + 1;

		Real e1 = dsolv1_[ l1 ];
		deriv = lk_hack_weight_ * ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );

		//std::cout << "dsolv1_ deriv: " << deriv << " numeric: " << ((solv1_[ l2 ] - solv1_[ l1 ] )*get_bins_per_A2_)*2*sqrt(d2) << std::endl;

		return deriv * one_over_d;
	} else {
		return 0.0;
	}

}

/// @details iterates across polar atoms, and increments the f1/f2 derivative vectors
/// a proportion of the derivative vectors each atom's pseudo-base atom had accumulated
/// into the derivatives of the atom's bonded neighbors
void
LK_hack::distribute_pseudo_base_atom_derivatives( pose::Pose const & pose ) const
{
	Size const total_residue = base_pseudo_atom_f1_f2s_.size();
	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		Size const ii_nheavy = base_pseudo_atom_f1_f2s_[ ii ].size();
		conformation::Residue const & ii_res( pose.residue( ii ) );
		for ( Size jj = 1; jj <= ii_nheavy; ++jj ) {
			if ( nneighbs_[ ii ][ jj ] == 0 ) continue;
			Size jj_nneighbs( nneighbs_[ ii ][ jj ] );
			Real divvy_proportion = 1.0 / jj_nneighbs;

			Vector const f1_to_divvy = divvy_proportion * base_pseudo_atom_f1_f2s_[ ii ][ jj ].first;
			Vector const f2_to_divvy = divvy_proportion * base_pseudo_atom_f1_f2s_[ ii ][ jj ].second;

			Size count_neighbors_found = 0;
			for (Size kk = 1; kk <= ii_res.bonded_neighbor(jj).size(); ++kk){
				Size neighbor_id = ii_res.bonded_neighbor(jj)[kk];
				if ( !  ii_res.atom_is_hydrogen(neighbor_id) ){
					++count_neighbors_found;
					atom_f1_f2s_[ ii ][ neighbor_id ].first  += f1_to_divvy;
					atom_f1_f2s_[ ii ][ neighbor_id ].second += f2_to_divvy;
				}
			}
			if ( ii_res.type().n_residue_connections_for_atom( jj ) > 0 ) {
				for ( Size kk = 1; kk <= ii_res.type().residue_connections_for_atom(jj).size(); ++kk ) {

					chemical::ResConnID const kk_conn = ii_res.connect_map( ii_res.type().residue_connections_for_atom( jj )[ kk ] );
					Size const neighbor_res_id( kk_conn.resid() );
					Size const nieghbor_atom_id( pose.residue( kk_conn.resid() ).residue_connection( kk_conn.connid() ).atomno() );
					if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
						atom_f1_f2s_[ neighbor_res_id ][ nieghbor_atom_id ].first  += f1_to_divvy;
						atom_f1_f2s_[ neighbor_res_id ][ nieghbor_atom_id ].second += f2_to_divvy;
						count_neighbors_found++;
					}
				}
			}
		debug_assert( count_neighbors_found == jj_nneighbs );
		}
	}

}


void
LK_hack::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,// sfxn,
	EnergyMap const & ,//weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( pose.residue( id.rsd() ).atom_is_hydrogen( id.atomno() ) ) return;

	F1 += atom_f1_f2s_[ id.rsd() ][ id.atomno() ].first;
	F2 += atom_f1_f2s_[ id.rsd() ][ id.atomno() ].second;
}


void
LK_hack::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}
core::Size
LK_hack::version() const
{
	return 1; // Initial versioning
}

}
}
}
