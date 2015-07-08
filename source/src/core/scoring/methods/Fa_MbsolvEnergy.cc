// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Fa_MbsolvEnergy.hh
/// @author Patrick Barth


// Unit headers
#include <core/scoring/methods/Fa_MbsolvEnergy.hh>
#include <core/scoring/methods/Fa_MbsolvEnergyCreator.hh>

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


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the Fa_MbsolvEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
Fa_MbsolvEnergyCreator::create_energy_method(
  methods::EnergyMethodOptions const & options
) const {
  return methods::EnergyMethodOP( new Fa_MbsolvEnergy(
  	*( ScoringManager::get_instance()->etable( options ).lock() ),
  	*( ScoringManager::get_instance()->memb_etable( options.etable_type() ).lock() )
  ) );
}

ScoreTypes
Fa_MbsolvEnergyCreator::score_types_for_method() const {
  ScoreTypes sts;
  sts.push_back( fa_mbsolv );
  return sts;
}

Fa_MbsolvEnergy::Fa_MbsolvEnergy( etable::Etable const & etable_in, etable::MembEtable const & memb_etable_in):
  parent( methods::EnergyMethodCreatorOP( new Fa_MbsolvEnergyCreator ) ),
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
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.get_bins_per_A2()),
  verbose_( false ),
  potential_( ScoringManager::get_instance()->get_Membrane_FAPotential() )
{}

Distance
Fa_MbsolvEnergy::atomic_interaction_cutoff() const
{
	return etable_.max_dis();
}

/// clone
EnergyMethodOP
Fa_MbsolvEnergy::clone() const
{
	return EnergyMethodOP( new Fa_MbsolvEnergy( *this ) );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
///
void
Fa_MbsolvEnergy::setup_for_scoring(
  pose::Pose & pose, ScoreFunction const &
) const
{
  potential_.compute_fa_projection( pose );
}


void
Fa_MbsolvEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	Real fa_mbsolv_score( 0.0 );

  get_residue_pair_energy( rsd1, rsd2, pose, fa_mbsolv_score);

	emap[ fa_mbsolv ] += fa_mbsolv_score;
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

  using namespace etable::count_pair;
  CountPairFunctionOP cpfxn( 0 );

  if( same_res ) {
    cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
  } else {
    cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
  }

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const heavy_atom_i( rsd1.xyz( i ) );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

      Real cp_weight = 1.0; Size path_dist( 0 );
      if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {

			  Vector const heavy_atom_j( rsd2.xyz( j ) );

			  Vector const d_ij = heavy_atom_j - heavy_atom_i;
			  Real const d2 = d_ij.length_squared();

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				Real dummy_deriv( 0.0 );

        //pbadebug WARNING
        bool debug(false);

				temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), d2, dummy_deriv,
           Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd1.seqpos(),i),
           Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd2.seqpos(),j),debug);

        if( same_res ) temp_score *= 0.5;

        fa_mbsolv_score += temp_score;

				/*if ( verbose_ && std::abs( fa_mbsolv_score ) > 0.1 )
					std::cout << "fa_mbsolv_score: rsd1 " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " rsd2 " << rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,3,fa_mbsolv_score) << ' ' << std::
endl;*/

			} // cp

		} // j
	} // i

}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
Fa_MbsolvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const & scfxn
) const
{
  potential_.compute_fa_projection( pose );
	pose.update_residue_neighbors();
  fa_mbsolv_weight_ = scfxn.weights()[ fa_mbsolv ];
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

	if ( ( d2 < safe_max_dis2_) && ( d2 != Real(0.0) ) ) {

		Real const d2_bin = d2 * get_bins_per_A2_;
		int	disbin = static_cast< int >( d2_bin ) + 1;
		Real	frac = d2_bin - ( disbin - 1 );

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

    if(debug) {
    std::cout << "f1 s1l1 mbs1l1 s1l2 mbs1l2 " << f1 << " " << solv1_[ l1 ] << " " << memb_solv1_[ l1 ] << " " <<
                  solv1_[ l2 ] << " " << memb_solv1_[ l2 ] << std::endl;
    std::cout << "f2 s2l1 mbs2l1 s2l2 mbs2l2 " << f2 << " " << solv2_[ l1 ] << " " << memb_solv2_[ l1 ] << " " <<
                  solv2_[ l2 ] << " " << memb_solv2_[ l2 ] << std::endl;
    }

		if (eval_deriv) {
			//			int const l1 = dsolv1_.index( disbin, atom2.type(), atom1.type() ),
			//				l2 = l1 + 1;
			//			Real e1 = dsolv1_[ l1 ];
			//			deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
      e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
      e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
      e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
      e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
      e1 = e11 + e21;
      e2 = e12 + e22;

      deriv = e1 + frac * ( e2 - e1 );
      deriv = deriv / std::sqrt( d2 );
		}

	}
	return temp_score;
}
////////////////////////////////////////////////////////////////////////////////////
Real
Fa_MbsolvEnergy::eval_dE_dR_over_r(
  conformation::Atom const & atom1,
  conformation::Atom const & atom2,
  EnergyMap const & /*weights*/,
  Vector & F1,
  Vector & F2,
  Real const & f1,
  Real const & f2
) const
{

  F1 = atom1.xyz().cross( atom2.xyz() );
  F2 = atom1.xyz() - atom2.xyz();
  Real d2,frac;
	
  d2 = atom1.xyz().distance_squared( atom2.xyz() );

  if ( ( d2 < safe_max_dis2_ ) && ( d2 != Real(0.0) ) ) {

     // bin by distance:
     Real const d2_bin = d2 * get_bins_per_A2_;
     int disbin = static_cast< int >( d2_bin ) + 1;
     frac = d2_bin - ( disbin - 1 );

    // l1 and l2 are FArray LINEAR INDICES for fast lookup:
    // [ l1 ] == (disbin  ,attype2,attype1)
    // [ l2 ] == (disbin+1,attype2,attype1)

    //Real deriv = 0.0;

    int const l1 = dsolv1_.index( disbin, atom1.type(), atom2.type()),
      l2 = l1 + 1;

    Real e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
    Real e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
    Real e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
    Real e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
    Real e1 = e11 + e21;
    Real e2 = e12 + e22;

    Real deriv = fa_mbsolv_weight_ * ( e1 + frac * ( e2 - e1 ) );

    return deriv / std::sqrt( d2 );
  } else {
    return 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////
void
Fa_MbsolvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights,
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
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	//	Real deriv( 0.0 );

	for ( graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );
    bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

		using namespace etable::count_pair;
    CountPairFunctionOP cpfxn( 0 );

    if( same_res ) {
      cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
    } else {
      cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
    }

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0; Size path_dist(0);

			if ( cpfxn->count(m, n, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();
				//				Real const d = std::sqrt( d2 );
				Vector const d_ij_norm = d_ij.normalized();

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;


				Vector f1( 0.0 ), f2( 0.0 );

        Real const dE_dR_over_r
           ( eval_dE_dR_over_r( rsd1.atom(m), rsd2.atom(n), weights, f1, f2,
             Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd1.seqpos(),m),
             Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd2.seqpos(),n) ) );
        if ( dE_dR_over_r != 0.0 ) {
          if( same_res ) {
            F1 += 0.5 * dE_dR_over_r * cp_weight * f1;
            F2 += 0.5 * dE_dR_over_r * cp_weight * f2;
          } else {
            F1 += dE_dR_over_r * cp_weight * f1;
            F2 += dE_dR_over_r * cp_weight * f2;
          }
        }
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
  ScoreFunction const &,
  EnergyMap & emap
) const
{

  Real fa_mbsolv_score( 0.0 );

  get_residue_pair_energy( rsd, rsd, pose, fa_mbsolv_score);

  emap[ fa_mbsolv ] += fa_mbsolv_score;

}

////////////////////////////////////////////////
void
Fa_MbsolvEnergy::finalize_total_energy(
	pose::Pose & /*pose*/,
	ScoreFunction const &,
	EnergyMap & /*emap*/
) const
{
 if (verbose_)	std::cout << "DONE SCORING" << std::endl;
}


/// @details Pose must already contain a cenlist object or this method will fail.
Membrane_FAEmbed const &
Fa_MbsolvEnergy::Membrane_FAEmbed_from_pose( pose::Pose const & pose ) const
{
  //using core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED;
  return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
}

MembraneTopology const &
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
}


