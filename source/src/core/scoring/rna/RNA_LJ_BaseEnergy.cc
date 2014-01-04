// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/methods/RNA_LJ_BaseEnergy.hh
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_LJ_BaseEnergy.hh>
#include <core/scoring/rna/RNA_LJ_BaseEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// June 2009. Quick hack to permit additional dispersional attraction between
//  RNA bases. NOTE: It will be worthwhile at some point to make this generically
//  compute rep and atr between all aromatic atoms, not just RNA bases!!
//
//
namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_LJ_BaseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_LJ_BaseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new RNA_LJ_BaseEnergy( *ScoringManager::get_instance()->etable( options.etable_type() ) );
}

ScoreTypes
RNA_LJ_BaseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_fa_atr_base );
	sts.push_back( rna_fa_rep_base );
	return sts;
}


RNA_LJ_BaseEnergy::RNA_LJ_BaseEnergy( etable::Etable const & etable_in ) :
	parent( new RNA_LJ_BaseEnergyCreator ),
	etable_( etable_in ),
	ljatr_( etable_in.ljatr() ),
	ljrep_( etable_in.ljrep() ),
	dljatr_( etable_in.dljatr() ),
	dljrep_( etable_in.dljrep() ),
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.get_bins_per_A2() ),
	verbose_( false )
{
	if ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] ){
		utility_exit_with_message(  "RNA_LJ_BaseEnergy not compatible with analytic_etable_evaluation yet -- rerun with flag -analytic_etable_evaluation false." );
	}

}

Distance
RNA_LJ_BaseEnergy::atomic_interaction_cutoff() const
{
	return etable_.max_dis();
}

// clone
methods::EnergyMethodOP
RNA_LJ_BaseEnergy::clone() const
{
	return new RNA_LJ_BaseEnergy( *this );
}

////////////////////////////////////////////////
RNA_LJ_BaseEnergy::RNA_LJ_BaseEnergy( RNA_LJ_BaseEnergy const & src ):
	parent( src ),
	etable_( src.etable_ ),
	ljatr_( src.ljatr_ ),
	ljrep_( src.ljrep_ ),
	dljatr_( src.dljatr_ ),
	dljrep_( src.dljrep_ ),
	safe_max_dis2_( src.safe_max_dis2_ ),
	get_bins_per_A2_( src.get_bins_per_A2_   ),
	verbose_( src.verbose_ )
{}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
RNA_LJ_BaseEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	//
	// Currently only works for RNA bases!!!
	// Could easily make it more general by checking for atoms that are *aromatic*
	// This information could be kept in the Residue (and ResidueType) class, e.g., is_aromatic( atom_index ).
	//
	if ( !rsd1.is_RNA() ) return ;
	if ( !rsd2.is_RNA() ) return ;

	Real fa_atr_score( 0.0 ), fa_rep_score( 0.0 );
	// Basically just re-calculating lennard jones terms.
	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	//Real lj_aro_score =  0.0;

	for ( Size i = rsd1.first_sidechain_atom() + 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const heavy_atom_i( rsd1.xyz( i ) );

		for ( Size j = rsd2.first_sidechain_atom() + 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( j ) );

				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();

				if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real( 0.0 ) ) ) continue;

				//Real dotprod( 1.0 );
				Real dummy_deriv( 0.0 );

				Real temp_atr_score( 0.0 ), temp_rep_score( 0.0 );

				eval_lj( rsd1.atom( i ), rsd2.atom( j ), d2,
								 temp_atr_score, temp_rep_score,
								 dummy_deriv, dummy_deriv );

				fa_atr_score += cp_weight * temp_atr_score;
				fa_rep_score += cp_weight * temp_rep_score;

			} // cp

		} // j
	} // i

	emap[ rna_fa_atr_base ] += fa_atr_score;
	emap[ rna_fa_rep_base ] += fa_rep_score;

}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
RNA_LJ_BaseEnergy::setup_for_derivatives(
																				 pose::Pose & pose,
																				 ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
}


////////////////////////////////////////////////
void
RNA_LJ_BaseEnergy::eval_lj(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const & d2,
	Real & fa_atr_score,
	Real & fa_rep_score,
	Real & deriv_atr,
	Real & deriv_rep
) const
{

	//	Real temp_score( 0.0 );
	deriv_atr = 0.0;
	deriv_rep = 0.0;

	//Make this an input option for efficiency
	bool const eval_deriv( true );

	if ( ( d2 < safe_max_dis2_ ) && ( d2 != Real( 0.0 ) ) ) {

		Real const d2_bin = d2 * get_bins_per_A2_;
		int	disbin = static_cast< int > ( d2_bin ) + 1;
		Real	frac = d2_bin - ( disbin - 1 );

		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		{
			int const l1 = ljatr_.index( disbin, atom2.type(), atom1.type() );
			int const l2 = l1 + 1;
			fa_atr_score = ( ( 1. - frac )* ljatr_[ l1 ] + frac * ljatr_[ l2 ] );
			if ( eval_deriv ) {
				deriv_atr = ( ljatr_[ l2 ] - ljatr_[ l1 ] ) * get_bins_per_A2_ * std::sqrt( d2 ) * 2;
			}
		}


		{
			int const l1 = ljrep_.index( disbin, atom2.type(), atom1.type() );
			int const l2 = l1 + 1;
			fa_rep_score = ( ( 1. - frac )* ljrep_[ l1 ] + frac * ljrep_[ l2 ] );
			if ( eval_deriv ) {
				deriv_rep = ( ljrep_[ l2 ] - ljrep_[ l1 ] ) * get_bins_per_A2_ * std::sqrt( d2 ) * 2;
			}
		}

	}

}

//////////////////////////////////////////////////////////////////////////////////////
void
RNA_LJ_BaseEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	// Currently extremely RNA specific.
	if ( !rsd1.is_RNA() ) return;
	if ( m > rsd1.nheavyatoms() ) return;
	if ( m < rsd1.first_sidechain_atom() + 1 ) return;

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

		Size const j( ( *iter )->get_other_ind( i ) );

		if ( pos1_fixed && domain_map( i ) == domain_map( j ) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );

		using namespace etable::count_pair;
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		for ( Size n = rsd2.first_sidechain_atom() + 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( m, n, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();
				//Real const d = std::sqrt( d2 );
				Vector const d_ij_norm = d_ij.normalized();

				if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real( 0.0 ) ) ) continue;

				Real dummy( 0.0 ), fa_atr_deriv( 0.0 ), fa_rep_deriv( 0.0 );
				eval_lj( rsd1.atom( m ), rsd2.atom( n ), d2, dummy, dummy, fa_atr_deriv, fa_rep_deriv );

				Vector const f2_fwd =   -1.0 * cp_weight * d_ij_norm * ( fa_atr_deriv * weights[ rna_fa_atr_base ] +
																										fa_rep_deriv * weights[ rna_fa_rep_base ] );
				Vector const f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

				F1 += f1_fwd;
				F2 += f2_fwd;

      }


    }
  }

}




//////////////////////////////////////////////////////////////////////////////////////
// hacky hacky.
Real
RNA_LJ_BaseEnergy::eval_atom_energy(
	id::AtomID const & atom_id,
	pose::Pose const & pose
) const
{

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	Real score( 0.0 );

	// Currently extremely RNA specific.
	if ( !rsd1.is_RNA() ) return score;
	if ( m > rsd1.nheavyatoms() ) return score;
	if ( m < rsd1.first_sidechain_atom() + 1 ) return score;

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	for ( Size j = 1; j < pose.total_residue(); j ++ ) {

		if ( i == j ) continue;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_RNA() ) continue;

		using namespace etable::count_pair;
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		for ( Size n = rsd2.first_sidechain_atom() + 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( m, n, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();
				//Real const d = std::sqrt( d2 );
				Vector const d_ij_norm = d_ij.normalized();

				if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real( 0.0 ) ) ) continue;

				Real dummy( 0.0 ), fa_atr( 0.0 ), fa_rep( 0.0 );
				eval_lj( rsd1.atom( m ), rsd2.atom( n ), d2, fa_atr, fa_rep, dummy, dummy );

				// In principle could pass in an emap and weight the components
				// by the Emap.
				score += cp_weight * ( fa_atr + fa_rep );
      }


    }
  }

	return score;
}


////////////////////////////////////////////////
void
RNA_LJ_BaseEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}

////////////////////////////////////////////////
void
RNA_LJ_BaseEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const
{
	if ( verbose_ )	std::cout << "DONE SCORING" << std::endl;
}
core::Size
RNA_LJ_BaseEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}



