// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_PolarNonPolarEnergy.hh
/// @author Rhiju Das


// Unit headers
#include <core/scoring/methods/LK_PolarNonPolarEnergy.hh>
#include <core/scoring/methods/LK_PolarNonPolarEnergyCreator.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/RNA_Util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

#include <ObjexxFCL/format.hh>

#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL::fmt;

/// @details This must return a fresh instance of the LK_PolarNonPolarEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
LK_PolarNonPolarEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new LK_PolarNonPolarEnergy( *( ScoringManager::get_instance()->etable( options.etable_type() )),
																options.analytic_etable_evaluation() );
}

ScoreTypes
LK_PolarNonPolarEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( lk_costheta );
	sts.push_back( lk_polar );
	sts.push_back( lk_nonpolar );
	sts.push_back( lk_polar_intra_RNA );
	sts.push_back( lk_nonpolar_intra_RNA );
	return sts;
}


LK_PolarNonPolarEnergy::LK_PolarNonPolarEnergy( etable::Etable const & etable_in, bool const analytic_etable_evaluation ):
	parent( new LK_PolarNonPolarEnergyCreator ),
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	max_dis_( etable_in.max_dis() ),
	verbose_( false )
{
	if ( analytic_etable_evaluation ) {
		etable_evaluator_ = new etable::AnalyticEtableEvaluator( etable_in );
	} else {
		etable_evaluator_ = new etable::TableLookupEvaluator( etable_in );
	}
}

////////////////////////////////////////////////
LK_PolarNonPolarEnergy::LK_PolarNonPolarEnergy( LK_PolarNonPolarEnergy const & src ):
	parent( src ),
	safe_max_dis2_( src.safe_max_dis2_ ),
	max_dis_( src.max_dis_ ),
	etable_evaluator_( src.etable_evaluator_ ),
	verbose_( src.verbose_ )
{}


void
LK_PolarNonPolarEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const{

	if(rsd.is_RNA()==false) return;

	Real lk_polar_intra_RNA_score, lk_nonpolar_intra_RNA_score, lk_costheta_intra_RNA_score;
	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );

	get_residue_energy_RNA_intra( rsd, pose, lk_polar_intra_RNA_score, lk_nonpolar_intra_RNA_score, lk_costheta_intra_RNA_score, compute_polar, compute_nonpolar );
	emap[ lk_polar_intra_RNA ]    += lk_polar_intra_RNA_score;
	emap[ lk_nonpolar_intra_RNA ] += lk_nonpolar_intra_RNA_score;

}

bool
LK_PolarNonPolarEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	bool method_1 = ( weights[lk_polar_intra_RNA] > 0.0 || weights[lk_nonpolar_intra_RNA] > 0.0 ) ? true : false;

	return method_1;
}


Distance
LK_PolarNonPolarEnergy::atomic_interaction_cutoff() const
{
	return max_dis_;
}

/// clone
EnergyMethodOP
LK_PolarNonPolarEnergy::clone() const
{
	return new LK_PolarNonPolarEnergy( *this );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	Real lk_polar_score, lk_nonpolar_score, lk_costheta_score;

	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );

	get_residue_pair_energy_one_way( rsd1, rsd2, pose, lk_polar_score, lk_nonpolar_score, lk_costheta_score,
																	 compute_polar, compute_nonpolar );
	emap[ lk_polar ]    += lk_polar_score;
	emap[ lk_nonpolar ] += lk_nonpolar_score;
	emap[ lk_costheta ] += lk_costheta_score;

	get_residue_pair_energy_one_way( rsd2, rsd1, pose, lk_polar_score, lk_nonpolar_score, lk_costheta_score,
																	 compute_polar, compute_nonpolar );
	emap[ lk_polar ]    += lk_polar_score;
	emap[ lk_nonpolar ] += lk_nonpolar_score;
	emap[ lk_costheta ] += lk_costheta_score;

}

////////////////////////////////////////////////
Vector
LK_PolarNonPolarEnergy::get_base_vector( conformation::Residue const & rsd1, Size const i, pose::Pose const & pose ) const
{
	// Use DB/APL prescription of looking through bonded neighbors.
	Size non_H_neighbors = 0;
	Vector  base_pseudo_atom(0);
	for (Size ii = 1; ii <=rsd1.bonded_neighbor(i).size(); ++ii){
		Size neighbor_id = rsd1.bonded_neighbor(i)[ii];
		if ( !  rsd1.atom_is_hydrogen(neighbor_id)){
			base_pseudo_atom += rsd1.xyz(neighbor_id);
			non_H_neighbors++;
		}
	}

	if ( rsd1.type().n_residue_connections_for_atom(i) > 0  ) {
		/// CONTEXT DEPENDENCY HERE -- e.g. if c_prev moves, rsd1 needs to be rescored.  Fortunately,
		/// if c_prev moves, the internal "psi" for rsd1 will be updated, and this residue will be rescored.

		for ( Size ii = 1; ii <= rsd1.type().residue_connections_for_atom(i).size(); ++ii ) {
			chemical::ResConnID const ii_conn = rsd1.connect_map( rsd1.type().residue_connections_for_atom(i)[ ii ] );
			Size const neighbor_res_id(  ii_conn.resid() );
			if (neighbor_res_id < 1 ) continue;
			Size const nieghbor_atom_id( pose.residue( ii_conn.resid() ).residue_connection( ii_conn.connid() ).atomno() );
			if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
				base_pseudo_atom += pose.residue( neighbor_res_id ).xyz( nieghbor_atom_id );
				non_H_neighbors++;
			}
		}
	}

	if (non_H_neighbors > 0 ) 	base_pseudo_atom /= non_H_neighbors;

	// Note -- probably should have a big WARNING show up
	// if we're trying to compute this scorefunction for water or something
	// where there is no base atom.  That's where non_H_neighbors is 0.

	Vector res1_base_vector =  rsd1.xyz(i) - base_pseudo_atom;
	Vector res1_base_vector_norm = (res1_base_vector).normalized();
	return res1_base_vector_norm;
}

////////////////////////////////////////////////
// Why is all this code copied? Argh. -- rhiju
void
LK_PolarNonPolarEnergy::get_residue_energy_RNA_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real & lk_polar_intra_RNA_score,
	Real & lk_nonpolar_intra_RNA_score,
	Real & lk_costheta_intra_RNA_score,
	bool const compute_polar,
	bool const compute_nonpolar
) const
{

	using namespace etable::count_pair;

	lk_polar_intra_RNA_score =  0.0;
	lk_nonpolar_intra_RNA_score =  0.0;
	lk_costheta_intra_RNA_score =  0.0;

	conformation::Residue const & rsd1=rsd; //Assume rsd1 contributes polar atoms.
	conformation::Residue const & rsd2=rsd; //Assume rsd2 contributes occluding atoms.

	CountPairFunctionOP cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd, CP_CROSSOVER_4); //intra_res version!

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const heavy_atom_i( rsd1.xyz( i ) );

		//Just compute cos(theta) for polars.
		bool is_polar( false );
		if ( rsd1.atom_type(i).is_acceptor() || rsd1.atom_type(i).is_donor()) is_polar = true;

		if ( is_polar  && !compute_polar    ) continue;
		if ( !is_polar && !compute_nonpolar ) continue;

		//Need to figure out "direction"...
		Vector const res1_base_vector_norm = is_polar ? get_base_vector( rsd1, i, pose ) : Vector( 0.0 );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {

				if( !Is_base_phosphate_atom_pair(rsd1, rsd2, i, j) ) continue;

				Vector const heavy_atom_j( rsd2.xyz( j ) );

				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				Real dotprod( 1.0 );
				Real dummy_deriv( 0.0 );
				Real temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), dummy_deriv);

				if ( is_polar ) {
					lk_polar_intra_RNA_score += temp_score;

					Vector const d_ij = heavy_atom_j - heavy_atom_i;
					Vector const d_ij_norm = d_ij.normalized();
					dotprod = dot( res1_base_vector_norm, d_ij_norm );
					temp_score *= dotprod;
					lk_costheta_intra_RNA_score += temp_score;

					if ( verbose_ && std::abs( temp_score ) > 0.1 ){
						std::cout << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
							rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,3,temp_score) << ' ' << F(8,3,dotprod) << std::endl;
					}

				} else {
					lk_nonpolar_intra_RNA_score += temp_score;
				}

			} // cp

		} // j
	} // i

}

////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::get_residue_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & lk_polar_score,
	Real & lk_nonpolar_score,
	Real & lk_costheta_score,
	bool const compute_polar,
	bool const compute_nonpolar
) const
{
	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	lk_polar_score =  0.0;
	lk_nonpolar_score =  0.0;
	lk_costheta_score =  0.0;

	//Assume rsd1 contributes polar atoms.
	//Assume rsd2 contributes occluding atoms.
	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const & heavy_atom_i( rsd1.xyz( i ) );

		//Just compute cos(theta) for polars.
		bool const is_polar = ( rsd1.atom_type(i).is_acceptor() || rsd1.atom_type(i).is_donor());

		if ( is_polar  && !compute_polar    ) continue;
		if ( !is_polar && !compute_nonpolar ) continue;

		//Need to figure out "direction"...
		Vector const res1_base_vector_norm = is_polar ? get_base_vector( rsd1, i, pose ) : Vector( 0.0 );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( j ) );

				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				Real dotprod( 1.0 );
				Real dummy_deriv( 0.0 );
				Real temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), dummy_deriv);

				if ( is_polar ) {
					lk_polar_score += temp_score;

					Vector const d_ij = heavy_atom_j - heavy_atom_i;
					Vector const d_ij_norm = d_ij.normalized();
					dotprod = dot( res1_base_vector_norm, d_ij_norm );
					temp_score *= dotprod;
					lk_costheta_score += temp_score;

					if ( verbose_ && std::abs( temp_score ) > 0.1 ){
						std::cout << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
							rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,3,temp_score) << ' ' << F(8,3,dotprod) << std::endl;
					}

				} else {
					lk_nonpolar_score += temp_score;
				}


			} // cp

		} // j
	} // i

}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::setup_for_derivatives(
																				 pose::Pose & pose,
																				 ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
}


////////////////////////////////////////////////
Real
LK_PolarNonPolarEnergy::eval_lk(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & deriv ) const
{

	Real temp_score( 0.0 );
	deriv = 0.0;
	bool const eval_deriv( true );

	etable_evaluator_->atom_pair_lk_energy_and_deriv_v( atom1, atom2, temp_score, deriv, eval_deriv );

	return temp_score;
}

//////////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::eval_atom_derivative_intra_RNA(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{


	bool do_eval_intra_RNA= (weights[lk_polar_intra_RNA]>0.0 || weights[lk_nonpolar_intra_RNA]>0.0) ? true : false;

	if(do_eval_intra_RNA==false) return; //early return.

	Size const i( atom_id.rsd() );
	Size const j( atom_id.rsd() );

	conformation::Residue const & rsd1( pose.residue( i ) );
	conformation::Residue const & rsd2( pose.residue( j ) );

	if(rsd1.is_RNA()==false) return;
	if(rsd2.is_RNA()==false) return; //no effect!

	Size const m( atom_id.atomno() );

	if ( m > rsd1.nheavyatoms() ) return;

	if(verbose_){
		std::cout << "Start LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res case, res= " << i << " atomno= " << m << "[" << rsd1.atom_name(m) <<  "]" << std::endl;
	}

	//bool const pos1_fixed( domain_map( i ) != 0 );

	//if( pos1_fixed && domain_map(i) == domain_map(j) ){ //MOD OUT ON July 19th, 2011...THIS MIGHT BE BUGGY!
		 //std::cout << "LK_PolarNonPolarEnergy::eval_atom_derivative early return since pos1_fixed && domain_map(i=" << i << ") == domain_map(j=" << j << ")" << std::endl;
	//	 return; //Fixed w.r.t. one another.
	//}

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	Real deriv( 0.0 );
	Vector const res1_base_vector_norm = get_base_vector( rsd1, m, pose );

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4); //intra_res version!

	bool atom1_is_polar( false );
	if (rsd1.atom_type(m).is_acceptor() || rsd1.atom_type(m).is_donor() ) atom1_is_polar = true;

	for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

		Real cp_weight = 1.0;
		Size path_dist( 0 );
		if ( cpfxn->count(m, n, cp_weight, path_dist ) ) {

			if(Is_base_phosphate_atom_pair(rsd1, rsd2, m, n)==false) continue;

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();
			Real const d = std::sqrt( d2 );
			Vector const d_ij_norm = d_ij.normalized();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			bool atom2_is_polar( false );
			if (rsd2.atom_type(n).is_acceptor() || rsd2.atom_type(n).is_donor() ) atom2_is_polar = true;

			Vector const res2_base_vector_norm = get_base_vector( rsd2, n, pose );

			// Real const dist_ij = d_ij.length();

			Vector f1_fwd( 0.0 ), f2_fwd( 0.0 ), f1_bkd( 0.0 ), f2_bkd( 0.0 );
			Real lk_score1( 0.0 ), lk_score2( 0.0 );
			Real dotprod_fwd( 1.0 ), dotprod_bkd( 1.0 );

			//Forward direction first.
			lk_score1 = cp_weight * eval_lk( rsd1.atom(m), rsd2.atom(n), deriv );

			f2_fwd =   -1.0 * cp_weight * deriv * d_ij_norm;
			f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

			if ( atom1_is_polar ) {

				F1 += weights[ lk_polar_intra_RNA ] * f1_fwd;
				F2 += weights[ lk_polar_intra_RNA ] * f2_fwd;

				dotprod_fwd = dot( res1_base_vector_norm, d_ij_norm );

				f2_fwd *= dotprod_fwd;
				f2_fwd -= lk_score1 * ( 1/d ) *  (res1_base_vector_norm  - dotprod_fwd * d_ij_norm );

				f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

				//F1 += weights[ lk_costheta ] * f1_fwd;
				//F2 += weights[ lk_costheta ] * f2_fwd;

				lk_score1 *= dotprod_fwd; //to check later (verbose)
			} else {

				F1 += weights[ lk_nonpolar_intra_RNA  ] * f1_fwd;
				F2 += weights[ lk_nonpolar_intra_RNA  ] * f2_fwd;

			}

			/////////////////////////////////
			// Backwards
			Vector d_ji_norm = -d_ij_norm;

			lk_score2 = cp_weight * eval_lk( rsd2.atom(n), rsd1.atom(m),  deriv );

			f2_bkd =   -1.0 * deriv * cp_weight * d_ji_norm;
			f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );


			if ( atom2_is_polar ){

				F1 -= weights[ lk_polar_intra_RNA ] * f1_bkd;
				F2 -= weights[ lk_polar_intra_RNA ] * f2_bkd;

				dotprod_bkd = dot( res2_base_vector_norm, d_ji_norm );

				f2_bkd *= dotprod_bkd;
				f2_bkd -= lk_score2 * ( 1/d ) *  (res2_base_vector_norm  - dotprod_bkd * d_ji_norm );

				f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );

				//F1 -= weights[ lk_costheta ] * f1_bkd;
				//F2 -= weights[ lk_costheta ] * f2_bkd;

				lk_score2 *= dotprod_bkd; //to check later (verbose)
			} else {

				F1 -= weights[ lk_nonpolar_intra_RNA ] * f1_bkd;
				F2 -= weights[ lk_nonpolar_intra_RNA  ] * f2_bkd;

			}

			if ( verbose_  &&  (std::abs( lk_score1 ) > 0.1 || std::abs( lk_score2 ) > 0.1 ) ) {
				std::cout << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( m ) << " covered by " <<
					rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(n) <<
					" " << F(8,3,lk_score1) << " " << F(8,3,lk_score2) <<
					" ==> " << " DERIV " <<
					F(8,6,f2_fwd( 1 ) ) <<  ' ' << F(8,6,f2_bkd(1) ) <<
					' ' << std::sqrt( d2 )	<< " " << cp_weight << " " <<
					F(8,3,dotprod_fwd) << " " << F(8,3,dotprod_bkd) << std::endl;
			}


    }


  }

	if(verbose_){
		std::cout << "LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res :";
	 	std::cout << " F1= " << F1[0] << " " << F1[1] << " " << F1[2];
	 	std::cout << " F2= " << F2[0] << " " << F2[1] << " " << F2[2] << std::endl;
		std::cout << "Finish LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res case, res= " << i << " atomno= " << m << "[" << rsd1.atom_name(m) <<  "]" << std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	eval_atom_derivative_intra_RNA(atom_id, pose, domain_map, weights, F1, F2);

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

	Real deriv( 0.0 );
	Vector const res1_base_vector_norm = get_base_vector( rsd1, m, pose );

	for ( graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );

		using namespace etable::count_pair;
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		bool atom1_is_polar( false );
		if (rsd1.atom_type(m).is_acceptor() || rsd1.atom_type(m).is_donor() ) atom1_is_polar = true;

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( cpfxn->count(m, n, cp_weight, path_dist ) ) {

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();
				Real const d = std::sqrt( d2 );
				Vector const d_ij_norm = d_ij.normalized();

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				bool atom2_is_polar( false );
				if (rsd2.atom_type(n).is_acceptor() || rsd2.atom_type(n).is_donor() ) atom2_is_polar = true;

				Vector const res2_base_vector_norm = get_base_vector( rsd2, n, pose );

				// Real const dist_ij = d_ij.length();

				Vector f1_fwd( 0.0 ), f2_fwd( 0.0 ), f1_bkd( 0.0 ), f2_bkd( 0.0 );
				Real lk_score1( 0.0 ), lk_score2( 0.0 );
				Real dotprod_fwd( 1.0 ), dotprod_bkd( 1.0 );

				//Forward direction first.
				lk_score1 = cp_weight * eval_lk( rsd1.atom(m), rsd2.atom(n), deriv );

				f2_fwd =   -1.0 * cp_weight * deriv * d_ij_norm;
				f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

				if ( atom1_is_polar ) {

					F1 += weights[ lk_polar ] * f1_fwd;
					F2 += weights[ lk_polar ] * f2_fwd;

					dotprod_fwd = dot( res1_base_vector_norm, d_ij_norm );

					f2_fwd *= dotprod_fwd;
					f2_fwd -= lk_score1 * ( 1/d ) *  (res1_base_vector_norm  - dotprod_fwd * d_ij_norm );

					f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

					F1 += weights[ lk_costheta ] * f1_fwd;
					F2 += weights[ lk_costheta ] * f2_fwd;

					lk_score1 *= dotprod_fwd; //to check later (verbose)
				} else {

					F1 += weights[ lk_nonpolar ] * f1_fwd;
					F2 += weights[ lk_nonpolar ] * f2_fwd;

				}

				/////////////////////////////////
				// Backwards
				Vector d_ji_norm = -d_ij_norm;

				lk_score2 = cp_weight * eval_lk( rsd2.atom(n), rsd1.atom(m),  deriv );

				f2_bkd =   -1.0 * deriv * cp_weight * d_ji_norm;
				f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );


				if ( atom2_is_polar ){

					F1 -= weights[ lk_polar ] * f1_bkd;
					F2 -= weights[ lk_polar ] * f2_bkd;

					dotprod_bkd = dot( res2_base_vector_norm, d_ji_norm );

					f2_bkd *= dotprod_bkd;
					f2_bkd -= lk_score2 * ( 1/d ) *  (res2_base_vector_norm  - dotprod_bkd * d_ji_norm );

					f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );

					F1 -= weights[ lk_costheta ] * f1_bkd;
					F2 -= weights[ lk_costheta ] * f2_bkd;

					lk_score2 *= dotprod_bkd; //to check later (verbose)
				} else {

					F1 -= weights[ lk_nonpolar ] * f1_bkd;
					F2 -= weights[ lk_nonpolar ] * f2_bkd;

				}

				if ( verbose_  &&  (std::abs( lk_score1 ) > 0.1 || std::abs( lk_score2 ) > 0.1 ) ) {
					std::cout << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( m ) << " covered by " <<
						rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(n) <<
						" " << F(8,3,lk_score1) << " " << F(8,3,lk_score2) <<
						" ==> " << " DERIV " <<
						F(8,6,f2_fwd( 1 ) ) <<  ' ' << F(8,6,f2_bkd(1) ) <<
						' ' << std::sqrt( d2 )	<< " " << cp_weight << " " <<
						F(8,3,dotprod_fwd) << " " << F(8,3,dotprod_bkd) << std::endl;
				}


      }


    }
  }

}


////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}

////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const
{
	if (verbose_)	std::cout << "DONE SCORING" << std::endl;
}
core::Size
LK_PolarNonPolarEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}



