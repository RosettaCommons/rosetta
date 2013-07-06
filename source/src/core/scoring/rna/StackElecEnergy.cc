// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/StackElecEnergy.cc
/// @brief  HackElec, but just 'perpendicular' to bases... trying to separate from H-bonds & geom_sol.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/StackElecEnergy.hh>
#include <core/scoring/rna/StackElecEnergyCreator.hh>

// Package headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::fmt;

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


////////////////////////////////////////////////////////////////////////////////////
// StackElecEnergy: Mash-up of HackElecEnergy & RNA_FullAtomStackingEnergy.
//
//  We need an approximate treatment of electrostatics, to handle several
//   puzzles in RNA modeling, including penalties  for 5'-CG-3' turner rule,
//   favored structures in the tandem GA motif, etc.
//
//  Unfortunately, just using hack_elec doesn't solve this -- severe convolution
//   with hydrogen bonds. This stack_elec function separates electrostatics
//   'perpendicular' to nucleobases away from hydrogen bonds & geom_sol which
//   have been balanced to account for in-plane interactions (base pairing).
//
//  Inspired by discussions with Fang-Chieh Chou.
//
//     -- Rhiju, July 2013
//
////////////////////////////////////////////////////////////////////////////////////

// C++
static basic::Tracer tr("core.scoring.rna.StackElecEnergy");

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the StackElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
StackElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new StackElecEnergy( options );
}

ScoreTypes
StackElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( stack_elec );
	sts.push_back( stack_elec_base_only );
	return sts;
}

/// c-tor
StackElecEnergy::StackElecEnergy( methods::EnergyMethodOptions const & options ) :
	parent( new StackElecEnergyCreator ),
	coulomb_( options ),
	base_base_only_( false ), //true will be faster computation but appears less accurate (and less physically consistent)
	verbose_( false )
{
	coulomb_.initialize();
}

////////////////////////////////////////////////////////////////////////////
StackElecEnergy::StackElecEnergy( StackElecEnergy const & src ):
	parent( src ),
	coulomb_( src.coulomb() ),
	base_base_only_( src.base_base_only_ ),
	verbose_( src.verbose_ )
{
	coulomb_.initialize();
}


//clone
methods::EnergyMethodOP
StackElecEnergy::clone() const
{
  return new StackElecEnergy( *this );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


///
void
StackElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
  pose.update_residue_neighbors();

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.update( pose );

}

void
StackElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
  pose.update_residue_neighbors();

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.update( pose );

}

//////////////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{


	Real score_base1( 0.0 ), score_base2( 0.0 );

	Real const score = residue_pair_energy_one_way( rsd1, rsd2, pose, score_base1 ) +
		residue_pair_energy_one_way( rsd2, rsd1, pose, score_base2 ) ;

  emap[ stack_elec ]       += score;

	emap[ stack_elec_base_only ]   += score_base1 + score_base2;

	if ( verbose_ && std::abs( score ) > 0.01 ) {
		tr.Info << "respair " << rsd1.name3()  << rsd1.seqpos() << " --- " << rsd2.name3() << rsd2.seqpos() << ": " << F( 8, 3,score) << std::endl;
	}

}


//////////////////////////////////////////////////////////////////////////////////////////
Real
StackElecEnergy::residue_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & score_base
) const
{

	score_base = 0.0;

	// currently defined only for RNA -- could easily generalize to DNA or proteins,
	//  as long as we precompute centroid & base-normals for rings.
	if ( !rsd1.is_RNA() ) return 0.0;

  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

  Real score( 0.0 );

  Size const i( rsd1.seqpos() );
  kinematics::Stub const & stub_i( base_stubs[i] );
  Matrix const M_i ( stub_i.M );

  // Loop over base heavy atoms.
  // If I want to generalize this to proteins, maybe could loop over "aromatic" atoms.
  for ( Size m = 1; m <= rsd1.natoms(); ++m )  {

		// following contains virtual check.
		if ( !is_rna_base( rsd1, m ) ) continue;

		Real const i_charge( rsd1.atomic_charge(m) );
		if ( i_charge == 0.0 ) continue;

    Vector const atom_i( rsd1.xyz( m ) );

    for ( Size n = 1; n <= rsd2.natoms(); ++n ) {

			if ( rsd2.is_virtual( n ) ) continue;
			if ( base_base_only_ && !is_rna_base( rsd2, n ) ) continue;

			Real const j_charge( rsd2.atomic_charge(n) );
			if ( j_charge == 0.0 ) continue;

      Vector const atom_j( rsd2.xyz( n ) );

			Real cos_kappa2( 0.0 ); // useful for output...
			Real const stack_elec_score = get_stack_elec_score( atom_i, atom_j, i_charge, j_charge, M_i, cos_kappa2 );
			score += stack_elec_score;

			if ( is_rna_base( rsd1, m) && is_rna_base( rsd2, n) ) score_base += stack_elec_score;

			//DEBUG
			if ( std::abs( stack_elec_score ) > 0.1 ){
				tr << rsd1.name1() << I( 2, rsd1.seqpos() ) ;
				tr << " " << rsd1.type().atom_name(m) << " " << A( 4, rsd1.atom_type(m).name() );
				tr << " [" << F(8, 3, i_charge) << "]";
				tr << " ---  ";
				tr << rsd2.name1() << I( 2, rsd2.seqpos() );
				tr << " " << rsd2.type().atom_name(n) << " " << A( 4, rsd2.atom_type(n).name() );
				tr << " [" << F(8, 3, j_charge) << "]";
				tr << ": " << F(8, 3, stack_elec_score);
				tr <<  " [ coskappa2 " << F(8, 3, cos_kappa2 ) << "]" << std::endl;
			}

    }
  }


	return score;
}


//////////////////////////////////////////////////////////////////////////////
bool
StackElecEnergy::is_rna_base(
	conformation::Residue const & rsd,
  Size const & m ) const {

	if ( !rsd.type().is_RNA() ) return false;

	//Need to be careful nheavyatoms count includes hydrogen! when the hydrogen is made virtual..
	if ( rsd.is_virtual( m ) ) return false;

	if ( rsd.atom_is_hydrogen( m ) ) return is_rna_base( rsd, rsd.atom_base( m ) );

	// Following is really really specific to RNA.
	return ( m > rsd.first_sidechain_atom() && m <= rsd.nheavyatoms() );

}

////////////////////////////////////////////////////////////////////////////// (Need to condition this? Parin Sep 2, 2009)
void
StackElecEnergy::eval_atom_derivative(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const
{

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	if ( rsd1.is_virtual( m ) ) return;

	Real const i_charge( rsd1.atomic_charge(m) );
	if ( i_charge == 0.0 ) return;

	if ( base_base_only_ && !is_rna_base( rsd1, m ) ) return;

  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	Vector const atom_i( rsd1.xyz( m ) );

	bool const pos1_fixed( domain_map( i ) != 0 );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );

		for ( Size n = 1; n <= rsd2.natoms(); ++n ) {

			if( rsd2.is_virtual(n) ) continue;

			Real const j_charge( rsd2.atomic_charge(n) );
			if ( j_charge == 0.0 ) continue;

			if ( base_base_only_ && !is_rna_base( rsd2, n ) ) continue;

			Vector const atom_j( rsd2.xyz( n ) );

			if ( is_rna_base( rsd1, m ) ) {

				kinematics::Stub const & stub_i( base_stubs[i] );
				Matrix const M_i ( stub_i.M );

				Vector const & deriv_vector_i = get_stack_elec_deriv( atom_i, atom_j, i_charge, j_charge, M_i );

				Vector force_vector_i = weights[ stack_elec ] * deriv_vector_i;

				if ( weights[ stack_elec_base_only ] != 0.0 && is_rna_base( rsd2, n ) ) force_vector_i += weights[ stack_elec_base_only ] * deriv_vector_i;

				//Force/torque with which occluding atom j acts on "dipole" i.
				F1 += -1.0 * cross( force_vector_i, atom_j );
				F2 += -1.0 * force_vector_i;
			}

			if ( is_rna_base( rsd2, n )  ){
					//Force/torque with which occluding atom i acts on "dipole" j.

				kinematics::Stub const & stub_j( base_stubs[j] );
				Matrix const M_j ( stub_j.M );

				Vector const & deriv_vector_j = get_stack_elec_deriv( atom_j, atom_i, j_charge, i_charge, M_j );

				Vector force_vector_j = weights[ stack_elec ] * deriv_vector_j;

				if ( weights[ stack_elec_base_only ] != 0.0 && is_rna_base( rsd1, m ) ) force_vector_j += weights[ stack_elec_base_only ] * deriv_vector_j;

				F1 += cross( force_vector_j, atom_i );
				F2 += force_vector_j;
			}

		}
	}

}


//////////////////////////////////////////////////////////////////////////////////////////
//get_stack_elec_deriv evaluates the stack_elec_score between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by residue_pair_energy_one_way
Real
StackElecEnergy::get_stack_elec_score( Vector const & r_i,
																			 Vector const & r_j,
																			 Real const & i_charge,
																			 Real const & j_charge,
																			 Matrix const & M_i,
																			 Real & cos_kappa2 ) const
{

  Vector const z_i = M_i.col_z();

	Vector const r_vec = r_j - r_i;
	Real const r = r_vec.length();
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

  Real score  = coulomb().eval_atom_atom_hack_elecE( r_i, i_charge, r_j, j_charge);

  //Orientation dependence
	cos_kappa2 = cos_kappa * cos_kappa;

  score *= cos_kappa2;

  return score;

}


////////////////////////////////////////////////////////////////////////////////
//get_stack_elec_deriv evaluates the stack_elec_deriv between a pair of atoms.
//This function is called by eval_atom_derivative
Vector
StackElecEnergy::get_stack_elec_deriv( Vector const & r_i,
																			 Vector const & r_j,
																			 Real const & i_charge,
																			 Real const & j_charge,
																			 Matrix const & M_i ) const
{

	//  Well, the energy is a function of the inter-atom distance and a special cos(angle)
	//            E = E ( r, cos(theta ) ).
	//  so
	//      dE/dx = dE/dr (x/r)   +   ( - x * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dy = dE/dr (y/r)   +   ( - y * z / r^3 ) ( dE/dcos(theta) )
	//      dE/dz = dE/dr (z/r)   +   ( (r^2 - z^2) / r^3 ) ( dE/dcos(theta) )
	//

  Vector const x_i = M_i.col_x();
  Vector const y_i = M_i.col_y();
  Vector const z_i = M_i.col_z();

	Vector const r_vec = r_j - r_i;
	Real const r = r_vec.length();
	Real const x = dot( x_i, r_vec );
	Real const y = dot( y_i, r_vec );
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	/////////////////////////////////
	//dE_dcoskappa
	/////////////////////////////////
	Real dE_dcoskappa = coulomb().eval_atom_atom_hack_elecE( r_i, i_charge, r_j, j_charge);
  dE_dcoskappa  *= 2 * cos_kappa;


	/////////////////////////////////
	//dE_dr
	/////////////////////////////////
	Real const r2 = r_vec.length_squared();
  Real dE_dr_over_r  = coulomb().eval_dhack_elecE_dr_over_r( r2, i_charge, j_charge );
  //Orientation dependence
  dE_dr_over_r *= cos_kappa * cos_kappa ;

	Real const dE_dx = ( dE_dr_over_r ) * x - (dE_dcoskappa) * (x * z)/ (r * r * r) ;
	Real const dE_dy = ( dE_dr_over_r ) * y - (dE_dcoskappa) * (y * z)/ (r * r * r) ;
	Real const dE_dz = ( dE_dr_over_r ) * z + (dE_dcoskappa) * (x * x  +  y * y)/ (r * r * r);

  return (dE_dx * x_i + dE_dy * y_i + dE_dz * z_i );

}


//////////////////////////////////////////////////////////////////////////////////////////
void
StackElecEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.calculated() = false;

}

/// @brief StackElecEnergy distance cutoff
Distance
StackElecEnergy::atomic_interaction_cutoff() const
{
  return 0.0; /// Uh, I don't know.
}

core::Size
StackElecEnergy::version() const
{
       return 1; // Initial versioning
}


}
}
}
