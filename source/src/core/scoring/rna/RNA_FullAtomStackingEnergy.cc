// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_FullAtomStacking.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_FullAtomStackingEnergy.hh>
#include <core/scoring/rna/RNA_FullAtomStackingEnergyCreator.hh>

// Package headers
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/kinematics/MinimizerMapBase.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>


// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>



// C++

static basic::Tracer tr("core.scoring.rna.RNA_FullAtomStackingEnergy");
using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_FullAtomStackingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_FullAtomStackingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_FullAtomStackingEnergy;
}

ScoreTypes
RNA_FullAtomStackingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_stack );
	sts.push_back( fa_stack_aro );
	return sts;
}


/// c-tor
RNA_FullAtomStackingEnergy::RNA_FullAtomStackingEnergy() :
	parent( new RNA_FullAtomStackingEnergyCreator ),
	//Parameters are totally arbitrary and made up!
	prefactor_ ( 0.2 ),
	full_stack_cutoff_ ( 4.0 ), //Encompass next base stack
	dist_cutoff_ ( 6.0 ), // Do not to next-nearest base stack!
	dist_cutoff2_ ( dist_cutoff_ * dist_cutoff_ ),
	base_base_only_( true )
{}

//clone
methods::EnergyMethodOP
RNA_FullAtomStackingEnergy::clone() const
{
  return new RNA_FullAtomStackingEnergy;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::setup_for_minimizing(
    pose::Pose & pose,
    ScoreFunction const & sfxn,
    kinematics::MinimizerMapBase const & min_map
) const
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    
    //set_nres_mono(pose);
    
    if ( pose.energies().use_nblist() ) {
        // stash our nblist inside the pose's energies object
        Energies & energies( pose.energies() );
        
        // setup the atom-atom nblist
        NeighborListOP nblist;
        Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
        Real const XX = dist_cutoff_ + 2 * tolerated_motion;
        nblist = new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX);
        if ( pose.energies().use_nblist_auto_update() ) {
            nblist->set_auto_update( tolerated_motion );
        }
        // this partially becomes the EtableEnergy classes's responsibility
        nblist->setup( pose, sfxn, *this);
        energies.set_nblist( EnergiesCacheableDataType::FA_STACK_NBLIST, nblist );
    }
}

///
void
RNA_FullAtomStackingEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn) const
{
  pose.update_residue_neighbors();

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.update( pose );

	pose.update_residue_neighbors();
	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::FA_STACK_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

void
RNA_FullAtomStackingEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & scfxn) const
{
  pose.update_residue_neighbors();

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.update( pose );

}

//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

  if ( !rsd1.is_RNA() ) return;
  if ( !rsd2.is_RNA() ) return;

	Real score_aro1( 0.0 ), score_aro2( 0.0 );

	Real const score = residue_pair_energy_one_way( rsd1, rsd2, pose, score_aro1 ) +
		residue_pair_energy_one_way( rsd2, rsd1, pose, score_aro2 ) ;

  emap[ fa_stack ]       += score;

  emap[ fa_stack_aro ]   += score_aro1 + score_aro2;

	if ( score <= -0.0001 ) {
		//		tr.Info << rsd1.name3()  << rsd1.seqpos() << "---" << rsd2.name3() << rsd2.seqpos() << ": " << score << std::endl;
	}

	//tr.Info << rsd1.name3()  << rsd1.seqpos() << "---" << rsd2.name3() << rsd2.seqpos() << ": " << score << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_FullAtomStackingEnergy::residue_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & score_aro
) const
{

	score_aro = 0.0;

  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

  Real score( 0.0 );

  Size const i( rsd1.seqpos() );
	//  Size const j( rsd2.seqpos() );

	//  Vector const & centroid_i( base_centroids[i] );
  kinematics::Stub const & stub_i( base_stubs[i] );
  Matrix const M_i ( stub_i.M );

  // Loop over base heavy atoms.
  // If I want to generalize this to proteins, maybe could loop over "aromatic" atoms.
  for ( Size m = rsd1.first_sidechain_atom(); m <= rsd1.nheavyatoms(); ++m ) {
	//Need to be careful nheavyatoms count includes hydrogen! when the hydrogen is made virtual..

		if(rsd1.is_virtual(m)) continue;

		if(m==rsd1.first_sidechain_atom()){
			//Consistency check
			if(rsd1.type().atom_name(m) !=" O2'") utility_exit_with_message( "m==rsd1.first_sidechain_atom() but rsd1.type().atom_name(m) !=\" O2'\" ");
			continue;
		}

    Vector const heavy_atom_i( rsd1.xyz( m ) );

		//Look for occlusion by other base atoms? Or all other heavy atoms?
		Size const atom_num_start = base_base_only_ ? rsd2.first_sidechain_atom() : 1 ;

    for ( Size n = atom_num_start; n <= rsd2.nheavyatoms(); ++n ) {

			if(rsd2.is_virtual(n)) continue;

			if(n==rsd2.first_sidechain_atom() && base_base_only_){
				//Consistency check
				if(rsd2.type().atom_name(n) !=" O2'") utility_exit_with_message( "n==rsd2.first_sidechain_atom() but rsd2.type().atom_name(n) !=\" O2'\" ");
				continue;
			}




      Vector const heavy_atom_j( rsd2.xyz( n ) );
      Vector r = heavy_atom_j - heavy_atom_i;
      Real const dist2 = r.length_squared();

      if ( dist2 < dist_cutoff2_ ) {

				//				Distance const dist = sqrt( dist2 );
				Real const fa_stack_score = get_fa_stack_score( r, M_i );
				score += fa_stack_score;

				if ( is_aro( rsd1, m) && is_aro( rsd2, n) ) score_aro += fa_stack_score;

				//				Real const cos_kappa_j = dot( r, z_j);
				//				score += get_score( dist, cos_kappa_j );

				//DEBUG
//				std::cout << "ONE  " << "res_name= " << rsd1.name() << " res_seqpos= " << rsd1.seqpos();
//				std::cout << " atom " << m  << " " << 	"name= " << rsd1.type().atom_name(m) << " type= " << rsd1.atom_type(m).name();

//				std::cout << "TWO  " << "res_name= " << rsd2.name() << " res_seqpos= " << rsd2.seqpos();
//				std::cout << " atom " << n  << " " << 	"name= " << rsd2.type().atom_name(n) << " type= " << rsd2.atom_type(n).name() << std::endl;

      }
    }
  }


	return score;
}


//////////////////////////////////////////////////////////////////////////////
// Depending on the option "base_base_only_" allow a stacking interaction
//  between an atom in a nucleobase and othe nucleobase atoms or *any*
//  other atom.
//
//if ( !check_base_base_OK( rsd1, rsd2, m, n ) && !check_base_base_OK( rsd2, rsd1, n, m ) ) continue;
//NEED BOTH TO BE NOT OK TO CONTINUE;
//rsd1 with m, rsd2 with n
//m needs to be a side chain atom no matter what
//IN base_base_only_mode, n needs to be a side chain atom

//So if base_base_only, then only consider interaction between base atoms of each base.
//if base_base_only==false, then consider interaction bewteen base atom of one base and all other atoms of the other base.

bool
RNA_FullAtomStackingEnergy::check_base_base_OK(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
  Size const & m, Size const & n ) const {

	if ( m < rsd1.first_sidechain_atom()  || m > rsd1.nheavyatoms() ) return false;

	if(m==rsd1.first_sidechain_atom()){
		//Consistency check
		if(rsd1.type().atom_name(m) !=" O2'") utility_exit_with_message( "m==rsd1.first_sidechain_atom() but rsd1.type().atom_name(m) !=\" O2'\" ");
		return false;
	}

	if ( n > rsd2.nheavyatoms() ) return false;

	if (base_base_only_ && n < rsd2.first_sidechain_atom() ) return false;

	if(n==rsd2.first_sidechain_atom() && base_base_only_){
		//Consistency check
		if(rsd2.type().atom_name(n) !=" O2'") utility_exit_with_message( "n==rsd2.first_sidechain_atom() but rsd2.type().atom_name(n) !=\" O2'\" ");
		return false;
	}

	return true;

}


//////////////////////////////////////////////////////////////////////////////
bool
RNA_FullAtomStackingEnergy::is_aro(
	conformation::Residue const & rsd1,
  Size const & m ) const {

	return ( rsd1.atom_type( m ).is_aromatic() );
	return true;

}

////////////////////////////////////////////////////////////////////////////// (Need to condition this? Parin Sep 2, 2009)
void
RNA_FullAtomStackingEnergy::eval_atom_derivative(
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

	if (rsd1.is_virtual(m)) return;
	if ( !rsd1.is_RNA() ) return;

	if ( m > rsd1.nheavyatoms() ) return;

  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	//  Vector const & centroid_i( base_centroids[i] );
  kinematics::Stub const & stub_i( base_stubs[i] );
  Matrix const M_i ( stub_i.M );

	Vector const heavy_atom_i( rsd1.xyz( m ) );

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

		if ( !rsd2.is_RNA() ) continue;

		//Look for occlusion by other base atoms? Or all other heavy atoms?
		//		Size const atom_num_start = base_base_only_ ? rsd2.first_sidechain_atom() : 1 ; //check_base_base_OK take care of this part.

		kinematics::Stub const & stub_j( base_stubs[j] );
		Matrix const M_j ( stub_j.M );

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {


			if(rsd2.is_virtual(n)) continue;

			if ( !check_base_base_OK( rsd1, rsd2, m, n ) && !check_base_base_OK( rsd2, rsd1, n, m ) ) continue; //This screen for nonbase atoms

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector r = heavy_atom_j - heavy_atom_i;
			Real const dist2 = r.length_squared();

			if ( dist2 < dist_cutoff2_ ) { //dist_cutoff2=dist_cutoff*dist_cutoff. Energy fade from max to 0 between  full_stack_cutoff_ and dist_cutoff_

				if ( check_base_base_OK( rsd1, rsd2, m, n ) ) {
					Vector force_vector_i = weights[ fa_stack ] * get_fa_stack_deriv( r, M_i );

					if ( is_aro( rsd1, m ) && is_aro( rsd2, n ) ) force_vector_i += weights[ fa_stack_aro ] * get_fa_stack_deriv( r, M_i );

					//Force/torque with which occluding atom j acts on "dipole" i.
					F1 += -1.0 * cross( force_vector_i, heavy_atom_j );
					F2 += -1.0 * force_vector_i;
				}

				if ( check_base_base_OK( rsd2, rsd1, n, m ) ) {
					//Force/torque with which occluding atom i acts on "dipole" j.
					// Note that this calculation is a repeat of some other call to this function.
					// Might make more sense to do some (alternative) bookkeeping.
					Vector force_vector_j = weights[ fa_stack ] * get_fa_stack_deriv( -r, M_j );

					if ( is_aro( rsd1, m ) && is_aro( rsd2, n ) ) force_vector_j += weights[ fa_stack_aro ] * get_fa_stack_deriv( -r, M_j );

					F1 += cross( force_vector_j, heavy_atom_i );
					F2 += force_vector_j;
				}

			}
		}
	}

}


//////////////////////////////////////////////////////////////////////////////////////////
//get_fa_stack_deriv evaluates the fa_stack_score between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by residue_pair_energy_one_way
Real
RNA_FullAtomStackingEnergy::get_fa_stack_score( Vector const r_vec, Matrix const M_i ) const
{

  Vector const z_i = M_i.col_z();

	Real const r = r_vec.length();
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

  Real score  = -1.0 * prefactor_;

  //Orientation dependence
  score *= cos_kappa * cos_kappa ;

  Distance b = r - full_stack_cutoff_;

  //Just use a simple cubic spline to fade from 1 to 0,
  // and to force derivatives to be continuous.
  if ( b > 0.0 ) {
    b /= (dist_cutoff_ - full_stack_cutoff_ );
    Real const b2 = b*b;
    Real const b3 = b2*b;
    score *= ( 2 * b3 - 3 * b2 + 1 );
  }

  return score;

}


//////////////////////////////////////////////////////////////////////////////////////////
//get_fa_stack_deriv evaluates the fa_stack_deriv between a pair of atoms. (r_vec is the vector between them, M_i is the coordinate matrix of one of the base)
//This function is called by eval_atom_derivative
Vector
RNA_FullAtomStackingEnergy::get_fa_stack_deriv( Vector const r_vec, Matrix const M_i ) const
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

	Real const r = r_vec.length();
	Real const x = dot( x_i, r_vec );
	Real const y = dot( y_i, r_vec );
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

	/////////////////////////////////
	//dE_dcoskappa
	/////////////////////////////////
  Real dE_dcoskappa  = -1.0 * prefactor_;

  //Orientation dependence
  dE_dcoskappa *= 2 * cos_kappa ;

  Distance b = r - full_stack_cutoff_;
	Distance distance_scale = (dist_cutoff_ - full_stack_cutoff_ );
	b /= distance_scale;

	Real b2( 0.0 ), b3( 0.0 );
  //Just use a simple cubic spline to fade from 1 to 0,
  // and to force derivatives to be continuous.
  if ( b > 0.0 && b < 1.0) {
    b2 = b*b;
    b3 = b2*b;
    dE_dcoskappa *= ( 2 * b3 - 3 * b2 + 1 );
  }

	/////////////////////////////////
	//dE_dr
	/////////////////////////////////
  Real dE_dr  = -1.0 * prefactor_;

  //Orientation dependence
  dE_dr *= cos_kappa * cos_kappa ;

  //Just use a simple cubic spline to fade from 1 to 0,
  // and to force derivatives to be continuous.
  if ( b > 0.0 && b < 1.0) {
    dE_dr *= ( 6 * b2 - 6 * b );
		dE_dr /= distance_scale;
  } else {
		dE_dr = 0.0;
	}

	Real const dE_dx = ( dE_dr ) * (x/r) - (dE_dcoskappa) * (x * z)/ (r*r*r) ;
	Real const dE_dy = ( dE_dr ) * (y/r) - (dE_dcoskappa) * (y * z)/ (r*r*r) ;
	Real const dE_dz = ( dE_dr ) * (z/r) + (dE_dcoskappa) * (x*x + y*y)/ (r*r*r);

	//	std::cout << " HELLO " << r << " " << dE_dr << " " << dE_dcoskappa << std::endl;

  return (dE_dx * x_i + dE_dy * y_i + dE_dz * z_i );

}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomStackingEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{

  rna::RNA_ScoringInfo  & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
  rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
  rna_centroid_info.calculated() = false;

}

/// @brief RNA_FullAtomStackingEnergy distance cutoff
Distance
RNA_FullAtomStackingEnergy::atomic_interaction_cutoff() const
{
  return dist_cutoff_ + 2 * chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH; //Similar to FA_ElecEnergy
}

core::Size
RNA_FullAtomStackingEnergy::version() const
{
       return 1; // Initial versioning
}

etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_intrares_countpair(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	utility_exit_with_message( "FA_ElecEnergy does not define intra-residue pair energies; do not call get_intrares_countpair()" );
	return 0;
}

etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	if ( res1 == res2 ) {
		return new CountPairNone;
	}

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}


etable::count_pair::CountPairFunctionCOP
RNA_FullAtomStackingEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace etable::count_pair;

	if ( !rsd1.is_RNA() ) return new CountPairNone;
	if ( !rsd2.is_RNA() ) return new CountPairNone;
	if ( rsd1.seqpos() == rsd2.seqpos() ) return new CountPairNone;
	return new CountPairAll;

}


}
}
}
