// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNATorsionPotential.cc
/// @brief  DNATorsionPotential potential class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Rhiju Das
/// @author Jim Havranek

// Unit Headers
#include <core/scoring/dna/DNATorsionPotential.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/dna/DNA_Util.hh>

// Project Headers
#include <core/pose/Pose.hh>
//#include <core/io/database/open.hh>
//#include <basic/options/option.hh>
#include <basic/basic.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AmberPeriodicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>
//#include <ObjexxFCL/fmt/formatted.o.hh>

#include <basic/Tracer.hh>
#include <numeric/conversions.hh>
#include <iostream>
static basic::Tracer tr( "core.scoring.dna.DNATorsionPotential" );

namespace core {
namespace scoring {
namespace dna {

enum{ WHATEVER, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI, NU0, NU1, NU2, NU3, NU4 };

DNATorsionPotential::DNATorsionPotential():
	////////////////////////////////////////////////////////////////////////////
	// Parameters for alpha,beta,gamma,delta, etc. torsion constraints...
	////////////////////////////////////////////////////////////////////////////
	DELTA_CUTOFF_( 115.0 ),
//	scale_dna_torsion_tether_( 0.05 ), // THIS IS A SCALING FACTOR FOR ALL CONSTRAINTS.
	scale_dna_torsion_tether_( 0.05 ), // THIS IS A SCALING FACTOR FOR ALL CONSTRAINTS.
	scale_dna_torsion_sd_( 1.0 / std::sqrt( scale_dna_torsion_tether_ ) ),
	////////////////////////////////////////////////////////////////////////////
	// Ribose closure weights
	////////////////////////////////////////////////////////////////////////////
	c2star_c3star_bond_length_( 1.526 ),
	c2star_c3star_sd_( 1.0/ sqrt( 310.0 ) ), // 310.0 is the value of k
	c2star_c3star_dist_harm_func_( new constraints::HarmonicFunc( c2star_c3star_bond_length_, scale_dna_torsion_sd_ * c2star_c3star_sd_ )),

	c4star_c3star_c2star_bond_angle_( numeric::conversions::radians( 109.50 ) ),
	c4star_c3star_c2star_angle_harm_func_(
		new constraints::HarmonicFunc( c4star_c3star_c2star_bond_angle_, scale_dna_torsion_sd_ * 1.0/sqrt( numeric::conversions::radians( 40.0 ) ) ) ),

	o3star_c3star_c2star_bond_angle_( numeric::conversions::radians( 109.50 ) ),
	o3star_c3star_c2star_angle_harm_func_(
		new constraints::HarmonicFunc( o3star_c3star_c2star_bond_angle_, scale_dna_torsion_sd_ * 1.0/sqrt( numeric::conversions::radians( 50.0 ) ) ) ),

	c3star_c2star_c1star_bond_angle_( numeric::conversions::radians( 109.50 ) ),
	c3star_c2star_c1star_angle_harm_func_(
		new constraints::HarmonicFunc( c3star_c2star_c1star_bond_angle_, scale_dna_torsion_sd_ * 1.0/sqrt( numeric::conversions::radians( 40.0 ) ) ) ),

	// Might also be good to have additional angle or torsional potentials
	// to preserve sugar geometry.
	////////////////////////////////////////////////////////////////////////////
	verbose_( false )
{
	init_dna_torsion_parameters();
}


/////////////////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::setup_constraints(
   pose::Pose & pose,
	 constraints::ConstraintSetOP & dna_torsion_constraints,
	 constraints::ConstraintSetOP & dna_sugar_close_constraints,
	 constraints::ConstraintSetOP & dna_base_distance_constraints) const
{

	// dna_torsion_constraints->clear()  ...   doesn't exist!
	//Constraints are atom-pair, angle, dihedral...
	dna_sugar_close_constraints = constraints::ConstraintSetOP( new constraints::ConstraintSet );
	add_sugar_ring_closure_constraints( pose, *dna_sugar_close_constraints );
	//	add_o2star_torsion_constraints(     pose, *dna_torsion_constraints );

	// Why can't these terms be "constraints", in dna_torsion_constraints above? Because
	//  some involve atoms that change types when residues are switched in and out (during design!).
	// Could either define a different sort of constraint ("TorsionConstraint")
	dna_torsion_constraints = constraints::ConstraintSetOP( new constraints::ConstraintSet );
 	add_dna_torsion_tethers(  pose, *dna_torsion_constraints );

	dna_base_distance_constraints = constraints::ConstraintSetOP( new constraints::ConstraintSet );
	add_dna_base_distance_constraints( pose, *dna_base_distance_constraints );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::add_sugar_ring_closure_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const {
	for (Size i = 1; i <= pose.total_residue(); i++ ){
		add_sugar_ring_closure_constraints( pose.residue( i ), cst_set );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const {
	if ( !rsd.is_DNA() ) return;

	Size const & i( rsd.seqpos() );

	Size const c1star_index = rsd.atom_index( "C1'" );
	Size const c2star_index = rsd.atom_index( "C2'" );
	Size const c3star_index = rsd.atom_index( "C3'" );
	Size const o3star_index = rsd.atom_index( "O3'" );
	Size const c4star_index = rsd.atom_index( "C4'" );

	cst_set.add_constraint( new constraints::AtomPairConstraint( id::AtomID( c2star_index, i),
																															 id::AtomID( c3star_index, i),
																															 c2star_c3star_dist_harm_func_,
																															 dna_sugar_close ) );

	constraints::ConstraintOP angle1 = new constraints::AngleConstraint( id::AtomID( c4star_index, i),
																																			 id::AtomID( c3star_index, i),
																																			 id::AtomID( c2star_index, i),
																																			 c4star_c3star_c2star_angle_harm_func_,
																																			 dna_sugar_close );
	cst_set.add_constraint( angle1 );

	constraints::ConstraintOP angle2 = new constraints::AngleConstraint( id::AtomID( o3star_index, i),
																																			 id::AtomID( c3star_index, i),
																																			 id::AtomID( c2star_index, i),
																																			 o3star_c3star_c2star_angle_harm_func_,
																																			 dna_sugar_close );
	cst_set.add_constraint( angle2 );

	constraints::ConstraintOP angle3 = new constraints::AngleConstraint( id::AtomID( c3star_index, i),
																																			 id::AtomID( c2star_index, i),
																																			 id::AtomID( c1star_index, i),
																																			 c3star_c2star_c1star_angle_harm_func_,
																																			 dna_sugar_close );
	cst_set.add_constraint( angle3 );

	// Need to add an improper dihedral to keep the hydrogens correct on C2'



}
///////////////////////////////////////
void
DNATorsionPotential::add_dna_base_distance_constraints(
		pose::Pose & pose,
		constraints::ConstraintSet & cst_set ) const
{
	using namespace core::chemical;
	Size const nres = pose.total_residue();
	for( Size i = 1; i < nres; ++i ){
		conformation::Residue const & rsd( pose.residue( i ) );
		conformation::Residue const & next_rsd( pose.residue( i + 1 ) );
		if( !rsd.is_DNA() || !next_rsd.is_DNA() || rsd.is_upper_terminus() ) continue; //job undone: need to add conditions when the rsd is not basepaired

		Size const H2star_index = pose.residue( i ).atom_index( "H21*" );
		Size const H1star_index = pose.residue( i ).atom_index( " H2'" );
		Size H68_index, next_H68_index;
		if( rsd.type().aa() == na_ade || rsd.type().aa() == na_gua )
			H68_index = pose.residue( i ).atom_index( "H8" );
		else
			H68_index = pose.residue( i ).atom_index( "H6" );

		if( next_rsd.type().aa() == na_ade || next_rsd.type().aa() == na_gua )
			next_H68_index = pose.residue( i + 1 ).atom_index( "H8" );
		else
			next_H68_index = pose.residue( i + 1 ).atom_index( "H6" );
//set up harmonic function
		Real angle_diff ( rsd.mainchain_torsion(5) - rsd.mainchain_torsion(6) );
		Real dist_1H = 0.0041 * angle_diff + 2.7092;
		constraints::HarmonicFuncOP H1_harm_func(	new constraints::HarmonicFunc( dist_1H, 0.307 ) );
		cst_set.add_constraint( new constraints::AtomPairConstraint( id::AtomID( H1star_index, i),
																															 id::AtomID( next_H68_index, i + 1),
																															 H1_harm_func,
																															 dna_base_distance ) );

		Real dist_2H = 0.0081 * angle_diff + 4.0213;
		constraints::HarmonicFuncOP H2_harm_func(	new constraints::HarmonicFunc( dist_2H, 0.381 ) );
		cst_set.add_constraint( new constraints::AtomPairConstraint( id::AtomID( H2star_index, i),
																															 id::AtomID( next_H68_index, i + 1),
																															 H2_harm_func,
																															 dna_base_distance ) );

		Real dist_H68 = 0.0068 * angle_diff + 5.4228;
		constraints::HarmonicFuncOP H68_harm_func(	new constraints::HarmonicFunc( dist_H68, 0.373 ) );
		cst_set.add_constraint( new constraints::AtomPairConstraint( id::AtomID( H68_index, i),
																															 id::AtomID( next_H68_index, i + 1),
																															 H68_harm_func,
																															 dna_base_distance ) );
//	std::cout << "TEST" << " angle " << angle_diff << " dist_1H " << dist_1H << std::endl;
}
}

///////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::add_dna_torsion_tethers(
	 pose::Pose & pose,
	 constraints::ConstraintSet & cst_set ) const
{
	using namespace numeric;

	Size const nres = pose.total_residue();

	for (Size i = 1; i <=nres; i++ ){

		conformation::Residue const & rsd( pose.residue( i ) );
		if (!rsd.is_DNA() ) continue;

		/////////////////////////////////////
		// alpha --
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, ALPHA, alpha_components_ );

		/////////////////////////////////////
		// beta -- 1 harmonic potential
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, BETA, beta_components_ );

		/////////////////////////////////////
		// gamma
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, GAMMA, gamma_components_ );


		/////////////////////////////////////
		// delta
		/////////////////////////////////////

		add_DNA_torsion_constraint( pose, i, cst_set, DELTA, delta_components_ );

		/////////////////////////////////////
		// epsilon
		/////////////////////////////////////

		add_DNA_torsion_constraint( pose, i, cst_set, EPSILON, epsilon_components_ );

		/////////////////////////////////////
		// zeta
		/////////////////////////////////////

		add_DNA_torsion_constraint( pose, i, cst_set, ZETA, zeta_components_ );

		/////////////////////////////////////
		// nu0
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU0, nu0_components_ );

		/////////////////////////////////////
		// nu1
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU1, nu1_components_ );

		/////////////////////////////////////
		// nu2
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU2, nu2_components_ );

		/////////////////////////////////////
		// nu4
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU4, nu4_components_ );

#ifdef NOTDEF

		/////////////////////////////////////
		// nu3
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU3, nu3_components_ );

		/////////////////////////////////////
		// chi
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, CHI, gaussian_parameter_set_chi_south_);

		/////////////////////////////////////
		// nu2
		/////////////////////////////////////
		add_DNA_torsion_constraint( pose, i, cst_set, NU2, gaussian_parameter_set_nu2_south_);

#endif

	}
}


///////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::add_DNA_torsion_constraint(
			 pose::Pose & pose,
			 Size const i,
			 constraints::ConstraintSet & cst_set,
			 Size const dna_torsion_number,
			 utility::vector1< constraints::AmberPeriodicFuncOP > const & torsion_components ) const
{

	conformation::Residue rsd( pose.residue( i ) );

	// Get the atoms involved
	id::AtomID id1,id2,id3,id4;
	bool fail = get_atom_ids_by_torsion( dna_torsion_number, pose, i, id1, id2, id3, id4 );
	if( fail ) {
//		tr << "Failed to get atom ids at residue " << i << " for torsion number " << dna_torsion_number << std::endl;
		return;
	}

	// Generate dihedral constraints for each term in the vector of Fourier components

	for( Size this_comp = 1 ; this_comp <= torsion_components.size() ; ++this_comp ) {
//		tr << "Adding torsion at residue " << i << " ids " << id1 << " " << id2 << " " << id3 << " " << id4 << std::endl;
		constraints::ConstraintOP dihedral = new constraints::DihedralConstraint( id1, id2, id3, id4,
											torsion_components[ this_comp ], dna_bb_torsion );
		cst_set.add_constraint( dihedral );

	}

}

///////////////////////////////////////////////////////////////////////////////
void
DNATorsionPotential::init_dna_torsion_parameters()
{

  // Parameters for DNA backbone torsions are taken from the Amber MM code -
	// Alpha and gamma are from the refitting of Perez et.al. in Biophys. J. (2007)
	// v92 pp. 3816-3829.
	// The rest are taken from the parm99 parameter set available with the Ambertools
	// online.
	// It is my understanding that these should also work for RNA.
	// -jjh

	// Note there are scaling factors to move Amber parameters into our Amber format
	alpha_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians( 31.79508), 0.185181, 1.0 ) );
	alpha_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(351.95960), 1.256531, 2.0 ) );
	alpha_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(357.24748), 0.354858, 3.0 ) );

	beta_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.150000, 3.0 ) );

	gamma_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(190.97653), 1.178040, 1.0 ) );
	gamma_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(295.63279), 0.092102, 2.0 ) );
	gamma_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(348.09535), 0.962830, 3.0 ) );

	delta_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.400000, 3.0 ) );

	epsilon_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.150000, 3.0 ) );

	zeta_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.200000, 2.0 ) );
	zeta_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 0.250000, 3.0 ) );

	// CT-CT-OS-CT
	nu0_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(180.00000), 0.100000, 2.0 ) );
	nu0_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 0.383000, 3.0 ) );

	// CT-CT-CT-OS
	nu1_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.400000, 3.0 ) );

	// CT-CT-CT-CT
	nu2_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(180.00000), 0.200000, 1.0 ) );
	nu2_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(180.00000), 0.250000, 2.0 ) );
	nu2_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 0.180000, 3.0 ) );

	// CT-CT-CT-OS
	nu3_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 1.400000, 3.0 ) );

	// CT-CT-OS-CT
	nu4_components_.push_back( new constraints::AmberPeriodicFunc( numeric::conversions::radians(180.00000), 0.100000, 2.0 ) );
	nu4_components_.push_back( new constraints::AmberPeriodicFunc( 0.00000, 0.383000, 3.0 ) );

	// Now the relevant atom names

	alpha_atom_names_.push_back( "O3'" );
	alpha_atom_names_.push_back( "P" );
	alpha_atom_names_.push_back( "O5'" );
	alpha_atom_names_.push_back( "C5'" );

	beta_atom_names_.push_back( "P" );
	beta_atom_names_.push_back( "O5'" );
	beta_atom_names_.push_back( "C5'" );
	beta_atom_names_.push_back( "C4'" );

	gamma_atom_names_.push_back( "O5'" );
	gamma_atom_names_.push_back( "C5'" );
	gamma_atom_names_.push_back( "C4'" );
	gamma_atom_names_.push_back( "C3'" );

	delta_atom_names_.push_back( "C5'" );
	delta_atom_names_.push_back( "C4'" );
	delta_atom_names_.push_back( "C3'" );
	delta_atom_names_.push_back( "O3'" );

	epsilon_atom_names_.push_back( "C4'" );
	epsilon_atom_names_.push_back( "C3'" );
	epsilon_atom_names_.push_back( "O3'" );
	epsilon_atom_names_.push_back( "P" );

	zeta_atom_names_.push_back( "C3'" );
	zeta_atom_names_.push_back( "O3'" );
	zeta_atom_names_.push_back( "P" );
	zeta_atom_names_.push_back( "O5'" );

	nu0_atom_names_.push_back( "C4'" );
	nu0_atom_names_.push_back( "O4'" );
	nu0_atom_names_.push_back( "C1'" );
	nu0_atom_names_.push_back( "C2'" );

	nu1_atom_names_.push_back( "O4'" );
	nu1_atom_names_.push_back( "C1'" );
	nu1_atom_names_.push_back( "C2'" );
	nu1_atom_names_.push_back( "C3'" );

	nu2_atom_names_.push_back( "C1'" );
	nu2_atom_names_.push_back( "C2'" );
	nu2_atom_names_.push_back( "C3'" );
	nu2_atom_names_.push_back( "C4'" );

	nu3_atom_names_.push_back( "C2'" );
	nu3_atom_names_.push_back( "C3'" );
	nu3_atom_names_.push_back( "C4'" );
	nu3_atom_names_.push_back( "O4'" );

	nu4_atom_names_.push_back( "C3'" );
	nu4_atom_names_.push_back( "C4'" );
	nu4_atom_names_.push_back( "O4'" );
	nu4_atom_names_.push_back( "C1'" );

}

	bool
	DNATorsionPotential::get_atom_ids_by_torsion(
		Size const dna_torsion_number,
		pose::Pose & pose,
		Size const resid,
		id::AtomID & id1,
		id::AtomID & id2,
		id::AtomID & id3,
		id::AtomID & id4 ) const
{

	// Note:  A return value of 'true' denotes failure / not applicable!

	if( ( dna_torsion_number == ALPHA ) &&
			( pose.residue_type( resid ).is_lower_terminus() ) ) {
		return true;
	}

	if( ( dna_torsion_number == EPSILON ) &&
			( pose.residue_type( resid ).is_upper_terminus() ) ) {
		return true;
	}

	if( ( dna_torsion_number == ZETA ) &&
			( pose.residue_type( resid ).is_upper_terminus() ) ) {
		return true;
	}


	switch ( dna_torsion_number ) {
	case ALPHA:
		id1 = id::AtomID( pose.residue( resid - 1 ).atom_index( alpha_atom_names_[1] ), resid - 1 );
		id2 = id::AtomID( pose.residue( resid ).atom_index( alpha_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( alpha_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( alpha_atom_names_[4] ), resid );
		return false;
	case BETA:
		id1 = id::AtomID( pose.residue( resid ).atom_index( beta_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( beta_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( beta_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( beta_atom_names_[4] ), resid );
		return false;
	case GAMMA:
		id1 = id::AtomID( pose.residue( resid ).atom_index( gamma_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( gamma_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( gamma_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( gamma_atom_names_[4] ), resid );
		return false;
	case DELTA:
		id1 = id::AtomID( pose.residue( resid ).atom_index( delta_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( delta_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( delta_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( delta_atom_names_[4] ), resid );
		return false;
	case EPSILON:
		id1 = id::AtomID( pose.residue( resid ).atom_index( epsilon_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( epsilon_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( epsilon_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid + 1 ).atom_index( epsilon_atom_names_[4] ), resid + 1);
		return false;
	case ZETA:
		id1 = id::AtomID( pose.residue( resid ).atom_index( zeta_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( zeta_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid + 1 ).atom_index( zeta_atom_names_[3] ), resid + 1);
		id4 = id::AtomID( pose.residue( resid + 1 ).atom_index( zeta_atom_names_[4] ), resid + 1);
		return false;
	case NU0:
		id1 = id::AtomID( pose.residue( resid ).atom_index( nu0_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( nu0_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( nu0_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( nu0_atom_names_[4] ), resid );
		return false;
	case NU1:
		id1 = id::AtomID( pose.residue( resid ).atom_index( nu1_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( nu1_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( nu1_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( nu1_atom_names_[4] ), resid );
		return false;
	case NU2:
		id1 = id::AtomID( pose.residue( resid ).atom_index( nu2_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( nu2_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( nu2_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( nu2_atom_names_[4] ), resid );
		return false;
	case NU3:
		id1 = id::AtomID( pose.residue( resid ).atom_index( nu3_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( nu3_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( nu3_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( nu3_atom_names_[4] ), resid );
		return false;
	case NU4:
		id1 = id::AtomID( pose.residue( resid ).atom_index( nu4_atom_names_[1] ), resid );
		id2 = id::AtomID( pose.residue( resid ).atom_index( nu4_atom_names_[2] ), resid );
		id3 = id::AtomID( pose.residue( resid ).atom_index( nu4_atom_names_[3] ), resid );
		id4 = id::AtomID( pose.residue( resid ).atom_index( nu4_atom_names_[4] ), resid );
		return false;
	default:
		utility_exit_with_message("bad dna torsion type for DNATorsionPotential: " );
	}

	return true;
}






}
}
}
