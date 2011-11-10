// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_TorsionPotential.cc
/// @brief  RNA_TorsionPotential potential class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_TorsionPotential.hh>

// Package Headers
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rna/RNA_Util.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/basic.hh>

#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/CharmmPeriodicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/scoring/constraints/ConstantFunc.hh>
#include <core/scoring/constraints/ConstantFunc.fwd.hh>

#include <core/chemical/AtomType.hh>  //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009
#include <core/chemical/AtomTypeSet.hh>  //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType Oct 14, 2009


// Numeric Headers
#include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <numeric/interpolation/periodic_range/half/interpolation.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2A.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/TorsionID.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <numeric/conversions.hh>

//Add on Oct 29, 2009
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
//using namespace basic::options::OptionKeys;
//using namespace basic::options;


static basic::Tracer tr( "core.scoring.rna.RNA_TorsionPotential" );

namespace core {
namespace scoring {
namespace rna {

using namespace ObjexxFCL::fmt;

RNA_TorsionPotential::RNA_TorsionPotential():
	////////////////////////////////////////////////////////////////////////////
	// Parameters for alpha,beta,gamma,delta, etc. torsion constraints...
	////////////////////////////////////////////////////////////////////////////
	// "Tight torsion" means use the width of the major peak in the torsion
	// histogram to determine how strong the harmonic constraint will be. Otherwise, be less restrictive,
	// and use width of the minor fatter peak that lies under the sharp peak.
	rna_tight_torsions_( true ),
	DELTA_CUTOFF_( 115.0 ),
	scale_rna_torsion_tether_( 0.05 ), // THIS IS A SCALING FACTOR FOR ALL CONSTRAINTS.
	scale_rna_torsion_sd_( 1.0 / std::sqrt( scale_rna_torsion_tether_ ) ),
	////////////////////////////////////////////////////////////////////////////
	// Ribose closure weights
	////////////////////////////////////////////////////////////////////////////
	o4star_index_( 7  ),
	c1star_index_( 10 ),
	c2star_index_( 11 ),
	c4star_index_( 6  ),
	o4star_c1star_bond_length_( 1.414 ),
	o4star_c1star_sd_( 0.01 ),
	o4star_c1star_dist_harm_func_( new constraints::HarmonicFunc( o4star_c1star_bond_length_, scale_rna_torsion_sd_ * o4star_c1star_sd_ )),
	angle_sd_( numeric::conversions::radians( 1.0 ) ),
	o4star_c1star_c2star_bond_angle_( numeric::conversions::radians( 106.39 ) ),
	o4star_c1star_c2star_angle_harm_func_(
		new constraints::HarmonicFunc( o4star_c1star_c2star_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) ),
	o4star_c1star_first_base_bond_angle_( numeric::conversions::radians( 108.2 ) ),
	o4star_c1star_first_base_angle_harm_func_(
		new constraints::HarmonicFunc( o4star_c1star_first_base_bond_angle_, angle_sd_ ) ),
	c4star_o4star_c1star_bond_angle_( numeric::conversions::radians( 110.4 ) ),
	c4star_o4star_c1star_angle_harm_func_(
		new constraints::HarmonicFunc( c4star_o4star_c1star_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) ),
	// Might also be good to have additional angle or torsional potentials
	// to preserve sugar geometry.
	////////////////////////////////////////////////////////////////////////////
	// 2'-OH proton chi potential. A 3-well potential. Could also add
	// a weak potential to break degeneracy.
	////////////////////////////////////////////////////////////////////////////
	o2star_potential_weight_( 2.7 ),
	o2star_best_torsion_( numeric::conversions::radians( -140.0 ) ),
	o2star_dihedral_constraint_func1_( new constraints::CharmmPeriodicFunc( o2star_best_torsion_, o2star_potential_weight_, 3 ) ),
 	o2star_dihedral_constraint_func2_( new constraints::CharmmPeriodicFunc( o2star_best_torsion_, 0.01, 1 ) ),
	verbose_( false )
{
	init_rna_torsion_gaussian_parameters();
}


/////////////////////////////////////////////////////////////////////////////////////////
// void
// RNA_TorsionPotential::update_constraints( pose::Pose & pose ) const{


// 	//	std::cout << "UPDATE CONSTRAINTS " << std::endl;

// 	constraints::ConstraintSetOP constraint_set_with_rna_torsions( (pose.constraint_set())->clone() );

// 	add_rna_torsion_score_constraints(  pose, *constraint_set_with_rna_torsions );
// 	add_sugar_ring_closure_constraints( pose, *constraint_set_with_rna_torsions );
// 	add_o2star_torsion_constraints(     pose, *constraint_set_with_rna_torsions );

// 	pose.constraint_set( constraint_set_with_rna_torsions);


// }

/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionPotential::setup_constraints(
   pose::Pose & pose,
	 constraints::ConstraintSetOP & rna_torsion_constraints,
	 constraints::ConstraintSetOP & rna_sugar_close_constraints,
	 rna::RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers ) const
{
	// rna_torsion_constraints->clear()  ...   doesn't exist!

	//Constraints are atom-pair, angle, dihedral...
	rna_sugar_close_constraints = constraints::ConstraintSetOP( new constraints::ConstraintSet );
	add_sugar_ring_closure_constraints( pose, *rna_sugar_close_constraints );
	//	add_o2star_torsion_constraints(     pose, *rna_torsion_constraints );

	// Why can't these terms be "constraints", in rna_torsion_constraints above? Because
	//  some involve atoms that change types when residues are switched in and out (during design!).
	// Could either define a different sort of constraint ("TorsionConstraint")
	rna_torsion_constraints = constraints::ConstraintSetOP( new constraints::ConstraintSet );
	rna_side_chain_torsion_tethers.clear();
 	add_rna_torsion_tethers(  pose, *rna_torsion_constraints, rna_side_chain_torsion_tethers );
// 	std::cout << "Setting up RNA torsion constraints" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionPotential::add_sugar_ring_closure_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const {
	for (Size i = 1; i <= pose.total_residue(); i++ ){
		add_sugar_ring_closure_constraints( pose.residue( i ), cst_set );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionPotential::add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const {

	if ( !rsd.is_RNA() ) return;

	Size const & i( rsd.seqpos() );

	cst_set.add_constraint( new constraints::AtomPairConstraint( id::AtomID( o4star_index_, i),
																															 id::AtomID( c1star_index_, i),
																															 o4star_c1star_dist_harm_func_,
																															 rna_sugar_close ) );

	constraints::ConstraintOP angle1 = new constraints::AngleConstraint( id::AtomID( o4star_index_, i),
																																			 id::AtomID( c1star_index_, i),
																																			 id::AtomID( c2star_index_, i),
																																			 o4star_c1star_c2star_angle_harm_func_,
																																			 rna_sugar_close );
	cst_set.add_constraint( angle1 );

	constraints::ConstraintOP angle2 = new constraints::AngleConstraint( id::AtomID( c4star_index_, i),
																																			 id::AtomID( o4star_index_, i),
																																			 id::AtomID( c1star_index_, i),
																																			 c4star_o4star_c1star_angle_harm_func_,
																																			 rna_sugar_close );
	cst_set.add_constraint( angle2 );

	Size const first_base_index = first_base_atom_index( rsd );
	constraints::ConstraintOP angle3 = new constraints::AngleConstraint( id::AtomID( o4star_index_, i),
																																			 id::AtomID( c1star_index_, i),
																																			 id::AtomID( first_base_index, i),
																																			 o4star_c1star_first_base_angle_harm_func_,
																																			 rna_sugar_close );
	//	std::cout << "CST FIRST BASE ATOM " << rsd.name() << " " << rsd.seqpos() << " " << rsd.atom_name( first_base_index ) << std::endl;
	cst_set.add_constraint( angle3 );

}

///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionPotential::add_rna_torsion_tethers(
	 pose::Pose & pose,
	 constraints::ConstraintSet & cst_set,
	 rna::RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers ) const
{
	using namespace numeric;

	Size const nres = pose.total_residue();

	for (Size i = 1; i <=nres; i++ ){

		conformation::Residue const & rsd( pose.residue( i ) );
		if (!rsd.is_RNA() ) continue;

		/////////////////////////////////////
		// alpha --  3 harmonic potentials
		/////////////////////////////////////
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, ALPHA, 	gaussian_parameter_set_alpha_);

		/////////////////////////////////////
		// beta -- 1 harmonic potential
		/////////////////////////////////////
		add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, BETA, gaussian_parameter_set_beta_);

		/////////////////////////////////////
		// gamma
		/////////////////////////////////////
		add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, GAMMA, gaussian_parameter_set_gamma_);

		/////////////////////////////////////
		// delta
		/////////////////////////////////////
		Real const & delta= numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA) ); //Parin July 25, account for the possibility that the torsion will be outside the [-180:180] range
/* //Oct_19_2009
		if (delta <= DELTA_CUTOFF_) { // North, or 3'-endo sugar pucker, favored by RNA.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, DELTA, gaussian_parameter_set_delta_north_);
		} else { // South, or 2'-endo sugar pucker.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, DELTA, gaussian_parameter_set_delta_south_);
		}
*/
		/////////////////////////////////////
		// epsilon
		/////////////////////////////////////
			if (delta <= DELTA_CUTOFF_) { // North, or 3'-endo sugar pucker, favored by RNA.
				add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, EPSILON, 	gaussian_parameter_set_epsilon_north_);
			} else { // South, or 2'-endo sugar pucker.
				add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, EPSILON, 	gaussian_parameter_set_epsilon_south_);
			}

		/////////////////////////////////////
		// zeta
		/////////////////////////////////////
			Real next_alpha( -60.0 );
			if ( i < nres && pose.residue(i+1).is_RNA() ) next_alpha = numeric::principal_angle_degrees(pose.residue( i+1 ).mainchain_torsion( ALPHA ));
			//Parin July 25, account for the possibility that the torsion will be outside the [-180:180] range

			if ( next_alpha > -120.0 && next_alpha <= 0.0 ) { //default A-form, alpha sc-
				add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, ZETA, gaussian_parameter_set_zeta_alpha_sc_minus_);
			} else if ( next_alpha > 0.0 && next_alpha < 100.0 ) { // alpha sc+
				add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, ZETA, gaussian_parameter_set_zeta_alpha_sc_plus_);
			} else { // alpha ap
				add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, ZETA, gaussian_parameter_set_zeta_alpha_ap_);
			}


		/////////////////////////////////////
		// chi
		/////////////////////////////////////
		if (delta <= DELTA_CUTOFF_) { // North, or 3'-endo sugar pucker, favored by RNA.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, CHI, gaussian_parameter_set_chi_north_);
		} else { // South, or 2'-endo sugar pucker.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, CHI, gaussian_parameter_set_chi_south_);
		}

/*	//Oct_19_2009
		/////////////////////////////////////
		// nu2
		/////////////////////////////////////
		if (delta <= DELTA_CUTOFF_) { // North, or 3'-endo sugar pucker, favored by RNA.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, NU2, gaussian_parameter_set_nu2_north_);
		} else { // South, or 2'-endo sugar pucker.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, NU2, gaussian_parameter_set_nu2_south_);
		}

		/////////////////////////////////////
		// nu1
		/////////////////////////////////////
		if (delta <= DELTA_CUTOFF_) { // North, or 3'-endo sugar pucker, favored by RNA.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, NU1, gaussian_parameter_set_nu1_north_);
		} else { // South, or 2'-endo sugar pucker.
			add_RNA_torsion_constraint( pose, i, cst_set, rna_side_chain_torsion_tethers, NU1, gaussian_parameter_set_nu1_south_);
		}
*/
	}

}


///////////////////////////////////////////////////////////////////////////////

void
RNA_TorsionPotential::add_RNA_torsion_constraint(
			 pose::Pose & pose,
			 Size const i,
			 constraints::ConstraintSet & cst_set,
			 RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers,
			 Size const rna_torsion_number,
			 Gaussian_parameter_set const & gaussian_parameter_set ) const
{

	conformation::Residue rsd( pose.residue( i ) );

	//All this converting between torsion number <--> mainchain/chi is perhaps unnecessary?
	id::TorsionID  torsion_id( i, id::BB, rna_torsion_number );
	if ( rna_torsion_number > rna::NUM_RNA_MAINCHAIN_TORSIONS ) {
		Size const chino( rna_torsion_number - rna::NUM_RNA_MAINCHAIN_TORSIONS );
		torsion_id = id::TorsionID( i, id::CHI, chino  );
	}

	id::AtomID id1,id2,id3,id4;
	bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	if (fail) return;


/////////////////Parin July 26, virtual atom test/////////////////////////////////////////////////////////////////////

	//Does this create a hard copy? If so, should it be optimized?
	conformation::Residue const & rsd_1=pose.residue(id1.rsd());
	conformation::Residue const & rsd_2=pose.residue(id2.rsd());
	conformation::Residue const & rsd_3=pose.residue(id3.rsd());
	conformation::Residue const & rsd_4=pose.residue(id4.rsd());

	bool Is_virtual_torsion=( rsd_1.is_virtual(id1.atomno()) || rsd_2.is_virtual(id2.atomno()) || rsd_3.is_virtual(id3.atomno()) || rsd_4.is_virtual(id4.atomno()) );

	bool Is_virtual_torsion_old=( std::abs( rsd_1.atomic_charge( id1.atomno() ) ) < 1e-3 || std::abs( rsd_2.atomic_charge( id2.atomno() ) ) < 1e-3 || std::abs( rsd_3.atomic_charge( id3.atomno() ) ) < 1e-3 || std::abs( rsd_4.atomic_charge( id4.atomno() ) ) < 1e-3);

	if(Is_virtual_torsion!=Is_virtual_torsion_old){ //This is precaution...will take out after I am sure the new screen condition works Sep3, 2009 Parin
			tr.Info << "Error with Is_virtual_torsion_condition" << std::endl;
			exit (1);
	}

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Oct_18_2009, fix eplison/gamma/alpha bug at chain break closure
			if(Is_virtual_torsion){
				tr.Info << "Is_virtual_torsion";
			}else{
				tr.Info << "Is_NOT_virtual_torsion";
			}
			tr.Info << "  torsion_id: " << torsion_id;
			tr.Info << "  atom_id: " << id1 << " " << id2 << " " << id3 << " " << id4;
			tr.Info << "  name: " << rsd_1.type().atom_name(id1.atomno()) << " " << rsd_2.type().atom_name(id2.atomno()) << " " << rsd_3.type().atom_name(id3.atomno()) << " " << rsd_4.type().atom_name(id4.atomno());
			tr.Info << "  type: " << rsd_1.atom_type(id1.atomno()).name() << " " << rsd_2.atom_type(id2.atomno()).name() << " " << rsd_3.atom_type(id3.atomno()).name() << " " << rsd_4.atom_type(id4.atomno()).name() << std::endl;

*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	if(Is_virtual_torsion){
		if (verbose_){
			tr.Info << "In RNA_TorsionPotential_parin, encounter torsion containing one or more virtual atom(s)" << std::endl;
			tr.Info << "  torsion_id: " << torsion_id;
			tr.Info << "  atom_id: " << id1 << " " << id2 << " " << id3 << " " << id4 << std::endl;
			tr.Info << "  name: " << rsd_1.type().atom_name(id1.atomno()) << " " << rsd_2.type().atom_name(id2.atomno()) << " " << rsd_3.type().atom_name(id3.atomno()) << " " << rsd_4.type().atom_name(id4.atomno()) << std::endl;
			tr.Info << "  type: " << rsd_1.atom_type(id1.atomno()).name() << " " << rsd_2.atom_type(id2.atomno()).name() << " " << rsd_3.atom_type(id3.atomno()).name() << " " << rsd_4.atom_type(id4.atomno()).name() << std::endl;
			tr.Info << "		atom_type_index: " << rsd_1.atom_type_index( id1.atomno()) << " " << rsd_2.atom_type_index( id2.atomno()) << " " << rsd_3.atom_type_index( id3.atomno())  << " " << rsd_4.atom_type_index( id4.atomno()) << std::endl;
			tr.Info << "		atomic_charge: " << rsd_1.atomic_charge( id1.atomno())	<< " " << rsd_2.atomic_charge( id2.atomno())	<< " " << rsd_3.atomic_charge( id3.atomno())	<< " " << rsd_4.atomic_charge( id4.atomno()) << std::endl;
		}

		return;
	}


	//Check if atom id is virtual and return if it is. NEED TO IMPLEMENT THIS Parin July 24.

//	Real torsion_value( pose.torsion( torsion_id ) );

	Real torsion_value = numeric::principal_angle_degrees( pose.torsion( torsion_id ) ); //Parin July 25, account for the possibility that the torsion will be outside the [-180:180] range
	Real principal_torsion_value=torsion_value;

	// If there are multiple harmonic tethers available choose the "closest" one.
	assert( gaussian_parameter_set.size() > 0 );

	if ( torsion_id.type() == id::CHI ) {
		//Keep a list of all possible harmonics -- this is because during packing, chi's can switch between anti and syn.
		utility::vector1< constraints::CircularHarmonicFuncOP > harm_funcs;
		for (Size n = 1; n <= gaussian_parameter_set.size(); n++ ){

			Real center = gaussian_parameter_set[n].center;
			Real const width = scale_rna_torsion_sd_ * gaussian_parameter_set[n].width;


//			constraints::CircularHarmonicFuncOP harm_func  (new constraints::CircularHarmonicFunc( numeric::conversions::radians( center ), numeric::conversions::radians( width ) ) );


/////////////Update on Oct 29 ,2009. SECOND attempt at implementing syn_chi_penalty, this time will make the function continuous (will still contain discontinuities)
			Real offset;
			if(n==1){ //anti_chi
				offset=0.0;
			}else if(n==2){ //syn_chi
				offset=basic::options::option[ basic::options::OptionKeys::score::syn_chi_penalty ]/2.9; //Calibrated with RNA_torsion weight= 2.9
			}else{
				tr.Info << " Error in RNA_TorsionPotential::add_RNA_torsion_constraint, torsion_id.type() == id::CHI " << std::endl;
			}
			constraints::CircularHarmonicFuncOP harm_func= new constraints::CircularHarmonicFunc( numeric::conversions::radians( center ), numeric::conversions::radians( width ), offset );
			harm_funcs.push_back( harm_func );
		}

		rna_side_chain_torsion_tethers[i].push_back( std::make_pair( torsion_id, harm_funcs ) );
		//if (verbose_) tr.Info << "SIDE-CHAIN " <<
		//											I( 3,i) << " " << I(3, rna_torsion_number) << F(8,3,best_center)  << " " << F(8,4,torsion_value) << ": " <<
		//											F(8,3, harm_func->func( numeric::conversions::radians( torsion_value ) ) ) << std::endl;


/*	//////////////////////////////////Oct 28, 2009 (FIRST attempt at implementing syn_chi_penalty)///////////////////////////////////////////////////////
			Real syn_chi_penalty=basic::options::option[ basic::options::OptionKeys::score::syn_chi_penalty ]/2.9; //Calibrated with RNA_torsion weight= 2.9
			if (verbose_) tr.Info << "syn_chi_penalty= " << syn_chi_penalty << std::endl;

		//	gaussian_parameter_set_chi_north_.push_back( Gaussian_parameter(228.92, 79.43, 11.43) );
		//	gaussian_parameter_set_chi_north_.push_back( Gaussian_parameter( 1.07, -50.76, 26.11) );
		//	gaussian_parameter_set_chi_south_.push_back( Gaussian_parameter(12.53, 116.60, 25.23) );
		//	gaussian_parameter_set_chi_south_.push_back( Gaussian_parameter( 0.64, -49.91, 28.10) );
		//Define syn chi as [-146.655:14.715] where
		//14.715 is the average of 70.43 (north anti minimum) and -50.76 (north syn minimum)
		//-146.555 is the average of -243.4 (116.60-360, south anti minimum) and -49.91 (south syn minimum)
		//These are Rhiju's original parameters. Parin syn chi minima parameter's is slightly difference. See get_full_rotamers_include_syn_chi() of parin_rna_test.cc

		if( (principal_torsion_value>-146.655 && principal_torsion_value<14.715)==true){ //Test if chi is syn conformation.
			constraints::FuncOP const_func =new constraints::ConstantFunc( syn_chi_penalty);
			//DOES THIS WORK FOR SIDE CHAIN (CHI) AS WELL?
			constraints::ConstraintOP dihedral = new constraints::DihedralConstraint( id1, id2, id3, id4, const_func, rna_torsion ); //Oct_19_2009 Where the hell is rna_torsion defined?
			cst_set.add_constraint( dihedral );
		}
*/	/////////////////////////////////////////////////////////////////////////////////////////////////////


	}	else {

		Real best_center( 0.0 ), best_width( 0.0 ), best_weight( 0.0 ), best_sigma2( 10000000.0 ), best_deviation( 10000000.0 );

		for (Size n = 1; n <= gaussian_parameter_set.size(); n++ ){

			Real center = gaussian_parameter_set[n].center;
			//scale_rna_torsion_sd_=4.472
			Real const width = scale_rna_torsion_sd_ * gaussian_parameter_set[n].width;
			Real const weight = 1.0/( width * width );
			Real deviation = numeric::principal_angle_degrees( torsion_value - center );
			Real const sigma2 = deviation * deviation * weight;

			if (sigma2 < best_sigma2 ){
				best_center = center;
				best_weight = weight;
				best_deviation = deviation;
				best_sigma2 = sigma2;
				best_width = width;
			}

		}

		assert( best_weight > 0.0 );

		// Perhaps this would be better as an "AmbiguousConstraint"? The only problem is that it will take longer to compute...
		constraints::FuncOP func  (new constraints::CircularHarmonicFunc( numeric::conversions::radians( best_center ), numeric::conversions::radians( best_width ) ) );


		if( torsion_id.torsion() == 5)  { //epsilon
			Real score=best_sigma2;
			Real max_score=1; //Calibrated with RNA_torsion weight= 2.9

			if(score> max_score){ //Impose constant maximum
				func = new constraints::ConstantFunc( max_score );
			}

		} else if( torsion_id.torsion() == 1 )  { //alpha

			if( (principal_torsion_value>-64 && principal_torsion_value<66)==false) {
				Real score=best_sigma2;
				Real max_score=0.25; //Calibrated with RNA_torsion weight= 2.9

				if(score> max_score){ //Impose constant maximum
					func = new constraints::ConstantFunc( max_score );
				}
			}

		} else if( torsion_id.torsion() == 6) { //zeta

			if( (gaussian_parameter_set[1].amplitude) > 0.9 && (gaussian_parameter_set[1].amplitude) < 1.1) { //Hacky thing to indicate that alpha is gauche_minus

				if((principal_torsion_value<-70)) {
						Real score=best_sigma2;
						Real max_score=0.172; //0.5 divide by 2.9
						if(score> max_score){ //Impose constant maximum
							func = new constraints::ConstantFunc( max_score );
						}
				} else if((principal_torsion_value>12.5 && principal_torsion_value<gaussian_parameter_set[2].center)) {
						Real score=best_sigma2;
						Real min_score=0.172; //0.5 divide by 2.9
						if(score < min_score){ //Impose constant maximum
							func = new constraints::ConstantFunc( min_score );
						}
				} else if(principal_torsion_value>=gaussian_parameter_set[2].center-2){ //Need the -2 so that potential is define exactly at mean[2]
							Real min_score=0.172; //0.5 divide by 2.9
							func = new constraints::ConstantFunc( min_score );
				}
			} else if( (gaussian_parameter_set[1].amplitude) > 1.9 && (gaussian_parameter_set[1].amplitude) < 2.1) { //Hacky thing to indicate that alpha is gauche_plus

				Real score=best_sigma2;
				Real min_score=0.172; //0.5 divide by 2.9

				if(score < min_score || principal_torsion_value<=gaussian_parameter_set[1].center || principal_torsion_value>=gaussian_parameter_set[2].center){ //Impose constant maximum
					func = new constraints::ConstantFunc( min_score );
				}

			} else if( (gaussian_parameter_set[1].amplitude) > 2.9 && (gaussian_parameter_set[1].amplitude) < 3.1) { //Hacky thing to indicate that alpha is trans
				if((principal_torsion_value<-70) || (principal_torsion_value>63)) {
						Real score=best_sigma2;
						Real max_score=0.3448; //1 divide by 2.9

						if(score> max_score){ //Impose constant maximum
							func = new constraints::ConstantFunc( max_score );
						}
				}
			}


		}



		////////////////////////////////////////////////////////////////////////////
		// This is handled by constraint machinery ...
		constraints::ConstraintOP dihedral = new constraints::DihedralConstraint( id1, id2, id3, id4, func, rna_torsion ); //Oct_19_2009 Where the hell is rna_torsion defined?
		cst_set.add_constraint( dihedral );
//		if (verbose_) tr.Info << "BB " << I( 3,i) << " " << I(3, rna_torsion_number) << F(8,3,best_center)  << " " << F(8,4,torsion_value) << ": " << F(8,3,dihedral->score( pose.conformation() ) ) << std::endl;

//Parin Apr 12, 2009
					if (verbose_) tr.Info << "BB " <<
										I( 3,i) << " " << I(3, rna_torsion_number) <<  " "<< F(8,4,torsion_value)  << " b_center=" << F(8,3,best_center)  << " b_width=" << F(8,3,best_width)  << " b_sigma2=" << F(8,3,best_sigma2)  << ": " <<
										F(8,3,dihedral->score( pose.conformation() ) ) << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Down with DOF constraints!!!
	// Same hacky thing as in Rosetta --- need to figure out offset of DOF value to torsion...
	////////////////////////////////////////////////////////////////////////////
	//	id::DOF_ID const dof_id(  pose.conformation().dof_id_from_torsion_id( torsion_id ) );

	//Real const torsion_offset = pose.conformation().dof( dof_id ) /*radians*/ - numeric::conversions::radians( torsion_value );

	//This could be done in a smarter way by making the Circular Harmonic Func's ahead of time.
	//	constraints::CircularHarmonicFuncOP harm_func  (new constraints::CircularHarmonicFunc( numeric::conversions::radians( best_center ) + torsion_offset,
	//																																													 numeric::conversions::radians( best_width ) ) );

	//	cst_set.add_dof_constraint( dof_id, harm_func, rna_torsion );


}

///////////////////////////////////////////////////////////////////////////////
// Should this stuff related to the 2'-OH proton be a different Energy Method and Potential?
void
RNA_TorsionPotential::add_o2star_torsion_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const
{
	using numeric::conversions::radians;
	using namespace constraints;

	for (Size i = 1; i<= pose.total_residue(); i++ ){

		if ( !pose.residue(i).is_RNA() ) continue;

	// 	///////////////////////////////////////////////////////////
// 		// Assumes chi torsion 4 is the 2'-OH torsion.
// 		id::TorsionID torsion_id = id::TorsionID( i, id::CHI, 4 );

// 		// A little check.
// 		{
// 			id::AtomID id1,id2,id3,id4;
// 			bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
// 			assert( !fail && pose.residue( id4.rsd() ).atom_name( id4.atomno() ) == "2HO*" );
// 		}

// 		id::DOF_ID const dof_id(  pose.conformation().dof_id_from_torsion_id( torsion_id ) );
// 		Real const torsion_offset = pose.conformation().dof( dof_id ) /*radians*/ - radians( pose.torsion( torsion_id) );

// 		FuncOP o2star_constraint_func1( new CharmmPeriodicFunc( o2star_best_torsion_ +torsion_offset, o2star_potential_weight_, 3 ) );
// 		FuncOP o2star_constraint_func2( new CharmmPeriodicFunc( o2star_best_torsion_ +torsion_offset, 0.01, 1 ) );

// 		//////////////////////////////////////////////////////////////////////////////
// 		//////////////////////////////////////////////////////////////////////////////
// 		// NOTE THAT BECAUSE THESE ARE DOF CONSTRAINTS, they won't actually
// 		// be calculated in the packer -- just the minimizer
// 		//////////////////////////////////////////////////////////////////////////////
// 		//////////////////////////////////////////////////////////////////////////////

// 		cst_set.add_dof_constraint( dof_id, o2star_constraint_func1, rna_torsion );
// 		cst_set.add_dof_constraint( dof_id, o2star_constraint_func2, rna_torsion );

		//Replace above with dihedral constraint, which will be calculation by packer.

 		id::TorsionID torsion_id = id::TorsionID( i, id::CHI, 4 );

		id::AtomID id1,id2,id3,id4;
		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		assert( !fail );
		assert( pose.residue( id4.rsd() ).atom_name( id4.atomno() ) == "2HO*" );

		constraints::ConstraintOP dihedral1 =
			new constraints::DihedralConstraint( id1, id2, id3, id4,
																					 o2star_dihedral_constraint_func1_,
																					 rna_torsion );
		cst_set.add_constraint( dihedral1 );

		constraints::ConstraintOP dihedral2 =
			new constraints::DihedralConstraint( id1, id2, id3, id4,
																					 o2star_dihedral_constraint_func2_,
																					 rna_torsion );
		cst_set.add_constraint( dihedral2 );

	}


}


///////////////////////////////////////////////////////////////////////////////
//WARNING!!: July 24th, the amplitude was not recalculated

void
RNA_TorsionPotential::init_rna_torsion_gaussian_parameters()
{
	// These numbers refer to amplitude, mean, and width of fitted Gaussians.

/*
	if ( rna_tight_torsions_ ) {
//		gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(222.03, -64.11,  9.64) );  By Parin July 17
		gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(222.03, -64.11,  20.0) ); By Parin July 24
	} else {
		gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(13.12, -73.74,  35.97) );
	}
	gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(10.51, 66.01,  18.09) );
	gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(18.40, 161.80,  18.12) );

	if ( rna_tight_torsions_ ) {
//		gaussian_parameter_set_beta_.push_back( Gaussian_parameter(181.33, 176.33,  11.54) );  By Parin July 17
		gaussian_parameter_set_beta_.push_back( Gaussian_parameter(181.33, 176.33,  30.0) ); //First trail change to 20....now make it wider to 45(too weak) change to 30.0
	} else {
		gaussian_parameter_set_beta_.push_back( Gaussian_parameter(32.30, 174.52,  43.56) );
	}


	if ( rna_tight_torsions_ ) {
//		gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(366.90, 53.08,  6.64) );  By Parin July 17
			gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(366.90, 53.08,  20.0) );
	} else {
		gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(18.26, 56.59,  20.57) );
	}
	gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(21.61, 178.19,  13.61) );
	gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(3.98, -64.02,  17.76) );

	gaussian_parameter_set_epsilon_north_.push_back( Gaussian_parameter(178.08, -150.17,  14.64) );
	gaussian_parameter_set_epsilon_north_.push_back( Gaussian_parameter(2.52, 68.28,  32.29) );
	gaussian_parameter_set_epsilon_south_.push_back( Gaussian_parameter(11.95, -98.45, 26.80) );
	gaussian_parameter_set_epsilon_south_.push_back( Gaussian_parameter( 0.58, 159.70, 103.86) );

	if ( rna_tight_torsions_ ) {
//		gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter( 143.97, -71.45, 7.91) ); By Parin July 17
//		gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter( 143.97, -71.45, 20.0) ); By Parin July 25
			gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter( 143.97, -71.45, 15.28) ); //Steal std from zeta_trans_alpha
	} else {
		gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter( 78.74, -68.60, 16.19) );
	}
	//	gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter(  2.43, 178.84, 114.82) );
	gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter(  2.43, 178.84, 40.00) );

	gaussian_parameter_set_zeta_alpha_sc_plus_.push_back( Gaussian_parameter(2.08, -137.28,  63.12) );
	gaussian_parameter_set_zeta_alpha_sc_plus_.push_back( Gaussian_parameter(3.37, 87.07,  32.69) );

	gaussian_parameter_set_zeta_alpha_ap_.push_back( Gaussian_parameter(13.65, -69.74,  15.28) );
	gaussian_parameter_set_zeta_alpha_ap_.push_back( Gaussian_parameter(2.69, 63.03,  33.61) );
*/

	gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(222.03, -69.51,  20.36) ); //Aform
	gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(10.51, 67.38,  18.12) );
	gaussian_parameter_set_alpha_.push_back( Gaussian_parameter(18.40, 168.3,  23.28) );


	gaussian_parameter_set_beta_.push_back( Gaussian_parameter(181.33, 179.5,  37.4) );

	gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(366.90, 54.20,  9.74) ); //Aform
	gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(3.98, -68.32	,  13.98) );
	gaussian_parameter_set_gamma_.push_back( Gaussian_parameter(21.61, 175.65,  12.82) );

	gaussian_parameter_set_epsilon_north_.push_back( Gaussian_parameter(178.08, -144.34,  25.02) );
	gaussian_parameter_set_epsilon_south_.push_back( Gaussian_parameter(11.95, -108.87, 25.31) );

//DO NOT CHANGE THE ORDER of push back for zeta, delta, chi, nu1 and nu2......many part of the code assumes this particular order!!!!

	////////////////////mean and std for zeta is choosen to get the desire potential shape///////////////////////////////
	//Design this so that there will be a maxima centered at zeta=12.5 degree,
	gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter( 1.0, -71.45, 22.5) );
	gaussian_parameter_set_zeta_alpha_sc_minus_.push_back( Gaussian_parameter(  1.0, 96, 22.5) );

	//Design this so that there will be a maxima centered at zeta=-25 degree,
	gaussian_parameter_set_zeta_alpha_sc_plus_.push_back( Gaussian_parameter(2.0, -125.00 ,  26.93) );
	gaussian_parameter_set_zeta_alpha_sc_plus_.push_back( Gaussian_parameter(2.0, 75.00,  26.93) );

	gaussian_parameter_set_zeta_alpha_ap_.push_back( Gaussian_parameter(3.0, -69.74,  15.28) );
	gaussian_parameter_set_zeta_alpha_ap_.push_back( Gaussian_parameter(3.0, 63.03,  33.61) );

//////For delta, chi, nu1 and nu2, use Rhiju's old parameter//////////////////////////////////////////////////////////////////////////

	gaussian_parameter_set_delta_north_.push_back( Gaussian_parameter(687.92, 82.90,  3.99) );
	gaussian_parameter_set_delta_south_.push_back( Gaussian_parameter(53.18, 145.25,  6.35) );

	gaussian_parameter_set_chi_north_.push_back( Gaussian_parameter(228.92, 79.43, 11.43) );
	gaussian_parameter_set_chi_north_.push_back( Gaussian_parameter( 1.07, -50.76, 26.11) );

	gaussian_parameter_set_chi_south_.push_back( Gaussian_parameter(12.53, 116.60, 25.23) );
	gaussian_parameter_set_chi_south_.push_back( Gaussian_parameter( 0.64, -49.91, 28.10) );

	gaussian_parameter_set_nu2_north_.push_back( Gaussian_parameter(1148.01, 38.82, -2.77) );
	gaussian_parameter_set_nu2_south_.push_back( Gaussian_parameter(173.15, -37.22, 3.00) );

	gaussian_parameter_set_nu1_north_.push_back( Gaussian_parameter(631.56, 95.34,  4.20) );
	gaussian_parameter_set_nu1_south_.push_back( Gaussian_parameter(57.04, 155.51,  6.00) );

}



}
}
}
