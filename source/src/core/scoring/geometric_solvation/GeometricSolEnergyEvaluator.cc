// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/chemical/rna/util.hh>

// Project headers
#include <core/chemical/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <ObjexxFCL/format.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

static basic::Tracer tr( "core.scoring.geometric_solvation.GeometricSolEnergyEvaluator" );

using namespace core::scoring::hbonds;
using namespace ObjexxFCL::format;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// TO DO --
// could allow cached calculation of derivatives.
///////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

///@brief copy c-tor
GeometricSolEnergyEvaluator::GeometricSolEnergyEvaluator( methods::EnergyMethodOptions const & opts ) :
	options_( new methods::EnergyMethodOptions( opts ) ),
	hb_database_( HBondDatabase::get_database( opts.hbond_options().params_database_tag() )),
	dist_cut2_( 27.0 ),   // 5.2*5.2
	geometric_sol_scale_( 0.4 * 1.17 / 0.65 ),
	correct_geom_sol_acceptor_base_( basic::options::option[ basic::options::OptionKeys::score::geom_sol_correct_acceptor_base ]() ),
	verbose_( false )
{
	//////////////////////////////
	// What's up with following?
	//////////////////////////////
	options_->exclude_DNA_DNA( false 	/*GEOMETRIC SOLVATION NOT COMPATIBLE WITH EXCLUDE_DNA_DNA FLAG YET*/ );
	options_->hbond_options().use_incorrect_deriv( true );  // override command line.
	options_->hbond_options().use_sp2_chi_penalty( false ); // override command line.
}

/// copy ctor
GeometricSolEnergyEvaluator::GeometricSolEnergyEvaluator( GeometricSolEnergyEvaluator const & src ) :
		utility::pointer::ReferenceCount(src),
		options_( new methods::EnergyMethodOptions( * src.options_ )),
		hb_database_( src.hb_database_ ),
		dist_cut2_( src.dist_cut2_ ),   // 5.2*5.2
		geometric_sol_scale_( src.geometric_sol_scale_ ),
		correct_geom_sol_acceptor_base_( src.correct_geom_sol_acceptor_base_ ),
		verbose_( src.verbose_ )
{}

//Destructor
GeometricSolEnergyEvaluator::~GeometricSolEnergyEvaluator()
{}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// Everything in here.
void
GeometricSolEnergyEvaluator::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap ) const
{


	if ( rsd1.seqpos() == rsd2.seqpos() ) return;  //Is this necessary?
	//	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	//////////////////////////////////
	// Yo, what about count_pair?
	//  Need to consider when including
	//  intra terms
	//////////////////////////////////

  Real geo_solE =
    res_res_geometric_sol_one_way( rsd1, rsd2, pose ) +
    res_res_geometric_sol_one_way( rsd2, rsd1, pose ) ;

  // store the energies
    emap[ geom_sol ] += geo_solE;



}


/////////////////////////////////////////////////////////////////////
// This is meant to be a reasonable clone of John Karanicolas'
//  Rosetta++ code. Some of the crazy options to do backbone only
//  or side chain only are not ported over, but they could
//  be implemented using the "atom_is_backbone" information.
//
// [Note to self (rhiju) : It might also make sense to precompute and cache this data
// in, say, a geometric solvation potential object, so that derivatives
// don't need to be computed over and over again, and code won't be
// copied. ]
//
Real
GeometricSolEnergyEvaluator::res_res_geometric_sol_one_way(
	conformation::Residue const & polar_rsd,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose ) const
{

	Real geo_solE =
		donorRes_occludingRes_geometric_sol_one_way(    polar_rsd, occ_rsd, pose ) +
		acceptorRes_occludingRes_geometric_sol_one_way( polar_rsd, occ_rsd, pose );

	return geo_solE;
}


//////////////////////////////////////////////////////////////////////////////////////
//optimization functions
//Added by Joseph Yesselman 9/5/2013
//////////////////////////////////////////////////////////////////////////////////////


Real
GeometricSolEnergyEvaluator::geometric_sol_one_way_sc(
                                                         conformation::Residue const & polar_rsd,
                                                         conformation::Residue const & occ_rsd,
                                                         pose::Pose const & pose ) const
{

  Real geo_solE =
  donorRes_occludingRes_geometric_sol_one_way_sc( polar_rsd, occ_rsd, pose ) +
  acceptorRes_occludingRes_geometric_sol_one_way_sc( polar_rsd, occ_rsd, pose );

  return geo_solE;


}


inline
Real
GeometricSolEnergyEvaluator::acceptorRes_occludingRes_geometric_sol_one_way_sc(
		conformation::Residue const & acc_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose) const
{
	Real res_solE( 0.0 ), energy( 0.0 );

	for ( chemical::AtomIndices::const_iterator
			anum  = acc_rsd.accpt_pos().begin(),
			anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atm( *anum );

		Size start (1);

		if(acc_rsd.atom_is_backbone(acc_atm)) { start = occ_rsd.first_sidechain_atom(); }

		for ( Size occ_atm = start; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			//Important NOTE. I originally had the code in the following function
			// written out inside this loop -- and packing was faster.
			// Perhaps something to do with inlining or compiler optimization.
			// I've left it this way for now, because it helps prevent copying too
			// much of the code shared between  residue pair scoring and for the derivatives.
			// However, if speed becomes important, here's a place to start.
			get_atom_atom_geometric_solvation_for_acceptor( acc_atm, acc_rsd, occ_atm, occ_rsd, pose, energy);
			res_solE += energy;
		}
	}

	return res_solE;
}

inline
Real
GeometricSolEnergyEvaluator::donorRes_occludingRes_geometric_sol_one_way_sc(
                                                                               conformation::Residue const & don_rsd,
                                                                               conformation::Residue const & occ_rsd,
                                                                               pose::Pose const & pose) const
{
  Real res_solE( 0.0 ), energy( 0.0 );

  // Here we go -- cycle through polar hydrogens in don_aa, everything heavy in occluding atom.
  for ( chemical::AtomIndices::const_iterator
       hnum  = don_rsd.Hpos_polar().begin(),
       hnume = don_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
    Size const don_h_atm( *hnum );


    Size start (1);

    if(don_rsd.atom_is_backbone(don_h_atm)) { start = occ_rsd.first_sidechain_atom(); }

    // if(assign_rna_atom_type(don_rsd, don_h_atm) != type_1) continue;

    for ( Size occ_atm = start; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

      //Important NOTE. I originally had the code in the following function
      // written out inside this loop -- and packing was faster.
      // Perhaps something to do with inlining or compiler optimization.
      // I've left it this way for now, because it helps prevent copying too
      // much of the code shared between  residue pair
      // scoring and for the derivatives.
      // However, if speed becomes important, here's a place to start.
      get_atom_atom_geometric_solvation_for_donor( don_h_atm, don_rsd,
                                                  occ_atm, occ_rsd, pose, energy );
      res_solE += energy;
    }
  }

  return res_solE;
}




///////////////////////////////////////////////////////////////////////////////////////

Real
GeometricSolEnergyEvaluator::geometric_sol_one_way_bb_bb(
  conformation::Residue const & polar_rsd,
  conformation::Residue const & occ_rsd,
  pose::Pose const & pose ) const
{

  Real geo_solE =
    donorRes_occludingRes_geometric_sol_one_way_bb_bb( polar_rsd, occ_rsd, pose ) +
    acceptorRes_occludingRes_geometric_sol_one_way_bb_bb( polar_rsd, occ_rsd, pose );

	return geo_solE;


}


Real
GeometricSolEnergyEvaluator::acceptorRes_occludingRes_geometric_sol_one_way_bb_bb(
  conformation::Residue const & acc_rsd,
  conformation::Residue const & occ_rsd,
  pose::Pose const & pose) const
{

  Real res_solE( 0.0 ), energy( 0.0 );

  for ( chemical::AtomIndices::const_iterator
       anum  = acc_rsd.accpt_pos().begin(),
       anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {

    Size const acc_atm( *anum );

    if(!acc_rsd.atom_is_backbone(acc_atm)) continue;


    //if(acc_atm > acc_rsd_last_backbone_num) continue;

    for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

      if(!occ_rsd.atom_is_backbone(occ_atm)) continue;

      //if(occ_atm > occ_rsd_last_backbone_num) continue;

      //Important NOTE. I originally had the code in the following function
      // written out inside this loop -- and packing was faster.
      // Perhaps something to do with inlining or compiler optimization.
      // I've left it this way for now, because it helps prevent copying too
      // much of the code shared between  residue pair
      // scoring and for the derivatives.
      // However, if speed becomes important, here's a place to start.
      get_atom_atom_geometric_solvation_for_acceptor( acc_atm, acc_rsd,
                                                     occ_atm, occ_rsd, pose, energy);
      res_solE += energy;
    }
  }

  return res_solE;

}


Real
GeometricSolEnergyEvaluator::donorRes_occludingRes_geometric_sol_one_way_bb_bb(
  conformation::Residue const & don_rsd,
  conformation::Residue const & occ_rsd,
  pose::Pose const & pose) const
{
	Real res_solE( 0.0 ), energy( 0.0 );

  // Here we go -- cycle through polar hydrogens in don_aa, everything heavy in occluding atom.
  for ( chemical::AtomIndices::const_iterator
       hnum  = don_rsd.Hpos_polar().begin(),
       hnume = don_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
    Size const don_h_atm( *hnum );

    //if(don_h_atm > don_rsd_last_backbone_num) continue;

    if(!don_rsd.atom_is_backbone(don_h_atm)) continue;


    for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

      if(!occ_rsd.atom_is_backbone(occ_atm)) continue;

      //if(occ_atm > occ_rsd_last_backbone_num) continue;

      //Important NOTE. I originally had the code in the following function
      // written out inside this loop -- and packing was faster.
      // Perhaps something to do with inlining or compiler optimization.
      // I've left it this way for now, because it helps prevent copying too
      // much of the code shared between  residue pair
      // scoring and for the derivatives.
      // However, if speed becomes important, here's a place to start.
      get_atom_atom_geometric_solvation_for_donor( don_h_atm, don_rsd,
                                                  occ_atm, occ_rsd, pose, energy );
      res_solE += energy;
    }
  }

  return res_solE;
}










//////////////////////////////////////////////////////////////////////////////////////
inline
Real
GeometricSolEnergyEvaluator::donorRes_occludingRes_geometric_sol_one_way(
	conformation::Residue const & don_rsd,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose ) const
{

	Real res_solE( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through polar hydrogens in don_aa, everything heavy in occluding atom.
	for ( chemical::AtomIndices::const_iterator
			hnum  = don_rsd.Hpos_polar().begin(),
			hnume = don_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			//Important NOTE. I originally had the code in the following function
			// written out inside this loop -- and packing was faster.
			// Perhaps something to do with inlining or compiler optimization.
			// I've left it this way for now, because it helps prevent copying too
			// much of the code shared between  residue pair
			// scoring and for the derivatives.
			// However, if speed becomes important, here's a place to start.
			get_atom_atom_geometric_solvation_for_donor( don_h_atm, don_rsd,
				occ_atm, occ_rsd, pose, energy );
			res_solE += energy;
		}
	}

	return res_solE;
}

///////////////////////////////////////////////////////////////////////////////////////
inline
Real
GeometricSolEnergyEvaluator::acceptorRes_occludingRes_geometric_sol_one_way(
	conformation::Residue const & acc_rsd,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose ) const
{

	Real res_solE( 0.0 ), energy( 0.0 );

	for ( chemical::AtomIndices::const_iterator
					anum  = acc_rsd.accpt_pos().begin(),
					anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {

		Size const acc_atm( *anum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			//Important NOTE. I originally had the code in the following function
			// written out inside this loop -- and packing was faster.
			// Perhaps something to do with inlining or compiler optimization.
			// I've left it this way for now, because it helps prevent copying too
			// much of the code shared between  residue pair
			// scoring and for the derivatives.
			// However, if speed becomes important, here's a place to start.
			get_atom_atom_geometric_solvation_for_acceptor( acc_atm, acc_rsd,
				occ_atm, occ_rsd, pose, energy);
			res_solE += energy;
		}
	}

	return res_solE;
}


//////////////////////////////////////////////////////////////////////////////
// Helper function for creating a mock H or O for the mock water.
// Assume it is in the same plane as that defined by the other
// atoms involved in the hydrogen bond.
// There are probably (better) examples of this function elsewhere in the code.
//
inline
void
GeometricSolEnergyEvaluator::set_water_base_atm(
  Vector const & base_v,
	Vector const & atom_v,
	Vector const & water_v,
	Vector & water_base_v,
	Real const & xH /*cos(theta)*/,
	Distance const & bond_length
) const
{
	Vector x,y,z, direction;

	//Define coordinate system.
	z = cross( (water_v - atom_v),  (atom_v - base_v) );
	z.normalize();
	y = water_v - atom_v;
	y.normalize();
	x = cross( y, z );

	// Plop the atom down
	direction = xH * y  + Real( std::sqrt( 1 - (xH * xH) ) ) * x;
	water_base_v = water_v + bond_length * direction;

}


/////////////////////////////////////////////////////////////////////////////////////////
inline
Real
GeometricSolEnergyEvaluator::occluded_water_hbond_penalty(
  bool const & is_donor,
	hbonds::HBEvalTuple const & hbond_eval_type,
	Vector const & polar_atm_xyz,
	Vector const & base_atm_xyz,
	Vector const & occluding_atm_xyz,
	Size const & polar_nb,
	Size const & occ_nb,
	bool const update_deriv /* = false*/,
	HBondDerivs & deriv /* = DUMMY_DERIV2D*/ ) const
{

	static bool const always_do_full_calculation( true );

	// jumpout criteria copied from hb_energy_deriv in hbonds.cc

	// Compute geometry
	Real AHdis( 0.0 ), xD( 0.0 ), xH( 0.0 ), energy( 0.0 );

	//Craziness... create an artifical atom to complete "water molecule".
	Vector water_base_atm;
	static Distance const water_O_H_distance( 0.958 );
	Real environment_weight( 1.0 );

	//Might be cleaner to separate this into two functions, one for donor, one for acceptor?
	if ( is_donor ) {

		// water is the acceptor, give it perfect geometry
		xH = 1./3.;  // perfect geometry is cos( 180 - 109.5 degrees), which is 1/3

		// compute the distance to the accepting water
		Real const AHdis2 = (polar_atm_xyz - occluding_atm_xyz).length_squared();
		if ( AHdis2 > MAX_R2 ) return 0.0;
		if ( AHdis2 < MIN_R2 ) return 0.0;
		AHdis = std::sqrt(AHdis2);

		// find the cosine of the base-proton-water angle (xD)
		xD = get_water_cos( base_atm_xyz, polar_atm_xyz, occluding_atm_xyz );
		if ( xD < MIN_xD ) return 0.0;
		if ( xD > MAX_xD ) return 0.0;

		//rhiju, testing alternative calculation that will give derivative.
		Vector occluding_base_atm_xyz( 0.0 );
		if ( update_deriv || always_do_full_calculation ) {
			set_water_base_atm( base_atm_xyz, polar_atm_xyz, occluding_atm_xyz, occluding_base_atm_xyz,
				xH, water_O_H_distance );
			hb_energy_deriv( *hb_database_, options_->hbond_options(),
				hbond_eval_type, base_atm_xyz, polar_atm_xyz,
				occluding_atm_xyz,
				occluding_base_atm_xyz,
				occluding_base_atm_xyz,
				energy, hbderiv_ABE_GO_NO_xH, deriv);
			if (verbose_) tr << "DERIV ENERGY DONOR:  " << energy;
		}

		if (options_->hbond_options().use_hb_env_dep() ) {
			environment_weight = hbonds::get_environment_dependent_weight( hbond_eval_type, polar_nb, occ_nb, options_->hbond_options() );
		}

	} else {

		// water is the donor, give it perfect geometry
		xD = 0.9999;

		// compute the distance to the accepting water proton
		// subtract the water's OH distance to get the AHdis,
		// since the distance computed was from the acceptor to the water oxygen
		// note: water proton lies on the line between the acceptor and the water oxygen
		AHdis = ( polar_atm_xyz - occluding_atm_xyz ).length();
		AHdis -= water_O_H_distance; // water O-H distance
		Real const AHdis2 = AHdis * AHdis;
		if ( AHdis2 > MAX_R2 ) return 0.;
		if ( AHdis2 < MIN_R2 ) return 0.;

		// find cosine of the base-acceptor-water_proton angle (xH)
		// note: this is the same as the base-acceptor-water_oxygen angle
		xH = get_water_cos( base_atm_xyz, polar_atm_xyz, occluding_atm_xyz );
		if ( xH < MIN_xH ) return 0.;
		if ( xH > MAX_xH ) return 0.;

		//rhiju, testing alternative calculation that will give derivative.
		if ( update_deriv || always_do_full_calculation ) {
			Vector occluding_base_atm_xyz( 0.0 );
 			set_water_base_atm( base_atm_xyz, polar_atm_xyz, occluding_atm_xyz, occluding_base_atm_xyz,
				-xD, water_O_H_distance );
 			hb_energy_deriv( *hb_database_, options_->hbond_options(), hbond_eval_type,
				occluding_atm_xyz, occluding_base_atm_xyz,
				polar_atm_xyz,	base_atm_xyz, base_atm_xyz,
				energy, hbderiv_ABE_GO_NO_xD, deriv);
			if (verbose_) tr << "DERIV ENERGY ACCPT: " << energy;

			if (options_->hbond_options().use_hb_env_dep() ) {
				environment_weight = hbonds::get_environment_dependent_weight( hbond_eval_type, polar_nb, occ_nb, options_->hbond_options() );
			}

		}


	}

	// Note that following should be a little faster and could be used if derivative is not necessary
	// However, use of hb_energy_deriv gives nearly exact match of analytical and numerical.
	// while every once in a while this does not...
	if ( !always_do_full_calculation ) {
		Real dummy_chi( 0.0 );
		assert( ! options_->hbond_options().use_sp2_chi_penalty() ); // APL avoid the new sp2 chi term.
		hbond_compute_energy(
			*hb_database_, options_->hbond_options(), hbond_eval_type,
			AHdis, xD, xH, dummy_chi, energy );
	}

	if (verbose_ ) tr << "  jk ENERGY: " << energy << std::endl;

	core::Real sol_penalty = -1.0 * energy;

	// jk THIS NEEDS TO BE FIT MORE RIGOROUSLY LATER...
	// Apply a scaling factor (effectively a weight), tying the weight of the
	// solvation term to the Hbond term (rather than to the LK weight)
	// Note: chose the bb-sc Hbond weight, because they're all about the same
	//	core::Real const sol_weight = geometric_sol_weight * pack_wts.Whbond_bb_sc() / pack_wts.Wsol();
	Real reweight = geometric_sol_scale_ * environment_weight;

	sol_penalty *= reweight;

	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//   WHERE DOES THE ENVIRONMENT WEIGHT COME IN?
	//  -- wasn't taken into account in rosetta++ --
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	if ( update_deriv ){
		if ( is_donor ) {
			deriv.h_deriv.f1()   *= -1.0 * reweight;
			deriv.h_deriv.f2()   *= -1.0 * reweight;
			deriv.don_deriv.f1() *= -1.0 * reweight;
			deriv.don_deriv.f2() *= -1.0 * reweight;
		} else {
			// had to flip reference frame to get HB energy and derivative.
			deriv.h_deriv.f1()   *= +1.0 * reweight;
			deriv.h_deriv.f2()   *= +1.0 * reweight;
			deriv.don_deriv.f1() *= +1.0 * reweight;
			deriv.don_deriv.f2() *= +1.0 * reweight;
		}

	}

	// this is a penalty, don't return a negative number
	if (sol_penalty < 0.) {
		if ( update_deriv ) deriv = ZERO_DERIV2D;
		return 0.;
	}

	//	if (sol_penalty > 0. ) {
	//		tr << F(7,3,polar_atm_xyz(1)) << " "  << F(7,3,base_atm_xyz(1)) <<  " " << F(7,3,occluding_atm_xyz(1) ) << std::endl;
		//		tr << "[AHdis " << F(7,3,AHdis) << "; xD " << F(7,3,xD) << ";  xH " << F(7,3,xH) << ", e " << F(7,3,energy)  <<  "]" ;
	//	}

	return sol_penalty; // return a positive number (this is a penalty)

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inline
void
GeometricSolEnergyEvaluator::get_atom_atom_geometric_solvation_for_donor(
	Size const & don_h_atm,
	conformation::Residue const & don_rsd,
	Size const  & occ_atm,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose,
	Real & energy,
	bool const update_deriv /*= false*/,
	HBondDerivs & deriv /* = DUMMY_DERIVS */
) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	//Why do we need to send in the pose and the residue stuff?
	// Well, the pose has info on backbone H-bonds.
	// and, during design, the residue type doesn't actually
	// have to match what's in the pose! Tricky!

	// In case of early return, initialize. Note that energy *does not* accumulate
	energy = 0.0;
	deriv = ZERO_DERIV2D;

	if ( occ_rsd.is_virtual( occ_atm ) ) return;
	if ( don_rsd.is_virtual( don_h_atm ) ) return;

	//Need to know about backbone/backbone H-bonds for proteins
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( HBOND_SET )));
	TenANeighborGraph const & tenA_neighbor_graph
		( pose.energies().tenA_neighbor_graph() );

	assert( atom_is_donor_h( don_rsd, don_h_atm ) );

	Size const don_base_atm( don_rsd.atom_base( don_h_atm ) );

	Vector const & don_h_atm_xyz( don_rsd.atom( don_h_atm ).xyz() );
	Vector const & don_base_atm_xyz( don_rsd.atom( don_base_atm ).xyz() );

	// the base atom isn't allowed to occlude solvent
	// Note: In current implementation, intraresidue pairs aren't checked...
	if ( ( don_rsd.seqpos() == occ_rsd.seqpos() ) && ( occ_atm == don_base_atm ) ) return;
	if (occ_rsd.is_protein() && occ_rsd.is_virtual(occ_atm)) return;
	if (don_rsd.is_protein() && don_rsd.is_virtual(don_h_atm)) return;

	// if a backbone donor participates in a backbone-backbone Hbond,
	// nothing is allowed to occlude solvent except a backbone acceptor
	bool const don_h_atm_is_protein_backbone
		( don_rsd.is_protein() && don_rsd.atom_is_backbone( don_h_atm ) );
	bool const occ_atm_is_protein_backbone_acceptor
		( occ_rsd.is_protein() && occ_rsd.atom_is_backbone( occ_atm ) && atom_is_acceptor(occ_rsd, occ_atm ) );
	if ( don_h_atm_is_protein_backbone && hbond_set.don_bbg_in_bb_bb_hbond( don_rsd.seqpos() ) &&
			 !occ_atm_is_protein_backbone_acceptor ) {
		return;
	}

	// if the distance is > 5.2 A, from the base atom, it doesn't occlude solvent
	Vector const & occ_atm_xyz( occ_rsd.atom( occ_atm ).xyz() );
	Real const base_dis2 = ( occ_atm_xyz - don_base_atm_xyz).length_squared();
	if ( base_dis2 > dist_cut2_ ) return;

	// if distance to base atom is greater than distance to donor, it doesn't occlude solvent
	Real const hdis2 = ( occ_atm_xyz - don_h_atm_xyz ).length_squared();
	if ( hdis2 > base_dis2 ) return;

	// For a backbone CO occluding a backbone NH, use the backbone-backbone (linear) geometry
	// to compute solvation penalty (only really matters for the CO acceptor's geometry, but
	// do it here as well for consistency)
	bool const potential_backbone_backbone_hbond =
		( don_h_atm_is_protein_backbone && occ_atm_is_protein_backbone_acceptor );
	HBEvalTuple hbe = potential_backbone_backbone_hbond ? ( HBEvalTuple( don_base_atm, don_rsd, occ_atm, occ_rsd) ) : HBEvalTuple( hbdon_HXL, hbacc_HXL, seq_sep_other ); // apl note: the false condition maps to hbe_dHXLaHXL

	Size don_nbrs( 0 ), occ_nbrs( 0 );
	if (options_->hbond_options().use_hb_env_dep() ) {
		don_nbrs = tenA_neighbor_graph.get_node( don_rsd.seqpos() )->num_neighbors_counting_self();
		occ_nbrs = tenA_neighbor_graph.get_node( occ_rsd.seqpos() )->num_neighbors_counting_self();
	}

	// jk Compute the Hbond energy as if this was a water
	// jk Add the water Hbond energy to a running total sum, as well as to the residue sum
	// mjo not sure what hbe should be if potential_backbone_backbone_hbond is false,
	// mjo I'm setting it to arbitrary type that previously was hbe_SP3SC just for consistency.
	// mjo Previously it ignores sequence separation.  why?
	energy = occluded_water_hbond_penalty(
		true /*is_donor*/, hbe,
		don_h_atm_xyz, don_base_atm_xyz, occ_atm_xyz,
		don_nbrs, occ_nbrs,
		update_deriv, deriv);

	if ( verbose_ && ( energy > 0.0 ) ) {
		tr <<"jk DON res "<< don_rsd.name1() << I(3,don_rsd.seqpos())<<
			" atom "<< don_rsd.atom_name( don_h_atm )<<" is occluded by occ_res " <<
			occ_rsd.name1()<< I(3, occ_rsd.seqpos()) <<
			" atom "<< occ_rsd.atom_name( occ_atm ) <<
			"  (HBEvalType " <<  I(2,hbe.eval_type()) << ") " <<
			" with energy "<< F(8,3,energy)<< std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Vector
GeometricSolEnergyEvaluator::get_acceptor_base_atm_xyz( conformation::Residue const & acc_rsd, Size const & acc_atm ) const{

	Vector base_atm_xyz;

	if ( correct_geom_sol_acceptor_base_ ) {
		// following handles the special case in which acceptor has two base residues -- occurs for
		//  N inside rings.
		// This matches machinery in hbonds_geom.cc.  That doesn't mean that the base atom is set
		//  totally correctly -- water (TIP3.params) and O4' in nucleic acids still have weird base atoms,
		//  but at least the hbonds and geom_sol match up.
		Vector dummy;
		chemical::Hybridization acc_hybrid( acc_rsd.atom_type( acc_atm ).hybridization());
		make_hbBasetoAcc_unitvector(
			options_->hbond_options(),
			acc_hybrid,
			acc_rsd.atom( acc_atm ).xyz(),
			acc_rsd.xyz( acc_rsd.atom_base( acc_atm ) ),
			acc_rsd.xyz( acc_rsd.abase2( acc_atm ) ),
			base_atm_xyz, dummy );
	} else {
		base_atm_xyz = acc_rsd.atom( acc_rsd.atom_base( acc_atm ) ).xyz();
	}
	return base_atm_xyz;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
void
GeometricSolEnergyEvaluator::get_atom_atom_geometric_solvation_for_acceptor(
	Size const & acc_atm,
	conformation::Residue const & acc_rsd,
	Size const & occ_atm,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose,
	Real & energy,
	bool const update_deriv /*= false*/,
	HBondDerivs & deriv /* = DUMMY_DERIV2D */) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	//Why do we need to send in the pose and the residue stuff?
	// Well, the pose has info on backbone H-bonds.
	// and, during design, the residue type doesn't actually
	// have to match what's in the pose! Tricky!

	// In case of early return, initialize. Note that energy *does not* accumulate
	energy = 0.0;
	deriv = ZERO_DERIV2D;

	if ( occ_rsd.is_virtual( occ_atm ) ) return;
	if ( acc_rsd.is_virtual( acc_atm ) ) return;

	// JUST FOR TESTING! DO NOT CHECK IN!
	//	if ( acc_rsd.atom_name( acc_atm ) !=  " O3'" ) return;
	//	if ( occ_rsd.atom_name( occ_atm ) !=  " C5'" ) return;

	//Need to know about backbone/backbone H-bonds for proteins
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( HBOND_SET )));
	TenANeighborGraph const & tenA_neighbor_graph
		( pose.energies().tenA_neighbor_graph() );

	assert( atom_is_acceptor( acc_rsd, acc_atm ) );

	Size const base_atm ( acc_rsd.atom_base( acc_atm ) );

	Vector const & acc_atm_xyz( acc_rsd.atom( acc_atm ).xyz() );

	Vector base_atm_xyz = get_acceptor_base_atm_xyz( acc_rsd, acc_atm );

	bool const acc_atm_is_protein_backbone
		( acc_rsd.is_protein() && acc_rsd.atom_is_backbone( acc_atm ) );

	// Virtual atom (e.g., Andrew Leaver-Fay's NV in proline) don't count.

	if ( occ_rsd.is_protein() && occ_rsd.is_virtual(occ_atm) ) return;
	if ( acc_rsd.is_protein() && acc_rsd.is_virtual(acc_atm) ) return;

	// the base atom isn't allowed to occlude solvent
	// Note: In current implementation, intraresidue pairs aren't checked...
	if ( ( acc_rsd.seqpos() == occ_rsd.seqpos() ) && ( occ_atm == base_atm ) ) return;

	// if a backbone acceptor participates in a backbone-backbone Hbond,
	// nothing is allowed to occlude solvent except a backbone donor
	bool const occ_atm_is_protein_backbone_donor
		( occ_rsd.is_protein() && occ_rsd.atom_is_backbone( occ_atm ) && atom_is_donor(occ_rsd, occ_atm ) );
	if ( acc_atm_is_protein_backbone && hbond_set.acc_bbg_in_bb_bb_hbond( acc_rsd.seqpos() ) &&
			 !occ_atm_is_protein_backbone_donor ) return;

	// an atom directly bound to the acceptor isn't allowed to occlude solvent
	// Note: In current implementation, intraresidue pairs aren't checked...
	if ( ( acc_rsd.seqpos() == occ_rsd.seqpos() ) && acc_rsd.path_distance( acc_atm, occ_atm ) < 2 ) return;

	// if the distance is > 5.2 A, from the acceptor, it doesn't occlude solvent
	Vector const & occ_atm_xyz( occ_rsd.atom( occ_atm ).xyz() );
	Real const acc_dis2 = ( occ_atm_xyz - acc_atm_xyz ).length_squared();
	if ( acc_dis2 > dist_cut2_ ) return;

	// if distance to base atom is greater than distance to donor, it doesn't occlude solvent
	Real const base_dis2 = ( occ_atm_xyz - base_atm_xyz ).length_squared();
	if ( acc_dis2 > base_dis2 ) return;

	// For a backbone NH occluding a backbone CO, use the backbone-backbone (linear) geometry
	// to compute solvation penalty
	bool const potential_backbone_backbone_hbond =
		( acc_atm_is_protein_backbone && occ_atm_is_protein_backbone_donor );
	// Following distinguished between helix/turn/etc. based on seq. separation. Is that what we want? Not in original jk code, so not for now.
	//			HBEvalType hbe = potential_backbone_backbone_hbond ? ( hbond_evaluation_type( occ_atm, occ_rsd, acc_atm, acc_rsd) ) : hbe_SP3SC;
	//			HBEvalType hbe = potential_backbone_backbone_hbond ? hbe_BSC : hbe_SP3SC;
 	HBEvalTuple hbe = potential_backbone_backbone_hbond ? ( HBEvalTuple( occ_atm, occ_rsd, acc_atm, acc_rsd ) ) :
 		HBEvalTuple( hbdon_H2O, get_hb_acc_chem_type( acc_atm, acc_rsd), seq_sep_other );

	Size acc_nbrs( 0 ), occ_nbrs( 0 );
	if (options_->hbond_options().use_hb_env_dep() ) {
		acc_nbrs = tenA_neighbor_graph.get_node( acc_rsd.seqpos() )->num_neighbors_counting_self();
		occ_nbrs = tenA_neighbor_graph.get_node( occ_rsd.seqpos() )->num_neighbors_counting_self();
	}

	// jk Compute the Hbond energy as if this was a water
	// jk Add the water Hbond energy to a running total sum, as well as to the residue sum
	energy = occluded_water_hbond_penalty(
		false /*is_donor*/, hbe,
		acc_atm_xyz, base_atm_xyz, occ_atm_xyz,
		acc_nbrs, occ_nbrs,
		update_deriv, deriv);

	if ( verbose_ && ( energy > 0.0 ) ) {
		tr<<"jk ACC res "<< acc_rsd.name1() << I(3, acc_rsd.seqpos())<<
			" atom "<< acc_rsd.atom_name( acc_atm )<<" is occluded by occ_res "<<
			occ_rsd.name1()<< I(3, occ_rsd.seqpos()) <<
			" atom "<< occ_rsd.atom_name( occ_atm ) <<
			"  (HBEvalType " <<  I(2,hbe.eval_type()) << ") " <<
			" with energy "<< F(8,3,energy)<<std::endl;

	}

}


///////////////////////////////////////////////////////////////////////////////
///    Compute the cosine required for calling water Hbond energies
inline
Real
GeometricSolEnergyEvaluator::get_water_cos( Vector const & base_atm_xyz,
																	 Vector const & polar_atm_xyz,
																	 Vector const & occluding_atm_xyz ) const
{
	return dot( (polar_atm_xyz - base_atm_xyz).normalize(),  (occluding_atm_xyz - polar_atm_xyz).normalize() );
}


//////////////////////////////////////////////////////////////////////////////////////
// Silly helper function
// These should probably live inside conformation::Residue.
//
bool
GeometricSolEnergyEvaluator::atom_is_donor( conformation::Residue const & rsd, Size const atm ) const
{
	for ( chemical::AtomIndices::const_iterator
			hnum  = rsd.Hpos_polar().begin(),
			hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		Size const don_base_atm( rsd.atom_base( don_h_atm ) );
		if ( don_base_atm == atm ) return true;
	}
	return false;
}
//////////////////////////////////////////////////////////////////////////////////////
// Silly helper function
// These should probably live inside conformation::Residue.
//
bool
GeometricSolEnergyEvaluator::atom_is_donor_h( conformation::Residue const & rsd, Size const atm ) const
{
	for ( chemical::AtomIndices::const_iterator
			hnum  = rsd.Hpos_polar().begin(),
			hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		if ( don_h_atm == atm ) return true;
	}
	return false;
}
//////////////////////////////////////////////////////////////////////////////
// Stupid helper function
// These should probably live inside conformation::Residue.
bool
GeometricSolEnergyEvaluator::atom_is_acceptor( conformation::Residue const & rsd, Size const atm ) const
{
	for ( chemical::AtomIndices::const_iterator
			anum  = rsd.accpt_pos().begin(),
			anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atm( *anum );
		if ( acc_atm == atm ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////
// Stupid helper function
// These should probably live inside conformation::Residue.
bool
GeometricSolEnergyEvaluator::atom_is_heavy( conformation::Residue const & rsd, Size const atm ) const
{
	//Could check if its hydrogen, but this is the same delineation used in the
	// residue-residue pair energy loop.
	return (atm <= rsd.nheavyatoms() && !rsd.is_virtual( atm ));
}

//////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergyEvaluator::eval_atom_derivative_intra_RNA(
		 id::AtomID const & atom_id,
		 pose::Pose const & pose,
		 EnergyMap const & weights,
		 Vector & F1,
		 Vector & F2
) const
{

	Size const i( atom_id.rsd() );

	conformation::Residue const & current_rsd( pose.residue( i ) );
	conformation::Residue const &   other_rsd( pose.residue( i ) );

	//Ok right now intrares energy is define only for the RNA case. Parin Sripakdeevong, June 26, 2011.
	if (current_rsd.is_RNA()==false) return;
	if (other_rsd.is_RNA()==false) return; //no effect!

	static bool const update_deriv( true );

	Real energy( 0.0 );
	hbonds::HBondDerivs deriv;

	Size const current_atm( atom_id.atomno() );

	if (verbose_) std::cout << "Start eval_atom_derivative, intra_res case, res= " << i << " atomno= " << current_atm << "[" << current_rsd.atom_name(current_atm) <<  "]" << std::endl;

	// If this atom is a donor, go over heavy atoms in other residue.
	if ( atom_is_donor_h( current_rsd, current_atm ) ) {
		for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){

			if (core::chemical::rna::Is_base_phosphate_atom_pair(current_rsd, other_rsd, current_atm, m)==false) continue;

			get_atom_atom_geometric_solvation_for_donor( current_atm, current_rsd, m, other_rsd, pose, energy, update_deriv, deriv );

			F1 += weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
			F2 += weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
		}
	}

	// If this atom is an acceptor, go over heavy atoms in other residue.
	if ( atom_is_acceptor( current_rsd, atom_id.atomno() ) ) {
		for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){

			if (core::chemical::rna::Is_base_phosphate_atom_pair(current_rsd, other_rsd, current_atm, m)==false) continue;

			get_atom_atom_geometric_solvation_for_acceptor( current_atm, current_rsd, m, other_rsd, pose, energy, update_deriv, deriv );

			F1 += weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
			F2 += weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
		}
	}

	//Treat atom as occluder if its heavy.
	if ( atom_is_heavy( current_rsd, atom_id.atomno() ) ) {
		//			Go over donors in other atom.
		for ( chemical::AtomIndices::const_iterator hnum  = other_rsd.Hpos_polar().begin(), hnume = other_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const don_h_atm( *hnum );

			if (core::chemical::rna::Is_base_phosphate_atom_pair(current_rsd, other_rsd, current_atm, don_h_atm)==false) continue;

			get_atom_atom_geometric_solvation_for_donor( don_h_atm, other_rsd, atom_id.atomno(), current_rsd, pose, energy, update_deriv, deriv );

			F1 -= weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
			F2 -= weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
		}

		// Go over acceptors in other atom.
		for ( chemical::AtomIndices::const_iterator anum  = other_rsd.accpt_pos().begin(), anume = other_rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const acc_atm ( *anum );

			if (core::chemical::rna::Is_base_phosphate_atom_pair(current_rsd, other_rsd, current_atm, acc_atm)==false) continue;

			get_atom_atom_geometric_solvation_for_acceptor( acc_atm, other_rsd, atom_id.atomno(), current_rsd, pose, energy, update_deriv, deriv );

			F1 -= weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
			F2 -= weights[ geom_sol_intra_RNA ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
		}
	}

	if (verbose_){
		std::cout << "eval_atom_derivative, intra_res :";
	 	std::cout << " F1= " << F1[0] << " " << F1[1] << " " << F1[2];
	 	std::cout << " F2= " << F2[0] << " " << F2[1] << " " << F2[2] << std::endl;
		std::cout << "Finish eval_atom_derivative, intra_res case, res= " << i << " atomno= " << current_atm << "[" << current_rsd.atom_name(current_atm) <<  "]" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Note that this computes every interaction *twice* -- three times if you
//  note that the score calculation above does most of the computation already.
// Oh well -- we currently assume derivative calculation doesn't happen too often!
//
// Also -- for context-independent calculations, should be able to use domain_map
//  to be smarter about what to calculate...
//
void
GeometricSolEnergyEvaluator::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	Real energy( 0.0 );
	hbonds::HBondDerivs deriv;

	eval_atom_derivative_intra_RNA(atom_id, pose, weights, F1, F2);

	conformation::Residue const & current_rsd( pose.residue( atom_id.rsd() ) );

	Size const i( atom_id.rsd() );
	Size const current_atm( atom_id.atomno() );

	//	Size const nres = pose.total_residue();
	static bool const update_deriv( true );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for( graph::Graph::EdgeListConstIter
			iter = energy_graph.get_node( i )->const_edge_list_begin();
			iter != energy_graph.get_node( i )->const_edge_list_end();
			++iter ){

		Size j( (*iter)->get_other_ind( i ) );

		if ( i == j ) continue; //Parin S. Already dealt with above for the intra_RNA case! Plus I think a edge doesn't point to itself!

		conformation::Residue const & other_rsd( pose.residue( j ) );

 		// If this atom is a donor, go over heavy atoms in other residue.
 		if ( atom_is_donor_h( current_rsd, current_atm ) ) {
 			for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){
				get_atom_atom_geometric_solvation_for_donor( current_atm, current_rsd,
					m, other_rsd,
					pose, energy, update_deriv, deriv );
				F1 += weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 += weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );

 			}
 		}

 		// If this atom is an acceptor, go over heavy atoms in other residue.
 		if ( atom_is_acceptor( current_rsd, atom_id.atomno() ) ) {
 			for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){
				get_atom_atom_geometric_solvation_for_acceptor( current_atm, current_rsd,
					m, other_rsd,
					pose, energy, update_deriv, deriv );
				F1 += weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 += weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );

 			}
 		}

		//Treat atom as occluder if its heavy.
 		if ( atom_is_heavy( current_rsd, atom_id.atomno() ) ) {
			//			Go over donors in other atom.
			for ( chemical::AtomIndices::const_iterator
							hnum  = other_rsd.Hpos_polar().begin(),
							hnume = other_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
				Size const don_h_atm( *hnum );
				get_atom_atom_geometric_solvation_for_donor( don_h_atm, other_rsd,
																										 atom_id.atomno(), current_rsd,
																										 pose, energy, update_deriv, deriv );
				F1 -= weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 -= weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );

			}

			// Go over acceptors in other atom.
			for ( chemical::AtomIndices::const_iterator
					anum  = other_rsd.accpt_pos().begin(),
					anume = other_rsd.accpt_pos().end(); anum != anume; ++anum ) {
				Size const acc_atm ( *anum );
				get_atom_atom_geometric_solvation_for_acceptor( acc_atm, other_rsd,
					atom_id.atomno(), current_rsd,
					pose, energy, update_deriv, deriv );
				F1 -= weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 -= weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );

			}
 		}

	}

}

////////////////////////////////////////////////////////////////////
// Only return energy for occluded polar atoms.
Real
GeometricSolEnergyEvaluator::eval_atom_energy(
	id::AtomID const & atom_id,
	pose::Pose const & pose
) const
{

	Real total_energy( 0.0 );

	conformation::Residue const & current_rsd( pose.residue( atom_id.rsd() ) );

	Size const i( atom_id.rsd() );
	Size const current_atm( atom_id.atomno() );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for( graph::Graph::EdgeListConstIter
				 iter = energy_graph.get_node( i )->const_edge_list_begin();
			 iter != energy_graph.get_node( i )->const_edge_list_end();
			 ++iter ){

		Size j( (*iter)->get_other_ind( i ) );

		//As above disallow contribution within a residue. Is this OK? Sure there shouldn't be an "intra" term?
		// anyway list of neighbors does not include i==j.
		if ( i == j ) continue;

		conformation::Residue const & other_rsd( pose.residue( j ) );

 		// If this atom is a donor, go over heavy atoms in other residue.
 		if ( atom_is_donor_h( current_rsd, current_atm ) ) {
 			for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){
				Real energy( 0.0 );
				get_atom_atom_geometric_solvation_for_donor( current_atm, current_rsd,
					m, other_rsd, pose, energy );
				total_energy += energy;
 			}
 		}

 		// If this atom is an acceptor, go over heavy atoms in other residue.
 		if ( atom_is_acceptor( current_rsd, atom_id.atomno() ) ) {
 			for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){
				Real energy( 0.0 );
				get_atom_atom_geometric_solvation_for_acceptor( current_atm, current_rsd,
					m, other_rsd, pose, energy );
				total_energy += energy;
 			}
 		}

	}

	return total_energy;
}


///////////////////////////////////////////////////////////
/// Note -- this is used in HBonds to
/// add up backbone/backbone energies,
/// while the residue_pair_energy does the
/// sidechain stuff. For now make this an
/// empty function, although it will come back into
/// the game if we decide to speed up the packer.
///////////////////////////////////////////////////////////
//void
//GeometricSolEnergyEvaluator::finalize_total_energy(
//	pose::Pose &,
//	ScoreFunction const &,
//	EnergyMap & totals
//) const
//{
//  using core::scoring::EnergiesCacheableDataType::HBOND_SET;
//
// 	hbonds::GeometricSolSet const & hbond_set
// 		( static_cast< hbonds::GeometricSolSet const & >
// 			( pose.energies().data().get( HBOND_SET )));

// 	Real lr_bbE( 0.0 ), sr_bbE( 0.0 ), bb_scE( 0.0 ), scE( 0.0 );

// 	get_hbond_energies( hbond_set, sr_bbE, lr_bbE, bb_scE, scE );

// 	// the current logic is that we fill the hbond set with backbone
// 	// hbonds only at the beginning of scoring. this is done to setup
// 	// the bb-bb hbond exclusion logic. so the hbondset should only
// 	// include bb-bb hbonds.
// 	// but see get_hb_don_chem_type in hbonds_geom.cc -- that only
// 	// classifies protein backbone donors as backbone, and the energy
// 	// accumulation by type is influenced by that via HBeval_lookup
// 	//
// 	// the important thing is that there's no double counting, which
// 	// is I think true since both fill_hbond_set and get_rsd-rsd-energy
// 	// use atom_is_backbone to check...
// 	//assert( std::abs( bb_scE ) < 1e-3 && std::abs( scE ) < 1e-3 );
// 	totals[ hbond_sr_bb ] += sr_bbE;
// 	totals[ hbond_lr_bb ] += lr_bbE;
//}


// COPIED OVER FROM HBondEnergy.cc ==> comment is not rhiju's!
///@brief HACK!  MAX_R defines the maximum donorH to acceptor distance.
// The atomic_interaction_cutoff method is meant to return the maximum distance
// between two *heavy atoms* for them to have a zero interaction energy.
// I am currently assuming a 1.35 A maximum distance between a hydrogen and the
// heavy atom it is bound to, stealing this number from the CYS.params file since
// the HG in CYS is much further from it's SG than aliphatic hydrogens are from their carbons.
// This is a bad idea.  Someone come up with a way to fix this!
//
// At 4.35 A interaction cutoff, the hbond energy function is incredibly short ranged!
Distance
GeometricSolEnergyEvaluator::atomic_interaction_cutoff() const
{
	return MAX_R + 1.35; // MAGIC NUMBER
}

///////////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergyEvaluator::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap & emap
) const{

	if (rsd.is_RNA()==false) return;

		Real geo_solE_intra_RNA =
			   donorRes_occludingRes_geometric_sol_RNA_intra( rsd, pose ) +
			acceptorRes_occludingRes_geometric_sol_RNA_intra( rsd, pose );

	// store the energies
	emap[ geom_sol_intra_RNA ] += geo_solE_intra_RNA;

}

///////////////////////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::donorRes_occludingRes_geometric_sol_RNA_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose ) const
{


	Real res_solE( 0.0 ), energy( 0.0 );

	conformation::Residue const & don_rsd=rsd;
	conformation::Residue const & occ_rsd=rsd;

	// Here we go -- cycle through polar hydrogens in don_aa, everything heavy in occluding atom.
	for ( chemical::AtomIndices::const_iterator hnum  = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			if (core::chemical::rna::Is_base_phosphate_atom_pair(rsd, rsd, occ_atm, don_h_atm)==false) continue;

			get_atom_atom_geometric_solvation_for_donor( don_h_atm, don_rsd, occ_atm, occ_rsd, pose, energy );
			res_solE += energy;
		}
	}

	return res_solE;
}

///////////////////////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::acceptorRes_occludingRes_geometric_sol_RNA_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose ) const
{

	conformation::Residue const & acc_rsd=rsd;
	conformation::Residue const & occ_rsd=rsd;


	Real res_solE( 0.0 ), energy( 0.0 );

	for ( chemical::AtomIndices::const_iterator anum  = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atm( *anum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			if (core::chemical::rna::Is_base_phosphate_atom_pair(rsd, rsd, occ_atm, acc_atm)==false) continue;

			get_atom_atom_geometric_solvation_for_acceptor( acc_atm, acc_rsd, occ_atm, occ_rsd, pose, energy);
			res_solE += energy;

		}
	}

	return res_solE;
}


} //geometric_solvation
} //scoring
} //core
