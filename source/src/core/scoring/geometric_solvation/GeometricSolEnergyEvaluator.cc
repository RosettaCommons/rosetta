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
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/chemical/rna/util.hh>

// Project headers
#include <core/chemical/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/types.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <numeric/xyz.io.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "core.scoring.geometric_solvation.GeometricSolEnergyEvaluator" );

using namespace core::scoring::hbonds;
using namespace ObjexxFCL::format;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
// When a heavy-atom occludes an acceptor or donor, estimate
//  the cost of displacing a water that was at that position.
//
// Actually involves placement of a pseudo-water atom with
//  perfect orientation at heavy-atom position, and a
//  snazzy reuse of the hbond derivative framework.
//
// Overhauled by the Das lab in 2013-2014, especially
// to fix derivatives to match new hbond system.
//
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

/// @brief copy c-tor
GeometricSolEnergyEvaluator::GeometricSolEnergyEvaluator( methods::EnergyMethodOptions const & opts ) :
	options_( opts ),
	hb_database_( HBondDatabase::get_database( opts.hbond_options().params_database_tag() )),
	dist_cut2_( 27.0 ),   // 5.2*5.2
	geometric_sol_scale_( 0.4 * 1.17 / 0.65 ),
	interres_path_distance_cutoff_( opts.geom_sol_interres_path_distance_cutoff() ), // 0, currently no cutoff here -- so this is counting way too many pairs, even for bonded atoms!
	intrares_path_distance_cutoff_( opts.geom_sol_intrares_path_distance_cutoff() ), // 6, note that this counts very few pairs
	verbose_( false )
{
	//runtime_assert( !options_.exclude_DNA_DNA() );	//GEOMETRIC SOLVATION NOT COMPATIBLE WITH EXCLUDE_DNA_DNA FLAG YET
}

/// copy ctor
GeometricSolEnergyEvaluator::GeometricSolEnergyEvaluator( GeometricSolEnergyEvaluator const & src ) :
		utility::pointer::ReferenceCount(src),
		options_( src.options_ ),
		hb_database_( src.hb_database_ ),
		dist_cut2_( src.dist_cut2_ ),   // 5.2*5.2
		geometric_sol_scale_( src.geometric_sol_scale_ ),
		interres_path_distance_cutoff_( src.interres_path_distance_cutoff_ ),
		intrares_path_distance_cutoff_( src.intrares_path_distance_cutoff_ ),
		verbose_( src.verbose_ )
{}

//Destructor
GeometricSolEnergyEvaluator::~GeometricSolEnergyEvaluator()
{}


/////////////////////////////////////////////////////////////////////////////
// This is meant to be a reasonable clone of John Karanicolas' original (pre-2008!)
//  Rosetta++ code. Over the years, the code has been optimized for packing/designing/minimizing
//  by members of the Das lab -- rhiju
//
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


//////////////////////////////////////////////////////////////////////////////////////
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
//optimization functions for packing/design.
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

//////////////////////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////////
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

    if(!don_rsd.atom_is_backbone(don_h_atm)) continue;

    for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

      if(!occ_rsd.atom_is_backbone(occ_atm)) continue;

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
	y = water_v - atom_v;
	y.normalize();

	z = cross( y,  (atom_v - base_v) );
	if ( z.length() == 0.0 )	z = cross( water_v - atom_v, water_v + atom_v ); // arbitrary
	runtime_assert( z.length() > 0.0 );
	z.normalize();

	x = cross( y, z );

	// Plop the atom down
	direction = xH * y  + Real( std::sqrt( 1 - (xH * xH) ) ) * x;
	water_base_v = water_v + bond_length * direction;

}


/////////////////////////////////////////////////////////////////////////////////////////
// jumpout criteria copied from hb_energy_deriv in hbonds.cc
inline
Real
GeometricSolEnergyEvaluator::occluded_water_hbond_penalty(
  bool const & is_donor,
	hbonds::HBEvalTuple const & hbond_eval_type,
	Vector const & polar_atm_xyz,
	Vector const & base_atm_xyz,
	Vector const & base2_atm_xyz,
	Vector const & occluding_atm_xyz,
	Size const & polar_nb,
	Size const & occ_nb,
	bool const update_deriv /* = false*/,
	HBondDerivs & deriv /* = DUMMY_DERIV2D*/ ) const
{

	// toggling this on turns on derivative calculations only when necessary
	// this actually does not help speed detectably. however, checking that
	// we get the same answer with true/false helps catch inconsistencies between
	// how hbonds are dealt with here vs. inside the hbonds_geom.cc machinery.
	static bool const always_do_full_calculation( false );

	// Compute geometry
	Real AHdis( 0.0 ), xD( 0.0 ), xH( 0.0 ), energy( 0.0 );

	// Create an artifical atom to make a complete pseudo "water molecule" at the occluding heavy atom.
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
		if ( update_deriv || always_do_full_calculation ) {
			Vector occluding_base_atm_xyz( 0.0 );
			set_water_base_atm( base_atm_xyz, polar_atm_xyz, occluding_atm_xyz, occluding_base_atm_xyz,
													xH, water_O_H_distance );
			hb_energy_deriv( *hb_database_, options_.hbond_options(),
				hbond_eval_type, base_atm_xyz, polar_atm_xyz,
				occluding_atm_xyz,
				occluding_base_atm_xyz,
				occluding_base_atm_xyz,
				energy, hbderiv_ABE_GO_GEOMSOL_OCC_DON, deriv);
		} else {
			Real chi( 0.0 );
			hbond_compute_energy( *hb_database_, options_.hbond_options(), hbond_eval_type, AHdis, xD, xH, chi, energy );
		}
	} else {

		// water is the donor, give it perfect geometry. This used to be 0.9999, but that caused imperfect derivatives.
		xD = 1.0;

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
		Vector pseudo_base_atm_xyz, dummy; // stolen from hbonds_geom.cc
		make_hbBasetoAcc_unitvector( options_.hbond_options(),
																 get_hbe_acc_hybrid( hbond_eval_type.eval_type() ),
																 polar_atm_xyz,
																 base_atm_xyz,
																 base2_atm_xyz,
																 pseudo_base_atm_xyz, dummy );

		xH = get_water_cos( pseudo_base_atm_xyz, polar_atm_xyz, occluding_atm_xyz );
		if ( xH < MIN_xH ) return 0.;
		if ( xH > MAX_xH ) return 0.;

		//rhiju, testing alternative calculation that will give derivative.
		Vector occluding_mock_hydrogen_atm_xyz( 0.0 );
		set_water_base_atm( pseudo_base_atm_xyz, polar_atm_xyz, occluding_atm_xyz, occluding_mock_hydrogen_atm_xyz,
												-xD, water_O_H_distance );
		if ( update_deriv || always_do_full_calculation || true ) {
 			hb_energy_deriv( *hb_database_, options_.hbond_options(), hbond_eval_type,
											 occluding_atm_xyz, occluding_mock_hydrogen_atm_xyz,
											 polar_atm_xyz, base_atm_xyz /*pseudo_base determined inside!*/, base2_atm_xyz,
											 energy, hbderiv_ABE_GO_GEOMSOL_OCC_ACC, deriv);
		} else {
			Real chi( 0.0 );
			// stolen from hbonds_geom.cc -- probably should make a shared function.
			if ( options_.hbond_options().use_sp2_chi_penalty() &&
					 get_hbe_acc_hybrid( hbond_eval_type.eval_type() ) == chemical::SP2_HYBRID &&
					 base2_atm_xyz != Vector(-1.0, -1.0, -1.0) ) {
				chi = numeric::dihedral_radians( occluding_mock_hydrogen_atm_xyz, polar_atm_xyz, pseudo_base_atm_xyz, base2_atm_xyz );
			} else if ( options_.hbond_options().measure_sp3acc_BAH_from_hvy() &&
									( hbond_eval_type.acc_type() == hbacc_AHX || hbond_eval_type.acc_type() == hbacc_HXL ) ) {
				/// Bxyz really is the heavy atom base and B2xyz really is the hydroxyl hydrogen
				/// this is guaranteed by the hbond_measure_sp3acc_BAH_from_hvy flag.
				chi = numeric::dihedral_radians( occluding_mock_hydrogen_atm_xyz, polar_atm_xyz, pseudo_base_atm_xyz, base2_atm_xyz );
			}
			hbond_compute_energy( *hb_database_, options_.hbond_options(), hbond_eval_type, AHdis, xD, xH, chi, energy );
		}
	}

	if (options_.hbond_options().use_hb_env_dep() ) {
		environment_weight = hbonds::get_environment_dependent_weight( hbond_eval_type, polar_nb, occ_nb, options_.hbond_options() );
	}

	if (verbose_ ) TR << "  jk ENERGY: " << energy << std::endl;
	core::Real sol_penalty = -1.0 * energy;

	// jk THIS NEEDS TO BE FIT MORE RIGOROUSLY LATER...
	// Apply a scaling factor (effectively a weight), tying the weight of the
	// solvation term to the Hbond term (rather than to the LK weight)
	Real reweight = geometric_sol_scale_ * environment_weight;
	sol_penalty *= reweight;
	deriv.h_deriv.f1()      *= reweight;
	deriv.h_deriv.f2()      *= reweight;
	deriv.acc_deriv.f1()    *= reweight;
	deriv.acc_deriv.f2()    *= reweight;
	deriv.abase_deriv.f1()  *= reweight;
	deriv.abase_deriv.f2()  *= reweight;
	deriv.abase2_deriv.f1() *= reweight;
	deriv.abase2_deriv.f2() *= reweight;
	deriv.don_deriv.f1()    *= reweight;
	deriv.don_deriv.f2()    *= reweight;

	// this is a penalty, don't return a negative number
	if (sol_penalty < 0.) {
		if ( update_deriv ) deriv = ZERO_DERIV2D;
		return 0.;
	}

	return sol_penalty; // return a positive number (this is a penalty)

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
GeometricSolEnergyEvaluator::check_path_distance(
	 conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 Size const & atm1,
	 Size const & atm2 ) const {

	path_distance_ = 0;

	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		return ( rsd1.path_distance( atm1, atm2 ) >= intrares_path_distance_cutoff_ );
	}

	if ( ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) && ( interres_path_distance_cutoff_ > 0 ) ) {
		etable::count_pair::CountPairGeneric count_pair( rsd1, rsd2 ); // this is inefficient... happens with each atom pair in rsd1, rsd2.
		path_distance_ = count_pair.path_distance( atm1, atm2 );
		return ( path_distance_ >= interres_path_distance_cutoff_ );
	}

	return true;
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
	HBondDerivs & deriv /* = DUMMY_DERIVS */,
	HBEvalTuple & hbe /* = HBEvalTuple() */
) const
{
	//Why do we need to send in the pose and the residue stuff?
	// Well, the pose has info on backbone H-bonds.
	// and, during design, the residue type doesn't actually
	// have to match what's in the pose! TRicky!

	// In case of early return, initialize. Note that energy *does not* accumulate
	energy = 0.0;
	deriv = ZERO_DERIV2D;

	if ( occ_rsd.is_virtual( occ_atm ) ) return;
	if ( don_rsd.is_virtual( don_h_atm ) ) return;

	debug_assert( don_rsd.atom_is_polar_hydrogen( don_h_atm ) );

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
		( occ_rsd.is_protein() && occ_rsd.atom_is_backbone( occ_atm ) && occ_rsd.heavyatom_is_an_acceptor( occ_atm ) );
	if ( don_h_atm_is_protein_backbone && !occ_atm_is_protein_backbone_acceptor ) {
		using core::scoring::EnergiesCacheableDataType::HBOND_SET;
		if ( pose.energies().data().has( HBOND_SET ) ){
			hbonds::HBondSet const & hbond_set
				( static_cast< hbonds::HBondSet const & >
					( pose.energies().data().get( HBOND_SET )));
			if ( hbond_set.don_bbg_in_bb_bb_hbond( don_rsd.seqpos() ) )	return;
		}
	}

	// an atom directly bound to the donor isn't allowed to occlude solvent
	if ( !check_path_distance( don_rsd, occ_rsd, don_h_atm, occ_atm ) ) return;

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
	hbe = potential_backbone_backbone_hbond ? ( HBEvalTuple( don_base_atm, don_rsd, occ_atm, occ_rsd) ) : HBEvalTuple( hbdon_HXL, hbacc_HXL, seq_sep_other ); // apl note: the false condition maps to hbe_dHXLaHXL

	Size don_nbrs( 0 ), occ_nbrs( 0 );
	if (options_.hbond_options().use_hb_env_dep() ) {
		//Need to know about backbone/backbone H-bonds for proteins
		TenANeighborGraph const & tenA_neighbor_graph
			( pose.energies().tenA_neighbor_graph() );
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
		don_h_atm_xyz, don_base_atm_xyz, don_base_atm_xyz /*base2, not in use*/,
		occ_atm_xyz,
		don_nbrs, occ_nbrs,
		update_deriv, deriv);

	if ( verbose_ && ( energy > 0.0 ) ) {
		TR << " DON res "<< don_rsd.name1() << I(3,don_rsd.seqpos())<<
			" atom "<< don_rsd.atom_name( don_h_atm )<<" is occluded by occ_res " <<
			occ_rsd.name1()<< I(3, occ_rsd.seqpos()) <<
			" atom "<< occ_rsd.atom_name( occ_atm ) <<
			" path-distance " << path_distance_ <<
			"  (HBEvalType " <<  I(2,hbe.eval_type()) << ") " <<
			" with energy "<< F(8,3,energy) << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// following handles the special case in which acceptor has two base residues -- occurs for
//  N inside rings.
// This matches machinery in hbonds_geom.cc.  That doesn't mean that the base atom is set
//  totally correctly -- water (TIP3.params) and O4' in nucleic acids still have weird base atoms,
//  but at least the hbonds and geom_sol match up.
inline
Vector
GeometricSolEnergyEvaluator::get_acceptor_base_atm_xyz( conformation::Residue const & acc_rsd, Size const & acc_atm,
																												hbonds::HBEvalTuple const & hbt ) const{

	Vector base_atm_xyz, dummy;
	chemical::Hybridization acc_hybrid( get_hbe_acc_hybrid( hbt.eval_type() ) ); //acc_rsd.atom_type( acc_atm ).hybridization());
	make_hbBasetoAcc_unitvector( options_.hbond_options(),
															 acc_hybrid,
															 acc_rsd.atom( acc_atm ).xyz(),
															 acc_rsd.xyz( acc_rsd.atom_base( acc_atm ) ),
															 acc_rsd.xyz( acc_rsd.abase2( acc_atm ) ),
															 base_atm_xyz, dummy );
	return base_atm_xyz;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergyEvaluator::get_atom_atom_geometric_solvation_for_acceptor(
	Size const & acc_atm,
	conformation::Residue const & acc_rsd,
	Size const & occ_atm,
	conformation::Residue const & occ_rsd,
	pose::Pose const & pose,
	Real & energy,
	bool const update_deriv /*= false*/,
	HBondDerivs & deriv /* = DUMMY_DERIV2D */,
	HBEvalTuple & hbe /* = HBEvalTuple() */
) const
{

	//Why do we need to send in the pose and the residue stuff?
	// Well, the pose has info on backbone H-bonds.
	// and, during design, the residue type doesn't actually
	// have to match what's in the pose! TRicky!

	// In case of early return, initialize. Note that energy *does not* accumulate
	energy = 0.0;
	deriv = ZERO_DERIV2D;

	if ( occ_rsd.is_virtual( occ_atm ) ) return;
	if ( acc_rsd.is_virtual( acc_atm ) ) return;

	debug_assert( acc_rsd.heavyatom_is_an_acceptor( acc_atm ) );

	bool const acc_atm_is_protein_backbone
		( acc_rsd.is_protein() && acc_rsd.atom_is_backbone( acc_atm ) );

	// if a backbone acceptor participates in a backbone-backbone Hbond,
	// nothing is allowed to occlude solvent except a backbone donor
	bool const occ_atm_is_protein_backbone_donor
		( occ_rsd.is_protein() && occ_rsd.atom_is_backbone( occ_atm ) &&
			occ_rsd.heavyatom_has_polar_hydrogens( occ_atm ) );
	if ( acc_atm_is_protein_backbone && !occ_atm_is_protein_backbone_donor) {
		//Need to know about backbone/backbone H-bonds for proteins
		using core::scoring::EnergiesCacheableDataType::HBOND_SET;
		if ( pose.energies().data().has( HBOND_SET ) ){
			hbonds::HBondSet const & hbond_set
				( static_cast< hbonds::HBondSet const & >
					( pose.energies().data().get( HBOND_SET )));
			if ( hbond_set.acc_bbg_in_bb_bb_hbond( acc_rsd.seqpos() ) ) return;
		}
	}

	// For a backbone NH occluding a backbone CO, use the backbone-backbone (linear) geometry
	// to compute solvation penalty
	bool const potential_backbone_backbone_hbond =
		( acc_atm_is_protein_backbone && occ_atm_is_protein_backbone_donor );
	// Following distinguished between helix/turn/etc. based on seq. separation. Is that what we want? Not in original jk code, so not for now.
	//			HBEvalType hbe = potential_backbone_backbone_hbond ? ( hbond_evaluation_type( occ_atm, occ_rsd, acc_atm, acc_rsd) ) : hbe_SP3SC;
	//			HBEvalType hbe = potential_backbone_backbone_hbond ? hbe_BSC : hbe_SP3SC;
	hbe = potential_backbone_backbone_hbond ? ( HBEvalTuple( occ_atm, occ_rsd, acc_atm, acc_rsd ) ) :
 		HBEvalTuple( hbdon_H2O, get_hb_acc_chem_type( acc_atm, acc_rsd), seq_sep_other );

	Size const base_atm ( acc_rsd.atom_base( acc_atm ) );

	Vector const & acc_atm_xyz( acc_rsd.atom( acc_atm ).xyz() );

	Vector base_atm_xyz  = acc_rsd.xyz( acc_rsd.atom_base(  acc_atm ) );//get_acceptor_base_atm_xyz( acc_rsd, acc_atm, hbe );
	Vector base2_atm_xyz = acc_rsd.xyz( acc_rsd.abase2( acc_atm ) );

	// Virtual atom (e.g., Andrew Leaver-Fay's NV in proline) don't count.
	if ( occ_rsd.is_protein() && occ_rsd.is_virtual(occ_atm) ) return;
	if ( acc_rsd.is_protein() && acc_rsd.is_virtual(acc_atm) ) return;

	// the base atom isn't allowed to occlude solvent
	// Note: In current implementation, intraresidue pairs aren't checked...
	if ( ( acc_rsd.seqpos() == occ_rsd.seqpos() ) && ( occ_atm == base_atm ) ) return;

	// an atom directly bound to the acceptor isn't allowed to occlude solvent
	if ( !check_path_distance( acc_rsd, occ_rsd, acc_atm, occ_atm ) ) return;

	// if the distance is > 5.2 A, from the acceptor, it doesn't occlude solvent
	Vector const & occ_atm_xyz( occ_rsd.atom( occ_atm ).xyz() );
	Real const acc_dis2 = ( occ_atm_xyz - acc_atm_xyz ).length_squared();
	if ( acc_dis2 > dist_cut2_ ) return;

	// if distance to base atom is greater than distance to donor, it doesn't occlude solvent
	Real const base_dis2 = ( occ_atm_xyz - base_atm_xyz ).length_squared();
	if ( acc_dis2 > base_dis2 ) return;

	Size acc_nbrs( 0 ), occ_nbrs( 0 );
	if (options_.hbond_options().use_hb_env_dep() ) {
		TenANeighborGraph const & tenA_neighbor_graph
			( pose.energies().tenA_neighbor_graph() );
		acc_nbrs = tenA_neighbor_graph.get_node( acc_rsd.seqpos() )->num_neighbors_counting_self();
		occ_nbrs = tenA_neighbor_graph.get_node( occ_rsd.seqpos() )->num_neighbors_counting_self();
	}

	// jk Compute the Hbond energy as if this was a water
	// jk Add the water Hbond energy to a running total sum, as well as to the residue sum
	energy = occluded_water_hbond_penalty(
		false /*is_donor*/, hbe,
		acc_atm_xyz, base_atm_xyz, base2_atm_xyz,
		occ_atm_xyz,
		acc_nbrs, occ_nbrs,
		update_deriv, deriv);

	if ( verbose_ && ( energy > 0.0 ) ) {
		TR << " ACC res "<< acc_rsd.name1() << I(3, acc_rsd.seqpos())<<
			" atom "<< acc_rsd.atom_name( acc_atm )<<" is occluded by occ_res "<<
			occ_rsd.name1()<< I(3, occ_rsd.seqpos()) <<
			" atom "<< occ_rsd.atom_name( occ_atm ) <<
			" path-distance " << path_distance_ <<
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

//////////////////////////////////////////////////////////////////////////////
// Silly helper function
// These should probably live inside conformation::Residue.
bool
GeometricSolEnergyEvaluator::atom_is_heavy( conformation::Residue const & rsd, Size const atm ) const
{
	//Could check if its hydrogen, but this is the same delineation used in the
	// residue-residue pair energy loop.
	return (atm <= rsd.nheavyatoms() && !rsd.is_virtual( atm ));
}

// COPIED OVER FROM HBondEnergy.cc ==> comment is not rhiju's!
/// @brief HACK!  MAX_R defines the maximum donorH to acceptor distance.
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
) const {

	Real geo_solE_intra_RNA =
		donorRes_occludingRes_geometric_sol_intra( rsd, pose, true /*just_RNA*/ ) +
		acceptorRes_occludingRes_geometric_sol_intra( rsd, pose, true /*just_RNA*/ );
	emap[ geom_sol_intra_RNA ] += geo_solE_intra_RNA;

	if ( options_.put_intra_into_total() ){
		Real geo_solE_intra =
			donorRes_occludingRes_geometric_sol_intra( rsd, pose, false /*just_RNA*/ ) +
			acceptorRes_occludingRes_geometric_sol_intra( rsd, pose, false /*just_RNA*/ );
		emap[ geom_sol ] += geo_solE_intra;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::donorRes_occludingRes_geometric_sol_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	bool const just_RNA /* legacy */ ) const
{
	Real res_solE( 0.0 ), energy( 0.0 );
	if ( !just_RNA && !calculate_intra_res_hbonds( rsd, options_.hbond_options() ) ) return res_solE;

	conformation::Residue const & don_rsd=rsd;
	conformation::Residue const & occ_rsd=rsd;

	// Here we go -- cycle through polar hydrogens in don_aa, everything heavy in occluding atom.
	for ( chemical::AtomIndices::const_iterator hnum  = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			if ( just_RNA && !core::chemical::rna::is_base_phosphate_atom_pair(rsd, rsd, occ_atm, don_h_atm) ) continue;

			get_atom_atom_geometric_solvation_for_donor( don_h_atm, don_rsd, occ_atm, occ_rsd, pose, energy );
			res_solE += energy;

		}
	}

	return res_solE;
}

///////////////////////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::acceptorRes_occludingRes_geometric_sol_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	bool const just_RNA /* legacy */ ) const
{

	Real res_solE( 0.0 ), energy( 0.0 );
	if ( !just_RNA && !calculate_intra_res_hbonds( rsd, options_.hbond_options() ) ) return res_solE;

	conformation::Residue const & acc_rsd=rsd;
	conformation::Residue const & occ_rsd=rsd;

	for ( chemical::AtomIndices::const_iterator anum  = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atm( *anum );
		for ( Size occ_atm = 1; occ_atm <= occ_rsd.nheavyatoms(); occ_atm++ ) {

			if ( just_RNA && !core::chemical::rna::is_base_phosphate_atom_pair(rsd, rsd, occ_atm, acc_atm) ) continue;

			get_atom_atom_geometric_solvation_for_acceptor( acc_atm, acc_rsd, occ_atm, occ_rsd, pose, energy);
			res_solE += energy;

		}
	}

	return res_solE;
}

////////////////////////////////////////////////////////////////////
// useful for coloring PDBs
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

		if ( i == j ) continue;

		conformation::Residue const & other_rsd( pose.residue( j ) );

 		// If this atom is a donor, go over heavy atoms in other residue.
 		if ( current_rsd.atom_is_polar_hydrogen( current_atm ) ) {
 			for (Size m = 1; m <= other_rsd.nheavyatoms(); m++ ){
				Real energy( 0.0 );
				get_atom_atom_geometric_solvation_for_donor( current_atm, current_rsd,
					m, other_rsd, pose, energy );
				total_energy += energy;
 			}
 		}

 		// If this atom is an acceptor, go over heavy atoms in other residue.
 		if (  current_rsd.heavyatom_is_an_acceptor( atom_id.atomno() ) ) {
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


//////////////////////////////////////////////////////////////////////////////////////////////
void
fill_atom_derivs_for_donor( hbonds::HBondDerivs const & deriv,
														core::conformation::Residue const & don_rsd,
														Size const don_hatm,
														Size const occatm,
														utility::vector1< DerivVectorPair > & don_atom_derivs,
														utility::vector1< DerivVectorPair > & occ_atom_derivs,
														hbonds::HBEvalTuple const & /*hbe_type*/,
														hbonds::HBondOptions const & /*hbond_options*/,
														Real const & weighted_energy )
{
	Size const donatm = don_rsd.atom_base( don_hatm );
	don_atom_derivs[ donatm ].f1() += weighted_energy * deriv.don_deriv.f1();
	don_atom_derivs[ donatm ].f2() += weighted_energy * deriv.don_deriv.f2();

	don_atom_derivs[ don_hatm ].f1() += weighted_energy * deriv.h_deriv.f1();
	don_atom_derivs[ don_hatm ].f2() += weighted_energy * deriv.h_deriv.f2();

	// the other side holds the 'acceptor' and its base atoms -- all get summed onto occluding atom.
	occ_atom_derivs[ occatm ].f1() += weighted_energy * deriv.acc_deriv.f1();
	occ_atom_derivs[ occatm ].f2() += weighted_energy * deriv.acc_deriv.f2();

	occ_atom_derivs[ occatm ].f1() += weighted_energy * deriv.abase_deriv.f1();
	occ_atom_derivs[ occatm ].f2() += weighted_energy * deriv.abase_deriv.f2();

	occ_atom_derivs[ occatm ].f1() += weighted_energy * deriv.abase2_deriv.f1();
	occ_atom_derivs[ occatm ].f2() += weighted_energy * deriv.abase2_deriv.f2();
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
fill_atom_derivs_for_acceptor( hbonds::HBondDerivs const & deriv,
															 core::conformation::Residue const & acc_rsd,
															 Size const aatm,
															 Size const occatm,
															 utility::vector1< DerivVectorPair > & acc_atom_derivs,
															 utility::vector1< DerivVectorPair > & occ_atom_derivs,
															 hbonds::HBEvalTuple const & hbe_type,
															 hbonds::HBondOptions const & hbond_options,
															 Real const & weighted_energy )
{

	acc_atom_derivs[ aatm ].f1() += weighted_energy * deriv.acc_deriv.f1();
	acc_atom_derivs[ aatm ].f2() += weighted_energy * deriv.acc_deriv.f2();

	// ring-acceptor derivative assignment logic is TRicky
	assign_abase_derivs( hbond_options, acc_rsd, aatm, hbe_type, deriv.abase_deriv, weighted_energy, acc_atom_derivs );

	int const base2( acc_rsd.abase2( aatm ) );
	acc_atom_derivs[ base2 ].f1() += weighted_energy * deriv.abase2_deriv.f1();
	acc_atom_derivs[ base2 ].f2() += weighted_energy * deriv.abase2_deriv.f2();

	// the other side holds both the 'donor' and the 'donor hydrogen'
	occ_atom_derivs[ occatm ].f1() += weighted_energy * deriv.don_deriv.f1();
	occ_atom_derivs[ occatm ].f2() += weighted_energy * deriv.don_deriv.f2();

	occ_atom_derivs[ occatm ].f1() += weighted_energy * deriv.h_deriv.f1();
	occ_atom_derivs[ occatm ].f2() += weighted_energy * deriv.h_deriv.f2();
}

//////////////////////////////////////////////////////////////////////////////
// copies some code from eval_residue_pair_derivatives, but that's necessary, I think
// if we want to retain Parin's base/phosphate check. However, I'd like to replace that
// with a CountPair call, eventually, in which case we should unify!
void
GeometricSolEnergyEvaluator::eval_intrares_derivatives(
		 conformation::Residue const & rsd,
		 pose::Pose const & pose,
		 Real const & geom_sol_intra_weight,
		 utility::vector1< DerivVectorPair > & atom_derivs,
		 bool const just_RNA
) const
{

	if ( just_RNA && !rsd.is_RNA() ) return;
	if ( !just_RNA && !calculate_intra_res_hbonds( rsd, options_.hbond_options() ) ) return;

	Real energy( 0.0 );
	hbonds::HBondDerivs deriv;
	hbonds::HBEvalTuple hbe_type;
	bool const update_deriv( true );
	Real const weighted_energy = -1.0 * geom_sol_intra_weight;

	for ( Size ii = 1; ii <= rsd.natoms(); ii++ ) {
		for ( Size jj = (ii+1); jj <= rsd.natoms(); jj++ ) {

			if ( just_RNA && !core::chemical::rna::is_base_phosphate_atom_pair( rsd, rsd, ii, jj) ) continue;

			if ( atom_is_heavy( rsd, jj ) ) {
				if ( rsd.atom_is_polar_hydrogen( ii ) ) {
					get_atom_atom_geometric_solvation_for_donor( ii, rsd, jj, rsd, pose, energy, update_deriv, deriv, hbe_type );
					fill_atom_derivs_for_donor( deriv, rsd, ii, jj, atom_derivs, atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
				} else if (  rsd.heavyatom_is_an_acceptor( ii ) ) {
					get_atom_atom_geometric_solvation_for_acceptor( ii, rsd, jj, rsd, pose, energy, update_deriv, deriv, hbe_type );
					fill_atom_derivs_for_acceptor( deriv, rsd, ii, jj, atom_derivs, atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
				}
			}
			if ( atom_is_heavy ( rsd, ii ) ) {
				if ( rsd.atom_is_polar_hydrogen( jj ) ) {
					get_atom_atom_geometric_solvation_for_donor( jj, rsd, ii, rsd, pose, energy, update_deriv, deriv, hbe_type );
					fill_atom_derivs_for_donor( deriv, rsd, jj, ii, atom_derivs, atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
				} else if (  rsd.heavyatom_is_an_acceptor( jj ) ) {
					get_atom_atom_geometric_solvation_for_acceptor( jj, rsd, ii, rsd, pose, energy, update_deriv, deriv, hbe_type );
					fill_atom_derivs_for_acceptor( deriv, rsd, jj, ii, atom_derivs, atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergyEvaluator::eval_residue_pair_derivatives(
	conformation::Residue const & ires,
	conformation::Residue const & jres,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	Real const geom_sol_weight,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{


	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( geom_solv_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	Real energy( 0.0 );
	hbonds::HBondDerivs deriv;
	hbonds::HBEvalTuple hbe_type;
	static bool const update_deriv( true );
	Real const weighted_energy = -1.0 * geom_sol_weight;

	for ( Size k = 1, kend = neighbs.size(); k <= kend; ++k ) {
		Size const ii = neighbs[ k ].atomno1();
		Size const jj = neighbs[ k ].atomno2();

		if ( atom_is_heavy( jres, jj ) ) {
			if ( ires.atom_is_polar_hydrogen( ii ) ) {
				get_atom_atom_geometric_solvation_for_donor( ii, ires, jj, jres, pose, energy, update_deriv, deriv, hbe_type );
				fill_atom_derivs_for_donor( deriv, ires, ii, jj, r1_atom_derivs, r2_atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );

			} else if (  ires.heavyatom_is_an_acceptor( ii ) ) {
				get_atom_atom_geometric_solvation_for_acceptor( ii, ires, jj, jres, pose, energy, update_deriv, deriv, hbe_type );
				fill_atom_derivs_for_acceptor( deriv, ires, ii, jj, r1_atom_derivs, r2_atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
			}
		}
		if ( atom_is_heavy ( ires, ii ) ) {
			if ( jres.atom_is_polar_hydrogen( jj ) ) {
				get_atom_atom_geometric_solvation_for_donor( jj, jres, ii, ires, pose, energy, update_deriv, deriv, hbe_type );
				fill_atom_derivs_for_donor( deriv, jres, jj, ii, r2_atom_derivs, r1_atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
			} else if (  jres.heavyatom_is_an_acceptor( jj ) ) {
				get_atom_atom_geometric_solvation_for_acceptor( jj, jres, ii, ires, pose, energy, update_deriv, deriv, hbe_type );
				fill_atom_derivs_for_acceptor( deriv, jres, jj, ii, r2_atom_derivs, r1_atom_derivs, hbe_type, options_.hbond_options(), weighted_energy );
			}
		}

	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::residue_pair_energy_ext(
											conformation::Residue const & rsd1,
											conformation::Residue const & rsd2,
											ResPairMinimizationData const & min_data,
											pose::Pose const & pose
											) const
{
	Real score( 0.0 );
	Real energy( 0.0 );

	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( geom_solv_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	Size m = 0;
	Size n = 0;

	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		m = neighbs[ ii ].atomno1();
		n = neighbs[ ii ].atomno2();

		if ( atom_is_heavy( rsd2, n ) ) {
			if ( rsd1.atom_is_polar_hydrogen( m ) ) {
				get_atom_atom_geometric_solvation_for_donor( m, rsd1, n, rsd2, pose, energy );
				score += energy;
			} else if (  rsd1.heavyatom_is_an_acceptor( m ) ) {
				get_atom_atom_geometric_solvation_for_acceptor( m, rsd1, n, rsd2, pose, energy );
				score += energy;
			}
		}
		if ( atom_is_heavy ( rsd1, m ) ) {
			if ( rsd2.atom_is_polar_hydrogen( n ) ) {
				get_atom_atom_geometric_solvation_for_donor( n, rsd2, m, rsd1, pose, energy );
				score += energy;
			} else if (  rsd2.heavyatom_is_an_acceptor( n ) ) {
				get_atom_atom_geometric_solvation_for_acceptor( n, rsd2, m, rsd1, pose, energy );
				score += energy;
			}
		}
	}
	return score;
}

//////////////////////////////////////////////////////////////////////////////
// helper functions that are usually in scorefiles, but now in this evaluator
// to prevent copying code between GeometricSolEnergy files.
bool
GeometricSolEnergyEvaluator::defines_score_for_residue_pair(
												   conformation::Residue const & rsd1,
												   conformation::Residue const & rsd2,
												   bool res_moving_wrt_eachother
												   ) const
{
	if ( rsd1.seqpos() == rsd2.seqpos() )	return false;
	return res_moving_wrt_eachother;
}

etable::count_pair::CountPairFunctionCOP
GeometricSolEnergyEvaluator::get_count_pair_function(
											Size const res1,
											Size const res2,
											pose::Pose const & pose
											) const
{
	using namespace etable::count_pair;
	if ( res1 == res2 )	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}

etable::count_pair::CountPairFunctionCOP
GeometricSolEnergyEvaluator::get_count_pair_function(
											conformation::Residue const & rsd1,
											conformation::Residue const & rsd2
											) const
{
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	// PUT THIS IN PROPERLY LATER.
	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		//		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairAll ) );
}

etable::count_pair::CountPairFunctionCOP
GeometricSolEnergyEvaluator::get_intrares_countpair(
			 conformation::Residue const & res
) const
{
	using namespace etable::count_pair;
	return CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_3 );
}


void
GeometricSolEnergyEvaluator::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData & pair_data
) const
{

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( geom_solv_pair_nblist ) ));
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( 5.2 + 2*tolerated_narrow_nblist_motion, 2 );

	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );

	pair_data.set_data( geom_solv_pair_nblist, nblist );
}


//////////////////////////////////////////////////////////////////////
Real
GeometricSolEnergyEvaluator::precalculate_bb_bb_energy_for_design(
 pose::Pose const & pose
) const {

	Real precalculated_bb_bb_energy( 0.0 );
  Size const total_residue = pose.total_residue();

  EnergyGraph const & energy_graph( pose.energies().energy_graph() );

  for (Size i = 1; i <= total_residue; i++ ){

    conformation::Residue const & res_i( pose.residue( i ) );

    for( graph::Graph::EdgeListConstIter
					 iter = energy_graph.get_node( i )->const_edge_list_begin();
				 iter != energy_graph.get_node( i )->const_edge_list_end();
				 ++iter ){

      Size j( (*iter)->get_other_ind( i ) );

      conformation::Residue const & res_j( pose.residue( j ) );

      //only need to do it one way since will sample the reverse when res_i is at j
      precalculated_bb_bb_energy += geometric_sol_one_way_bb_bb(res_i, res_j, pose);

    }
  }

	return precalculated_bb_bb_energy;
}


} //geometric_solvation
} //scoring
} //core
