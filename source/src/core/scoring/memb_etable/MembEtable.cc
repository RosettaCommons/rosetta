// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/scoring/etable/Etable.cc
/// @brief
/// @author

// Unit Headers
#include <core/scoring/memb_etable/MembEtable.hh>
#include <core/scoring/etable/Etable.hh>

#include <basic/options/option.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility Headers

// C++ Headers
#include <iostream>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <fstream>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace ObjexxFCL;

static basic::Tracer TR("core.scoring.etable");

namespace core {
namespace scoring {
namespace etable {

using namespace basic::options;
using namespace basic::options::OptionKeys;

///  constructor
MembEtable::MembEtable(
	chemical::AtomTypeSetCAP atom_set_in, // like etable namespace
	EtableOptions const & options,
	std::string const alternate_parameter_set // = ""
):
  //pba
  Etable( atom_set_in, options, alternate_parameter_set ),

	// from atop_props_in:
	atom_set_                 ( atom_set_in ),
	n_atomtypes               ( atom_set_in->n_atomtypes() ),

	// from options
	max_dis_                  ( options.max_dis ),
	bins_per_A2               ( options.bins_per_A2 ),
	Wradius                   ( options.Wradius ), // global mod to radii
	lj_switch_dis2sigma       ( options.lj_switch_dis2sigma ),
	max_dis2                  ( max_dis_*max_dis_ ),
	etable_disbins            ( static_cast< int >( max_dis2 * bins_per_A2)+1),

	// hard-coded for now
	lj_use_lj_deriv_slope     ( true ),
	lj_slope_intercept        ( 0.0 ),
	lj_use_hbond_radii        ( true ),
	lj_hbond_dis              ( 3.0 ),
	lj_hbond_hdis             ( 1.95 ),
	lj_hbond_accOch_dis       ( 2.80 ), // unused
	lj_hbond_accOch_hdis      ( 1.75 ), // unused
	lj_use_water_radii        ( true ),
	lj_water_dis              ( 3.0 ),
	lj_water_hdis             ( 1.95 ),
	lk_min_dis2sigma          ( 0.89 ),
	min_dis                   ( 0.01 ),
	min_dis2                  ( min_dis * min_dis ),
	add_long_range_damping    ( true ),
	long_range_damping_length ( 0.5 ),
	epsilon                   ( 0.0001 ),
	safe_max_dis2             ( max_dis2 - epsilon ),
	hydrogen_interaction_cutoff2_( option[ score::fa_Hatr ] ?
		std::pow( max_dis_ + 2*chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH, 2 ) :
		std::pow(5.0,2) ),
	nblist_dis2_cutoff_XX_    ((options.max_dis+1.5)          *  (options.max_dis+1.5)          ),
	nblist_dis2_cutoff_XH_    ( option[ score::fa_Hatr ] ?
															nblist_dis2_cutoff_XX_ :
															((options.max_dis+0.5)/2.0+2.2) * ((options.max_dis+0.5)/2.0+2.2) ),
	nblist_dis2_cutoff_HH_    ( option[ score::fa_Hatr ] ?
															nblist_dis2_cutoff_XX_ :
															4.4 * 4.4 ),
	max_non_hydrogen_lj_radius_( 0.0 ),
	max_hydrogen_lj_radius_( 0.0 )

{

	// size the arrays
	//ljatr_.dimension(  etable_disbins, n_atomtypes, n_atomtypes );
	//ljrep_.dimension(  etable_disbins, n_atomtypes, n_atomtypes );
	solv1_.dimension(  etable_disbins, n_atomtypes, n_atomtypes );
	solv2_.dimension(  etable_disbins, n_atomtypes, n_atomtypes );
  memb_solv1_.dimension(  etable_disbins, n_atomtypes, n_atomtypes ); //pba
  memb_solv2_.dimension(  etable_disbins, n_atomtypes, n_atomtypes ); //pba
	//dljatr_.dimension( etable_disbins, n_atomtypes, n_atomtypes );
	//dljrep_.dimension( etable_disbins, n_atomtypes, n_atomtypes );
	dsolv2_.dimension(  etable_disbins, n_atomtypes, n_atomtypes );
	dsolv1_.dimension( etable_disbins, n_atomtypes, n_atomtypes );
  memb_dsolv2_.dimension(  etable_disbins, n_atomtypes, n_atomtypes ); //pba
  memb_dsolv1_.dimension( etable_disbins, n_atomtypes, n_atomtypes ); //pba

	lj_radius_.resize( n_atomtypes, 0.0 );
	//lj_wdepth_.resize( n_atomtypes, 0.0 );
	lk_dgfree_.resize( n_atomtypes, 0.0 );
	lk_lambda_.resize( n_atomtypes, 0.0 );
	lk_volume_.resize( n_atomtypes, 0.0 );
  memb_lk_dgfree_.resize( n_atomtypes, 0.0 );  //pba
  //lk_dgrefce_.resize( n_atomtypes, 0.0 );       //pba
  //memb_lk_dgrefce_.resize( n_atomtypes, 0.0 ); //pba
  lk_dgrefce_.dimension( n_atomtypes );
  memb_lk_dgrefce_.dimension( n_atomtypes );


	for ( int i=1; i<= n_atomtypes; ++i ) {
		lj_radius_[i] = (*atom_set_in)[i].lj_radius();
		//lj_wdepth_[i] = (*atom_set_in)[i].lj_wdepth();
		//lk_dgfree_[i] = (*atom_set_in)[i].lk_dgfree();
		lk_lambda_[i] = (*atom_set_in)[i].lk_lambda();
		lk_volume_[i] = (*atom_set_in)[i].lk_volume();
    //memb_lk_dgfree_[i] = (*atom_set_in)[i].memb_lk_dgfree();   //pba
    //lk_dgrefce_[i] = (*atom_set_in)[i].lk_dgrefce();             //pba
    //memb_lk_dgrefce_[i] = (*atom_set_in)[i].memb_lk_dgrefce(); //pba
    //lk_dgrefce_(i) = (*atom_set_in)[i].lk_dgrefce();             //pba
    //memb_lk_dgrefce_(i) = (*atom_set_in)[i].memb_lk_dgrefce(); //pba


		if ( (*atom_set_in)[i].is_hydrogen() ) {
			if ( lj_radius_[i] > max_hydrogen_lj_radius_ ) max_hydrogen_lj_radius_ = lj_radius_[i];
		} else {
			if ( lj_radius_[i] > max_non_hydrogen_lj_radius_ ) max_non_hydrogen_lj_radius_ = lj_radius_[i];
		}
	}

	//if ( alternate_parameter_set.size() ) {
		/// uses alternate paramers
		std::string param_name;
   /*
		param_name = "LJ_RADIUS_"+alternate_parameter_set;
		if ( atom_set_in->has_extra_parameter( param_name ) ) {
			TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
			Size const index( atom_set_in->extra_parameter_index( param_name ) );
			for ( int i=1; i<= n_atomtypes; ++i ) lj_radius_[i] = (*atom_set_in)[i].extra_parameter( index );
		}

		param_name = "LJ_WDEPTH_"+alternate_parameter_set;
		if ( atom_set_in->has_extra_parameter( param_name ) ) {
			TR << "Using alternate parameters: " << param_name << " in Etable construction."<<std::endl;
			Size const index( atom_set_in->extra_parameter_index( param_name ) );
			for ( int i=1; i<= n_atomtypes; ++i ) lj_wdepth_[i] = (*atom_set_in)[i].extra_parameter( index );
		}

		param_name = "LK_DGFREE_"+alternate_parameter_set;
		if ( atom_set_in->has_extra_parameter( param_name ) ) {
			TR << "Using alternate parameters: " << param_name << " in Etable construction."<<std::endl;
			Size const index( atom_set_in->extra_parameter_index( param_name ) );
			for ( int i=1; i<= n_atomtypes; ++i ) lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( index );
		}

		param_name = "LK_LAMBDA_"+alternate_parameter_set;
		if ( atom_set_in->has_extra_parameter( param_name ) ) {
			TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
			Size const index( atom_set_in->extra_parameter_index( param_name ) );
			for ( int i=1; i<= n_atomtypes; ++i ) lk_lambda_[i] = (*atom_set_in)[i].extra_parameter( index );
		}

		param_name = "LK_VOLUME_"+alternate_parameter_set;
		if ( atom_set_in->has_extra_parameter( param_name ) ) {
			TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
			Size const index( atom_set_in->extra_parameter_index( param_name ) );
			for ( int i=1; i<= n_atomtypes; ++i ) lk_volume_[i] = (*atom_set_in)[i].extra_parameter( index );
		}
    */
    //param_name = "WAT_LK_DGFREE_"+alternate_parameter_set; //pba
    param_name = "LK_DGFREE"; //pba
    //if ( atom_set_in->has_extra_parameter( param_name ) ) {
      TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
      Size const lkfree_index( atom_set_in->extra_parameter_index( param_name ) );
      for ( int i=1; i<= n_atomtypes; ++i ) lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( lkfree_index );
    //}

    //param_name = "MEMB_LK_DGFREE_"+alternate_parameter_set; //pba
    param_name = "MEMB_LK_DGFREE"; //pba
    //if ( atom_set_in->has_extra_parameter( param_name ) ) {
      TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
      Size const mb_lkfree_index( atom_set_in->extra_parameter_index( param_name ) );
      for ( int i=1; i<= n_atomtypes; ++i ) memb_lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( mb_lkfree_index );
    //}

    //param_name = "LK_DGREFCE_"+alternate_parameter_set; //pba
    param_name = "LK_DGREFCE"; //pba
    //if ( atom_set_in->has_extra_parameter( param_name ) ) {
      TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
      Size const lkref_index( atom_set_in->extra_parameter_index( param_name ) );
      //for ( int i=1; i<= n_atomtypes; ++i ) lk_dgrefce_[i] = (*atom_set_in)[i].extra_parameter( index );
      for ( int i=1; i<= n_atomtypes; ++i ) lk_dgrefce_(i) = (*atom_set_in)[i].extra_parameter( lkref_index );
    //}

    //param_name = "MEMB_LK_DGREFCE_"+alternate_parameter_set; //pba
    param_name = "MEMB_LK_DGREFCE"; //pba
    //if ( atom_set_in->has_extra_parameter( param_name ) ) {
      TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
      Size const mb_lkref_index( atom_set_in->extra_parameter_index( param_name ) );
      //for ( int i=1; i<= n_atomtypes; ++i ) memb_lk_dgrefce_[i] = (*atom_set_in)[i].extra_parameter( index );
      for ( int i=1; i<= n_atomtypes; ++i ) memb_lk_dgrefce_(i) = (*atom_set_in)[i].extra_parameter( mb_lkref_index );
    //}
	//} if alternate param set


    //pbadebug
  /*for ( int i=1; i<= n_atomtypes; ++i ) {
     std::cout << "atom " << i << " lj_rad " << lj_radius_[i] << " lk_vol " << lk_volume_[i] <<
                 " lk_free " << lk_dgfree_[i] << " lk_ref " << lk_dgrefce_(i) << " lambda " << lk_lambda_[i] <<
                 " mb_lk_ref " << memb_lk_dgrefce_(i) << " mb_lk_free " << memb_lk_dgfree_[i] << "\n";
  }*/

//	Real const MAX_H_HEAVY_DISTANCE = 1.35; // FIX THIS SULFUR to hydrogen bond length.

/*  Real max_lj_rep_for_h = std::max(
		2 * max_hydrogen_lj_radius_ + 2 * MAX_H_HEAVY_DISTANCE,
		max_hydrogen_lj_radius_ + MAX_H_HEAVY_DISTANCE + max_non_hydrogen_lj_radius_ );


	hydrogen_interaction_cutoff2_ = basic::options::option[ score::fa_Hatr ] ?
		std::pow( 2 * MAX_H_HEAVY_DISTANCE + max_dis_, 2 ) : max_lj_rep_for_h * max_lj_rep_for_h;
*/
	make_pairenergy_table();
}

////////////////////////////////////////////////////////////////////////////////
/*void
MembEtable::copy_from( Etable const * source )
{
  assert( dynamic_cast< TenANeighborNode const * > (source) );
  assert( static_cast< TenANeighborNode const * > (source) );
  TenANeighborNode const * tAsource( static_cast< TenANeighborNode const * > (source) );
  neighbor_mass_           = tAsource->neighbor_mass_;
  sum_of_neighbors_masses_ = tAsource->sum_of_neighbors_masses_;
  since_last_sonm_update_  = tAsource->since_last_sonm_update_;
}*/

////////////////////////////////////////////////////////////////////////////////
/// @begin MembEtable::make_pairenergy_table
///
/// @brief calculate fast lookup arrays for vdw and solvation energy
///
/// @detailed
///
/// Several energies are precomputed when fullatom mode is initialized and
/// stored in lookup tables to speed fullatom calculation. Currently
/// pre-computed values are the Lennard-Jones van der Waals approximation
/// (lj) and the Lazaridis-Karplus implicit solvation function (lk). For
/// each of these energies the derivative w.r.t. atom pair separation distance is
/// calculated and stored as well. Note that the lj energy is artificially
/// divided into atractive and repulsive components.
///
/// @global_read
///
/// pdbstatistics_pack.h: several threshold distances
/// energy.h: just about everything should be used in the tree of
///           etable functions
///
/// @global_write
///
/// the etable: ljatr,dljatr,ljrep,dljrep,solvE
///
/// @remarks
///
/// @references
///
/// @authors ctsa 10-2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::make_pairenergy_table()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//  locals
	Real dis2,dis,dis2_step;
	//Real atrE,d_atrE,repE,d_repE;
  Real solvE1,solvE2,dsolvE1,dsolvE2;
  Real memb_solvE1,memb_solvE2,memb_dsolvE1,memb_dsolvE2; //pba


	//  pre-calculated coefficients for high-speed potential calculation
	//    computed in initialization function
	//  Note: +1 atom type is a disulfide-bonded cysteine
	//        (only if you are using the old disulfide hack)
	FArray2D< Real > lj_sigma( n_atomtypes + 1, n_atomtypes + 1 );
	//FArray2D< Real > lj_r6_coeff( n_atomtypes + 1, n_atomtypes + 1 );
	//FArray2D< Real > lj_r12_coeff( n_atomtypes + 1, n_atomtypes + 1 );
	//FArray2D< Real > lj_switch_intercept( n_atomtypes + 1, n_atomtypes + 1 );
	//FArray2D< Real > lj_switch_slope( n_atomtypes + 1, n_atomtypes + 1 );
	FArray1D< Real > lk_inv_lambda2( n_atomtypes );
	FArray2D< Real > lk_coeff( n_atomtypes, n_atomtypes );
	FArray2D< Real > lk_min_dis2sigma_value( n_atomtypes + 1, n_atomtypes + 1 );
  FArray2D< Real > memb_lk_coeff( n_atomtypes, n_atomtypes ); //pba
  FArray2D< Real > memb_lk_min_dis2sigma_value( n_atomtypes + 1, n_atomtypes + 1 ); //pba

	// parameters to calculate the damping at max_dis range
	Real damping_thresh_dis2;
	int damping_disbins, normal_disbins;
	//Real dljatr_damp, dljrep_damp, dsolv1_damp, dsolv2_damp;
	//Real intercept_ljatr_damp, intercept_ljrep_damp;
	Real dsolv1_damp, dsolv2_damp, intercept_solv1_damp, intercept_solv2_damp;
  Real memb_dsolv1_damp, memb_dsolv2_damp, intercept_memb_solv1_damp, intercept_memb_solv2_damp; //pba

	//  etable parameters
	TR << "Starting membrane specific energy table calculation" << std::endl;

// 	std::TR << "Energy table parameter set: " << etable_label << ' ' <<
// 	 etable_revision.substr( 1, etable_revision.length() - 3 ) << ' ' <<
// 	 etable_date.substr( 1, etable_date.length() - 3 ) << std::endl;

	// ctsa - values precomputed from etable parameters:
	//   these were all originally constants but need
	//   to be calculated here b/c the switch value:
	//   lj_switch_dis2sigma
	//   is set at runtime
	//lj_switch_sigma2dis = 1.0/lj_switch_dis2sigma;

	// ctsa - value of the lennard-jones potential at the linear
	//   switch point divided by the wdepth (note that this
	//   coefficient is independent of atomtype)
	//lj_switch_value2wdepth = std::pow( lj_switch_sigma2dis, 12 ) -
	// 2.0 * std::pow( lj_switch_sigma2dis, 6 );

	// ctsa - slope of the lennard-jones potential at the linear
	//   switch point times sigma divided by wdepth (note that this
	//   coefficient is independent of atomtype)
	//lj_switch_slope_sigma2wdepth = -12.0 * (
	// std::pow( lj_switch_sigma2dis, 13 ) -
	// std::pow( lj_switch_sigma2dis, 7 ) );

	//  initialize non-distance dependent coefficients
	/*precalc_etable_coefficients(lj_sigma, lj_r6_coeff, lj_r12_coeff,
	 lj_switch_intercept, lj_switch_slope, lk_inv_lambda2, lk_coeff,
	 memb_lk_coeff, lk_min_dis2sigma_value, memb_lk_min_dis2sigma_value );*/

  precalc_etable_coefficients(lj_sigma,lk_inv_lambda2,lk_coeff,memb_lk_coeff,lk_min_dis2sigma_value,memb_lk_min_dis2sigma_value); //pba

	//  calc distance**2 step per bin
	dis2_step = 1.0 / bins_per_A2;

	//  get the number of damping disbins
	if ( add_long_range_damping ) {
		Real const dif = max_dis_ - long_range_damping_length;
		damping_thresh_dis2 = max_dis2 - ( dif * dif );
		damping_disbins = static_cast< int >( damping_thresh_dis2*bins_per_A2 );
		normal_disbins = etable_disbins-damping_disbins;
	} else {
		normal_disbins = etable_disbins;
	}

	//  ctsa - step through distance**2 bins and calculate potential
	for ( int atype1 = 1, atype_end = n_atomtypes; atype1 <= atype_end; ++atype1 ) {
		for ( int atype2 = 1; atype2 <= atype_end; ++atype2 ) {
			//  ctsa - normal bins have their lj and lk values
			//    calculated analytically
			for ( int disbin = 1; disbin <= normal_disbins; ++disbin ) {
				dis2 = ( disbin - 1 ) * dis2_step;

			  /*calc_etable_value(dis2,atype1,atype2,atrE,d_atrE,repE,d_repE,solvE1,
				 solvE2,dsolvE1,dsolvE2,lj_sigma, lj_r6_coeff,lj_r12_coeff,
				 lj_switch_intercept,lj_switch_slope, lk_inv_lambda2,
				 lk_coeff, lk_min_dis2sigma_value, memb_solvE1, memb_solvE2,
         memb_lk_coeff, memb_lk_min_dis2sigma_value, memb_dsolvE); //pba */

        calc_etable_value(dis2,atype1,atype2,solvE1,solvE2,dsolvE1,dsolvE2,lj_sigma,lk_inv_lambda2,
         lk_coeff, lk_min_dis2sigma_value,memb_solvE1, memb_solvE2,
         memb_lk_coeff, memb_lk_min_dis2sigma_value, memb_dsolvE1, memb_dsolvE2); //pba

				//ljatr_(disbin,atype2,atype1) = atrE;
				//dljatr_(disbin,atype2,atype1) = d_atrE;
				//ljrep_(disbin,atype2,atype1) = repE;
				//dljrep_(disbin,atype2,atype1) = d_repE;
				solv1_(disbin,atype2,atype1) = solvE1;
				solv2_(disbin,atype2,atype1) = solvE2;
				dsolv2_(disbin,atype2,atype1) = dsolvE2;
				dsolv1_(disbin,atype2,atype1) = dsolvE1;
        memb_solv1_(disbin,atype2,atype1) = memb_solvE1; //pba
        memb_solv2_(disbin,atype2,atype1) = memb_solvE2; //pba
        memb_dsolv2_(disbin,atype2,atype1) = memb_dsolvE2; //pba
        memb_dsolv1_(disbin,atype2,atype1) = memb_dsolvE1; //pba
    //pbadebug WARNING
    //if(atype1==18 && atype2==19)
    //  std::cout << "bin s1 s2 " << disbin << " " << solv1_(disbin,atype2,atype1) << " " << solv2_(disbin,atype2,atype1) << std::endl;

			}

			if ( add_long_range_damping ) {

				// ctsa - remaining bins damp to 0. on a linear path
				/*dljatr_damp = -ljatr_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;
				dljrep_damp = -ljrep_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;*/
				dsolv1_damp = -solv1_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;
				dsolv2_damp = -solv2_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;
				//pba
				memb_dsolv1_damp = -memb_solv1_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;
				memb_dsolv2_damp = -memb_solv2_(normal_disbins,atype2,atype1) /
				 long_range_damping_length;

				/*intercept_ljatr_damp = -dljatr_damp*max_dis_;
				intercept_ljrep_damp = -dljrep_damp*max_dis_;*/
				intercept_solv1_damp = -dsolv1_damp*max_dis_;
				intercept_solv2_damp = -dsolv2_damp*max_dis_;
				//pba
				intercept_memb_solv1_damp = -memb_dsolv1_damp*max_dis_;
				intercept_memb_solv2_damp = -memb_dsolv2_damp*max_dis_;

				for ( int disbin = normal_disbins+1; disbin <= etable_disbins; ++disbin ) {
					dis2 = ( disbin - 1 ) * dis2_step;
					dis = std::sqrt(dis2);

					/*ljatr_(disbin,atype2,atype1) = intercept_ljatr_damp + dis *dljatr_damp;
					ljrep_(disbin,atype2,atype1) = intercept_ljrep_damp + dis *dljrep_damp;*/
					solv1_(disbin,atype2,atype1) = intercept_solv1_damp + dis *dsolv1_damp;
					solv2_(disbin,atype2,atype1) = intercept_solv2_damp + dis *dsolv2_damp;
					//pba
					memb_solv1_(disbin,atype2,atype1) = intercept_memb_solv1_damp + dis * memb_dsolv1_damp;
					memb_solv2_(disbin,atype2,atype1) = intercept_memb_solv2_damp + dis * memb_dsolv2_damp;

					/*dljatr_(disbin,atype2,atype1) = dljatr_damp;
					dljrep_(disbin,atype2,atype1) = dljrep_damp;*/
					dsolv1_(disbin,atype2,atype1)  = dsolv1_damp;
					dsolv2_(disbin,atype2,atype1) = dsolv2_damp;
					//pba
					memb_dsolv2_(disbin,atype2,atype1) = memb_dsolv2_damp;
					memb_dsolv1_(disbin,atype2,atype1) = memb_dsolv1_damp;

					//pbadebug WARNING
					//if(atype1==18 && atype2==19)
					//  std::cout << "bin s1 s2 " << disbin << " " << solv1_(disbin,atype2,atype1) << " " << solv2_(disbin,atype2,atype1) << std::endl;
				}
			}

			//  ctsa - set last bin of all values to zero
			/*ljatr_(etable_disbins,atype2,atype1) = 0.0;
			ljrep_(etable_disbins,atype2,atype1) = 0.0;*/
			solv1_(etable_disbins,atype2,atype1) = 0.0;
			solv2_(etable_disbins,atype2,atype1) = 0.0;
      memb_solv1_(etable_disbins,atype2,atype1) = 0.0; //pba
      memb_solv2_(etable_disbins,atype2,atype1) = 0.0; //pba
		}
	}

	//db  the following function call modifies the potential in three ways:
	//db     (1) the solvation energy for nonpolar atoms is held constant below
	//db     4.2A to avoid shifting the minimum in the LJ potential.
	//db     (2) a short range repulsion is added between backbone oxygens which are
	//db     otherwise brought too close together by the LJatr.
	//db     (3) the range of the repulsive interaction between non polar hydrogens is
	//db     increased slightly.  (this is currently commented out because the effects
	//db     on design have not been tested)

	//db  all three modifications are based on inspection of the atom pair distributions
	//db  after extensive refinement.
  //pbadebug WARNING
	modify_pot();

	if( !option[ score::no_smooth_etables ] ) smooth_etables();

	// sheffler changed this to happen after modify_pot so that
	// etable smoothing of Hydrogens and Waters will work properly
	//if ( !option[ score::fa_Hatr ] ) zero_hydrogen_and_water_ljatr();


	////////////////////////////////////////////////////////
	// etable I/O stuff
	////////////////////////////////////////////////////////
	using namespace std;
	using namespace ObjexxFCL;

	// just for convenience in input/output_etables
	map< string, FArray3D<Real>* > etables;
	/*etables[ "ljatr"] = &  ljatr_;	etables[ "ljrep"] = &  ljrep_;
	etables["dljatr"] = & dljatr_;	etables["dljrep"] = & dljrep_;*/
	etables[ "solv1"] = &  solv1_;	etables[ "solv2"] = &  solv2_;
	etables["dsolv1"]  = & dsolv1_;
  etables["dsolv2"]  = & dsolv2_;
  etables[ "memb_solv1"] = &  memb_solv1_;  etables[ "memb_solv2"] = &  memb_solv2_; //pba
  etables["memb_dsolv2"]  = & memb_dsolv2_; etables["memb_dsolv1"]  = & memb_dsolv1_;  //pba

	if ( option[ score::input_etables ].user() ) {
		string tag = option[ score::input_etables ];
		TR << "INPUT ETABLES " << tag << std::endl;
		for (map<string,FArray3D<Real>*>::iterator i = etables.begin(); i != etables.end(); i++) {
			string ename = i->first;
			string fname = tag+"."+ename+".etable";
			std::ifstream input( fname.c_str() ); // TODO sheffler: figure out how to do this the right way
			input_etable(*(i->second),ename,input);
			input.close();
		}
	}

	if ( option[ score::output_etables ].user() ) {
		string header = option[ score::output_etables ];
		TR << "OUTPUT ETABLES " << header << std::endl;
		for (map<string,FArray3D<Real>*>::iterator i = etables.begin(); i != etables.end(); i++) {
			string ename = i->first;
			string fname = header+"."+ename+".etable";
			TR << "output_etable: writing etable: " << ename << " to " << fname << std::endl;
			ofstream out(fname.c_str());
			output_etable(*(i->second),ename,out);
			out.close();
		}
	}

	TR << "Finished calculating membrane specific energy tables." << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin modify_pot
///
/// @brief modify Etable to better treat 0-0, C-C, and H-H interactions
///
/// @detailed
///$$$ the Etables are  modified in three ways:
///$$$ (1) the LK solvation energy is set to a constant below 4.2A to avoid shifting the position
///$$$ of the minimum on the LJatr potential.  in refined decoys the peak in the C C pair
///$$$ distribution function shifts to around 3.8 A from ~4.0A in native structures; the LJatr
///$$$ potential has a minimum at 4.0A but this shifts towards smaller values because of the LK
///$$$ solvation term, which became increasingly favorable at shorter distances
///$$$ (2) the backbone carbonyl oxygen-carbonyl oxygen LJrep term has been modified to become
///$$$ moderately repulsive at distances less than 3.6A.  this is to counteract the favorable
///$$$ LJatr between the atoms (which have radii of ~1.4A and so a minimum at ~2.8A; very few
///$$$ counts are observed in the pdb until around 3.2A) which leads to a significant shift in the
///$$$ O O pair distribution function towards smaller values in refined decoys.  the repulsion is
///$$$ a temporary proxy for the lack of explicit electrostatic repulsion in the current force
///$$$ field.
///$$$ (3) a third soft repulsion between non polar hydrogens that was also based on comparison of
///$$$ refined decoy to native pdf's is currently commented out as the effects on packing have not
///$$$ been tested.  it was observed that the protons tend to pile up at just beyond the point
///$$$ where the repulsion becomes strong, perhaps due to a general tendency to overcontraction
///$$$ because of long range LJatr interactions not compensated by interactions with solvent
///$$$ (which are of course missing)
///
/// @global_read
/// pdbstatistics_pack.h:  ljatr,dljatr,ljrep, dljrep, solv1,solv2,dsolv
///
/// @global_write
/// pdbstatistics_pack.h:  ljatr,dljatr,ljrep, dljrep, solv1,solv2,dsolv
///
/// @remarks
///
/// @references
///
/// @authors ctsa 10-2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::modify_pot()
{

	using namespace std;

	bool mod_hhrep      = false; //truefalseoption("mod_hhrep");
	Real const h = 0.0;// fullatom_setup_ns::mod_hhrep_height;
	Real const c = 0.0;// fullatom_setup_ns::mod_hhrep_center;
	Real const w = 0.0;// fullatom_setup_ns::mod_hhrep_width;
	Real const e = 0.0;// fullatom_setup_ns::mod_hhrep_exponent;


	if(mod_hhrep) {
		TR << "fullatom_setup: modifying h-h repulsion "
							<< " hhrep center "   << c
							<< " hhrep height "   << h
							<< " hhrep width "    << w
							<< " hhrep exponent " << e
							<< std::endl;
	}

	{

			Real const bin = ( 4.2 * 4.2 / .05 ) + 1.0; //SGM Off-by-1 bug fix: Add + 1.0: Index 1 is at distance^2==0
		int const ibin( static_cast< int >( bin ) );
 		for ( int k = 1; k <= etable_disbins; ++k ) {
 			Real const dis = std::sqrt( ( k - 1 ) * .05f ); //SGM Off-by-1 bug fix: k -> ( k - 1 )


			// if( !SMOOTH_ETABLES ) {
				utility::vector1<int> carbon_types;
				carbon_types.push_back( atom_set_->atom_type_index("CH1") );
				carbon_types.push_back( atom_set_->atom_type_index("CH2") );
				carbon_types.push_back( atom_set_->atom_type_index("CH3") );
				carbon_types.push_back( atom_set_->atom_type_index("aroC") );
				if ( dis < 4.2 ) {
					for ( int i = 1, i_end = carbon_types.size(); i <= i_end; ++i ) {
						for ( int j = 1, j_end = carbon_types.size(); j <= j_end; ++j ) {
							int const ii = carbon_types[i];
							int const jj = carbon_types[j];
							solv1_(k,jj,ii) = solv1_(ibin,jj,ii);
							solv1_(k,ii,jj) = solv1_(ibin,ii,jj); // Why is this duplicated?
							solv2_(k,jj,ii) = solv2_(ibin,jj,ii);
							solv2_(k,ii,jj) = solv2_(ibin,ii,jj); // Why is this duplicated?
							dsolv2_(k,jj,ii) = 0.0;
							dsolv2_(k,ii,jj) = 0.0; // Why is this duplicated?
							dsolv1_(k,ii,jj) = 0.0;
							dsolv1_(k,jj,ii) = 0.0; // Why is this duplicated?
              memb_solv1_(k,jj,ii) = memb_solv1_(ibin,jj,ii);
              memb_solv1_(k,ii,jj) = memb_solv1_(ibin,ii,jj); // Why is this duplicated?
              memb_solv2_(k,jj,ii) = memb_solv2_(ibin,jj,ii);
              memb_solv2_(k,ii,jj) = memb_solv2_(ibin,ii,jj); // Why is this duplicated?
              memb_dsolv2_(k,jj,ii) = 0.0;
              memb_dsolv2_(k,ii,jj) = 0.0; // Why is this duplicated?
              memb_dsolv1_(k,ii,jj) = 0.0;
              memb_dsolv1_(k,jj,ii) = 0.0; // Why is this duplicated?
						}
					}
				}
			// }

//   push carbonyl oxygens (in beta sheets) apart.  a proxy for the missing
//   electrostatic repulsion needed to counteract the LJatr which pulls the oxyens
//   together
			//mjo commenting out 'OCbb_idx' because it is unused and causes a warning
			//int const OCbb_idx = atom_set_->atom_type_index("OCbb");
			if ( dis <= 3.6 ) {
				//mjo commenting out 'fac' because it is unused and causes a warning
				//Real const fac = std::max( dis - 3.6, -1.5 );
// 				Real const fac = std::max( dis - 3.6f, -1.5f );
//				ljrep_(k,OCbb_idx,OCbb_idx) += 2 * ( fac * fac );
//				dljrep_(k,OCbb_idx,OCbb_idx) += 4 * fac;
			}
//  the following gives peak at 2.4 in 22 and 23 interactions. maybe push out a
//  bit further.  (this is commented out because effects on design have not been
//  tested)

			// use one half of a single term polynomial
			// as the repulsive term for apolar hydrogens (types 23 & 24)

// 			if( mod_hhrep ) {
// 				if( dis < c ) {
// 					for ( int j = 23; j <= 24; ++j ) {
// 						for ( int kk = 23; kk <= 24; ++kk ) {
// 							ljrep_(k,j,kk)  = h * pow( min(0.0f, (dis-c)/w ), e );//  +
// 							dljrep_(k,j,kk) = h * e / w * pow( min(0.0f, dis-c)/w, e-1 ) ;//  +
// 						}
// 					}
// 				}
// 			}

		} // end  for ( int k = 1; k <= etable_disbins; ++k ) {


	} //Objexx:SGM Extra {} scope is a VC++ work-around


}

void
MembEtable::smooth_etables()
{
	using namespace numeric::interpolation;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	FArray1D<Real> dis(etable_disbins);
	for(int i = 1; i <= etable_disbins; ++i) {
		dis(i) = sqrt( ( i - 1 ) * 1.0 / bins_per_A2 );
	}
	FArray1D<int> bin_of_dis(100);
  for(int i = 1; i <= 100; ++i) {
    float d = ((float)i)/10.0;
    bin_of_dis(i) = (int)( d*d * bins_per_A2 + 1 );
  }
  /*
	//////////////////////////////////////////////////////////////////
	// 1) change ljatr / ljrep split so that scaling is smooth
	//////////////////////////////////////////////////////////////////
	TR << "smooth_etable: changing atr/rep split to bottom of energy well" << std::endl;
	for( int at1 = 1; at1 <= n_atomtypes; ++at1 ) {
		for( int at2 = 1; at2 <= n_atomtypes; ++at2 ) {

			// find ljatr min
			core::Real min_atr = 0.0;
			int which_min = -1;
			for(int i = 1; i <= etable_disbins; ++i) {
				if ( ljatr_(i,at1,at2) < min_atr ) {
					min_atr = ljatr_(i,at1,at2);
					which_min = i;
				}
			}
			// for dis below ljatr min, transfer ljatr to ljrep
			if( min_atr < 0.0 ) {
				for( int i = 1; i <= which_min; ++i ) {
					 ljrep_(i,at1,at2) += ljatr_(i,at1,at2) - min_atr;
					 ljatr_(i,at1,at2) -= ljatr_(i,at1,at2) - min_atr; // = min_atr;
					dljrep_(i,at1,at2) += dljatr_(i,at1,at2);
					dljatr_(i,at1,at2) -= dljatr_(i,at1,at2); // = 0;
				}
			}

		}
	}

	using namespace numeric::interpolation::spline;
	///////////////////////////////////////////////////////////////////
	// 2) spline smooth ljatr/ljrep at fa_max_dis cut
	//////////////////////////////////////////////////////////////////
 	TR << "smooth_etable: spline smoothing lj etables (maxdis = " << max_dis_ << ")" << std::endl;
 	for( int at1 = 1; at1 <= n_atomtypes; ++at1 ) {
 		for( int at2 = 1; at2 <= n_atomtypes; ++at2 ) {

 			int start = bin_of_dis( (int)((max_dis_-1.5)*10.0) );//arg_max_first(dljatr_(1,at1,at2),1,600);

			Real lbx  =     dis(start);
			Real lby  =  ljatr_(start,at1,at2);
			Real lbdy = dljatr_(start,at1,at2);
			Real ubx  = dis(etable_disbins);
			Real uby  = 0.0;
			Real ubdy = 0.0;
			SplineGenerator gen( lbx, lby, lbdy, ubx, uby, ubdy );

		 	if( option[ score::etable_lr ].user() ) {
				Real modx = option[ score::etable_lr ]();
				Real mody = lbx * 0.5;
				gen.add_known_value( modx, mody );
			}

			InterpolatorOP interp( gen.get_interpolator() );
			for( int i = start; i <= etable_disbins; ++i ) {
				interp->interpolate( dis(i), ljatr_(i,at1,at2), dljatr_(i,at1,at2) );
			}

 		}
 	}
  */
	/////////////////////////////////////////////////////////////////////////////
	// 2) spline smooth solv1/solv2 at fa_max_dis cut && and first switchpoint
	/////////////////////////////////////////////////////////////////////////////
  using namespace numeric::interpolation::spline;
	TR << "smooth_etable: spline smoothing solvation etables (max_dis = " << max_dis_ << ")" << std::endl;
	for( int at1 = 1; at1 <= n_atomtypes; ++at1 ) {
		for( int at2 = 1; at2 <= n_atomtypes; ++at2 ) {

			int SWTCH = 1;
			for( SWTCH=1; SWTCH <= etable_disbins; ++SWTCH) {
				// std::cerr << SWTCH << "," << at1 << "," << at2 << "," << solv1_(SWTCH,at1,at2) << " ";
				if( (solv1_(SWTCH,at1,at2) != solv1_(1,at1,at2)) ||
						(solv2_(SWTCH,at1,at2) != solv2_(1,at1,at2)) ||
            (memb_solv1_(SWTCH,at1,at2) != memb_solv1_(1,at1,at2)) ||
            (memb_solv2_(SWTCH,at1,at2) != memb_solv2_(1,at1,at2)) ) {
					break;
				}
			}
			if ( SWTCH > etable_disbins ) continue;

			int const S1 = std::max(1,SWTCH - 30);
			int const E1 = std::min(SWTCH + 20,406);
			int const S2 = bin_of_dis( (int)((max_dis_-1.5)*10.0) );
			int const E2 = etable_disbins;
			// std::cerr << "smooth solv " << S1 << " " << E1 << " " << S2 << " " << E2 << std::endl;

			Real dsolv1e1 = (solv1_(E1+1,at1,at2)-solv1_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real dsolv1s2 = (solv1_(S2  ,at1,at2)-solv1_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));
			Real dsolv2e1 = (solv2_(E1+1,at1,at2)-solv2_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real dsolv2s2 = (solv2_(S2  ,at1,at2)-solv2_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));

      Real mdsolv1e1 = (memb_solv1_(E1+1,at1,at2)-memb_solv1_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
      Real mdsolv1s2 = (memb_solv1_(S2  ,at1,at2)-memb_solv1_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));
      Real mdsolv2e1 = (memb_solv2_(E1+1,at1,at2)-memb_solv2_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
      Real mdsolv2s2 = (memb_solv2_(S2  ,at1,at2)-memb_solv2_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));

			SplineGenerator gen11( dis(S1), solv1_(S1,at1,at2), 0.0     ,	dis(E1), solv1_(E1,at1,at2), dsolv1e1 );
			SplineGenerator gen21( dis(S1), solv2_(S1,at1,at2), 0.0     ,	dis(E1), solv2_(E1,at1,at2), dsolv2e1 );
			SplineGenerator gen12( dis(S2), solv1_(S2,at1,at2), dsolv1s2,	dis(E2), solv1_(E2,at1,at2), 0.0      );
			SplineGenerator gen22( dis(S2), solv2_(S2,at1,at2), dsolv2s2,	dis(E2), solv2_(E2,at1,at2), 0.0      );

      SplineGenerator mgen11( dis(S1), memb_solv1_(S1,at1,at2), 0.0     , dis(E1), memb_solv1_(E1,at1,at2), mdsolv1e1 );
      SplineGenerator mgen21( dis(S1), memb_solv2_(S1,at1,at2), 0.0     , dis(E1), memb_solv2_(E1,at1,at2), mdsolv2e1 );
      SplineGenerator mgen12( dis(S2), memb_solv1_(S2,at1,at2), mdsolv1s2, dis(E2), memb_solv1_(E2,at1,at2), 0.0      );
      SplineGenerator mgen22( dis(S2), memb_solv2_(S2,at1,at2), mdsolv2s2, dis(E2), memb_solv2_(E2,at1,at2), 0.0      );

			InterpolatorOP interp11( gen11.get_interpolator() );
			InterpolatorOP interp21( gen21.get_interpolator() );

      InterpolatorOP minterp11( mgen11.get_interpolator() );
      InterpolatorOP minterp21( mgen21.get_interpolator() );

			for( int i = S1; i <= E1; ++i ) {
				Real d1,d2;
				interp11->interpolate( dis(i), solv1_(i,at1,at2), d1 );
				interp21->interpolate( dis(i), solv2_(i,at1,at2), d2 );
				dsolv2_(i,at1,at2) = d2;
				dsolv1_(i,at1,at2) = d1;
        minterp11->interpolate( dis(i), memb_solv1_(i,at1,at2), d1 );
        minterp21->interpolate( dis(i), memb_solv2_(i,at1,at2), d2 );
        memb_dsolv2_(i,at1,at2) = d2;
        memb_dsolv1_(i,at1,at2) = d1;
			}

			InterpolatorOP interp12( gen12.get_interpolator() );
			InterpolatorOP interp22( gen22.get_interpolator() );

      InterpolatorOP minterp12( mgen12.get_interpolator() );
      InterpolatorOP minterp22( mgen22.get_interpolator() );

			for( int i = S2; i <= E2; ++i ) {
				Real d1,d2;
				interp12->interpolate( dis(i), solv1_(i,at1,at2), d1 );
				interp22->interpolate( dis(i), solv2_(i,at1,at2), d2 );
				dsolv2_(i,at1,at2) = d2;
				dsolv1_(i,at1,at2) = d1;
        minterp12->interpolate( dis(i), memb_solv1_(i,at1,at2), d1 );
        minterp22->interpolate( dis(i), memb_solv2_(i,at1,at2), d2 );
        memb_dsolv2_(i,at1,at2) = d2;
        memb_dsolv1_(i,at1,at2) = d1;
			}

		}
	}

}


////////////////////////////////////////////////////////////////////////////////
/// @begin output_etable
///
/// @brief output an etable data file in the same format used in input_etable
///
/// @detailed
///$$$ file first line is <etable> <etable_disbins>
///$$$ other lines are <atom type 1> <atomtype 1> <eval bin 1> <eval bin 2>...
///
/// @global_read
/// pdbstatistics_pack.h:  ljatr, dljatr, ljrep, dljrep, solv1,solv2,dsolv
///
///
/// @remarks
///
/// @references
///
/// @authors sheffler mar 19 2006
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::output_etable(
	ObjexxFCL::FArray3D<Real> & etable,
	std::string label,
	std::ostream & out
) {
	using namespace std;
	using namespace ObjexxFCL;

	out << label << " " << etable_disbins << endl;
	for(int at1 = 1; at1 <= n_atomtypes; at1++) {
		for(int at2 = 1; at2 <= n_atomtypes; at2++) {
			out << at1 << " "
					<< at2 << " ";
			for(int bin = 1; bin <= etable_disbins; bin++) {
				float evalue = etable(bin,at1,at2);
				out << evalue << ' ';
			}
			out << endl;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
/// @begin input_etable
///
/// @brief read in etable from a datafile
///
/// @detailed
///$$$ file first line is <etable> <etable_disbins>
///$$$ other lines are <atom type 1> <atomtype 1> <eval bin 1> <eval bin 2>...
///
/// @global_read
/// pdbstatistics_pack.h:  ljatr, dljatr, ljrep, dljrep, solv1,solv2,dsolv
///
/// @global_write
/// pdbstatistics_pack.h:  ljatr, dljatr, ljrep, dljrep, solv1,solv2,dsolv
///
/// @remarks
///
/// @references
///
/// @authors sheffler mar 19 2006
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::input_etable(
	ObjexxFCL::FArray3D<Real> & etable,
	const std::string label,
	std::istream & in
) {
	using namespace std;

	TR << "input_etable: reading etable... " << label << endl;
	istringstream intmp;
	string lblin;
	int numdisbin;
	char buf[100000];

	// check header
	in.getline(buf,100000);
	intmp.str(buf);
	if( !(intmp >> lblin >> numdisbin) ) {
		TR << "input_etable: WARNING bad etable header " << buf << endl;
		return;
	}
	if( lblin != label || etable_disbins != numdisbin ) {
		TR << "input_etable: WARNING etable types don't match! "<< endl;
		TR << "              expected " << label << "," << etable_disbins
			 << " got " << lblin << ',' << numdisbin<< endl;
		utility_exit_with_message( "input_etable: WARNING etable types don't match! " );
	} else {
		TR << "input_etable expected etable " << label << " of size " << etable_disbins
				 << ", got " << lblin << ',' << numdisbin<< endl;

	}

	// read in etable
	int at1,at2,count=0,scount=0;
	float evalue;
	while( in.getline(buf,100000) ) {
		//TR << "ASDFF buf  " << buf << endl;
		count++;
		intmp.clear();
		intmp.str(buf);
		if( ! (intmp >> at1 >> at2 ) ) {
			TR << "input_etable: error reading etable line: " << buf << endl;
		} else {
			for(int bin = 1; bin <=numdisbin; bin++) {
				if( !(intmp >> evalue) ) {
					TR << "input_etable: not enough bins on etable line: " << buf << endl;
					utility_exit_with_message( "input_etable: not enough bins on etable line: " );
				}
				//TR << "ASDFF read " << lblin << ' ' << at1 << ' ' << at2 << ' ' << bin << ' ' << evalue << endl;
				//TR << "ASDFF " << lblin << at1 << at2 << evalue;
				etable(bin,at1,at2) = evalue;
				etable(bin,at2,at1) = evalue;
			}
			scount++;
		}
	}
	TR << "              read " << scount << " of " << count << " lines" << endl;
}




////////////////////////////////////////////////////////////////////////////////
/// @begin precalc_etable_coefficients
///
/// @brief precalculate non-distance dependent coefficients of energy functions
///
/// @detailed
///
/// @param[out]   lj_sigma - out - for atomtypes i and j: (radius_i+radius_j)
/// @param[out]   lj_r6_coeff - out - precalced coefficient on the (1/dis)**6 term in lj
/// @param[out]   lj_r12_coeff - out - precalced coefficient on the (1/dis)**12 term in lj
/// @param[out]   lj_switch_intercept - out -
///            for close contacts calculate lj from a line with this intercept
/// @param[out]   lj_switch_slope - out -
///            for close contacts calculate lj from a line with this slope
/// @param[out]   lk_inv_lambda2 - out -
///            surprise! it's the 1/(lambda)**2 term in the lk equation
/// @param[out]   lk_coeff - out - precalculation of all non-distance dependent terms
///            outside of the exponential in the lk equation
/// @param[out]   lk_min_dis2sigma_value - out - below the min dis2sigma ratio for lk,
///            this value is assigned to the solvation
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors ctsa 10-2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::precalc_etable_coefficients(
	FArray2< Real > & lj_sigma,
	/*FArray2< Real > & lj_r6_coeff,
	FArray2< Real > & lj_r12_coeff,
	FArray2< Real > & lj_switch_intercept,
	FArray2< Real > & lj_switch_slope,*/
	FArray1< Real > & lk_inv_lambda2,
	FArray2< Real > & lk_coeff,
  FArray2< Real > & memb_lk_coeff, //pba
	FArray2< Real > & lk_min_dis2sigma_value,
  FArray2< Real > & memb_lk_min_dis2sigma_value //pba
)
{

// 	using namespace etable;
// 	using namespace param;
// 	using namespace water;
//   using namespace hbonds;

	// locals
	Real sigma; //sigma6,sigma12,wdepth;
	Real inv_lambda;
	FArray1D< Real > lk_coeff_tmp( n_atomtypes );
  FArray1D< Real > memb_lk_coeff_tmp( n_atomtypes ); //pba
	Real thresh_dis,inv_thresh_dis2,x_thresh;
	//int dstype;
	//int const atype_sulfur = { 16 };
	Real const inv_neg2_tms_pi_sqrt_pi = { -0.089793561062583294 };
	// coefficient for lk solvation

	// include follows locals so that data statements can initialize included arrays
	for ( int i = 1, e = n_atomtypes; i <= e; ++i ) {
		inv_lambda = 1.0/lk_lambda(i);
		//inv_lambda = 1.0/atom_type(i).lk_lambda();
		lk_inv_lambda2(i) = inv_lambda * inv_lambda;
		lk_coeff_tmp(i) = inv_neg2_tms_pi_sqrt_pi * lk_dgfree(i) * inv_lambda;
    memb_lk_coeff_tmp(i) = inv_neg2_tms_pi_sqrt_pi * memb_lk_dgfree(i) * inv_lambda; //pba
	}

	for ( int i = 1, e = n_atomtypes; i <= e; ++i ) {
		for ( int j = i; j <= e; ++j ) {

			sigma = Wradius * ( lj_radius(i) + lj_radius(j) );
			//sigma = Wradius * ( atom_type(i).lj_radius() + atom_type(j).lj_radius() );
//jjh temporary fix to prevent division by zero below
			sigma = ( sigma < 1.0e-9 ? 1.0e-9 : sigma );

			// (bk) modify sigma for hbond donors and acceptors
			//   ctsa -- should these scale down by Wradius as well?

			// pb specific sigma correction for pairs between charged oxygen acceptors (15)
			// pb and hydroxyl oxygen donors (13). sigma correction for all polar H and charged oxygen
			// pb acceptor. Combinations of these corrections allow better prediction of both
			// pb hydroxyl O donor/charged O acceptor and charged NH donor/charged O acceptor
			// pb distances.

			if ( lj_use_hbond_radii ) {
				if ( ( atom_type(i).is_acceptor() && atom_type(j).is_donor() ) ||
						 ( atom_type(i).is_donor() && atom_type(j).is_acceptor() ) ) {
					sigma = lj_hbond_dis;
//           if (tight_hb && ((i == 15 && j == 13) || (j == 15  && i == 13)))
// 						sigma = lj_hbond_accOch_dis;
				} else if ( ( atom_type(i).is_acceptor() && atom_type(j).is_polar_hydrogen() ) ||
										( atom_type(i).is_polar_hydrogen() && atom_type(j).is_acceptor() ) ) {
					sigma = lj_hbond_hdis;
//           if (tight_hb && ((i == 15 && j == 22) || (j == 15 && i == 22)))
// 						sigma = lj_hbond_accOch_hdis;
				}
			}

//lin   modify sigma for water and hbond donors/acceptors
			if ( lj_use_water_radii ) {
				if ( ( ( atom_type(i).is_acceptor() ||
								 atom_type(i).is_donor() ) &&
							 atom_type(j).is_h2o() ) ||
						 ( ( atom_type(j).is_acceptor() ||
								 atom_type(j).is_donor() ) &&
							 atom_type(i).is_h2o() ) ) {
					sigma = lj_water_dis;
				} else if ( ( atom_type(i).is_polar_hydrogen() &&
											atom_type(j).is_h2o() ) ||
										( atom_type(j).is_polar_hydrogen() &&
											atom_type(i).is_h2o() ) ) {
					sigma = lj_water_hdis;
				}
			}
/*
			sigma6  = std::pow( sigma, 6 );
			sigma12 = sigma6 * sigma6;
			wdepth = std::sqrt(lj_wdepth(i)*lj_wdepth(j));
			//wdepth = std::sqrt(atom_type(i).lj_wdepth()*atom_type(j).lj_wdepth());
*/
			lj_sigma(i,j) = sigma;
			lj_sigma(j,i) = lj_sigma(i,j);
/*
			lj_r6_coeff(i,j) = -2. * wdepth * sigma6;
			lj_r6_coeff(j,i) = lj_r6_coeff(i,j);

			lj_r12_coeff(i,j) = wdepth * sigma12;
			lj_r12_coeff(j,i) = lj_r12_coeff(i,j);

			// ctsa - create coefficients for linear projection of lj repulsive used
			//  for low distance values
			if ( lj_use_lj_deriv_slope ) {

				// ctsa - use the slope of the true lj repulsive at the
				//  linear switch point to create a linear projection of
				//  lj for low distances

				//  slope = wdepth/sigma *
				//          (slope@switch_point*sigma/wdepth)
				lj_switch_slope(i,j) = (wdepth/sigma)*
				 lj_switch_slope_sigma2wdepth;
				lj_switch_slope(j,i) = lj_switch_slope(i,j);

				// intercept = wdepth*(lj@switch_point/wdepth)
				//             - slope*switch_point_distance
				lj_switch_intercept(i,j) = wdepth*lj_switch_value2wdepth -
				 lj_switch_slope(i,j)*sigma*lj_switch_dis2sigma;
				lj_switch_intercept(j,i) = lj_switch_intercept(i,j);
			} else {

				// ctsa - create a linear projection of lj for low distances which
				//  is defined by a constant y intercept and the true lj repulsive
				//  value at the linear switch point
				lj_switch_slope(i,j) = -(1./sigma)*lj_switch_sigma2dis*
				 (lj_slope_intercept-wdepth*lj_switch_value2wdepth);
				lj_switch_slope(j,i) = lj_switch_slope(i,j);

				lj_switch_intercept(i,j) = lj_slope_intercept;
				lj_switch_intercept(j,i) = lj_switch_intercept(i,j);
			}
*/
			// ctsa - precalculated lk solvation coefficients
			lk_coeff(i,j) = lk_coeff_tmp(i) * lk_volume(j);
			//lk_coeff(i,j) = lk_coeff_tmp(i) * atom_type(j).lk_volume();
			lk_coeff(j,i) = lk_coeff_tmp(j) * lk_volume(i);
			//lk_coeff(j,i) = lk_coeff_tmp(j) * atom_type(i).lk_volume();
      memb_lk_coeff(i,j) = memb_lk_coeff_tmp(i) * lk_volume(j); //pba
      memb_lk_coeff(j,i) = memb_lk_coeff_tmp(j) * lk_volume(i); //pba

			// ctsa - when dis/sigma drops below lk_min_dis2sigma,
			//   a constant lk solvation value equal to the value at the
			//   switchover point is used. That switchover-point value
			//   is calculated here and stored in lk_min_dis2sigma_value
			thresh_dis = lk_min_dis2sigma*sigma;
			inv_thresh_dis2 = 1./( thresh_dis * thresh_dis );
			Real dis_rad = thresh_dis - lj_radius(i);
			//Real dis_rad = thresh_dis - atom_type(i).lj_radius();
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(i);
			lk_min_dis2sigma_value(i,j) = std::exp(-x_thresh) * lk_coeff(i,j) *
			 inv_thresh_dis2;
      memb_lk_min_dis2sigma_value(i,j) = std::exp(-x_thresh) * memb_lk_coeff(i,j) *
       inv_thresh_dis2; //pba

			dis_rad = thresh_dis - lj_radius(j);
			//dis_rad = thresh_dis - atom_type(j).lj_radius();
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(j);
			lk_min_dis2sigma_value(j,i) = std::exp(-x_thresh) * lk_coeff(j,i) *
			 inv_thresh_dis2;
      memb_lk_min_dis2sigma_value(j,i) = std::exp(-x_thresh) * memb_lk_coeff(j,i) *
       inv_thresh_dis2; //pba
    //pbadebug WARNING
    //if(i==18 && j==19)
    //std::cout << "lk_min_ij lk_min_ji " << lk_min_dis2sigma_value(i,j) << " " << lk_min_dis2sigma_value(j,i) << std::endl;

		}
	}

	// ctsa - calculate disulfide coefficients
	// pb -- killed this

}


////////////////////////////////////////////////////////////////////////////////
/// @begin calc_etable_value
///
/// @brief calc all etable values given a distance and atom-type pair
///
/// @detailed
///
/// given a pair of atom types and the squared inter-atomic separation
/// distance (and a whole bunch of pre-computed coeffecients), this returns
/// the value of the lennard-jones and lk solvation potentials and
/// their derivatives w.r.t. the separation distance
///
///
/// @param[in]   dis2 - in - atomic separation distance squared
/// @param[in]   atype1 - in - chemical type of atom 1
/// @param[in]   atype2 - in - chemical type of atom 2
/// @param[out]   atrE - out - atractive lj energy
/// @param[out]   d_atrE - out - d(atrE)/d(dis)
/// @param[out]   repE - out - repulsive lj energy
/// @param[out]   d_repE - out - d(repE)/d(dis)
/// @param[out]   solvE1 - out - lk solvation energy
/// @param[out]   solvE2 - out - lk solvation energy
/// @param[out]   dsolvE - out - d(solvE1+solvE2)/d(dis)
/// @param[in]   lj_sigma - in - for atomtypes i and j: (radius_i+radius_j)
/// @param[in]   lj_r6_coeff - in - precalced coefficient on the (1/dis)**6 term in lj
/// @param[in]   lj_r12_coeff - in - precalced coefficient on the (1/dis)**12 term in lj
/// @param[in]   lj_switch_intercept - in -
///            for close contacts calculate lj from a line with this intercept
/// @param[in]   lj_switch_slope - in -
///            for close contacts calculate lj from a line with this slope
/// @param[in]   lk_inv_lambda2 - in -
///            surprise! it's the 1/(lambda)**2 term in the lk equation
/// @param[in]   lk_coeff - in - precalculation of all non-distance dependent terms
///            outside of the exponential in the lk equation
/// @param[in]   lk_min_dis2sigma_value - in - below the min dis2sigma ratio for lk,
///            this value is assigned to the solvation
///
/// @global_read
///
/// pdbstatistics_pack.h
/// etable.h
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors ctsa 10-2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::calc_etable_value(
	Real & dis2,
	int & atype1,
	int & atype2,
	/*Real & atrE,
	Real & d_atrE,
	Real & repE,
	Real & d_repE,*/
	Real & solvE1,
	Real & solvE2,
	Real & dsolvE1,
	Real & dsolvE2,
	FArray2< Real > & lj_sigma,
	/*FArray2< Real > & lj_r6_coeff,
	FArray2< Real > & lj_r12_coeff,
	FArray2< Real > & lj_switch_intercept,
	FArray2< Real > & lj_switch_slope,*/
	FArray1< Real > & lk_inv_lambda2,
	FArray2< Real > & lk_coeff,
	FArray2< Real > & lk_min_dis2sigma_value,
  Real & memb_solvE1, //pba
  Real & memb_solvE2,
  FArray2< Real > & memb_lk_coeff,
  FArray2< Real > & memb_lk_min_dis2sigma_value,
  Real & memb_dsolvE1,
  Real & memb_dsolvE2
)
{
// 	using namespace etable;
// 	using namespace param;
// 	using namespace pdbstatistics_pack;

	// locals
	//Real ljE,d_ljE,
  Real x1,x2;
	Real dis;
	Real inv_dis,inv_dis2; //,inv_dis6,inv_dis7,inv_dis12,inv_dis13;
	Real dis2sigma;
	int xtra_atype1,xtra_atype2;
	//int const atype_sulfur = { 16 };

	// include after local variables to allow data statements to initialize
	/*atrE = 0.;
	d_atrE = 0.;
	repE = 0.;
	d_repE = 0.;*/
	solvE1 = 0.;
	solvE2 = 0.;
	dsolvE1 = 0.;
	dsolvE2 = 0.;
  memb_solvE1 = 0.; //pba
  memb_solvE2 = 0.;
  memb_dsolvE1 = 0.;
  memb_dsolvE2 = 0.;

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2 + epsilon ) return;

	if ( dis2 < min_dis2 ) dis2 = min_dis2;

	dis = std::sqrt(dis2);
	inv_dis = 1.0/dis;
	inv_dis2 = inv_dis * inv_dis;



	//  ctsa - switch to disulfide bonded atom types
	//    when conditions are met
// 	if ( ( atype1 == atype_sulfur && atype2 == atype_sulfur ) &&
// 	 dis < disulfide_dis_thresh ) {
// 		xtra_atype1 = n_atomtypes + 1;
// 		xtra_atype2 = n_atomtypes + 1;
// 	} else {
		xtra_atype1 = atype1;
		xtra_atype2 = atype2;
		//}


	dis2sigma = dis / lj_sigma(xtra_atype1,xtra_atype2);

/*
	if ( dis2sigma < lj_switch_dis2sigma ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma
		//    ratio drops below theshold
		d_ljE = lj_switch_slope(xtra_atype1,xtra_atype2);
		ljE = dis*d_ljE + lj_switch_intercept(xtra_atype1,xtra_atype2);
	} else {
		//  ctsa - calc regular lennard-jones
		inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
		inv_dis7  = inv_dis6 * inv_dis;
		inv_dis12 = inv_dis6 * inv_dis6;
		inv_dis13 = inv_dis12 * inv_dis;

		ljE = lj_r12_coeff(xtra_atype1,xtra_atype2) * inv_dis12 +
		 lj_r6_coeff(xtra_atype1,xtra_atype2) * inv_dis6;

		d_ljE = -12.*lj_r12_coeff(xtra_atype1,xtra_atype2) * inv_dis13-6. *
		 lj_r6_coeff(xtra_atype1,xtra_atype2) * inv_dis7;
	}

	if ( ljE < 0. ) {
		atrE = ljE;
		d_atrE = d_ljE;
	} else {
		repE = ljE;
		d_repE = d_ljE;
	}
*/
    //pbadebug WARNING
    //if(atype1==18 && atype2==19)
    //  std::cout << "d lj d2s lkmin " << dis << " " << lj_sigma(xtra_atype1,xtra_atype2) << " " << dis2sigma << " " << lk_min_dis2sigma << std::endl;
	// ctsa - calc lk
	if ( dis2sigma < lk_min_dis2sigma ) {
		// ctsa - solvation is constant when the dis/sigma ratio
		//   falls below minimum threshold
		solvE1 = lk_min_dis2sigma_value(xtra_atype1,xtra_atype2);
		solvE2 = lk_min_dis2sigma_value(xtra_atype2,xtra_atype1);
		dsolvE1 = 0.0;
    dsolvE2 = 0.0;
    memb_solvE1 = memb_lk_min_dis2sigma_value(xtra_atype1,xtra_atype2); //pba
    memb_solvE2 = memb_lk_min_dis2sigma_value(xtra_atype2,xtra_atype1);
    memb_dsolvE1 = memb_dsolvE2 = 0.0;

	} else {

		Real dis_rad = dis - lj_radius(atype1);
		//Real dis_rad = dis - atom_type(atype1).lj_radius();
		x1 = ( dis_rad * dis_rad ) * lk_inv_lambda2(atype1);
		dis_rad = dis - lj_radius(atype2);
		//dis_rad = dis - atom_type(atype2).lj_radius();
		x2 = ( dis_rad * dis_rad ) * lk_inv_lambda2(atype2);

		solvE1 = std::exp(-x1) * lk_coeff(atype1,atype2) * inv_dis2;
		solvE2 = std::exp(-x2) * lk_coeff(atype2,atype1) * inv_dis2;
    memb_solvE1 = std::exp(-x1) * memb_lk_coeff(atype1,atype2) * inv_dis2; //pba
    memb_solvE2 = std::exp(-x2) * memb_lk_coeff(atype2,atype1) * inv_dis2;

    //pbadebug WARNING
    //if(atype1==18 && atype2==19)
    //  std::cout << "d x1 x2 lk1 lk2 s1 s2 " << dis << " " << x1 << " " << x2 << " " << lk_coeff(atype1,atype2) << " " << lk_coeff(atype2,atype1) << " " << solvE1 << " " << solvE2 << std::endl;

		// ctsa - get d(lk_E)/dr
		dsolvE1 = -2.0 * solvE1 *
		 (((dis-lj_radius(atype1))*lk_inv_lambda2(atype1))+inv_dis);
		//(((dis-atom_type(atype1).lj_radius())*lk_inv_lambda2(atype1))+inv_dis);
		dsolvE2 = -2.0 * solvE2 *
		 (((dis-lj_radius(atype2))*lk_inv_lambda2(atype2))+inv_dis);
		//(((dis-atom_type(atype2).lj_radius())*lk_inv_lambda2(atype2))+inv_dis);

    // pba - get d(memblk_E)/dr
    memb_dsolvE1 = -2.0 * memb_solvE1 *
     (((dis-lj_radius(atype1))*lk_inv_lambda2(atype1))+inv_dis);
    memb_dsolvE2 = -2.0 * memb_solvE2 *
     (((dis-lj_radius(atype2))*lk_inv_lambda2(atype2))+inv_dis);

	}

}

/*
void
Etable::zero_hydrogen_and_water_ljatr()
{
	int const HOH = atom_set_->atom_type_index("HOH");
	for( int at1 = 1; at1 <= n_atomtypes; ++at1 ) {
		for( int at2 = 1; at2 <= n_atomtypes; ++at2 ) {

			// cbk  don't give hydrogens or water attractive lennard-jones
			// cbk  this is so very short range cut-offs can be used
			//if ( ( at1 >= 22 && at2 <= 26 ) || ( at2 >= 22 && at2 <= 26 ) ) { // BUG!  introduced in r19802
			if( atom_type(at1).is_hydrogen() || atom_type(at2).is_hydrogen() || at1 == HOH || at2 == HOH ) {
				for( int i = 1; i <= etable_disbins; ++i ) {
		 			 ljatr_(i,at1,at2) = 0.0;
					dljatr_(i,at1,at2) = 0.0;
				}
			}

		}
	}

}
*/
/// @brief Returns the maximum lj radius for any non-hydrogen
/// atom as defined by the atom-type-set used to create this Etable.
Real
MembEtable::max_non_hydrogen_lj_radius() const
{
	return max_non_hydrogen_lj_radius_;
}

/// @brief Returns the maximum lj radius for any hydrogen atom as
/// defined by the input atom-type-set used to create this Etable.
Real
MembEtable::max_hydrogen_lj_radius() const
{
	return max_hydrogen_lj_radius_;
}



} // etable
} // scoring
} // core
