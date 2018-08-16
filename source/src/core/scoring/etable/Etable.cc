// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for generating the table for fa_atr/rep and fa_sol
///
/// @details
/// This class is called upon by the ScoringManager. Since actual calculating of the LJ potential
/// is time consuming if done multiple times, this class precomputes and discritizes the potential
/// (meaning that the potential is broken down into bins). Once the bins have been created, it will
/// smooth out the bins, for better interpolation.
///
///
/// @author
/// I dont know?
/// Steven Combs - comments and skipping of virtual atoms
///
/////////////////////////////////////////////////////////////////////////


// Unit Headers
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/EnergyMap.hh>

// Package headers
#include <core/chemical/util.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>

#include <basic/options/option.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/cubic_polynomial.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility Headers

// C++ Headers
#include <iostream>
#include <fstream>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <core/chemical/AtomType.hh>

using basic::Error;
using basic::Warning;
using namespace ObjexxFCL;

static basic::Tracer TR( "core.scoring.etable" );

namespace core {
namespace scoring {
namespace etable {

/// @details Auto-generated virtual destructor
Etable::~Etable() = default;

using namespace basic::options;
using namespace basic::options::OptionKeys;

///  Note from rhiju & andrew leaver-day: we should get rid of command-line option calls here
///  in favor of setting options from EtableOptions object. Would allow user to have multiple
///  etables around!

//  constructor
Etable::Etable(
	chemical::AtomTypeSetCAP atom_set_in, // like etable namespace
	EtableOptions const & options,
	std::string const & alternate_parameter_set // = ""
) :
	// from atop_props_in:
	atom_set_                 ( atom_set_in ),
	n_atomtypes_              ( atom_set_in.lock()->n_atomtypes() ),

	// from options
	max_dis_                  ( options.max_dis ),
	bins_per_A2_              ( options.bins_per_A2 ),
	Wradius_                  ( options.Wradius ), // global mod to radii
	lj_switch_dis2sigma_       ( options.lj_switch_dis2sigma ),
	max_dis2_                  ( max_dis_*max_dis_ ),
	etable_disbins_            ( static_cast< int >( max_dis2_ * bins_per_A2_)+1),
	lj_hbond_OH_donor_dis_     ( options.lj_hbond_OH_donor_dis ),
	lj_hbond_dis_              ( 3.0 ),

	// hard-coded for now
	lj_use_lj_deriv_slope_     ( true ),
	lj_slope_intercept_        ( 0.0 ),
	lj_use_hbond_radii_        ( true ),
	lj_hbond_hdis_             ( options.lj_hbond_hdis ),
	lj_use_water_radii_        ( true ),
	lj_water_dis_              ( 3.0 ),
	lj_water_hdis_             ( 1.95 ),
	enlarge_h_lj_wdepth_      ( options.enlarge_h_lj_wdepth ),
	lk_min_dis2sigma_          ( 0.89 ),
	no_lk_polar_desolvation_  ( options.no_lk_polar_desolvation ),
	proline_N_is_lk_nonpolar_ ( options.proline_N_is_lk_nonpolar ),
	min_dis_                   ( 0.01 ),
	min_dis2_                  ( min_dis_ * min_dis_ ),
	add_long_range_damping_    ( true ),
	long_range_damping_length_ ( 0.5 ),
	epsilon_                   ( 0.0001 ),
	safe_max_dis2_             ( max_dis2_ - epsilon_ ),
	hydrogen_interaction_cutoff2_( options.fa_hatr ?
	std::pow( max_dis_ + 2*chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH, 2 ) :
	std::pow(5.0,2) ),
	max_non_hydrogen_lj_radius_( 0.0 ),
	max_hydrogen_lj_radius_( 0.0 ),
	fa_hatr_( options.fa_hatr ),
	slim_( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] ) // command-line option -- bad form; see note above.
{
	dimension_etable_arrays();
	initialize_from_input_atomset( atom_set_in );
	calculate_nblist_distance_thresholds( options );
	read_alternate_parameter_set( atom_set_in, alternate_parameter_set );
	calculate_hydrogen_atom_reach();
	initialize_carbontypes_to_linearize_fasol();
	make_pairenergy_table();
}

void
Etable::dimension_etable_arrays()
{
	// size the arrays
	if ( !slim_ ) {
		ljatr_.dimension(  etable_disbins_, n_atomtypes_, n_atomtypes_ );
		ljrep_.dimension(  etable_disbins_, n_atomtypes_, n_atomtypes_ );
		solv1_.dimension(  etable_disbins_, n_atomtypes_, n_atomtypes_ );
		solv2_.dimension(  etable_disbins_, n_atomtypes_, n_atomtypes_ );
		dljatr_.dimension( etable_disbins_, n_atomtypes_, n_atomtypes_ );
		dljrep_.dimension( etable_disbins_, n_atomtypes_, n_atomtypes_ );
		dsolv_.dimension(  etable_disbins_, n_atomtypes_, n_atomtypes_ );
		dsolv1_.dimension( etable_disbins_, n_atomtypes_, n_atomtypes_ );
	}

	lj_radius_.resize( n_atomtypes_, 0.0 );
	lj_wdepth_.resize( n_atomtypes_, 0.0 );
	lk_dgfree_.resize( n_atomtypes_, 0.0 );
	lk_lambda_.resize( n_atomtypes_, 0.0 );
	lk_volume_.resize( n_atomtypes_, 0.0 );

	lj_sigma_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );
	lj_r6_coeff_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );
	lj_r12_coeff_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );
	lj_switch_intercept_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );
	lj_switch_slope_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );
	lk_inv_lambda2_.dimension( n_atomtypes_ );
	lk_coeff_.dimension( n_atomtypes_, n_atomtypes_ );
	lk_min_dis2sigma_value_.dimension( n_atomtypes_ + 1, n_atomtypes_ + 1 );

	Size n_unique_pair_types = n_atomtypes_ * n_atomtypes_ - ( n_atomtypes_ * (n_atomtypes_ - 1 ) / 2 );
	analytic_parameters_.resize( n_unique_pair_types );

}

void
Etable::initialize_from_input_atomset(
	chemical::AtomTypeSetCAP atom_set_in_ap
)
{
	chemical::AtomTypeSetCOP atom_set_in( atom_set_in_ap );
	for ( int i=1; i<= n_atomtypes_; ++i ) {
		lj_radius_[i] = (*atom_set_in)[i].lj_radius();
		lj_wdepth_[i] = (*atom_set_in)[i].lj_wdepth();
		lk_dgfree_[i] = (*atom_set_in)[i].lk_dgfree();
		lk_lambda_[i] = (*atom_set_in)[i].lk_lambda();
		lk_volume_[i] = (*atom_set_in)[i].lk_volume();
		if ( (*atom_set_in)[i].is_hydrogen() ) {
			if ( lj_radius_[i] > max_hydrogen_lj_radius_ ) max_hydrogen_lj_radius_ = lj_radius_[i];
		} else {
			if ( lj_radius_[i] > max_non_hydrogen_lj_radius_ ) max_non_hydrogen_lj_radius_ = lj_radius_[i];
		}
	}

	// rhiju/fang -- Use larger LJ_WDEPTH for protons to avoid clashes in RNA
	if ( enlarge_h_lj_wdepth_ ) core::chemical::enlarge_h_lj_wdepth( lj_wdepth_, *atom_set_in );

	// APL -- hydrophobic desolvation only; turn off hydrophilic desolvation penalty
	if ( no_lk_polar_desolvation_ ) {
		chemical::AtomTypeSetCOP atom_set( atom_set_ );
		for ( int i=1; i<= n_atomtypes_; ++i ) {
			if ( (*atom_set)[i].is_acceptor() || (*atom_set)[i].is_donor() ) {
				if ( !proline_N_is_lk_nonpolar_ || (*atom_set)[i].atom_type_name() != "Npro" ) { // bazzoli
					lk_dgfree_[i] = 0.0;
				}
			}
		}
	}

}

void
Etable::calculate_nblist_distance_thresholds(
	EtableOptions const & options
)
{
	Real const max_tolerated_movement = 0.75;
	max_heavy_hydrogen_cutoff_    = fa_hatr_ ? max_dis_ : max_non_hydrogen_lj_radius_ + max_hydrogen_lj_radius_;
	max_hydrogen_hydrogen_cutoff_ = fa_hatr_ ? max_dis_ : 2 * max_hydrogen_lj_radius_;
	nblist_dis2_cutoff_XX_        = std::pow(options.max_dis + 2*max_tolerated_movement, 2 );
	nblist_dis2_cutoff_XH_        = std::pow(max_heavy_hydrogen_cutoff_ + 2*max_tolerated_movement, 2 );
	nblist_dis2_cutoff_HH_        = std::pow(max_hydrogen_hydrogen_cutoff_+ 2*max_tolerated_movement, 2 );
}

void
Etable::read_alternate_parameter_set(
	chemical::AtomTypeSetCAP atom_set_in_ap,
	std::string const & alternate_parameter_set
)
{
	if ( ! alternate_parameter_set.size() )  return;

	// uses alternate paramers
	std::string param_name;
	chemical::AtomTypeSetCOP atom_set_in( atom_set_in_ap );

	param_name = "LJ_RADIUS_"+alternate_parameter_set;
	if ( atom_set_in->has_extra_parameter( param_name ) ) {
		TR << "Using alternate parameters: " << param_name << " in Etable construction." << std::endl;
		Size const index( atom_set_in->extra_parameter_index( param_name ) );
		for ( int i=1; i<= n_atomtypes_; ++i ) lj_radius_[i] = (*atom_set_in)[i].extra_parameter( index );
	}

	param_name = "LJ_WDEPTH_"+alternate_parameter_set;
	if ( atom_set_in->has_extra_parameter( param_name ) ) {
		TR << "Using alternate parameters: " << param_name << " in Etable construction."<<std::endl;
		Size const index( atom_set_in->extra_parameter_index( param_name ) );
		for ( int i=1; i<= n_atomtypes_; ++i ) lj_wdepth_[i] = (*atom_set_in)[i].extra_parameter( index );
	}

	param_name = "LK_DGFREE_"+alternate_parameter_set;
	if ( atom_set_in->has_extra_parameter( param_name ) ) {
		TR << "Using alternate parameters: " << param_name << " in Etable construction."<<std::endl;
		Size const index( atom_set_in->extra_parameter_index( param_name ) );
		for ( int i=1; i<= n_atomtypes_; ++i ) lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( index );
	}

	param_name = "LK_LAMBDA_"+alternate_parameter_set;
	if ( atom_set_in->has_extra_parameter( param_name ) ) {
		TR << "Using alternate parameters: " << param_name << " in Etable construction." << std::endl;
		Size const index( atom_set_in->extra_parameter_index( param_name ) );
		for ( int i=1; i<= n_atomtypes_; ++i ) lk_lambda_[i] = (*atom_set_in)[i].extra_parameter( index );
	}

	param_name = "LK_VOLUME_"+alternate_parameter_set;
	if ( atom_set_in->has_extra_parameter( param_name ) ) {
		TR << "Using alternate parameters: " << param_name << " in Etable construction." << std::endl;
		Size const index( atom_set_in->extra_parameter_index( param_name ) );
		for ( int i=1; i<= n_atomtypes_; ++i ) lk_volume_[i] = (*atom_set_in)[i].extra_parameter( index );
	}
}

void
Etable::calculate_hydrogen_atom_reach()
{
	Real const MAX_H_HEAVY_DISTANCE = chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH;

	Real max_lj_rep_for_h = std::max(
		2 * max_hydrogen_lj_radius_ + 2 * MAX_H_HEAVY_DISTANCE,
		max_hydrogen_lj_radius_ + MAX_H_HEAVY_DISTANCE + max_non_hydrogen_lj_radius_ );

	hydrogen_interaction_cutoff2_ = fa_hatr_ ?
		std::pow( 2 * MAX_H_HEAVY_DISTANCE + max_dis_, 2 ) : max_lj_rep_for_h * max_lj_rep_for_h;
}

void
Etable::initialize_carbontypes_to_linearize_fasol()
{
	chemical::AtomTypeSetCOP atom_set( atom_set_);
	carbon_types_.push_back( atom_set->atom_type_index("CH1") );
	carbon_types_.push_back( atom_set->atom_type_index("CH2") );
	carbon_types_.push_back( atom_set->atom_type_index("CH3") );
	carbon_types_.push_back( atom_set->atom_type_index("aroC") );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculate fast lookup arrays for vdw and solvation energy
///
/// @details
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
/// @author ctsa 10-2003
///
/////////////////////////////////////////////////////////////////////////////////
void
Etable::make_pairenergy_table()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//  locals

	ObjexxFCL::FArray1D< Real > ljatr( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > dljatr( etable_disbins_);
	ObjexxFCL::FArray1D< Real > ljrep( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > dljrep( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > fasol1( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > fasol2( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > dfasol( etable_disbins_ );
	ObjexxFCL::FArray1D< Real > dfasol1( etable_disbins_ );

	// The index defining the range [1..normal_disbins] for which the
	// energy function is calculated analytically.
	int const normal_disbins = calculate_normal_disbins();

	//  etable parameters
	TR << "Starting energy table calculation" << std::endl;


	//  std::TR << "Energy table parameter set: " << etable_label << ' ' <<
	//   etable_revision.substr( 1, etable_revision.length() - 3 ) << ' ' <<
	//   etable_date.substr( 1, etable_date.length() - 3 ) << std::endl;

	// ctsa - values precomputed from etable parameters:
	//   these were all originally constants but need
	//   to be calculated here b/c the switch value:
	//   lj_switch_dis2sigma
	//   is set at runtime
	lj_switch_sigma2dis_ = 1.0/lj_switch_dis2sigma_;

	// ctsa - value of the lennard-jones potential at the linear
	//   switch point divided by the wdepth (note that this
	//   coefficient is independent of atomtype)
	lj_switch_value2wdepth_ = std::pow( lj_switch_sigma2dis_, 12 ) -
		2.0 * std::pow( lj_switch_sigma2dis_, 6 );

	// ctsa - slope of the lennard-jones potential at the linear
	//   switch point times sigma divided by wdepth (note that this
	//   coefficient is independent of atomtype)
	lj_switch_slope_sigma2wdepth_ = -12.0 * (
		std::pow( lj_switch_sigma2dis_, 13 ) -
		std::pow( lj_switch_sigma2dis_, 7 ) );

	//  initialize non-distance dependent coefficients
	precalc_etable_coefficients( lj_sigma_, lj_r6_coeff_, lj_r12_coeff_,
		lj_switch_intercept_, lj_switch_slope_, lk_inv_lambda2_, lk_coeff_,
		lk_min_dis2sigma_value_ );

	//  calc distance**2 step per bin
	Real const dis2_step = 1.0 / bins_per_A2_;


	bool const only_save_one_way = basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ];

	//  ctsa - step through distance**2 bins and calculate potential
	for ( int atype1 = 1, atype_end = n_atomtypes_; atype1 <= atype_end; ++atype1 ) {
		//check to make sure that atype1 is not virtual. This is important later
		bool atype1_virtual(atom_type(atype1).is_virtual());
		for ( int atype2 = only_save_one_way ? atype1 : 1 ; atype2 <= atype_end; ++atype2 ) {
			bool atype2_virtual(atom_type(atype2).is_virtual());

			//  ctsa - normal bins have their lj and lk values
			//  calculated analytically
			for ( int disbin = 1; disbin <= normal_disbins; ++disbin ) {
				Real const dis2 = ( disbin - 1 ) * dis2_step;
				if ( atype1_virtual || atype2_virtual ) {
					ljatr(  disbin) = 0;
					dljatr( disbin) = 0;
					ljrep(  disbin) = 0;
					dljrep( disbin) = 0;
					fasol1( disbin) = 0;
					fasol2( disbin) = 0;
					dfasol( disbin) = 0;
					dfasol1(disbin) = 0;
				} else {
					Real atrE,d_atrE,repE,d_repE,solvE1,solvE2,dsolvE1,dsolvE2;
					calc_etable_value(dis2,atype1,atype2,atrE,d_atrE,repE,d_repE,solvE1,
						solvE2,dsolvE1,dsolvE2 );

					ljatr(  disbin) = atrE;
					dljatr( disbin) = d_atrE;
					ljrep(  disbin) = repE;
					dljrep( disbin) = d_repE;
					fasol1( disbin) = solvE1;
					fasol2( disbin) = solvE2;
					dfasol( disbin) = dsolvE1 + dsolvE2;
					dfasol1(disbin) = dsolvE1;
				}
			}
			// save these parameters for the analytic evaluation of the etable energy
			if ( atype1 <= atype2 ) {
				EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
				p.maxd2 = safe_max_dis2_;
				p.ljrep_linear_ramp_d2_cutoff = std::pow( lj_switch_dis2sigma_ * lj_sigma_( atype1, atype2 ), 2); // ie dis / lj_sigma < lj_switch_d2sigma
				p.lj_r6_coeff = lj_r6_coeff_( atype1, atype2 );
				p.lj_r12_coeff = lj_r12_coeff_( atype1, atype2 );
				p.lj_switch_intercept = lj_switch_intercept_( atype1, atype2 );
				p.lj_switch_slope = lj_switch_slope_( atype1, atype2 );
				p.lk_coeff1 = lk_coeff_( atype1, atype2 );
				p.lk_coeff2 = lk_coeff_( atype2, atype1 );
				p.lk_min_dis2sigma_value = lk_min_dis2sigma_value_( atype2, atype1 );
			}

			if ( add_long_range_damping_ ) {
				damp_long_range( normal_disbins,
					ljatr, dljatr, ljrep, dljrep,
					fasol1, fasol2, dfasol, dfasol1 );
			}

			//  ctsa - set last bin of all values to zero
			ljatr( etable_disbins_) = 0.0;
			ljrep( etable_disbins_) = 0.0;
			fasol1(etable_disbins_) = 0.0;
			fasol2(etable_disbins_) = 0.0;

			// throwing this all into one function
			modify_pot_one_pair(
				atype1,atype2,
				ljatr, dljatr, ljrep, dljrep,
				fasol1, fasol2, dfasol, dfasol1 );

			if ( !option[ score::no_smooth_etables ] ) {
				smooth_etables_one_pair(
					atype1, atype2,
					ljatr, dljatr, ljrep, dljrep,
					fasol1, fasol2, dfasol, dfasol1 );
			}

			if ( ! fa_hatr_ ) {
				zero_hydrogen_and_water_ljatr_one_pair(
					atype1, atype2,
					ljrep, ljatr, dljatr,
					fasol1, fasol2, dfasol, dfasol1 );
			}

			if ( ! slim_ ) {
				assign_parameters_to_full_etables(
					atype1, atype2,
					ljatr, dljatr, ljrep, dljrep,
					fasol1, fasol2, dfasol, dfasol1 );
			}
		}
	}


	////////////////////////////////////////////////////////
	// etable I/O stuff
	////////////////////////////////////////////////////////
	using namespace std;
	using namespace ObjexxFCL;

	// just for convenience in input/output_etables
	map< string, FArray3D<Real>* > etables;
	etables[ "ljatr"] = &  ljatr_; etables[ "ljrep"] = &  ljrep_;
	etables["dljatr"] = & dljatr_; etables["dljrep"] = & dljrep_;
	etables[ "solv1"] = &  solv1_; etables[ "solv2"] = &  solv2_;
	etables["dsolv"]  = & dsolv_; etables["dsolv1"]  = & dsolv1_;


	if ( option[ score::input_etables ].user() ) {
		string tag = option[ score::input_etables ];
		TR << "INPUT ETABLES " << tag << std::endl;
		for ( auto & etable : etables ) {
			string ename = etable.first;
			string fname = tag+"."+ename+".etable";
			std::ifstream input( fname.c_str() ); // TODO sheffler: figure out how to do this the right way
			input_etable(*(etable.second),ename,input);
			input.close();
		}
	}

	if ( option[ score::output_etables ].user() ) {
		string header = option[ score::output_etables ];
		TR << "OUTPUT ETABLES " << header << std::endl;
		for ( auto & etable : etables ) {
			string ename = etable.first;
			string fname = header+"."+ename+".etable";
			TR << "output_etable: writing etable: " << ename << " to " << fname << std::endl;
			ofstream out(fname.c_str());
			output_etable(*(etable.second),ename,out);
			out.close();
		}
	}

	TR << "Finished calculating energy tables." << std::endl;
}

int
Etable::calculate_normal_disbins() const
{
	// parameters to calculate the damping at max_dis range
	Real damping_thresh_dis2;
	int normal_disbins;

	//  get the number of damping disbins
	if ( add_long_range_damping_ ) {
		Real const dif = max_dis_ - long_range_damping_length_;
		damping_thresh_dis2 = max_dis2_ - ( dif * dif );
		auto damping_disbins = static_cast< int >( damping_thresh_dis2*bins_per_A2_ );
		normal_disbins = etable_disbins_-damping_disbins;
	} else {
		normal_disbins = etable_disbins_;
	}
	return normal_disbins;
}

void
Etable::damp_long_range(
	int const normal_disbins,
	ObjexxFCL::FArray1A< Real > ljatr,
	ObjexxFCL::FArray1A< Real > dljatr,
	ObjexxFCL::FArray1A< Real > ljrep,
	ObjexxFCL::FArray1A< Real > dljrep,
	ObjexxFCL::FArray1A< Real > fasol1,
	ObjexxFCL::FArray1A< Real > fasol2,
	ObjexxFCL::FArray1A< Real > dfasol,
	ObjexxFCL::FArray1A< Real > dfasol1
)
{
	ljatr.dimension(   etable_disbins_ );
	dljatr.dimension(  etable_disbins_ );
	ljrep.dimension(   etable_disbins_ );
	dljrep.dimension(  etable_disbins_ );
	fasol1.dimension(  etable_disbins_ );
	fasol2.dimension(  etable_disbins_ );
	dfasol.dimension(  etable_disbins_ );
	dfasol1.dimension( etable_disbins_ );

	// ctsa - remaining bins damp to 0. on a linear path
	Real const dljatr_damp = -ljatr( normal_disbins) / long_range_damping_length_;
	Real const dljrep_damp = -ljrep( normal_disbins) / long_range_damping_length_;
	Real const dsolv1_damp = -fasol1(normal_disbins) / long_range_damping_length_;
	Real const dsolv2_damp = -fasol2(normal_disbins) / long_range_damping_length_;

	Real const intercept_ljatr_damp = -dljatr_damp*max_dis_;
	Real const intercept_ljrep_damp = -dljrep_damp*max_dis_;
	Real const intercept_solv1_damp = -dsolv1_damp*max_dis_;
	Real const intercept_solv2_damp = -dsolv2_damp*max_dis_;

	Real const dis2_step = 1.0 / bins_per_A2_;

	for ( int disbin = normal_disbins+1; disbin <= etable_disbins_; ++disbin ) {
		Real const dis2 = ( disbin - 1 ) * dis2_step;
		Real const dis = std::sqrt(dis2);

		ljatr( disbin) = intercept_ljatr_damp + dis *dljatr_damp;
		ljrep( disbin) = intercept_ljrep_damp + dis *dljrep_damp;
		fasol1(disbin) = intercept_solv1_damp + dis *dsolv1_damp;
		fasol2(disbin) = intercept_solv2_damp + dis *dsolv2_damp;

		dljatr( disbin) = dljatr_damp;
		dljrep( disbin) = dljrep_damp;
		dfasol( disbin) = dsolv1_damp + dsolv2_damp;
		dfasol1(disbin) = dsolv1_damp;
	}

}

void
Etable::assign_parameters_to_full_etables(
	Size atype1,
	Size atype2,
	ObjexxFCL::FArray1A< Real > ljatr,
	ObjexxFCL::FArray1A< Real > dljatr,
	ObjexxFCL::FArray1A< Real > ljrep,
	ObjexxFCL::FArray1A< Real > dljrep,
	ObjexxFCL::FArray1A< Real > fasol1,
	ObjexxFCL::FArray1A< Real > fasol2,
	ObjexxFCL::FArray1A< Real > dfasol,
	ObjexxFCL::FArray1A< Real > dfasol1
)
{
	debug_assert( ljatr_.size() != 0 );

	ljatr.dimension(   etable_disbins_ );
	dljatr.dimension(  etable_disbins_ );
	ljrep.dimension(   etable_disbins_ );
	dljrep.dimension(  etable_disbins_ );
	fasol1.dimension(  etable_disbins_ );
	fasol2.dimension(  etable_disbins_ );
	dfasol.dimension(  etable_disbins_ );
	dfasol1.dimension( etable_disbins_ );

	// Take slices of the large arrays
	ObjexxFCL::FArray1A< Real > ljatr_full(   ljatr_(  1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > dljatr_full(  dljatr_( 1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > ljrep_full(   ljrep_(  1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > dljrep_full(  dljrep_( 1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > fasol1_full(  solv1_(  1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > fasol2_full(  solv2_(  1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > dfasol_full(  dsolv_(  1, atype2, atype1 ) );
	ObjexxFCL::FArray1A< Real > dfasol1_full( dsolv1_( 1, atype2, atype1 ) );

	ljatr_full.dimension(   etable_disbins_ );
	dljatr_full.dimension(  etable_disbins_ );
	ljrep_full.dimension(   etable_disbins_ );
	dljrep_full.dimension(  etable_disbins_ );
	fasol1_full.dimension(  etable_disbins_ );
	fasol2_full.dimension(  etable_disbins_ );
	dfasol_full.dimension(  etable_disbins_ );
	dfasol1_full.dimension( etable_disbins_ );

	// Slice assignment
	ljatr_full   = ljatr;
	dljatr_full  = dljatr;
	ljrep_full   = ljrep;
	dljrep_full  = dljrep;
	fasol1_full  = fasol1;
	fasol2_full  = fasol2;
	dfasol_full  = dfasol;
	dfasol1_full = dfasol1;


}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief modify Etable to better treat 0-0, C-C, and H-H interactions
///
/// @details
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
/// Comments from DB:
///db  the following function call modifies the potential in three ways:
///db     (1) the solvation energy for nonpolar atoms is held constant below
///db     4.2A to avoid shifting the minimum in the LJ potential.
///db     (2) a short range repulsion is added between backbone oxygens which are
///db     otherwise brought too close together by the LJatr.
///db     (3) the range of the repulsive interaction between non polar hydrogens is
///db     increased slightly.  (this is currently commented out because the effects
///db     on design have not been tested)
//db  all three modifications are based on inspection of the atom pair distributions
//db  after extensive refinement.


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
/// @author ctsa 10-2003
///
void
Etable::modify_pot_one_pair(
	Size const atype1,
	Size const atype2,
	ObjexxFCL::FArray1A< Real > ljatr,
	ObjexxFCL::FArray1A< Real > dljatr,
	ObjexxFCL::FArray1A< Real > ljrep,
	ObjexxFCL::FArray1A< Real > dljrep,
	ObjexxFCL::FArray1A< Real > fasol1,
	ObjexxFCL::FArray1A< Real > fasol2,
	ObjexxFCL::FArray1A< Real > dfasol,
	ObjexxFCL::FArray1A< Real > dfasol1
)
{

	ljatr.dimension(   etable_disbins_ );
	dljatr.dimension(  etable_disbins_ );
	ljrep.dimension(   etable_disbins_ );
	dljrep.dimension(  etable_disbins_ );
	fasol1.dimension(  etable_disbins_ );
	fasol2.dimension(  etable_disbins_ );
	dfasol.dimension(  etable_disbins_ );
	dfasol1.dimension( etable_disbins_ );

	using namespace std;

	bool mod_hhrep      = false; //truefalseoption("mod_hhrep");
	Real const h = 0.0;// fullatom_setup_ns::mod_hhrep_height;
	Real const c = 0.0;// fullatom_setup_ns::mod_hhrep_center;
	Real const w = 0.0;// fullatom_setup_ns::mod_hhrep_width;
	Real const e = 0.0;// fullatom_setup_ns::mod_hhrep_exponent;

	bool const skip_mod_OCbb_OCbb_rep = basic::options::option[ basic::options::OptionKeys::score::unmodifypot ];

	if ( mod_hhrep ) {
		TR << "fullatom_setup: modifying h-h repulsion "
			<< " hhrep center "   << c
			<< " hhrep height "   << h
			<< " hhrep width "    << w
			<< " hhrep exponent " << e
			<< std::endl;
	}

	{
		chemical::AtomTypeSetCOP atom_set( atom_set_ );
		Size const OCbb_idx = atom_set->atom_type_index("OCbb");
		if ( atype1 == OCbb_idx && atype2 == OCbb_idx && ! skip_mod_OCbb_OCbb_rep ) {
			ExtraQuadraticRepulsion OCbb_OCbb_exrep;
			OCbb_OCbb_exrep.xlo = 2.1; OCbb_OCbb_exrep.xhi = 3.6; OCbb_OCbb_exrep.slope = 2;
			//OCbb_OCbb_exrep.extrapolated_slope = 0; OCbb_OCbb_exrep.ylo = 4.5; // 2*(1.5)^2;
			OCbb_OCbb_exrep.extrapolated_slope = -6; /* d/dx 2*( 3.6 - x )^2 at x=2.1 ==> -2*2*1.5 */ OCbb_OCbb_exrep.ylo = 4.5; // 2*(1.5)^2;
			//ljrep_extra_repulsion( OCbb_idx, OCbb_idx ) = OCbb_OCbb_exrep;
			EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
			p.ljrep_extra_repulsion = OCbb_OCbb_exrep;
		}


		Real const bin = ( 4.2 * 4.2 / .05 ) + 1.0; //SGM Off-by-1 bug fix: Add + 1.0: Index 1 is at distance^2==0
		auto const ibin( static_cast< int >( bin ) );
		for ( int k = 1; k <= etable_disbins_; ++k ) {
			Real const dis = std::sqrt( ( k - 1 ) * .05f ); //SGM Off-by-1 bug fix: k -> ( k - 1 )

			// if( !SMOOTH_ETABLES ) {
			if ( dis < 4.2 ) {
				bool at1iscarbon(false), at2iscarbon(false);
				for ( Size i = 1; i <= carbon_types_.size(); ++i ) { at1iscarbon |= (atype1 == carbon_types_[i] ); }
				for ( Size i = 1; i <= carbon_types_.size(); ++i ) { at2iscarbon |= (atype2 == carbon_types_[i] ); }
				if ( at1iscarbon && at2iscarbon ) {
					fasol1(  k ) = fasol1( ibin );
					fasol2(  k ) = fasol2( ibin );
					dfasol(  k ) = 0.0;
					dfasol1( k ) = 0.0;
				}
			}
			//   push carbonyl oxygens (in beta sheets) apart.  a proxy for the missing
			//   electrostatic repulsion needed to counteract the LJatr which pulls the oxyens
			//   together
			if ( dis > 3.6 )  continue;

			if ( atype1 == OCbb_idx && atype2 == OCbb_idx && ! skip_mod_OCbb_OCbb_rep ) {
				Real const fac = std::max( dis - 3.6, -1.5 );
				//std::cout << "adding extra repulsion " << dis << " " << 2*fac*fac << " + " << ljrep_(k,OCbb_idx,OCbb_idx) <<  std::endl;
				//Real const fac = std::max( dis - 3.6f, -1.5f );
				ljrep(  k ) += 2 * ( fac * fac );
				dljrep( k ) += 4 * fac;
			}
		} // end  for ( int k = 1; k <= etable_disbins_; ++k ) {


	} //Objexx:SGM Extra {} scope is a VC++ work-around

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// zero out the attractive and solvation energies for the REPLS and HREPS atom types
	//
	// TR << "setting up REPLONLY residues to zero" << std::endl;
	chemical::AtomTypeSetCOP atom_set( atom_set_ );
	Size repls_idx = atom_set->atom_type_index( "REPLS" );
	Size hreps_idx = atom_set->atom_type_index( "HREPS" );
	if ( atype1 == repls_idx || atype2 == repls_idx || atype1 == hreps_idx || atype2 == hreps_idx ) {
		Size last_rep_bin( 0 );
		for ( int i = 1; i <= etable_disbins_; i++ ) {
			ljatr(   i ) = 0.0;
			dljatr(  i ) = 0.0;
			fasol1(  i ) = 0.0;
			fasol2(  i ) = 0.0;
			dfasol(  i ) = 0.0;
			dfasol1( i ) = 0.0;
			if ( last_rep_bin == 0 && ljrep( i ) == 0.0 ) {
				last_rep_bin = i;
			}
		}
		//std::cout << "zeroing ljatr_final_weight " << (*atom_set_)[at1].name() << " " << (*atom_set_)[at2].name() << std::endl;

		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.ljrep_from_negcrossing = true;
		// record the minimum value as the bin after the last value at which the repulsive energy is positive.
		p.lj_minimum = std::sqrt( ( last_rep_bin - 1 ) * 1.0 / bins_per_A2_ );
		p.maxd2 = p.lj_minimum*p.lj_minimum; // also set this value as the maximum distance for which the energy should be evaluated
		p.ljatr_final_weight = 0.0;
		p.fasol_final_weight = 0.0;
	}
}


void
Etable::smooth_etables_one_pair(
	Size const atype1,
	Size const atype2,
	ObjexxFCL::FArray1A< Real > ljatr,
	ObjexxFCL::FArray1A< Real > dljatr,
	ObjexxFCL::FArray1A< Real > ljrep,
	ObjexxFCL::FArray1A< Real > dljrep,
	ObjexxFCL::FArray1A< Real > fasol1,
	ObjexxFCL::FArray1A< Real > fasol2,
	ObjexxFCL::FArray1A< Real > dfasol,
	ObjexxFCL::FArray1A< Real > dfasol1
)
{
	using namespace numeric::interpolation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ljatr.dimension(   etable_disbins_ );
	dljatr.dimension(  etable_disbins_ );
	ljrep.dimension(   etable_disbins_ );
	dljrep.dimension(  etable_disbins_ );
	fasol1.dimension(  etable_disbins_ );
	fasol2.dimension(  etable_disbins_ );
	dfasol.dimension(  etable_disbins_ );
	dfasol1.dimension( etable_disbins_ );

	FArray1D<Real> dis(etable_disbins_);
	for ( int i = 1; i <= etable_disbins_; ++i ) {
		dis(i) = sqrt( ( i - 1 ) * 1.0 / bins_per_A2_ );
	}
	FArray1D<int> bin_of_dis(etable_disbins_);
	for ( int i = 1; i <= etable_disbins_; ++i ) {
		float d = ((float)i)/10.0;
		bin_of_dis(i) = (int)( d*d * bins_per_A2_ + 1 );
	}

	int minima_bin_index = 0;

	//////////////////////////////////////////////////////////////////
	// 1) change ljatr / ljrep split so that scaling is smooth
	//////////////////////////////////////////////////////////////////
	if ( atype1 == 1 && atype2 == 1 ) {
		TR << "smooth_etable: changing atr/rep split to bottom of energy well" << std::endl;
	}

	// find ljatr min
	core::Real min_atr = 0;
	int which_min = -1;
	for ( int i = 1; i <= etable_disbins_; ++i ) {
		if ( ljatr(i) < min_atr ) {
			min_atr = ljatr(i);
			which_min = i;
		}
	}
	if ( which_min != -1 ) {
		minima_bin_index = which_min;
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );

		// Correct version
		p.lj_minimum = lj_sigma_(atype1,atype2);
		p.lj_val_at_minimum = -1 * std::sqrt(lj_wdepth_[atype1]*lj_wdepth_[atype2]);

	} else {
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.lj_minimum = max_dis_;
		p.lj_val_at_minimum = 0;
	}
	// for dis below ljatr min, transfer ljatr to ljrep
	if ( min_atr < 0.0 ) {
		for ( int i = 1; i <= which_min; ++i ) {
			ljrep(i)  += ljatr(i) - min_atr;
			ljatr(i)  -= ljatr(i) - min_atr; // = min_atr;
			dljrep(i) += dljatr(i);
			dljatr(i) -= dljatr(i); // = 0;
		}
	}

	using namespace numeric::interpolation::spline;
	///////////////////////////////////////////////////////////////////
	// 2) spline smooth ljatr/ljrep at fa_max_dis cut
	//////////////////////////////////////////////////////////////////
	if ( atype1 == 1 && atype2 == 1 ) {
		TR << "smooth_etable: spline smoothing lj etables (maxdis = " << max_dis_ << ")" << std::endl;
	}

	// APL - 2012/6/14 - Take the bin maximum of the LJ-radii sum and (w/famxd=6) 4.5 A.  This change
	// effects the interaction of Br/Br and Br/I.
	int start = std::max( bin_of_dis( (int)( (max_dis_-1.5) * 10.0) ), minima_bin_index+1 );
	Real lbx  = dis(start);
	Real ubx  = dis(etable_disbins_);
	Real lby  =  ljatr(start);
	Real lbdy = dljatr(start);
	Real uby  = 0.0;
	Real ubdy = 0.0;
	{
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.ljatr_cubic_poly_xlo = lbx;
		p.ljatr_cubic_poly_xhi = ubx;
	}
	//std::cout << (*atom_set_)[at1].name() << " " << (*atom_set_)[at2].name() << " lj_minima: " << lj_minima(at1,at2) << " lj_vals_at_minima " << lj_vals_at_minima(at1,at2);
	//std::cout << " lbx " << lbx << " ubx " << ubx << std::endl;

	SplineGenerator gen( lbx, lby, lbdy, ubx, uby, ubdy );

	// APL -- Disabling this behavior since I think it is unused and if it is used, then it would
	// Prevent the analytic etable evaluation code I'm working on

	InterpolatorOP interp( gen.get_interpolator() );
	for ( int i = start; i <= etable_disbins_; ++i ) {
		interp->interpolate( dis(i), ljatr(i), dljatr(i) );
	}

	// Save the spline parameters for the analytic evaluation
	SimpleInterpolatorOP sinterp = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp );
	if ( ! sinterp ) {
		utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
	}
	{
		numeric::SplineParameters sparams;
		sparams.ylo  = sinterp->y()[ 1 ];
		sparams.yhi  = sinterp->y()[ 2 ];
		sparams.y2lo = sinterp->ddy()[ 1 ];
		sparams.y2hi = sinterp->ddy()[ 2 ];
		//ljatr_spline_parameters( atype1, atype2 ) = sparams;
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.ljatr_cubic_poly_parameters = cubic_polynomial_from_spline( p.ljatr_cubic_poly_xlo, p.ljatr_cubic_poly_xhi, sparams );
	}

	/////////////////////////////////////////////////////////////////////////////
	// 3) spline smooth solv1/solv2 at fa_max_dis cut && and first switchpoint
	/////////////////////////////////////////////////////////////////////////////
	if ( atype1 == 1 && atype2 == 1 ) {
		TR << "smooth_etable: spline smoothing solvation etables (max_dis = " << max_dis_ << ")" << std::endl;
	}

	// These two values do not depend on at1 or at2
	int const S2 = bin_of_dis( (int)((max_dis_-1.5)*10.0) ); // start of spline 2
	int const E2 = etable_disbins_; // end os spline 2

	// Save these values for the splines
	fasol_cubic_poly_far_xlo_ = dis( S2 );
	fasol_cubic_poly_far_xhi_ = dis( E2 );

	int SWTCH = 1;
	for ( SWTCH=1; SWTCH <= etable_disbins_; ++SWTCH ) {
		// std::cerr << SWTCH << "," << at1 << "," << at2 << "," << solv1_(SWTCH,at1,at2) << " ";
		if ( (fasol1(SWTCH) != fasol1(1)) ||
				(fasol2(SWTCH) != fasol2(1)) ) {
			break;
		}
	}
	if ( SWTCH > etable_disbins_ ) {
		// treat the range [0,famaxdis] as a constant (of zero) -- everything below
		// the distance fasol_spline_close_start_end(at1,at2).first is evaluated as
		// the constant fasol_spline_close(at1,at2).ylo.
		numeric::SplineParameters sp; sp.ylo=0; sp.yhi=0; sp.y2lo = 0; sp.y2hi = 0;
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.fasol_cubic_poly_close_start = fasol_cubic_poly_far_xhi_;
		p.fasol_cubic_poly_close_end   = fasol_cubic_poly_far_xhi_ + 1.0;
		p.fasol_cubic_poly_close_flat  = 0;
		p.fasol_cubic_poly_close = cubic_polynomial_from_spline( p.fasol_cubic_poly_close_start, p.fasol_cubic_poly_close_end, sp );
		return;
	}

	int const S1 = std::max(1,SWTCH - 30);
	int const E1 = std::min(SWTCH + 20,406);
	//std::cout << (*atom_set_)[atype1].name() << " " << (*atom_set_)[atype2].name() << " smooth solv " << S1 << " " << E1 << " " << S2 << " " << E2 << std::endl;

	Real dsolv1e1 = (fasol1(E1+1)-fasol1(E1  ))/(dis(E1+1)-dis(E1  ));
	Real dsolv1s2 = (fasol1(S2  )-fasol1(S2-1))/(dis(S2  )-dis(S2-1));
	Real dsolv2e1 = (fasol2(E1+1)-fasol2(E1  ))/(dis(E1+1)-dis(E1  ));
	Real dsolv2s2 = (fasol2(S2  )-fasol2(S2-1))/(dis(S2  )-dis(S2-1));

	SplineGenerator gen11( dis(S1), fasol1(S1), 0.0     , dis(E1), fasol1(E1), dsolv1e1 );
	SplineGenerator gen21( dis(S1), fasol2(S1), 0.0     , dis(E1), fasol2(E1), dsolv2e1 );
	SplineGenerator gen12( dis(S2), fasol1(S2), dsolv1s2, dis(E2), fasol1(E2), 0.0      );
	SplineGenerator gen22( dis(S2), fasol2(S2), dsolv2s2, dis(E2), fasol2(E2), 0.0      );

	InterpolatorOP interp11( gen11.get_interpolator() );
	InterpolatorOP interp21( gen21.get_interpolator() );
	for ( int i = S1; i <= E1; ++i ) {
		Real d1,d2;
		interp11->interpolate( dis(i), fasol1(i), d1 );
		interp21->interpolate( dis(i), fasol2(i), d2 );
		dfasol(i) = d1+d2;
		dfasol1(i) = d1;
	}

	InterpolatorOP interp12( gen12.get_interpolator() );
	InterpolatorOP interp22( gen22.get_interpolator() );
	for ( int i = S2; i <= E2; ++i ) {
		Real d1,d2;
		interp12->interpolate( dis(i), fasol1(i), d1 );
		interp22->interpolate( dis(i), fasol2(i), d2 );
		dfasol(i)  = d1+d2;
		dfasol1(i) = d1;
	}

	// Save the spline parameters
	// Represent the energy as simply the sum of the two splines: the polynomials may simply be added together.
	{
		SplineGenerator genclose( dis(S1), fasol1(S1) + fasol2(S1), 0, dis(E1), fasol1(E1)+fasol2(E1), dfasol(E1) );
		InterpolatorOP interp_close( genclose.get_interpolator() );

		SimpleInterpolatorOP sinterp_close = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_close );

		if ( ! sinterp_close ) {
			utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
		}
		numeric::SplineParameters sparams;
		sparams.ylo  = sinterp_close->y()[ 1 ];
		sparams.yhi  = sinterp_close->y()[ 2 ];
		sparams.y2lo = sinterp_close->ddy()[ 1 ];
		sparams.y2hi = sinterp_close->ddy()[ 2 ];

		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.fasol_cubic_poly_close_start = dis(S1);
		p.fasol_cubic_poly_close_end = dis(E1);
		p.fasol_cubic_poly_close_flat  = sparams.ylo;
		p.fasol_cubic_poly_close = cubic_polynomial_from_spline( p.fasol_cubic_poly_close_start, p.fasol_cubic_poly_close_end, sparams );
	}

	// Now compute the splines for the individual desolvation of atoms 1 and 2 between S1 and E1
	{
		Real dummy1, dummy2, dsolvE1, dsolvE2;
		lk_solv_energy_and_deriv( atype1, atype2, dis(E1), dummy1, dummy2, dsolvE1, dsolvE2 );

		// Create a close spline for the desolvation of atom 1 by atom 2
		{

			SplineGenerator genclose( dis(S1), fasol1(S1), 0, dis(E1), fasol1(E1), dsolvE1 );
			InterpolatorOP interp_close( genclose.get_interpolator() );

			SimpleInterpolatorOP sinterp_close = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_close );

			if ( ! sinterp_close ) {
				utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
			}
			numeric::SplineParameters sparams;
			sparams.ylo  = sinterp_close->y()[ 1 ];
			sparams.yhi  = sinterp_close->y()[ 2 ];
			sparams.y2lo = sinterp_close->ddy()[ 1 ];
			sparams.y2hi = sinterp_close->ddy()[ 2 ];

			EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
			p.fasol_cubic_poly1_close_flat  = sparams.ylo;
			p.fasol_cubic_poly1_close = cubic_polynomial_from_spline( p.fasol_cubic_poly_close_start, p.fasol_cubic_poly_close_end, sparams );
		}

		// Create a close spline for the desolvation of atom 2 by atom 1
		{
			SplineGenerator genclose( dis(S1), fasol2(S1), 0, dis(E1), fasol2(E1), dsolvE2 );
			InterpolatorOP interp_close( genclose.get_interpolator() );

			SimpleInterpolatorOP sinterp_close = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_close );

			if ( ! sinterp_close ) {
				utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
			}
			numeric::SplineParameters sparams;
			sparams.ylo  = sinterp_close->y()[ 1 ];
			sparams.yhi  = sinterp_close->y()[ 2 ];
			sparams.y2lo = sinterp_close->ddy()[ 1 ];
			sparams.y2hi = sinterp_close->ddy()[ 2 ];

			EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
			p.fasol_cubic_poly2_close_flat  = sparams.ylo;
			p.fasol_cubic_poly2_close = cubic_polynomial_from_spline( p.fasol_cubic_poly_close_start, p.fasol_cubic_poly_close_end, sparams );
		}
	}

	{
		// 2. far spline
		// APL Smoother derivatives -- create a spline where the derivatives are actually correct
		SplineGenerator genfar( dis(S2), fasol1(S2) + fasol2(S2), dfasol(S2), dis(E2), 0.0, 0.0 );
		InterpolatorOP interp_far( genfar.get_interpolator() );

		SimpleInterpolatorOP sinterp_far = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_far );
		if ( ! sinterp_far ) {
			utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
		}
		numeric::SplineParameters sparams;
		sparams.ylo  = sinterp_far->y()[ 1 ];
		sparams.yhi  = sinterp_far->y()[ 2 ];
		sparams.y2lo = sinterp_far->ddy()[ 1 ];
		sparams.y2hi = sinterp_far->ddy()[ 2 ];
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.fasol_cubic_poly_far = cubic_polynomial_from_spline( fasol_cubic_poly_far_xlo_, fasol_cubic_poly_far_xhi_,  sparams );
	}

	// Now compute the splines for the individual desolvation of atoms 1 and 2 between S2 and E2
	{
		Real dummy1, dummy2, dsolvE1, dsolvE2;
		lk_solv_energy_and_deriv( atype1, atype2, dis(S2), dummy1, dummy2, dsolvE1, dsolvE2 );

		{
			// Create a spline for the desolvation of atom 1 by atom 2
			SplineGenerator genfar( dis(S2), fasol1(S2) , dsolvE1, dis(E2), 0.0, 0.0 );
			InterpolatorOP interp_far( genfar.get_interpolator() );

			SimpleInterpolatorOP sinterp_far = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_far );
			if ( ! sinterp_far ) {
				utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
			}
			numeric::SplineParameters sparams;
			sparams.ylo  = sinterp_far->y()[ 1 ];
			sparams.yhi  = sinterp_far->y()[ 2 ];
			sparams.y2lo = sinterp_far->ddy()[ 1 ];
			sparams.y2hi = sinterp_far->ddy()[ 2 ];
			EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
			p.fasol_cubic_poly1_far = cubic_polynomial_from_spline( fasol_cubic_poly_far_xlo_, fasol_cubic_poly_far_xhi_, sparams );
		}

		// Create a spline for the desolvation of atom 2 by atom 1
		{
			SplineGenerator genfar( dis(S2), fasol2(S2), dsolvE2, dis(E2), 0.0, 0.0 );
			InterpolatorOP interp_far( genfar.get_interpolator() );

			SimpleInterpolatorOP sinterp_far = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_far );
			if ( ! sinterp_far ) {
				utility_exit_with_message( "Etable created non-simple-interpolator in smooth_etables()" );
			}
			numeric::SplineParameters sparams;
			sparams.ylo  = sinterp_far->y()[ 1 ];
			sparams.yhi  = sinterp_far->y()[ 2 ];
			sparams.y2lo = sinterp_far->ddy()[ 1 ];
			sparams.y2hi = sinterp_far->ddy()[ 2 ];
			EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
			p.fasol_cubic_poly2_far = cubic_polynomial_from_spline( fasol_cubic_poly_far_xlo_, fasol_cubic_poly_far_xhi_, sparams );
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief output an etable data file in the same format used in input_etable
///
/// @details
///$$$ file first line is <etable> <etable_disbins_>
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
/// @author sheffler mar 19 2006
///
/////////////////////////////////////////////////////////////////////////////////
void
Etable::output_etable(
	ObjexxFCL::FArray3D<Real> & etable,
	std::string label,
	std::ostream & out
) {
	using namespace std;
	using namespace ObjexxFCL;

	out << label << " " << etable_disbins_ << endl;
	for ( int at1 = 1; at1 <= n_atomtypes_; at1++ ) {
		for ( int at2 = 1; at2 <= n_atomtypes_; at2++ ) {
			out << at1 << " "
				<< at2 << " ";
			for ( int bin = 1; bin <= etable_disbins_; bin++ ) {
				float evalue = etable(bin,at1,at2);
				out << evalue << ' ';
			}
			out << endl;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief read in etable from a datafile
///
/// @details
///$$$ file first line is <etable> <etable_disbins_>
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
/// @author sheffler mar 19 2006
///
/////////////////////////////////////////////////////////////////////////////////
void
Etable::input_etable(
	ObjexxFCL::FArray3D<Real> & etable,
	std::string const & label,
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
	if ( !(intmp >> lblin >> numdisbin) ) {
		TR << "input_etable: WARNING bad etable header " << buf << endl;
		return;
	}
	if ( lblin != label || etable_disbins_ != numdisbin ) {
		TR << "input_etable: WARNING etable types don't match! "<< endl;
		TR << "              expected " << label << "," << etable_disbins_
			<< " got " << lblin << ',' << numdisbin<< endl;
		utility_exit_with_message( "input_etable: WARNING etable types don't match! " );
	} else {
		TR << "input_etable expected etable " << label << " of size " << etable_disbins_
			<< ", got " << lblin << ',' << numdisbin<< endl;

	}

	// read in etable
	int at1,at2,count=0,scount=0;
	float evalue;
	while ( in.getline(buf,100000) ) {
		//TR << "ASDFF buf  " << buf << endl;
		count++;
		intmp.clear();
		intmp.str(buf);
		if ( ! (intmp >> at1 >> at2 ) ) {
			TR << "input_etable: error reading etable line: " << buf << endl;
		} else {
			for ( int bin = 1; bin <=numdisbin; bin++ ) {
				if ( !(intmp >> evalue) ) {
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
///
/// @brief precalculate non-distance dependent coefficients of energy functions
///
/// @details
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
/// @author ctsa 10-2003
///
/////////////////////////////////////////////////////////////////////////////////
void
Etable::precalc_etable_coefficients(
	FArray2< Real > & lj_sigma,
	FArray2< Real > & lj_r6_coeff,
	FArray2< Real > & lj_r12_coeff,
	FArray2< Real > & lj_switch_intercept,
	FArray2< Real > & lj_switch_slope,
	FArray1< Real > & lk_inv_lambda2,
	FArray2< Real > & lk_coeff,
	FArray2< Real > & lk_min_dis2sigma_value
)
{

	//  using namespace etable;
	//  using namespace param;
	//  using namespace water;
	//   using namespace hbonds;

	// locals
	Real sigma,sigma6,sigma12,wdepth;
	Real inv_lambda;
	FArray1D< Real > lk_coeff_tmp( n_atomtypes_ );
	Real thresh_dis,inv_thresh_dis2,x_thresh;
	//int dstype;
	//int const atype_sulfur = { 16 };
	Real const inv_neg2_tms_pi_sqrt_pi = { -0.089793561062583294 };
	// coefficient for lk solvation

	// include follows locals so that data statements can initialize included arrays
	for ( int i = 1, e = n_atomtypes_; i <= e; ++i ) {
		inv_lambda = 1.0/lk_lambda(i);
		//inv_lambda = 1.0/atom_type(i).lk_lambda();
		lk_inv_lambda2(i) = inv_lambda * inv_lambda;
		lk_coeff_tmp(i) = inv_neg2_tms_pi_sqrt_pi *
			lk_dgfree(i) * inv_lambda;
		//atom_type(i).lk_dgfree() * inv_lambda;
	}

	for ( int i = 1, e = n_atomtypes_; i <= e; ++i ) {
		for ( int j = i; j <= e; ++j ) {

			sigma = Wradius_ * ( lj_radius(i) + lj_radius(j) );
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

			if ( lj_use_hbond_radii_ ) {
				if ( ( atom_type(i).is_acceptor() && atom_type(j).is_donor() ) ||
						( atom_type(i).is_donor() && atom_type(j).is_acceptor() ) ) {
					// fpd this is bad ... can this be a property instead?
					// apl -- yep, it sure is
					if (
							( atom_type(j).is_donor() && atom_type(j).name().substr(0,2) == "OH" ) ||
							( atom_type(i).is_donor() && atom_type(i).name().substr(0,2) == "OH" ) ||
							( atom_type(j).is_donor() && atom_type(j).name().substr(0,2) == "OW" ) ||
							( atom_type(i).is_donor() && atom_type(i).name().substr(0,2) == "OW" ) ||
							( atom_type(j).is_donor() && atom_type(j).name() == "Oet3" ) || // "Oet3" is ester atom (used for nucleic acid O4', used to be "OH")
							( atom_type(i).is_donor() && atom_type(i).name() == "Oet3" ) ) {
						sigma = lj_hbond_OH_donor_dis_;
					} else {
						sigma = lj_hbond_dis_;
						//           if (tight_hb && ((i == 15 && j == 13) || (j == 15  && i == 13)))
						//       sigma = lj_hbond_accOch_dis;
					}
				} else if ( ( atom_type(i).is_acceptor() && atom_type(j).is_polar_hydrogen() ) ||
						( atom_type(i).is_polar_hydrogen() && atom_type(j).is_acceptor() ) ) {
					sigma = lj_hbond_hdis_;
				} else if ( ( atom_type(i).is_acceptor() && atom_type(j).name() == "MG2p" ) ||
						( atom_type(i).name() == "MG2p" && atom_type(j).is_acceptor() ) ) {
					sigma = lj_hbond_hdis_;
					//           if (tight_hb && ((i == 15 && j == 22) || (j == 15 && i == 22)))
					//       sigma = lj_hbond_accOch_hdis;
				}
			}

			//lin   modify sigma for water and hbond donors/acceptors
			if ( lj_use_water_radii_ ) {
				if ( ( ( atom_type(i).is_acceptor() ||
						atom_type(i).is_donor() ) &&
						atom_type(j).is_h2o() ) ||
						( ( atom_type(j).is_acceptor() ||
						atom_type(j).is_donor() ) &&
						atom_type(i).is_h2o() ) ) {
					sigma = lj_water_dis_;
				} else if ( ( atom_type(i).is_polar_hydrogen() &&
						atom_type(j).is_h2o() ) ||
						( atom_type(j).is_polar_hydrogen() &&
						atom_type(i).is_h2o() ) ) {
					sigma = lj_water_hdis_;
				}
			}

			sigma6  = std::pow( sigma, 6 );
			sigma12 = sigma6 * sigma6;
			wdepth = std::sqrt(lj_wdepth(i)*lj_wdepth(j));
			//wdepth = std::sqrt(atom_type(i).lj_wdepth()*atom_type(j).lj_wdepth());

			lj_sigma(i,j) = sigma;
			lj_sigma(j,i) = lj_sigma(i,j);

			lj_r6_coeff(i,j) = -2. * wdepth * sigma6;
			lj_r6_coeff(j,i) = lj_r6_coeff(i,j);

			lj_r12_coeff(i,j) = wdepth * sigma12;
			lj_r12_coeff(j,i) = lj_r12_coeff(i,j);

			// ctsa - create coefficients for linear projection of lj repulsive used
			//  for low distance values
			if ( lj_use_lj_deriv_slope_ ) {

				// ctsa - use the slope of the true lj repulsive at the
				//  linear switch point to create a linear projection of
				//  lj for low distances

				//  slope = wdepth/sigma *
				//          (slope@switch_point*sigma/wdepth)
				lj_switch_slope(i,j) = (wdepth/sigma)*
					lj_switch_slope_sigma2wdepth_;
				lj_switch_slope(j,i) = lj_switch_slope(i,j);

				// intercept = wdepth*(lj@switch_point/wdepth)
				//             - slope*switch_point_distance
				lj_switch_intercept(i,j) = wdepth*lj_switch_value2wdepth_ -
					lj_switch_slope(i,j)*sigma*lj_switch_dis2sigma_;
				lj_switch_intercept(j,i) = lj_switch_intercept(i,j);
			} else {

				// ctsa - create a linear projection of lj for low distances which
				//  is defined by a constant y intercept and the true lj repulsive
				//  value at the linear switch point
				lj_switch_slope(i,j) = -(1./sigma)*lj_switch_sigma2dis_*
					(lj_slope_intercept_-wdepth*lj_switch_value2wdepth_);
				lj_switch_slope(j,i) = lj_switch_slope(i,j);

				lj_switch_intercept(i,j) = lj_slope_intercept_;
				lj_switch_intercept(j,i) = lj_switch_intercept(i,j);
			}

			// ctsa - precalculated lk solvation coefficients
			lk_coeff(i,j) = lk_coeff_tmp(i) * lk_volume(j);
			lk_coeff(j,i) = lk_coeff_tmp(j) * lk_volume(i);

			// ctsa - when dis/sigma drops below lk_min_dis2sigma,
			//   a constant lk solvation value equal to the value at the
			//   switchover point is used. That switchover-point value
			//   is calculated here and stored in lk_min_dis2sigma_value
			thresh_dis = lk_min_dis2sigma_*sigma;
			inv_thresh_dis2 = 1./( thresh_dis * thresh_dis );
			Real dis_rad = thresh_dis - lj_radius(i);
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(i);
			lk_min_dis2sigma_value(i,j) = std::exp(-x_thresh) * lk_coeff(i,j) *
				inv_thresh_dis2;

			dis_rad = thresh_dis - lj_radius(j);
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(j);
			lk_min_dis2sigma_value(j,i) = std::exp(-x_thresh) * lk_coeff(j,i) *
				inv_thresh_dis2;

		}
	}
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief calc all etable values given a distance and atom-type pair
///
/// @details
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
/// @author ctsa 10-2003
///
/////////////////////////////////////////////////////////////////////////////////
void
Etable::calc_etable_value(
	Real dis2,
	int const atype1,
	int const atype2,
	Real & atrE,
	Real & d_atrE,
	Real & repE,
	Real & d_repE,
	Real & solvE1,
	Real & solvE2,
	Real & dsolvE1,
	Real & dsolvE2
) const
{
	//  using namespace etable;
	//  using namespace param;
	//  using namespace pdbstatistics_pack;

	// locals
	Real ljE,d_ljE;
	Real dis;
	Real inv_dis,inv_dis2,inv_dis6,inv_dis7,inv_dis12,inv_dis13;
	Real dis2sigma;

	// include after local variables to allow data statements to initialize
	atrE = 0.;
	d_atrE = 0.;
	repE = 0.;
	d_repE = 0.;
	solvE1 = 0.;
	solvE2 = 0.;
	dsolvE1 = 0.;
	dsolvE2 = 0.;

	//  ctsa - epsilon_ allows final bin value to be calculated
	if ( dis2 > max_dis2_ + epsilon_ ) return;

	if ( dis2 < min_dis2_ ) dis2 = min_dis2_;

	dis = std::sqrt(dis2);
	inv_dis = 1.0/dis;
	inv_dis2 = inv_dis * inv_dis;

	dis2sigma = dis / lj_sigma_(atype1,atype2);

	if ( dis2sigma < lj_switch_dis2sigma_ ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma
		//    ratio drops below theshold
		d_ljE = lj_switch_slope_(atype1,atype2);
		ljE = dis*d_ljE + lj_switch_intercept_(atype1,atype2);
	} else {
		//  ctsa - calc regular lennard-jones
		inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
		inv_dis7  = inv_dis6 * inv_dis;
		inv_dis12 = inv_dis6 * inv_dis6;
		inv_dis13 = inv_dis12 * inv_dis;

		ljE = lj_r12_coeff_(atype1,atype2) * inv_dis12 +
			lj_r6_coeff_(atype1,atype2) * inv_dis6;

		d_ljE = -12.*lj_r12_coeff_(atype1,atype2) * inv_dis13-6. *
			lj_r6_coeff_(atype1,atype2) * inv_dis7;
	}

	if ( ljE < 0. ) {
		atrE = ljE;
		d_atrE = d_ljE;
	} else {
		repE = ljE;
		d_repE = d_ljE;
	}


	// ctsa - calc lk
	if ( dis2sigma < lk_min_dis2sigma_ ) {
		// ctsa - solvation is constant when the dis/sigma ratio
		//   falls below minimum threshold
		solvE1 = lk_min_dis2sigma_value_(atype1,atype2);
		solvE2 = lk_min_dis2sigma_value_(atype2,atype1);
		dsolvE1 = dsolvE2 = 0.0;

	} else {
		lk_solv_energy_and_deriv( atype1, atype2, dis, solvE1, solvE2, dsolvE1, dsolvE2 );
	}

}

void
Etable::lk_solv_energy_and_deriv(
	int const atype1,
	int const atype2,
	Real const dis,
	Real & solvE1,
	Real & solvE2,
	Real & dsolvE1,
	Real & dsolvE2
) const
{
	Real const inv_dis = 1.0/dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	Real const dis_rad1 = dis - lj_radius( atype1 );
	Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);

	Real const dis_rad2 = dis - lj_radius( atype2 );
	Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

	solvE1 = std::exp(-x1) * lk_coeff_(atype1,atype2) * inv_dis2;
	solvE2 = std::exp(-x2) * lk_coeff_(atype2,atype1) * inv_dis2;

	// ctsa - get d(lk_E)/dr
	dsolvE1 = -2.0 * solvE1 *
		(((dis-lj_radius(atype1))*lk_inv_lambda2_(atype1))+inv_dis);
	dsolvE2 = -2.0 * solvE2 *
		(((dis-lj_radius(atype2))*lk_inv_lambda2_(atype2))+inv_dis);
}

void
Etable::zero_hydrogen_and_water_ljatr_one_pair(
	Size const atype1,
	Size const atype2,
	ObjexxFCL::FArray1A< Real > ljrep,
	ObjexxFCL::FArray1A< Real > ljatr,
	ObjexxFCL::FArray1A< Real > dljatr,
	ObjexxFCL::FArray1A< Real > fasol1,
	ObjexxFCL::FArray1A< Real > fasol2,
	ObjexxFCL::FArray1A< Real > dfasol,
	ObjexxFCL::FArray1A< Real > dfasol1
)
{
	ljrep.dimension(  etable_disbins_ );
	ljatr.dimension(  etable_disbins_ );
	dljatr.dimension( etable_disbins_ );
	fasol1.dimension( etable_disbins_ );
	fasol2.dimension( etable_disbins_ );
	dfasol.dimension( etable_disbins_ );
	dfasol1.dimension( etable_disbins_ );

	chemical::AtomTypeSetCOP atom_set( atom_set_ );
	Size const HOH = atom_set->atom_type_index("HOH");
	if ( atype1 == HOH || atype2 == HOH || atom_type(atype1).is_hydrogen() || atom_type(atype2).is_hydrogen() ) {

		// cbk  don't give hydrogens or water attractive lennard-jones
		// cbk  this is so very short range cut-offs can be used

		Size first_zero_ljrep = 0;
		for ( int i = 1; i <= etable_disbins_; ++i ) {
			ljatr(i) = 0.0;
			dljatr(i) = 0.0;
			fasol1(i) = fasol2(i) = dfasol(i) = dfasol1(i) = 0.0; // APL TEMP.  Disable solvation term for hydrogen atoms
			if ( first_zero_ljrep == 0 && ljrep(i) == 0.0 ) {
				first_zero_ljrep = i;
			}
		}
		EtableParamsOnePair & p = analytic_params_for_pair( atype1, atype2 );
		p.ljatr_final_weight = 0.0;
		p.fasol_final_weight = 0.0;
		p.maxd2 = lj_sigma_(atype1,atype2) * lj_sigma_(atype1,atype2); // ( first_zero_ljrep - 1 ) * 1.0 / bins_per_A2;
		p.hydrogen_interaction = true;
		// Disable fasol for hydrogens? Technically, fasol doesn't get disabled for hydrogens! fasol_final_weight(at1,at2) = 0.0;
		// This is surely a bug.
	}
}

/// @brief Returns the maximum lj radius for any non-hydrogen
/// atom as defined by the atom-type-set used to create this Etable.
Real
Etable::max_non_hydrogen_lj_radius() const
{
	return max_non_hydrogen_lj_radius_;
}

/// @brief Returns the maximum lj radius for any hydrogen atom as
/// defined by the input atom-type-set used to create this Etable.
Real
Etable::max_hydrogen_lj_radius() const
{
	return max_hydrogen_lj_radius_;
}

void Etable::interpolated_analytic_etable_evaluation(
	core::conformation::Atom const & at1,
	core::conformation::Atom const & at2,
	core::Real & lj_atrE,
	core::Real & lj_repE,
	core::Real & fa_solE,
	Real & dis2
) const
{
	Real ljatrE_lo, ljrepE_lo, fasolE_lo;
	Real ljatrE_hi, ljrepE_hi, fasolE_hi;
	core::conformation::Atom at2_lo( at2 ), at2_hi( at2 );
	dis2 = at1.xyz().distance_squared( at2.xyz() );
	Real d2dummy;
	auto dis2bin = (int) ( dis2 * bins_per_A2_ );

	Real dis2lo = Real(dis2bin) / bins_per_A2_;
	Real dis2hi = Real(dis2bin  + 1 )/bins_per_A2_;
	at2_lo.xyz( Vector( std::sqrt( dis2lo ), 0, 0 ));
	at2_hi.xyz( Vector( std::sqrt( dis2hi ), 0, 0 ));

	analytic_etable_evaluation( at1, at2_lo, ljatrE_lo, ljrepE_lo, fasolE_lo, d2dummy );
	analytic_etable_evaluation( at1, at2_hi, ljatrE_hi, ljrepE_hi, fasolE_hi, d2dummy );

	Real alpha = (dis2*bins_per_A2_ - dis2bin);
	//std::cout << "debug " << dis2 << " " << dis2lo << " " << dis2hi << " " << alpha << " " << 1-alpha << " " << ljatrE_lo << " " << ljatrE_hi << std::endl;

	lj_atrE = alpha * ljatrE_hi + (1-alpha) * ljatrE_lo;
	lj_repE = alpha * ljrepE_hi + (1-alpha) * ljrepE_lo;
	fa_solE = alpha * fasolE_hi + (1-alpha) * fasolE_lo;
}


} // etable
} // scoring
} // core
