// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/memb_etable/MembEtable.cc
/// @brief Table of pre-computed LK membrane solvation energies
/// @author Patrick Barth (original)
/// @author Rebecca Alford (rfalford12@gmail.com)

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

// C++ Headers
#include <iostream>
#include <fstream>

#include <basic/Tracer.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "core.scoring.etable" );

namespace core {
namespace scoring {
namespace etable {

using namespace basic::options;
using namespace basic::options::OptionKeys;

// Constructors & Setup
MembEtable::MembEtable(
	chemical::AtomTypeSetCAP atom_set_in_ap,
	EtableOptions const & options,
	std::string const alternate_parameter_set
) : Etable( atom_set_in_ap, options, alternate_parameter_set ),
	max_non_hydrogen_lj_radius_( 0.0 ),
	max_hydrogen_lj_radius_( 0.0 ),
	lj_radius_(),
	lk_dgfree_(),
	lk_volume_(),
	lk_lambda_(),
	memb_lk_dgfree_(),
	lk_dgrefce_(),
	memb_lk_dgrefce_(),
	solv1_(),
	solv2_(),
	dsolv1_(),
	dsolv2_(),
	memb_solv1_(),
	memb_solv2_(),
	memb_dsolv1_(),
	memb_dsolv2_()
{

	dimension_memb_etable_arrays();
	initialize_from_atomset( atom_set_in_ap );
	make_pairenergy_table();
}

MembEtable::MembEtable( MembEtable const & src ) :
	Etable( src ),
	max_non_hydrogen_lj_radius_( src.max_hydrogen_lj_radius_ ),
	max_hydrogen_lj_radius_( src.max_hydrogen_lj_radius_ ),
	lj_radius_( src.lj_radius_ ),
	lk_dgfree_( src.lk_dgfree_ ),
	lk_volume_( src.lk_volume_ ),
	lk_lambda_( src.lk_lambda_ ),
	memb_lk_dgfree_( src.memb_lk_dgfree_ ),
	lk_dgrefce_( src.lk_dgrefce_ ),
	memb_lk_dgrefce_( src.memb_lk_dgrefce_ ),
	solv1_( src.solv1_ ),
	solv2_( src.solv2_ ),
	dsolv1_( src.dsolv1_ ),
	dsolv2_( src.dsolv2_ ),
	memb_solv1_( src.memb_solv1_ ),
	memb_solv2_( src.memb_solv2_ ),
	memb_dsolv1_( src.memb_dsolv1_ ),
	memb_dsolv2_( src.memb_dsolv2_ )
{}

MembEtable::~MembEtable() {}

// Initialization Functions

void
MembEtable::dimension_memb_etable_arrays() {

	// size the arrays
	solv1_.dimension( etable_disbins(), n_atomtypes(), n_atomtypes() );
	solv2_.dimension( etable_disbins(), n_atomtypes(), n_atomtypes() );
	memb_solv1_.dimension(  etable_disbins(), n_atomtypes(), n_atomtypes() );
	memb_solv2_.dimension(  etable_disbins(), n_atomtypes(), n_atomtypes() );
	dsolv2_.dimension(  etable_disbins(), n_atomtypes(), n_atomtypes() );
	dsolv1_.dimension( etable_disbins(), n_atomtypes(), n_atomtypes() );
	memb_dsolv2_.dimension(  etable_disbins(), n_atomtypes(), n_atomtypes() );
	memb_dsolv1_.dimension( etable_disbins(), n_atomtypes(), n_atomtypes() );

	lj_radius_.resize( n_atomtypes(), 0.0 );
	lk_dgfree_.resize( n_atomtypes(), 0.0 );
	lk_lambda_.resize( n_atomtypes(), 0.0 );
	lk_volume_.resize( n_atomtypes(), 0.0 );
	memb_lk_dgfree_.resize( n_atomtypes(), 0.0 );
	lk_dgrefce_.dimension( n_atomtypes() );
	memb_lk_dgrefce_.dimension( n_atomtypes() );
}

void
MembEtable::initialize_from_atomset( chemical::AtomTypeSetCAP atom_set_in_ap ) {

	chemical::AtomTypeSetCOP atom_set_in( atom_set_in_ap );

	for ( int i=1; i<= n_atomtypes(); ++i ) {
		lj_radius_[i] = (*atom_set_in)[i].lj_radius();
		lk_lambda_[i] = (*atom_set_in)[i].lk_lambda();
		lk_volume_[i] = (*atom_set_in)[i].lk_volume();

		if ( (*atom_set_in)[i].is_hydrogen() ) {
			if ( lj_radius_[i] > max_hydrogen_lj_radius_ ) max_hydrogen_lj_radius_ = lj_radius_[i];
		} else {
			if ( lj_radius_[i] > max_non_hydrogen_lj_radius_ ) max_non_hydrogen_lj_radius_ = lj_radius_[i];
		}
	}

	std::string param_name;

	param_name = "LK_DGFREE";
	TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
	Size const lkfree_index( atom_set_in->extra_parameter_index( param_name ) );
	for ( int i=1; i<= n_atomtypes(); ++i ) lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( lkfree_index );


	param_name = "MEMB_LK_DGFREE";
	TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
	Size const mb_lkfree_index( atom_set_in->extra_parameter_index( param_name ) );
	for ( int i=1; i<= n_atomtypes(); ++i ) memb_lk_dgfree_[i] = (*atom_set_in)[i].extra_parameter( mb_lkfree_index );

	param_name = "LK_DGREFCE";
	TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
	Size const lkref_index( atom_set_in->extra_parameter_index( param_name ) );
	for ( int i=1; i<= n_atomtypes(); ++i ) lk_dgrefce_(i) = (*atom_set_in)[i].extra_parameter( lkref_index );


	param_name = "MEMB_LK_DGREFCE";
	TR << "Using alternate parameters: " << param_name << " in MembEtable construction." << std::endl;
	Size const mb_lkref_index( atom_set_in->extra_parameter_index( param_name ) );
	for ( int i=1; i<= n_atomtypes(); ++i ) memb_lk_dgrefce_(i) = (*atom_set_in)[i].extra_parameter( mb_lkref_index );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Make Pair Energy Table: Compute Fast Lookup for vdw and solvation energy
/// @details Several energies are precomputed when fullatom mode is initialized and
/// stored in lookup tables to speed fullatom calculation. Currently
/// pre-computed values are the Lennard-Jones van der Waals approximation
/// (lj) and the Lazaridis-Karplus implicit solvation function (lk). For
/// each of these energies the derivative w.r.t. atom pair separation distance is
/// calculated and stored as well. Note that the lj energy is artificially
/// divided into atractive and repulsive components.
///
/// @author ctsa (last modified: 10-2003)
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::make_pairenergy_table()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real dis2,dis,dis2_step;
	Real solvE1,solvE2,dsolvE1,dsolvE2;
	Real memb_solvE1,memb_solvE2,memb_dsolvE1,memb_dsolvE2;

	// Pre-Computed coefficients for high-speed potential calculation computed
	// in intiailization. Note: +1 atom type is a disulfide bonded cyzteine
	// (only if using the old dissulfide hack)
	FArray2D< Real > lj_sigma( n_atomtypes() + 1, n_atomtypes() + 1 );
	FArray1D< Real > lk_inv_lambda2( n_atomtypes() );
	FArray2D< Real > lk_coeff( n_atomtypes(), n_atomtypes() );
	FArray2D< Real > lk_min_dis2sigma_value( n_atomtypes() + 1, n_atomtypes() + 1 );
	FArray2D< Real > memb_lk_coeff( n_atomtypes(), n_atomtypes() );
	FArray2D< Real > memb_lk_min_dis2sigma_value( n_atomtypes() + 1, n_atomtypes() + 1 );

	// Parameters to compute damping at max_dis range
	Real damping_thresh_dis2;
	Real dsolv1_damp, dsolv2_damp, intercept_solv1_damp, intercept_solv2_damp;
	Real memb_dsolv1_damp, memb_dsolv2_damp, intercept_memb_solv1_damp, intercept_memb_solv2_damp;

	TR << "Starting membrane specific energy table calculation" << std::endl;

	precalc_etable_coefficients( lj_sigma, lk_inv_lambda2, lk_coeff, memb_lk_coeff, lk_min_dis2sigma_value, memb_lk_min_dis2sigma_value);

	// calc distance**2 step per bin
	dis2_step = 1.0 / get_bins_per_A2();

	int normal_disbins;
	// get the number of damping disbins
	if ( add_long_range_damping() ) {
		Real const dif = max_dis() - long_range_damping_length();
		damping_thresh_dis2 = max_dis2() - ( dif * dif );
		int damping_disbins = static_cast< int >( damping_thresh_dis2*get_bins_per_A2() );
		normal_disbins = etable_disbins()-damping_disbins;
	} else {
		normal_disbins = etable_disbins();
	}

	// ctsa - step through distance**2 bins and calculate potential
	for ( int atype1 = 1, atype_end = n_atomtypes(); atype1 <= atype_end; ++atype1 ) {
		for ( int atype2 = 1; atype2 <= atype_end; ++atype2 ) {

			// ctsa - normal bins have their lj and lk values
			// calculated analytically
			for ( int disbin = 1; disbin <= normal_disbins; ++disbin ) {
				dis2 = ( disbin - 1 ) * dis2_step;

				calc_etable_value(dis2,atype1,atype2,solvE1,solvE2,dsolvE1,dsolvE2,lj_sigma,lk_inv_lambda2,
					lk_coeff, lk_min_dis2sigma_value,memb_solvE1, memb_solvE2,
					memb_lk_coeff, memb_lk_min_dis2sigma_value, memb_dsolvE1, memb_dsolvE2);

				solv1_(disbin,atype2,atype1) = solvE1;
				solv2_(disbin,atype2,atype1) = solvE2;
				dsolv2_(disbin,atype2,atype1) = dsolvE2;
				dsolv1_(disbin,atype2,atype1) = dsolvE1;
				memb_solv1_(disbin,atype2,atype1) = memb_solvE1;
				memb_solv2_(disbin,atype2,atype1) = memb_solvE2;
				memb_dsolv2_(disbin,atype2,atype1) = memb_dsolvE2;
				memb_dsolv1_(disbin,atype2,atype1) = memb_dsolvE1;
			}

			if ( add_long_range_damping() ) {

				dsolv1_damp = -solv1_(normal_disbins,atype2,atype1) /
					long_range_damping_length();
				dsolv2_damp = -solv2_(normal_disbins,atype2,atype1) /
					long_range_damping_length();

				memb_dsolv1_damp = -memb_solv1_(normal_disbins,atype2,atype1) /
					long_range_damping_length();
				memb_dsolv2_damp = -memb_solv2_(normal_disbins,atype2,atype1) /
					long_range_damping_length();


				intercept_solv1_damp = -dsolv1_damp*max_dis();
				intercept_solv2_damp = -dsolv2_damp*max_dis();

				intercept_memb_solv1_damp = -memb_dsolv1_damp*max_dis();
				intercept_memb_solv2_damp = -memb_dsolv2_damp*max_dis();

				for ( int disbin = normal_disbins+1; disbin <= etable_disbins(); ++disbin ) {
					dis2 = ( disbin - 1 ) * dis2_step;
					dis = std::sqrt(dis2);

					solv1_(disbin,atype2,atype1) = intercept_solv1_damp + dis *dsolv1_damp;
					solv2_(disbin,atype2,atype1) = intercept_solv2_damp + dis *dsolv2_damp;
					memb_solv1_(disbin,atype2,atype1) = intercept_memb_solv1_damp + dis * memb_dsolv1_damp;
					memb_solv2_(disbin,atype2,atype1) = intercept_memb_solv2_damp + dis * memb_dsolv2_damp;

					dsolv1_(disbin,atype2,atype1)  = dsolv1_damp;
					dsolv2_(disbin,atype2,atype1) = dsolv2_damp;
					memb_dsolv2_(disbin,atype2,atype1) = memb_dsolv2_damp;
					memb_dsolv1_(disbin,atype2,atype1) = memb_dsolv1_damp;

				}
			}

			solv1_(etable_disbins(),atype2,atype1) = 0.0;
			solv2_(etable_disbins(),atype2,atype1) = 0.0;
			memb_solv1_(etable_disbins(),atype2,atype1) = 0.0;
			memb_solv2_(etable_disbins(),atype2,atype1) = 0.0;
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
	modify_pot();

	if ( !option[ score::no_smooth_etables ] ) smooth_etables();

	////////////////////////////////////////////////////////
	// etable I/O stuff
	////////////////////////////////////////////////////////
	using namespace std;
	using namespace ObjexxFCL;

	// just for convenience in input/output_etables
	map< string, FArray3D<Real>* > etables;

	etables[ "solv1"] = &  solv1_; etables[ "solv2"] = &  solv2_;
	etables["dsolv1"]  = & dsolv1_;
	etables["dsolv2"]  = & dsolv2_;
	etables[ "memb_solv1"] = &  memb_solv1_;  etables[ "memb_solv2"] = &  memb_solv2_;
	etables["memb_dsolv2"]  = & memb_dsolv2_; etables["memb_dsolv1"]  = & memb_dsolv1_;

	if ( option[ score::input_etables ].user() ) {
		string tag = option[ score::input_etables ];
		TR << "INPUT ETABLES " << tag << std::endl;
		for ( map<string,FArray3D<Real>*>::iterator i = etables.begin(), end = etables.end(); i != end; ++i ) {
			string ename = i->first;
			string fname = tag+"."+ename+".etable";
			std::ifstream input( fname.c_str() );
			input_etable(*(i->second),ename,input);
			input.close();
		}
	}

	if ( option[ score::output_etables ].user() ) {
		string header = option[ score::output_etables ];
		TR << "OUTPUT ETABLES " << header << std::endl;
		for ( map<string,FArray3D<Real>*>::iterator i = etables.begin(), end = etables.end(); i != end; ++i ) {
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
/// @brief Modify Potential: Modify Etable to better treat 0-0, C-C, and H-H interactions
///
/// @detailes
/// The Etables are  modified in three ways:
///    (1) the LK solvation energy is set to a constant below 4.2A to avoid shifting the position
///    of the minimum on the LJatr potential.  in refined decoys the peak in the C C pair
///    distribution function shifts to around 3.8 A from ~4.0A in native structures; the LJatr
///    potential has a minimum at 4.0A but this shifts towards smaller values because of the LK
///    solvation term, which became increasingly favorable at shorter distances
///    (2) the backbone carbonyl oxygen-carbonyl oxygen LJrep term has been modified to become
///    moderately repulsive at distances less than 3.6A.  this is to counteract the favorable
///    LJatr between the atoms (which have radii of ~1.4A and so a minimum at ~2.8A; very few
///    counts are observed in the pdb until around 3.2A) which leads to a significant shift in the
///    O O pair distribution function towards smaller values in refined decoys.  the repulsion is
///    a temporary proxy for the lack of explicit electrostatic repulsion in the current force
///    field.
///    (3) a third soft repulsion between non polar hydrogens that was also based on comparison of
///    refined decoy to native pdf's is currently commented out as the effects on packing have not
///    been tested.  it was observed that the protons tend to pile up at just beyond the point
///    where the repulsion becomes strong, perhaps due to a general tendency to overcontraction
///    because of long range LJatr interactions not compensated by interactions with solvent
///    (which are of course missing)
///
/// @author ctsa (last modified: 10-2003)
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::modify_pot()
{
	using namespace std;

	bool mod_hhrep = false;
	Real const h = 0.0;
	Real const c = 0.0;
	Real const w = 0.0;
	Real const e = 0.0;

	if ( mod_hhrep ) {
		TR << "fullatom_setup: modifying h-h repulsion "
			<< " hhrep center "   << c
			<< " hhrep height "   << h
			<< " hhrep width "    << w
			<< " hhrep exponent " << e
			<< std::endl;
	} {

		Real const bin = ( 4.2 * 4.2 / .05 ) + 1.0;
		Real dis;
		int const ibin( static_cast< int >( bin ) );
		chemical::AtomTypeSetCOP atom_set_ac( atom_set() );
		for ( int k = 1; k <= etable_disbins(); ++k ) {
			dis = std::sqrt( ( k - 1 ) * .05f );

			utility::vector1<int> carbon_types;
			carbon_types.push_back( atom_set_ac->atom_type_index("CH1") );
			carbon_types.push_back( atom_set_ac->atom_type_index("CH2") );
			carbon_types.push_back( atom_set_ac->atom_type_index("CH3") );
			carbon_types.push_back( atom_set_ac->atom_type_index("aroC") );

			if ( dis >= 4.2 )  continue;
			
			for ( int i = 1, i_end = carbon_types.size(); i <= i_end; ++i ) {
				for ( int j = 1, j_end = carbon_types.size(); j <= j_end; ++j ) {
					
					int const ii = carbon_types[i];
					int const jj = carbon_types[j];
					
					solv1_(k,jj,ii) = solv1_(ibin,jj,ii);
					solv1_(k,ii,jj) = solv1_(ibin,ii,jj);
					solv2_(k,jj,ii) = solv2_(ibin,jj,ii);
					solv2_(k,ii,jj) = solv2_(ibin,ii,jj);
					dsolv2_(k,jj,ii) = 0.0;
					dsolv2_(k,ii,jj) = 0.0;
					dsolv1_(k,ii,jj) = 0.0;
					dsolv1_(k,jj,ii) = 0.0;
					memb_solv1_(k,jj,ii) = memb_solv1_(ibin,jj,ii);
					memb_solv1_(k,ii,jj) = memb_solv1_(ibin,ii,jj);
					memb_solv2_(k,jj,ii) = memb_solv2_(ibin,jj,ii);
					memb_solv2_(k,ii,jj) = memb_solv2_(ibin,ii,jj);
					memb_dsolv2_(k,jj,ii) = 0.0;
					memb_dsolv2_(k,ii,jj) = 0.0;
					memb_dsolv1_(k,ii,jj) = 0.0;
					memb_dsolv1_(k,jj,ii) = 0.0;
				}
			}
		}
	}
}

void
MembEtable::smooth_etables()
{
	using namespace numeric::interpolation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	FArray1D<Real> dis(etable_disbins());
	for ( int i = 1; i <= etable_disbins(); ++i ) {
		dis(i) = sqrt( ( i - 1 ) * 1.0 / get_bins_per_A2() );
	}

	FArray1D<int> bin_of_dis(100);

	for ( int i = 1; i <= 100; ++i ) {
		float d = ((float)i)/10.0;
		bin_of_dis(i) = (int)( d*d * get_bins_per_A2() + 1 );
	}

	/////////////////////////////////////////////////////////////////////////////
	// 2) spline smooth solv1/solv2 at fa_max_dis cut && and first switchpoint
	/////////////////////////////////////////////////////////////////////////////

	using namespace numeric::interpolation::spline;
	TR << "smooth_etable: spline smoothing solvation etables (max_dis = " << max_dis() << ")" << std::endl;

	for ( int at1 = 1; at1 <= n_atomtypes(); ++at1 ) {
		for ( int at2 = 1; at2 <= n_atomtypes(); ++at2 ) {

			int SWTCH = 1;
			for ( SWTCH=1; SWTCH <= etable_disbins(); ++SWTCH ) {
				if ( (solv1_(SWTCH,at1,at2) != solv1_(1,at1,at2)) ||
						(solv2_(SWTCH,at1,at2) != solv2_(1,at1,at2)) ||
						(memb_solv1_(SWTCH,at1,at2) != memb_solv1_(1,at1,at2)) ||
						(memb_solv2_(SWTCH,at1,at2) != memb_solv2_(1,at1,at2)) ) {
					break;
				}
			}
			if ( SWTCH > etable_disbins() ) continue;

			int const S1 = std::max(1,SWTCH - 30);
			int const E1 = std::min(SWTCH + 20,406);
			int const S2 = bin_of_dis( (int)((max_dis()-1.5)*10.0) );
			int const E2 = etable_disbins();

			Real dsolv1e1 = (solv1_(E1+1,at1,at2)-solv1_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real dsolv1s2 = (solv1_(S2  ,at1,at2)-solv1_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));
			Real dsolv2e1 = (solv2_(E1+1,at1,at2)-solv2_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real dsolv2s2 = (solv2_(S2  ,at1,at2)-solv2_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));

			Real mdsolv1e1 = (memb_solv1_(E1+1,at1,at2)-memb_solv1_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real mdsolv1s2 = (memb_solv1_(S2  ,at1,at2)-memb_solv1_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));
			Real mdsolv2e1 = (memb_solv2_(E1+1,at1,at2)-memb_solv2_(E1  ,at1,at2))/(dis(E1+1)-dis(E1  ));
			Real mdsolv2s2 = (memb_solv2_(S2  ,at1,at2)-memb_solv2_(S2-1,at1,at2))/(dis(S2  )-dis(S2-1));

			SplineGenerator gen11( dis(S1), solv1_(S1,at1,at2), 0.0     , dis(E1), solv1_(E1,at1,at2), dsolv1e1 );
			SplineGenerator gen21( dis(S1), solv2_(S1,at1,at2), 0.0     , dis(E1), solv2_(E1,at1,at2), dsolv2e1 );
			SplineGenerator gen12( dis(S2), solv1_(S2,at1,at2), dsolv1s2, dis(E2), solv1_(E2,at1,at2), 0.0      );
			SplineGenerator gen22( dis(S2), solv2_(S2,at1,at2), dsolv2s2, dis(E2), solv2_(E2,at1,at2), 0.0      );

			SplineGenerator mgen11( dis(S1), memb_solv1_(S1,at1,at2), 0.0     , dis(E1), memb_solv1_(E1,at1,at2), mdsolv1e1 );
			SplineGenerator mgen21( dis(S1), memb_solv2_(S1,at1,at2), 0.0     , dis(E1), memb_solv2_(E1,at1,at2), mdsolv2e1 );
			SplineGenerator mgen12( dis(S2), memb_solv1_(S2,at1,at2), mdsolv1s2, dis(E2), memb_solv1_(E2,at1,at2), 0.0      );
			SplineGenerator mgen22( dis(S2), memb_solv2_(S2,at1,at2), mdsolv2s2, dis(E2), memb_solv2_(E2,at1,at2), 0.0      );

			InterpolatorOP interp11( gen11.get_interpolator() );
			InterpolatorOP interp21( gen21.get_interpolator() );

			InterpolatorOP minterp11( mgen11.get_interpolator() );
			InterpolatorOP minterp21( mgen21.get_interpolator() );

			for ( int i = S1; i <= E1; ++i ) {
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

			for ( int i = S2; i <= E2; ++i ) {
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
/// @brief Output Etable: Output etable data in a file in the same format used in input
/// @details File first line is <etable> <etable_disbins()> and other lines are <atom type 1>
/// <atomtype 1> <eval bin 1> <eval bin 2>...
///
/// @author sheffler (last modified: mar 19 2006)
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::output_etable(
	ObjexxFCL::FArray3D<Real> & etable,
	std::string label,
	std::ostream & out
) {

	using namespace std;
	using namespace ObjexxFCL;

	out << label << " " << etable_disbins() << endl;
	for ( int at1 = 1; at1 <= n_atomtypes(); at1++ ) {
		for ( int at2 = 1; at2 <= n_atomtypes(); at2++ ) {
			out << at1 << " "
				<< at2 << " ";
			for ( int bin = 1; bin <= etable_disbins(); bin++ ) {
				float evalue = etable(bin,at1,at2);
				out << evalue << ' ';
			}
			out << endl;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @brief: Input Etable: read in etable from a datafile
/// @details file first line is <etable> <etable_disbins()> and other lines are
/// <atom type 1> <atomtype 1> <eval bin 1> <eval bin 2>...
///
/// @author sheffler (last modified: mar 19 2006)
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
	if ( !(intmp >> lblin >> numdisbin) ) {
		TR << "input_etable: WARNING bad etable header " << buf << endl;
		return;
	}
	if ( lblin != label || etable_disbins() != numdisbin ) {
		TR << "input_etable: WARNING etable types don't match! "<< endl;
		TR << "              expected " << label << "," << etable_disbins()
			<< " got " << lblin << ',' << numdisbin<< endl;
		utility_exit_with_message( "input_etable: WARNING etable types don't match! " );
	} else {
		TR << "input_etable expected etable " << label << " of size " << etable_disbins()
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
				etable(bin,at1,at2) = evalue;
				etable(bin,at2,at1) = evalue;
			}
			scount++;
		}
	}
	TR << "              read " << scount << " of " << count << " lines" << endl;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief precalculate non-distance dependent coefficients of energy functions
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
/// @author ctsa (last modified: 10-2003)
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::precalc_etable_coefficients(
	FArray2< Real > & lj_sigma,
	FArray1< Real > & lk_inv_lambda2,
	FArray2< Real > & lk_coeff,
	FArray2< Real > & memb_lk_coeff,
	FArray2< Real > & lk_min_dis2sigma_value,
	FArray2< Real > & memb_lk_min_dis2sigma_value
)
{

	// locals
	Real sigma; //sigma6,sigma12,wdepth;
	Real inv_lambda;
	FArray1D< Real > lk_coeff_tmp( n_atomtypes() );
	FArray1D< Real > memb_lk_coeff_tmp( n_atomtypes() );
	Real thresh_dis,inv_thresh_dis2,x_thresh;
	Real const inv_neg2_tms_pi_sqrt_pi = { -0.089793561062583294 };
	// coefficient for lk solvation

	// include follows locals so that data statements can initialize included arrays
	for ( int i = 1, e = n_atomtypes(); i <= e; ++i ) {
		inv_lambda = 1.0/lk_lambda(i);
		lk_inv_lambda2(i) = inv_lambda * inv_lambda;
		lk_coeff_tmp(i) = inv_neg2_tms_pi_sqrt_pi * lk_dgfree(i) * inv_lambda;
		memb_lk_coeff_tmp(i) = inv_neg2_tms_pi_sqrt_pi * memb_lk_dgfree(i) * inv_lambda;
	}

	for ( int i = 1, e = n_atomtypes(); i <= e; ++i ) {
		for ( int j = i; j <= e; ++j ) {

			sigma = Wradius() * ( lj_radius(i) + lj_radius(j) );
			sigma = ( sigma < 1.0e-9 ? 1.0e-9 : sigma );

			// (bk) modify sigma for hbond donors and acceptors
			//   ctsa -- should these scale down by Wradius as well?

			// pb specific sigma correction for pairs between charged oxygen acceptors (15)
			// pb and hydroxyl oxygen donors (13). sigma correction for all polar H and charged oxygen
			// pb acceptor. Combinations of these corrections allow better prediction of both
			// pb hydroxyl O donor/charged O acceptor and charged NH donor/charged O acceptor
			// pb distances.

			if ( lj_use_hbond_radii() ) {
				if ( ( atom_type(i).is_acceptor() && atom_type(j).is_donor() ) ||
						( atom_type(i).is_donor() && atom_type(j).is_acceptor() ) ) {
					sigma = get_lj_hbond_dis();
				} else if ( ( atom_type(i).is_acceptor() && atom_type(j).is_polar_hydrogen() ) ||
						( atom_type(i).is_polar_hydrogen() && atom_type(j).is_acceptor() ) ) {
					sigma = get_lj_hbond_hdis();

				}
			}

			//lin   modify sigma for water and hbond donors/acceptors
			if ( lj_use_water_radii() ) {
				if ( ( ( atom_type(i).is_acceptor() ||
						atom_type(i).is_donor() ) &&
						atom_type(j).is_h2o() ) ||
						( ( atom_type(j).is_acceptor() ||
						atom_type(j).is_donor() ) &&
						atom_type(i).is_h2o() ) ) {
					sigma = lj_water_dis();
				} else if ( ( atom_type(i).is_polar_hydrogen() &&
						atom_type(j).is_h2o() ) ||
						( atom_type(j).is_polar_hydrogen() &&
						atom_type(i).is_h2o() ) ) {
					sigma = lj_water_hdis();
				}
			}

			lj_sigma(i,j) = sigma;
			lj_sigma(j,i) = lj_sigma(i,j);

			// ctsa - precalculated lk solvation coefficients
			lk_coeff(i,j) = lk_coeff_tmp(i) * lk_volume(j);
			lk_coeff(j,i) = lk_coeff_tmp(j) * lk_volume(i);
			memb_lk_coeff(i,j) = memb_lk_coeff_tmp(i) * lk_volume(j);
			memb_lk_coeff(j,i) = memb_lk_coeff_tmp(j) * lk_volume(i);

			// ctsa - when dis/sigma drops below lk_min_dis2sigma,
			//   a constant lk solvation value equal to the value at the
			//   switchover point is used. That switchover-point value
			//   is calculated here and stored in lk_min_dis2sigma_value
			thresh_dis = lk_min_dis2sigma()*sigma;
			inv_thresh_dis2 = 1./( thresh_dis * thresh_dis );
			Real dis_rad = thresh_dis - lj_radius(i);
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(i);
			lk_min_dis2sigma_value(i,j) = std::exp(-x_thresh) * lk_coeff(i,j) *
				inv_thresh_dis2;
			memb_lk_min_dis2sigma_value(i,j) = std::exp(-x_thresh) * memb_lk_coeff(i,j) *
				inv_thresh_dis2;

			dis_rad = thresh_dis - lj_radius(j);
			x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2(j);
			lk_min_dis2sigma_value(j,i) = std::exp(-x_thresh) * lk_coeff(j,i) *
				inv_thresh_dis2;
			memb_lk_min_dis2sigma_value(j,i) = std::exp(-x_thresh) * memb_lk_coeff(j,i) *
				inv_thresh_dis2;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @brief calc all etable values given a distance and atom-type pair
///
/// @details
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
/// @author ctsa (last modified: 10-2003)
/////////////////////////////////////////////////////////////////////////////////
void
MembEtable::calc_etable_value(
	Real & dis2,
	int & atype1,
	int & atype2,
	Real & solvE1,
	Real & solvE2,
	Real & dsolvE1,
	Real & dsolvE2,
	FArray2< Real > & lj_sigma,
	FArray1< Real > & lk_inv_lambda2,
	FArray2< Real > & lk_coeff,
	FArray2< Real > & lk_min_dis2sigma_value,
	Real & memb_solvE1,
	Real & memb_solvE2,
	FArray2< Real > & memb_lk_coeff,
	FArray2< Real > & memb_lk_min_dis2sigma_value,
	Real & memb_dsolvE1,
	Real & memb_dsolvE2
) {
	// locals
	Real x1,x2;
	Real dis;
	Real inv_dis,inv_dis2;
	Real dis2sigma;
	int xtra_atype1,xtra_atype2;

	// include after local variables to allow data statements to initialize
	solvE1 = 0.;
	solvE2 = 0.;
	dsolvE1 = 0.;
	dsolvE2 = 0.;
	memb_solvE1 = 0.;
	memb_solvE2 = 0.;
	memb_dsolvE1 = 0.;
	memb_dsolvE2 = 0.;

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2() + epsilon() ) return;
	if ( dis2 < min_dis2() ) dis2 = min_dis2();

	dis = std::sqrt(dis2);
	inv_dis = 1.0/dis;
	inv_dis2 = inv_dis * inv_dis;

	xtra_atype1 = atype1;
	xtra_atype2 = atype2;

	dis2sigma = dis / lj_sigma(xtra_atype1,xtra_atype2);

	if ( dis2sigma < lk_min_dis2sigma() ) {
		solvE1 = lk_min_dis2sigma_value(xtra_atype1,xtra_atype2);
		solvE2 = lk_min_dis2sigma_value(xtra_atype2,xtra_atype1);
		dsolvE1 = 0.0;
		dsolvE2 = 0.0;
		memb_solvE1 = memb_lk_min_dis2sigma_value(xtra_atype1,xtra_atype2);
		memb_solvE2 = memb_lk_min_dis2sigma_value(xtra_atype2,xtra_atype1);
		memb_dsolvE1 = memb_dsolvE2 = 0.0;

	} else {

		Real dis_rad = dis - lj_radius(atype1);
		x1 = ( dis_rad * dis_rad ) * lk_inv_lambda2(atype1);
		dis_rad = dis - lj_radius(atype2);
		x2 = ( dis_rad * dis_rad ) * lk_inv_lambda2(atype2);

		solvE1 = std::exp(-x1) * lk_coeff(atype1,atype2) * inv_dis2;
		solvE2 = std::exp(-x2) * lk_coeff(atype2,atype1) * inv_dis2;
		memb_solvE1 = std::exp(-x1) * memb_lk_coeff(atype1,atype2) * inv_dis2;
		memb_solvE2 = std::exp(-x2) * memb_lk_coeff(atype2,atype1) * inv_dis2;

		dsolvE1 = -2.0 * solvE1 *
			(((dis-lj_radius(atype1))*lk_inv_lambda2(atype1))+inv_dis);
		dsolvE2 = -2.0 * solvE2 *
			(((dis-lj_radius(atype2))*lk_inv_lambda2(atype2))+inv_dis);

		memb_dsolvE1 = -2.0 * memb_solvE1 *
			(((dis-lj_radius(atype1))*lk_inv_lambda2(atype1))+inv_dis);
		memb_dsolvE2 = -2.0 * memb_solvE2 *
			(((dis-lj_radius(atype2))*lk_inv_lambda2(atype2))+inv_dis);
	}
}

} // etable
} // scoring
} // core
