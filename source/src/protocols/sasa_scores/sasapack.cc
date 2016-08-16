// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/phil/test1.cc
/// @brief  Some simple examples of how to use basic functionality + some DNA functionality
/// @author Phil Bradley (pbradley@fhcrc.org)

// libRosetta headers
#include <protocols/sasa_scores/sasapack.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/VariantType.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sasa_scores.sasapack" );


namespace protocols {
namespace sasa_scores {

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::chemical;
using namespace core::scoring;
using utility::vector1;
//using namespace basic::options;

typedef utility::vector1< Real > Reals;

using namespace std;
using namespace ObjexxFCL::format;

class Poly {
public:
	Real
	operator()( Real const x ) const;

	Size
	degree() const { return coeffs_.size() -1; }

	friend
	std::istream & operator >>( std::istream & is, Poly & p );

private:
	Reals coeffs_;


};

///////
Real
Poly::operator()( Real const x ) const
{
	runtime_assert( !coeffs_.empty() );

	Size const degree( coeffs_.size() - 1 );

	Real polyval(0.0);
	Real x_raised( 1.0 );
	for ( Size i=0; i<= degree; ++i ) {
		polyval += coeffs_[ degree+1-i ] * x_raised;
		x_raised *= x;
	}
	return polyval;
}

std::istream & operator >>( std::istream & is, Poly & p )
{
	Size degree(0);
	string tmp1,tmp2,tmp3;
	is >> tmp1 >> tmp2 >> degree >> tmp3;
	if ( is.fail() || tmp1 != "POLY" || tmp2 != "DEGREE" || tmp3 != "COEFFS_HI2LO" || degree > 10000 ) {
		is.setstate( std::ios_base::failbit );
		p.coeffs_.clear();
		return is;
	}
	p.coeffs_.resize( degree+1 );
	for ( Size i=1; i<=degree+1; ++i ) {
		is >> p.coeffs_[i];
	}
	if ( is.fail() ) {
		p.coeffs_.clear();
		return is;
	}
	return is;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PPoly {
public:
	Real
	operator()( Real const x ) const;

	friend
	std::istream & operator >>( std::istream & is, PPoly & p );

	Size
	max_degree() const;

private:
	vector1< Poly > polys_;
	Reals xmins_;
	Reals xmaxs_;

};

Real
PPoly::operator()( Real const x ) const
{
	for ( Size ii=1; ii<= xmins_.size(); ++ii ) {
		if ( x >= xmins_[ii] && x <= xmaxs_[ii] ) return polys_[ii](x);
	}
	utility_exit_with_message("PPoly::operator() x out of range: "+ObjexxFCL::string_of(x)); // probably a nicer way to handle this...
	return 0.0;
}

Size
PPoly::max_degree() const
{
	Size md(0);
	for ( Size ii=1; ii<= polys_.size(); ++ii ) {
		md = max( md, polys_[ii].degree() );
	}
	return md;
}

std::istream & operator >>( std::istream & is, PPoly & pp )
{
	string tmp1,tmp2;
	Size npoly;
	is >> tmp1 >> tmp2 >> npoly;
	if ( is.fail() || tmp1 != "PPOLY" || tmp2 != "NPOLY" ) {
		is.setstate( std::ios_base::failbit );
		pp.polys_.clear();
		return is;
	}
	for ( Size i=1; i<= npoly; ++i ) {
		Poly p;
		Real xmin,xmax;
		is >> tmp1 >> xmin >> tmp2 >> xmax >> p;
		if ( is.fail() || tmp1 != "XMIN" || tmp2 != "XMAX" ) break;
		pp.polys_.push_back(p);
		pp.xmins_.push_back( xmin );
		pp.xmaxs_.push_back( xmax );
	}
	if ( is.fail() || pp.polys_.size() != npoly ) {
		is.setstate( std::ios_base::failbit );
		pp.polys_.clear(); pp.xmins_.clear(); pp.xmaxs_.clear();
		return is;
	}
	return is;
}

void
help_load_data(
	utility::vector1< PPoly > & polys,
	Reals & avg_sasa14s,
	std::string datafile,
	std::string tag1,
	std::string tag2,
	std::string function_name
);

///////////////////////////////////////////////////////////////////////////////
void
load_sasapack_polynomial_coefficients(
	vector1< PPoly > & polys,
	Reals & avg_sasa14s
) {
	string const datafile( "scoring/sasa_scores/sasapack_datafile_v1.txt" );// TO DO: make this an option
	help_load_data( polys, avg_sasa14s, datafile, "PPOLY_SASAPACK", "PPOLYVAL_SASAPACK", "load_sasapack_polynomial_coefficients" );
}


///////////////////////////////////////////////////////////////////////////////
void
load_avge_polynomial_coefficients(
	vector1< PPoly > & polys,
	Reals & avg_sasa14s
) {
	string const datafile( "scoring/sasa_scores/avge_datafile_score12prime_v1.txt" );// TO DO: make this an option
	help_load_data( polys, avg_sasa14s, datafile, "PPOLY_NORME", "PPOLYVAL_NORME", "load_avge_polynomial_coefficients" );
}

///////////////////////////////////////////////////////////////////////////////
/// helper function to open a file and respond to appropriate flags
void
help_load_data(
	vector1< PPoly > & polys,
	Reals & avg_sasa14s,
	string datafile,
	string tag1,
	string tag2,
	string function_name
) {
	polys.clear();
	polys.resize( num_canonical_aas );
	avg_sasa14s.clear();
	avg_sasa14s.resize( num_canonical_aas );

	utility::io::izstream data;
	//string const datafile( "scoring/sasa_scores/avge_datafile_score12prime_v1.txt" );// TO DO: make this an option
	basic::database::open( data, datafile );

	map< std::pair< AA, Size >, vector1< std::pair< Real, Real > > > all_polyvals;

	string line;
	bool found_avg_sasa( false );
	while ( getline( data,line ) ) {
		string linetag, name1;
		Size degree;
		istringstream is(line );
		is >> linetag;
		if ( linetag == tag1 ) {
			PPoly p;
			is >> name1 >> degree >> p;
			runtime_assert( !is.fail() );
			AA const aa( aa_from_oneletter_code( name1[0] ) );
			runtime_assert( aa >= aa_ala && aa <= num_canonical_aas );
			polys[ aa ] = p;
		} else if ( linetag == tag2 ) {
			Real x,y;
			is >> name1 >> degree >> x >> y;
			std::pair< AA, Size > const p( make_pair( aa_from_oneletter_code( name1[0] ), degree ) );
			all_polyvals[p].push_back( make_pair(x,y) );
		} else if ( linetag == "AVG_SASA" ) {
			Real mean;
			is >> name1 >> linetag >> mean;
			found_avg_sasa = true;
			AA const aa( aa_from_oneletter_code( name1[0] ) );
			runtime_assert( aa >= aa_ala && aa <= num_canonical_aas );
			avg_sasa14s[ aa ] = mean;
		}
	}
	runtime_assert( found_avg_sasa );

	for ( Size i=1; i<= 20; ++i ) {
		AA const aa = AA(i);
		PPoly const & p( polys[i] );
		Size const degree( p.max_degree() );
		cout << function_name << ": using degree " << degree << " polynomial for " << aa <<
			" datafile: " << datafile << endl; /// NOTE: cout
		vector1< std::pair< Real, Real > > const polyvals( all_polyvals[ make_pair( aa, degree ) ] );
		Size count(0);
		Real err( 0.0 );
		for ( Size ii=1; ii<= polyvals.size(); ++ii ) {
			Real const x( polyvals[ii].first ), expected_y( polyvals[ii].second ), recomputed_y( p( x ) );
			err += ( expected_y - recomputed_y ) * ( expected_y - recomputed_y );
			++count;
		}
		if ( count ) err /= count;
		TR.Trace << "polyval_err: " << aa << ' ' << degree << " npoints: " << I(4,count) << " err-per-point: " <<
			F(9,3,err) << endl;
		runtime_assert( err < 1e-3 );
	}
	data.close();
}


///////////////////////////////////////////////////////////////////////////////
/// NOTE: this does not include the probe radius in the sasa value, hence somewhat specialized for sasapack
void
compute_residue_sasas_for_sasa_scores(
	Real const probe_radius,
	Pose const & pose,
	Reals & rsd_sasa
)
{
	bool const use_big_polar_H( false ), use_naccess_sasa_radii( false ), expand_polar_radii( false ),
		include_probe_radius_in_atom_radii( false ), use_lj_radii( true ); // NOTE -- LJ RADII
	Real const polar_expansion_radius_unused( 1.0 );

	id::AtomID_Map< Real > atom_sasa;
	id::AtomID_Map< bool > atom_subset;
	atom_subset.clear(); // unnecessary?
	core::pose::initialize_atomid_map( atom_subset, pose, true );

	calc_per_atom_sasa( pose,   atom_sasa,   rsd_sasa, probe_radius, use_big_polar_H,
		atom_subset, use_naccess_sasa_radii, expand_polar_radii, polar_expansion_radius_unused,
		include_probe_radius_in_atom_radii, use_lj_radii );

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// NOTE: for this to make any sense it's critical that the same scorefxn be used in scoring decoys and natives
///
void
compute_avge_scores(
	Pose const & pose_in,
	Reals & residue_avge,
	Reals & residue_normsasa,
	Real & average_avge, // over all residues that were counted
	Real & average_normsasa
)
{
	bool const ignore_gly_paa( true ), ignore_pro_close( true ), ignore_omega( true ), ignore_fa_dun( true ),
		ignore_fa_rep( true );

	static ScoreFunctionOP fa_scorefxn( 0 );

	if ( !fa_scorefxn ) {
		fa_scorefxn = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
		if ( pose::symmetry::is_symmetric( *fa_scorefxn ) ) {
			/// seems like residue energies may be messed up in the symmetric case, for intra-monomer interactions
			///
			fa_scorefxn = scoring::symmetry::asymmetrize_scorefunction( *fa_scorefxn );
			runtime_assert( !pose::symmetry::is_symmetric( *fa_scorefxn ) );
		}
	}


	residue_avge.clear(); residue_avge.resize( pose_in.total_residue(), 0.0 );
	residue_normsasa.clear(); residue_normsasa.resize( pose_in.total_residue(), 0.0 );

	static vector1< PPoly > polys;
	static Reals avg_sasa14s;
	if ( polys.empty() ) {
		load_avge_polynomial_coefficients( polys, avg_sasa14s );
	}

	Pose pose( pose_in );

	if ( pose::symmetry::is_symmetric( pose ) ) pose::symmetry::make_asymmetric_pose( pose );

	/// need to compute sasa
	Real const probe_radius( 1.4 );
	Reals rsd_sasa;
	compute_residue_sasas_for_sasa_scores( probe_radius, pose, rsd_sasa );

	(*fa_scorefxn)( pose ); // rescore with the same energy function used in fitting coeffs

	EnergyMap const & wts( fa_scorefxn->weights() );

	average_avge = average_normsasa = 0.0;

	Size count(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_protein() ) continue;
		if ( rsd.is_lower_terminus() || rsd.is_upper_terminus() ) continue;
		if ( rsd.has_variant_type( chemical::DISULFIDE ) ) continue; // avge ==> 0.0 for disulfs
		EnergyMap const & rsd_energies( pose.energies().residue_total_energies(i) );
		Real const total_energy( rsd_energies.dot( fa_scorefxn->weights() ) );

		Real normE( total_energy );
		if ( ignore_gly_paa   ) normE -= rsd_energies[ p_aa_pp   ] * wts[ p_aa_pp   ];
		if ( ignore_pro_close ) normE -= rsd_energies[ pro_close ] * wts[ pro_close ];
		if ( ignore_omega     ) normE -= rsd_energies[ omega     ] * wts[ omega     ];
		if ( ignore_fa_dun    ) normE -= rsd_energies[ fa_dun    ] * wts[ fa_dun    ];
		if ( ignore_fa_rep    ) normE -= rsd_energies[ fa_rep    ] * wts[ fa_rep    ];

		Real const sasa14( rsd_sasa[i] ), expected_normE( polys[ rsd.aa() ]( sasa14 ) );

		residue_avge    [i] = normE - expected_normE;
		residue_normsasa[i] = sasa14 - avg_sasa14s[ rsd.aa() ];

		average_avge     += residue_avge    [ i ];
		average_normsasa += residue_normsasa[ i ];
		++count;
	}

	if ( count ) {
		average_avge     /= count;
		average_normsasa /= count;
	}
}


///////////////////////////////////////////////////////////////////////////////
void
compute_sasapack_scores(
	Pose const & pose,
	Reals & residue_sasapack,
	Reals & residue_normsasa,
	Real & average_sasapack,
	Real & average_normsasa
)
{
	static vector1< PPoly > polys;
	static Reals avg_sasa14s;
	if ( polys.empty() ) {
		load_sasapack_polynomial_coefficients( polys, avg_sasa14s );
	}

	//// need to compute sasas
	Real const big_probe_radius( 1.4 ), small_probe_radius( 0.5 );
	Reals rsd_sasa_big_probe, rsd_sasa_small_probe;
	compute_residue_sasas_for_sasa_scores(   big_probe_radius, pose, rsd_sasa_big_probe );
	compute_residue_sasas_for_sasa_scores( small_probe_radius, pose, rsd_sasa_small_probe );

	residue_sasapack.clear(); residue_sasapack.resize( pose.total_residue(), 0. );
	residue_normsasa.clear(); residue_normsasa.resize( pose.total_residue(), 0. );

	average_sasapack = average_normsasa = 0.0;

	Size count(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_protein() ) continue;
		if ( rsd.is_lower_terminus() || rsd.is_upper_terminus() ) continue;
		if ( rsd.has_variant_type( chemical::DISULFIDE ) ) continue; // sasapack ==> 0.0 for disulfs

		Real const actual_sasa14( rsd_sasa_big_probe[i] ), actual_sasa5( rsd_sasa_small_probe[i] );
		Real const expected_sasa5( polys[ rsd.aa() ]( actual_sasa14 ) );
		residue_sasapack[ i ] = actual_sasa5 - expected_sasa5;
		residue_normsasa[ i ] = actual_sasa14 - avg_sasa14s[ rsd.aa() ];

		average_sasapack += residue_sasapack[i];
		average_normsasa += residue_normsasa[i];
		++count;
	}
	if ( count ) {
		average_sasapack /= count;
		average_normsasa /= count;
	}
}


} // ns sasa_scores
} // ns protocols
