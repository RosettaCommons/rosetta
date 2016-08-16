// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/dna/DNA_DihedralPotential.cc
/// @brief  dna scoring
/// @author Phil Bradley

#include <core/scoring/dna/DNA_DihedralPotential.hh>
//#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/base_geometry.hh>
//#include <core/scoring/dna/DNA_BasePotential.hh>

#include <core/scoring/ScoringManager.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/izstream.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/format.hh>

//using namespace ObjexxFCL;

namespace core {
namespace scoring {
namespace dna {

static basic::Tracer TR("core.scoring.dna.DNA_DihedralPotential" );
using utility::vector1;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;

// Size const
// DNA_DihedralPotential::chi_torsion_index_( 7 );

/// iname1 goes from 1-4 : acgt
inline
Size
get_iname1_from_name1( char const name1 )
{
	static std::string const bases( "acgt" );
	if ( bases.find( name1 ) == std::string::npos ) {
		utility_exit_with_message("unrecognized name1: "+std::string(1,name1) );
	}
	return bases.find( name1 ) + 1;
}

DNA_DihedralPotential::DNA_DihedralPotential( std::string const & filename )
{
	read_dna_geometry_log_file( filename );
}


DNA_DihedralPotential::DNA_DihedralPotential()
{
	read_dna_geometry_log_file_from_database( "scoring/dna/bound_dna_dihedrals.txt" );
	//read_dna_geometry_log_file( "bound_dna.stats" );
}

// [0,360)
void
put_angle_in_0to360( Real & angle )
{
	while ( angle <0 ) angle += 360.0;
	while ( angle >= 360.0 ) angle -= 360.0;
}

Real
get_triple_bin( Real angle )
{
	put_angle_in_0to360( angle );
	if ( angle < 120.0 ) return 1;
	else if ( angle < 240.0 ) return 2;
	else return 3;
}

Real
get_b1b2_bin( Real epsilon, Real zeta )
{
	Real const dev( basic::subtract_degree_angles( epsilon, zeta ) );
	if ( dev < 0 ) return 1;
	else return 2;
}


void
get_mean_median_and_sdev(
	vector1< Real > vals,
	Real & mean,
	Real & median,
	Real & sdev
)
{
	mean = median = sdev = 0.0;

	if ( vals.empty() ) return;
	std::sort( vals.begin(), vals.end() );
	median = vals[ (vals.size()+1)/2 ];

	for ( vector1< Real >::const_iterator val= vals.begin(); val != vals.end(); ++val ) mean += *val;
	mean /= vals.size();
	for ( vector1< Real >::const_iterator val= vals.begin(); val != vals.end(); ++val ) {
		sdev += ( *val - mean ) * ( *val - mean );
	}
	sdev = std::sqrt( sdev / vals.size() );
}


bool
DNA_DihedralPotential::skip_torsion( conformation::Residue const & rsd, Size const tor ) const
{
	return ( ( !rsd.is_DNA() ) ||
		( rsd.is_lower_terminus() && tor == 1 ) ||
		( rsd.is_upper_terminus() && ( tor == 5 || tor == 6 ) ) );
}




/// 1- delta
/// 2- "chi2"
/// 3- "chi3"
/// 4- "chi4"
///
void
DNA_DihedralPotential::eval_sugar_torsion_score_and_deriv(
	Real const torsion,
	Size const tor,
	conformation::Residue const & rsd,
	Size const pucker,
	Real & score,
	Real & dscore_dtor
) const
{
	static Real const sdev(basic::options::option[basic::options::OptionKeys::dna::specificity::dna_sugar_torsion_sdev]);
	assert( pucker<= 9 );
	assert( tor >=1 && tor<= 4 );

	Size const iname( get_iname1_from_name1( rsd.name1() ) );

	Real const mean( mean_sugar_torsion_[ iname ][ pucker ][ tor ] );
	Real const dev( basic::subtract_degree_angles( torsion, mean ) );

	score = numeric::square( dev / sdev );
	dscore_dtor = ( 2* dev ) / (sdev*sdev);

}

void
DNA_DihedralPotential::get_sugar_torsion_mean_and_sdev(
	Size const tor,
	conformation::Residue const & rsd,
	Size const pucker,
	Real & mean,
	Real & sdev
) const
{
	static Real const sdev_(basic::options::option[basic::options::OptionKeys::dna::specificity::dna_sugar_torsion_sdev]);

	assert( pucker<= 9 );
	assert( tor >=1 && tor<= 4 );

	Size const iname( get_iname1_from_name1( rsd.name1() ) );

	mean = mean_sugar_torsion_[ iname ][ pucker ][ tor ];
	sdev = sdev_;
}

void
DNA_DihedralPotential::eval_harmonic_backbone_torsion_score_and_deriv(
	Size const tor,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real & score,
	Real & dscore_dtor
) const
{
	using namespace basic::options;
	score = dscore_dtor = 0.0;
	if ( tor == 4 ) return;

	static utility::vector1< Real > const sdev_backbone_torsion
		( option[ OptionKeys::dna::specificity::dna_backbone_torsion_sdevs ].user() ?
		vector1< Real >( option[ OptionKeys::dna::specificity::dna_backbone_torsion_sdevs ] ) :
		vector1< Real >( utility::tools::make_vector1( 17.0, 30.0, 15.0, 0.0, 20.0, 30.0 ) ) );

	runtime_assert( sdev_backbone_torsion.size() == 6 ); // in case from cmdline

	Size const seqpos( rsd.seqpos() );
	assert( tor <= 6 );
	Real const angle( rsd.mainchain_torsion( tor ) );

	Real mean;

	if ( tor == 1 || tor == 3 ) {
		// alpha and gamma
		Size const bin( get_triple_bin( angle ) );
		mean = mean_backbone_torsion_[ tor ][ bin ];
	} else if ( tor == 2 ) {
		if ( rsd.is_lower_terminus() || rsd.seqpos() == 1 || !pose.residue( rsd.seqpos() - 1 ).is_DNA() ) {
			mean = mean_backbone_torsion_[ tor ][ 1 ]; // use the B1 value
		} else {
			Real const epsilon( pose.residue( seqpos-1 ).mainchain_torsion( 5 ) ),
				zeta( pose.residue( seqpos-1 ).mainchain_torsion( 6 ) );
			mean = mean_backbone_torsion_[ tor ][ get_b1b2_bin( epsilon, zeta ) ];
		}
	} else {
		assert( tor == 5 || tor == 6 );
		Real const epsilon( rsd.mainchain_torsion( 5 ) ), zeta( rsd.mainchain_torsion( 6 ) );
		mean = mean_backbone_torsion_[ tor ][ get_b1b2_bin( epsilon, zeta ) ];
	}

	Real const dev( basic::subtract_degree_angles( angle, mean ) );
	Real const sdev( sdev_backbone_torsion[ tor ] );

	score = numeric::square( dev / sdev );
	dscore_dtor = ( 2* dev ) / (sdev*sdev);
}

void
DNA_DihedralPotential::eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv
(
	conformation::Residue const & rsd,
	pose::Pose const &,// pose,
	Size const pucker,
	Real & score,
	Real & dscore_dchi
) const
{
	static Real const sdev(basic::options::option[basic::options::OptionKeys::dna::specificity::dna_chi_torsion_sdev]);
	Real angle( rsd.chi(1) );
	put_angle_in_0to360( angle );

	Real mean(0.0);
	if ( angle < 100.0 && angle > 20.0 ) {
		// syn !!
		mean = 50.0;
	} else {
		mean = mean_sugar_torsion_[ get_iname1_from_name1( rsd.name1() ) ][ pucker ][ 5 ];
	}

	Real const dev( basic::subtract_degree_angles( rsd.chi(1), mean ) );

	score = numeric::square( dev / sdev );
	dscore_dchi = ( 2* dev ) / (sdev*sdev);

}

Real
DNA_DihedralPotential::get_mean_sugar_pucker_dependent_chi( conformation::Residue const & rsd ) const
{
	std::pair< std::string, int > pucker;
	get_sugar_pucker( rsd, pucker );
	Size const pucker_index( pucker.second );

	return mean_sugar_torsion_[ get_iname1_from_name1( rsd.name1() ) ][ pucker_index ][ 5 ];
}


void
filter_torsions_by_iname(
	Size const iname,
	vector1< Size > const & torsion_inames,
	vector1< Real > & torsions
)
{
	assert( torsions.size() == torsion_inames.size() );

	for ( Size i= torsions.size(); i>= 1; --i ) {
		if ( torsion_inames[i] != iname ) torsions.erase( torsions.begin() + i-1 );
	}

}



/// should only do this once
void
DNA_DihedralPotential::parse_dna_geometry_log(
	std::istream & data
)
{
	// now we are also going to get means and sdevs for the sugar torsions
	//
	// indexed by sugar_torsions[ iname ][ pucker-integer ][ which-angle ]
	//
	// iname:             1-4
	// pucker-integer:    0-9
	// which-angle:       1-5
	//
	// where which-angle =
	// 1 - delta
	// 2 - chi2
	// 3 - chi3
	// 4 - chi4
	// 5 - chi1 (not really a sugar torsion!)
	//

	utility::vector1< utility::vector1< Real > > all_torsions_; // tmp legacy underscore
	//utility::vector1< Size > all_iname1_;
	utility::vector1< bool > is_lower_terminus_;
	utility::vector1< bool > is_upper_terminus_;

	utility::vector0< utility::vector1< utility::vector1< Real > > > sugar_torsions;
	utility::vector0< utility::vector1< Size > > sugar_inames;
	sugar_torsions.resize( 10 );
	sugar_inames.resize( 10 );
	for ( Size p=0; p<10; ++p ) {
		sugar_torsions[p].resize( 5 ); // delta chi2 chi3 chi4 chi
	}

	// open the file from the database
	utility::vector1< Real > tor( 7 ), params( 6 );
	std::string line, tag;
	utility::vector1< std::string > tags( 6 );
	while ( getline( data, line ) ) {
		std::istringstream l( line );
		l >> tag;
		if ( tag == "DNA_DIHEDRALS" ) {
			// DNA_DIHEDRALS   85 a C2* endo  6 -- -- B1  -9    0.0  178.6   53.6  143.5  178.8  274.5   60.0  122.5  339.5  153.9 dna_dihedral_energy:    15.806 1a1f.pdb
			// DNA_DIHEDRALS  299 g C2* endo  6 g- g+ B2   7  270.9  187.6   51.3  141.0  244.4  173.9  260.7  138.3  320.5  165.8 dna_dihedral_energy:    25.891 anii.pdb
			// DNA_DIHEDRALS  300 c C1* exxo  5 g- g+ --   0  327.4  131.9   21.3  133.6    0.0    0.0  251.9  142.0  322.5  162.5 dna_dihedral_energy:    16.872 anii.pdb
			char name1;
			std::string agtag1, agtag2, b12tag;
			Size pucker_int;
			Real chi2, chi3, chi4;
			l >> tag >> name1 >> tag >> tag >> pucker_int >> agtag1 >> agtag2 >> b12tag >> tag >>
				tor[1] >> tor[2] >> tor[3] >> tor[4] >> tor[5] >> tor[6] >> tor[7] >> chi2 >> chi3 >> chi4;
			if ( l.fail() || pucker_int>9  ) utility_exit_with_message("bad line! "+line );

			bool const is_lower_terminus( ( agtag1 == "--" ) );
			bool const is_upper_terminus( ( b12tag == "--" ) );

			if ( !is_lower_terminus && !is_upper_terminus ) {
				sugar_torsions[ pucker_int ][ 1 ].push_back( tor[4] );
				sugar_torsions[ pucker_int ][ 2 ].push_back( chi2 );
				sugar_torsions[ pucker_int ][ 3 ].push_back( chi3 );
				sugar_torsions[ pucker_int ][ 4 ].push_back( chi4 );
				sugar_torsions[ pucker_int ][ 5 ].push_back( tor[7] ); // chi, ie glycosidic dihedral
				sugar_inames[ pucker_int ].push_back( get_iname1_from_name1( name1 ) ); // for filtering
			}

			all_torsions_.push_back( tor );
			is_lower_terminus_.push_back( is_lower_terminus );
			is_upper_terminus_.push_back( is_upper_terminus );
			//all_iname1_.push_back( get_iname1_from_name1( name1 ) );
		}
	}

	{ // get means for sugar torsions and chi
		//
		Real const min_chi_torsion( 120.0 ); // dont include syn conformations
		Size const min_torsions( 8 );
		using basic::subtract_degree_angles;
		mean_sugar_torsion_.clear(); mean_sugar_torsion_.resize( 4 );
		for ( Size iname=1; iname<= 4; ++iname ) { // acgt
			mean_sugar_torsion_[ iname ].resize( 10 );
			for ( Size p=0; p<= 9; ++p ) {
				mean_sugar_torsion_[ iname ][ p ].resize( 5 ); // delta, "chi2", "chi3", "chi4", chi
				for ( Size d=1; d<= 5; ++d ) {
					bool const sequence_dependent( d == 5 ); // chi is sequence-dependent
					utility::vector1< Real > torsions( sugar_torsions[p][d] );
					if ( sequence_dependent ) filter_torsions_by_iname( iname, sugar_inames[p], torsions );

					if ( torsions.size() < min_torsions ) {
						int window(0);
						while ( torsions.size() < min_torsions ) {
							++window;
							for ( int offset= -window; offset<= window; ++offset ) {
								if ( abs( offset ) < window ) continue; // dont use torsions we already got last time
								if ( int(p) + offset < 0 || int(p) + offset>9 ) continue;
								utility::vector1< Real > other_torsions( sugar_torsions[p+offset][d] );
								if ( sequence_dependent ) filter_torsions_by_iname( iname, sugar_inames[p+offset], other_torsions );
								Size const needed( min_torsions - torsions.size() );
								for ( Size i=1; i<= needed && i<= other_torsions.size(); ++i ) {
									//TR.Trace << "steal torsion: " << p << ' ' << p+offset << ' ' << i << std::endl;
									torsions.push_back( other_torsions[i] );
								}
							}
						}
					}

					std::sort( torsions.begin(), torsions.end() );
					if ( d == 5 ) {
						// no torsions less than min_chi_torsion
						while ( torsions.size() > min_torsions-2 && torsions.front() < min_chi_torsion ) {
							TR.Trace << "skip-syn-chi: " << iname << ' ' << p << ' ' << d << ' ' << torsions.front() << std::endl;
							torsions.erase( torsions.begin() );
						}
						TR.Trace << "still-syn-chi? " << iname << ' ' << p << ' ' << d << ' ' << torsions.front() << std::endl;
					}
					Size const ntorsions( torsions.size() );
					Real const md( torsions[ ( ntorsions+1)/2 ] );
					Real mn( 0.0 );
					for ( Size i=1; i<= ntorsions; ++i ) mn += ( subtract_degree_angles( torsions[i], md ) );
					mn /= ntorsions;
					mn += md;
					Real sdev( 0.0 );
					for ( Size i=1; i<= ntorsions; ++i ) sdev += numeric::square( subtract_degree_angles( torsions[i], mn ) );
					sdev = std::sqrt( sdev / ntorsions );
					TR.Trace << "sugar_torsion: " << iname << ' ' << p << ' ' << d <<
						" ntor: " << ObjexxFCL::format::I(5,ntorsions) <<
						" median: " << ObjexxFCL::format::F(9,3,md) <<
						" mean: "   << ObjexxFCL::format::F(9,3,mn ) <<
						" sdev: "   << ObjexxFCL::format::F(9,3,sdev) << std::endl;
					mean_sugar_torsion_[ iname ][ p ][ d ] = mn;
				} // d=1,4
			} // p=0,9
		} // iname=1,4
	} // scope

	Size const max_bb_bins( 3 );
	mean_backbone_torsion_.clear(); mean_backbone_torsion_.resize( 6 );
	for ( Size t=1; t<= 6; ++t ) mean_backbone_torsion_[t].resize( max_bb_bins, 0.0 );


	{ // compute alpha and gamma mean values
		// alpha and gamma parameterization depends on the bin (g+/g-/t)
		// bin is [0:120][120:240][240:360]
		//
		for ( Size t=1; t<=3; ++t ) {
			if ( t==2 ) continue; //not beta
			for ( Size bin=1; bin<= 3; ++bin ) {
				vector1< Real > torsions;
				for ( Size i=1; i<= all_torsions_.size(); ++i ) {
					if ( is_lower_terminus_[i] || is_upper_terminus_[i] ) continue;
					Real const angle( all_torsions_[i][t] );
					if ( get_triple_bin( angle ) != bin ) continue;
					torsions.push_back( angle );
				}
				Real mean, median, sdev;
				get_mean_median_and_sdev( torsions, mean, median, sdev );
				mean_backbone_torsion_[ t ][ bin ] = mean;
				TR.Trace << "mean_backbone_torsion: " << t << ' ' << bin << " ntor: " << I(6,torsions.size() ) <<
					" mean: " <<F(9,3,mean) << " median: " << F(9,3,median) << " sdev: " << F(9,3,sdev) << std::endl;
			}
		}
	}


	{ // compute beta mean values
		// beta depends on the b1/b2 bin for the previous rsd
		Size const t( 2 );
		for ( Size bin=1; bin<= 2; ++bin ) { // B1, B2
			vector1< Real > torsions;
			for ( Size i=1; i<= all_torsions_.size(); ++i ) {
				if ( is_lower_terminus_[i] || is_lower_terminus_[i-1] || is_upper_terminus_[i] ) continue;
				Real const angle( all_torsions_[i][t] ), epsilon( all_torsions_[i-1][5] ), zeta( all_torsions_[i-1][6] );
				if ( get_b1b2_bin( epsilon, zeta ) != bin ) continue;
				torsions.push_back( angle );
			}
			Real mean, median, sdev;
			get_mean_median_and_sdev( torsions, mean, median, sdev );
			mean_backbone_torsion_[ t ][ bin ] = mean;
			TR.Trace << "mean_backbone_torsion: " << t << ' ' << bin << " ntor: " << I(6,torsions.size() ) <<
				" mean: " <<F(9,3,mean) << " median: " << F(9,3,median) << " sdev: " << F(9,3,sdev) << std::endl;
		}
	}


	{ // compute epsilon, zeta mean values
		// they depend on b1/b2 bin of the same residue

		for ( Size t=5; t<= 6; ++t ) {
			for ( Size bin=1; bin<= 2; ++bin ) { // B1, B2
				vector1< Real > torsions;
				for ( Size i=1; i<= all_torsions_.size(); ++i ) {
					if ( is_lower_terminus_[i] || is_upper_terminus_[i] ) continue;
					Real const angle( all_torsions_[i][t] ), epsilon( all_torsions_[i][5] ), zeta( all_torsions_[i][6] );
					if ( get_b1b2_bin( epsilon, zeta ) != bin ) continue;
					torsions.push_back( angle );
				}
				Real mean, median, sdev;
				get_mean_median_and_sdev( torsions, mean, median, sdev );
				mean_backbone_torsion_[ t ][ bin ] = mean;
				TR.Trace << "mean_backbone_torsion: " << t << ' ' << bin << " ntor: " << I(6,torsions.size() ) <<
					" mean: " <<F(9,3,mean) << " median: " << F(9,3,median) << " sdev: " << F(9,3,sdev) << std::endl;
			}
		}
	}

}

void
DNA_DihedralPotential::read_dna_geometry_log_file_from_database(
	std::string const & database_file
)
{
	utility::io::izstream data;
	basic::database::open( data, database_file );

	parse_dna_geometry_log( data );

	data.close();

} // namespace dna

void
DNA_DihedralPotential::read_dna_geometry_log_file(
	std::string const & filename
)
{
	utility::io::izstream data( filename );

	parse_dna_geometry_log( data );

	data.close();
}


} // namespace dna
}} // scoring core
