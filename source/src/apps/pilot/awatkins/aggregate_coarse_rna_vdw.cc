// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Andy Watkins

#include <utility/json_utilities.hh>

// libRosetta headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

#include <core/pose/Pose.hh>

#include <utility/io/ozstream.hh>
#include <utility/tools/make_vector.hh>

#include <protocols/rna/movers/RNA_Coarsify.hh>
#include <numeric/MathNTensor.hh>
#include <numeric/MathNTensor_io.hh>

// C++ headers
#include <iostream>
#include <string>
#include <algorithm>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/options/keys/OptionKeyList.hh>

using namespace core;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::pose;
using namespace core::import_pose;
using namespace core::import_pose::pose_stream;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.awatkins.aggregate_coarse_rna_vdw" );

core::Real
calc_IQR(
	utility::vector1< core::Real > const & values
) {
	//values = std::sort( values.begin(), values.end() );

	core::Size q1 = values.size()/4;
	core::Size q3 = 3*values.size()/4;

	return values[ q3 ] - values[ q1 ];
}

core::Real
fdrule(
	utility::vector1< core::Real > const & values
) {
	core::Real IQR = calc_IQR( values );
	core::Real denominator = pow( values.size(), 1.0/3.0 );

	return 2 * IQR / denominator;
}

class NonceHistogram {
public:

	NonceHistogram( utility::vector1< core::Real > & values /*, core::Size n_bins*/ ) {
		// Figure out bins, using whatever that algorithm is.
		//  Freedman-Diaconis rule: bin-width is 2x IQR / n ^ 1/3
		//  nbins is max-min/binw
		std::sort( values.begin(), values.end() );
		binw_ = fdrule( values );
		core::Size n_bins = core::Size( ( values[ values.size() ] - values[ 1 ] ) / binw_ );

		// One fewer than nbins
		utility::vector1< core::Real > bin_boundaries;
		for ( core::Size ii = 1; ii < n_bins; ++ii ) {
			bin_boundaries.push_back( values[ 1 ] + binw_ * ii );
		}
		// TR << "bin bounds" << bin_boundaries << std::endl;
		// TR << "max value " << values[ values.size() ] << std::endl;
		counts_.resize( n_bins );
		distribution_.resize( n_bins );

		// Put stuff into bins
		for ( auto const & value : values ) {
			bool bin_found = false;
			for ( core::Size ii = 1; ii <= bin_boundaries.size(); ++ii ) {
				auto const bound = bin_boundaries[ ii ];
				if ( value < bound ) {
					++counts_[ ii ];
					bin_found = true;
					break;
				}
			}
			if ( !bin_found ) {
				++counts_[ counts_.size() ];
			}
		}

		// normalize
		core::Size total = values.size();
		for ( core::Size ii = 1; ii <= n_bins; ++ii ) {
			distribution_[ ii ] = Real( counts_[ ii ] ) / Real( total );
		}

	}

	utility::vector1< core::Real > distribution() const { return distribution_; }
	utility::vector1< core::Size > counts() const { return counts_; }
	core::Real binw() const { return binw_; }

private:

	utility::vector1< core::Size > counts_;
	utility::vector1< core::Real > distribution_;
	core::Real binw_ = 0;
};

numeric::MathNTensor< core::Real, 1 >
mathntensor_from_dist( utility::vector1< Real > const & dist ) {
	auto tens = numeric::MathNTensor< core::Real, 1 >( utility::fixedsizearray1< core::Size, 1 >{ dist.size() } );
	// TR << "dist is " << dist << std::endl;
	// find max value of negln dist
	core::Real maxval = 0;
	for ( auto const & d : dist ) {
		core::Real trial = -1.0 * std::log( d );
		if ( trial > maxval ) {
			maxval = trial;
		}
	}
	if ( maxval > 6 ) maxval = 6;
	TR << std::endl << "Maximum statistical potential value is " << maxval << std::endl;

	for ( core::Size ii = 1; ii <= dist.size(); ++ii ) {
		tens( utility::fixedsizearray1< core::Size, 1 >{ ii - 1 } ) = -1.0 * std::log( dist[ ii ] ) - maxval;
		if ( -1.0 * std::log( dist[ ii ] ) - maxval > 0 ) {
			tens( utility::fixedsizearray1< core::Size, 1 >{ ii - 1 } ) = 0;
		}
		TR << tens( utility::fixedsizearray1< core::Size, 1 >{ ii - 1 } ) << " ";
	}
	TR << std::endl;

	return tens;
}

///////////////////////////////////////////////////////////////
void
my_main()
{
	using namespace basic::options;

	utility::io::ozstream out1( "vdw.txt" );
	utility::io::ozstream out2( "vdw_1_1.txt" );
	utility::io::ozstream out3( "vdw_1_2.txt" );
	utility::io::ozstream out4( "vdw_1_3.txt" );
	utility::io::ozstream out5( "vdw_2_2.txt" );
	utility::io::ozstream out6( "vdw_2_3.txt" );
	utility::io::ozstream out7( "vdw_3_3.txt" );

	utility::vector1< core::Real > values_P_P;
	utility::vector1< core::Real > values_P_S;
	utility::vector1< core::Real > values_P_CEN;
	utility::vector1< core::Real > values_S_S;
	utility::vector1< core::Real > values_S_CEN;
	utility::vector1< core::Real > values_CEN_CEN;

	// Read poses in and coarsify them all.
	// setup input stream
	PoseInputStreamOP input;
	if ( option[ in::file::tags ].user() ) {
		input = utility::pointer::make_shared< SilentFilePoseInputStream >(
			option[ in::file::silent ](),
			option[ in::file::tags ]()
		);
	} else if ( option[ in::file::silent ].user() ) {
		input = utility::pointer::make_shared< SilentFilePoseInputStream >(
			option[ in::file::silent ]()
		);
	} else {
		input = utility::pointer::make_shared< PDBPoseInputStream >( option[ in::file::s ]() );
	}

	auto rna_coarsify = utility::pointer::make_shared< protocols::rna::movers::RNA_Coarsify >();
	// iterate over input stream
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	while ( input->has_another_pose() ) {
		Pose pose;
		input->fill_pose( pose, *rsd_set, false );
		Real const CUTOFF( 12.0 );

		Size const NATOMS( 3 ); //P, S, CEN
		for ( Size i = 1; i <= pose.size(); i++ ) {
			Residue const & rsd1 = pose.residue( i );

			for ( Size j = i+2; j <= pose.size(); j++ ) {
				Residue const & rsd2 = pose.residue( j );

				utility::vector1< Real > dists;
				Real min_dist( CUTOFF );
				for ( Size k = 1; k <= NATOMS; k++ ) {
					for ( Size m = 1; m <= NATOMS; m++ ) {

						Real const dist = ( rsd1.xyz( k ) - rsd2.xyz( m ) ).length();
						dists.push_back( dist );
						if ( dist < min_dist ) min_dist = dist;

					}
				}

				if ( min_dist < CUTOFF ) {
					for ( Size k = 1; k <= dists.size(); k++ ) {
						out1 << ' ' << dists[ k ];
					}
					out1 << std::endl;

					// 1-1, 1-2, 1-3, 2-1, 2-2, 2-3, 3-1, 3-2, 3-3
					out2 << dists[ 1 ] << std::endl;
					values_P_P.push_back( dists[ 1 ] );
					out3 << dists[ 2 ] << ' ' << dists[ 4 ] << std::endl;
					values_P_S.push_back( dists[ 2 ] );
					values_P_S.push_back( dists[ 4 ] );
					out4 << dists[ 3 ] << ' ' << dists[ 7 ] << std::endl;
					values_P_CEN.push_back( dists[ 3 ] );
					values_P_CEN.push_back( dists[ 7 ] );
					out5 << dists[ 5 ] << std::endl;
					values_S_S.push_back( dists[ 5 ] );
					out6 << dists[ 6 ] << ' ' << dists[ 8 ] << std::endl;
					values_S_CEN.push_back( dists[ 6 ] );
					values_S_CEN.push_back( dists[ 8 ] );
					out7 << dists[ 9 ] << std::endl;
					values_CEN_CEN.push_back( dists[ 9 ] );

				}

			}
		}
	}

	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	out6.close();
	out7.close();

	auto P_P_hist = NonceHistogram( values_P_P );
	auto P_S_hist = NonceHistogram( values_P_S );
	auto P_CEN_hist = NonceHistogram( values_P_CEN );
	auto S_S_hist = NonceHistogram( values_S_S );
	auto S_CEN_hist = NonceHistogram( values_S_CEN );
	auto CEN_CEN_hist = NonceHistogram( values_CEN_CEN );

	auto P_P_dist = P_P_hist.distribution();
	auto P_S_dist = P_S_hist.distribution();
	auto P_CEN_dist = P_CEN_hist.distribution();
	auto S_S_dist = S_S_hist.distribution();
	auto S_CEN_dist = S_CEN_hist.distribution();
	auto CEN_CEN_dist = CEN_CEN_hist.distribution();

	// TR << "P_P counts     " << P_P_hist.counts() << std::endl;
	// TR << "P_S counts     " << P_S_hist.counts() << std::endl;
	// TR << "P_CEN counts   " << P_CEN_hist.counts() << std::endl;
	// TR << "S_S counts     " << S_S_hist.counts() << std::endl;
	// TR << "S_CEN counts   " << S_CEN_hist.counts() << std::endl;
	// TR << "CEN_CEN counts " << CEN_CEN_hist.counts() << std::endl;

	auto P_P_tensor = mathntensor_from_dist( P_P_dist );
	auto P_S_tensor = mathntensor_from_dist( P_S_dist );
	auto P_CEN_tensor = mathntensor_from_dist( P_CEN_dist );
	auto S_S_tensor = mathntensor_from_dist( S_S_dist );
	auto S_CEN_tensor = mathntensor_from_dist( S_CEN_dist );
	auto CEN_CEN_tensor = mathntensor_from_dist( CEN_CEN_dist );

	using namespace utility::json_spirit;
	using namespace utility::tools;
	Value P_P_val = make_vector( Pair( "binwidth", Array({P_P_hist.binw()} )), Pair( "minval", Array({values_P_P[1]} )), Pair( "maxval", Array({values_P_P[values_P_P.size()]} )), Pair( "n_bins", Array({P_P_dist.size()+1}) ) );
	Value P_S_val = make_vector( Pair( "binwidth", Array({P_S_hist.binw()} )), Pair( "minval", Array({values_P_S[1]} )), Pair( "maxval", Array({values_P_S[values_P_S.size()]} )), Pair( "n_bins", Array({P_S_dist.size()+1}) ) );
	Value P_CEN_val = make_vector( Pair( "binwidth", Array({P_CEN_hist.binw()} )), Pair( "minval", Array({values_P_CEN[1]} )), Pair( "maxval", Array({values_P_CEN[values_P_CEN.size()]} )), Pair( "n_bins", Array({P_CEN_dist.size()+1}) ) );
	Value S_S_val = make_vector( Pair( "binwidth", Array({S_S_hist.binw()} )), Pair( "minval", Array({values_S_S[1]} )), Pair( "maxval", Array({values_S_S[values_S_S.size()]} )), Pair( "n_bins", Array({S_S_dist.size()+1}) ) );
	Value S_CEN_val = make_vector( Pair( "binwidth", Array({S_CEN_hist.binw()} )), Pair( "minval", Array({values_S_CEN[1]} )), Pair( "maxval", Array({values_S_CEN[values_S_CEN.size()]} )), Pair( "n_bins", Array({S_CEN_dist.size()+1}) ) );
	Value CEN_CEN_val = make_vector( Pair( "binwidth", Array({CEN_CEN_hist.binw()} )), Pair( "minval", Array({values_CEN_CEN[1]} )), Pair( "maxval", Array({values_CEN_CEN[values_CEN_CEN.size()]} )), Pair( "n_bins", Array({CEN_CEN_dist.size()+1}) ) );

	numeric::write_tensor_to_file( "dist_P_P.bin", P_P_tensor, P_P_val );
	numeric::write_tensor_to_file( "dist_P_S.bin", P_S_tensor, P_S_val );
	numeric::write_tensor_to_file( "dist_P_CEN.bin", P_CEN_tensor, P_CEN_val );
	numeric::write_tensor_to_file( "dist_S_S.bin", S_S_tensor, S_S_val );
	numeric::write_tensor_to_file( "dist_S_CEN.bin", S_CEN_tensor, S_CEN_val );
	numeric::write_tensor_to_file( "dist_CEN_CEN.bin", CEN_CEN_tensor, CEN_CEN_val );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence>  [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::native );
		option.add_relevant( in::file::input_res );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( out::overwrite );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		my_main();

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}

