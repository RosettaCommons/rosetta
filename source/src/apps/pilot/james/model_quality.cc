// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <apps/pilot/james/james_util.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintForest.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/FuncFactory.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <utility/excn/Exceptions.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

int
main( int argc, char * argv [] ) {
	try {

	using core::Size;
	using core::Real;
	using utility::vector1;
	using std::string;
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	devel::init( argc, argv );

	string const atom_name( "CA" );
	int output_width = 16;
	int precision = 4;
	Real dist_upper( 10.0 );
	Size const min_seqsep( 3 );

	utility::vector1< core::Real > dist_thresholds;
	utility::vector1< core::Real > torsion_thresholds;
	if ( option[ james::dist_thresholds ].user() ) {
		dist_thresholds = option[ james::dist_thresholds ]();
	} else {
		dist_thresholds.clear();
		dist_thresholds.push_back(  1.0 );
		dist_thresholds.push_back(  2.0 );
		dist_thresholds.push_back(  5.0 );
	}
	if ( option[ james::torsion_thresholds ].user() ) {
		torsion_thresholds = option[ james::torsion_thresholds ]();
	} else {
		torsion_thresholds.clear();
		torsion_thresholds.push_back(  15.0 );
		torsion_thresholds.push_back(  30.0 );
		torsion_thresholds.push_back(  60.0 );
		torsion_thresholds.push_back(  90.0 );
	}

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set(
			option[ in::file::residue_type_set ]()
		);

	// read in a native pose
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
		);
	} else {
		utility_exit_with_message( "Error: must provide in::file::native!" );
	}
	string outfile_prefix = option[ out::file::silent ]();
	string dist_outfile = outfile_prefix + ".dist_features";
	string torsion_outfile = outfile_prefix + ".torsion_features";
	std::ofstream dist_output( dist_outfile.c_str() );
	std::ofstream torsion_output( torsion_outfile.c_str() );
	if ( ! dist_output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + dist_outfile + '\n' );
	}
	if ( ! torsion_output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + torsion_outfile + '\n' );
	}

	TorsionList native_torsions( native_pose );
	DistanceMatrix native_dist( native_pose, atom_name );
	Size const nres( native_pose.total_residue() );
	core::import_pose::pose_stream::MetaPoseInputStream input
		= core::import_pose::pose_stream::streams_from_cmd_line();
	Size const n_torsions( native_torsions.n_torsions() );

	// torsion_counts[ thresholds ][ n_torsions ][ nres ]
	vector1< vector1< vector1< Size > > > torsion_counts(
		torsion_thresholds.size(), vector1< vector1< Size > > ( n_torsions, vector1< Size >( nres, 0 ) )
	);

	// distance_counts[ thresholds ][ nres ][ nres ]
	vector1< vector1< vector1< Size > > > distance_counts(
		dist_thresholds.size(), vector1< vector1< Size > >( nres, vector1< Size >( nres, 0 ) )
	);
	Size decoy_count( 0 );
	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		runtime_assert( pose.sequence() == native_pose.sequence() );

		TorsionList pose_torsions( pose );
		DistanceMatrix pose_dist( pose, atom_name );
		for ( Size ii = 1; ii <= nres; ++ii ) {
			for ( Size jj = ii + 1; jj <= nres; ++jj ) {
				if ( native_dist.distance( ii, jj ) >= dist_upper ) continue;
				if ( (jj - ii) < min_seqsep ) continue;
				Real const dist_delta(
					std::abs(
						native_dist.distance( ii, jj ) - pose_dist.distance( ii, jj )
					)
				);
				for ( Size kk = 1; kk <= dist_thresholds.size(); ++kk ) {
					if ( dist_delta <= dist_thresholds[kk] ) distance_counts[kk][ii][jj]++;
				}
			} // jj

			for ( Size kk = 1; kk <= n_torsions; ++kk ) {
				Real const torsion_delta( std::abs(
					native_torsions.torsion( ii, kk ) - pose_torsions.torsion( ii, kk )
				) );

				for ( Size ll = 1; ll <= torsion_thresholds.size(); ++ll ) {
					if ( torsion_delta <= torsion_thresholds[ll] )
					++torsion_counts[ll][kk][ii];
				}
			}
		} // ii
		++decoy_count;
	} // while ( input.has_another_pose() )

	dist_output << A( output_width, "resi" ) << A( output_width, "resj" )
		<< A( output_width, "native_dist" );
	for ( Size ii = 1; ii <= dist_thresholds.size(); ++ii ) {
		dist_output << A( output_width, "ndec_" + string_of( dist_thresholds[ii] ) );
	}
	dist_output << std::endl;

	for ( Size ii = 1; ii <= nres; ++ii ) {
		for ( Size jj = ii + 1; jj <= nres; ++jj ) {
			if ( native_dist.distance( ii, jj ) >= dist_upper ) continue;
			if ( (jj - ii) < min_seqsep ) continue;
			dist_output << I( output_width, ii ) << I( output_width, jj )
				<< F( output_width, precision, native_dist.distance( ii, jj ) );
			for ( Size kk = 1; kk <= dist_thresholds.size(); ++kk ) {
				dist_output << I( output_width, distance_counts[kk][ii][jj] );
			}
			dist_output << std::endl;
		}
	}

	torsion_output << A( output_width, "resi" ) << A( output_width, "tors_idx" )
		<< A( output_width, "nat_tors" );
 	for ( Size tt = 1; tt <= torsion_thresholds.size(); ++tt ) {
		torsion_output << A( output_width, "tors_" + string_of( torsion_thresholds[tt] ) );
	}
	torsion_output << std::endl;
	for ( Size ii = 1; ii <= nres; ++ii ) {
		for ( Size kk = 1; kk <= n_torsions; ++kk ) {
			torsion_output << I( output_width, ii )  << I( output_width, kk )
				<< F( output_width, precision, native_torsions.torsion( ii, kk ) ) ;

		 	for ( Size tt = 1; tt <= torsion_thresholds.size(); ++tt ) {
				torsion_output << I( output_width, torsion_counts[tt][kk][ii] );
			}

			torsion_output << std::endl;
		}
	}

	dist_output 	<< "# total of " << decoy_count << " decoys." << std::endl;
	torsion_output 	<< "# total of " << decoy_count << " decoys." << std::endl;

	dist_output.close();
	torsion_output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // int main
