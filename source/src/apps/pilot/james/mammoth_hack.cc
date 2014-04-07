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

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SWAligner.hh>

#include <utility/vector1.hh>

#include <ObjexxFCL/string.functions.hh>

#include <numeric/model_quality/util.hh>
#include <numeric/model_quality/rms.hh>

using core::Size;
using core::Real;
using utility::vector1;

#include <apps/pilot/james/james_util.hh>

// C++ headers
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
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
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <utility/stream_util.hh>
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
#include <utility/io/izstream.fwd.hh>
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
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/RmsData.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/CArray.fwd.hh>
#include <ObjexxFCL/ChunkVector.fwd.hh>
#include <ObjexxFCL/Cstring.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.all.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.all.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3P.fwd.hh>
#include <ObjexxFCL/FArray4.all.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4P.fwd.hh>
#include <ObjexxFCL/FArray5.all.fwd.hh>
#include <ObjexxFCL/FArray5.fwd.hh>
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArray5P.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/ObjexxFCL.Project.hh>
#include <ObjexxFCL/ObjexxFCL.fwd.hh>
#include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/sbyte.fwd.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
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


class SequenceCoords;
typedef utility::pointer::owning_ptr< SequenceCoords > SequenceCoordsOP;

class SequenceCoords : public core::sequence::Sequence {
public:
	SequenceCoords(
		core::pose::Pose const & pose
	) {
		std::string const atom_name("CA");
		coords_.clear();

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			core::id::NamedAtomID atom_ii( atom_name, ii );
			coords_.push_back( pose.xyz( atom_ii ) );
			//std::cerr << "idx " << ii << ": "
			//	<< coords_[ii].x() << "," << coords_[ii].y() << "," << coords_[ii].z() << std::endl;
		}
		//std::cerr << "have " << coords_.size() << " coords." << std::endl;
		sequence( pose.sequence() );
	}

	core::PointPosition coord( Size const idx ) const {
		runtime_assert( idx <= coords_.size() );
		//std::cerr << "idx " << idx << " returning "
		//	<< coords_[idx].x() << "," << coords_[idx].y() << "," << coords_[idx].z()
		//	<< std::endl;
		return coords_[ idx ];
	}

	utility::vector1< core::PointPosition> coords() const {
		return coords_;
	}

	void coords( utility::vector1< core::PointPosition > coords ) {
		coords_ = coords;
	}

	void print() const {
		std::cerr << "SEQUENCE: " << sequence() << std::endl;
		for ( Size ii = 1; ii <= coords_.size(); ++ii ) {
			std::cerr << "idx " << ii << ": "
				<< coords_[ii].x() << "," << coords_[ii].y() << "," << coords_[ii].z()
					<< std::endl;
		}
	}

	virtual core::sequence::SequenceOP clone() const {
		return new SequenceCoords( *this );
	}

	  /// @brief copy ctor
	SequenceCoords( SequenceCoords const & src ):
		Sequence()
	{
		*this = src;
	}

	/// @brief assignment operator.
	SequenceCoords & operator = ( SequenceCoords const & rhs ) {
		if ( this == &rhs ) return *this;
		sequence( rhs.sequence() );
		gap_char( rhs.gap_char() );
		start   ( rhs.start() );
		coords  ( rhs.coords() );
		id      ( rhs.id() );
		return *this;
	}

private:
	utility::vector1< core::PointPosition > coords_;
}; // class SequenceCoords

class RMS_ScoringScheme : public core::sequence::ScoringScheme {
public:
	RMS_ScoringScheme( Size const len ) : length_( len ) {}
	core::sequence::ScoringSchemeOP clone() const {
		return new RMS_ScoringScheme(*this);
	}

	Real score(
		core::sequence::SequenceOP seq1,
		core::sequence::SequenceOP seq2,
		Size pos1,
		Size pos2
	) {
		SequenceCoordsOP coords1 = SequenceCoordsOP(
			static_cast < SequenceCoords * > ( seq1() )
		);
		SequenceCoordsOP coords2 = SequenceCoordsOP(
			static_cast < SequenceCoords * > ( seq2() )
		);

		if ( length() + pos1 > seq1->length() || length() + pos2 > seq2->length() ) {
			//std::cerr << "returning zero!" << std::endl;
			return 0;
		}

		utility::vector1< core::PointPosition > p1, p2;
		for ( Size ii = 1; ii <= length(); ++ii ) {
			p1.push_back( coords1->coord( pos1 + ii - 1 ) );
			p2.push_back( coords2->coord( pos2 + ii - 1 ) );
		}

		//for ( Size ii = 1; ii <= length(); ++ii ) {
		//	std::cerr << "p1 " << p1[ii].x() << "," << p1[ii].y() << "," << p1[ii].z() << std::endl;
		//	std::cerr << "p2 " << p2[ii].x() << "," << p2[ii].y() << "," << p2[ii].z() << std::endl;
		//}

		using std::min;
		using std::sqrt;
		//static Real const random_rms( sqrt( 1.0 - 2.84 / sqrt( length() ) ) );
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		static Real const random_rms( 6 );
		Real const rms( numeric::model_quality::calc_rms( p1, p2 ) );
		Real score( 0.0 );
		if ( random_rms > rms ) {
			score = ( random_rms - rms ) / random_rms * 10;
		}

		//Real const score( std::min( 10.0, 10.0 - rms ) );
		//std::cout << "rms = " << rms << ", random_rms = " << random_rms << std::endl;
		//std::cout << "score(" << pos1 << "," << pos2 << ") = " << score << std::endl;
		return score;
	}

	Size length() const {
		return length_;
	}

private:
	Size length_;
}; // RMS_ScoringScheme

utility::vector1< core::PointPosition >
forward_in_time(
	FArray2D< numeric::Real > p1,
	int const natoms
) {
	//FArray2D< numeric::Real > retval( 3, p1_coords.size() );
	utility::vector1< core::PointPosition > retval;
	for ( int ii = 1; ii <= natoms; ++ii ) {
		core::PointPosition temp;
		temp.x( p1(ii,1) );
		temp.y( p1(ii,2) );
		temp.z( p1(ii,3) );
		retval.push_back( temp );
	} // i
	return retval;
}

FArray2D< core::Real > back_in_time_2d(
	utility::vector1< core::PointPosition > p1
) {

	FArray2D< numeric::Real > retval( 3, p1.size() );
	for ( Size ii = 1; ii <= p1.size(); ++ii ) {
		for ( Size kk = 1; kk <= 3; ++kk ) { // k = X, Y and Z
			retval(kk,ii) = p1[ii][kk];
		} // k
	} // i

	return retval;
}

FArray1D< core::Real > back_in_time_1d(
	utility::vector1< numeric::Real > w
) {
	int const natoms( w.size() );
	FArray1D< numeric::Real > retval( natoms, 1.0 );
	for ( int ii = 1; ii <= natoms; ++ii ) {
		retval(ii) = w[ii];
	}

	return retval;
}

core::id::SequenceMapping extend_mapping_by_threshold(
	utility::vector1< core::PointPosition > & /*p1*/,
	utility::vector1< core::PointPosition > & /*p2*/,
	core::id::SequenceMapping const & map
) {
	for ( Size ii = 1; ii <= map.size1(); ++ii ) {

	}
	core::id::SequenceMapping extended_map( map );
	return extended_map;
}

void apply_rotation(
	utility::vector1< core::PointPosition > & p2,
	FArray2D< numeric::Real > uu
) {
	// rotate p2f onto p1f using p2f x uu
	// DEFINITELY CHECK THIS! WRITTEN AT 4AM 4/22/09
	for ( Size ii = 1; ii <= p2.size(); ++ii ) {
		p2[ii].x( p2[ii].x() * uu(1,1) + p2[ii].y() * uu(2,1) + p2[ii].z() * uu(3,1) );
		p2[ii].y( p2[ii].x() * uu(1,2) + p2[ii].y() * uu(2,2) + p2[ii].z() * uu(3,2) );
		p2[ii].z( p2[ii].x() * uu(1,3) + p2[ii].y() * uu(2,3) + p2[ii].z() * uu(3,3) );
	}
}

/// @brief Superimposes p2 onto p1, using the SequenceMapping in map.
FArray2D< numeric::Real > optimal_rotation(
	utility::vector1< core::PointPosition > & p1,
	utility::vector1< core::PointPosition > & p2
) {
	FArray2D< numeric::Real > p1f = back_in_time_2d( p1 );
	FArray2D< numeric::Real > p2f = back_in_time_2d( p2 );
	int const natoms( p1.size() );
	//FArray1D< numeric::Real > ww = back_in_time_1d( weights );
	FArray1D< numeric::Real > ww( natoms, 1.0 );
	FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	numeric::model_quality::findUU( p1f, p2f, ww, natoms, uu, ctx );
	return uu;
}

/// @brief Aligns the points in p1 with the points in p2. The alignment
/// in the provided SequenceMapping is used as to seed the algorithm, and
/// the locally optimal alignment between p1 and p2 is returned.
core::id::SequenceMapping maxsub_align(
	utility::vector1< core::PointPosition > & p1,
	utility::vector1< core::PointPosition > & p2,
	core::id::SequenceMapping const & map
) {
	core::id::SequenceMapping new_map ( map );
	core::id::SequenceMapping last_map( map );
	numeric::model_quality::center_atoms_at_origin< numeric::Real > ( p1 );
	numeric::model_quality::center_atoms_at_origin< numeric::Real > ( p2 );

	while( !( new_map == last_map ) ) {
		// create copies of p1 and p2 based on map
		utility::vector1< core::PointPosition > current_p1;
		utility::vector1< core::PointPosition > current_p2;
		for ( Size ii = 1; ii <= map.size1(); ++ii ) {
			Size const jj( new_map[ii] );
			if ( jj != 0 ) {
				current_p1.push_back( p1[ ii ] );
				current_p2.push_back( p2[ jj ] );
			}
		}

		// calculate optimal superposition of p1 onto p2
		FArray2D< numeric::Real > uu = optimal_rotation( current_p1, current_p2 );
		apply_rotation( p2, uu );

		// update last_map and new_map
		last_map = new_map;
		new_map = extend_mapping_by_threshold( p1, p2, map );
	}

	return new_map;
}

int
main( int argc, char * argv [] ) {
	try {

	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::import_pose::pose_stream;

	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using std::string;

	devel::init( argc, argv );

	// setup residue types
	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	// read in native and template poses
	core::pose::Pose native_pose, template_pose;
	core::import_pose::pose_from_pdb(
		native_pose, *rsd_set, option[ in::file::native ]()
	);
	core::import_pose::pose_from_pdb(
		template_pose, *rsd_set, option[ in::file::template_pdb ]()[1]
	);

	protocols::jumping::assign_ss_dssp( template_pose );
	protocols::jumping::assign_ss_dssp( native_pose   );

	SequenceCoordsOP template_seq( new SequenceCoords( template_pose ) );
	template_seq->id( option[ in::file::template_pdb ]()[1] );
	SequenceCoordsOP native_seq( new SequenceCoords( native_pose ) );
	native_seq->id( option[ in::file::native ]() );

	std::cout << "template " << template_pose.secstruct() << std::endl;
	std::cout << "native   " << native_pose.secstruct() << std::endl;
	std::cout << "template " << template_pose.sequence() << std::endl;
	std::cout << "native   " << native_pose.sequence() << std::endl;
	SequenceOP template_ss(
		new Sequence( template_pose.secstruct(), template_seq->id(), template_seq->start() )
	);
	SequenceOP native_ss(
		new Sequence( native_pose.secstruct(), native_seq->id(), native_seq->start() )
	);

	std::cout << template_pose.secstruct() << std::endl;
	std::cout << native_pose.secstruct() << std::endl;

	SWAligner sw_align;
	ScoringSchemeOP ss( new RMS_ScoringScheme( 7 ) );
	ss->gap_open  ( -7.00 );
	ss->gap_extend( -0.35 );
	SequenceAlignment current_align(
		sw_align.align( native_seq, template_seq, ss )
	);

	// iteratively update alignment using Maxsub algorithm

	//string outfile( option[ out::file::alignment ]() );
	//std::ofstream output( outfile.c_str() );
	//if ( !output.is_open() ) {
	//	utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	//}

	//output << current_align << std::endl;
	std::cout << current_align << std::endl;
	//using namespace core::sequence;
	//SequenceAlignment ss_align(
	//	mapping_to_alignment(
	//		current_align.sequence_mapping(1,2),
	//		native_ss,
	//		template_ss
	//	)
	//);
	//std::cout << ss_align << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
