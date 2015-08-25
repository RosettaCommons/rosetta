// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh
/// @brief   Implementation of semi-rotameric rotamer libraries from Jun08
/// @author  Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_tmpl_hh
#define INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_tmpl_hh

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

// Unit Headers
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>

// Project Headers
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>

#include <core/conformation/Residue.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Basic Headers
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.functions.hh>
#include <utility/backtrace.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/MathTensor.hh>
#include <numeric/interpolation/spline/PolycubicSpline.hh>
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/util.hh>

// Boost Headers
#include <boost/cstdint.hpp>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/LexicographicalIterator.fwd.hh>
#include <utility/LexicographicalIterator.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/mpistream.ipp>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
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
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.all.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray2.all.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3P.fwd.hh>
#include <ObjexxFCL/FArray4.all.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4.hh>
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/FArray4P.fwd.hh>
#include <ObjexxFCL/FArray5.all.fwd.hh>
#include <ObjexxFCL/FArray5.fwd.hh>
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArray5P.fwd.hh>
#include <ObjexxFCL/FArray6.all.fwd.hh>
#include <ObjexxFCL/FArray6.fwd.hh>
#include <ObjexxFCL/FArray6A.fwd.hh>
#include <ObjexxFCL/FArray6D.fwd.hh>
#include <ObjexxFCL/FArray6P.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/KeyFArray4D.fwd.hh>
#include <ObjexxFCL/KeyFArray5D.fwd.hh>
#include <ObjexxFCL/KeyFArray6D.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/interpolate.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

namespace core {
namespace pack {
namespace dunbrack {

template < class P >
bool
BBDepNRChiSample< P >::operator==( BBDepNRChiSample<P> const & other ) const {
	static thread_local basic::Tracer TR( "core.pack.dunbrack.BBDepNRChiSample" );

	core::Real const & ANGLE_DELTA = core::pack::dunbrack::SingleResidueDunbrackLibrary::ANGLE_DELTA;
	core::Real const & PROB_DELTA = core::pack::dunbrack::SingleResidueDunbrackLibrary::PROB_DELTA;

	bool equal( true );

	if ( packed_rotno_ != other.packed_rotno_ ) {
		TR.Debug << "Comparision error, packed_rotno - " << packed_rotno_ << " vs. " << other.packed_rotno_ << std::endl;
		equal = false;
	}
	if ( nrchi_bin_ != other.nrchi_bin_ ) {
		TR.Debug << "Comparision error, nrchi_bin - " << nrchi_bin_ << " vs. " << other.nrchi_bin_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( nrchi_mean_, other.nrchi_mean_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparision error, nrchi_mean - " << nrchi_mean_ << " vs. " << other.nrchi_mean_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( nrchi_sd_, other.nrchi_sd_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparision error, nrchi_sd - " << nrchi_sd_ << " vs. " << other.nrchi_sd_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( prob_, other.prob_, PROB_DELTA ) ) {
		TR.Debug << "Comparision error, prob - " << prob_ << " vs. " << other.prob_ << std::endl;
		equal = false;
	}

	return equal;
}

template < class P >
bool
BBIndNRChiSample< P >::operator==( BBIndNRChiSample<P> const & other ) const {
	static thread_local basic::Tracer TR( "core.pack.dunbrack.BBIndNRChiSample" );

	core::Real const & ANGLE_DELTA = core::pack::dunbrack::SingleResidueDunbrackLibrary::ANGLE_DELTA;
	core::Real const & PROB_DELTA = core::pack::dunbrack::SingleResidueDunbrackLibrary::PROB_DELTA;

	bool equal( true );

	if ( ! numeric::equal_by_epsilon( left_, other.left_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparision error, left - " << left_ << " vs. " << other.left_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( median_, other.median_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparision error, median - " << median_ << " vs. " << other.median_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( right_, other.right_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparision error, right - " << right_ << " vs. " << other.right_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( prob_, other.prob_, PROB_DELTA ) ) {
		TR.Debug << "Comparision error, prob - " << prob_ << " vs. " << other.prob_ << std::endl;
		equal = false;
	}

	return equal;
}

template < Size N >
inline
bool
BBDepScoreInterpData< N >::operator==( BBDepScoreInterpData< N > const & other ) const {
	static thread_local basic::Tracer TR( "core.pack.dunbrack.BBDepScoreInterpData" );

	core::Real const & COEF_DELTA = core::pack::dunbrack::SingleResidueDunbrackLibrary::COEF_DELTA;

	bool equal( true );

	for ( Size i = 1; i <= ( 1 << ( N + 1 ) ); ++i ) {
		if ( ! numeric::equal_by_epsilon( n_derivs_[ i ], other.n_derivs_[ i ], COEF_DELTA ) ) {
			TR.Debug << "Comparision error, deriv " << i << " - " << n_derivs_[ i ] << " vs. " << other.n_derivs_[ i ] << std::endl;
			equal = false;
		}
	}

	return equal;
}


template < Size T, Size N >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::SemiRotamericSingleResidueDunbrackLibrary(
	chemical::AA const aa_in,
	bool const backbone_independent_scoring,   // true uses less memory
	bool const backbone_independent_rotamer_sampling // true uses less memory
) :
parent( aa_in, false /*dun02*/ ), // false, since the semi rotameric library is new in 2008
bbind_nrchi_scoring_( backbone_independent_scoring ),
bbind_nrchi_sampling_( backbone_independent_rotamer_sampling ),
nrchi_periodicity_( 0.0 ),
nrchi_lower_angle_( 0.0 ),
bbdep_nrchi_nbins_( 0 ),
bbdep_nrchi_binsize_( 0 ),
bbind_nrchi_nbins_( 0 ),
bbind_nrchi_binsize_( 0 ),
n_nrchi_sample_bins_( 0 )
{}

template < Size T, Size N >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::~SemiRotamericSingleResidueDunbrackLibrary() throw () {}

template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	return rotamer_energy_deriv_bbdep( rsd, scratch, false );
}

template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	return rotamer_energy_deriv_bbdep( rsd, scratch, true );
}

template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::rotamer_energy_deriv_bbdep(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const
{
	debug_assert( ! bbind_nrchi_scoring_ );

	//Multiplier for D-amino acids:
	const core::Real d_multiplier = rsd.has_property( "D_AA" ) ? -1.0 : 1.0;

	parent::eval_rotameric_energy_deriv( rsd, scratch, eval_deriv ); //This has been updated to allow scoring of D-amino acids.

	Real chidevpen_score( 0.0 );
	for ( Size ii = 1; ii <= T; ++ii ) chidevpen_score += scratch.chidevpen()[ ii ];

	Real nrchi_score, dnrchiscore_dchi;
	utility::fixedsizearray1< Real, N > dnrchiscore_dbb;
	nrchi_score = bbdep_nrchi_score( rsd, scratch, dnrchiscore_dchi, dnrchiscore_dbb ); //Handles D-amino acids; derivatives WRT phi, psi, and chi are inverted iff D-amino acid.

	if ( nrchi_score != nrchi_score ) { //NaN check
		nrchi_score = 0;
		std::cerr << "NaN in SemiRot: " << rsd.seqpos() << " " << rsd.name() << std::endl;
	}

	// Corrections for Shanon Entropy
	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
		nrchi_score += scratch.entropy();
	}

	scratch.fa_dun_tot() = chidevpen_score + nrchi_score;
	scratch.fa_dun_rot() = 0;
	scratch.fa_dun_semi() = nrchi_score;

	//scratch.fa_dun_dev() = chidevpensum; // keep original chidev score...

	//std::cout << "SRSRDL bbdep score: " << rsd.name() << " " << rsd.mainchain_torsion( 1 ) << " " << rsd.mainchain_torsion( 2)  << " " << rsd.chi()[ T + 1 ] << ": " << nrchi_score << " " << chidevpen_score << " " << std::endl;

	if ( ! eval_deriv ) return chidevpen_score + nrchi_score;

	/// sum derivatives.
	Real5 & dE_dbb(  scratch.dE_dbb() );
	Real5 & dE_dbb_dev(  scratch.dE_dbb_dev() );
	Real5 & dE_dbb_rot(  scratch.dE_dbb_rot() );
	Real5 & dE_dbb_semi(  scratch.dE_dbb_semi() );
	Real4 & dE_dchi( scratch.dE_dchi() );
	Real4 & dE_dchi_dev( scratch.dE_dchi_dev() );
	Real4 & dE_dchi_semi( scratch.dE_dchi_semi() );
	std::fill( dE_dchi.begin(), dE_dchi.end(), 0.0 );
	std::fill( dE_dchi_dev.begin(), dE_dchi_dev.end(), 0.0 );
	std::fill( dE_dchi_semi.begin(), dE_dchi_semi.end(), 0.0 );

	for ( Size initi = 1; initi <= N; ++initi ) {
		dE_dbb[ initi ]      = 0.0;
		dE_dbb_dev[ initi ]  = 0.0;
		dE_dbb_rot[ initi ]  = 0.0;
		dE_dbb_semi[ initi ] = 0.0;
	}

	for ( Size bbi = 1; bbi <= N; ++bbi ) {
		dE_dbb[ bbi ]     = d_multiplier * scratch.dchidevpen_dbb()[ bbi ];
		dE_dbb_dev[ bbi ] = d_multiplier * scratch.dchidevpen_dbb()[ bbi ];

		// Correction for entropy
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
			dE_dbb[ bbi ]    += d_multiplier * scratch.dentropy_dbb()[ bbi ];
			dE_dbb_semi[ bbi ] += d_multiplier * scratch.dentropy_dbb()[ bbi ];
		}
	}
	for ( Size bbi = 1; bbi <= N; ++bbi ) {
		dE_dbb[ bbi ]    += d_multiplier * dnrchiscore_dbb[ bbi ];
		dE_dbb_semi[ bbi ] += d_multiplier * dnrchiscore_dbb[ bbi ];
	}

	for ( Size i = 1; i <= T; ++i ) {
		dE_dchi[ i ]  = d_multiplier * scratch.dchidevpen_dchi()[ i ];
		dE_dchi_dev[ i ] = d_multiplier * scratch.dchidevpen_dchi()[ i ];
	}
	dE_dchi[ T + 1 ]   = d_multiplier * dnrchiscore_dchi;
	dE_dchi_semi[ T + 1 ] = d_multiplier * dnrchiscore_dchi;

	parent::correct_termini_derivatives( rsd, scratch );

	return chidevpen_score + nrchi_score;
}

/// @brief Trilinear interpolation.  Derivatives discontinuous at edge planes.
/// Correction: it seems that this was updated at some point to use tricubic interpolation.
/// Note: This function has been updated to handle D-amino acids properly.
template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::bbdep_nrchi_score(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	Real & dnrchi_score_dnrchi,
	utility::fixedsizearray1< Real, N > & dnrchi_score_dbb
) const
{
	debug_assert( ! bbind_nrchi_scoring_ );

	Size const packed_rotno( parent::rotwell_2_packed_rotno( scratch.rotwell() ));

	Real nrchi = rsd.chi( T + 1 );
	core::Real d_multiplier = 1.0; //A multiplier for D-amino acid angles.  -1.0 if D-amino acid, 1.0 otherwise.

	if ( rsd.has_property( "D_AA" ) ) {
		nrchi *= -1.0; //Invert chi if this is a D-amino acid
		d_multiplier = -1.0;
	}
	Size nrchi_bin, nrchi_bin_next;
	Real nrchi_alpha;
	get_bbdep_nrchi_bin( nrchi, nrchi_bin, nrchi_bin_next, nrchi_alpha );

	utility::fixedsizearray1< Real, N > bbs = parent::get_bbs_from_rsd( rsd );
	for ( Size bbi = 1; bbi <= N; ++bbi ) bbs[ bbi ] *= d_multiplier;

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );
	for ( Size bbi = 1; bbi <= N; ++bbi ) {
		debug_assert( bb_bin[ bbi ] >= 1 && bb_bin[ bbi ] <= parent::N_PHIPSI_BINS[bbi] );
	}

	utility::fixedsizearray1< Size, ( 1 << N ) > interp_indices;
	for ( Size ii = 1; ii <= ( 1 << N ); ++ii ) {
		interp_indices[ ii ] = make_conditional_index( N, parent::N_PHIPSI_BINS, ii, bb_bin_next, bb_bin );
	}

	utility::vector1< BBDepScoreInterpData< N > > bbdep_interp_data( ( 1 << ( N + 1 ) ), BBDepScoreInterpData<N>() );

	for ( Size dati = 1; dati <= ( 1 << ( N + 1 ) ); ++dati ) {
		Size i = ( dati + 1 ) / 2;
		// this seems backwards to me too but it is 100% correct
		bbdep_interp_data[ dati ] = bbdep_nrc_interpdata_[ packed_rotno ]( interp_indices[ i ],
			( ( dati % 2 ) ? nrchi_bin : nrchi_bin_next ) );
	}

	Real interpolated_energy( 0.0 );
	//amw test consistency of new formulation
	utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << ( N + 1 ) ) >, ( 1 << ( N + 1 ) ) > n_derivs;
	for ( Size deriv_i = 1; deriv_i <= ( 1 << ( N + 1 ) ); ++deriv_i ) {
		for ( Size dati = 1; dati <= ( 1 << ( N + 1 ) ); ++dati ) {
			n_derivs[ deriv_i ][ dati ] = bbdep_interp_data[ dati ].n_derivs_[ deriv_i ];
		}
	}

	utility::fixedsizearray1< Real, N+1 > dbbp;
	utility::fixedsizearray1< Real, N+1 > binwbb;
	utility::fixedsizearray1< Real, N+1 > dvaldbb;
	for ( Size i = 1; i <= N; ++i ) {
		dbbp[ i ] = bb_alpha[ i ];
		binwbb[ i ] = 10;
	}
	dbbp[ N + 1 ]   = nrchi_alpha;
	binwbb[ N + 1 ] = bbdep_nrchi_binsize_;

	alternate_tricubic_interpolation( n_derivs, dbbp, binwbb, interpolated_energy, dvaldbb );
	for ( Size bbi = 1; bbi <= N; ++bbi ) dnrchi_score_dbb[ bbi ] = dvaldbb[ bbi ];
	dnrchi_score_dnrchi = dvaldbb[ N + 1 ];

	return interpolated_energy;
}

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::best_rotamer_energy(
	conformation::Residue const & rsd,
	bool curr_rotamer_only,
	RotamerLibraryScratchSpace & scratch
) const
{
	debug_assert( rsd.nchi() == T+1 );
	Real nrchi_score( 0 );
	if ( curr_rotamer_only ) {
		Real dnrchiscore_dchi;
		utility::fixedsizearray1< Real, N > dnrchiscore_dbb;
		parent::eval_rotameric_energy_deriv( rsd, scratch, false);

		nrchi_score = bbdep_nrchi_score( rsd, scratch, dnrchiscore_dchi, dnrchiscore_dbb );
		core::conformation::Residue rsd_copy ( rsd );
		utility::vector1< Real > rsd_chi = rsd.chi();

		for ( Size jj = 0; jj <= bbdep_nrchi_nbins_; ++jj ) {
			rsd_chi[ rsd_copy.nchi() ] = nrchi_lower_angle_ + bbdep_nrchi_binsize_ * jj;
			rsd_copy.chi( rsd_chi );
			parent::eval_rotameric_energy_deriv( rsd_copy, scratch, false);
			Real tmp_nrchi_score = bbdep_nrchi_score( rsd_copy, scratch, dnrchiscore_dchi, dnrchiscore_dbb );
			if ( tmp_nrchi_score < nrchi_score ) nrchi_score = tmp_nrchi_score;
		}
	} else {

		// AMW
		utility::fixedsizearray1< Real, N > const bbs2( parent::get_bbs_from_rsd( rsd ) );
		utility::fixedsizearray1< Real, 5 > bbs( 0 );
		for ( Size i = 1; i <= N; ++i ) bbs[ i ] = bbs2[ i ];

		utility::vector1< DunbrackRotamerSampleData > rotamer_samples=/*dunlib->*/get_all_rotamer_samples( bbs );
		//this could be smarter since the T+1 position of the sc_torsions are not used

		Real dnrchiscore_dchi;
		utility::fixedsizearray1< Real, N > dnrchiscore_dbb;
		parent::eval_rotameric_energy_deriv( rsd, scratch, false);
		nrchi_score = bbdep_nrchi_score( rsd, scratch, dnrchiscore_dchi, dnrchiscore_dbb );
		Real tmp_nrchi_score;
		//search the space of terminal chiT
		core::conformation::Residue rsd_copy ( rsd );
		utility::vector1< Real > rsd_chi = rsd.chi();

		for ( Size jj = 1; jj <= rotamer_samples.size(); ++jj ) {
			for ( Size ii = 1; ii <= T; ++ii ) rsd_chi[ ii ] = rotamer_samples[ jj ].chi_mean()[ ii ];

			for ( Size kk = 0; kk <= bbdep_nrchi_nbins_; ++kk ) {
				rsd_chi[ rsd_copy.nchi() ] = nrchi_lower_angle_ + bbdep_nrchi_binsize_ * kk;
				rsd_copy.chi( rsd_chi );
				parent::eval_rotameric_energy_deriv( rsd_copy, scratch, false);
				tmp_nrchi_score = bbdep_nrchi_score( rsd_copy, scratch, dnrchiscore_dchi, dnrchiscore_dbb );
				if ( tmp_nrchi_score < nrchi_score ) nrchi_score = tmp_nrchi_score;
			}
		}
	}
	return nrchi_score;
}


template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::assign_random_rotamer_with_bias(
	conformation::Residue const & rsd,
	pose::Pose const & /*pose*/,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	assign_random_rotamer_with_bias_bbdep( rsd, scratch, RG, new_chi_angles, perturb_from_rotamer_center );
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::assign_random_rotamer_with_bias_bbdep(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	Real random_prob = RG.uniform();

	utility::fixedsizearray1< Real, N > bbs = parent::get_bbs_from_rsd( rsd );
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	BBDepNRChiSample< Real > interpolated_nrchi_sample;
	Size count( 0 );
	while ( random_prob > 0 ) {

		Size index00 = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

		BBDepNRChiSample<> const & nrchi_sample_00( bbdep_rotamers_to_sample_( index00, ++count ) );

		interpolated_nrchi_sample = interpolate_bbdep_nrchi_sample(
			nrchi_sample_00.packed_rotno_, nrchi_sample_00.nrchi_bin_,
			bb_bin, bb_bin_next, bb_alpha );

		random_prob -= interpolated_nrchi_sample.prob_;
		if ( count == bbdep_rotamers_to_sample_.size2() ) break;
	}

	/// Get chimean and chisdves for rotameric chi
	Size packed_rotno = interpolated_nrchi_sample.packed_rotno_;
	PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
	parent::interpolate_rotamers( scratch, packed_rotno, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );

	this->assign_chi_for_interpolated_rotamer( interpolated_rotamer, rsd, RG, new_chi_angles, perturb_from_rotamer_center );

	if ( ! perturb_from_rotamer_center ) {
		new_chi_angles[ T + 1 ] = interpolated_nrchi_sample.nrchi_mean_;
	} else {
		new_chi_angles[ T + 1 ] = interpolated_nrchi_sample.nrchi_mean_ +
			RG.gaussian() * interpolated_nrchi_sample.nrchi_sd_;
	}
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::interpolate_nrchi_values(
	utility::fixedsizearray1< Size, N > bb_bin,
	utility::fixedsizearray1< Size, N > bb_bin_next,
	utility::fixedsizearray1< Real, N > bb_alpha,
	Size packed_rotno,
	utility::vector1< Real > & interpolated_nrchi_distribution
) const {
	// collect 2^N values
	utility::fixedsizearray1< Size, ( 1 << N ) > indices;
	for ( Size ii = 1; ii <= ( 1 << N ); ++ii ) {
		indices[ ii ] = make_conditional_index( N, parent::N_PHIPSI_BINS, ii, bb_bin, bb_bin_next );
	}

	utility::fixedsizearray1< Real, ( 1 << N ) > vals;
	utility::vector1< Real > energies( 0, bbdep_nrchi_nbins_ );
	utility::fixedsizearray1< Real, N > dummy;
	for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {
		for ( Size ii = 1; ii <= ( 1 << N ); ++ii ) {
			utility::fixedsizearray1< Size, ( 1 << N ) > vals;
			vals[ ii ] = bbdep_nrc_interpdata_[ packed_rotno ]( indices[ ii ], jj ).n_derivs_[ 1 ];
		}
		interpolate_polylinear_by_value( vals, bb_alpha, parent::PHIPSI_BINRANGE, false, energies[ jj ], dummy );
	}

	for ( Size ii = 1; ii <= ( 1 << N ); ++ii ) {
		interpolated_nrchi_distribution[ ii ] = std::exp( -energies[ ii ] );
	}
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	rotamers::RotamerVector & rotamers
) const
{
	fill_rotamer_vector_bbdep( pose, scorefxn, task, packer_neighbor_graph,
		concrete_residue, existing_residue, extra_chi_steps, buried, rotamers );
}


template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::fill_rotamer_vector_bbdep(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	rotamers::RotamerVector & rotamers
) const
{
	RotamerLibraryScratchSpace scratch;
	utility::vector1< Size > rotamer_has_been_interpolated( grandparent::n_packed_rots(), 0 );
	utility::vector1< PackedDunbrackRotamer< T, N, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	//Determine whether this is a D-amino acid:
	core::Real d_multiplier = 1.0;
	if ( existing_residue.has_property( "D_AA" ) ) d_multiplier = -1.0;

	/// Save backbone interpolation data for reuse
	utility::fixedsizearray1< Real, N > bbs = parent::get_bbs_from_rsd( existing_residue );
	for ( Size bbi = 1; bbi <= N; ++bbi ) bbs[ bbi ] *= d_multiplier;
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Real const requisit_probability = buried ? 0.95 : 0.87;
	//grandparent::probability_to_accumulate_while_building_rotamers( buried ); -- 98/95 split generates too many samples
	Real accumulated_probability( 0.0 );

	Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	Size count_rotamers_built = 0;
	while ( accumulated_probability < requisit_probability ) {
		++count_rotamers_built;
		Size index00 = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

		utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > nrchi_sample( ( BBDepNRChiSample<>() ) );

		nrchi_sample[ 1 ] = bbdep_rotamers_to_sample_( index00, count_rotamers_built );
		Size const packed_rotno00 = nrchi_sample[ 1 ].packed_rotno_;
		Size const nrchi_bin = nrchi_sample[ 1 ].nrchi_bin_;

		utility::fixedsizearray1< Size, ( 1 << N ) > ind;
		for ( Size indi = 2; indi <= ( 1 << N ); ++indi ) {
			Size index = make_conditional_index( N, parent::N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
			ind[ indi ] = bbdep_rotsample_sorted_order_( index, packed_rotno00, nrchi_bin );
			nrchi_sample[ indi ] = bbdep_rotamers_to_sample_( index, ind[ indi ] );
		}

		BBDepNRChiSample< Real > const interpolated_nrchi_sample = interpolate_bbdep_nrchi_sample( nrchi_sample, bb_alpha );

		if ( rotamer_has_been_interpolated[ packed_rotno00 ] == 0 ) {
			/// interpolate the rotameric chi at most once
			rotamer_has_been_interpolated[ packed_rotno00 ] = 1;
			parent::interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamers[ packed_rotno00 ] );
		}

		build_bbdep_rotamers(
			pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried,
			interpolated_rotamers[ packed_rotno00 ], interpolated_nrchi_sample,
			rotamers );

		accumulated_probability += interpolated_nrchi_sample.prob_;

		if ( count_rotamers_built == max_rots_that_can_be_built ) break;
	}
}

template < Size T, Size N >
utility::vector1< DunbrackRotamerSampleData >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_all_rotamer_samples(
	utility::fixedsizearray1< Real, 5 > bbs
) const
{
	return get_all_rotamer_samples_bbdep( bbs );
}

template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_probability_for_rotamer(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	return get_probability_for_rotamer_bbdep( bbs, rot_ind );
}

template < Size T, Size N >
DunbrackRotamerSampleData
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	return get_rotamer_bbdep( bbs, rot_ind );
}

template < Size T, Size N >
Size
SemiRotamericSingleResidueDunbrackLibrary< T, N >::nchi() const
{
	return T + 1;
}

template < Size T, Size N >
Size
SemiRotamericSingleResidueDunbrackLibrary< T, N >::nbb() const
{
	return N;
}

template < Size T, Size N >
Size
SemiRotamericSingleResidueDunbrackLibrary< T, N >::n_rotamer_bins() const
{
	return grandparent::n_possible_rots() * n_nrchi_sample_bins_;
}

template < Size T, Size N >
utility::vector1< DunbrackRotamerSampleData >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_all_rotamer_samples_bbdep(
	utility::fixedsizearray1< Real, 5 > bbs2
) const
{
	RotamerLibraryScratchSpace scratch;
	utility::vector1< Size > rotamer_has_been_interpolated( grandparent::n_packed_rots(), 0 );
	utility::vector1< PackedDunbrackRotamer< T, N, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	utility::fixedsizearray1< Real, N > bbs;
	for ( Size i = 1 ; i <= N; ++i ) bbs[ i ] = bbs2[ i ];

	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Size const n_rots = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	for ( Size ii = 1; ii <= n_rots; ++ii ) {

		Size index00 = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

		utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > nrchi_sample( ( BBDepNRChiSample<>() ) );
		nrchi_sample[ 1 ] = bbdep_rotamers_to_sample_( index00, ii );
		Size const packed_rotno00 = nrchi_sample[1].packed_rotno_;
		Size const nrchi_bin = nrchi_sample[1].nrchi_bin_;

		utility::fixedsizearray1< Size, ( 1 << N ) > ind;
		for ( Size indi = 2; indi <= ( 1 << N ); ++indi ) {
			Size index = make_conditional_index( N, parent::N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
			ind[ indi ] = bbdep_rotsample_sorted_order_( index, packed_rotno00, nrchi_bin );
			nrchi_sample[ indi ] = bbdep_rotamers_to_sample_( index, ind[ indi ] );
		}

		BBDepNRChiSample< Real > const interpolated_nrchi_sample = interpolate_bbdep_nrchi_sample( nrchi_sample, bb_alpha );

		if ( rotamer_has_been_interpolated[ packed_rotno00 ] == 0 ) {
			/// interpolate the rotameric chi at most once
			rotamer_has_been_interpolated[ packed_rotno00 ] = 1;
			parent::interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamers[ packed_rotno00 ] );
		}

		PackedDunbrackRotamer< T, N, Real > nextrot( interpolated_rotamers[ packed_rotno00 ] );
		DunbrackRotamerSampleData sample( true );
		sample.set_nchi( T + 1 );
		sample.set_rotwell( parent::packed_rotno_2_rotwell( packed_rotno00 ));
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, nextrot.chi_mean( jj ) );
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, nextrot.chi_sd( jj ) );
		sample.set_rotwell( T + 1, nrchi_bin );
		sample.set_chi_mean( T + 1, interpolated_nrchi_sample.nrchi_mean_ );
		sample.set_chi_sd( T + 1, interpolated_nrchi_sample.nrchi_sd_ ); // bogus value
		sample.set_prob( interpolated_nrchi_sample.prob_ );
		sample.set_nrchi_lower_boundary( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno00 ).left_ );
		sample.set_nrchi_upper_boundary( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno00 ).right_ );
		// trick: the interpolated rotamer has the sum of the probabilities for each of the nrchi rotamers
		// so that the probability of getting this nrchi bin given that you're in this bin for the rotameric
		// residues is the quotient of the unconditioned probability for this rotamer and the unconditioned
		// probability of the rotameric chi being in this bin.
		sample.set_nrchi_probability( interpolated_nrchi_sample.prob_ / nextrot.rotamer_probability() );

		all_rots.push_back( sample );
	}
	return all_rots;
}

template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_probability_for_rotamer_bbdep(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Size index00 = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
	utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > nrchi_samples( ( BBDepNRChiSample<>() ) );

	BBDepNRChiSample<> const & nrchi_sample_00( bbdep_rotamers_to_sample_( index00, rot_ind ) );
	nrchi_samples[ 1 ] = nrchi_sample_00;

	Size const packed_rotno00 = nrchi_sample_00.packed_rotno_;
	Size const nrchi_bin = nrchi_sample_00.nrchi_bin_;

	utility::fixedsizearray1< Size, ( 1 << N ) > ind;
	ind[ 1 ] = index00;
	utility::fixedsizearray1< Real, ( 1 << N ) > nrchi_sample_probs;
	nrchi_sample_probs[ 1 ] = static_cast< Real >( nrchi_samples[ 1 ].prob_ );
	for ( Size indi = 2; indi <= ( 1 << N ); ++indi ) {
		Size index = make_conditional_index( N, parent::N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
		ind[ indi ] = bbdep_rotsample_sorted_order_( index, packed_rotno00, nrchi_bin );
		nrchi_samples[ indi ] = bbdep_rotamers_to_sample_( index, ind[ indi ] );
		nrchi_sample_probs[ indi ] = static_cast< Real >( nrchi_samples[ indi ].prob_ );
	}

	Real nrchi_prob;
	utility::fixedsizearray1< Real, N > dummy_dprob;
	utility::fixedsizearray1< Real, N > binw( parent::PHIPSI_BINRANGE );
	interpolate_polylinear_by_value( nrchi_sample_probs, bb_alpha, binw, false /*treat_as_angles*/, nrchi_prob, dummy_dprob );

	return nrchi_prob;
}

template < Size T, Size N >
DunbrackRotamerSampleData
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_bbdep(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	RotamerLibraryScratchSpace scratch;
	Size const n_packed_rots = 1 << N;

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	parent::get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > nrchi_sample( ( BBDepNRChiSample<>() ) );
	Size index00 = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
	nrchi_sample[ 1 ] = bbdep_rotamers_to_sample_( index00, rot_ind );
	Size const packed_rotno00 = nrchi_sample[1].packed_rotno_;
	Size const nrchi_bin = nrchi_sample[1].nrchi_bin_;

	utility::fixedsizearray1< Size, ( 1 << N ) > ind;
	for ( Size indi = 2; indi <= n_packed_rots; ++indi ) {
		Size index = make_conditional_index( N, parent::N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
		ind[ indi ] = bbdep_rotsample_sorted_order_( index, packed_rotno00, nrchi_bin );
		nrchi_sample[ indi ] = bbdep_rotamers_to_sample_( index, ind[ indi ] );
	}

	BBDepNRChiSample< Real > const interpolated_nrchi_sample = interpolate_bbdep_nrchi_sample( nrchi_sample, bb_alpha );

	PackedDunbrackRotamer< T, N, Real > rot;
	parent::interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, rot );

	DunbrackRotamerSampleData sample( true );
	sample.set_nchi( T + 1 );
	sample.set_rotwell( parent::packed_rotno_2_rotwell( packed_rotno00 ));
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, rot.chi_mean( jj ) );
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, rot.chi_sd( jj ) );
	sample.set_rotwell( T + 1, nrchi_bin );
	sample.set_chi_mean( T + 1, interpolated_nrchi_sample.nrchi_mean_ );
	sample.set_chi_sd( T + 1, interpolated_nrchi_sample.nrchi_sd_ ); // bogus value
	sample.set_prob( interpolated_nrchi_sample.prob_ );
	sample.set_nrchi_lower_boundary( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno00 ).left_ );
	sample.set_nrchi_upper_boundary( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno00 ).right_ );
	// trick: the interpolated rotamer has the sum of the probabilities for each of the nrchi rotamers
	// so that the probability of getting this nrchi bin given that you're in this bin for the rotameric
	// residues is the quotient of the unconditioned probability for this rotamer and the unconditioned
	// probability of the rotameric chi being in this bin.
	sample.set_nrchi_probability( interpolated_nrchi_sample.prob_ / rot.rotamer_probability() );

	return sample;
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::build_bbdep_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer,
	BBDepNRChiSample< Real > const interpolated_sample,
	rotamers::RotamerVector & rotamers
) const
{
	DunbrackRotamer< T, N, Real > interpolated_rot( parent::packed_rotamer_2_regular_rotamer( interpolated_rotamer ));
	interpolated_rot.rotamer_probability() = interpolated_rotamer.rotamer_probability();
	BBDepSemiRotamericData< T, N > bbdep_rotamer_building_data( interpolated_rot, interpolated_sample );

	// now build the chi sets derived from this base rotamer
	utility::vector1< ChiSetOP > chi_set_vector;
	parent::enumerate_chi_sets( *concrete_residue, task, existing_residue.seqpos(), buried, bbdep_rotamer_building_data, extra_chi_steps, chi_set_vector );

	parent::create_rotamers_from_chisets( pose, scorefxn, task, packer_neighbor_graph, concrete_residue, existing_residue, chi_set_vector, rotamers );
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::chisamples_for_rotamer_and_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T, N > const & rotamer_data,
	utility::vector1< Real > const & extra_steps,
	utility::vector1< Real > & total_chi,
	utility::vector1< int  > & total_rot,
	utility::vector1< Real > & total_ex_steps,
	utility::vector1< Real > & rotsample_prob
) const
{
	if ( chi_index != T + 1 ) {
		parent::chisamples_for_rotamer_and_chi( rsd_type, rtask, buried, chi_index, rotamer_data, extra_steps,
			total_chi, total_rot, total_ex_steps, rotsample_prob );
	} else {
		bbdep_chisamples_for_rotamer_chi( rsd_type, rtask, buried, chi_index, rotamer_data, extra_steps,
			total_chi, total_rot, total_ex_steps, rotsample_prob );
	}
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::bbdep_chisamples_for_rotamer_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T, N > const & rotamer_data,
	utility::vector1< Real > const & extra_steps,
	utility::vector1< Real > & total_chi,
	utility::vector1< int  > & total_rot,
	utility::vector1< Real > & total_ex_steps,
	utility::vector1< Real > & chisample_prob
) const
{
	using namespace pack::task;
	debug_assert( chi_index == T + 1 );

	BBDepSemiRotamericData< T, N > const & bbdep_rotameric_data(
		static_cast< BBDepSemiRotamericData< T, N > const & > ( rotamer_data ) );

	BBDepNRChiSample< Real > const & nrsample = bbdep_rotameric_data.bbdep_nrchi_sample();

	total_chi.push_back( nrsample.nrchi_mean_ );
	total_ex_steps.push_back( 0.0 );
	total_rot.push_back( nrsample.nrchi_bin_ );
	chisample_prob.push_back( nrsample.prob_ );

	ExtraRotSample ex_samp_level = rtask.extrachi_sample_level( buried, chi_index, rsd_type );

	if ( ex_samp_level == NO_EXTRA_CHI_SAMPLES ) return;

	Real const min_extrachi_sd = 0.0; // What's an appropriate value here?
	if ( nrsample.nrchi_sd_ <= min_extrachi_sd ) return;

	for ( Size k=1; k<= extra_steps.size(); ++k ) {
		total_chi.push_back( nrsample.nrchi_mean_  + extra_steps[ k ] * nrsample.nrchi_sd_  );
		total_ex_steps.push_back( extra_steps[ k ] );
		total_rot.push_back( nrsample.nrchi_bin_ );
		chisample_prob.push_back( nrsample.prob_  ); // no "penalty" for off-mean samples.
	}
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	utility_exit_with_message("Unimplemented!");
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::write_to_binary( utility::io::ozstream & out ) const
{
	using namespace boost;
	parent::write_to_binary( out );
	// For safety, write out const data so that when reading in const data later, we can verify we
	// have a valid input.

	// 1. bbind_nrchi_scoring_
	{
		boost::int32_t bbind_nrchi_scoring = bbind_nrchi_scoring_;
		out.write( ( char* ) &bbind_nrchi_scoring, sizeof( boost::int32_t ) );
	}
	// 2. bbind_nrchi_sampling_
	{
		boost::int32_t bbind_nrchi_sampling = bbind_nrchi_sampling_;
		out.write( ( char* ) &bbind_nrchi_sampling, sizeof( boost::int32_t ) );
	}

	// 3a. bbdep_non_rotameric_chi_scores_
	Size const n_bbdep_scores = grandparent::n_packed_rots() * product( parent::N_PHIPSI_BINS ) * bbdep_nrchi_nbins_;
	BBDepScoreInterpData<N> * bbdep_nrc_interpdata = new BBDepScoreInterpData<N>[ n_bbdep_scores ];

	Size count( 0 );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				bbdep_nrc_interpdata[ count ] = bbdep_nrc_interpdata_[ ii ]( make_index( N, parent::N_PHIPSI_BINS, bb_bin ), jj );
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {//-1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}
	}
	out.write( (char*) bbdep_nrc_interpdata, n_bbdep_scores * sizeof( BBDepScoreInterpData<N> ) );
	delete [] bbdep_nrc_interpdata;

	/// 4.n_nrchi_sample_bins_
	{
		boost::int32_t n_nrchi_sample_bins =  n_nrchi_sample_bins_;
		out.write( ( char* ) & n_nrchi_sample_bins, sizeof( boost::int32_t ) );
	}

	// 5. bbdep_rotamers_to_sample_, bbdep_rotsample_sorted_order_
	Size const n_bbdep_nrchi_samples =
		grandparent::n_packed_rots() * product( parent::N_PHIPSI_BINS ) * n_nrchi_sample_bins_;

	boost::int32_t * bbdep_nrchi_packed_rotnos = new boost::int32_t[ n_bbdep_nrchi_samples ];
	boost::int32_t * bbdep_nrchi_bin     = new boost::int32_t[ n_bbdep_nrchi_samples ];
	DunbrackReal   * bbdep_nrchi_means     = new DunbrackReal[   n_bbdep_nrchi_samples ];
	DunbrackReal   * bbdep_nrchi_sdevs     = new DunbrackReal[   n_bbdep_nrchi_samples ];
	DunbrackReal   * bbdep_nrchi_probs     = new DunbrackReal[   n_bbdep_nrchi_samples ];

	boost::int32_t * bbdep_rotsample_sorted_order = new boost::int32_t[ n_bbdep_nrchi_samples ];

	count = 0;
	Size count_iijj( 1 ); // 3d array sorted by frequency instead of 4d array divided first by rotameric rotno.
	for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
		for ( Size jj = 1; jj <= grandparent::n_packed_rots(); ++jj ) {

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size llkk = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

				bbdep_nrchi_packed_rotnos[    count ] = bbdep_rotamers_to_sample_( llkk, count_iijj ).packed_rotno_;
				bbdep_nrchi_bin[              count ] = bbdep_rotamers_to_sample_( llkk, count_iijj ).nrchi_bin_;
				bbdep_nrchi_means[            count ] = bbdep_rotamers_to_sample_( llkk, count_iijj ).nrchi_mean_;
				bbdep_nrchi_sdevs[            count ] = bbdep_rotamers_to_sample_( llkk, count_iijj ).nrchi_sd_;
				bbdep_nrchi_probs[            count ] = bbdep_rotamers_to_sample_( llkk, count_iijj ).prob_;
				bbdep_rotsample_sorted_order[ count ] = bbdep_rotsample_sorted_order_( llkk, jj, ii );
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[ p ] + 1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[ p ] + 1 ) p = 1;
				}
			}
			++count_iijj;
		}
	}
	out.write( (char*) bbdep_nrchi_packed_rotnos, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
	out.write( (char*) bbdep_nrchi_bin, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
	out.write( (char*) bbdep_nrchi_means, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbdep_nrchi_sdevs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbdep_nrchi_probs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbdep_rotsample_sorted_order, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );

	delete [] bbdep_nrchi_packed_rotnos; bbdep_nrchi_packed_rotnos = 0;
	delete [] bbdep_nrchi_bin;       bbdep_nrchi_bin = 0;
	delete [] bbdep_nrchi_means;   bbdep_nrchi_means = 0;
	delete [] bbdep_nrchi_sdevs;   bbdep_nrchi_sdevs = 0;
	delete [] bbdep_nrchi_probs;   bbdep_nrchi_probs = 0;
	delete [] bbdep_rotsample_sorted_order; bbdep_rotsample_sorted_order = 0;

	{ // Save the bbind rotamer definitions reguardless of bbind_nrchi_sampling_'s status.

		// 6. bbind_rotamers_to_sample_, bbind_rotamers_sorted_by_probability_
		Size const n_bbind_nrchi_samples = n_nrchi_sample_bins_ * grandparent::n_packed_rots();

		DunbrackReal * bbind_nrchi_lefts   = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_medians = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_rights  = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_probs   = new DunbrackReal[ n_bbind_nrchi_samples ];

		boost::int32_t * bbind_rotsample_sorted_order( 0 );

		count = 0;
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
				bbind_nrchi_lefts[   count ] = bbind_rotamers_to_sample_( jj, ii ).left_;
				bbind_nrchi_medians[ count ] = bbind_rotamers_to_sample_( jj, ii ).median_;
				bbind_nrchi_rights[  count ] = bbind_rotamers_to_sample_( jj, ii ).right_;
				bbind_nrchi_probs[   count ] = bbind_rotamers_to_sample_( jj, ii ).prob_;
				++count;
			}
		}
		out.write( (char*) bbind_nrchi_lefts,   n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		out.write( (char*) bbind_nrchi_medians, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		out.write( (char*) bbind_nrchi_rights,  n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		out.write( (char*) bbind_nrchi_probs,   n_bbind_nrchi_samples * sizeof( DunbrackReal ) );

		delete [] bbind_nrchi_lefts;
		delete [] bbind_nrchi_medians;
		delete [] bbind_nrchi_rights;
		delete [] bbind_nrchi_probs;
		delete [] bbind_rotsample_sorted_order;
	}

}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::read_from_binary( utility::io::izstream & in )
{
	using namespace boost;
	parent::read_from_binary( in );
	Size const num_rot_bin = product( parent::N_PHIPSI_BINS );
	// Check that constant member data initialized in the ctor is consistent with this input file
	// 1. bbind_nrchi_scoring_
	{
		boost::int32_t bbind_nrchi_scoring( 0 );
		in.read( (char*) &bbind_nrchi_scoring, sizeof( boost::int32_t ));
		if ( bbind_nrchi_scoring_ != static_cast< bool > ( bbind_nrchi_scoring ) ) {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone dependent non-rotameric chi scoring,\nbut the instantiated Semirotamer library"
				<< " was created with backbone independent non-rotameric chi scoring.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
			utility_exit();
		}
	}

	// 2. bbind_nrchi_sampling_
	{
		boost::int32_t bbind_nrchi_sampling( 0 );
		in.read( (char*) &bbind_nrchi_sampling, sizeof( boost::int32_t ));
		if ( bbind_nrchi_sampling_ != static_cast< bool > (  bbind_nrchi_sampling ) ) {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone dependent non-rotameric chi sampling,\nbut the instantiated Semirotamer library"
				<< " was created with backbone independent non-rotameric chi sampling.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
			utility_exit();
		}
	}

	// 3a. bbdep_non_rotameric_chi_scores_
	Size const n_bbdep_scores = grandparent::n_packed_rots() * num_rot_bin * bbdep_nrchi_nbins_;
	BBDepScoreInterpData<N> * bbdep_nrc_interpdata = new BBDepScoreInterpData<N>[ n_bbdep_scores ];
	in.read( (char*) bbdep_nrc_interpdata, n_bbdep_scores * sizeof( BBDepScoreInterpData<N> ) );

	Size count( 0 );
	bbdep_nrc_interpdata_.resize( grandparent::n_packed_rots() );//, far2d );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		bbdep_nrc_interpdata_[ ii ].dimension( num_rot_bin, bbdep_nrchi_nbins_, BBDepScoreInterpData<N>() );
		for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {
			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size llkk = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
				bbdep_nrc_interpdata_[ ii ]( llkk, jj ) = bbdep_nrc_interpdata[ count ];//BBDepScoreInterpData<N>();
				bbdep_nrc_interpdata_[ ii ]( llkk, jj ).n_derivs_ = bbdep_nrc_interpdata[ count ].n_derivs_;
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[ p ] + 1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[ p ] + 1 ) p = 1;
				}
			}
		}
	}
	delete [] bbdep_nrc_interpdata;

	/// 4.n_nrchi_sample_bins_
	{
		boost::int32_t n_nrchi_sample_bins( 0 );
		in.read( (char*) & n_nrchi_sample_bins, sizeof( boost::int32_t ) );
		n_nrchi_sample_bins_ = n_nrchi_sample_bins;
	}

	// 5. bbdep_rotamers_to_sample_, bbdep_rotsample_sorted_order_
	Size const n_bbdep_nrchi_samples =
		grandparent::n_packed_rots() * num_rot_bin * n_nrchi_sample_bins_;

	boost::int32_t * bbdep_nrchi_packed_rotnos = new boost::int32_t[ n_bbdep_nrchi_samples ];
	boost::int32_t * bbdep_nrchi_bin     = new boost::int32_t[ n_bbdep_nrchi_samples ];
	DunbrackReal *   bbdep_nrchi_means     = new DunbrackReal[   n_bbdep_nrchi_samples ];
	DunbrackReal *   bbdep_nrchi_sdevs     = new DunbrackReal[   n_bbdep_nrchi_samples ];
	DunbrackReal *   bbdep_nrchi_probs     = new DunbrackReal[   n_bbdep_nrchi_samples ];

	boost::int32_t * bbdep_rotsample_sorted_order = new boost::int32_t[ n_bbdep_nrchi_samples ];

	in.read( (char*) bbdep_nrchi_packed_rotnos, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
	in.read( (char*) bbdep_nrchi_bin, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
	in.read( (char*) bbdep_nrchi_means, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbdep_nrchi_sdevs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbdep_nrchi_probs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbdep_rotsample_sorted_order, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );

	bbdep_rotamers_to_sample_.dimension( num_rot_bin, grandparent::n_packed_rots()*n_nrchi_sample_bins_ );
	bbdep_rotsample_sorted_order_.dimension( num_rot_bin, grandparent::n_packed_rots(), n_nrchi_sample_bins_ );

	count = 0;
	Size count_iijj( 1 ); // 3d array sorted by frequency instead of 4d array divided first by rotameric rotno.
	for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
		for ( Size jj = 1; jj <= grandparent::n_packed_rots(); ++jj ) {
			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size rot_index = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

				bbdep_rotamers_to_sample_( rot_index, count_iijj ).packed_rotno_ = bbdep_nrchi_packed_rotnos[ count ];
				bbdep_rotamers_to_sample_( rot_index, count_iijj ).nrchi_bin_ = bbdep_nrchi_bin[ count ];
				bbdep_rotamers_to_sample_( rot_index, count_iijj ).nrchi_mean_ = bbdep_nrchi_means[ count ];
				bbdep_rotamers_to_sample_( rot_index, count_iijj ).nrchi_sd_ = bbdep_nrchi_sdevs[ count ];
				bbdep_rotamers_to_sample_( rot_index, count_iijj ).prob_ = bbdep_nrchi_probs[ count ];
				bbdep_rotsample_sorted_order_( rot_index, jj, ii ) = bbdep_rotsample_sorted_order[count];
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[ p ] + 1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[ p ] + 1 ) p = 1;
				}
			}
			++count_iijj;
		}
	}

	delete [] bbdep_nrchi_packed_rotnos; bbdep_nrchi_packed_rotnos = 0;
	delete [] bbdep_nrchi_bin;       bbdep_nrchi_bin = 0;
	delete [] bbdep_nrchi_means;   bbdep_nrchi_means = 0;
	delete [] bbdep_nrchi_sdevs;   bbdep_nrchi_sdevs = 0;
	delete [] bbdep_nrchi_probs;   bbdep_nrchi_probs = 0;
	delete [] bbdep_rotsample_sorted_order; bbdep_rotsample_sorted_order = 0;

	{ // Save the bbind rotamer definitions reguardless of bbind_nrchi_sampling_'s status.

		// 6. bbind_rotamers_to_sample_, bbind_rotamers_sorted_by_probability_
		Size const n_bbind_nrchi_samples = n_nrchi_sample_bins_ * grandparent::n_packed_rots();

		DunbrackReal * bbind_nrchi_lefts   = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_medians = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_rights  = new DunbrackReal[ n_bbind_nrchi_samples ];
		DunbrackReal * bbind_nrchi_probs   = new DunbrackReal[ n_bbind_nrchi_samples ];
		boost::int32_t * bbind_rotsample_sorted_order = new boost::int32_t[  n_bbind_nrchi_samples ];

		in.read( (char*) bbind_nrchi_lefts, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbind_nrchi_medians, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbind_nrchi_rights, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbind_nrchi_probs, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );

		bbind_rotamers_to_sample_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );

		Size count( 0 );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
				bbind_rotamers_to_sample_( jj, ii ).left_   = bbind_nrchi_lefts[ count ];
				bbind_rotamers_to_sample_( jj, ii ).median_ = bbind_nrchi_medians[ count ];
				bbind_rotamers_to_sample_( jj, ii ).right_  = bbind_nrchi_rights[  count ];
				bbind_rotamers_to_sample_( jj, ii ).prob_   = bbind_nrchi_probs[   count ];
				++count;
			}
		}

		delete [] bbind_nrchi_lefts;
		delete [] bbind_nrchi_medians;
		delete [] bbind_nrchi_rights;
		delete [] bbind_nrchi_probs;
		delete [] bbind_rotsample_sorted_order;
	}

}

/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
/// Values tested should parallel those used in the read_from_binary() function.
template < Size T, Size N >
bool
SemiRotamericSingleResidueDunbrackLibrary< T, N >::operator ==( rotamers::SingleResidueRotamerLibrary const & rhs) const {
	static thread_local basic::Tracer TR( "core.pack.dunbrack.SemiRotamericSingleResidueDunbrackLibrary" );

	// Raw pointer okay, we're just using it to check for conversion
	SemiRotamericSingleResidueDunbrackLibrary< T, N> const * ptr( dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< T, N > const * > ( &rhs ) );
	if ( ptr == 0 ) {
		TR << "In comparison operator: right-hand side is not a matching SemiRotamericSingleResidueDunbrackLibrary." << std::endl;
		return false;
	}
	SemiRotamericSingleResidueDunbrackLibrary< T, N > const & other( dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< T, N > const & > ( rhs ) );

	bool equal( true );

	if ( ! parent::operator==( rhs ) ) {
		//TR.Debug << "In SemiRotamericSingleResidueDunbrackLibrary< T >::operator== : Parent comparsion returns false, stopping equality test." << std::endl;
		//return false;
		equal = false;
	}

	core::Real const & ANGLE_DELTA( grandparent::ANGLE_DELTA );

	// 0. Non-binary-loaded data
	if ( ! numeric::equal_by_epsilon( nrchi_periodicity_, other.nrchi_periodicity_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " nrchi_periodicity_ " << nrchi_periodicity_ << " vs. " << other.nrchi_periodicity_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( nrchi_lower_angle_, other.nrchi_lower_angle_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " nrchi_lower_angle_ " << nrchi_lower_angle_ << " vs. " << other.nrchi_lower_angle_ << std::endl;
		equal = false;
	}
	if ( bbdep_nrchi_nbins_ != other.bbdep_nrchi_nbins_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbdep_nrchi_nbins_ " << bbdep_nrchi_nbins_ << " vs. " << other.bbdep_nrchi_nbins_ << std::endl;
		return false; // Major data consistency issues - don't bother with the rest
	}
	if ( ! numeric::equal_by_epsilon( bbdep_nrchi_binsize_, other.bbdep_nrchi_binsize_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbdep_nrchi_binsize_ " << bbdep_nrchi_binsize_ << " vs. " << other.bbdep_nrchi_binsize_ << std::endl;
		equal = false;
	}
	if ( bbind_nrchi_nbins_ != other.bbind_nrchi_nbins_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbind_nrchi_nbins_ " << bbind_nrchi_nbins_ << " vs. " << other.bbind_nrchi_nbins_ << std::endl;
		return false; // Major data consistency issues - don't bother with the rest
	}
	if ( ! numeric::equal_by_epsilon( bbind_nrchi_binsize_, other.bbind_nrchi_binsize_, ANGLE_DELTA ) ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbind_nrchi_binsize_ " << bbind_nrchi_binsize_ << " vs. " << other.bbind_nrchi_binsize_ << std::endl;
		equal = false;
	}


	// 1. bbind_nrchi_scoring_
	if ( bbind_nrchi_scoring_ != other.bbind_nrchi_scoring_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbind_nrchi_scoring_ " << bbind_nrchi_scoring_ << " vs. " << other.bbind_nrchi_scoring_ << std::endl;
		return false; // Major data consistency issues - don't bother with the rest
	}

	// 2. bbind_nrchi_sampling_
	if ( bbind_nrchi_sampling_ != other.bbind_nrchi_sampling_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " bbind_nrchi_sampling_ " << bbind_nrchi_sampling_ << " vs. " << other.bbind_nrchi_sampling_ << std::endl;
		return false; // Major data consistency issues - don't bother with the rest
	}

	if ( grandparent::n_packed_rots() != other.n_packed_rots() ) {
		return false;
	}

	// 3a. bbdep_non_rotameric_chi_scores_
	assert( bbdep_nrchi_nbins_ == other.bbdep_nrchi_nbins_ ); // Assumed the same?
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {
			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {
				Size bb_rot_index = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
				if ( !( bbdep_nrc_interpdata_[ ii ]( bb_rot_index, jj ) == other.bbdep_nrc_interpdata_[ ii ]( bb_rot_index, jj ) ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
						<< " bbdep_nrc_interpdata:   " << ii << " " << bb_bin[1] << " " << bb_bin[2] << " " << jj << std::endl;
					equal = false;
				}

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[ p ] + 1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[ p ] + 1 ) p = 1;
				}
			}
		}
	}

	/// 4.n_nrchi_sample_bins_
	if ( n_nrchi_sample_bins_ != other.n_nrchi_sample_bins_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
			<< " n_nrchi_sample_bins_ " << n_nrchi_sample_bins_ << " vs. " << other.n_nrchi_sample_bins_ << std::endl;
		return false; // Major data consistency issues - don't bother with the rest
	}

	// 5. bbdep_rotamers_to_sample_, bbdep_rotsample_sorted_order_
	Size count_iijj( 1 ); // 3d array sorted by frequency instead of 4d array divided first by rotameric rotno.
	for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
		for ( Size jj = 1; jj <= grandparent::n_packed_rots(); ++jj ) {

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size bb_rot_index = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
				if ( !( bbdep_rotamers_to_sample_( bb_rot_index, count_iijj ) == other.bbdep_rotamers_to_sample_( bb_rot_index, count_iijj ) ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
						<< " bbdep_rotamers_to_sample: " << ii << " " << jj << " " << bb_bin[1] << " " << bb_bin[2] << std::endl;
					equal = false;
				}
				if ( bbdep_rotsample_sorted_order_( bb_rot_index, jj, ii ) != other.bbdep_rotsample_sorted_order_( bb_rot_index, jj, ii ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
						<< " bbdep_rotsample_sorted_order " << ii << " " << jj << " " << bb_bin[1] << " " << bb_bin[2] << " - "
						<< bbdep_rotsample_sorted_order_( bb_rot_index, jj, ii ) << " vs. " << other.bbdep_rotsample_sorted_order_( bb_rot_index, jj, ii ) << std::endl;
					equal = false;
				}
				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[ p ] + 1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[ p ] + 1 ) p = 1;
				}
			}
			++count_iijj;
		}
	}

	// 6. bbind_rotamers_to_sample_, bbind_rotamers_sorted_by_probability_
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
			if ( !( bbind_rotamers_to_sample_( jj, ii ) == other.bbind_rotamers_to_sample_( jj, ii ) ) ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
					<< " bbind_rotamers_to_sample: " << ii << " " << jj << std::endl;
				equal = false;
			}
			if ( bbind_nrchi_sampling_ ) {
				if ( bbind_rotamers_sorted_by_probability_( jj, ii ) != other.bbind_rotamers_sorted_by_probability_( jj, ii ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( grandparent::aa() )
						<< " bbind_rotamers_sorted_by_probability: " << ii << " " << jj << " - "
						<< bbind_rotamers_sorted_by_probability_( jj, ii ) << " vs. " << bbind_rotamers_sorted_by_probability_( jj, ii ) << std::endl;
					equal = false;
				}
			}
		}
	}

	return equal;
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi(
	ChiVector const & chi,
	RotVector & rot
) const
{
	parent::get_rotamer_from_chi( chi, rot );
	Size const packed_rotno( grandparent::rotwell_2_packed_rotno( rot ) );
	Real const nrchi = clip_to_nrchi_range( chi[ T + 1 ] );
	/// find bracketting chi's as defined in the bbind rotamer definitions file.
	for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
		Real const left = bbind_rotamers_to_sample_( ii, packed_rotno ).left_;
		Real const right = bbind_rotamers_to_sample_( ii, packed_rotno ).right_;
		if ( left < nrchi_lower_angle_ ) {
			/// nrchi bin spans periodicity gap...
			if ( ! ( left <= nrchi && nrchi <= right ) &&
					! ( left + nrchi_periodicity_ <= nrchi && nrchi <= right + nrchi_periodicity_ ) ) {
				continue;
			}
		} else if ( right > nrchi_lower_angle_ + nrchi_periodicity_ ) {
			if ( ! ( left <= nrchi && nrchi <= right ) &&
					! ( left - nrchi_periodicity_ <= nrchi && nrchi <= right - nrchi_periodicity_ ) ) {
				continue;
			}
		} else if ( left > nrchi ) {
			continue;
		} else if (  nrchi > right ) {
			continue;
		}

		rot[ T + 1 ] = ii;
		break;
	}

	if ( rot[ T + 1 ] == 0 ) {
		std::cerr << "Failed to find bracketing bin for : " << nrchi << T;
		for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
			std::cerr << "( " << ii << " : " << bbind_rotamers_to_sample_( ii, packed_rotno ).left_ << ", " << bbind_rotamers_to_sample_(ii, packed_rotno ).right_ << ") ";
		}
		std::cerr << std::endl;
		utility_exit();
	}
}

/// @brief For the non rotameric chi, how many asymmetric degrees are there?  (e.g. 360 for asn, 180 for asp)
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::set_nrchi_periodicity(
	Real angle_in_degrees
)
{
	nrchi_periodicity_ = angle_in_degrees;
}

/// @brief What angle do the interpolation data start from?
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::set_nonrotameric_chi_start_angle(
	Real angle_in_degrees
)
{
	nrchi_lower_angle_ = angle_in_degrees;
}

/// @brief What is the angular step size of the bbdep score?
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::set_nonrotameric_chi_bbdep_scoring_step_size(
	Real step_size_in_degrees
)
{
	debug_assert( nrchi_periodicity_ != 0.0 );
	debug_assert( step_size_in_degrees != 0.0 );
	bbdep_nrchi_binsize_ = step_size_in_degrees;
	bbdep_nrchi_nbins_ = static_cast< Size > (nrchi_periodicity_ / step_size_in_degrees);
}

/// @brief What is the angular step size of the bbind score?
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::set_nonrotameric_chi_bbind_scoring_step_size(
	Real step_size_in_degrees
)
{
	debug_assert( nrchi_periodicity_ != 0.0 );
	debug_assert( step_size_in_degrees != 0.0 );
	bbind_nrchi_binsize_ = step_size_in_degrees;
	bbind_nrchi_nbins_ = static_cast< Size > (nrchi_periodicity_ / step_size_in_degrees);
}


template < Size T, Size N >
Size
SemiRotamericSingleResidueDunbrackLibrary< T, N >::memory_usage_dynamic() const
{
	Size total_memory = parent::memory_usage_dynamic();

	/// for bbdep nrchi scoring
	for ( Size ii = 1; ii <= bbdep_nrc_interpdata_.size(); ++ii ) {
		total_memory += bbdep_nrc_interpdata_[ ii ].size() * sizeof( BBDepScoreInterpData<N> );
	}
	total_memory += bbdep_nrc_interpdata_.size() * sizeof( ObjexxFCL::FArray3D< BBDepScoreInterpData<N> > );

	/// for bbind nrchi scoring
	total_memory += bbind_non_rotameric_chi_scores_.size() * sizeof( Real );

	/// for bbind nrchi sampling
	total_memory += bbind_rotamers_to_sample_.size() * sizeof( BBIndNRChiSample<> );
	total_memory += bbind_rotamers_sorted_by_probability_.size() * sizeof( Size );

	/// for bbdep nrchi sampling
	total_memory += bbdep_rotamers_to_sample_.size() * sizeof( BBDepNRChiSample<> );
	total_memory += bbdep_rotsample_sorted_order_.size() * sizeof( Size );
	/// parental cost
	total_memory += parent::memory_usage_dynamic();

	return total_memory;
}


template < Size T, Size N >
Size SemiRotamericSingleResidueDunbrackLibrary< T, N >::memory_usage_static() const
{
	return sizeof( SemiRotamericSingleResidueDunbrackLibrary< T, N > );
}

/// @breif read-from-files that does not require the in_continmin_bbind file
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::read_from_files(
	utility::io::izstream & in_rotdef,
	utility::io::izstream & in_rotameric,
	utility::io::izstream & in_continmin_bbdep
)
{
	debug_assert( ! bbind_nrchi_scoring_ ); // make sure you're not actually trying to use backbone-independent scoring.

	read_rotamer_definitions( in_rotdef );
	read_rotameric_data( in_rotameric );
	read_bbdep_continuous_minimization_data( in_continmin_bbdep );
}

/// @details The rotamer definition file has already been read before this
/// function is called, so all of the rotameric rotamer wells have been encountered
/// and the packed_rotno enumeration may be used instead of the regular rotno.
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::read_rotameric_data(
	utility::io::izstream & in_rotameric
)
{
	std::string three_letter_code;
	Size const num_bb_bins = product( parent::N_PHIPSI_BINS );

	utility::fixedsizearray1< Real, N > bbs;
	Size count;
	utility::vector1< Size > rotwell( DUNBRACK_MAX_SCTOR, 0 );
	DunbrackReal prob;
	//MaximCode
	//The natural logarithm is the base-e logarithm: the inverse of the natural exponential function (exp).
	DunbrackReal minusLogProbability(-1e-10);
	utility::vector1< DunbrackReal > mean( DUNBRACK_MAX_SCTOR, 0.0 );
	utility::vector1< DunbrackReal > stdev( DUNBRACK_MAX_SCTOR, 0.0 );

	ObjexxFCL::FArray1D< Size > rotameric_count( num_bb_bins, Size( 0 ) );
	ObjexxFCL::FArray1D< Size > semi_rotameric_count( num_bb_bins, Size( 0 ));

	while ( in_rotameric ) {
		/// 1. peek at the line; if it starts with #, skip to the next line.
		char first_char = in_rotameric.peek();
		if ( first_char == '#' ) {
			std::string line;
			in_rotameric.getline( line );
			continue;
		}

		/// 2. Line contains:
		/// a. AA three-letter code
		/// b. phi (degrees)
		/// c. psi (degrees)
		/// d. observed count
		/// e-h. chi-(1-4) wells (integers)
		/// i. rotamer probability
		/// j-m. chi-(1-4) mean (degrees)
		/// n-o. chi-(1-4) stdev (degrees );

		in_rotameric >> three_letter_code;
		for ( Size ii = 1; ii <= N; ++ii ) in_rotameric >> bbs[ ii ];
		in_rotameric >> count;
		in_rotameric >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		in_rotameric >> prob;
		//MaximCode
		if ( basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
				&& basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable] ) {
			in_rotameric >> minusLogProbability;
		}
		in_rotameric >> mean[ 1 ] >> mean[ 2 ] >> mean[ 3 ] >> mean[ 4 ];
		in_rotameric >> stdev[ 1 ] >> stdev[ 2 ] >> stdev[ 3 ] >> stdev[ 4 ];

		if ( ! in_rotameric ) break; // we've read past the end of the file...

		bool cont = false;
		// skip -180 because some of our libraries omit -180
		for ( Size ii = 1; ii <= N; ++ii ) if ( bbs[ ii ] == -180 ) cont = true; // duplicated data...
		if ( cont ) continue;

		/// AVOID inf and NaN by correcting the database
		for ( Size ii = 1; ii <= T; ++ii ) if ( stdev[ ii ] == 0.0 ) stdev[ ii ] = 5; // bogus!

		if ( prob <= 1e-6 ) {
			prob = 1e-6;
			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
			/// AMW -- Changing this to be a <= criterion and to 1e-6 to actually match the resolution
			/// of the library
		}

		utility::fixedsizearray1< Size, N > bb_bin;
		parent::get_bb_bins( bbs, bb_bin );

		// calc index
		Size ind = make_index( N, parent::N_PHIPSI_BINS, bb_bin );
		++semi_rotameric_count( ind );

		/// Store first T elements in rotamers_, and, if we're going to build bbdep
		/// rotamers during packing, store data for the last chi in the bbdep_rotamers_to_sample_
		/// array.  If this rotamer has already been encountered, increment its probability.
		/// Because the probabilities of the rotameric chi are accumulated over all the
		/// non-rotameric chi samples, but the rotamers do not change order once they're found, they will not be
		/// in guaranteed sorted order at the end of reading this file.  The "sorted" index
		/// is somewhat of a misnomer -- it simply means the sorted order in which the rotamer
		/// was first encountered.  Fortunately, none of the rotamer building subroutines
		/// require that the rotamers_ tables be sorted by probability.

		Size const packed_rotno( grandparent::rotwell_2_packed_rotno( rotwell ) );
		if ( parent::packed_rotno_2_sorted_rotno()( ind, packed_rotno ) == 0 ) {
			// first time this rotwell has been encountered.
			Size nencountered = ++rotameric_count( ind );
			parent::rotamers()( ind, nencountered ).packed_rotno() = packed_rotno;
			parent::packed_rotno_2_sorted_rotno()( ind, packed_rotno ) = nencountered;
			for ( Size ii = 1; ii <= T; ++ii ) {
				parent::rotamers()( ind, nencountered ).chi_mean( ii, mean[  ii ] );
				parent::rotamers()( ind, nencountered ).chi_sd(   ii, stdev[ ii ] );
			}
			parent::rotamers()( ind, nencountered ).rotamer_probability() = prob;
			for ( Size di = 1; di <= ( 1 << N ); ++di ) parent::rotamers()( ind, nencountered ).n_derivs()[ di ] = 0;
		} else {
			Size const sorted_rotno = parent::packed_rotno_2_sorted_rotno()( ind, packed_rotno );
			parent::rotamers()( ind, sorted_rotno ).packed_rotno() = packed_rotno;
			parent::rotamers()( ind, sorted_rotno ).rotamer_probability() += prob;
			// verify that the data file is as expected: chimean and chisd do not change for
			// rotameric chi across different entries of the nrchi.
			for ( Size ii = 1; ii <= T; ++ii ) {
				if ( std::abs( basic::periodic_range( parent::rotamers()( ind, sorted_rotno ).chi_mean( ii ) - mean[  ii ], 360 )) > 1e-5 ) {
					std::cerr << "ERROR: semi-rotameric input invalid -- chi means differ." << std::endl;
					for ( Size bbi = 1; bbi <= N; ++bbi ) std::cerr << "bb" << bbi << ": " << bbs[ bbi ];
					std::cerr << " original chi_" << ii << ": ";
					std::cerr << parent::rotamers()( ind, sorted_rotno ).chi_mean( ii );
					std::cerr << " later chi_" << ii << ": "<< mean[  ii ] << std::endl;
					utility_exit();
				}
				if ( std::abs( parent::rotamers()( ind, sorted_rotno ).chi_sd( ii ) - stdev[ ii ]) > 1e-5 ) {
					std::cerr << "ERROR: semi-rotameric input invalid -- chi stdevs differ." << std::endl;
					for ( Size bbi = 1; bbi <= N; ++bbi ) std::cerr << "bb" << bbi << ": " << bbs[ bbi ];
					std::cerr << " original chi_" << ii << ": ";
					std::cerr << parent::rotamers()( ind, sorted_rotno ).chi_sd( ii );
					std::cerr << " later chi_" << ii << ": "<< stdev[  ii ] << std::endl;
					utility_exit();
				}
			}
		}

		Size const which_rotamer = semi_rotameric_count( ind );
		BBDepNRChiSample<> & sample( bbdep_rotamers_to_sample_( ind, which_rotamer ) );
		sample.packed_rotno_ = packed_rotno;
		sample.nrchi_bin_ = rotwell[ T + 1 ];
		sample.nrchi_mean_ = mean[ T + 1 ];
		sample.nrchi_sd_ = stdev[ T + 1 ];
		sample.prob_ = prob;
		bbdep_rotsample_sorted_order_( ind, packed_rotno, rotwell[ T + 1 ] ) = which_rotamer;
	}
	parent::initialize_bicubic_splines();
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::read_bbdep_continuous_minimization_data(
	utility::io::izstream & in_continmin
)
{
	Size const n_bin_pow = product( parent::N_PHIPSI_BINS );
	if ( bbind_nrchi_scoring_ ) return;

	std::string three_letter_code;
	utility::fixedsizearray1< Real, N > bbs;
	Real base_prob;
	//MaximCode
	//The natural logarithm is the base-e logarithm: the inverse of the natural exponential function (exp).
	Real minusLogProbability(-1e-10);
	Size count;
	utility::vector1< Size > rotwell( T );
	utility::vector1< Real > chimean( T );
	utility::vector1< Real > chisd( T );
	utility::vector1< Real > nrchi_probs( bbdep_nrchi_nbins_, 0.0 );

	utility::vector1< ObjexxFCL::FArray2D< Real > > bbdep_non_rotameric_chi_scores;
	bbdep_non_rotameric_chi_scores.resize( grandparent::n_packed_rots() );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		bbdep_non_rotameric_chi_scores[ ii ].dimension( n_bin_pow, bbdep_nrchi_nbins_, Real(0.0) );
	}

	while ( in_continmin ) {
		char first_char = in_continmin.peek();
		if ( first_char == '#' ) {
			std::string line; in_continmin.getline( line );
			continue;
		}

		// The structure of an input line is:
		// a. Three-letter code
		// b. phi
		// c. psi
		// d. observation count
		// e... rotamer well indices for rotameric residues (1 or 2)
		// f. probability of this well being occupied
		// g... rotameric chi means
		// h... rotameric chi sdevs
		// i... non-rotameric chi probabilities binned -- sum to 1.
		in_continmin >> three_letter_code;
		for ( Size ii = 1; ii <= N; ++ii )                  in_continmin >> bbs[ ii ];
		in_continmin >> count;
		for ( Size ii = 1; ii <= T; ++ii )                  in_continmin >> rotwell[ ii ];
		in_continmin >> base_prob;
		//MaximCode
		if ( basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
				&& basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable] ) {
			in_continmin >> minusLogProbability;
		}
		for ( Size ii = 1; ii <= T; ++ii )                  in_continmin >> chimean[ ii ];
		for ( Size ii = 1; ii <= T; ++ii )                  in_continmin >> chisd[ ii ];
		for ( Size ii = 1; ii <= bbdep_nrchi_nbins_; ++ii ) in_continmin >> nrchi_probs[ ii ];

		if ( ! in_continmin ) break; // we've read past the end of the file...
		bool cont = false;
		for ( Size ii = 1; ii <= N; ++ii ) if ( bbs[ ii ] == -180 ) cont = true; // duplicated data...
		if ( cont ) continue;

		/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
		/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
		/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
		/// AMW -- Changing this to be a <= criterion and to 1e-6 to actually match the resolution
		/// of the library
		if ( base_prob <= 1e-6 ) base_prob = 1e-6; // correct probability 0 events.

		/// Now, convert the probabilities into scores
		utility::fixedsizearray1< Size, N > bb_bin;
		parent::get_bb_bins( bbs, bb_bin );
		Size ind = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

		Size const rotno = grandparent::rotwell_2_rotno( rotwell );
		for ( Size ii = 1; ii <= bbdep_nrchi_nbins_; ++ii ) {
			/// input data has probability 0 events; assign a tiny probability for these cases.

			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
			/// AMW -- Changing this to be a <= criterion and to 1e-6 to actually match the resolution
			/// of the library
			Real const prob = base_prob * ( nrchi_probs[ ii ] <= 1e-6 ? 1e-6 : nrchi_probs[ ii ] );
			bbdep_non_rotameric_chi_scores[ rotno ]( ind, ii ) = -std::log( prob );
		}
	}

	/// Now create the tricubic interpolation data and store that data in the
	using namespace numeric;
	using namespace numeric::interpolation::spline;
	utility::vector1< BorderFlag > border( N, e_Periodic );
	border.push_back( e_Periodic );
	utility::vector1< Real > start( N, -180.0 );
	start.push_back( nrchi_lower_angle_ );
	utility::vector1< Real > delta( N, 10.0 );
	delta.push_back( bbdep_nrchi_binsize_ );
	utility::vector1< bool > lin_cont( N, true );
	lin_cont.push_back( true );
	utility::vector1< std::pair< double, double>  > first_be( N, std::pair< double, double > ( 10, 10 ) );
	first_be.push_back( std::pair< double, double>( 10, 10 ) );

	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {

		utility::vector1< Size > dimensions( N );
		for ( Size i = 1; i <= N; ++i ) dimensions[ i ] = parent::N_PHIPSI_BINS[i];
		dimensions.push_back( bbdep_nrchi_nbins_ );
		MathNTensor< Real > data( dimensions, Real(0.0) );
		utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
		utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
		for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];

		Size p = 1;
		while ( bb_bin[ N + 1 ] == 1 ) {

			utility::vector1< Size > position( N );
			for ( Size qq = 1; qq <= N; ++qq ) position[ qq ] = bb_bin[ qq ];
			position.push_back( 1 );
			utility::vector1< Size > indices;
			for ( Size i = 1; i <= N; ++i ) indices.push_back( bb_bin[ i ] - 1 );

			Size bbdepindex = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

			for ( Size ll = 1; ll <= bbdep_nrchi_nbins_; ++ll ) {
				position[ N + 1 ] = ll;
				data( position ) = bbdep_non_rotameric_chi_scores[ ii ]( bbdepindex, ll );
			}

			bb_bin[ 1 ]++;
			while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
				bb_bin[ p ] = 1;
				bb_bin[ ++p ]++;
				if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
			}
		}

		PolycubicSpline spline( N + 1 ); // extra dimension because of the nrchi
		spline.train( border, start, delta, data, lin_cont, first_be );

		for ( Size bbi = 1; bbi <= ( N+1 ); ++bbi ) bb_bin[ bbi ] = 1;
		for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = parent::N_PHIPSI_BINS[i];
		p = 1;
		while ( bb_bin[ N + 1 ] == 1 ) {

			utility::vector1< Size > position( N );
			for ( Size qq = 1; qq <= N; ++qq ) position[ qq ] = bb_bin[ qq ];
			position.push_back( 1 );
			Size bbdepindex = make_index( N, parent::N_PHIPSI_BINS, bb_bin );

			for ( Size ll = 1; ll <= bbdep_nrchi_nbins_; ++ll ) {
				position[ N + 1 ] = ll;
				bbdep_nrc_interpdata_[ ii ]( bbdepindex, ll ).n_derivs_[ 1 ]   = ( DunbrackReal ) data( position );
				for ( Size di = 2; di <= ( 1 << ( N + 1 ) ); ++di ) {
					bbdep_nrc_interpdata_[ ii ]( bbdepindex, ll ).n_derivs_[ di ] = ( DunbrackReal ) spline.get_deriv( di )( position );
				}
			}

			bb_bin[ 1 ]++;
			while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
				bb_bin[ p ] = 1;
				bb_bin[ ++p ]++;
				if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
			}
		}
	}
}

/// @details the rotamer definition file must be the first file read.
/// Reading this file gives the number of pseudorotamers used for
/// both the bbdep and bbind rotamer building.  The bbind rotamer data
/// is saved iff the bbind_rotamer_building_ flag is true.
template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::read_rotamer_definitions(
	utility::io::izstream & in_rotdef
)
{
	utility::vector1< Size > rotwell( 4, 0 );
	Real prob, lnprob, left, median, right;
	Size num_bins = product( parent::N_PHIPSI_BINS );

	utility::vector1< utility::vector1< Real > > lefts(   grandparent::n_possible_rots() );
	utility::vector1< utility::vector1< Real > > medians( grandparent::n_possible_rots() );
	utility::vector1< utility::vector1< Real > > rights(  grandparent::n_possible_rots() );
	utility::vector1< utility::vector1< Real > > probs(   grandparent::n_possible_rots() );

	while ( in_rotdef ) {
		char first_char = in_rotdef.peek();
		if ( first_char == '#' ) {
			std::string line; in_rotdef.getline( line );
			continue;
		}

		/// Input line is formatted:
		/// a.b.c.d. r1, r2, r3, r4
		/// e. prob
		/// f. lnprob
		/// g. non-rotameric chi pseudorotamer left range
		/// h. non-rotameric chi pseudorotamer median (mode?)
		/// i. non-rotameric chi pseudorotamer right range
		in_rotdef >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		in_rotdef >> prob >> lnprob >> left >> median >> right;

		if ( ! in_rotdef ) break; // we've read past the end of the file...

		Size const rotno( grandparent::rotwell_2_rotno( rotwell ) );
		//std::cout << "amw rotwell[ T + 1 ] = rotwell[ " << (T+1) << " ] = " << rotwell[ T+1 ] << std::endl;
		if ( rotwell[ T + 1 ] == 1 ) grandparent::mark_rotwell_exists( rotwell );
		lefts[   rotno ].push_back( left );
		medians[ rotno ].push_back( median );
		rights[  rotno ].push_back( right );
		probs[   rotno ].push_back( prob );
	}

	for ( Size ii = 1; ii <= grandparent::n_possible_rots(); ++ii ) {
		/// Verify that there are the same number of nrchi rotamers for each rotno.
		if ( lefts[ ii ].size() != 0 ) {
			/// find the first rotno with data.
			for ( Size jj = ii+1; jj <= grandparent::n_possible_rots(); ++jj ) {
				if ( lefts[ jj ].size() != 0 && lefts[ jj ].size() != lefts[ ii ].size() ) {
					std::cerr << "ERROR: Error in reading rotamer definition file" << std::endl;
					std::cerr << "ERROR: rotno's " << ii << " and " << jj;
					std::cerr << " have each defined a different number of pseudorots: ";
					std::cerr << lefts[ ii ].size() << " != " << lefts[ jj ].size();
					std::cerr << std::endl;
					utility_exit();
				}
			}
			n_nrchi_sample_bins_ = lefts[ ii ].size();
			//std::cout << "Setting n_nrchi_sample_bins_ to lefts[ " << ii << " ].size() = " << lefts[ ii ].size() << std::endl;
			break;
		}
	}

	/// After reading this file, we have seen all the rotameric chi that we are going to.
	/// Tell the base class that the rotno to packedrotno conversion is now appropriate.
	/// This triggers a call to convert rotno indices to packedrotno indices
	grandparent::declare_all_existing_rotwells_encountered();

	/// Allocate space for rotamers and rotamer-sorted-order mapping
	parent::rotamers().dimension( num_bins, grandparent::n_packed_rots(), PackedDunbrackRotamer< T, N >() );
	parent::packed_rotno_2_sorted_rotno().dimension( num_bins, grandparent::n_packed_rots() );
	parent::packed_rotno_2_sorted_rotno() = 0;

	/// save this data for rotamer creation and for binning nrchi values regardless of whether
	/// we're using bbind rotamer sampling.
	bbind_rotamers_to_sample_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );

	Size count_packed_rots( 0 );
	for ( Size ii = 1; ii <= grandparent::n_possible_rots(); ++ii ) {
		if ( lefts[ ii ].size() == 0 ) continue;
		++count_packed_rots;
		for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
			bbind_rotamers_to_sample_( jj, count_packed_rots ).left_   = lefts[ ii ][ jj ];
			bbind_rotamers_to_sample_( jj, count_packed_rots ).median_ = medians[ ii ][ jj ];
			bbind_rotamers_to_sample_( jj, count_packed_rots ).right_  = rights[ ii ][ jj ];
			bbind_rotamers_to_sample_( jj, count_packed_rots ).prob_   = probs[ ii ][ jj ];
		}
	}

	/// Allocate space for bbdep rotamer sampling for the nrchi.
	bbdep_rotamers_to_sample_.dimension( num_bins, grandparent::n_packed_rots() * n_nrchi_sample_bins_ );
	bbdep_rotsample_sorted_order_.dimension( num_bins, grandparent::n_packed_rots(), n_nrchi_sample_bins_ );
	bbdep_rotsample_sorted_order_ = 0;

	debug_assert( bbdep_nrchi_nbins_ != 0 );
	bbdep_nrc_interpdata_.resize( grandparent::n_packed_rots() );//, far2d );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		bbdep_nrc_interpdata_[ ii ].dimension( num_bins, bbdep_nrchi_nbins_ );
	}
}

/// @details Clips to range [ nrchi_lower_angle_, nrchi_lower_angle_ + nrchi_periodicity )
/// Note: this isn't really "clipping" (i.e. values outside of this range are not set to
/// the minimum or maximum of the range).  Instead, it is properly returning the value
/// within the range corresponding to a value outside of the range.
template < Size T, Size N >
Real
SemiRotamericSingleResidueDunbrackLibrary< T, N >::clip_to_nrchi_range( Real chi ) const
{
	Real chip = basic::periodic_range( chi, nrchi_periodicity_ );

	if ( chip >= nrchi_lower_angle_ + nrchi_periodicity_ ) {
		while ( chip >= nrchi_lower_angle_ + nrchi_periodicity_ ) chip -= nrchi_periodicity_;
	} else if ( chip < nrchi_lower_angle_ ) {
		while ( chip < nrchi_lower_angle_ ) chip += nrchi_periodicity_;
	}
	return chip;
}

template < Size T, Size N >
void
SemiRotamericSingleResidueDunbrackLibrary< T, N >::get_bbdep_nrchi_bin(
	Real nrchi,
	Size & bin_lower,
	Size & bin_upper,
	Real & nrchi_alpha
) const
{
	Real clipped_nrchi( clip_to_nrchi_range( nrchi ));
	grandparent::bin_angle( nrchi_lower_angle_, bbdep_nrchi_binsize_, nrchi_periodicity_,
		bbdep_nrchi_nbins_, clipped_nrchi, bin_lower, bin_upper, nrchi_alpha );
}

template < Size T, Size N >
BBDepNRChiSample< Real >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::interpolate_bbdep_nrchi_sample(
	Size const packed_rotno,
	Size const nrchi_bin,
	utility::fixedsizearray1< Size, N > const bb_bin,
	utility::fixedsizearray1< Size, N > const bb_bin_next,
	utility::fixedsizearray1< Real, N > const bb_alpha
) const
{
	utility::fixedsizearray1< Size, ( 1 << N ) > ind;
	utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > nrchi_sample( ( BBDepNRChiSample<>() ) );
	for ( Size indi = 1; indi <= ( 1 << N ); ++indi ) {
		Size index = make_conditional_index( N, parent::N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
		ind[ indi ] = bbdep_rotsample_sorted_order_( index, packed_rotno, nrchi_bin );
		nrchi_sample[ indi ] = ( bbdep_rotamers_to_sample_( index, ind[ indi ] ) );
	}

	return interpolate_bbdep_nrchi_sample( nrchi_sample, bb_alpha );
}

template < Size T, Size N >
BBDepNRChiSample< Real >
SemiRotamericSingleResidueDunbrackLibrary< T, N >::interpolate_bbdep_nrchi_sample(
	utility::fixedsizearray1< BBDepNRChiSample<>, ( 1 << N ) > const & nrchi_sample,
	utility::fixedsizearray1< Real, N > const bb_alpha
) const
{
	BBDepNRChiSample< Real > interpolated_sample;
	interpolated_sample.packed_rotno_ = nrchi_sample[ 1 ].packed_rotno_;
	interpolated_sample.nrchi_bin_ = nrchi_sample[ 1 ].nrchi_bin_;

	utility::fixedsizearray1< Real, N > dummy_dprob, dummy_dmean, dummy_dsd;
	utility::fixedsizearray1< Real, ( 1 << N ) > rot_prob, nrchi_mean, nrchi_sd;
	utility::fixedsizearray1< Real, N > binw( parent::PHIPSI_BINRANGE );

	for ( Size i = 1; i <= ( 1 << N ) /*nrchi_sample.size()*/; ++i ) {
		rot_prob[ i ] = static_cast< Real >   ( nrchi_sample[ i ].prob_ );
		nrchi_mean[ i ] = static_cast< Real > ( nrchi_sample[ i ].nrchi_mean_ );
		nrchi_sd[ i ] = static_cast< Real >   ( nrchi_sample[ i ].nrchi_sd_ );
	}

	interpolate_polylinear_by_value( rot_prob,   bb_alpha, binw, false /*treat_as_angles*/, interpolated_sample.prob_,       dummy_dprob );
	interpolate_polylinear_by_value( nrchi_mean, bb_alpha, binw, true  /*treat_as_angles*/, interpolated_sample.nrchi_mean_, dummy_dmean );
	interpolate_polylinear_by_value( nrchi_sd,   bb_alpha, binw, false /*treat_as_angles*/, interpolated_sample.nrchi_sd_,   dummy_dsd   );


	return interpolated_sample;
}

template < Size T, Size N >
void
initialize_and_read_srsrdl(
	SemiRotamericSingleResidueDunbrackLibrary< T, N > & srsrdl,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle,
	utility::io::izstream & rotamer_definitions,
	utility::io::izstream & regular_library,
	utility::io::izstream & continuous_minimization_bbdep
)
{
	initialize_srsrdl( srsrdl, nrchi_is_symmetric, nrchi_start_angle );
	srsrdl.read_from_files( rotamer_definitions, regular_library, continuous_minimization_bbdep );
}

template < Size T, Size N >
void
initialize_srsrdl(
	SemiRotamericSingleResidueDunbrackLibrary< T, N > & srsrdl,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle
)
{
	Real const nrchi_periodicity  = nrchi_is_symmetric ? 180.0 : 360.0;
	Real const nrchi_bbdep_step_size = nrchi_is_symmetric ?   5.0 :  10.0;
	Real const nrchi_bbind_step_size = nrchi_is_symmetric ?   0.5 :   1.0;

	srsrdl.set_nrchi_periodicity( nrchi_periodicity );
	srsrdl.set_nonrotameric_chi_start_angle( nrchi_start_angle );
	srsrdl.set_nonrotameric_chi_bbdep_scoring_step_size( nrchi_bbdep_step_size );
	srsrdl.set_nonrotameric_chi_bbind_scoring_step_size( nrchi_bbind_step_size );
}

} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_TMPL_HH
