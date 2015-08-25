// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/dunbrack/RotamericSingleResidueDunbrackLibrary.hh
/// @brief   Declaration of rotameric rotamer libraries from Jun08 (as opposed to semi-rotameric)
/// @author  Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_tmpl_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_tmpl_hh

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

// Unit Headers
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/ChiSet.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <core/pose/Pose.hh>
#include <basic/basic.hh>
#include <basic/interpolate.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/backtrace.hh>
#include <utility/io/ozstream.hh>
#include <utility/LexicographicalIterator.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/MathNTensor.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/interpolation/spline/PolycubicSpline.hh>
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
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
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
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
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
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
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
#include <utility/io/mpistream.hh>
#include <utility/io/mpistream.ipp>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
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
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
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
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <basic/Tracer.hh>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

namespace core {
namespace pack {
namespace dunbrack {

template < Size T, Size N >
RotamericSingleResidueDunbrackLibrary< T, N >::RotamericSingleResidueDunbrackLibrary(
	AA const aa_in,
	bool dun02
) :
parent( aa_in, T, dun02 )
{
	// This is horrible but temporarily necessarily until we figure out how to distribute full beta rotlibs
	// Betas automatically (if 3 BB) have phipsi binrange of 30
	if ( N > 2 ) {
		for ( Size i = 1; i <= N; ++i ) {
			N_PHIPSI_BINS[ i ] = 12;
			PHIPSI_BINRANGE[ i ] = 360 / N_PHIPSI_BINS[ i ];
		}
	} else {
		for ( Size i = 1; i <= N; ++i ) {
			N_PHIPSI_BINS[ i ] = 36;
			PHIPSI_BINRANGE[ i ] = 360 / N_PHIPSI_BINS[ i ];
		}
	}
}

template < Size T, Size N >
RotamericSingleResidueDunbrackLibrary< T, N >::~RotamericSingleResidueDunbrackLibrary() throw()
{}


template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	return eval_rotameric_energy_deriv( rsd, scratch, false );
}


template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	/// most of the work is done in this call
	Real score = eval_rotameric_energy_deriv( rsd, scratch, true );

	//Multiplier for D-amino acids:
	const core::Real d_multiplier = rsd.has_property( "D_AA" ) ? -1.0 : 1.0;

	if ( score != score ) { // NaN check
		score = 0;
		std::cerr << "NaN at residue rsd: " << rsd.seqpos() << " " << rsd.name() << std::endl;
	}
	if ( score > 1e16 ) { // inf check
		std::cerr << "inf at residue rsd: " << rsd.seqpos() << " " << rsd.name() <<  " " << score << std::endl;
		score = 0;
	}

	/// sum derivatives.
	Real5 & dE_dbb(  scratch.dE_dbb() );
	Real5 & dE_dbb_dev(  scratch.dE_dbb_dev() );
	Real5 & dE_dbb_rot(  scratch.dE_dbb_rot() );
	Real4 & dE_dchi( scratch.dE_dchi() );
	Real4 & dE_dchi_dev( scratch.dE_dchi_dev() );

	// p0 - the base probability -- not modified by the chi-dev penalty
	Real const rotprob( scratch.rotprob() );
	Real const invp( ( rotprob == Real( 0.0 ) ) ? 0.0 : -1.0 / rotprob );

	for ( Size bbi = 1; bbi <= N; ++bbi ) {
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {
			dE_dbb[ bbi ]   = d_multiplier * ( scratch.dneglnrotprob_dbb()[ bbi ] + scratch.dchidevpen_dbb()[ bbi ] );
			dE_dbb_dev[ bbi ] = d_multiplier * scratch.dchidevpen_dbb()[ bbi ];
			dE_dbb_rot[ bbi ] = d_multiplier * scratch.dneglnrotprob_dbb()[ bbi ];
		} else {
			dE_dbb[ bbi ]   = d_multiplier * ( invp * scratch.drotprob_dbb()[ bbi ] + scratch.dchidevpen_dbb()[ bbi ] );
			dE_dbb_dev[ bbi ] = d_multiplier * scratch.dchidevpen_dbb()[ bbi ];
			dE_dbb_rot[ bbi ] = d_multiplier * invp * scratch.drotprob_dbb()[ bbi ];
		}

		// Correction for entropy
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
			dE_dbb[ bbi ]   += d_multiplier * scratch.dentropy_dbb()[ bbi ];
			dE_dbb_rot[ bbi ] += d_multiplier * scratch.dentropy_dbb()[ bbi ];
		}
	}

	for ( Size i = 1; i <= T; ++i ) {
		dE_dchi[ i ]     = d_multiplier * scratch.dchidevpen_dchi()[ i ];
		dE_dchi_dev[ i ] = d_multiplier * scratch.dchidevpen_dchi()[ i ];
	}

	correct_termini_derivatives( rsd, scratch );

	return score;
}

template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::eval_rotameric_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const
{
	//There's probably a better way to check this.
	debug_assert( rsd.type().backbone_aa() == aa() ||
		( core::chemical::is_canonical_D_aa( rsd.aa() )
		&& core::chemical::get_L_equivalent( rsd.aa() ) == aa() ) );

	ChiVector chi ( rsd.chi() );
	//Invert if we're dealing with a D-amino acid.
	if ( core::chemical::is_canonical_D_aa( rsd.aa() ) ) for ( core::Size i = 1; i <= chi.size(); i++ ) chi[ i ] *= -1.0;

	Real4 & chimean    ( scratch.chimean() );
	Real4 & chisd    ( scratch.chisd()   );
	Real4 & chidev    ( scratch.chidev()  );
	Real4 & chidevpen   ( scratch.chidevpen() );
	Real4 & dchidevpen_dchi  ( scratch.dchidevpen_dchi() );
	Real5 & drotprob_dbb  ( scratch.drotprob_dbb() );
	Real5 & dneglnrotprob_dbb ( scratch.dneglnrotprob_dbb() ); // for bicubic interpolation
	Real5 & dchidevpen_dbb  ( scratch.dchidevpen_dbb()  );
	FiveReal4 & dchimean_dbb ( scratch.dchimean_dbb()  );
	FiveReal4 & dchisd_dbb  ( scratch.dchisd_dbb()  );

	std::fill( chimean.begin(),        chimean.end(),           0.0 );
	std::fill( chisd.begin(),        chisd.end(),             0.0 );
	std::fill( chidev.begin(),        chidev.end(),            0.0 );
	std::fill( chidevpen.begin(),       chidevpen.end(),         0.0 );
	std::fill( dchidevpen_dchi.begin(),   dchidevpen_dchi.end(),   0.0 );
	std::fill( drotprob_dbb.begin(),      drotprob_dbb.end(),      0.0 );
	std::fill( dneglnrotprob_dbb.begin(), dneglnrotprob_dbb.end(), 0.0 );
	std::fill( dchidevpen_dbb.begin(),    dchidevpen_dbb.end(),    0.0 );
	for ( Size bbi = 1; bbi <= N; ++bbi ) {
		std::fill( dchimean_dbb[ bbi ].begin(), dchimean_dbb[ bbi ].end(), 0.0 );
		std::fill( dchisd_dbb[ bbi ].begin(),  dchisd_dbb[ bbi ].end(), 0.0 );
	}

	scratch.fa_dun_tot() = 0;
	scratch.fa_dun_rot() = 0;
	scratch.fa_dun_semi() = 0;
	scratch.fa_dun_dev() = 0;

	// compute rotamer number from chi
	Size4 & rotwell( scratch.rotwell() );

	/// Don't use derived class's version of this function.
	RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi_static( chi, rotwell );

	Size packed_rotno( rotwell_2_packed_rotno( rotwell ) );

	// panic!  Extremely unlikely rotamer found.  Find another rotamer that has at least some probability,
	// and move this rotamer toward it -- do so in a predictable manner so that the score function is continuous
	// as it tries to move from this rotamer to another.
	if ( packed_rotno == 0 ) packed_rotno = find_another_representative_for_unlikely_rotamer( rsd, rotwell );

	PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
	interpolate_rotamers( rsd, scratch, packed_rotno, interpolated_rotamer );

	if ( dun02() ) {
		for ( Size ii = 1; ii <= T; ++ii ) chidev[ ii ] = subtract_chi_angles( chi[ ii ], chimean[ ii ], aa(), ii );
	} else {
		for ( Size ii = 1; ii <= T; ++ii ) chidev[ ii ] = basic::periodic_range( chi[ ii ] - chimean[ ii ], 360 );
	}

	for ( Size ii = 1; ii <= T; ++ii ) {
		/// chidev penalty: define a gaussian with a height of 1 and a standard deviation of chisd[ ii ];
		/// exp( -1 * (chi_i - chi_mean_i )**2 / ( 2 * chisd_ii**2 ) )
		/// Pretend this gaussian defines a probability for having an angle at a certain deviation from the mean.
		/// Convert that probability into an energy:
		/// -ln p = -ln exp( -1* (chi_i - chi_mean_i)**2/(2 chisd_i**2) ) = (chi_i - chi_mean_i)**2/(2 chisd_i**2)
		chidevpen[ ii ] = chidev[ ii ] * chidev[ ii ] / ( 2 * chisd[ ii ] * chisd[ ii ] );
		/// Add in gaussian height normalization
		/// p = prod( i, 1/(stdv_i*sqrt(2*pi)) * exp( -1* (chi_i - chi_mean_i)**2/(2 chisd_i**2) ) )
		/// -ln(p) == sum( i, (chi_i - chi_mean_i)**2/(2 chisd_i**2) + log(stdv_i*sqrt(2*pi)) )
		///  == sum( i, (chi_i - chi_mean_i)**2/(2 chisd_i**2) + log(stdv_i) + log(sqrt(2*pi)))  <-- where sum(i,sqrt2pi) is just a constant per amino acid
		///     and can be treated as part of the reference energies.
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_normsd ] ) {
			chidevpen[ ii ] += std::log( chisd[ ii ] );
		}
	}

	Real chidevpensum( 0.0 );
	for ( Size ii = 1; ii <= T; ++ii ) chidevpensum += chidevpen[ ii ];


	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {
		scratch.fa_dun_tot() = scratch.negln_rotprob() + chidevpensum;
		scratch.fa_dun_rot() = scratch.negln_rotprob();
		scratch.fa_dun_dev() = chidevpensum;
	} else {
		scratch.fa_dun_rot() = -std::log(scratch.rotprob());
		scratch.fa_dun_tot() = scratch.fa_dun_rot() + chidevpensum;
		scratch.fa_dun_dev() = chidevpensum;
	}

	// Corrections for Shanon Entropy
	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
		scratch.fa_dun_rot() += scratch.entropy();
		scratch.fa_dun_tot() += scratch.entropy();
	}

	Real const score( scratch.fa_dun_tot() );

	if ( ! eval_deriv ) return score;

	for ( Size bbi = 1; bbi <= N; ++bbi ) dchidevpen_dbb[ bbi ] = 0.0;

	for ( Size ii = 1; ii <= T; ++ii ) {

		/// Backbone derivatives for chi-dev penalty.
		/// Xmean_i and sd_i both depend on phi and psi.
		/// Let: f = (X_i-Xmean_i)**2
		/// Let: g = 2 sd_i**2
		/// Then, chidevpen = f/g
		/// and, dchidevpen = (f'g - fg')/(gg)
		Real const f   = chidev[ ii ] * chidev[ ii ];
		Real const fprime = -2 * chidev[ ii ];
		Real const g   = 2 * chisd[ ii ] * chisd[ ii ];
		Real const gprime = 4 * chisd[ ii ];
		Real const invgg  = 1 / ( g * g );

		/// also need to add in derivatives for log(stdv) height normalization
		/// dE/dphi = (above) + 1/stdv * dstdv/dphi
		for ( Size bbi = 1; bbi <= N; ++bbi ) {
			scratch.dE_dbb_dev_perchi()[ bbi ][ ii ] = ( g * fprime * dchimean_dbb[ bbi ][ ii ]
				- f * gprime *   dchisd_dbb[ bbi ][ ii ] ) * invgg;

			dchidevpen_dbb[ bbi ] += scratch.dE_dbb_dev_perchi()[ bbi ][ ii ];

			// Derivatives for the change in the Gaussian height normalization due to sd changing as a function of phi
			if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_normsd ] ) {
				dchidevpen_dbb[ bbi ]      += 1.0 / chisd[ ii ] * dchisd_dbb[ bbi ][ ii ];
				scratch.dE_dbb_dev_perchi()[ bbi ][ ii ] += 1.0 / chisd[ ii ] * dchisd_dbb[ bbi ][ ii ];
			}
		}
		dchidevpen_dchi[ ii ] = chidev[ ii ] / ( chisd[ ii ] * chisd[ ii ] );
	}
	return score;
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::assign_random_rotamer_with_bias(
	conformation::Residue const & rsd,
	pose::Pose const & /*pose*/,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	Size packed_rotno( 0 );
	assign_random_rotamer( rsd, scratch, RG, new_chi_angles, perturb_from_rotamer_center, packed_rotno );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::assign_random_rotamer(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center,
	Size & packed_rotno
) const
{
	Real random_prob = RG.uniform();

	utility::fixedsizearray1< Real, N > const bbs( get_bbs_from_rsd( rsd ) );
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;

	/// Go through rotamers in decreasing order by probability and stop when the
	Size count = 0;
	while ( random_prob > 0 ) {
		Size index = make_index( N, N_PHIPSI_BINS, bb_bin );
		packed_rotno = rotamers_( index, ++count ).packed_rotno();
		interpolate_rotamers( scratch, packed_rotno, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );
		random_prob -= interpolated_rotamer.rotamer_probability();
		//loop condition might end up satisfied even if we've walked through all possible rotamers
		// if the chosen random number was nearly 1
		// (and interpolation introduced a tiny bit of numerical noise).
		if ( count == rotamers_.size2() ) break;
	}
	assign_chi_for_interpolated_rotamer( interpolated_rotamer, rsd, RG, new_chi_angles, perturb_from_rotamer_center );

}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::assign_chi_for_interpolated_rotamer(
	PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer,
	conformation::Residue const & rsd,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	new_chi_angles.resize( rsd.nchi() );
	if ( ! perturb_from_rotamer_center ) {
		for ( Size ii = 1; ii <= T; ++ii ) new_chi_angles[ ii ] = interpolated_rotamer.chi_mean( ii );
	} else {
		for ( Size ii = 1; ii <= T; ++ii ) new_chi_angles[ ii ] = interpolated_rotamer.chi_mean( ii ) + RG.gaussian() * interpolated_rotamer.chi_sd(ii);
	}

	/// Set any remaining chi uniformly? ( proton chi)
	for ( Size ii = T + 1; ii <= rsd.nchi(); ++ii ) new_chi_angles[ ii ] = RG.uniform() * 360.0 - 180.0;
}


/// @details The new rotamer library represents only 75 of 81 possible arginine rotamers,
/// and 75 of 81 lysine rotamers. In the unlikely event that a rotamer is encountered
/// that's not represented in the library, find another rotamer to represent it.
/// Ideally, this stand-in rotamer would be closest to the input rotamer in physical geometry.
/// (sidechain atom rms, e.g.)
/// The following code instead first looks through all possible rotamers with a Hamming distance
/// of one from input rotamer trying to find one that works, looking first for rotamers from
/// the furthest chi toward the closest chi.  If no such rotamer may be found, it gives up and
/// returns the most-probable rotamer for a phi-psi bin.
///
/// This function modifies the "rotwell" assigned to this rotamer so that later code that relies
/// on the consistency of the rotwell and packed_rotno information will behave correctly.
template < Size T, Size N >
Size
RotamericSingleResidueDunbrackLibrary< T, N >::find_another_representative_for_unlikely_rotamer(
	conformation::Residue const & rsd,
	Size4 & rotwell
) const
{
	/// Start from the furthest chi and work inwards trying to find another rotamer that could work
	for ( Size ii = T; ii >= 1; --ii ) {
		Size ii_orig_value = rotwell[ ii ];
		for ( Size jj = 1; jj <= 3; ++jj ) { /// 3 == NUMBER OF ROTAMERIC CHI BINS
			if ( jj == ii_orig_value ) continue;
			rotwell[ ii ] = jj;
			Size new_packed_rotno = rotwell_2_packed_rotno( rotwell );
			if ( new_packed_rotno != 0 ) return new_packed_rotno;
		}
		/// didn't find one for chi_ii, reset rotwell to original value
		rotwell[ ii ] = ii_orig_value;
	}

	// At this point, we've examined all the rotamer wells with a Hamming distance of 1 from
	// the input rotamer well and have come up empty.
	// Just go for the most likely rotamer in this well -- no guarantee that as the rotamer swings from the
	// current rotamer well to the target rotamer well that it doesn't first swing through some other

	utility::fixedsizearray1< Real, N > bbs( get_bbs_from_rsd( rsd ) );
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Size index = make_index( N, N_PHIPSI_BINS, bb_bin );
	Size packed_rotno = rotamers_( index, 1 ).packed_rotno();
	packed_rotno_2_rotwell( packed_rotno, rotwell );
	return packed_rotno;
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::correct_termini_derivatives(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	// mt: for the termini, these derivatives are 0, because these "pseudo"
	// mt: torsions are kept fixed.
	// amw: assume that we are more interested in 1 and N than PHI and PSI per se; this reduces in the 1,2 case
	if ( rsd.is_lower_terminus() ) {
		// amw: now 1 is actually better here than RotamerLibraryScratchSpace::AA_PHI_INDEX because the right element to alter is always element 1 but not necessarily has the meaning of "PHI"
		scratch.dE_dbb()[ 1 ] = 0;
		scratch.dE_dbb_dev()[ 1 ] = 0;
		scratch.dE_dbb_rot()[ 1 ] = 0;
		scratch.dE_dbb_semi()[ 1 ] = 0;
	}
	if ( rsd.is_upper_terminus() ) {
		scratch.dE_dbb()[ N ] = 0;
		scratch.dE_dbb_dev()[ N ] = 0;
		scratch.dE_dbb_rot()[ N ] = 0;
		scratch.dE_dbb_semi()[ N ] = 0;
	}


}

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::best_rotamer_energy(
	conformation::Residue const & rsd,
	bool curr_rotamer_only,
	RotamerLibraryScratchSpace & scratch
) const
{
	Size const num_packed_rots = ( 1 << N );
	Real maxprob( 0 );
	if ( curr_rotamer_only ) {
		Size4 & rotwell( scratch.rotwell() );
		RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi_static( rsd.chi(), rotwell );
		Size const packed_rotno = rotwell_2_packed_rotno( rotwell );
		PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
		interpolate_rotamers( rsd, scratch, packed_rotno, interpolated_rotamer );
		maxprob = interpolated_rotamer.rotamer_probability();

	} else {
		utility::fixedsizearray1< Real, N > const bbs = get_bbs_from_rsd( rsd );
		utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
		utility::fixedsizearray1< Real, N > bb_alpha;
		get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

		utility::fixedsizearray1< Size, (1 << N ) > packed_rotnos;
		for ( Size indi = 1; indi <= num_packed_rots; ++indi ) {
			Size index = make_conditional_index( N, N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
			packed_rotnos[ indi ] = rotamers_( index, 1 ).packed_rotno();
		}

		PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
		for ( Size ii = 1; ii <= num_packed_rots; ++ii ) {
			interpolate_rotamers( rsd, scratch, packed_rotnos[ ii ], interpolated_rotamer );
			maxprob = ( maxprob < interpolated_rotamer.rotamer_probability() ?
				interpolated_rotamer.rotamer_probability() : maxprob );
		}
	}
	return -1 * std::log( maxprob );
}

/// @details The Dunbrack library's phi/psi data for a grid point
/// (x,y) collects data in the neighborhood of (x,y).  As an interpolation
/// point p moves toward a grid point (x,y), the (x,y) share in the
/// interpolated value goes to 1.  This is distinct from having interpolation wells
/// where an interpolation in the center of a well produces the maximum contribution
/// from that well.  Most of the basic::interpolation code is designed for
/// the second interpretation of interpolation, and so CTSA's Dunbrack library code did
/// funky thinks like shift by 5 degrees so that the basic::interpolation code could
/// shift it back by 5 degrees again.  The code below makes no such shift.
///
/// The alpha fraction is the distance along each axis that the interpolation point
/// has progressed from the lower grid point toward the upper grid point; it ranges
/// from 0 to 1.
template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::get_bb_bins(
	utility::fixedsizearray1< Real, N > bbs,
	utility::fixedsizearray1< Size, N > & bb_bin,
	utility::fixedsizearray1< Size, N > & bb_bin_next,
	utility::fixedsizearray1< Real, N > & bb_alpha
) const
{
	for ( Size ii = 1; ii <= N; ++ii ) {
		parent::bin_angle( -180.0, PHIPSI_BINRANGE[ ii ], 360.0, N_PHIPSI_BINS[ ii ], basic::periodic_range( bbs[ ii ], 360 ),
			bb_bin[ ii ], bb_bin_next[ ii ], bb_alpha[ ii ] );
	}

	verify_bb_bins( bbs, bb_bin, bb_bin_next );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::get_bb_bins(
	utility::fixedsizearray1< Real, N > bbs,
	utility::fixedsizearray1< Size, N > & bb_bin
) const
{
	utility::fixedsizearray1< Size, N > bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::interpolate_rotamers(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	Size const packed_rotno,
	PackedDunbrackRotamer< T, N, Real > & interpolated_rotamer
) const
{
	utility::fixedsizearray1< Real, N > bbs = get_bbs_from_rsd( rsd );
	//For D-amino acids, invert phi and psi:
	if ( core::chemical::is_canonical_D_aa( rsd.aa() ) ) for ( Size bbi = 1; bbi <= N; ++bbi ) bbs[ bbi ] *= -1.0;

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	interpolate_rotamers( scratch, packed_rotno, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::interpolate_rotamers(
	RotamerLibraryScratchSpace & scratch,
	Size const packed_rotno,
	utility::fixedsizearray1< Size, N > const bb_bin,
	utility::fixedsizearray1< Size, N > const bb_bin_next,
	utility::fixedsizearray1< Real, N > const bb_alpha,
	PackedDunbrackRotamer< T, N, Real > & interpolated_rotamer
) const
{
	using namespace basic;

	interpolated_rotamer.packed_rotno() = packed_rotno;
	utility::fixedsizearray1< Size, ( 1 << N ) > sorted_rotno;
	utility::vector1< PackedDunbrackRotamer< T, N > > rot( 1 << N );
	utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << N ) >, ( 1 << N ) > derivs_vec;
	utility::fixedsizearray1< Real, ( 1 << N ) > rotprob;

	for ( Size sri = 1; sri <= ( 1 << N ); ++sri ) {
		Size index = make_conditional_index( N, N_PHIPSI_BINS, sri, bb_bin_next, bb_bin );
		sorted_rotno[ sri ] = packed_rotno_2_sorted_rotno_( index, packed_rotno );
		rot.push_back( rotamers_( index, sorted_rotno[ sri ] ) );
		for ( Size i = 1; i <= T; ++i ) rot[ sri ].chi_mean( i ) = rotamers_( index, sorted_rotno[ sri ] ).chi_mean( i );
		for ( Size i = 1; i <= T; ++i ) rot[ sri ].chi_sd(   i ) = rotamers_( index, sorted_rotno[ sri ] ).chi_sd(   i );
		rot[ sri ].rotamer_probability() = rotamers_( index, sorted_rotno[ sri ] ).rotamer_probability();
		for ( Size di = 1; di <= ( 1 << N ); ++di ) {
			derivs_vec[ sri ][ di ] = static_cast< Real >( rotamers_( index, sorted_rotno[ sri ] ).n_derivs()[ di ] );
		}
		rotprob[ sri ] = static_cast< Real >( rot[ sri ].rotamer_probability() );
		if ( rotprob[ sri ] <= 1e-6 ) rotprob[ sri ] = 1e-6;
	}

	utility::fixedsizearray1< Real, N > binw( PHIPSI_BINRANGE );
	utility::fixedsizearray1< Real, N > scratch_drotprob_dbb;

	/// Note: with bicubic splines to interpolate the energy for a rotamer, bilinear
	/// interpolation of the rotamer probabilities no longer serves a useful purpose
	/// for the 2002 library.  However, the 2008 and 2010 libraries are not yet up to
	/// speed with cubic interpolation (requiring tricubic splines), so keep this value
	/// around for the moment.  Also being preserved for the sake of backwards compatibility.
	interpolate_polylinear_by_value( rotprob, bb_alpha, binw, false, // treat as angles
		scratch.rotprob(), scratch_drotprob_dbb );

	for ( Size i = 1; i <= N; ++i ) scratch.drotprob_dbb()[ i ] = scratch_drotprob_dbb[ i ];

	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {

		// populate derivatives vector
		utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << N ) >, ( 1 << N ) > n_derivs;
		for ( Size i = 1; i <= ( 1 << N ); ++i ) {
			for ( Size j = 1; j <= ( 1 << N ); ++j ) n_derivs[i][j] = derivs_vec[ j ][ i ];
		}
		utility::fixedsizearray1< Real, N > binw( PHIPSI_BINRANGE );
		utility::fixedsizearray1< Real, N > scratch_dneglnrotprob_dbb;
		alternate_tricubic_interpolation( n_derivs, bb_alpha, binw, scratch.negln_rotprob(), scratch_dneglnrotprob_dbb );
		for ( Size i = 1; i <= N; ++i ) scratch.dneglnrotprob_dbb()[ i ] = scratch_dneglnrotprob_dbb[ i ];
		interpolated_rotamer.rotamer_probability() = std::exp( -scratch.negln_rotprob() );

	} else {
		interpolated_rotamer.rotamer_probability() = scratch.rotprob();
	}

	// Entropy correction
	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
		// populate S_n_derivs
		utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << N ) >, ( 1 << N ) > S_n_derivs;
		for ( Size i = 1; i <= ( 1 << N ); ++i ) {
			for ( Size j = 1; j <= ( 1 << N ); ++j ) S_n_derivs[i][j] = ShanonEntropy_n_derivs_[ i ][ j ];
		}
		utility::fixedsizearray1< Real, N > binw( PHIPSI_BINRANGE );

		utility::fixedsizearray1< Real, N > scratch_dentropy_dbb;
		alternate_tricubic_interpolation( S_n_derivs, bb_alpha, binw, scratch.entropy(), scratch_dentropy_dbb );
		for ( Size i = 1; i <= N; ++i ) scratch.dentropy_dbb()[ i ] = scratch_dentropy_dbb[ i ];
	}

	for ( Size ii = 1; ii <= T; ++ii ) {

		utility::fixedsizearray1< Real, ( 1 << N ) > chi_mean;
		utility::fixedsizearray1< Real, ( 1 << N ) > chi_sd;
		for ( Size roti = 1; roti <= ( 1 << N ); ++roti ) {
			chi_mean[ roti ] = static_cast< Real >( rot[ roti ].chi_mean( ii ) );
			chi_sd[ roti ]   = static_cast< Real >( rot[ roti ].chi_sd( ii ) );
		}

		utility::fixedsizearray1< Real, N > scratch_dchi_mean;
		utility::fixedsizearray1< Real, N > scratch_dchi_sd;

		for ( Size bbi = 1; bbi <= N; ++bbi ) {
			scratch_dchi_mean[ bbi ] = scratch.dchimean_dbb()[ bbi ][ ii ];
			scratch_dchi_sd[ bbi ]   = scratch.dchisd_dbb()[   bbi ][ ii ];
		}
		utility::fixedsizearray1< Real, N > binw( PHIPSI_BINRANGE );

		interpolate_polylinear_by_value( chi_mean, bb_alpha, binw, true, scratch.chimean()[ ii ], scratch_dchi_mean );
		interpolated_rotamer.chi_mean( ii ) = scratch.chimean()[ ii ];
		interpolate_polylinear_by_value( chi_sd, bb_alpha, binw, false, scratch.chisd()[ ii ], scratch_dchi_sd );
		interpolated_rotamer.chi_sd( ii ) = scratch.chisd()[ ii ];

		for ( Size bbi = 1; bbi <= N; ++bbi ) {
			scratch.dchimean_dbb()[ bbi ][ ii ] = scratch_dchi_mean[ bbi ];
			scratch.dchisd_dbb()[ bbi ][ ii ]   = scratch_dchi_sd[ bbi ];
		}
	}
}

/// @details Handle lower-term residues by returning a "neutral" phi value
template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::get_phi_from_rsd(
	conformation::Residue const & rsd
) const
{
	debug_assert( rsd.is_protein() || rsd.is_peptoid() );
	if ( rsd.is_lower_terminus() ) return parent::NEUTRAL_PHI;
	else return rsd.mainchain_torsion( 1 );
}

/// @details Handle upper-term residues by returning a "neutral" psi value
template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::get_psi_from_rsd(
	conformation::Residue const & rsd
) const
{
	debug_assert( rsd.is_protein() || rsd.is_peptoid() );
	if ( rsd.is_upper_terminus() ) return parent::NEUTRAL_PSI;
	else return rsd.mainchain_torsion( 2 );
}

/// @details Handle lower-term residues by returning a "neutral" phi value
template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::get_bb_from_rsd(
	Size bbn,
	conformation::Residue const & rsd
) const
{
	debug_assert( rsd.is_protein() || rsd.is_peptoid() );
	if      ( rsd.is_lower_terminus() && bbn == 1 ) return parent::NEUTRAL_PHI;
	else if ( rsd.is_upper_terminus() && bbn == N ) return parent::NEUTRAL_PSI;
	else return rsd.mainchain_torsion( bbn );
}

/// @details Handle upper-term residues by returning a "neutral" psi value
template < Size T, Size N >
utility::fixedsizearray1< Real, N >
RotamericSingleResidueDunbrackLibrary< T, N >::get_bbs_from_rsd(
	conformation::Residue const & rsd
) const
{
	debug_assert( rsd.is_protein() || rsd.is_peptoid() );

	utility::fixedsizearray1< Real, N > tors;
	for ( Size ii = 1; ii <= N; ++ii ) {
		if      ( rsd.is_lower_terminus() && ii == 1 ) tors[ ii ] = parent::NEUTRAL_PHI;
		else if ( rsd.is_upper_terminus() && ii == N ) tors[ ii ] = parent::NEUTRAL_PSI;
		else tors[ ii ] = rsd.mainchain_torsion( ii );
	}
	return tors;
}


template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::fill_rotamer_vector(
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

	///Determine whether this is a D-amino acid:
	const core::Real d_multiplier = existing_residue.has_property( "D_AA" ) ? -1.0 : 1.0;

	/// Save backbone interpolation data for reuse
	utility::fixedsizearray1< Real, N > bbs = get_bbs_from_rsd( existing_residue );
	for ( Size bbi = 1; bbi <= N; ++bbi ) bbs[ bbi ] *= d_multiplier;
	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Real const requisit_probability = probability_to_accumulate_while_building_rotamers( buried ); // ( buried  ? 0.98 : 0.95 )
	Real accumulated_probability( 0.0 );

	Size const max_rots_that_can_be_built = n_packed_rots();
	Size count_rotamers_built = 0;
	while ( accumulated_probability < requisit_probability ) {
		// Iterate through rotamaers in decreasing order of probabilities, stopping once the requisit probility is hit.
		++count_rotamers_built;

		Size index = make_index( N, N_PHIPSI_BINS, bb_bin );
		Size const packed_rotno00 = rotamers_( index, count_rotamers_built ).packed_rotno();
		PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
		interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );

		build_rotamers( pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried, rotamers,
			interpolated_rotamer );

		accumulated_probability += interpolated_rotamer.rotamer_probability();
		if ( count_rotamers_built == max_rots_that_can_be_built ) break; // this shouldn't happen...
	}
}

template < Size T, Size N >
utility::vector1< DunbrackRotamerSampleData >
RotamericSingleResidueDunbrackLibrary< T, N >::get_all_rotamer_samples(
	utility::fixedsizearray1< Real, 5 > bbs2
) const
{
	RotamerLibraryScratchSpace scratch;

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;

	utility::fixedsizearray1< Real, N > bbs;
	for ( Size i = 1 ; i <= N; ++i ) bbs[ i ] = bbs2[ i ];
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Size const n_rots = n_packed_rots();
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	for ( Size ii = 1; ii <= n_rots; ++ii ) {
		// Iterate through rotamaers in decreasing order of probabilities
		Size index = make_index( N, N_PHIPSI_BINS, bb_bin );
		Size const packed_rotno00 = rotamers_( index, ii ).packed_rotno();
		PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
		interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );

		DunbrackRotamerSampleData sample( false );
		sample.set_nchi( T );
		sample.set_rotwell( parent::packed_rotno_2_rotwell( interpolated_rotamer.packed_rotno() ));
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, interpolated_rotamer.chi_mean( jj ) );
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, interpolated_rotamer.chi_sd( jj ) );
		sample.set_prob( interpolated_rotamer.rotamer_probability() );

		all_rots.push_back( sample );
	}
	return all_rots;
}

template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::get_probability_for_rotamer(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	Size const n_rot = 1 << N;

	utility::fixedsizearray1< Size, N >  bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N >  bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	utility::vector1< PackedDunbrackRotamer< T, N > > rot;
	utility::fixedsizearray1< Size, ( 1 << N ) > sorted_rotno;

	Size ind00 = make_index( N, N_PHIPSI_BINS, bb_bin );
	utility::fixedsizearray1< Real, ( 1 << N ) > interp_probs;
	rot.push_back( rotamers_( ind00, rot_ind ) );
	rot[ 1 ].rotamer_probability() = rotamers_( ind00, rot_ind ).rotamer_probability();
	interp_probs[ 1 ] = static_cast< Real >( rot[ 1 ].rotamer_probability() );
	Size const packed_rotno00 = rot[ 1 ].packed_rotno();

	for ( Size indi = 2; indi <= n_rot; ++indi ) {
		Size index = make_conditional_index( N, N_PHIPSI_BINS, indi, bb_bin_next, bb_bin );
		sorted_rotno[ indi ] = packed_rotno_2_sorted_rotno_( index, packed_rotno00 );
		rot.push_back( rotamers_( index, sorted_rotno[ indi ] ) );
		rot[ indi ].rotamer_probability() = rotamers_( index, sorted_rotno[ indi ] ).rotamer_probability();
		interp_probs[ indi ] = static_cast< Real >( rot[ indi ].rotamer_probability() );
	}

	Real rot_prob;
	utility::fixedsizearray1< Real, N > dummy;
	utility::fixedsizearray1< Real, N > binw( PHIPSI_BINRANGE );

	interpolate_polylinear_by_value( interp_probs, bb_alpha, binw, false, rot_prob, dummy );
	return rot_prob;
}

template < Size T, Size N >
Real
RotamericSingleResidueDunbrackLibrary< T, N >::get_probability_for_rotamer(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	utility::fixedsizearray1< Real, N > bbs;
	bbs[1] = phi;
	bbs[2] = psi;
	return get_probability_for_rotamer( bbs, rot_ind );
}

template < Size T, Size N >
DunbrackRotamerSampleData
RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer(
	utility::fixedsizearray1< Real, N > bbs,
	Size rot_ind
) const
{
	RotamerLibraryScratchSpace scratch;

	utility::fixedsizearray1< Size, N > bb_bin, bb_bin_next;
	utility::fixedsizearray1< Real, N > bb_alpha;
	get_bb_bins( bbs, bb_bin, bb_bin_next, bb_alpha );

	Size index = make_index( N, N_PHIPSI_BINS, bb_bin );
	PackedDunbrackRotamer< T, N > const & rot00( rotamers_( index, rot_ind ) );
	Size const packed_rotno00 = rot00.packed_rotno();
	PackedDunbrackRotamer< T, N, Real > interpolated_rotamer;
	interpolate_rotamers( scratch, packed_rotno00, bb_bin, bb_bin_next, bb_alpha, interpolated_rotamer );

	DunbrackRotamerSampleData sample( false );
	sample.set_nchi( T );
	sample.set_rotwell( parent::packed_rotno_2_rotwell( interpolated_rotamer.packed_rotno() ));
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, interpolated_rotamer.chi_mean( jj ) );
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, interpolated_rotamer.chi_sd( jj ) );
	sample.set_prob( interpolated_rotamer.rotamer_probability() );

	return sample;
}


template < Size T, Size N >
DunbrackRotamerSampleData
RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	utility::fixedsizearray1< Real, N > bbs;
	bbs[1] = phi;
	bbs[2] = psi;
	return get_rotamer( bbs, rot_ind );
}

template < Size T, Size N >
Size
RotamericSingleResidueDunbrackLibrary< T, N >::nchi() const
{
	return T;
}

template < Size T, Size N >
Size
RotamericSingleResidueDunbrackLibrary< T, N >::nbb() const
{
	return N;
}

template < Size T, Size N >
Size
RotamericSingleResidueDunbrackLibrary< T, N >::n_rotamer_bins() const
{
	return parent::n_possible_rots();
}


/// @details Load interpolated rotamer data into a RotamericData object on the stack
/// so that it can be handed into the enumerate_chi_sets method, which itself relies
/// on the "template method" chisamples_for_rotamer_and_chi.  Enumerate the chi samples,
/// and build rotamers from these chi samples.
///
/// "template method" is a design pattern where a base class calls a polymorphic method
/// that can be overloaded by a derived class, usually in the middle of a function that
/// does a lot of work.  See "Design Patterns," Gamma et al.
template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::build_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	rotamers::RotamerVector & rotamers,
	PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer
) const
{
	DunbrackRotamer< T, N, Real > interpolated_rot( packed_rotamer_2_regular_rotamer( interpolated_rotamer ) );
	interpolated_rot.rotamer_probability() = interpolated_rotamer.rotamer_probability();
	RotamericData< T, N > rotameric_rotamer_building_data( interpolated_rot );

	// now build the chi sets derived from this base rotamer
	utility::vector1< ChiSetOP > chi_set_vector;
	enumerate_chi_sets(
		*concrete_residue, task, existing_residue.seqpos(), buried,
		rotameric_rotamer_building_data,
		extra_chi_steps, chi_set_vector );

	create_rotamers_from_chisets(
		pose, scorefxn, task,
		packer_neighbor_graph, concrete_residue, existing_residue,
		chi_set_vector, rotamers );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::verify_bb_bins(
	utility::fixedsizearray1< Real, N > bbs,
	utility::fixedsizearray1< Size, N > const bb_bin,
	utility::fixedsizearray1< Size, N > const bb_bin_next
) const
{
	for ( Size ii = 1; ii <= N; ++ii ) {
		if ( ( bb_bin[ ii ] < 1 || bb_bin[ ii ] > 36 ) || ( bb_bin_next[ ii ] < 1 || bb_bin_next[ ii ] > 36 ) ) {
			std::cerr << "ERROR: bb " << ii << " bin out of range: " <<
				aa() << " " << bbs[ ii ] << " " << bb_bin[ ii ] << " " << bb_bin_next[ ii ] <<  std::endl;
			utility_exit();
		}
	}
}

/// @details once a list of chi samples has been enumerated, this function
/// instantiates Residue objectes and give them the correct geometry.
template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::create_rotamers_from_chisets(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< ChiSetOP > const & chi_set_vector,
	rotamers::RotamerVector & rotamers
) const
{
	using namespace utility;
	pack::task::ResidueLevelTask const & rtask( task.residue_task( existing_residue.seqpos() ) );

	// construct real rotamers
	for ( vector1< ChiSetOP >::const_iterator chi_set( chi_set_vector.begin() );
			chi_set != chi_set_vector.end(); ++chi_set ) {

		conformation::ResidueOP rotamer = conformation::ResidueFactory::create_residue(
			*concrete_residue, existing_residue, pose.conformation(), rtask.preserve_c_beta() );
		for ( Size jj = 1; jj <= (*chi_set)->chi.size(); ++jj ) rotamer->set_chi( jj, (*chi_set)->chi[ jj ] );

		// apply an operation (or a filter) to this rotamer at build time
		bool reject(false);
		for ( pack::rotamer_set::RotamerOperations::const_iterator
				op( rtask.rotamer_operations().begin() );
				op != rtask.rotamer_operations().end(); ++op ) {
			reject |= ! (**op)( rotamer, pose, scorefxn, rtask, packer_neighbor_graph, *chi_set );
		}
		if ( !reject ) rotamers.push_back( rotamer );
	}
}

template< Size T, Size N >
template< class P >
DunbrackRotamer< T, N, P >
RotamericSingleResidueDunbrackLibrary< T, N >::packed_rotamer_2_regular_rotamer(
	PackedDunbrackRotamer< T, N, P > const & packedrot
) const
{
	DunbrackRotamer< T, N, P > dunrot;
	for ( Size ii = 1; ii <= T; ++ii ) {
		dunrot.chi_mean( ii ) = packedrot.chi_mean( ii );
		dunrot.chi_sd(   ii ) = packedrot.chi_sd( ii );
		dunrot.rotwell(  ii ) = packed_rotno_2_rotwell( packedrot.packed_rotno() )[ ii ];
	}
	dunrot.rotamer_probability() = packedrot.rotamer_probability();
	return dunrot;
}


template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::enumerate_chi_sets(
	chemical::ResidueType const & rsd_type,
	pack::task::PackerTask const & task,
	Size const seqpos,
	bool buried,
	RotamericData< T, N > const & rotamer_data,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	utility::vector1< ChiSetOP > & chi_set_vector
) const
{
	Real const base_probability( rotamer_data.rotamer().rotamer_probability() );
	using namespace utility;
	Size const nchi( rsd_type.nchi() );

	vector1< vector1< Real > > total_chi( nchi );
	vector1< vector1< int  > > total_rot( nchi );
	vector1< vector1< Real > > total_ex_steps( nchi ); // remember which extra chi steps were used for each rotamer
	vector1< vector1< Real > > chisample_prob( nchi ); // let derived class define sample probability

	pack::task::ResidueLevelTask const & rtask( task.residue_task( seqpos ) );
	/// List all chi samples for each chi
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		chisamples_for_rotamer_and_chi(
			rsd_type, rtask, buried, ii, rotamer_data, extra_chi_steps[ ii ],
			total_chi[ ii ], total_rot[ ii ], total_ex_steps[ ii ], chisample_prob[ ii ]
		);
	}

	// now convert from total_chi to chisets....

	Size exchi_product = 1;
	utility::vector1< Size > lex_sizes( total_chi.size() );
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		lex_sizes[ ii ] = total_chi[ ii ].size();
		exchi_product *= total_chi[ ii ].size();
	}
	chi_set_vector.reserve( exchi_product );

	// previously named min_extrachi_rot_prob in rosetta++
	Real const minimum_probability_for_extra_rotamers = 0.001; // make this a virtual function call...

	bool first_rotamer( true );
	for ( utility::LexicographicalIterator lex( lex_sizes ); ! lex.at_end(); ++lex ) {
		Real chi_set_prob = base_probability;

		runtime_assert( nchi <=  lex.size() );
		runtime_assert( nchi <=  chisample_prob.size() );

		for ( Size ii = 1; ii <= nchi; ++ii ) {
			if ( chisample_prob[ ii ].size() >= lex[ ii ] ) chi_set_prob *= chisample_prob[ ii ][ lex[ ii ] ];
		}

		if ( ! first_rotamer && chi_set_prob < minimum_probability_for_extra_rotamers ) continue;

		first_rotamer = false;

		ChiSetOP chi_set( new ChiSet( nchi ) );
		for ( Size ii = 1; ii <= nchi; ++ii ) {
			chi_set->chi[ ii ] = total_chi[ ii ][ lex[ ii ] ];
			chi_set->ex_chi_steps[ ii ] = total_ex_steps[ ii ][ lex[ ii ] ];
			chi_set->rot[ ii ] = total_rot[ ii ][ lex[ ii ] ];
		}
		chi_set_vector.push_back( chi_set );
	}
}


template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::chisamples_for_rotamer_and_chi(
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

	//debug_assert( dynamic_cast< RotamericData const & > ( rotamer_data ) );
	//RotamericData const & rotameric_data( static_cast< RotamericData const & > ( rotamer_data ) );

	// setting this value to zero disables the small-standdev filter -- is that what we want?
	DunbrackReal const min_extrachi_sd( 0.0 );

	if ( chi_index > T ) {
		// chi angle not present in dunbrack library, eg. ser,thr,tyr Hchi

		utility::vector1< std::pair< Real, Real > > const & chi_rotamers( rsd_type.chi_rotamers( chi_index ) );
		if ( ! chi_rotamers.empty() ) {
			// use the rotamer info encoded in the residue
			for ( Size j = 1; j <= chi_rotamers.size(); ++j ) {
				Real const rot_chi_mean( chi_rotamers[ j ].first );
				Real const rot_chi_sdev( chi_rotamers[ j ].second );
				total_chi.push_back( rot_chi_mean );
				total_ex_steps.push_back( 0. );
				total_rot.push_back( j );
				for ( Size k = 1; k <= extra_steps.size(); ++k ) {
					total_chi.push_back( rot_chi_mean + extra_steps[ k ] * rot_chi_sdev );
					total_ex_steps.push_back( extra_steps[ k ] );
					total_rot.push_back( j );
					chisample_prob.push_back( 1.0 ); // assume perfectly likely rotamer sample?
				}
			}
		} else if ( rsd_type.is_proton_chi( chi_index ) ) {
			pack::task::ExtraRotSample ex_samp_level = rtask.extrachi_sample_level( buried, chi_index, rsd_type );
			bool const skip_extra_proton_chi_samples = ( ex_samp_level == pack::task::NO_EXTRA_CHI_SAMPLES );

			Size const proton_chi_index = rsd_type.chi_2_proton_chi( chi_index );
			utility::vector1< Real > const & samples = rsd_type.proton_chi_samples( proton_chi_index );
			utility::vector1< Real > const & extra_samples = rsd_type.proton_chi_extra_samples( proton_chi_index );
			for ( Size ii = 1; ii <= samples.size(); ++ii ) {
				total_chi.push_back( samples[ ii ] );
				total_ex_steps.push_back( 0.0 );
				total_rot.push_back( 1 );
				chisample_prob.push_back( 1.0 );

				if ( skip_extra_proton_chi_samples ) continue;

				for ( Size jj = 1; jj <= extra_samples.size(); ++jj ) {
					total_chi.push_back( samples[ ii ] + extra_samples[ jj ] );
					total_ex_steps.push_back( 0.0 );
					total_rot.push_back( 1 );
					chisample_prob.push_back( 1.0 );

					total_chi.push_back( samples[ ii ] - extra_samples[ jj ] );
					total_ex_steps.push_back( 0.0 );
					total_rot.push_back( 1 );
					chisample_prob.push_back( 1.0 );
				}
			}

		} else {
			// just use the chi of the platonic residue -- this could be made faster, but chi is not part of rsdtype intrfc
			// Is this code ever executed?  It's not going to be executed for any amino acid...
			Real const icoor_chi( numeric::dihedral(
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[ 1 ] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[ 2 ] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[ 3 ] ).ideal_xyz(),
				rsd_type.atom( rsd_type.chi_atoms( chi_index )[ 4 ] ).ideal_xyz()
				));
			total_chi.push_back( icoor_chi );
			total_ex_steps.push_back( 0. );
			total_rot.push_back( 1 );
			chisample_prob.push_back( 1.0 );
		}
	} else {
		// use the dunbrack data
		total_chi.push_back( rotamer_data.rotamer().chi_mean( chi_index ) );
		total_ex_steps.push_back( 0. );
		total_rot.push_back( rotamer_data.rotamer().rotwell( chi_index ) );
		chisample_prob.push_back( 1.0 );

		// ctsa - eliminate extra chi angles with insignificant sd's
		if ( rotamer_data.rotamer().chi_sd( chi_index ) <= min_extrachi_sd ) return;

		for ( Size k = 1; k <= extra_steps.size(); ++k ) {
			total_chi.push_back( rotamer_data.rotamer().chi_mean( chi_index )
				+ extra_steps[k] * rotamer_data.rotamer().chi_sd( chi_index ) );
			total_ex_steps.push_back( extra_steps[k] );
			total_rot.push_back( rotamer_data.rotamer().rotwell( chi_index ) );
			/// prob = exp( -( x - xmean )**2 / 2sd**2). Since we're sampling at intrvals of sd,
			/// prob = exp( -( step * sd )**2 / 2sd**2 ) = exp( - step**2/2 ).
			/// of course, a true gaussian probability would be scaled by 1/sd(2p)**0.5...
			chisample_prob.push_back( std::exp( -0.5 * (extra_steps[ k ]*extra_steps[ k ])) );
		}
	}
}

//XRW_B_T1
/*
template < Size T, Size N >
SingleResidueRotamerLibraryOP
RotamericSingleResidueDunbrackLibrary< T, N >::coarsify(coarse::Translator const & map) const
{
utility_exit_with_message("Unimplemented!");
return 0;
}
*/
//XRW_E_T1

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	utility_exit_with_message("Unimplemented!");
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::write_to_binary( utility::io::ozstream & out ) const
{
	using namespace boost;

	parent::write_to_binary( out );

	/// 1. rotamers_
	{
		Size const ntotalrot = product( N_PHIPSI_BINS ) * parent::n_packed_rots();
		Size const ntotalchi = ntotalrot * T;
		///a. means
		DunbrackReal * rotamer_means = new DunbrackReal[ ntotalchi ];
		// b. standard deviations
		DunbrackReal * rotamer_stdvs = new DunbrackReal[ ntotalchi ];
		// c. rotamer probabilities
		DunbrackReal * rotamer_probs = new DunbrackReal[ ntotalrot ];
		// d. bicubic polynomial parameters
		utility::vector1< DunbrackReal * > rotamer_negln_n_derivs( 1 << N );
		for ( Size initi = 1; initi <= ( 1 << N ); ++initi ) {
			rotamer_negln_n_derivs[ initi ] = new DunbrackReal[ ntotalrot ];
		}

		// e. packed rotamer numbers
		DunbrackReal * packed_rotnos = new DunbrackReal[ ntotalrot ];

		Size count_chi( 0 ), count_rots( 0 );
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );

				for ( Size ll = 1; ll <= T; ++ll ) {
					rotamer_means[ count_chi ] = rotamers_( bb_rot_index, ii ).chi_mean( ll );
					rotamer_stdvs[ count_chi ] = rotamers_( bb_rot_index, ii ).chi_sd( ll );
					++count_chi;
				}
				rotamer_probs[ count_rots ] = rotamers_( bb_rot_index, ii ).rotamer_probability();
				for ( Size deriv_i = 1; deriv_i <= ( 1 << N ); ++deriv_i ) {
					rotamer_negln_n_derivs[ deriv_i ][ count_rots ] = rotamers_( bb_rot_index, ii ).n_derivs()[ deriv_i ];
				}

				packed_rotnos[ count_rots ] = rotamers_( bb_rot_index, ii ).packed_rotno();
				++count_rots;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}
		out.write( (char*) rotamer_means, ntotalchi * sizeof( DunbrackReal ));
		out.write( (char*) rotamer_stdvs, ntotalchi * sizeof( DunbrackReal ));
		out.write( (char*) rotamer_probs, ntotalrot * sizeof( DunbrackReal ));
		for ( Size writei = 1; writei <= ( 1 << N ); ++writei ) {
			out.write( (char*) rotamer_negln_n_derivs[ writei ], ntotalrot * sizeof( DunbrackReal ));
		}

		out.write( (char*) packed_rotnos, ntotalrot * sizeof( DunbrackReal ));
		delete [] rotamer_means;
		delete [] rotamer_stdvs;
		delete [] rotamer_probs;
		for ( Size writei = 1; writei <= ( 1 << N ); ++writei ) {
			delete[] rotamer_negln_n_derivs[ writei ];
		}
		delete [] packed_rotnos;
	}

	/// 2. packed_rotno_to_sorted_rotno_
	{
		Size const ntotalpackedrots = product( N_PHIPSI_BINS ) * parent::n_packed_rots();
		boost::int32_t * packed_rotno_2_sorted_rotno = new boost::int32_t[ ntotalpackedrots ];
		Size count( 0 );
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				packed_rotno_2_sorted_rotno[ count ] = packed_rotno_2_sorted_rotno_( make_index( N, N_PHIPSI_BINS, bb_bin ), ii );
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}
		out.write( (char*) packed_rotno_2_sorted_rotno, ntotalpackedrots * sizeof( boost::int32_t ) );
		delete [] packed_rotno_2_sorted_rotno;
	}
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::read_from_binary( utility::io::izstream & in )
{
	parent::read_from_binary( in );
	Size n_rot_bins = product( N_PHIPSI_BINS );
	/// 1. rotamers_
	{
		Size const ntotalrot = n_rot_bins * parent::n_packed_rots();
		Size const ntotalchi = ntotalrot * T;
		rotamers_.dimension( n_rot_bins, parent::n_packed_rots(), PackedDunbrackRotamer< T, N >() );

		///a. means
		DunbrackReal * rotamer_means = new DunbrackReal[ ntotalchi ];
		// b. standard deviations
		DunbrackReal * rotamer_stdvs = new DunbrackReal[ ntotalchi ];
		// c. rotamer probabilities
		DunbrackReal * rotamer_probs = new DunbrackReal[ ntotalrot ];
		// d. bicubic polynomial parameters
		utility::vector1< DunbrackReal * > rotamer_negln_n_derivs( 1 << N );

		in.read( (char*) rotamer_means, ntotalchi * sizeof( DunbrackReal ));
		in.read( (char*) rotamer_stdvs, ntotalchi * sizeof( DunbrackReal ));
		in.read( (char*) rotamer_probs, ntotalrot * sizeof( DunbrackReal ));
		for ( Size initi = 1; initi <= ( 1 << N ); ++initi ) {
			rotamer_negln_n_derivs[ initi ] = new DunbrackReal[ ntotalrot ];
			in.read( (char*) rotamer_negln_n_derivs[ initi ], ntotalrot * sizeof( DunbrackReal ));
		}
		// e. packed rotamer numbers
		DunbrackReal * packed_rotnos = new DunbrackReal[ ntotalrot ];

		in.read( (char*) packed_rotnos, ntotalrot * sizeof( DunbrackReal ));

		Size count_chi( 0 ), count_rots( 0 );
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );
				for ( Size ll = 1; ll <= T; ++ll ) {
					rotamers_( bb_rot_index, ii ).chi_mean( ll ) = rotamer_means[ count_chi ];
					rotamers_( bb_rot_index, ii ).chi_sd(  ll ) = rotamer_stdvs[ count_chi ];
					++count_chi;
				}
				rotamers_( bb_rot_index, ii ).rotamer_probability() = rotamer_probs[ count_rots ];
				for ( Size deriv_i = 1; deriv_i <= ( 1 << N ); ++deriv_i ) {
					rotamers_( bb_rot_index, ii ).n_derivs()[ deriv_i ] = rotamer_negln_n_derivs[ deriv_i ][ count_rots ];
				}
				rotamers_( bb_rot_index, ii ).packed_rotno( ( Size ) packed_rotnos[ count_rots ] );

				++count_rots;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}

		delete [] rotamer_means;
		delete [] rotamer_stdvs;
		delete [] rotamer_probs;
		for ( Size deli = 1; deli <= ( 1 << N ); ++deli ) {
			delete [] rotamer_negln_n_derivs[ deli ];
		}
		delete [] packed_rotnos;
	}

	/// 2. packed_rotno_to_sorted_rotno_
	{
		Size const ntotalpackedrots = n_rot_bins * parent::n_packed_rots();
		boost::int32_t * packed_rotno_2_sorted_rotno = new boost::int32_t[ ntotalpackedrots ];
		in.read( (char*) packed_rotno_2_sorted_rotno, ntotalpackedrots * sizeof( boost::int32_t ) );
		Size count( 0 );
		packed_rotno_2_sorted_rotno_.dimension( n_rot_bins, parent::n_packed_rots() );
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				packed_rotno_2_sorted_rotno_( make_index( N, N_PHIPSI_BINS, bb_bin ), ii ) = packed_rotno_2_sorted_rotno[ count ];
				++count;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}
		delete [] packed_rotno_2_sorted_rotno;
	}

	/// Entropy setup once reading is finished
	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
		setup_entropy_correction();
	}
}

/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
/// Values tested should parallel those used in the read_from_binary() function.
template < Size T, Size N >
bool
RotamericSingleResidueDunbrackLibrary< T, N >::operator ==( SingleResidueRotamerLibrary const & rhs) const {
	static thread_local basic::Tracer TR( "core.pack.dunbrack.RotamericSingleResidueDunbrackLibrary" );

	// Raw pointer okay, we're just using it to check for conversion
	RotamericSingleResidueDunbrackLibrary< T, N > const * ptr( dynamic_cast< RotamericSingleResidueDunbrackLibrary< T, N > const * > ( &rhs ) );
	if ( ptr == 0 ) {
		TR << "In comparison operator: right-hand side is not a matching RotamericSingleResidueDunbrackLibrary." << std::endl;
		return false;
	}
	RotamericSingleResidueDunbrackLibrary< T, N > const & other( dynamic_cast< RotamericSingleResidueDunbrackLibrary< T, N > const & > ( rhs ) );

	bool equal( true );

	if ( ! parent::operator==( rhs ) ) {
		//TR.Debug << "In RotamericSingleResidueDunbrackLibrary< T >::operator== : Parent comparsion returns false, stopping equality test." << std::endl;
		//return false;
		equal = false;
	}


	/// 1. rotamers_
	if (  parent::n_packed_rots() != other.n_packed_rots() ) {
		return false;
	}
	debug_assert( parent::n_packed_rots() == other.n_packed_rots() );
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
		utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
		for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
		Size p = 1;
		while ( bb_bin[ N + 1 ] == 1 ) {
			Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );

			PackedDunbrackRotamer< T, N > const & this_rot( rotamers_( bb_rot_index, ii ) );
			PackedDunbrackRotamer< T, N > const & other_rot( other.rotamers_( bb_rot_index, ii ) );
			for ( Size ll = 1; ll <= T; ++ll ) {
				if ( ! numeric::equal_by_epsilon( this_rot.chi_mean( ll ), other_rot.chi_mean( ll ), ANGLE_DELTA ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
						<< " rotamers chi mean " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " " << ll << " - "
						<< this_rot.chi_mean( ll ) << " vs. " <<  other_rot.chi_mean( ll ) << std::endl;
					equal = false;
				}
				if ( ! numeric::equal_by_epsilon( this_rot.chi_sd( ll ), other_rot.chi_sd( ll ), ANGLE_DELTA ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
						<< " rotamers chi sd " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " " << ll << " - "
						<< this_rot.chi_sd( ll ) << " vs. " <<  other_rot.chi_sd( ll ) << std::endl;
					equal = false;
				}
			}
			if ( ! numeric::equal_by_epsilon( this_rot.rotamer_probability(), other_rot.rotamer_probability(), PROB_DELTA ) ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
					<< " rotamer probability " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " - "
					<< this_rot.rotamer_probability() << " vs. " <<  other_rot.rotamer_probability() << std::endl;
				equal = false;
			}
			for ( Size n_deriv_index = 1; n_deriv_index <= (1 << N); ++n_deriv_index ) {
				if ( ! numeric::equal_by_epsilon( this_rot.n_derivs()[n_deriv_index], other_rot.n_derivs()[n_deriv_index], ENERGY_DELTA ) ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
						<< " deriv " << n_deriv_index << " " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " - "
						<< this_rot.n_derivs()[n_deriv_index] << " vs. " <<  other_rot.n_derivs()[n_deriv_index] << std::endl;
					equal = false;
				}
			}
			if ( this_rot.packed_rotno() != other_rot.packed_rotno() ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
					<< " packed rotno " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " - "
					<< this_rot.packed_rotno() << " vs. " <<  other_rot.packed_rotno() << std::endl;
				equal = false;
			}

			bb_bin[ 1 ]++;
			while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
				bb_bin[ p ] = 1;
				bb_bin[ ++p ]++;
				if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
			}
		}
	}

	/// 2. packed_rotno_to_sorted_rotno_
	for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
		utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
		utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
		for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
		Size p = 1;
		while ( bb_bin[ N + 1 ] == 1 ) {
			Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );
			if ( packed_rotno_2_sorted_rotno_( bb_rot_index, ii ) != packed_rotno_2_sorted_rotno_( bb_rot_index, ii ) ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa( aa() )
					<< " packed_rotno_2_sorted_rotno " << bb_bin[1] << " " << bb_bin[2] << " " << ii << " - "
					<< packed_rotno_2_sorted_rotno_( bb_rot_index, ii ) << " vs. " << other.packed_rotno_2_sorted_rotno_( bb_rot_index, ii ) << std::endl;
				equal = false;
			}

			bb_bin[ 1 ]++;
			while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
				bb_bin[ p ] = 1;
				bb_bin[ ++p ]++;
				if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
			}
		}
	}

	/// Other variables ( ShanonEntropy_, S_dsecophi_, S_dsecopsi_, S_dsecophipsi_) are set up from the loaded data
	/// automatically - they should be equivalent if all other data is equivalent.

	return equal;
}

// Entropy
template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::setup_entropy_correction()
{
	// Iter again in order to reference ShanonEntropy
	for ( Size dimi = 1; dimi <= ( 1 << N ); ++dimi ) {
		ShanonEntropy_n_derivs_[ dimi ].dimension( product( N_PHIPSI_BINS ) );
	}

	// Get entropy values first
	utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
	utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
	for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
	Size p = 1;
	while ( bb_bin[ N + 1 ] == 1 ) {

		Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );

		// Initialize
		for ( Size deriv_i = 1; deriv_i <= ( 1 << N ); ++deriv_i ) {
			ShanonEntropy_n_derivs_[ deriv_i ]( bb_rot_index ) = 0.0;
		}

		// Divide probability by psum in order to make probability normalized
		// especially for semi-rotameric amino acids
		Real psum( 0.0 );
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
			psum += rotamers_( bb_rot_index, ii ).rotamer_probability();
		}

		// The values actually stored are positive sign ( == negative entropy )
		// which corresponds to Free energy contribution by Entropy
		for ( Size ii = 1; ii <= parent::n_packed_rots(); ++ii ) {
			ShanonEntropy_n_derivs_[ 1 ]( bb_rot_index ) += -rotamers_( bb_rot_index, ii ).n_derivs()[ 1 ] * rotamers_( bb_rot_index, ii ).rotamer_probability() / psum;
		}

		bb_bin[ 1 ]++;
		while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
			bb_bin[ p ] = 1;
			bb_bin[ ++p ]++;
			if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
		}
	}

	/*
	std::cout << std::setw(3) << "AA";
	std::cout << " " << std::setw(6) << "Phi" << " " << std::setw(6) << "Psi";
	std::cout << " " << std::setw(8) << "Prot1" << " " << std::setw(8) << "Prot2" << " " << std::setw(8) << "Prot3";
	std::cout << " " << std::setw(8) << "S" << " " << std::setw(8) << "d2Sdphi2";
	std::cout << " " << std::setw(8) << "d2Sdpsi2" << " " << std::setw(8) << "d4Sdp2p2";
	std::cout << std::endl;
	*/

	// Derivatives

	for ( Size i = 1; i <= N; ++i ) bb_bin       = 1;
	for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
	p = 1;
	while ( bb_bin[ N + 1 ] == 1 ) {

		// AMW
		// Entropy derivative setup for each bb bin (used to be over jj, kk
		// phi/psi indices, but now you construct bb_bin from one plus the
		// first N indices and call make_index)

		// Set up denominators for each bb
		utility::fixedsizearray1< Real, N > dbb2( 0 );
		for ( Size ii = 1; ii <= N; ++ii ) {
			dbb2[ ii ] = 4.0 * PHIPSI_BINRANGE[ ii ] * PHIPSI_BINRANGE[ ii ];
		}

		// Make "this index"
		Size bin_index = make_index( N, N_PHIPSI_BINS, bb_bin );

		// Make an array of "next indices" for each BB
		utility::fixedsizearray1< Real, N > next_bb_index( 0 );
		for ( Size ii = 1; ii <= N; ++ii ) {
			Size old_number = bb_bin[ ii ];
			if ( bb_bin[ ii ] == N_PHIPSI_BINS[ ii ] ) {
				bb_bin[ ii ] = 1;
			} else {
				++bb_bin[ ii ];
			}
			next_bb_index[ ii ] = make_index( N, N_PHIPSI_BINS, bb_bin );
			bb_bin[ ii ] = old_number;
		}

		// Make an array of "previous indices" for each BB
		utility::fixedsizearray1< Real, N > prev_bb_index( 0 );
		for ( Size ii = 1; ii <= N; ++ii ) {
			Size old_number = bb_bin[ ii ];
			if ( bb_bin[ ii ] == 1 ) {
				bb_bin[ ii ] = N_PHIPSI_BINS[ ii ];
			} else {
				--bb_bin[ ii ];
			}
			prev_bb_index[ ii ] = make_index( N, N_PHIPSI_BINS, bb_bin );
			bb_bin[ ii ] = old_number;
		}

		// Indices all acquired... now to calculate.
		// Note that ShanonEntropy_n_derivs_ has already been set for 1
		// because that is the value
		for ( Size deriv_i = 2; deriv_i <= ( 1 << N ); ++deriv_i ) {
			ShanonEntropy_n_derivs_[ deriv_i ]( bin_index ) = 1;
		}

		for ( Size deriv_i = 2; deriv_i <= ( 1 << N ); ++deriv_i ) {
			// No clever pre-computation available
			// (In the 2D case, you just wrote down dsecophi and dsecopsi
			// and dsecophipsi was the product)
			// No ordering guarantee in 3+ dimensions

			for ( Size bbi = 1; bbi <= N; ++bbi ) {
				// You get a factor of (next_val - prev_val) / dbb2[bbi] if
				// the derivative bit is set for this bb angle
				// Otherwise, factor of 1
				if ( ( deriv_i - 1 ) & ( 1 << ( N - bbi ) ) ) {
					ShanonEntropy_n_derivs_[ deriv_i ]( bin_index ) *=
						( ShanonEntropy_n_derivs_[ 1 ]( next_bb_index[ bbi ] )
						- ShanonEntropy_n_derivs_[ 1 ]( prev_bb_index[ bbi ] ) )
						/ dbb2[ bbi ];
				}
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

/// @details Returns the three letter string of the next amino acid specified in the
/// input library.
template < Size T, Size N >
std::string
RotamericSingleResidueDunbrackLibrary< T, N >::read_from_file(
	utility::io::izstream & infile,
	bool first_line_three_letter_code_already_read
)
{
	Size const num_bb_bins = product( N_PHIPSI_BINS );

	std::string const my_name( chemical::name_from_aa( aa() ) );
	std::string next_name; // empty string to start

	/// Read all the data from the first phi/psi bin and keep it temporarily.
	/// Note all the rotamer wells encountered along the way; then declare_all_rotwells_encountered,
	/// allocate the big tables, transfer the data in the temporary arrays into the big table,
	typename utility::vector1< DunbrackRotamer< T, N > > first_phipsibin_data;
	first_phipsibin_data.reserve( n_possible_rots() );

	utility::fixedsizearray1< Real, N > bb;
	DunbrackReal probability(0.0);
	//MaximCode
	//The natural logarithm is the base-e logarithm: the inverse of the natural exponential function (exp).
	DunbrackReal minusLogProbability(-1e-10);

	utility::fixedsizearray1< Size, N > bb_bin;
	utility::fixedsizearray1< Size, N > last_bb_bin;
	Size count( 0 );
	bool very_first_rotamer( true );
	bool finished_first_phipsi_bin( false );
	Size count_in_this_phipsi_bin( 1 );
	utility::vector1< Size > rotwell( DUNBRACK_MAX_SCTOR, 0 );
	utility::vector1< DunbrackReal > chimean( DUNBRACK_MAX_SCTOR, 0.0 );
	utility::vector1< DunbrackReal > chisd( DUNBRACK_MAX_SCTOR, 0.0 );
	std::string three_letter_code;
	Size count_read_rotamers(0);

	while ( infile ) { // if the file should suddenly end, we'd be in trouble.

		/// 1. peek at the line; if it starts with #, skip to the next line.
		char first_char = infile.peek();
		if ( first_char == '#' ) {
			std::string line;
			infile.getline( line );
			continue;
		}

		/// 2. Read the line.  Format is:
		/// a.   three-letter-code,
		/// b.   phi,
		/// c.   psi,
		/// d.   count
		/// e..h r1..r4
		/// i.   probability
		/// j..m chimean1..chimean4
		/// n..q chisd1..chisd4

		if ( first_line_three_letter_code_already_read ) {
			/// The last library to read from this file already ate my three letter code;
			/// skip directly to reading phi/psi values...
			first_line_three_letter_code_already_read = false;
		} else {
			infile >> three_letter_code;
			if ( three_letter_code != my_name ) {
				next_name = three_letter_code;

				initialize_bicubic_splines();
				/// Entropy setup once reading is finished
				if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
					setup_entropy_correction();
				}

				return next_name;
			}/// else, we're still reading data for the intended amino acid, so let's continue...
		}

		for ( Size ii = 1; ii <= N; ++ii ) {
			infile >> bb[ ii ];
		}

		infile >> count;
		infile >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		infile >> probability;
		//MaximCode
		if ( basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
				&& basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable] ) {
			infile >> minusLogProbability;
		}
		infile >> chimean[ 1 ] >> chimean[ 2 ] >> chimean[ 3 ] >> chimean[ 4 ];
		infile >> chisd[ 1 ] >> chisd[ 2 ] >> chisd[ 3 ] >> chisd[ 4 ];

		bool cont = false;
		for ( Size ii = 1; ii <= N; ++ii ) if ( bb[ ii ] == -180 ) cont = true; // duplicated data...
		if ( cont ) continue;

		++count_read_rotamers;

		/// AVOID INF!
		for ( Size ii = 1; ii <= T; ++ii ) if ( chisd[ ii ] == 0.0 ) chisd[ ii ] = 5; // bogus!

		if ( probability <= 1e-6 ) {
			probability = 1e-6;
			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
			/// AMW -- Changing this to be a <= criterion and to 1e-6 to actually match the resolution
			/// of the library
		}

		get_bb_bins( bb, bb_bin );

		if ( finished_first_phipsi_bin ) {
			Size const packed_rotno = rotwell_2_packed_rotno( rotwell );

			if ( packed_rotno < 1 || packed_rotno > parent::n_packed_rots() ) {
				std::cerr << "ERROR in converting rotwell to packed rotno: ";
				for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwell[ ii ];
				std::cerr << " " << packed_rotno << std::endl;
				utility_exit();
			}
			PackedDunbrackRotamer< T, N > rotamer( chimean, chisd, probability, packed_rotno );

			for ( Size ii = 1; ii <= N; ++ii ) {
				if ( bb_bin[ ii ] != last_bb_bin[ ii ] ) {
					count_in_this_phipsi_bin = 1; break;
				}
			}
			if ( count_in_this_phipsi_bin > parent::n_packed_rots() ) break;

			Size bb_bin_index = make_index( N, N_PHIPSI_BINS, bb_bin );
			rotamers_( bb_bin_index, count_in_this_phipsi_bin ) = rotamer;
			packed_rotno_2_sorted_rotno_( bb_bin_index, packed_rotno ) = count_in_this_phipsi_bin;
			++count_in_this_phipsi_bin;
			for ( Size ii = 1; ii <= N; ++ii ) last_bb_bin[ ii ] = bb_bin[ ii ];

		} else {
			bool notdone = false;
			for ( Size ii = 1; ii <= N; ++ii ) {
				if ( bb_bin[ ii ] != last_bb_bin[ ii ] ) {
					notdone = true; break;
				}
			}
			if ( !very_first_rotamer && notdone ) {
				// We have now read all rotamers from this phi/psi bin.
				// 1. Declare all existing rotwells encountered
				// 2. Allocate space for rotamer data
				// 3. Store data from first rotwell.
				// 4. Store the data from the rotamer we just read.

				// 1.
				declare_all_existing_rotwells_encountered();
				// 2.
				rotamers_.dimension( num_bb_bins, n_packed_rots() );//, PackedDunbrackRotamer< T, N >() );
				packed_rotno_2_sorted_rotno_.dimension( num_bb_bins, n_packed_rots() );
				// 3.
				utility::vector1< Size > first_bin_rotwell( DUNBRACK_MAX_SCTOR, 0 );
				for ( Size ii = 1; ii <= first_phipsibin_data.size(); ++ii ) {
					for ( Size jj = 1; jj <= T; ++jj ) first_bin_rotwell[ jj ] = first_phipsibin_data[ ii ].rotwell( jj );
					Size const packed_rotno = rotwell_2_packed_rotno( first_bin_rotwell );
					if ( packed_rotno < 1 || packed_rotno > parent::n_packed_rots() ) {
						std::cerr << "ERROR in converting rotwell to packed rotno: ";
						for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwell[ ii ];
						std::cerr << " " << packed_rotno << std::endl;
						utility_exit();
					}

					Size bb_bin_index = make_index( N, N_PHIPSI_BINS, last_bb_bin );
					rotamers_( bb_bin_index, ii ) = PackedDunbrackRotamer< T, N >( first_phipsibin_data[ ii ], packed_rotno );
					rotamers_( bb_bin_index, ii ).rotamer_probability() = first_phipsibin_data[ ii ].rotamer_probability();

					packed_rotno_2_sorted_rotno_( bb_bin_index, packed_rotno ) = ii;
				}

				// 4.
				debug_assert( count_in_this_phipsi_bin == 1 );
				Size const packed_rotno = rotwell_2_packed_rotno( rotwell );
				PackedDunbrackRotamer< T, N > rotamer( chimean, chisd, probability, packed_rotno );
				Size bb_bin_index = make_index( N, N_PHIPSI_BINS, bb_bin );

				rotamers_( bb_bin_index, count_in_this_phipsi_bin ) = rotamer;
				packed_rotno_2_sorted_rotno_( bb_bin_index, packed_rotno ) = count_in_this_phipsi_bin;
				++count_in_this_phipsi_bin;
				for ( Size bbi = 1; bbi <= N; ++bbi ) last_bb_bin[ bbi ] = bb_bin[ bbi ];

				finished_first_phipsi_bin = true;
			} else {
				very_first_rotamer = false;
				mark_rotwell_exists( rotwell );
				first_phipsibin_data.push_back( DunbrackRotamer< T, N >( chimean, chisd, probability, rotwell ) );
				first_phipsibin_data[first_phipsibin_data.size()].rotamer_probability() = probability;
				for ( Size bbi = 1; bbi <= N; ++bbi ) last_bb_bin[ bbi ] = bb_bin[ bbi ];
			}
		}
	}

	initialize_bicubic_splines();
	/// Entropy setup once reading is finished
	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ) {
		setup_entropy_correction();
	}

	return next_name; // if we arived here, we got to the end of the file, so return the empty string.
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::initialize_bicubic_splines()
{
	// AMW: we have to treat the n_bb = 2 base case differently
	// because of how differently BicubicSpline works as an object
	// i.e. higher order splines solve problems by reducing to BicubicSpline as a base case
	// there may be a more elegant method to integrate the two but I don't know it.
	if ( N == 2 ) {
		for ( Size ii = 1; ii <= n_packed_rots(); ++ii ) {
			using namespace numeric;
			using namespace numeric::interpolation::spline;

			BicubicSpline rotEspline;
			MathMatrix< Real > energy_vals( N_PHIPSI_BINS[1], N_PHIPSI_BINS[2] );
			for ( Size jj = 0; jj < N_PHIPSI_BINS[1]; ++jj ) {
				for ( Size kk = 0; kk < N_PHIPSI_BINS[2]; ++kk ) {
					Size jp1kp1 = jj * N_PHIPSI_BINS[2] + kk+1;
					Size sortedrotno = packed_rotno_2_sorted_rotno_( jp1kp1, ii );
					PackedDunbrackRotamer< T, N > & rot = rotamers_( jp1kp1, sortedrotno );
					DunbrackReal rotamer_energy = -std::log( rot.rotamer_probability() );
					rot.n_derivs()[ 1 ] = rotamer_energy;
					energy_vals( jj, kk ) = rotamer_energy;
				}
			}
			BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
			Real start_vals[2] = {0.0, 0.0}; // grid is placed on the ten-degree marks
			Real deltas[2] = {10.0, 10.0}; // grid is 10 degrees wide
			bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
			std::pair< Real, Real > unused[2];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			rotEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
			for ( Size jj = 0; jj < N_PHIPSI_BINS[1]; ++jj ) {
				for ( Size kk = 0; kk < N_PHIPSI_BINS[2]; ++kk ) {
					Size jp1kp1 = jj * N_PHIPSI_BINS[2] + kk+1;
					Size sortedrotno = packed_rotno_2_sorted_rotno_( jp1kp1, ii );
					PackedDunbrackRotamer< T, N > & rot = rotamers_( jp1kp1, sortedrotno );
					rot.n_derivs()[ 3 ] = rotEspline.get_dsecox()(  jj, kk );
					rot.n_derivs()[ 2 ] = rotEspline.get_dsecoy()(  jj, kk );
					rot.n_derivs()[ 4 ] = rotEspline.get_dsecoxy()( jj, kk );
				}
			}
		}
		// might as well hard-code case 3 to use tricubic splines
	} else if ( N == 3 ) {
		for ( Size ii = 1; ii <= n_packed_rots(); ++ii ) {
			using namespace numeric;
			using namespace numeric::interpolation::spline;

			TricubicSpline rotEspline;
			MathTensor< Real > energy_vals( N_PHIPSI_BINS[1], N_PHIPSI_BINS[2], N_PHIPSI_BINS[3] );
			for ( Size jj = 0; jj < N_PHIPSI_BINS[1]; ++jj ) {
				for ( Size kk = 0; kk < N_PHIPSI_BINS[2]; ++kk ) {
					for ( Size ll = 0; ll < N_PHIPSI_BINS[3]; ++ll ) {
						Size jp1kp1lp1 = jj * N_PHIPSI_BINS[3] * N_PHIPSI_BINS[2] + kk * N_PHIPSI_BINS[3] + ll+1;
						Size sortedrotno = packed_rotno_2_sorted_rotno_( jp1kp1lp1, ii );
						PackedDunbrackRotamer< T, N > & rot = rotamers_( jp1kp1lp1, sortedrotno );
						rot.rotamer_probability() = rotamers_( jp1kp1lp1, sortedrotno ).rotamer_probability();
						//if ( rot.rotamer_probability() <= 1e-6 ) rot.rotamer_probability() = 1e-6;
						DunbrackReal rotamer_energy = -std::log( rot.rotamer_probability() );
						rot.n_derivs()[ 1 ] = rotamer_energy;
						energy_vals( jj, kk, ll ) = rotamer_energy;
					}
				}
			}
			BorderFlag periodic_boundary[3] = { e_Periodic, e_Periodic, e_Periodic };
			Real start_vals[3] = {0.0, 0.0, 0.0}; // grid is placed on the ten-degree marks
			Real deltas[3] = {PHIPSI_BINRANGE[1], PHIPSI_BINRANGE[2], PHIPSI_BINRANGE[3]}; // grid is 10 degrees wide
			bool lincont[3] = {false,false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
			std::pair< Real, Real > unused[3];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			unused[2] = std::make_pair( 0.0, 0.0 );
			rotEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
			for ( Size jj = 0; jj < N_PHIPSI_BINS[1]; ++jj ) {
				for ( Size kk = 0; kk < N_PHIPSI_BINS[2]; ++kk ) {
					for ( Size ll = 0; ll < N_PHIPSI_BINS[3]; ++ll ) {
						Size jp1kp1lp1 = jj * N_PHIPSI_BINS[3] * N_PHIPSI_BINS[2] + kk * N_PHIPSI_BINS[3] + ll+1;
						Size sortedrotno = packed_rotno_2_sorted_rotno_( jp1kp1lp1, ii );
						PackedDunbrackRotamer< T, N > & rot = rotamers_( jp1kp1lp1, sortedrotno );
						rot.n_derivs()[ 2 ] = rotEspline.get_dsecoz()(  jj, kk, ll  );
						rot.n_derivs()[ 3 ] = rotEspline.get_dsecoy()(  jj, kk, ll );
						rot.n_derivs()[ 4 ] = rotEspline.get_dsecoyz()( jj, kk, ll  );
						rot.n_derivs()[ 5 ] = rotEspline.get_dsecox()(  jj, kk, ll );
						rot.n_derivs()[ 6 ] = rotEspline.get_dsecoxz()(  jj, kk, ll  );
						rot.n_derivs()[ 7 ] = rotEspline.get_dsecoxy()(  jj, kk, ll );
						rot.n_derivs()[ 8 ] = rotEspline.get_dsecoxyz()( jj, kk, ll  );
					}
				}
			}
		}
	} else {
		for ( Size ii = 1; ii <= n_packed_rots(); ++ii ) {
			using namespace numeric;
			using namespace numeric::interpolation::spline;

			PolycubicSpline rotEspline( N );
			utility::vector1< Size > n_dims( N );
			for ( Size i = 1; i <= N; ++i ) n_dims[ i ] = N_PHIPSI_BINS[i];
			MathNTensor< Real > energy_vals( n_dims, 0.0 );

			utility::fixedsizearray1< Size, (N+1) > bb_bin( 1 );
			utility::fixedsizearray1< Size, (N+1) > bb_bin_maxes( 1 );
			for ( Size i = 1; i <= N; ++i ) bb_bin_maxes = N_PHIPSI_BINS[i];
			Size p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );
				Size sortedrotno = packed_rotno_2_sorted_rotno_( bb_rot_index, ii );
				PackedDunbrackRotamer< T, N > & rot = rotamers_( bb_rot_index, sortedrotno );
				DunbrackReal rotamer_energy = -std::log( rot.rotamer_probability() );
				rot.n_derivs()[ 1 ] = rotamer_energy;
				utility::vector1< Size > indices;
				for ( Size i = 1; i <= N; ++i ) indices.push_back( bb_bin[ i ] - 1 );
				energy_vals( indices ) = rotamer_energy;

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}

			utility::vector1< BorderFlag > periodic_boundary( N, e_Periodic );
			utility::vector1< Real > start_vals( N, 0.0 );
			utility::vector1< Real > deltas( N, 10.0 );
			for ( Size s = 1; s <= N; ++s ) deltas[ s ] = PHIPSI_BINRANGE[ s ];
			utility::vector1< bool > lincont( N, false );
			utility::vector1< std::pair< Real, Real > > unused( N, std::make_pair( 0.0, 0.0 ) );
			rotEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );

			for ( Size clear_i = 1; clear_i <= N+1; ++clear_i ) bb_bin[ clear_i ] = 1;
			p = 1;
			while ( bb_bin[ N + 1 ] == 1 ) {

				Size bb_rot_index = make_index( N, N_PHIPSI_BINS, bb_bin );
				Size sortedrotno = packed_rotno_2_sorted_rotno_( bb_rot_index, ii );
				PackedDunbrackRotamer< T, N > & rot = rotamers_( bb_rot_index, sortedrotno );

				utility::vector1< Size > indices;
				for ( Size i = 1; i <= N; ++i ) indices.push_back( bb_bin[ i ] - 1 );
				for ( Size i = 1; i <= ( 1 << N ); ++i ) rot.n_derivs()[ i ] = rotEspline.get_deriv( i )( indices );

				bb_bin[ 1 ]++;
				while ( bb_bin[ p ] == bb_bin_maxes[p]+1 ) {
					bb_bin[ p ] = 1;
					bb_bin[ ++p ]++;
					if ( bb_bin[ p ] != bb_bin_maxes[p]+1 ) p = 1;
				}
			}
		}
	}
}

/// @details called only if the library is actually an RSRDL<T> object.  Derived classes
/// should not call this function or recurse.  Accounts for the statically allocated data
/// that's part of this class.
template < Size T, Size N >
Size RotamericSingleResidueDunbrackLibrary< T, N >::memory_usage_static() const
{
	return sizeof( RotamericSingleResidueDunbrackLibrary< T, N > );
}

/// @details Measures the amount of dynamically allocated data in this class.  Must recurse to parent
/// to count parent's dynamically allocated data.
template < Size T, Size N >
Size RotamericSingleResidueDunbrackLibrary< T, N >::memory_usage_dynamic() const
{
	Size total = parent::memory_usage_dynamic(); // recurse to parent.
	total += rotamers_.size() * sizeof( PackedDunbrackRotamer< T, N > );
	total += packed_rotno_2_sorted_rotno_.size() * sizeof( Size ); // could make these shorts or chars!
	//total += max_rotprob_.size() * sizeof( DunbrackReal );
	return total;
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi(
	ChiVector const & chi,
	RotVector & rot
) const
{
	Size4 rot4;
	get_rotamer_from_chi_static( chi, rot4 );
	rot.resize( chi.size() );
	for ( Size ii = 1; ii <= T; ++ii ) rot[ ii ] = rot4[ ii ];
	for ( Size ii = T+1; ii <= chi.size(); ++ii ) rot[ ii ] = 0;
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi_static(
	ChiVector const & chi,
	Size4 & rot
) const
{
	Real4 chi4;
	for ( Size ii = 1; ii <= T; ++ii ) chi4[ ii ] = chi[ ii ];
	get_rotamer_from_chi_static( chi4, rot );
}

template < Size T, Size N >
void
RotamericSingleResidueDunbrackLibrary< T, N >::get_rotamer_from_chi_static(
	Real4 const & chi,
	Size4 & rot
) const
{
	core::chemical::AA aa2=aa();
	if ( core::chemical::is_canonical_D_aa( aa2 ) ) aa2 = core::chemical::get_L_equivalent( aa2 );

	if ( dun02() ) { rotamer_from_chi_02( chi, aa2, T, rot ); return; }

	debug_assert( chi.size() >= T );

	/// compiler will unroll this loop
	for ( Size ii = 1; ii <= T; ++ii ) rot[ ii ] = bin_rotameric_chi( basic::periodic_range( chi[ ii ], 360 ), ii );
}


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_TMPL_HH
