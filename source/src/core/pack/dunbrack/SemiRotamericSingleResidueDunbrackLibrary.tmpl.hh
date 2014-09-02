// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh
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
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/conformation/Residue.hh>
#include <basic/basic.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.functions.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/MathTensor.hh>
#include <numeric/interpolation/spline/TricubicSpline.hh>

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
#include <cassert>
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

//using namespace ObjexxFCL;

template < Size T >
SemiRotamericSingleResidueDunbrackLibrary< T >::SemiRotamericSingleResidueDunbrackLibrary(
	chemical::AA const aa_in,
	bool const backbone_independent_scoring,         // true uses less memory
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
{

}

template < Size T >
SemiRotamericSingleResidueDunbrackLibrary< T >::~SemiRotamericSingleResidueDunbrackLibrary() {}

template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	if ( bbind_nrchi_scoring_ ) {
		return rotamer_energy_deriv_bbind( rsd, scratch, false );
	} else {
		return rotamer_energy_deriv_bbdep( rsd, scratch, false );
	}
}


template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	if ( bbind_nrchi_scoring_ ) {
		return rotamer_energy_deriv_bbind( rsd, scratch, true );
	} else {
		return rotamer_energy_deriv_bbdep( rsd, scratch, true );
	}
}

template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::rotamer_energy_deriv_bbdep(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const
{
	assert( ! bbind_nrchi_scoring_ );

	//Multiplier for D-amino acids:
	const core::Real d_multiplier = (core::chemical::is_canonical_D_aa(rsd.aa()) ) ? -1.0 : 1.0;

	parent::eval_rotameric_energy_deriv( rsd, scratch, eval_deriv ); //This has been updated to allow scoring of D-amino acids.

	Real chidevpen_score( 0.0 );
	for ( Size ii = 1; ii <= T; ++ii ) {
		chidevpen_score += scratch.chidevpen()[ ii ];
	}

	Real nrchi_score, dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi;
	nrchi_score = bbdep_nrchi_score( rsd, scratch,
		dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi ); //Handles D-amino acids; derivatives WRT phi, psi, and chi are inverted iff D-amino acid.

	if ( nrchi_score != nrchi_score ) { //NaN check
		nrchi_score = 0;
		std::cerr << "NaN in SemiRot: " << rsd.seqpos() << " " << rsd.name() << std::endl;
	}

	// Corrections for Shanon Entropy
	if( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ){
		nrchi_score  += scratch.entropy();
	}

	scratch.fa_dun_tot() = chidevpen_score + nrchi_score;
	scratch.fa_dun_rot() = 0;
	scratch.fa_dun_semi() = nrchi_score;

	//scratch.fa_dun_dev() = chidevpensum; // keep original chidev score...

	//std::cout << "SRSRDL bbdep score: " << rsd.name() << " " <<
	//	rsd.mainchain_torsion( 1 ) << " " << rsd.mainchain_torsion( 2 ) <<
	//	" " << rsd.chi()[ T + 1 ] << ": " << nrchi_score << " " << chidevpen_score << " " << std::endl;

	if ( ! eval_deriv ) return chidevpen_score + nrchi_score;

	/// sum derivatives.
	Real3 & dE_dbb(  scratch.dE_dbb() );
	Real3 & dE_dbb_dev(  scratch.dE_dbb_dev() );
	Real3 & dE_dbb_rot(  scratch.dE_dbb_rot() );
	Real3 & dE_dbb_semi(  scratch.dE_dbb_semi() );
	Real4 & dE_dchi( scratch.dE_dchi() );
	Real4 & dE_dchi_dev( scratch.dE_dchi_dev() );
	Real4 & dE_dchi_semi( scratch.dE_dchi_semi() );

	std::fill( dE_dchi.begin(), dE_dchi.end(), 0.0 );
	std::fill( dE_dchi_dev.begin(), dE_dchi_dev.end(), 0.0 );
	std::fill( dE_dchi_semi.begin(), dE_dchi_semi.end(), 0.0 );
	std::fill( dE_dbb.begin(), dE_dbb.end(), 0.0 );
	std::fill( dE_dbb_dev.begin(), dE_dbb_dev.end(), 0.0 );
	std::fill( dE_dbb_rot.begin(), dE_dbb_rot.end(), 0.0 );
	std::fill( dE_dbb_semi.begin(), dE_dbb_semi.end(), 0.0 );

	for ( Size i=1; i<= DUNBRACK_MAX_BBTOR; ++i ) {
		dE_dbb[ i ] = d_multiplier*scratch.dchidevpen_dbb()[ i ];
		dE_dbb_dev[ i ] = d_multiplier*scratch.dchidevpen_dbb()[ i ];
		//dE_dbb_rot[ i ] = scratch.dneglnrotprob_dbb()[ i ];

		// Correction for entropy
		if( basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ] ){
			dE_dbb[ i ]      += d_multiplier * scratch.dentropy_dbb()[ i ];
			dE_dbb_semi[ i ] += d_multiplier * scratch.dentropy_dbb()[ i ];
		}
	}
	dE_dbb[ RotamerLibraryScratchSpace::AA_PHI_INDEX ] += d_multiplier*dnrchiscore_dphi;
	dE_dbb[ RotamerLibraryScratchSpace::AA_PSI_INDEX ] += d_multiplier*dnrchiscore_dpsi;
	dE_dbb_semi[ RotamerLibraryScratchSpace::AA_PHI_INDEX ] = d_multiplier*dnrchiscore_dphi;
	dE_dbb_semi[ RotamerLibraryScratchSpace::AA_PSI_INDEX ] = d_multiplier*dnrchiscore_dpsi;

	for ( Size i=1; i <= T; ++i ) {
		dE_dchi[ i ] = d_multiplier*scratch.dchidevpen_dchi()[ i ];
		dE_dchi_dev[ i ] = d_multiplier*scratch.dchidevpen_dchi()[ i ];
	}
	dE_dchi[ T + 1 ] = d_multiplier*dnrchiscore_dchi;
	dE_dchi_semi[ T + 1 ] = d_multiplier*dnrchiscore_dchi;

	parent::correct_termini_derivatives( rsd, scratch );

	return chidevpen_score + nrchi_score;

}

template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::rotamer_energy_deriv_bbind(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const
{
	assert( bbind_nrchi_scoring_ );

	Real rotameric_score = parent::eval_rotameric_energy_deriv( rsd, scratch, eval_deriv );

	/// The true score is -ln ( p0 ) + nrchi_score + chidev_penalty.
	Real nrchi_score, dnrchiscore_dnrchi;
	nrchi_score = bbind_nrchi_score( rsd, scratch, dnrchiscore_dnrchi );

	scratch.fa_dun_tot() = rotameric_score + nrchi_score;
	//scratch.fa_dun_rot() = 0; // keep original rot score...
	scratch.fa_dun_semi() = nrchi_score;
	//scratch.fa_dun_dev() = chidevpensum; // keep original chidevpen score...

	if ( ! eval_deriv ) return rotameric_score + nrchi_score;

	/// sum derivatives.
	Real3 & dE_dbb(  scratch.dE_dbb() );
	Real3 & dE_dbb_dev(  scratch.dE_dbb_dev() );
	Real3 & dE_dbb_rot(  scratch.dE_dbb_rot() );
	//Real3 & dE_dbb_semi(  scratch.dE_dbb_semi() );
	Real4 & dE_dchi( scratch.dE_dchi() );
	Real4 & dE_dchi_dev( scratch.dE_dchi_dev() );
	Real4 & dE_dchi_semi( scratch.dE_dchi_semi() );
	std::fill( dE_dchi.begin(), dE_dchi.end(), 0.0 );
	std::fill( dE_dchi_dev.begin(), dE_dchi_dev.end(), 0.0 );
	std::fill( dE_dchi_semi.begin(), dE_dchi_semi.end(), 0.0 );

	// p0 - the base probability -- not modified by the chi-dev penalty
	Real const rotprob( scratch.rotprob() );

	/*Size const nbb ( rsd.mainchain_torsions().size() );
	Size const nchi( rsd.nchi() );

	dE_dbb.resize ( nbb );
	dE_dchi.resize( nchi );*/

	Real const invp( ( rotprob == Real( 0.0 ) ) ? 0.0 : -1.0 / rotprob );

	for ( Size i=1; i<= DUNBRACK_MAX_BBTOR; ++i ) {
		dE_dbb[ i ] = invp * scratch.drotprob_dbb()[ i ] + scratch.dchidevpen_dbb()[ i ];
		dE_dbb_dev[ i ] = scratch.dchidevpen_dbb()[ i ];
		dE_dbb_rot[ i ] = invp * scratch.drotprob_dbb()[ i ];
	}

	for ( Size i=1; i<= T; ++i ) {
		dE_dchi[ i ] = scratch.dchidevpen_dchi()[ i ];
		dE_dchi_dev[ i ] = scratch.dchidevpen_dchi()[ i ];
	}
	dE_dchi[ T + 1 ] = dnrchiscore_dnrchi;
	dE_dchi_semi[ T + 1 ] = dnrchiscore_dnrchi;

	parent::correct_termini_derivatives( rsd, scratch );

	return rotameric_score + nrchi_score;
}

/// @details simple interpolation; derivatives discontinuous at grid points
template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::bbind_nrchi_score(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	Real & dnrchi_score_dnrchi
) const
{
	assert( bbind_nrchi_scoring_ );

	Size const packed_rotno( grandparent::rotwell_2_packed_rotno( scratch.rotwell() ));

	Real nrchi = rsd.chi( T + 1 );
	Size nrchi_bin, nrchi_bin_next;
	Real nrchi_alpha;

	get_bbind_nrchi_bin( nrchi, nrchi_bin, nrchi_bin_next, nrchi_alpha );

	Real const nrchi_beta = 1 - nrchi_alpha;
	Real const score_lower_bin = bbind_non_rotameric_chi_scores_( nrchi_bin,      packed_rotno );
	Real const score_upper_bin = bbind_non_rotameric_chi_scores_( nrchi_bin_next, packed_rotno );

	dnrchi_score_dnrchi = ( score_upper_bin - score_lower_bin ) / bbind_nrchi_binsize_;
	return nrchi_beta * score_lower_bin + nrchi_alpha * score_upper_bin;
}


/// @brief Trilinear interpolation.  Derivatives discontinuous at edge planes.
/// Correction: it seems that this was updated at some point to use tricubic interpolation.
/// Note: This function has been updated to handle D-amino acids properly.
template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::bbdep_nrchi_score(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	Real & dnrchi_score_dnrchi,
	Real & dnrchi_score_dphi,
	Real & dnrchi_score_dpsi
) const
{
	assert( ! bbind_nrchi_scoring_ );

	Size const packed_rotno( parent::rotwell_2_packed_rotno( scratch.rotwell() ));

	Real nrchi = rsd.chi( T + 1 );
	core::Real d_multiplier = 1.0; //A multiplier for D-amino acid angles.  -1.0 if D-amino acid, 1.0 otherwise.
	if(core::chemical::is_canonical_D_aa(rsd.aa())) {
		nrchi*=-1.0; //Invert chi if this is a D-amino acid
		d_multiplier = -1.0;
	}
	Size nrchi_bin, nrchi_bin_next;
	Real nrchi_alpha;
	get_bbdep_nrchi_bin( nrchi, nrchi_bin, nrchi_bin_next, nrchi_alpha );

	Real phi( d_multiplier * parent::get_phi_from_rsd( rsd ) ); //Inverted iff D-amino acid
	Real psi( d_multiplier * parent::get_psi_from_rsd( rsd ) ); //Inverted iff D-amino acid

	// get phi/psi bins
	Size phibin, phibin_next, psibin, psibin_next;
	Real phi_alpha, psi_alpha;
	// the alpha fraction being the fraction ( range [0,1)) of the way from the
	// lower grid point to the upper grid point
	parent::get_phipsi_bins(
		phi, psi, phibin, psibin,
		phibin_next, psibin_next,
		phi_alpha, psi_alpha );

	assert( phibin >= 1 && phibin <= parent::N_PHIPSI_BINS );
	assert( psibin >= 1 && psibin <= parent::N_PHIPSI_BINS );

	BBDepScoreInterpData const & d000( bbdep_nrc_interpdata_[ packed_rotno ]( phibin     , psibin     , nrchi_bin      ));
	BBDepScoreInterpData const & d001( bbdep_nrc_interpdata_[ packed_rotno ]( phibin     , psibin     , nrchi_bin_next ));
	BBDepScoreInterpData const & d010( bbdep_nrc_interpdata_[ packed_rotno ]( phibin     , psibin_next, nrchi_bin      ));
	BBDepScoreInterpData const & d011( bbdep_nrc_interpdata_[ packed_rotno ]( phibin     , psibin_next, nrchi_bin_next ));
	BBDepScoreInterpData const & d100( bbdep_nrc_interpdata_[ packed_rotno ]( phibin_next, psibin     , nrchi_bin      ));
	BBDepScoreInterpData const & d101( bbdep_nrc_interpdata_[ packed_rotno ]( phibin_next, psibin     , nrchi_bin_next ));
	BBDepScoreInterpData const & d110( bbdep_nrc_interpdata_[ packed_rotno ]( phibin_next, psibin_next, nrchi_bin      ));
	BBDepScoreInterpData const & d111( bbdep_nrc_interpdata_[ packed_rotno ]( phibin_next, psibin_next, nrchi_bin_next ));

	Real interpolated_energy(0.0);

	tricubic_interpolation(
		d000.value_, d000.dsecox_, d000.dsecoy_, d000.dsecoz_, d000.dsecoxy_, d000.dsecoxz_, d000.dsecoyz_, d000.dsecoxyz_,
		d001.value_, d001.dsecox_, d001.dsecoy_, d001.dsecoz_, d001.dsecoxy_, d001.dsecoxz_, d001.dsecoyz_, d001.dsecoxyz_,
		d010.value_, d010.dsecox_, d010.dsecoy_, d010.dsecoz_, d010.dsecoxy_, d010.dsecoxz_, d010.dsecoyz_, d010.dsecoxyz_,
		d011.value_, d011.dsecox_, d011.dsecoy_, d011.dsecoz_, d011.dsecoxy_, d011.dsecoxz_, d011.dsecoyz_, d011.dsecoxyz_,
		d100.value_, d100.dsecox_, d100.dsecoy_, d100.dsecoz_, d100.dsecoxy_, d100.dsecoxz_, d100.dsecoyz_, d100.dsecoxyz_,
		d101.value_, d101.dsecox_, d101.dsecoy_, d101.dsecoz_, d101.dsecoxy_, d101.dsecoxz_, d101.dsecoyz_, d101.dsecoxyz_,
		d110.value_, d110.dsecox_, d110.dsecoy_, d110.dsecoz_, d110.dsecoxy_, d110.dsecoxz_, d110.dsecoyz_, d110.dsecoxyz_,
		d111.value_, d111.dsecox_, d111.dsecoy_, d111.dsecoz_, d111.dsecoxy_, d111.dsecoxz_, d111.dsecoyz_, d111.dsecoxyz_,
		phi_alpha, psi_alpha, nrchi_alpha,
		10, 10, bbdep_nrchi_binsize_,
		interpolated_energy,
		dnrchi_score_dphi,
		dnrchi_score_dpsi,
		dnrchi_score_dnrchi );

	return interpolated_energy;

	/* Trilinear interpolation scheme below
	/// I would have defined this trilinear interpolation in the numeric library but for
	/// two problems: the bin width is not uniform for all three dimensions, and the
	/// data structure "F" that the interpolation library wants to read from must be indexed
	/// from 0, which, while possible with FArrays is distinctly not fun.  Indices that work
	/// for the other 1-based tables will not work for a 0-based table; code to compute indices
	/// needs to be nearly duplicated.  The better solution would be to adapt the numeric
	/// libraries to use 1-based indexing, since nothing in Rosetta uses 0-based indexing.

	Real const
		flll( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin     , psibin     , nrchi_bin      )),
		fllu( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin     , psibin     , nrchi_bin_next )),
		flul( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin     , psibin_next, nrchi_bin      )),
		fluu( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin     , psibin_next, nrchi_bin_next )),
		full( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin_next, psibin     , nrchi_bin      )),
		fulu( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin_next, psibin     , nrchi_bin_next )),
		fuul( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin_next, psibin_next, nrchi_bin      )),
		fuuu( bbdep_non_rotameric_chi_scores_[ packed_rotno ]( phibin_next, psibin_next, nrchi_bin_next ));

  //if ( flll > 1e16 ) {std::cout << "inf: flll " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( fllu > 1e16 ) {std::cout << "inf: fllu " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( flul > 1e16 ) {std::cout << "inf: flul " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( fluu > 1e16 ) {std::cout << "inf: fluu " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( full > 1e16 ) {std::cout << "inf: full " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( fulu > 1e16 ) {std::cout << "inf: fulu " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( fuul > 1e16 ) {std::cout << "inf: fuul " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }
  //if ( fuuu > 1e16 ) {std::cout << "inf: fuuu " <<  phibin << " " << psibin << " " << phibin_next << " " << psibin_next << " " << nrchi_bin << " " << nrchi_bin_next << std::endl; }


	Real const a1 = phi_alpha, a2 = psi_alpha, a3 = nrchi_alpha;
	Real const b1 = 1 - a1,    b2 = 1 - a2,    b3 = 1 - a3;

	dnrchi_score_dphi = (
		( b2 * b3 * ( full - flll ) ) + ( a2 * b3 * ( fuul - flul ) ) +
		( b2 * a3 * ( fulu - fllu ) ) + ( a2 * a3 * ( fuuu - fluu ) ) ) / parent::PHIPSI_BINRANGE;
	dnrchi_score_dpsi = (
		( b1 * b3 * ( flul - flll ) ) + ( a1 * b3 * ( fuul - full ) ) +
		( b1 * a3 * ( fluu - fllu ) ) + ( a1 * a3 * ( fuuu - fulu ) ) ) / parent::PHIPSI_BINRANGE;
	dnrchi_score_dnrchi = (
		( b1 * b2 * ( fllu - flll ) ) + ( a1 * b2 * ( fulu - full ) ) +
		( b1 * a2 * ( fluu - flul ) ) + ( a1 * a2 * ( fuuu - fuul ) ) ) / bbdep_nrchi_binsize_;

	return
	 ( b1 * b2 * b3 * flll ) +
	 ( a1 * b2 * b3 * full ) +
	 ( b1 * a2 * b3 * flul ) +
	 ( a1 * a2 * b3 * fuul ) +
	 ( b1 * b2 * a3 * fllu ) +
	 ( a1 * b2 * a3 * fulu ) +
	 ( b1 * a2 * a3 * fluu ) +
	 ( a1 * a2 * a3 * fuuu ); */

}

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::best_rotamer_energy(
	conformation::Residue const & rsd,
	bool curr_rotamer_only,
	RotamerLibraryScratchSpace & scratch
) const
{
		assert( rsd.nchi() == T+1 );
    Real nrchi_score( 0 );
  if ( curr_rotamer_only ) {
	      Real dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi;
			// Unused variable, but since thescratch is a non-const reference, it seems wise to continue calling the method
			//Real rotameric_score =
			parent::eval_rotameric_energy_deriv( rsd, scratch, false);

        nrchi_score = bbdep_nrchi_score( rsd, scratch, dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi );
				core::conformation::Residue rsd_copy (rsd);
				utility::vector1< Real > rsd_chi=rsd.chi();

				for ( Size jj = 0; jj <= bbdep_nrchi_nbins_; ++jj ) {
							rsd_chi[rsd_copy.nchi()]=nrchi_lower_angle_+bbdep_nrchi_binsize_*jj;
							rsd_copy.chi(rsd_chi);
							parent::eval_rotameric_energy_deriv( rsd_copy, scratch, false);
							Real tmp_nrchi_score=bbdep_nrchi_score( rsd_copy, scratch, dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi );
							if ( tmp_nrchi_score < nrchi_score)
										nrchi_score=tmp_nrchi_score;
					}

  } else {
				core::pack::dunbrack::SingleResidueRotamerLibraryCAP rotlib = RotamerLibrary::get_instance().get_rsd_library( rsd.type() );
				core::pack::dunbrack::SingleResidueDunbrackLibraryCAP dunlib( static_cast< SingleResidueDunbrackLibrary const * > ( rotlib() ));

				Real const phi( parent::get_phi_from_rsd( rsd ) );
    		Real const psi( parent::get_psi_from_rsd( rsd ) );

				utility::vector1< DunbrackRotamerSampleData > rotamer_samples=dunlib->get_all_rotamer_samples( phi, psi);
				//this could be smarter since the T+1 position of the sc_torsions are not used

				Real dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi;
				parent::eval_rotameric_energy_deriv( rsd, scratch, false);
        nrchi_score = bbdep_nrchi_score( rsd, scratch, dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi );
				Real tmp_nrchi_score;
				//search the space of terminal chiT
				core::conformation::Residue rsd_copy (rsd);
				utility::vector1< Real > rsd_chi=rsd.chi();

				for ( Size jj = 1; jj <= rotamer_samples.size(); ++jj ) {

						for ( Size ii = 1; ii <= T; ++ii ) {
              		rsd_chi[ii]=rotamer_samples[jj].chi_mean()[ii];
						}

        		for ( Size kk = 0; kk <= bbdep_nrchi_nbins_; ++kk ) {
              		rsd_chi[rsd_copy.nchi()]=nrchi_lower_angle_+bbdep_nrchi_binsize_*kk;
              		rsd_copy.chi(rsd_chi);
              		parent::eval_rotameric_energy_deriv( rsd_copy, scratch, false);
              		tmp_nrchi_score=bbdep_nrchi_score( rsd_copy, scratch, dnrchiscore_dchi, dnrchiscore_dphi, dnrchiscore_dpsi );
              		if ( tmp_nrchi_score < nrchi_score)
                    		nrchi_score=tmp_nrchi_score;
          		}

				}
	}

	return nrchi_score;
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::assign_random_rotamer_with_bias(
	conformation::Residue const & rsd,
	pose::Pose const & /*pose*/,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	if ( bbind_nrchi_sampling_ ) {
		assign_random_rotamer_with_bias_bbind( rsd, scratch, RG, new_chi_angles, perturb_from_rotamer_center );
	} else {
		assign_random_rotamer_with_bias_bbdep( rsd, scratch, RG, new_chi_angles, perturb_from_rotamer_center );
	}

}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::assign_random_rotamer_with_bias_bbind(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	// Parent will pick a rotwell for the rotameric chi.
	Size packed_rotno( 0 );
	parent::assign_random_rotamer( rsd, scratch, RG, new_chi_angles, perturb_from_rotamer_center, packed_rotno );

	/// Now, given the packed_rotno for the rotameric chi, pick a non-rotameric chi sample.
	Real random_prob = RG.uniform();
	Size count( 0 ), nrchi_bin( 0 );
	while ( random_prob > 0 ) {
		nrchi_bin = bbind_rotamers_sorted_by_probability_( ++count, packed_rotno );
		random_prob -= bbind_rotamers_to_sample_( nrchi_bin, packed_rotno ).prob_;
		if ( count == n_nrchi_sample_bins_ ) break;
	}

	if ( perturb_from_rotamer_center ) {
		/// Could be a lot of work to bias the sample in this well by the continuous energy distribution.
		/// why not just sample uniformly in this well?

		Real nrchi_prob = RG.uniform();
		Real left( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno ).left_ );
		Real right( bbind_rotamers_to_sample_( nrchi_bin, packed_rotno ).right_ );
		assert( left < right ); // I'm 99% this is true, that bins don't wrap around the center divide.

		new_chi_angles[ T + 1 ] = left + (left-right) * nrchi_prob;
	} else {
		new_chi_angles[ T + 1 ] = bbind_rotamers_to_sample_( nrchi_bin, packed_rotno ).median_;
	}
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::assign_random_rotamer_with_bias_bbdep(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const
{
	Real random_prob = RG.uniform();

	Real const phi( parent::get_phi_from_rsd( rsd ) );
	Real const psi( parent::get_psi_from_rsd( rsd ) );

	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next,phi_alpha, psi_alpha );

	BBDepNRChiSample< Real > interpolated_nrchi_sample;
	Size count( 0 );
	while ( random_prob > 0 ) {
		BBDepNRChiSample<> const & nrchi_sample_00(
			bbdep_rotamers_to_sample_( phibin, psibin, ++count ));

		interpolated_nrchi_sample = interpolate_bbdep_nrchi_sample(
			nrchi_sample_00.packed_rotno_, nrchi_sample_00.nrchi_bin_,
			phibin, psibin, phibin_next, psibin_next,
			phi_alpha, psi_alpha
		);

		random_prob -= interpolated_nrchi_sample.prob_;
		if ( count == bbdep_rotamers_to_sample_.size3() ) break;
	}

	/// Get chimean and chisdves for rotameric chi
	Size packed_rotno = interpolated_nrchi_sample.packed_rotno_;
	PackedDunbrackRotamer< T, Real > interpolated_rotamer;
	parent::interpolate_rotamers( scratch, packed_rotno,
		phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha,
		interpolated_rotamer );

	this->assign_chi_for_interpolated_rotamer( interpolated_rotamer, rsd, RG, new_chi_angles, perturb_from_rotamer_center );

	if ( ! perturb_from_rotamer_center ) {
		new_chi_angles[ T + 1 ] = interpolated_nrchi_sample.nrchi_mean_;
	} else {
		new_chi_angles[ T + 1 ] = interpolated_nrchi_sample.nrchi_mean_ +
			RG.gaussian() * interpolated_nrchi_sample.nrchi_sd_;
	}
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	RotamerVector & rotamers
) const
{
	if ( bbind_nrchi_sampling_ ) {
		fill_rotamer_vector_bbind( pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried, rotamers );
	} else {
		fill_rotamer_vector_bbdep( pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried, rotamers );
	}
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::fill_rotamer_vector_bbind(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	RotamerVector & rotamers
) const
{
	RotamerLibraryScratchSpace scratch;

	//Determine whether this is a D-amino acid:
	core::Real d_multiplier = 1.0;
	if(core::chemical::is_canonical_D_aa( existing_residue.aa() ) ) d_multiplier=-1.0;

	/// Save backbone interpolation data for reuse
	Real phi( d_multiplier * parent::get_phi_from_rsd( existing_residue ) );
	Real psi( d_multiplier * parent::get_psi_from_rsd( existing_residue ) );
	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	utility::vector1< Real > probs_of_next_rotamers( grandparent::n_packed_rots() );
	utility::vector1< Size > next_nrchi_rotamer_to_take( grandparent::n_packed_rots(), 1 );
	typename utility::vector1< PackedDunbrackRotamer< T, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		parent::interpolate_rotamers(
			scratch, ii,
			phibin, psibin,
			phibin_next, psibin_next,phi_alpha, psi_alpha,
			interpolated_rotamers[ ii ] );

		probs_of_next_rotamers[ ii ] = interpolated_rotamers[ ii ].rotamer_probability() *
			bbind_rotamers_to_sample_( bbind_rotamers_sorted_by_probability_( 1, ii ), ii ).prob_;
	}

	Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;

	Real const requisit_probability = buried ? 0.87 : 0.95;
	//grandparent::probability_to_accumulate_while_building_rotamers( buried ); -- 98/95 split generates too many samples
	Real accumulated_probability( 0.0 );

	Size count_rotamers_built = 0;
	while ( accumulated_probability < requisit_probability ) {
		/// Build the most probable rotamer of those remaining
		Size const which_packedrotno_to_build = utility::arg_max( probs_of_next_rotamers );

		Size const count = next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];
		Size const next_count = ++next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];

		Size const nrchi_rotno = bbind_rotamers_sorted_by_probability_( count, which_packedrotno_to_build );
		Size const next_nrchi_rotno = count < n_nrchi_sample_bins_ ?
			bbind_rotamers_sorted_by_probability_( next_count, which_packedrotno_to_build ) : 0 ;

		accumulated_probability += probs_of_next_rotamers[ which_packedrotno_to_build ];
		++count_rotamers_built;

		if ( probs_of_next_rotamers[ which_packedrotno_to_build ] <= 0 ) break; // this should never happen.

		if ( next_nrchi_rotno == 0 ) {
			probs_of_next_rotamers[ which_packedrotno_to_build ] = 0;
		} else {
			probs_of_next_rotamers[ which_packedrotno_to_build ] =
				interpolated_rotamers[ which_packedrotno_to_build ].rotamer_probability() *
				bbind_rotamers_to_sample_(
				next_nrchi_rotno,
				which_packedrotno_to_build ).prob_;
		}

		build_bbind_rotamers(
			pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried,
			interpolated_rotamers[ which_packedrotno_to_build ], nrchi_rotno,
			bbind_rotamers_to_sample_( nrchi_rotno, which_packedrotno_to_build ),
			rotamers);

		if ( count_rotamers_built == max_rots_that_can_be_built ) break;
	}
	//If this is a D-amino acid, rotamers will subsequently need to be inverted.
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::fill_rotamer_vector_bbdep(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	RotamerVector & rotamers
) const
{
	RotamerLibraryScratchSpace scratch;
	utility::vector1< Size > rotamer_has_been_interpolated( grandparent::n_packed_rots(), 0 );
	typename utility::vector1< PackedDunbrackRotamer< T, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	//Determine whether this is a D-amino acid:
	core::Real d_multiplier = 1.0;
	if(core::chemical::is_canonical_D_aa( existing_residue.aa() ) ) d_multiplier=-1.0;

	/// Save backbone interpolation data for reuse
	Real phi( d_multiplier * parent::get_phi_from_rsd( existing_residue ) );
	Real psi( d_multiplier * parent::get_psi_from_rsd( existing_residue ) );
	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	Real const requisit_probability = buried ? 0.95 : 0.87;
	//grandparent::probability_to_accumulate_while_building_rotamers( buried ); -- 98/95 split generates too many samples
	Real accumulated_probability( 0.0 );

	Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	Size count_rotamers_built = 0;
	while ( accumulated_probability < requisit_probability ) {
		++count_rotamers_built;

		BBDepNRChiSample<> const & nrchi_sample_00(
			bbdep_rotamers_to_sample_( phibin, psibin, count_rotamers_built ));

		Size const packed_rotno00 = nrchi_sample_00.packed_rotno_;
		Size const nrchi_bin = nrchi_sample_00.nrchi_bin_;

		Size const
		ind01( bbdep_rotsample_sorted_order_( phibin     , psibin_next, packed_rotno00, nrchi_bin )),
		ind10( bbdep_rotsample_sorted_order_( phibin_next, psibin     , packed_rotno00, nrchi_bin )),
		ind11( bbdep_rotsample_sorted_order_( phibin_next, psibin_next, packed_rotno00, nrchi_bin ));

		BBDepNRChiSample<> const & nrchi_sample_01( bbdep_rotamers_to_sample_( phibin     , psibin_next, ind01 ));
		BBDepNRChiSample<> const & nrchi_sample_10( bbdep_rotamers_to_sample_( phibin_next, psibin     , ind10 ));
		BBDepNRChiSample<> const & nrchi_sample_11( bbdep_rotamers_to_sample_( phibin_next, psibin_next, ind11 ));

		BBDepNRChiSample< Real > const interpolated_nrchi_sample =
			interpolate_bbdep_nrchi_sample(
			nrchi_sample_00, nrchi_sample_01,
			nrchi_sample_10, nrchi_sample_11,
			phi_alpha, psi_alpha
		);

		if ( rotamer_has_been_interpolated[ packed_rotno00 ] == 0 ) {
			/// interpolate the rotameric chi at most once
			rotamer_has_been_interpolated[ packed_rotno00 ] = 1;
			parent::interpolate_rotamers(
				scratch, packed_rotno00,
				phibin, psibin,
				phibin_next, psibin_next,phi_alpha, psi_alpha,
				interpolated_rotamers[ packed_rotno00 ] );
		}

		build_bbdep_rotamers(
			pose, scorefxn, task, packer_neighbor_graph,
			concrete_residue, existing_residue, extra_chi_steps, buried,
			interpolated_rotamers[ packed_rotno00 ], interpolated_nrchi_sample,
			rotamers );

		accumulated_probability += interpolated_nrchi_sample.prob_;

		if ( count_rotamers_built == max_rots_that_can_be_built ) break;
	}
	//If this is a D-amino acid, rotamers will subsequently need to be inverted.
}

template < Size T >
utility::vector1< DunbrackRotamerSampleData >
SemiRotamericSingleResidueDunbrackLibrary< T >::get_all_rotamer_samples(
	Real phi,
	Real psi
) const
{
	if ( bbind_nrchi_sampling_ ) {
		return get_all_rotamer_samples_bbind( phi, psi );
	} else {
		return get_all_rotamer_samples_bbdep( phi, psi );
	}

}

template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::get_probability_for_rotamer(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	if ( bbind_nrchi_sampling_ ) {
		return get_probability_for_rotamer_bbind( phi, psi, rot_ind );
	} else {
		return get_probability_for_rotamer_bbdep( phi, psi, rot_ind );
	}
}

template < Size T >
DunbrackRotamerSampleData
SemiRotamericSingleResidueDunbrackLibrary< T >::get_rotamer(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	if ( bbind_nrchi_sampling_ ) {
		return get_rotamer_bbind( phi, psi, rot_ind );
	} else {
		return get_rotamer_bbdep( phi, psi, rot_ind );
	}
}


template < Size T >
Size
SemiRotamericSingleResidueDunbrackLibrary< T >::nchi() const
{
	return T + 1;
}

template < Size T >
Size
SemiRotamericSingleResidueDunbrackLibrary< T >::n_rotamer_bins() const
{
	return grandparent::n_possible_rots() * n_nrchi_sample_bins_;
}


template < Size T >
utility::vector1< DunbrackRotamerSampleData >
SemiRotamericSingleResidueDunbrackLibrary< T >::get_all_rotamer_samples_bbind(
	Real phi,
	Real psi
) const
{
	RotamerLibraryScratchSpace scratch;

	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	utility::vector1< Real > probs_of_next_rotamers( grandparent::n_packed_rots() );
	utility::vector1< Size > next_nrchi_rotamer_to_take( grandparent::n_packed_rots(), 1 );
	typename utility::vector1< PackedDunbrackRotamer< T, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		parent::interpolate_rotamers(
			scratch, ii,
			phibin, psibin,
			phibin_next, psibin_next,phi_alpha, psi_alpha,
			interpolated_rotamers[ ii ] );

		probs_of_next_rotamers[ ii ] = interpolated_rotamers[ ii ].rotamer_probability() *
			bbind_rotamers_to_sample_( bbind_rotamers_sorted_by_probability_( 1, ii ), ii ).prob_;
	}

	//Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	Size const n_rots = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	for ( Size ii = 1; ii <= n_rots; ++ii ) {
		/// Build the most probable rotamer of those remaining
		Size const which_packedrotno_to_build = utility::arg_max( probs_of_next_rotamers );

		PackedDunbrackRotamer< T, Real > nextrot( interpolated_rotamers[ which_packedrotno_to_build ] );

		Size const count = next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];
		Size const next_count = ++next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];

		Size const nrchi_rotno = bbind_rotamers_sorted_by_probability_( count, which_packedrotno_to_build );
		Size const next_nrchi_rotno = count < n_nrchi_sample_bins_ ?
			bbind_rotamers_sorted_by_probability_( next_count, which_packedrotno_to_build ) : 0 ;

		DunbrackRotamerSampleData sample( true );
		sample.set_nchi( T + 1 );
		sample.set_rotwell( parent::packed_rotno_2_rotwell(nextrot.packed_rotno()) );
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, nextrot.chi_mean( jj ) );
		for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, nextrot.chi_sd( jj ) );
		sample.set_rotwell( T + 1, nrchi_rotno );
		sample.set_chi_mean( T + 1, bbind_rotamers_to_sample_( nrchi_rotno, which_packedrotno_to_build ).median_ );
		sample.set_chi_sd( T + 1, 0.0 ); // bogus value
		sample.set_prob( probs_of_next_rotamers[ which_packedrotno_to_build ] );
		sample.set_nrchi_lower_boundary( bbind_rotamers_to_sample_( nrchi_rotno, which_packedrotno_to_build ).left_ );
		sample.set_nrchi_upper_boundary( bbind_rotamers_to_sample_( nrchi_rotno, which_packedrotno_to_build ).right_ );
		sample.set_nrchi_probability( bbind_rotamers_to_sample_( nrchi_rotno, which_packedrotno_to_build ).prob_ );


		if ( probs_of_next_rotamers[ which_packedrotno_to_build ] <= 0 ) break; // this should never happen.

		if ( next_nrchi_rotno == 0 ) {
			probs_of_next_rotamers[ which_packedrotno_to_build ] = 0;
		} else {
			probs_of_next_rotamers[ which_packedrotno_to_build ] =
				interpolated_rotamers[ which_packedrotno_to_build ].rotamer_probability() *
				bbind_rotamers_to_sample_( next_nrchi_rotno, which_packedrotno_to_build ).prob_;
		}
		all_rots.push_back( sample );

	}
	return all_rots;
}


template < Size T >
utility::vector1< DunbrackRotamerSampleData >
SemiRotamericSingleResidueDunbrackLibrary< T >::get_all_rotamer_samples_bbdep(
	Real phi,
	Real psi
) const
{
	RotamerLibraryScratchSpace scratch;
	utility::vector1< Size > rotamer_has_been_interpolated( grandparent::n_packed_rots(), 0 );
	typename utility::vector1< PackedDunbrackRotamer< T, Real > > interpolated_rotamers( grandparent::n_packed_rots() );

	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	Size const n_rots = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	for ( Size ii = 1; ii <= n_rots; ++ii ) {
		BBDepNRChiSample<> const & nrchi_sample_00(
			bbdep_rotamers_to_sample_( phibin, psibin, ii ));

		Size const packed_rotno00 = nrchi_sample_00.packed_rotno_;
		Size const nrchi_bin = nrchi_sample_00.nrchi_bin_;

		Size const
		ind01( bbdep_rotsample_sorted_order_( phibin     , psibin_next, packed_rotno00, nrchi_bin )),
		ind10( bbdep_rotsample_sorted_order_( phibin_next, psibin     , packed_rotno00, nrchi_bin )),
		ind11( bbdep_rotsample_sorted_order_( phibin_next, psibin_next, packed_rotno00, nrchi_bin ));

		BBDepNRChiSample<> const & nrchi_sample_01( bbdep_rotamers_to_sample_( phibin     , psibin_next, ind01 ));
		BBDepNRChiSample<> const & nrchi_sample_10( bbdep_rotamers_to_sample_( phibin_next, psibin     , ind10 ));
		BBDepNRChiSample<> const & nrchi_sample_11( bbdep_rotamers_to_sample_( phibin_next, psibin_next, ind11 ));

		BBDepNRChiSample< Real > const interpolated_nrchi_sample =
			interpolate_bbdep_nrchi_sample(
			nrchi_sample_00, nrchi_sample_01,
			nrchi_sample_10, nrchi_sample_11,
			phi_alpha, psi_alpha
		);

		if ( rotamer_has_been_interpolated[ packed_rotno00 ] == 0 ) {
			/// interpolate the rotameric chi at most once
			rotamer_has_been_interpolated[ packed_rotno00 ] = 1;
			parent::interpolate_rotamers(
				scratch, packed_rotno00,
				phibin, psibin,
				phibin_next, psibin_next,phi_alpha, psi_alpha,
				interpolated_rotamers[ packed_rotno00 ] );
		}

		PackedDunbrackRotamer< T, Real > nextrot( interpolated_rotamers[ packed_rotno00 ] );

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


/// @details Lookup for bbind is significantly slower than for bbdep
template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::get_probability_for_rotamer_bbind(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	utility::vector1< Real > probs_of_next_rotamers( grandparent::n_packed_rots() );
	utility::vector1< Size > next_nrchi_rotamer_to_take( grandparent::n_packed_rots(), 1 );
	utility::vector1< Real > probability_of_rotameric_chi( grandparent::n_packed_rots() );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {

		PackedDunbrackRotamer< T > const & rot00( parent::rotamers()( phibin, psibin, ii ) );
		Size const packed_rotno = rot00.packed_rotno();

		Size const
			sorted_rotno_01( parent::packed_rotno_2_sorted_rotno()( phibin     , psibin_next, packed_rotno )),
			sorted_rotno_10( parent::packed_rotno_2_sorted_rotno()( phibin_next, psibin     , packed_rotno )),
			sorted_rotno_11( parent::packed_rotno_2_sorted_rotno()( phibin_next, psibin_next, packed_rotno ));

		PackedDunbrackRotamer< T > const &
			rot01( parent::rotamers()( phibin     , psibin_next, sorted_rotno_01 ) ),
			rot10( parent::rotamers()( phibin_next, psibin     , sorted_rotno_10 ) ),
			rot11( parent::rotamers()( phibin_next, psibin_next, sorted_rotno_11 ) );

		Real interpolated_probability, dummy1, dummy2;
		basic::interpolate_bilinear_by_value(
			static_cast< Real >  ( rot00.rotamer_probability()),
			static_cast< Real >  ( rot10.rotamer_probability()),
			static_cast< Real >  ( rot01.rotamer_probability()),
			static_cast< Real >  ( rot11.rotamer_probability()),
			phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, false /*treat_as_angles*/,
			interpolated_probability,
			dummy1, dummy2
		);

		probability_of_rotameric_chi[ ii ]  = interpolated_probability;
		probs_of_next_rotamers[ ii ] = probability_of_rotameric_chi[ ii ] *
			bbind_rotamers_to_sample_( bbind_rotamers_sorted_by_probability_( 1, ii ), ii ).prob_;
	}

	//Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	Size const n_rots = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	Real last_prob( 0 );
	/// Recover the rotamers in order, stop at the rot_ind rotamer.
	for ( Size ii = 1; ii <= rot_ind; ++ii ) {
		/// Build the most probable rotamer of those remaining
		Size const which_packedrotno_to_build = utility::arg_max( probs_of_next_rotamers );

		Size const count = next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];
		Size const next_count = ++next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];

		//Size const nrchi_rotno = bbind_rotamers_sorted_by_probability_( count, which_packedrotno_to_build );
		Size const next_nrchi_rotno = count < n_nrchi_sample_bins_ ?
			bbind_rotamers_sorted_by_probability_( next_count, which_packedrotno_to_build ) : 0 ;

		last_prob = probs_of_next_rotamers[ which_packedrotno_to_build ];

		if ( probs_of_next_rotamers[ which_packedrotno_to_build ] <= 0 ) break; // this should never happen.

		if ( next_nrchi_rotno == 0 ) {
			probs_of_next_rotamers[ which_packedrotno_to_build ] = 0;
		} else {
			probs_of_next_rotamers[ which_packedrotno_to_build ] =
				probability_of_rotameric_chi[ which_packedrotno_to_build ] *
				bbind_rotamers_to_sample_( next_nrchi_rotno, which_packedrotno_to_build ).prob_;
		}

	}
	return last_prob;

}

template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::get_probability_for_rotamer_bbdep(
	Real phi,
	Real psi,
	Size rot_ind
) const
{

	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	BBDepNRChiSample<> const & nrchi_sample_00(
		bbdep_rotamers_to_sample_( phibin, psibin, rot_ind ));

	Size const packed_rotno00 = nrchi_sample_00.packed_rotno_;
	Size const nrchi_bin = nrchi_sample_00.nrchi_bin_;

	Size const
	ind01( bbdep_rotsample_sorted_order_( phibin     , psibin_next, packed_rotno00, nrchi_bin )),
	ind10( bbdep_rotsample_sorted_order_( phibin_next, psibin     , packed_rotno00, nrchi_bin )),
	ind11( bbdep_rotsample_sorted_order_( phibin_next, psibin_next, packed_rotno00, nrchi_bin ));

	BBDepNRChiSample<> const & nrchi_sample_01( bbdep_rotamers_to_sample_( phibin     , psibin_next, ind01 ));
	BBDepNRChiSample<> const & nrchi_sample_10( bbdep_rotamers_to_sample_( phibin_next, psibin     , ind10 ));
	BBDepNRChiSample<> const & nrchi_sample_11( bbdep_rotamers_to_sample_( phibin_next, psibin_next, ind11 ));

	Real nrchi_prob, dummy_dprob_1, dummy_dprob_2;

	basic::interpolate_bilinear_by_value(
		static_cast< Real >  ( nrchi_sample_00.prob_),
		static_cast< Real >  ( nrchi_sample_10.prob_),
		static_cast< Real >  ( nrchi_sample_01.prob_),
		static_cast< Real >  ( nrchi_sample_11.prob_),
		phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, false /*treat_as_angles*/,
		nrchi_prob,
		dummy_dprob_1, dummy_dprob_2
	);

	return nrchi_prob;
}

template < Size T >
DunbrackRotamerSampleData
SemiRotamericSingleResidueDunbrackLibrary< T >::get_rotamer_bbind(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	utility::vector1< Real > probs_of_next_rotamers( grandparent::n_packed_rots() );
	utility::vector1< Size > next_nrchi_rotamer_to_take( grandparent::n_packed_rots(), 1 );
	utility::vector1< Real > probability_of_rotameric_chi( grandparent::n_packed_rots() );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {

		PackedDunbrackRotamer< T > const & rot00( parent::rotamers()( phibin, psibin, ii ) );
		Size const packed_rotno = rot00.packed_rotno();

		Size const
			sorted_rotno_01( parent::packed_rotno_2_sorted_rotno()( phibin     , psibin_next, packed_rotno )),
			sorted_rotno_10( parent::packed_rotno_2_sorted_rotno()( phibin_next, psibin     , packed_rotno )),
			sorted_rotno_11( parent::packed_rotno_2_sorted_rotno()( phibin_next, psibin_next, packed_rotno ));

		PackedDunbrackRotamer< T > const &
			rot01( parent::rotamers()( phibin     , psibin_next, sorted_rotno_01 ) ),
			rot10( parent::rotamers()( phibin_next, psibin     , sorted_rotno_10 ) ),
			rot11( parent::rotamers()( phibin_next, psibin_next, sorted_rotno_11 ) );

		Real interpolated_probability, dummy1, dummy2;
		basic::interpolate_bilinear_by_value(
			static_cast< Real >  ( rot00.rotamer_probability()),
			static_cast< Real >  ( rot10.rotamer_probability()),
			static_cast< Real >  ( rot01.rotamer_probability()),
			static_cast< Real >  ( rot11.rotamer_probability()),
			phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, false /*treat_as_angles*/,
			interpolated_probability,
			dummy1, dummy2
		);

		probability_of_rotameric_chi[ ii ]  = interpolated_probability;
		probs_of_next_rotamers[ ii ] = probability_of_rotameric_chi[ ii ] *
			bbind_rotamers_to_sample_( bbind_rotamers_sorted_by_probability_( 1, ii ), ii ).prob_;
	}

	//Size const max_rots_that_can_be_built = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	Size const n_rots = grandparent::n_packed_rots() * n_nrchi_sample_bins_;
	utility::vector1< DunbrackRotamerSampleData > all_rots;
	all_rots.reserve( n_rots );

	Real last_prob( 0 );
	Size last_nrchi_rotno( 0 );
	Size last_rchi_rotno( 0 );
	/// Recover the rotamers in order, stop at the rot_ind rotamer.
	for ( Size ii = 1; ii <= rot_ind; ++ii ) {
		/// Build the most probable rotamer of those remaining
		Size const which_packedrotno_to_build = utility::arg_max( probs_of_next_rotamers );
		last_prob = probs_of_next_rotamers[ which_packedrotno_to_build ];
		last_rchi_rotno = which_packedrotno_to_build;

		Size const count = next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];
		Size const next_count = ++next_nrchi_rotamer_to_take[ which_packedrotno_to_build ];

		Size const nrchi_rotno = bbind_rotamers_sorted_by_probability_( count, which_packedrotno_to_build );
		last_nrchi_rotno = nrchi_rotno;
		Size const next_nrchi_rotno = count < n_nrchi_sample_bins_ ?
			bbind_rotamers_sorted_by_probability_( next_count, which_packedrotno_to_build ) : 0 ;

		last_prob = probs_of_next_rotamers[ which_packedrotno_to_build ];

		if ( probs_of_next_rotamers[ which_packedrotno_to_build ] <= 0 ) break; // this should never happen.

		if ( next_nrchi_rotno == 0 ) {
			probs_of_next_rotamers[ which_packedrotno_to_build ] = 0;
		} else {
			probs_of_next_rotamers[ which_packedrotno_to_build ] =
				probability_of_rotameric_chi[ which_packedrotno_to_build ] *
				bbind_rotamers_to_sample_( next_nrchi_rotno, which_packedrotno_to_build ).prob_;
		}

	}

	RotamerLibraryScratchSpace scratch;

	PackedDunbrackRotamer< T, Real > rchi_rotamer;
	parent::interpolate_rotamers(
		scratch, last_rchi_rotno,
		phibin, psibin,
		phibin_next, psibin_next,phi_alpha, psi_alpha,
		rchi_rotamer );

	DunbrackRotamerSampleData sample( true );
	sample.set_nchi( T + 1 );
	sample.set_rotwell( parent::packed_rotno_2_rotwell(rchi_rotamer.packed_rotno()) );
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_mean( jj, rchi_rotamer.chi_mean( jj ) );
	for ( Size jj = 1; jj <= T; ++jj ) sample.set_chi_sd( jj, rchi_rotamer.chi_sd( jj ) );
	sample.set_rotwell( T + 1, last_nrchi_rotno );
	sample.set_chi_mean( T + 1, bbind_rotamers_to_sample_( last_nrchi_rotno, last_rchi_rotno ).median_ );
	sample.set_chi_sd( T + 1, 0.0 ); // bogus value
	sample.set_prob( last_prob );
	sample.set_nrchi_lower_boundary( bbind_rotamers_to_sample_( last_nrchi_rotno, last_rchi_rotno ).left_ );
	sample.set_nrchi_upper_boundary( bbind_rotamers_to_sample_( last_nrchi_rotno, last_rchi_rotno ).right_ );
	sample.set_nrchi_probability( bbind_rotamers_to_sample_(    last_nrchi_rotno, last_rchi_rotno ).prob_ );

	return sample;
}

template < Size T >
DunbrackRotamerSampleData
SemiRotamericSingleResidueDunbrackLibrary< T >::get_rotamer_bbdep(
	Real phi,
	Real psi,
	Size rot_ind
) const
{
	RotamerLibraryScratchSpace scratch;

	Size phibin, psibin, phibin_next, psibin_next;
	Real phi_alpha, psi_alpha;
	parent::get_phipsi_bins( phi, psi, phibin, psibin, phibin_next, psibin_next, phi_alpha, psi_alpha );

	BBDepNRChiSample<> const & nrchi_sample_00(
		bbdep_rotamers_to_sample_( phibin, psibin, rot_ind ));

	Size const packed_rotno00 = nrchi_sample_00.packed_rotno_;
	Size const nrchi_bin = nrchi_sample_00.nrchi_bin_;

	Size const
	ind01( bbdep_rotsample_sorted_order_( phibin     , psibin_next, packed_rotno00, nrchi_bin )),
	ind10( bbdep_rotsample_sorted_order_( phibin_next, psibin     , packed_rotno00, nrchi_bin )),
	ind11( bbdep_rotsample_sorted_order_( phibin_next, psibin_next, packed_rotno00, nrchi_bin ));

	BBDepNRChiSample<> const & nrchi_sample_01( bbdep_rotamers_to_sample_( phibin     , psibin_next, ind01 ));
	BBDepNRChiSample<> const & nrchi_sample_10( bbdep_rotamers_to_sample_( phibin_next, psibin     , ind10 ));
	BBDepNRChiSample<> const & nrchi_sample_11( bbdep_rotamers_to_sample_( phibin_next, psibin_next, ind11 ));

	BBDepNRChiSample< Real > const interpolated_nrchi_sample =
		interpolate_bbdep_nrchi_sample(
		nrchi_sample_00, nrchi_sample_01,
		nrchi_sample_10, nrchi_sample_11,
		phi_alpha, psi_alpha
	);

	PackedDunbrackRotamer< T, Real > rot;
	parent::interpolate_rotamers(
		scratch, packed_rotno00,
		phibin, psibin,
		phibin_next, psibin_next,phi_alpha, psi_alpha,
		rot );

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


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::build_bbdep_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	PackedDunbrackRotamer< T, Real > const & interpolated_rotamer,
	BBDepNRChiSample< Real > const interpolated_sample,
	RotamerVector & rotamers
) const
{
	DunbrackRotamer< T, Real > interpolated_rot( parent::packed_rotamer_2_regular_rotamer( interpolated_rotamer ));
	BBDepSemiRotamericData< T > bbdep_rotamer_building_data( interpolated_rot, interpolated_sample );

	// now build the chi sets derived from this base rotamer
	utility::vector1< ChiSetOP > chi_set_vector;
	parent::enumerate_chi_sets(
		*concrete_residue, task, existing_residue.seqpos(), buried,
		bbdep_rotamer_building_data,
		extra_chi_steps, chi_set_vector );

	parent::create_rotamers_from_chisets(
		pose, scorefxn, task,
		packer_neighbor_graph, concrete_residue, existing_residue,
		chi_set_vector, rotamers );
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::build_bbind_rotamers(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	PackedDunbrackRotamer< T, Real > const & interpolated_rotamer,
	Size const nrchi_rotno,
	BBIndNRChiSample<> const & interpolated_sample,
	RotamerVector & rotamers
) const
{
	DunbrackRotamer< T, Real > interpolated_rot( parent::packed_rotamer_2_regular_rotamer( interpolated_rotamer ));
	BBIndSemiRotamericData< T > bbind_rotamer_building_data( interpolated_rot, interpolated_sample, nrchi_rotno );

	// now build the chi sets derived from this base rotamer
	utility::vector1< ChiSetOP > chi_set_vector;
	parent::enumerate_chi_sets(
		*concrete_residue, task, existing_residue.seqpos(), buried,
		bbind_rotamer_building_data,
		extra_chi_steps, chi_set_vector );

	parent::create_rotamers_from_chisets(
		pose, scorefxn, task,
		packer_neighbor_graph, concrete_residue, existing_residue,
		chi_set_vector, rotamers );

}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::chisamples_for_rotamer_and_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T > const & rotamer_data,
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
	} else if ( bbind_nrchi_sampling_ ) {
		bbind_chisamples_for_rotamer_chi( rsd_type, rtask, buried, chi_index, rotamer_data, extra_steps,
			total_chi, total_rot, total_ex_steps, rotsample_prob );
	} else {
		bbdep_chisamples_for_rotamer_chi( rsd_type, rtask, buried, chi_index, rotamer_data, extra_steps,
			total_chi, total_rot, total_ex_steps, rotsample_prob );
	}
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::bbind_chisamples_for_rotamer_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T > const & rotamer_data,
	utility::vector1< Real > const & /*extra_steps*/,
	utility::vector1< Real > & total_chi,
	utility::vector1< int  > & total_rot,
	utility::vector1< Real > & total_ex_steps,
	utility::vector1< Real > & chisample_prob
) const
{
	using namespace pack::task;
	assert( chi_index == T + 1 );

	assert( dynamic_cast< BBIndSemiRotamericData< T > const * > ( & rotamer_data ) );
	BBIndSemiRotamericData< T > const & bbind_rotameric_data(
		static_cast< BBIndSemiRotamericData< T > const & > ( rotamer_data ) );

	BBIndNRChiSample<> const & nrsample = bbind_rotameric_data.bbind_nrchi_sample();

	Real const minimum_angular_distance_for_extra_bbind_rotamer_sample( 10.0 );

	total_chi.push_back( nrsample.median_ );
	total_ex_steps.push_back( 0.0 );
	total_rot.push_back( bbind_rotameric_data.nrchi_bin_id() );
	chisample_prob.push_back( nrsample.prob_ );
	//Real const min_extrachi_sd = 0.0; // ???

	ExtraRotSample ex_samp_level = rtask.extrachi_sample_level( buried, chi_index, &rsd_type );

	if ( ex_samp_level == NO_EXTRA_CHI_SAMPLES ) return;


	if ( std::abs( basic::periodic_range( nrsample.median_ - nrsample.left_, nrchi_periodicity_ )) >
		minimum_angular_distance_for_extra_bbind_rotamer_sample ) {

		/// halfway between the median and the left boundary
		total_chi.push_back( nrsample.median_  - 0.5 * (  nrsample.median_ - nrsample.left_ )    );
		total_ex_steps.push_back( 1.0 ); // not a measure of # stdevs...
		total_rot.push_back( bbind_rotameric_data.nrchi_bin_id() );
		chisample_prob.push_back( nrsample.prob_ );
	}

	if ( std::abs( basic::periodic_range( nrsample.median_ - nrsample.right_, nrchi_periodicity_ )) >
		minimum_angular_distance_for_extra_bbind_rotamer_sample ) {

		/// halfway between the median and the right boundary
		total_chi.push_back( nrsample.median_  + 0.5 * (  nrsample.right_ - nrsample.median_ )    );
		total_ex_steps.push_back( 1.0 ); // not a measure of # stdevs...
		total_rot.push_back( bbind_rotameric_data.nrchi_bin_id() );
		chisample_prob.push_back( nrsample.prob_ );
	}
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::bbdep_chisamples_for_rotamer_chi(
	chemical::ResidueType const & rsd_type,
	pack::task::ResidueLevelTask const & rtask,
	bool buried,
	Size const chi_index,
	RotamericData< T > const & rotamer_data,
	utility::vector1< Real > const & extra_steps,
	utility::vector1< Real > & total_chi,
	utility::vector1< int  > & total_rot,
	utility::vector1< Real > & total_ex_steps,
	utility::vector1< Real > & chisample_prob
) const
{
	using namespace pack::task;
	assert( chi_index == T + 1 );

	assert( dynamic_cast< BBDepSemiRotamericData< T > const * > ( & rotamer_data ) );
	BBDepSemiRotamericData< T > const & bbdep_rotameric_data(
		static_cast< BBDepSemiRotamericData< T > const & > ( rotamer_data ) );

	BBDepNRChiSample< Real > const & nrsample = bbdep_rotameric_data.bbdep_nrchi_sample();

	total_chi.push_back( nrsample.nrchi_mean_ );
	total_ex_steps.push_back( 0. );
	total_rot.push_back( nrsample.nrchi_bin_ );
	chisample_prob.push_back( nrsample.prob_ );

	ExtraRotSample ex_samp_level = rtask.extrachi_sample_level( buried, chi_index, &rsd_type );

	if ( ex_samp_level == NO_EXTRA_CHI_SAMPLES ) return;

	Real const min_extrachi_sd = 0.0; // What's an appropriate value here?
	if ( nrsample.nrchi_sd_ <= min_extrachi_sd ) return;

	for ( Size k=1; k<= extra_steps.size(); ++k ) {
		total_chi.push_back( nrsample.nrchi_mean_  + extra_steps[k] * nrsample.nrchi_sd_  );
		total_ex_steps.push_back( extra_steps[k] );
		total_rot.push_back( nrsample.nrchi_bin_ );
		chisample_prob.push_back( nrsample.prob_  ); // no "penalty" for off-mean samples.
	}

}

//XRW_B_T1
/*
template < Size T >
SingleResidueRotamerLibraryOP
SemiRotamericSingleResidueDunbrackLibrary< T >::coarsify(coarse::Translator const & map ) const
{
	return 0;
}
*/
//XRW_E_T1

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	utility_exit_with_message("Unimplemented!");
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::write_to_binary( utility::io::ozstream & out ) const
{
	using namespace boost;
	parent::write_to_binary( out );
	// For safety, write out const data so that when reading in const data later, we can verify we
	// have a valid input.

	// 1. bbind_nrchi_scoring_
	{
	boost::int32_t bbind_nrchi_scoring = bbind_nrchi_scoring_;
	out.write( (char*) &bbind_nrchi_scoring, sizeof( boost::int32_t ));
	}
	// 2. bbind_nrchi_sampling_
	{
	boost::int32_t bbind_nrchi_sampling = bbind_nrchi_sampling_;
	out.write( (char*) &bbind_nrchi_sampling, sizeof( boost::int32_t ));
	}

	if ( ! bbind_nrchi_scoring_ ) {
		// 3a. bbdep_non_rotameric_chi_scores_
		Size const n_bbdep_scores = grandparent::n_packed_rots() * parent::N_PHIPSI_BINS * parent::N_PHIPSI_BINS * bbdep_nrchi_nbins_;
		BBDepScoreInterpData * bbdep_nrc_interpdata = new BBDepScoreInterpData[ n_bbdep_scores ];
		Size count( 0 );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= parent::N_PHIPSI_BINS; ++ll ) {
						bbdep_nrc_interpdata[ count ] = bbdep_nrc_interpdata_[ ii ]( ll, kk, jj );
						++count;
					}
				}
			}
		}
		out.write( (char*) bbdep_nrc_interpdata, n_bbdep_scores * sizeof( BBDepScoreInterpData ) );
		delete [] bbdep_nrc_interpdata;
	} else {
		// 3b. bbind_non_rotameric_chi_scores_
		Size const n_bbind_scores = bbind_non_rotameric_chi_scores_.size();
		Real * bbind_non_rotameric_chi_scores = new Real[ n_bbind_scores ];
		Size count( 0 );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= bbind_nrchi_nbins_; ++jj ) {
				bbind_non_rotameric_chi_scores[ count ] = bbind_non_rotameric_chi_scores_( jj, ii );
				++count;
			}
		}
		out.write( (char*) bbind_non_rotameric_chi_scores, n_bbind_scores * sizeof( Real ));
		delete [] bbind_non_rotameric_chi_scores;
	}

	/// 4.n_nrchi_sample_bins_
	{
	boost::int32_t n_nrchi_sample_bins =  n_nrchi_sample_bins_;
	out.write( (char*) & n_nrchi_sample_bins, sizeof( boost::int32_t ) );
	}

	if ( ! bbind_nrchi_sampling_ ) {
		// 5. bbdep_rotamers_to_sample_, bbdep_rotsample_sorted_order_
		Size const n_bbdep_nrchi_samples =
			grandparent::n_packed_rots() * parent::N_PHIPSI_BINS * parent::N_PHIPSI_BINS * n_nrchi_sample_bins_;

		boost::int32_t * bbdep_nrchi_packed_rotnos = new boost::int32_t[ n_bbdep_nrchi_samples ];
		boost::int32_t * bbdep_nrchi_bin           = new boost::int32_t[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_means    = new DunbrackReal[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_sdevs    = new DunbrackReal[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_probs    = new DunbrackReal[ n_bbdep_nrchi_samples ];

		boost::int32_t * bbdep_rotsample_sorted_order = new boost::int32_t[ n_bbdep_nrchi_samples ];

		Size count( 0 );
		Size count_iijj( 1 ); // 3d array sorted by frequency instead of 4d array divided first by rotameric rotno.
		for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
			for ( Size jj = 1; jj <= grandparent::n_packed_rots(); ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= parent::N_PHIPSI_BINS; ++ll ) {
						bbdep_nrchi_packed_rotnos[ count ] = bbdep_rotamers_to_sample_( ll, kk, count_iijj ).packed_rotno_;
						bbdep_nrchi_bin[           count ] = bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_bin_;
						bbdep_nrchi_means[         count ] = bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_mean_;
						bbdep_nrchi_sdevs[         count ] = bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_sd_;
						bbdep_nrchi_probs[         count ] = bbdep_rotamers_to_sample_( ll, kk, count_iijj ).prob_;
						bbdep_rotsample_sorted_order[count]= bbdep_rotsample_sorted_order_( ll, kk, jj, ii );
						++count;
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
		delete [] bbdep_nrchi_bin;           bbdep_nrchi_bin = 0;
		delete [] bbdep_nrchi_means;         bbdep_nrchi_means = 0;
		delete [] bbdep_nrchi_sdevs;         bbdep_nrchi_sdevs = 0;
		delete [] bbdep_nrchi_probs;         bbdep_nrchi_probs = 0;
		delete [] bbdep_rotsample_sorted_order; bbdep_rotsample_sorted_order = 0;

	}

	{ // Save the bbind rotamer definitions reguardless of bbind_nrchi_sampling_'s status.

	// 6. bbind_rotamers_to_sample_, bbind_rotamers_sorted_by_probability_
	Size const n_bbind_nrchi_samples = n_nrchi_sample_bins_ * grandparent::n_packed_rots();

	DunbrackReal * bbind_nrchi_lefts   = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_medians = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_rights  = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_probs  = new DunbrackReal[  n_bbind_nrchi_samples ];

	boost::int32_t * bbind_rotsample_sorted_order( 0 );
	if ( bbind_nrchi_sampling_ ) {
		bbind_rotsample_sorted_order = new boost::int32_t[  n_bbind_nrchi_samples ];
	}

	Size count( 0 );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
			bbind_nrchi_lefts[   count ] = bbind_rotamers_to_sample_( jj, ii ).left_;
			bbind_nrchi_medians[ count ] = bbind_rotamers_to_sample_( jj, ii ).median_;
			bbind_nrchi_rights[  count ] = bbind_rotamers_to_sample_( jj, ii ).right_;
			bbind_nrchi_probs[   count ] = bbind_rotamers_to_sample_( jj, ii ).prob_;

			if ( bbind_nrchi_sampling_ ) {
				bbind_rotsample_sorted_order[ count ] = bbind_rotamers_sorted_by_probability_( jj, ii );
			}
			++count;
		}
	}
	out.write( (char*) bbind_nrchi_lefts,   n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbind_nrchi_medians, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbind_nrchi_rights,  n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	out.write( (char*) bbind_nrchi_probs,   n_bbind_nrchi_samples * sizeof( DunbrackReal ) );

	if ( bbind_nrchi_sampling_ ) {
		out.write( (char*) bbind_rotsample_sorted_order, n_bbind_nrchi_samples * sizeof( boost::int32_t ) );
	}

	delete [] bbind_nrchi_lefts;
	delete [] bbind_nrchi_medians;
	delete [] bbind_nrchi_rights;
	delete [] bbind_nrchi_probs;
	delete [] bbind_rotsample_sorted_order;
	}

}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_from_binary( utility::io::izstream & in )
{
	using namespace boost;
	parent::read_from_binary( in );

	// Check that constant member data initialized in the ctor is consistent with this input file
	// 1. bbind_nrchi_scoring_
	{
	boost::int32_t bbind_nrchi_scoring( 0 );
	in.read( (char*) &bbind_nrchi_scoring, sizeof( boost::int32_t ));
	if ( bbind_nrchi_scoring_ != static_cast< bool > ( bbind_nrchi_scoring ) ) {
		if ( bbind_nrchi_scoring ) {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone independent non-rotameric chi scoring,\nbut the instantiated Semirotamer library"
				<< " was created with backbone dependent non-rotameric chi scoring.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
		} else {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone dependent non-rotameric chi scoring,\nbut the instantiated Semirotamer library"
				<< " was created with backbone independent non-rotameric chi scoring.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
		}
		utility_exit();
	}
	}

	// 2. bbind_nrchi_sampling_
	{
	boost::int32_t bbind_nrchi_sampling( 0 );
	in.read( (char*) &bbind_nrchi_sampling, sizeof( boost::int32_t ));
	if ( bbind_nrchi_sampling_ != static_cast< bool > (  bbind_nrchi_sampling ) ) {
		if ( bbind_nrchi_sampling ) {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone independent non-rotameric chi sampling,\nbut the instantiated Semirotamer library"
				<< " was created with backbone dependent non-rotameric chi sampling.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
		} else {
			std::cerr << "ERROR: binary file description for " << grandparent::aa() << " contains data for"
				<< " backbone dependent non-rotameric chi sampling,\nbut the instantiated Semirotamer library"
				<< " was created with backbone independent non-rotameric chi sampling.\nCheck"
				<< " semirotameric parameter definitions in RotamerLibrary::initialize_dun10_aa_parameters.";
		}
		utility_exit();
	}
	}

	if ( ! bbind_nrchi_scoring_ ) {
		// 3a. bbdep_non_rotameric_chi_scores_
		Size const n_bbdep_scores = grandparent::n_packed_rots() * parent::N_PHIPSI_BINS * parent::N_PHIPSI_BINS * bbdep_nrchi_nbins_;
		BBDepScoreInterpData * bbdep_nrc_interpdata = new BBDepScoreInterpData[ n_bbdep_scores ];
		in.read( (char*) bbdep_nrc_interpdata, n_bbdep_scores * sizeof( BBDepScoreInterpData ) );

		Size count( 0 );
		bbdep_nrc_interpdata_.resize( grandparent::n_packed_rots() );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			bbdep_nrc_interpdata_[ ii ].dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, bbdep_nrchi_nbins_ );
			for ( Size jj = 1; jj <= bbdep_nrchi_nbins_; ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= parent::N_PHIPSI_BINS; ++ll ) {
						bbdep_nrc_interpdata_[ ii ]( ll, kk, jj ) = bbdep_nrc_interpdata[ count ];
						++count;
					}
				}
			}
		}
		delete [] bbdep_nrc_interpdata;
	} else {
		// 3b. bbind_non_rotameric_chi_scores_
		Size const n_bbind_scores = grandparent::n_packed_rots() * bbind_nrchi_nbins_;
		Real * bbind_non_rotameric_chi_scores = new Real[ n_bbind_scores ];
		in.read( (char*) bbind_non_rotameric_chi_scores, n_bbind_scores * sizeof( Real ));
		bbind_non_rotameric_chi_scores_.dimension( bbind_nrchi_nbins_, grandparent::n_packed_rots() );

		Size count( 0 );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= bbind_nrchi_nbins_; ++jj ) {
				bbind_non_rotameric_chi_scores_( jj, ii ) = bbind_non_rotameric_chi_scores[ count ];
				++count;
			}
		}
		delete [] bbind_non_rotameric_chi_scores;
	}

	/// 4.n_nrchi_sample_bins_
	{
	boost::int32_t n_nrchi_sample_bins( 0 );
	in.read( (char*) & n_nrchi_sample_bins, sizeof( boost::int32_t ) );
	n_nrchi_sample_bins_ = n_nrchi_sample_bins;
	}


	if ( ! bbind_nrchi_sampling_ ) {
		// 5. bbdep_rotamers_to_sample_, bbdep_rotsample_sorted_order_
		Size const n_bbdep_nrchi_samples =
			grandparent::n_packed_rots() * parent::N_PHIPSI_BINS * parent::N_PHIPSI_BINS * n_nrchi_sample_bins_;

		boost::int32_t * bbdep_nrchi_packed_rotnos = new boost::int32_t[ n_bbdep_nrchi_samples ];
		boost::int32_t * bbdep_nrchi_bin           = new boost::int32_t[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_means    = new DunbrackReal[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_sdevs    = new DunbrackReal[ n_bbdep_nrchi_samples ];
		DunbrackReal * bbdep_nrchi_probs    = new DunbrackReal[ n_bbdep_nrchi_samples ];

		boost::int32_t * bbdep_rotsample_sorted_order = new boost::int32_t[ n_bbdep_nrchi_samples ];

		in.read( (char*) bbdep_nrchi_packed_rotnos, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
		in.read( (char*) bbdep_nrchi_bin, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );
		in.read( (char*) bbdep_nrchi_means, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbdep_nrchi_sdevs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbdep_nrchi_probs, n_bbdep_nrchi_samples * sizeof( DunbrackReal ) );
		in.read( (char*) bbdep_rotsample_sorted_order, n_bbdep_nrchi_samples * sizeof( boost::int32_t ) );

		bbdep_rotamers_to_sample_.dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, grandparent::n_packed_rots()*n_nrchi_sample_bins_ );
		bbdep_rotsample_sorted_order_.dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, grandparent::n_packed_rots(), n_nrchi_sample_bins_ );

		Size count( 0 );
		Size count_iijj( 1 ); // 3d array sorted by frequency instead of 4d array divided first by rotameric rotno.
		for ( Size ii = 1; ii <= n_nrchi_sample_bins_; ++ii ) {
			for ( Size jj = 1; jj <= grandparent::n_packed_rots(); ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= parent::N_PHIPSI_BINS; ++ll ) {
						bbdep_rotamers_to_sample_( ll, kk, count_iijj ).packed_rotno_ = bbdep_nrchi_packed_rotnos[ count ];
						bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_bin_ = bbdep_nrchi_bin[ count ];
						bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_mean_ = bbdep_nrchi_means[ count ];
						bbdep_rotamers_to_sample_( ll, kk, count_iijj ).nrchi_sd_ = bbdep_nrchi_sdevs[ count ];
						bbdep_rotamers_to_sample_( ll, kk, count_iijj ).prob_ = bbdep_nrchi_probs[ count ];
						bbdep_rotsample_sorted_order_( ll, kk, jj, ii ) = bbdep_rotsample_sorted_order[count];
						++count;
					}
				}
			++count_iijj;
			}
		}

		delete [] bbdep_nrchi_packed_rotnos; bbdep_nrchi_packed_rotnos = 0;
		delete [] bbdep_nrchi_bin;           bbdep_nrchi_bin = 0;
		delete [] bbdep_nrchi_means;         bbdep_nrchi_means = 0;
		delete [] bbdep_nrchi_sdevs;         bbdep_nrchi_sdevs = 0;
		delete [] bbdep_nrchi_probs;         bbdep_nrchi_probs = 0;
		delete [] bbdep_rotsample_sorted_order; bbdep_rotsample_sorted_order = 0;

	}

	{ // Save the bbind rotamer definitions reguardless of bbind_nrchi_sampling_'s status.

	// 6. bbind_rotamers_to_sample_, bbind_rotamers_sorted_by_probability_
	Size const n_bbind_nrchi_samples = n_nrchi_sample_bins_ * grandparent::n_packed_rots();

	DunbrackReal * bbind_nrchi_lefts   = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_medians = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_rights  = new DunbrackReal[ n_bbind_nrchi_samples ];
	DunbrackReal * bbind_nrchi_probs  = new DunbrackReal[  n_bbind_nrchi_samples ];
	boost::int32_t * bbind_rotsample_sorted_order = new boost::int32_t[  n_bbind_nrchi_samples ];

	in.read( (char*) bbind_nrchi_lefts, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbind_nrchi_medians, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbind_nrchi_rights, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	in.read( (char*) bbind_nrchi_probs, n_bbind_nrchi_samples * sizeof( DunbrackReal ) );
	if ( bbind_nrchi_sampling_ ) {
		in.read( (char*) bbind_rotsample_sorted_order, n_bbind_nrchi_samples * sizeof( boost::int32_t ) );
	}

	bbind_rotamers_to_sample_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );
	if ( bbind_nrchi_sampling_ ) {
		bbind_rotamers_sorted_by_probability_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );
	}

	Size count( 0 );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
			bbind_rotamers_to_sample_( jj, ii ).left_   = bbind_nrchi_lefts[ count ];
			bbind_rotamers_to_sample_( jj, ii ).median_ = bbind_nrchi_medians[ count ];
			bbind_rotamers_to_sample_( jj, ii ).right_  = bbind_nrchi_rights[  count ];
			bbind_rotamers_to_sample_( jj, ii ).prob_   = bbind_nrchi_probs[   count ];

			if ( bbind_nrchi_sampling_ ) {
				bbind_rotamers_sorted_by_probability_( jj, ii ) = bbind_rotsample_sorted_order[ count ];
			}
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


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::get_rotamer_from_chi(
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
			if (
					! ( left <= nrchi && nrchi <= right ) &&
					! ( left + nrchi_periodicity_ <= nrchi && nrchi <= right + nrchi_periodicity_ ) ) {
				continue;
			}
		} else if ( right > nrchi_lower_angle_ + nrchi_periodicity_ ) {
			if (
					! ( left <= nrchi && nrchi <= right ) &&
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
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::set_nrchi_periodicity(
	Real angle_in_degrees
)
{
	nrchi_periodicity_ = angle_in_degrees;
}

/// @brief What angle do the interpolation data start from?
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::set_nonrotameric_chi_start_angle(
	Real angle_in_degrees
)
{
	nrchi_lower_angle_ = angle_in_degrees;
}

/// @brief What is the angular step size of the bbdep score?
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::set_nonrotameric_chi_bbdep_scoring_step_size(
	Real step_size_in_degrees
)
{
	assert( nrchi_periodicity_ != 0.0 );
	assert( step_size_in_degrees != 0.0 );
	bbdep_nrchi_binsize_ = step_size_in_degrees;
	bbdep_nrchi_nbins_ = static_cast< Size > (nrchi_periodicity_ / step_size_in_degrees);
}

/// @brief What is the angular step size of the bbind score?
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::set_nonrotameric_chi_bbind_scoring_step_size(
	Real step_size_in_degrees
)
{
	assert( nrchi_periodicity_ != 0.0 );
	assert( step_size_in_degrees != 0.0 );
	bbind_nrchi_binsize_ = step_size_in_degrees;
	bbind_nrchi_nbins_ = static_cast< Size > (nrchi_periodicity_ / step_size_in_degrees);
}


template < Size T >
Size
SemiRotamericSingleResidueDunbrackLibrary< T >::memory_usage_dynamic() const
{
	Size total_memory = parent::memory_usage_dynamic();

	/// for bbdep nrchi scoring
	for ( Size ii = 1; ii <= bbdep_nrc_interpdata_.size(); ++ii ) {
		total_memory += bbdep_nrc_interpdata_[ ii ].size() * sizeof( BBDepScoreInterpData );
	}
	total_memory += bbdep_nrc_interpdata_.size() * sizeof( ObjexxFCL::FArray3D< BBDepScoreInterpData > );

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


template < Size T >
Size SemiRotamericSingleResidueDunbrackLibrary< T >::memory_usage_static() const
{
	return sizeof( SemiRotamericSingleResidueDunbrackLibrary< T > );
}


template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_from_files(
	utility::io::izstream & in_rotdef,
	utility::io::izstream & in_rotameric,
	utility::io::izstream & in_continmin_bbdep,
	utility::io::izstream & in_continmin_bbind
)
{
	read_rotamer_definitions( in_rotdef );
	read_rotameric_data( in_rotameric );
	read_bbdep_continuous_minimization_data( in_continmin_bbdep );
	read_bbind_continuous_minimization_data( in_continmin_bbind );
}

/// @breif read-from-files that does not require the in_continmin_bbind file
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_from_files(
	utility::io::izstream & in_rotdef,
	utility::io::izstream & in_rotameric,
	utility::io::izstream & in_continmin_bbdep
)
{
	assert( ! bbind_nrchi_scoring_ ); // make sure you're not actually trying to use backbone-independent scoring.

	read_rotamer_definitions( in_rotdef );
	read_rotameric_data( in_rotameric );
	read_bbdep_continuous_minimization_data( in_continmin_bbdep );
}




/// @details The rotamer definition file has already been read before this
/// function is called, so all of the rotameric rotamer wells have been encountered
/// and the packed_rotno enumeration may be used instead of the regular rotno.
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_rotameric_data(
	utility::io::izstream & in_rotameric
)
{
	std::string three_letter_code;
	Real phi, psi;
	Size count;
	utility::vector1< Size > rotwell( DUNBRACK_MAX_SCTOR, 0 );
	DunbrackReal prob;
	utility::vector1< DunbrackReal > mean( DUNBRACK_MAX_SCTOR, 0.0 );
	utility::vector1< DunbrackReal > stdev( DUNBRACK_MAX_SCTOR, 0.0 );

	ObjexxFCL::FArray2D< Size > rotameric_count( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, Size( 0 ) );
	ObjexxFCL::FArray2D< Size > semi_rotameric_count( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, Size( 0 ));

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

		in_rotameric >> three_letter_code >> phi >> psi >> count;
		in_rotameric >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		in_rotameric >> prob;
		in_rotameric >> mean[ 1 ] >> mean[ 2 ] >> mean[ 3 ] >> mean[ 4 ];
		in_rotameric >> stdev[ 1 ] >> stdev[ 2 ] >> stdev[ 3 ] >> stdev[ 4 ];

		if ( ! in_rotameric ) break; // we've read past the end of the file...

		if ( phi == 180 || psi == 180 ) continue; // duplicated data...

		/// AVOID inf and NaN by correcting the database
		for ( Size ii = 1; ii <= T; ++ii ) {
			if ( stdev[ ii ] == 0.0 ) {
				stdev[ ii ] = 5; // bogus!
			}
		}
		if ( prob == 0.0 ) {
			prob = 1e-4;
			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
		}

		Size phibin, psibin;
		parent::get_phipsi_bins( phi, psi, phibin, psibin );

		++semi_rotameric_count( phibin, psibin );

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
		if ( parent::packed_rotno_2_sorted_rotno()( phibin, psibin, packed_rotno ) == 0 ) {
			// first time this rotwell has been encountered.
			Size nencountered = ++rotameric_count( phibin, psibin );
			parent::packed_rotno_2_sorted_rotno()( phibin, psibin, packed_rotno ) = nencountered;
			for ( Size ii = 1; ii <= T; ++ii ) {
				parent::rotamers()( phibin, psibin, nencountered ).chi_mean( ii, mean[  ii ] );
				parent::rotamers()( phibin, psibin, nencountered ).chi_sd(   ii, stdev[ ii ] );
			}
			parent::rotamers()( phibin, psibin, nencountered ).rotamer_probability() = prob;
			parent::rotamers()( phibin, psibin, nencountered ).rotE() = 0; // temp -- fill in sentinel values for unused variables
			parent::rotamers()( phibin, psibin, nencountered ).rotE_dsecophi() = 0;
			parent::rotamers()( phibin, psibin, nencountered ).rotE_dsecopsi() = 0;
			parent::rotamers()( phibin, psibin, nencountered ).rotE_dsecophipsi() = 0;
		} else {
			Size const sorted_rotno = parent::packed_rotno_2_sorted_rotno()( phibin, psibin, packed_rotno );
			parent::rotamers()( phibin, psibin, sorted_rotno ).rotamer_probability() += prob;
			// verify that the data file is as expected: chimean and chisd do not change for
			// rotameric chi across different entries of the nrchi.
			for ( Size ii = 1; ii <= T; ++ii ) {
				if ( std::abs( basic::periodic_range( parent::rotamers()( phibin, psibin, sorted_rotno ).chi_mean( ii ) - mean[  ii ], 360 )) > 1e-5 ) {
					std::cerr << "ERROR: semi-rotameric input invalid -- chi means differ." << std::endl;
					std::cerr << "Phi: " << phi << " Psi: " << psi << " original chi_" << ii << ": ";
					std::cerr << parent::rotamers()( phibin, psibin, sorted_rotno ).chi_mean( ii );
					std::cerr << " later chi_" << ii << ": "<< mean[  ii ] << std::endl;
					utility_exit();
				}
				if ( std::abs( parent::rotamers()( phibin, psibin, sorted_rotno ).chi_sd( ii ) - stdev[ ii ]) > 1e-5 ) {
					std::cerr << "ERROR: semi-rotameric input invalid -- chi stdevs differ." << std::endl;
					std::cerr << "Phi: " << phi << " Psi: " << psi << " original chi_" << ii << ": ";
					std::cerr << parent::rotamers()( phibin, psibin, sorted_rotno ).chi_sd( ii );
					std::cerr << " later chi_" << ii << ": "<< stdev[  ii ] << std::endl;
					utility_exit();
				}
			}
		}

		if ( ! bbind_nrchi_sampling_ ) {
			Size const which_rotamer = semi_rotameric_count( phibin, psibin );
			BBDepNRChiSample<> & sample( bbdep_rotamers_to_sample_( phibin, psibin, which_rotamer ) );
			sample.packed_rotno_ = packed_rotno;
			sample.nrchi_bin_ = rotwell[ T + 1 ];
			sample.nrchi_mean_ = mean[ T + 1 ];
			sample.nrchi_sd_ = stdev[ T + 1 ];
			sample.prob_ = prob;
			bbdep_rotsample_sorted_order_( phibin, psibin, packed_rotno, rotwell[ T + 1 ]) = which_rotamer;
		}

	}
	parent::initialize_bicubic_splines();


}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_bbdep_continuous_minimization_data(
	utility::io::izstream & in_continmin
)
{

	if ( bbind_nrchi_scoring_ ) return;

	std::string three_letter_code;
	Real phi, psi, base_prob;
	Size count;
	utility::vector1< Size > rotwell( T );
	utility::vector1< Real > chimean( T );
	utility::vector1< Real > chisd( T );
	utility::vector1< Real > nrchi_probs( bbdep_nrchi_nbins_, 0.0 );

	utility::vector1< ObjexxFCL::FArray3D< Real > > bbdep_non_rotameric_chi_scores;
	bbdep_non_rotameric_chi_scores.resize( grandparent::n_packed_rots() );
	for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
		bbdep_non_rotameric_chi_scores[ ii ].dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, bbdep_nrchi_nbins_ );
		bbdep_non_rotameric_chi_scores[ ii ] = 0;
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
		in_continmin >> three_letter_code >> phi >> psi >> count;
		for ( Size ii = 1; ii <= T; ++ii ) {
			in_continmin >> rotwell[ ii ];
		}
		in_continmin >> base_prob;
		for ( Size ii = 1; ii <= T; ++ii ) {
			in_continmin >> chimean[ ii ];
		}
		for ( Size ii = 1; ii <= T; ++ii ) {
			in_continmin >> chisd[ ii ];
		}
		for ( Size ii = 1; ii <= bbdep_nrchi_nbins_; ++ii ) {
			in_continmin >> nrchi_probs[ ii ];
		}

		if ( ! in_continmin ) break; // we've read past the end of the file...
		if ( phi == 180 || psi == 180 ) continue; // duplicated data...

		/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
		/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
		/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
		if (base_prob == 0.0 ) { base_prob = 1e-4; } // correct probability 0 events.

		/// Now, convert the probabilities into scores
		Size phibin, psibin;
		parent::get_phipsi_bins( phi, psi, phibin, psibin );
		Size const rotno = grandparent::rotwell_2_rotno( rotwell );
		for ( Size ii = 1; ii <= bbdep_nrchi_nbins_; ++ii ) {
			/// input data has probability 0 events; assign a tiny probability for these cases.

			/// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			/// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			/// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
			Real const prob = base_prob * ( nrchi_probs[ ii ] == 0.0 ? 1e-4 : nrchi_probs[ ii ] );
			bbdep_non_rotameric_chi_scores[ rotno ]( phibin, psibin, ii ) = -std::log( prob );
		}
	}

	/// Now create the tricubic interpolation data and store that data in the
	if ( true ) {
		//std::cout << "Creating tricubic splines for the backbone-dependent non-rotameric chi scores" << std::endl;
		using namespace numeric;
		using namespace numeric::interpolation::spline;
		BorderFlag border[3] = { e_Periodic, e_Periodic, e_Periodic};
		Real const start[3] = { -180.0, -180.0, nrchi_lower_angle_};
		Real const delta[3] = { 10.0, 10.0, bbdep_nrchi_binsize_ };
		bool const lin_cont[3] = { true, true, true}; // not used since we're imposing periodic boundary conditions
		std::pair< double, double> const first_be[3] = // also not used since we're imposing periodic boundary conditions
		{
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10)
		};

		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			MathTensor< Real > data( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, bbdep_nrchi_nbins_ );
			for ( Size jj = 1; jj <= parent::N_PHIPSI_BINS; ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= bbdep_nrchi_nbins_; ++ll ) {
						data( jj-1, kk-1, ll-1 ) = bbdep_non_rotameric_chi_scores[ ii ]( jj, kk, ll );
					}
				}
			}
			TricubicSpline spline;
			spline.train( border, start, delta, data, lin_cont, first_be );
			for ( Size jj = 1; jj <= parent::N_PHIPSI_BINS; ++jj ) {
				for ( Size kk = 1; kk <= parent::N_PHIPSI_BINS; ++kk ) {
					for ( Size ll = 1; ll <= bbdep_nrchi_nbins_; ++ll ) {
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).value_    = (DunbrackReal) data( jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecox_   = (DunbrackReal) spline.get_dsecox()(   jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoy_   = (DunbrackReal) spline.get_dsecoy()(   jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoz_   = (DunbrackReal) spline.get_dsecoz()(   jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoxy_  = (DunbrackReal) spline.get_dsecoxy()(  jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoxz_  = (DunbrackReal) spline.get_dsecoxz()(  jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoyz_  = (DunbrackReal) spline.get_dsecoyz()(  jj-1, kk-1, ll-1 );
						bbdep_nrc_interpdata_[ ii ]( jj, kk, ll ).dsecoxyz_ = (DunbrackReal) spline.get_dsecoxyz()( jj-1, kk-1, ll-1 );
					}
				}
			}

		}
		//std::cout << "Finished creating tricubic splines" << std::endl;
	}
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_bbind_continuous_minimization_data(
	utility::io::izstream & in_continmin
)
{

	if ( ! bbind_nrchi_scoring_ ) return;

	std::string three_letter_code;
	utility::vector1< Size > rotwell( 4, 0 );
	Real nrchi_val, rho, prob, neglnprob;

	while ( in_continmin ) {
		char first_char = in_continmin.peek();
		if ( first_char == '#' ) {
			std::string line; in_continmin.getline( line );
			continue;
		}

		/// Input line format is
		/// a. Three letter code
		/// b.c.d.e r1, r2, r3, r4
		/// f. nrchi grid point
		/// g. rho
		/// h. prob
		/// i. -lnprob
		in_continmin >> three_letter_code >> rotwell[ 1 ] >> rotwell[ 2 ] >> rotwell[ 3 ] >> rotwell[ 4 ];
		in_continmin >> nrchi_val >> rho >> prob >> neglnprob;

		if ( ! in_continmin ) break; // we've read past the end of the file...

		Size const packed_rotno( grandparent::rotwell_2_packed_rotno( rotwell ) );
		Size bbind_nrchi_bin, dummy; Real dummy_alpha;
		get_bbind_nrchi_bin( nrchi_val, bbind_nrchi_bin, dummy, dummy_alpha );
		bbind_non_rotameric_chi_scores_( bbind_nrchi_bin, packed_rotno ) = neglnprob;
	}
}

/// @details the rotamer definition file must be the first file read.
/// Reading this file gives the number of pseudorotamers used for
/// both the bbdep and bbind rotamer building.  The bbind rotamer data
/// is saved iff the bbind_rotamer_building_ flag is true.
template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::read_rotamer_definitions(
	utility::io::izstream & in_rotdef
)
{
	utility::vector1< Size > rotwell( 4, 0 );
	Real prob, lnprob, left, median, right;

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
		if ( rotwell[ T + 1 ] == 1 ) {
			grandparent::mark_rotwell_exists( rotwell );
		}
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

			/// DONE -- quit this loop
			break;
		}
	}


	/// After reading this file, we have seen all the rotameric chi that we are going to.
	/// Tell the base class that the rotno to packedrotno conversion is now appropriate.
	/// This triggers a call to convert rotno indices to packedrotno indices
	grandparent::declare_all_existing_rotwells_encountered();

	/// Allocate space for rotamers and rotamer-sorted-order mapping
	parent::rotamers().dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, grandparent::n_packed_rots() );
	parent::packed_rotno_2_sorted_rotno().dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, grandparent::n_packed_rots() );
	parent::packed_rotno_2_sorted_rotno() = 0;

	/// save this data for rotamer creation and for binning nrchi values regardless of whether
	/// we're using bbind rotamer sampling.
	bbind_rotamers_to_sample_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );
	// demension this mapping array iff bbind_nrchi_sampling_.
	//bbind_rotamers_sorted_by_probability_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );
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
	if ( ! bbind_nrchi_sampling_ ) {
		/// Allocate space for bbdep rotamer sampling for the nrchi.
		bbdep_rotamers_to_sample_.dimension(
			parent::N_PHIPSI_BINS,
			parent::N_PHIPSI_BINS,
			grandparent::n_packed_rots() * n_nrchi_sample_bins_ );
		bbdep_rotsample_sorted_order_.dimension(
			parent::N_PHIPSI_BINS,
			parent::N_PHIPSI_BINS,
			grandparent::n_packed_rots(),
			n_nrchi_sample_bins_ );
		bbdep_rotsample_sorted_order_ = 0;
	}

	if ( bbind_nrchi_scoring_ ) {
		assert( bbind_nrchi_nbins_ != 0 );
		bbind_non_rotameric_chi_scores_.dimension( bbind_nrchi_nbins_, grandparent::n_packed_rots() );
		bbind_non_rotameric_chi_scores_ = 0;
	} else {
		assert( bbdep_nrchi_nbins_ != 0 );
		bbdep_nrc_interpdata_.resize( grandparent::n_packed_rots() );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			bbdep_nrc_interpdata_[ ii ].dimension( parent::N_PHIPSI_BINS, parent::N_PHIPSI_BINS, bbdep_nrchi_nbins_ );
		}
	}

	if ( bbind_nrchi_sampling_ ) {
		bbind_rotamers_sorted_by_probability_.dimension( n_nrchi_sample_bins_, grandparent::n_packed_rots() );
		utility::vector1< ProbSortClass > rot_probs_and_inds( n_nrchi_sample_bins_ );
		for ( Size ii = 1; ii <= grandparent::n_packed_rots(); ++ii ) {
			for ( Size jj = 1; jj <= n_nrchi_sample_bins_; ++jj ) {
				rot_probs_and_inds[ jj ].probability_ = bbind_rotamers_to_sample_( jj, ii ).prob_;
				rot_probs_and_inds[ jj ].index_ = jj;
			}
			/// sort into ascending order by probability
			std::sort( rot_probs_and_inds.begin(), rot_probs_and_inds.end(), psc_compare );
			for ( Size jj = 1, jj_down = n_nrchi_sample_bins_; jj_down >= 1; --jj_down, ++jj ) {
				/// store in descending order; most likely to least likely
				bbind_rotamers_sorted_by_probability_( jj, ii ) = rot_probs_and_inds[ jj_down ].index_;
			}
		}

	}
}

/// @details Clips to range [ nrchi_lower_angle_, nrchi_lower_angle_ + nrchi_periodicity )
/// Note: this isn't really "clipping" (i.e. values outside of this range are not set to
/// the minimum or maximum of the range).  Instead, it is properly returning the value
/// within the range corresponding to a value outside of the range.
template < Size T >
Real
SemiRotamericSingleResidueDunbrackLibrary< T >::clip_to_nrchi_range( Real chi ) const
{
	Real chip = basic::periodic_range( chi, nrchi_periodicity_ );

	if ( chip >= nrchi_lower_angle_ + nrchi_periodicity_ ) {
		while ( chip >= nrchi_lower_angle_ + nrchi_periodicity_ ) {
			chip -= nrchi_periodicity_;
		}
	} else if ( chip < nrchi_lower_angle_ ) {
		while ( chip < nrchi_lower_angle_ ) {
			chip += nrchi_periodicity_;
		}
	}
	return chip;
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::get_bbdep_nrchi_bin(
	Real nrchi,
	Size & bin_lower,
	Size & bin_upper,
	Real & nrchi_alpha
) const
{
	Real clipped_nrchi( clip_to_nrchi_range( nrchi ));
	grandparent::bin_angle(
		nrchi_lower_angle_, bbdep_nrchi_binsize_, nrchi_periodicity_, bbdep_nrchi_nbins_,
		clipped_nrchi, bin_lower, bin_upper, nrchi_alpha
	);
}

template < Size T >
void
SemiRotamericSingleResidueDunbrackLibrary< T >::get_bbind_nrchi_bin(
	Real nrchi,
	Size & bin_lower,
	Size & bin_upper,
	Real & nrchi_alpha
) const
{
	Real clipped_nrchi( clip_to_nrchi_range( nrchi ));
	grandparent::bin_angle(
		nrchi_lower_angle_, bbind_nrchi_binsize_, nrchi_periodicity_, bbind_nrchi_nbins_,
		clipped_nrchi, bin_lower, bin_upper, nrchi_alpha
	);
}


template < Size T >
BBDepNRChiSample< Real >
SemiRotamericSingleResidueDunbrackLibrary< T >::interpolate_bbdep_nrchi_sample(
	Size const packed_rotno,
	Size const nrchi_bin,
	Size const phibin,
	Size const psibin,
	Size const phibin_next,
	Size const psibin_next,
	Real const phi_alpha,
	Real const psi_alpha
) const
{
	Size const
	ind00( bbdep_rotsample_sorted_order_( phibin     , psibin     , packed_rotno, nrchi_bin )),
	ind01( bbdep_rotsample_sorted_order_( phibin     , psibin_next, packed_rotno, nrchi_bin )),
	ind10( bbdep_rotsample_sorted_order_( phibin_next, psibin     , packed_rotno, nrchi_bin )),
	ind11( bbdep_rotsample_sorted_order_( phibin_next, psibin_next, packed_rotno, nrchi_bin ));

	BBDepNRChiSample<> const & nrchi_sample_00( bbdep_rotamers_to_sample_( phibin     , psibin     , ind00 ));
	BBDepNRChiSample<> const & nrchi_sample_01( bbdep_rotamers_to_sample_( phibin     , psibin_next, ind01 ));
	BBDepNRChiSample<> const & nrchi_sample_10( bbdep_rotamers_to_sample_( phibin_next, psibin     , ind10 ));
	BBDepNRChiSample<> const & nrchi_sample_11( bbdep_rotamers_to_sample_( phibin_next, psibin_next, ind11 ));

	return interpolate_bbdep_nrchi_sample(
		nrchi_sample_00, nrchi_sample_01,
		nrchi_sample_10, nrchi_sample_11,
		phi_alpha, psi_alpha
	);
}

template < Size T >
BBDepNRChiSample< Real >
SemiRotamericSingleResidueDunbrackLibrary< T >::interpolate_bbdep_nrchi_sample(
	BBDepNRChiSample<> const & nrchi_sample_00,
	BBDepNRChiSample<> const & nrchi_sample_01,
	BBDepNRChiSample<> const & nrchi_sample_10,
	BBDepNRChiSample<> const & nrchi_sample_11,
	Real const phi_alpha,
	Real const psi_alpha
) const
{
	BBDepNRChiSample< Real > interpolated_sample;
	interpolated_sample.packed_rotno_ = nrchi_sample_00.packed_rotno_;
	interpolated_sample.nrchi_bin_    = nrchi_sample_00.nrchi_bin_;

	Real dummy_dprob_1, dummy_dprob_2;
	basic::interpolate_bilinear_by_value(
		static_cast< Real >  ( nrchi_sample_00.prob_),
		static_cast< Real >  ( nrchi_sample_10.prob_),
		static_cast< Real >  ( nrchi_sample_01.prob_),
		static_cast< Real >  ( nrchi_sample_11.prob_),
		phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, false /*treat_as_angles*/,
		interpolated_sample.prob_,
		dummy_dprob_1, dummy_dprob_2
	);
	Real dummy_dmean_1, dummy_dmean_2;
	basic::interpolate_bilinear_by_value(
		static_cast< Real >  ( nrchi_sample_00.nrchi_mean_),
		static_cast< Real >  ( nrchi_sample_10.nrchi_mean_),
		static_cast< Real >  ( nrchi_sample_01.nrchi_mean_),
		static_cast< Real >  ( nrchi_sample_11.nrchi_mean_),
		phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, true /*treat_as_angles*/,
		interpolated_sample.nrchi_mean_,
		dummy_dmean_1, dummy_dmean_2
	);
	Real dummy_dsd_1, dummy_dsd_2;
	basic::interpolate_bilinear_by_value(
		static_cast< Real >  ( nrchi_sample_00.nrchi_sd_),
		static_cast< Real >  ( nrchi_sample_10.nrchi_sd_),
		static_cast< Real >  ( nrchi_sample_01.nrchi_sd_),
		static_cast< Real >  ( nrchi_sample_11.nrchi_sd_),
		phi_alpha, psi_alpha, parent::PHIPSI_BINRANGE, false /*treat_as_angles*/,
		interpolated_sample.nrchi_sd_,
		dummy_dsd_1, dummy_dsd_2
	);
	return interpolated_sample;
}

/// @details never call this ctor!
//template < Size T >
//SemiRotamericSingleResidueDunbrackLibrary< T >::SemiRotamericSingleResidueDunbrackLibrary()
//	:
//{
//	utility_exit();
//}



} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_TMPL_HH
