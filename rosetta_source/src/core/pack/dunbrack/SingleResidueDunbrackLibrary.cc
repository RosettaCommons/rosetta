// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

// Package headers
// AUTO-REMOVED #include <core/pack/dunbrack/SingleResidueDunbrackLibrary.tmpl.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/ChiSet.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/CoarseRotamer.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/CoarseSingleResidueLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>

// Project headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

// AUTO-REMOVED #include <core/coarse/Translator.hh>
#include <core/graph/Graph.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pack/task/PackerTask_.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <basic/basic.hh>
// AUTO-REMOVED #include <basic/interpolate.hh>
#include <basic/Tracer.hh>

#include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSetOperation.hh>

// Utility headers
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/LexicographicalIterator.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

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
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
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
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/PackerTask_.fwd.hh>
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
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
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
#include <utility/vector1.functions.hh>
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
#include <numeric/random/random.hh>
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
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
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
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/interpolate.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers



// Numeric headers
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

namespace core {
namespace pack {
namespace dunbrack {

static basic::Tracer SRDL_TR("core.pack.dunbrack");
static const Real MIN_ROT_PROB = 1.e-8;

Real const SingleResidueDunbrackLibrary::NEUTRAL_PHI = -90; // R++ value.  Roland Dunbrack suggests -60.
Real const SingleResidueDunbrackLibrary::NEUTRAL_PSI = 130; // R++ value.  Roland Dunbrack suggests  60.

SingleResidueDunbrackLibrary::SingleResidueDunbrackLibrary(
	AA const aa,
	Size const n_rotameric_chi,
	bool dun02
) :
	dun02_( dun02 ),
	aa_( aa ),
	n_rotameric_chi_( n_rotameric_chi ),
	n_chi_bins_( n_rotameric_chi, 0 ),
	n_chi_products_( n_rotameric_chi, 1),
	n_packed_rots_( 0 ),
	n_possible_rots_( 0 ),
	prob_to_accumulate_buried_(0.98),
	prob_to_accumulate_nonburied_(0.95),
	packed_rotno_conversion_data_current_( false )
{
	// this builds on the hard coded hack bellow
	// since NCAAs are aa_unk we cannot hard code the information
	// alternativly it is added to the residue type paramater files
	read_options();
	if (aa_ != chemical::aa_unk) {
		n_rotameric_bins_for_aa( aa_, n_chi_bins_, dun02 );
		for ( Size ii = n_rotameric_chi_; ii > 1; --ii ) {
			n_chi_products_[ ii - 1 ] = n_chi_products_[ ii ] * n_chi_bins_[ ii ];
		}
		n_possible_rots_ = n_rotameric_chi_ == 0 ? 0 : n_chi_products_[ 1 ] * n_chi_bins_[ 1 ];
		rotwell_exists_.resize( n_possible_rots_ );
		rotno_2_packed_rotno_.resize( n_possible_rots_ );
		std::fill( rotwell_exists_.begin(), rotwell_exists_.end(), false );
		std::fill( rotno_2_packed_rotno_.begin(), rotno_2_packed_rotno_.end(), 0 );
	}
}

/// @details Sets the number of bins for a particular chi angle, used for the NCAAs, info for CAAs is hardcoded bellow
void
SingleResidueDunbrackLibrary::set_n_chi_bins( utility::vector1< Size > const & n_chi_bins )
{
	// load rot size vector
	n_chi_bins_.resize( n_chi_bins.size() );
	for ( Size i = 1; i <= n_chi_bins.size(); ++i ) {
		n_chi_bins_[i] = n_chi_bins[i];
	}

	// basically the same as the ctor
	for ( Size ii = n_rotameric_chi_; ii > 1; --ii ) {
		n_chi_products_[ ii - 1 ] = n_chi_products_[ ii ] * n_chi_bins_[ ii ];
	}
	n_possible_rots_ = n_rotameric_chi_ == 0 ? 0 : n_chi_products_[ 1 ] * n_chi_bins_[ 1 ];
	rotwell_exists_.resize( n_possible_rots_ );
	rotno_2_packed_rotno_.resize( n_possible_rots_ );
	std::fill( rotwell_exists_.begin(), rotwell_exists_.end(), false );
	std::fill( rotno_2_packed_rotno_.begin(), rotno_2_packed_rotno_.end(), 0 );
}

SingleResidueDunbrackLibrary::~SingleResidueDunbrackLibrary() {}

void
SingleResidueDunbrackLibrary::read_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::packing::dunbrack_prob_buried ].user() )
		prob_to_accumulate_buried( option[ OptionKeys::packing::dunbrack_prob_buried ]() );
	if ( option[ OptionKeys::packing::dunbrack_prob_nonburied ].user() )
		prob_to_accumulate_nonburied( option[ OptionKeys::packing::dunbrack_prob_nonburied ]() );
}

/// @details the number of wells defined for the 08 library; includes
/// the number of wells for the semi-rotameric chi (arbitrarily chosen).
void
SingleResidueDunbrackLibrary::n_rotamer_bins_for_aa(
	chemical::AA const aa,
	RotVector & rot
)
{
	/// HARD CODED HACK -- SHOULD MOVE THIS TO BECOME A VIRTUAL FUNCTION CALL
	using namespace chemical;
	switch ( aa ) {
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 6; break;
		case aa_glu: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 6; break;
		case aa_phe: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 12; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_tyr: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		default:
			rot.resize( 0 );
	}
}

/// @details For the rotameric chi -- not all chi are rotameric
void
SingleResidueDunbrackLibrary::n_rotameric_bins_for_aa(
	chemical::AA const aa,
	RotVector & rot,
	bool dun02
)
{
	if ( dun02 ) {
		n_rotamer_bins_for_aa_02( aa, rot );
		return;
	}


	using namespace chemical;
	switch ( aa ) {
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_glu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_phe: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_tyr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		default:
			rot.resize( 0 );
	}
}

/// @details To continue supporting the 2002 Dunbrack library, we
/// need to preserve the rotamer well definitions that it made.
void
SingleResidueDunbrackLibrary::n_rotamer_bins_for_aa_02(
	chemical::AA const aa,
	RotVector & rot
)
{
	using namespace chemical;
	switch ( aa ) {
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_glu: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_phe: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 4; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_tyr: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
		default:
			rot.resize( 0 );
	}
}


/// 2002 Library hard code symmetry information.
Real
subtract_chi_angles(
	Real chi1,
	Real chi2,
	chemical::AA const & aa,
	int chino
)
{
	using namespace chemical;
	using basic::periodic_range;


	if ( chino == 2 ) {   // handle symmetry cases for chi2
		//  -- not for HIS in new dunbrack version
		if ( aa == aa_phe || aa == aa_tyr ) {
			return periodic_range((chi1-chi2),180.);
		} else if ( aa == aa_asp ) {
			return periodic_range((chi1-chi2),180.);
		}
	} else if ( chino == 3 ) {   // handle symmetry cases for chi3
		if ( aa == aa_glu ) {
			return periodic_range((chi1-chi2),180.);
		}
	}
	return periodic_range((chi1-chi2),360.);
}


/// @brief The base class needs to be informed about which rotamer wells
/// exist in order to create the rotwell to packed rot conversion data.
/// set_chi_nbins must be called first.
void
SingleResidueDunbrackLibrary::mark_rotwell_exists(
	utility::vector1< Size > const & rotwell
)
{
	assert( ! packed_rotno_conversion_data_current_ );
	Size const rotno( rotwell_2_rotno( rotwell ) );
	rotwell_exists_[ rotno ] = true;
}

/// @brief After the derived class has marked all the rotwells that do exist,
/// the base class will create the rotwell to packerot conversion data.
void
SingleResidueDunbrackLibrary::declare_all_existing_rotwells_encountered()
{
	Size count = 0;
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if ( rotwell_exists_[ ii ] ) ++count;
	}
	n_packed_rots_ = count;
	packed_rotno_2_rotno_.resize( n_packed_rots_ );
	packed_rotno_2_rotwell_.resize( n_packed_rots_ );
	count = 0;
	utility::vector1< Size > rotwell;
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if ( rotwell_exists_[ ii ] ) {
			++count;
			packed_rotno_2_rotno_[ count ] = ii;
			rotno_2_packed_rotno_[ ii    ] = count;
			rotno_2_rotwell( ii, rotwell );
			packed_rotno_2_rotwell_[ count ] = rotwell;
		} else {
			rotno_2_packed_rotno_[ ii ] = 0;
		}
	}
	packed_rotno_conversion_data_current_ = true;

	// force deallocatation of rotamer_exists_ data
	// neither resize() nor clear() guarantee emptying; swap does.
	utility::vector1< bool > empty_vector;
	rotwell_exists_.swap( empty_vector );

}

/// Conversion functions

/// @brief Convert from the rotamer bin indices for each chi to the
/// (non-compact) "rotamer number"
Size
SingleResidueDunbrackLibrary::rotwell_2_rotno(
	utility::vector1< Size > const & rotwell
) const
{
	assert( n_chi_products_.size() <= rotwell.size() );
	Size runsum = 1;
	for ( Size ii = 1; ii <= n_chi_products_.size(); ++ii ) {
		runsum += n_chi_products_[ ii ]*( rotwell[ ii ] - 1 );
	}
	return runsum;
}

/// @brief Convert from the rotamer bin indices for each chi to the
/// (non-compact) "rotamer number"
Size
SingleResidueDunbrackLibrary::rotwell_2_rotno( Size4 const & rotwell ) const
{
	assert( n_chi_products_.size() <= rotwell.size() );
	Size runsum = 1;
	for ( Size ii = 1; ii <= n_chi_products_.size(); ++ii ) {
		runsum += n_chi_products_[ ii ]*( rotwell[ ii ] - 1 );
	}
	return runsum;

}


/// @brief Convert from the rotamer number to the compacted
/// "packed rotamer number".  Returns 0 if rotno has no corresponding packed rotno.
Size
SingleResidueDunbrackLibrary::rotno_2_packed_rotno( Size const rotno ) const
{
	assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno_[ rotno ];
}


Size
SingleResidueDunbrackLibrary::rotwell_2_packed_rotno(
	utility::vector1< Size > const & rotwell
) const
{
	assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}

/// @brief Convert from the rotamer bin indices for each chi to the
/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
Size
SingleResidueDunbrackLibrary::rotwell_2_packed_rotno( Size4 const & rotwell ) const
{
	assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}


void
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	utility::vector1< Size > & rotwell
) const
{
	rotwell = packed_rotno_2_rotwell_[ packed_rotno ];
}

void
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	Size4 & rotwell
) const
{
	assert( packed_rotno_2_rotwell_[ packed_rotno ].size() <= rotwell.size() );
	std::copy( packed_rotno_2_rotwell_[ packed_rotno ].begin(), packed_rotno_2_rotwell_[ packed_rotno ].end(), rotwell.begin() );
}

utility::vector1< Size > const &
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell( Size const packed_rotno ) const
{
	return packed_rotno_2_rotwell_[ packed_rotno ];
}


void
SingleResidueDunbrackLibrary::write_to_binary( utility::io::ozstream & out ) const
{
	using namespace boost;
	/// 1. n_packed_rots_
	{
	boost::int32_t n_packed_rots( n_packed_rots_ );
	out.write( (char*) & n_packed_rots, sizeof( boost::int32_t ));
	}


	/// 2. rotno_2_packed_rotno_
	{
	boost::int32_t * rotno_2_packed_rotno = new boost::int32_t[ n_possible_rots_ ];
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) rotno_2_packed_rotno[ ii - 1 ] = rotno_2_packed_rotno_[ ii ];
	out.write( (char*) rotno_2_packed_rotno, n_possible_rots_ * sizeof( boost::int32_t ) );
	delete [] rotno_2_packed_rotno; rotno_2_packed_rotno = 0;
	}

	/// 3. packed_rotno_2_rotno_
	{
	boost::int32_t * packed_rotno_2_rotno = new boost::int32_t[ n_packed_rots_ ];
	for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) packed_rotno_2_rotno[ ii - 1 ] = packed_rotno_2_rotno_[ ii ];
	out.write( (char*) packed_rotno_2_rotno, n_packed_rots_ * sizeof( boost::int32_t ) );
	delete [] packed_rotno_2_rotno; packed_rotno_2_rotno = 0;
	}

	/// 4. packed_rotno_2_rotwell_
	{
	boost::int32_t * packed_rotno_2_rotwell = new boost::int32_t[ n_packed_rots_ * n_rotameric_chi_ ];
	Size count( 0 );
	for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
		for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
			packed_rotno_2_rotwell[ count ] = packed_rotno_2_rotwell_[ ii ][ jj ];
			++count;
		}
	}
	out.write( (char*) packed_rotno_2_rotwell, n_packed_rots_ * n_rotameric_chi_ * sizeof( boost::int32_t ) );
	delete [] packed_rotno_2_rotwell; packed_rotno_2_rotwell = 0;
	}

}

void
SingleResidueDunbrackLibrary::read_from_binary( utility::io::izstream & in )
{

	/// 1. n_packed_rots_
	{
	boost::int32_t n_packed_rots( 0 );
	in.read( (char*) & n_packed_rots, sizeof( boost::int32_t ));
	n_packed_rots_ = n_packed_rots;
	}

	/// 2. rotno_2_packed_rotno_
	{
	boost::int32_t * rotno_2_packed_rotno = new boost::int32_t[ n_possible_rots_ ];
	in.read( (char*) rotno_2_packed_rotno, n_possible_rots_ * sizeof( boost::int32_t ) );
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) rotno_2_packed_rotno_[ ii ] = rotno_2_packed_rotno[ ii - 1 ];
	delete [] rotno_2_packed_rotno; rotno_2_packed_rotno = 0;
	}

	/// 3. packed_rotno_2_rotno_
	{
	boost::int32_t * packed_rotno_2_rotno = new boost::int32_t[ n_packed_rots_ ];
	in.read( (char*) packed_rotno_2_rotno, n_packed_rots_ * sizeof( boost::int32_t ) );
	packed_rotno_2_rotno_.resize( n_packed_rots_ );
	for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) packed_rotno_2_rotno_[ ii ] = packed_rotno_2_rotno[ ii - 1 ];
	delete [] packed_rotno_2_rotno; packed_rotno_2_rotno = 0;
	}

	/// 4. packed_rotno_2_rotwell_
	{
	boost::int32_t * packed_rotno_2_rotwell = new boost::int32_t[ n_packed_rots_ * n_rotameric_chi_ ];
	in.read( (char*) packed_rotno_2_rotwell, n_packed_rots_ * n_rotameric_chi_ * sizeof( boost::int32_t ) );
	packed_rotno_2_rotwell_.resize( n_packed_rots_ );
	Size count( 0 );
	for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
		packed_rotno_2_rotwell_[ ii ].resize( n_rotameric_chi_ );
		for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
			packed_rotno_2_rotwell_[ ii ][ jj ] = packed_rotno_2_rotwell[ count ];
			++count;
		}
	}
	delete [] packed_rotno_2_rotwell;
	}
	packed_rotno_conversion_data_current_ = true;
}

/// @details not as fast as going through the packed_rotno_2_rotwell_
/// lookup table, but does the modulo converstion from a 1-based index
/// to a lexicographical index ordering.
///
/// if there are 3 chi, and 3 rotamer bins per chi, then
/// 21 would represent (3-1) * 3**2 + (1-1) * 3**1 + (3-1) * 3**0  + 1 = [ 3, 1, 3 ];
void
SingleResidueDunbrackLibrary::rotno_2_rotwell(
	Size const rotno,
	utility::vector1< Size > & rotwell
) const
{
	Size remainder = rotno - 1;
	rotwell.resize( n_rotameric_chi_ );
	std::fill( rotwell.begin(), rotwell.end(), 1 );
	for ( Size ii = 1; ii <= n_rotameric_chi_; /* no increment */ ) {
		if ( remainder > n_chi_products_[ ii ] ) {
			remainder -= n_chi_products_[ ii ];
			++rotwell[ ii ];
		} else {
			++ii;
		}
	}
}

Real
SingleResidueDunbrackLibrary::probability_to_accumulate_while_building_rotamers(
	bool buried
) const
{
	return ( buried ? prob_to_accumulate_buried_ : prob_to_accumulate_nonburied_ );
}

/// @brief setters for accumulation probability cutoff (to support externally-controlled option dependence)
void
SingleResidueDunbrackLibrary::prob_to_accumulate( Real buried, Real nonburied )
{
	prob_to_accumulate_buried( buried );
	prob_to_accumulate_nonburied( nonburied );
}
void
SingleResidueDunbrackLibrary::prob_to_accumulate_buried( Real buried )
{
	if ( buried <= 0. || buried > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_buried_ = buried;
}
void
SingleResidueDunbrackLibrary::prob_to_accumulate_nonburied( Real nonburied )
{
	if ( nonburied <= 0. || nonburied > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_nonburied_ = nonburied;
}

Size SingleResidueDunbrackLibrary::memory_usage_in_bytes() const
{
	return memory_usage_static() + memory_usage_dynamic();
}

Size SingleResidueDunbrackLibrary::memory_usage_dynamic() const
{
	Size total = 0;
	total += n_chi_bins_.size() * sizeof( Size );
	total += n_chi_products_.size() * sizeof( Size );
	total += rotwell_exists_.size() * sizeof( Size );
	total += rotno_2_packed_rotno_.size() * sizeof( Size );
	total += packed_rotno_2_rotno_.size() * sizeof( Size );
	for ( Size ii = 1; ii <= packed_rotno_2_rotwell_.size(); ++ii ) {
		total += packed_rotno_2_rotwell_[ ii ].size() * sizeof( Size );
	}
	total += packed_rotno_2_rotwell_.size() * sizeof( utility::vector1< Size > );
	return total;
}


/// @details forces instantiation of virtual functions for templated
/// derived classes... part of the uglyness of mixing templates and
/// polymorphism.  Never invoke this function.
void
SingleResidueDunbrackLibrary::hokey_template_workaround()
{
	assert( false );
	utility_exit_with_message(
	"ERROR: SingleResidueDunbrackLibrary::hokey_template_workaround should never be called!");

	RotamericSingleResidueDunbrackLibrary< ONE >   rsrdl_1( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO >   rsrdl_2( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE > rsrdl_3( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR >  rsrdl_4( chemical::aa_ala, false );

	SemiRotamericSingleResidueDunbrackLibrary< ONE >   srsrdl_1( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO >   srsrdl_2( chemical::aa_ala, true, true ); // e.g. glu
	// three and four do not exist... they could in the future if needed.

	chemical::ResidueType blah( 0, 0, 0, 0 ); ///TODO shouldnt this variable be named something else?
	conformation::Residue rsd( blah, true );
	RotamerLibraryScratchSpace scratch;

	rsrdl_1.nchi();
	rsrdl_2.nchi();
	rsrdl_3.nchi();
	rsrdl_4.nchi();
	srsrdl_1.nchi();
	srsrdl_2.nchi();

	rsrdl_1.n_rotamer_bins();
	rsrdl_2.n_rotamer_bins();
	rsrdl_3.n_rotamer_bins();
	rsrdl_4.n_rotamer_bins();
	srsrdl_1.n_rotamer_bins();
	srsrdl_2.n_rotamer_bins();


	rsrdl_1.rotamer_energy( rsd, scratch );
	rsrdl_2.rotamer_energy( rsd, scratch );
	rsrdl_3.rotamer_energy( rsd, scratch );
	rsrdl_4.rotamer_energy( rsd, scratch );
	srsrdl_1.rotamer_energy( rsd, scratch );
	srsrdl_2.rotamer_energy( rsd, scratch );

	rsrdl_1.rotamer_energy_deriv( rsd, scratch );
	rsrdl_2.rotamer_energy_deriv( rsd, scratch );
	rsrdl_3.rotamer_energy_deriv( rsd, scratch );
	rsrdl_4.rotamer_energy_deriv( rsd, scratch );
	srsrdl_1.rotamer_energy_deriv( rsd, scratch );
	srsrdl_2.rotamer_energy_deriv( rsd, scratch );

	rsrdl_1.best_rotamer_energy( rsd, true, scratch );
	rsrdl_2.best_rotamer_energy( rsd, true, scratch );
	rsrdl_3.best_rotamer_energy( rsd, true, scratch );
	rsrdl_4.best_rotamer_energy( rsd, true, scratch );
	srsrdl_1.best_rotamer_energy( rsd, true, scratch );
	srsrdl_2.best_rotamer_energy( rsd, true, scratch );

	pose::Pose pose; scoring::ScoreFunction sfxn; pack::task::PackerTask_ task( pose );
	core::graph::GraphOP png;
	chemical::ResidueTypeCAP cr;
	utility::vector1< utility::vector1< Real > > ecs;
	RotamerVector rv;

	rsrdl_1.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_2.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_3.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_4.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_1.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_2.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );

	utility::io::ozstream ozs;
	rsrdl_1.write_to_file( ozs );
	rsrdl_2.write_to_file( ozs );
	rsrdl_3.write_to_file( ozs );
	rsrdl_4.write_to_file( ozs );
	srsrdl_1.write_to_file( ozs );
	srsrdl_2.write_to_file( ozs );

	utility::io::ozstream os;
	rsrdl_1.write_to_binary( os );
	rsrdl_2.write_to_binary( os );
	rsrdl_3.write_to_binary( os );
	rsrdl_4.write_to_binary( os );
	srsrdl_1.write_to_binary( os );
	srsrdl_2.write_to_binary( os );

	utility::io::izstream is;
	rsrdl_1.read_from_binary( is );
	rsrdl_2.read_from_binary( is );
	rsrdl_3.read_from_binary( is );
	rsrdl_4.read_from_binary( is );
	srsrdl_1.read_from_binary( is );
	srsrdl_2.read_from_binary( is );

	rsrdl_1.memory_usage_in_bytes();
	rsrdl_2.memory_usage_in_bytes();
	rsrdl_3.memory_usage_in_bytes();
	rsrdl_4.memory_usage_in_bytes();
	srsrdl_1.memory_usage_in_bytes();
	srsrdl_2.memory_usage_in_bytes();

	utility::vector1< Real > chi; utility::vector1< Size > rot;
	rsrdl_1.get_rotamer_from_chi( chi, rot );
	rsrdl_2.get_rotamer_from_chi( chi, rot );
	rsrdl_3.get_rotamer_from_chi( chi, rot );
	rsrdl_4.get_rotamer_from_chi( chi, rot );
	srsrdl_1.get_rotamer_from_chi( chi, rot );
	srsrdl_2.get_rotamer_from_chi( chi, rot );

	rsrdl_1.memory_usage_dynamic();
	rsrdl_2.memory_usage_dynamic();
	rsrdl_3.memory_usage_dynamic();
	rsrdl_4.memory_usage_dynamic();
	srsrdl_1.memory_usage_dynamic();
	srsrdl_2.memory_usage_dynamic();

	rsrdl_1.memory_usage_static();
	rsrdl_2.memory_usage_static();
	rsrdl_3.memory_usage_static();
	rsrdl_4.memory_usage_static();
	srsrdl_1.memory_usage_static();
	srsrdl_2.memory_usage_static();

	rsrdl_1.get_all_rotamer_samples( 0.0, 0.0 );
	rsrdl_2.get_all_rotamer_samples( 0.0, 0.0 );
	rsrdl_3.get_all_rotamer_samples( 0.0, 0.0 );
	rsrdl_4.get_all_rotamer_samples( 0.0, 0.0 );
	srsrdl_1.get_all_rotamer_samples( 0.0, 0.0 );
	srsrdl_2.get_all_rotamer_samples( 0.0, 0.0 );


	numeric::random::RandomGenerator RG(1);
	ChiVector chiv;

	rsrdl_1.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );
	rsrdl_2.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );
	rsrdl_3.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );
	rsrdl_4.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );
	srsrdl_1.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );
	srsrdl_2.assign_random_rotamer_with_bias( rsd, scratch, RG, chiv, true );

	rsrdl_1.get_probability_for_rotamer( 0.0, 0.0, 1 );
	rsrdl_2.get_probability_for_rotamer( 0.0, 0.0, 1 );
	rsrdl_3.get_probability_for_rotamer( 0.0, 0.0, 1 );
	rsrdl_4.get_probability_for_rotamer( 0.0, 0.0, 1 );
	srsrdl_1.get_probability_for_rotamer( 0.0, 0.0, 1 );
	srsrdl_2.get_probability_for_rotamer( 0.0, 0.0, 1 );

	rsrdl_1.get_rotamer( 0.0, 0.0, 1 );
	rsrdl_2.get_rotamer( 0.0, 0.0, 1 );
	rsrdl_3.get_rotamer( 0.0, 0.0, 1 );
	rsrdl_4.get_rotamer( 0.0, 0.0, 1 );
	srsrdl_1.get_rotamer( 0.0, 0.0, 1 );
	srsrdl_2.get_rotamer( 0.0, 0.0, 1 );

	/*
	rsrdl_1
	rsrdl_2
	rsrdl_3
	rsrdl_4
	srsrdl_1
	srsrdl_2
	*/

}


} // namespace dunbrack
} // namespace scoring
} // namespace core

