// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/enzdes/DesignSilentStruct.cc
///
/// @brief protein silent-file structures for designs, also contains functionality to query
/// @brief a given pose and extract some more data to print
/// @author Florian Richter



// C++ Headers
// AUTO-REMOVED #include <fstream>
#include <iostream>
// AUTO-REMOVED #include <utility>
#include <vector>
#include <map>

#include <core/pose/Pose.hh>

// mini headers
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <devel/enzdes/DesignSilentStruct.hh>
#include <basic/MetricValue.hh>
// Auto-header: duplicate removed #include <core/pose/Pose.hh>

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
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
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
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0_bool.hh>
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
#include <utility/io/mpistream.hh>
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
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
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
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
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
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/pool/poolfwd.hpp>

//#include <devel/PoseMetricCalculators/InterfaceDeltaEnergeticsCalculator.hh> //needed for type checking

using namespace core;
using namespace core::io::silent;
using namespace ObjexxFCL::fmt;

namespace devel {
namespace enzdes {

static basic::Tracer tr("devel.enzdes.DesignSilentStruct");


DesignSilentStruct::DesignSilentStruct(
	core::pose::Pose const & pose,
	std::string tag,
  bool const add_in,
  bool const onlyadd_in )
{
	decoy_tag( tag );
	sequence( pose.sequence() );
	nres( pose.total_residue() );

  print_additional_ = add_in;
  print_only_additional_ = onlyadd_in;
}


DesignSilentStruct::DesignSilentStruct(
  pose::Pose const & pose,
  std::string tag, // = "empty_tag",
  utility::vector1<Size> const & spec_res_in,
  utility::vector1< std::string > const & rel_score_in,
  bool const add_in,
  bool const onlyadd_in )
{

  std::map< Size, utility::vector1< std::pair< std::string, std::string > > > dummycalcs;
  dummycalcs.clear();
  DesignSilentStruct( pose, tag, spec_res_in, rel_score_in, dummycalcs, add_in, onlyadd_in);
}


DesignSilentStruct::DesignSilentStruct(
  pose::Pose const & pose,
  std::string tag,
  utility::vector1<Size> const & spec_res_in,
  utility::vector1< std::string > const & rel_score_in,
  std::map< Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators,
  bool const add_in,
  bool const onlyadd_in
) {

  if( (add_in == false) && (onlyadd_in == true) ){
    std::cerr << "do you want a silent structure or not?" << std::endl;
    utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
  }

  decoy_tag( tag );
  sequence( pose.sequence() );
  nres( pose.total_residue() );

  print_additional_ = add_in;
  print_only_additional_ = onlyadd_in;

  if( (add_in == true) || (onlyadd_in == true) ){
    calculate_additional_info( pose, spec_res_in, rel_score_in, calculators );
  }

}

void
DesignSilentStruct::print_header( std::ostream& out ) const {

  out << LJ( 50, "#Design_name");
  for ( utility::vector1< SilentEnergy >::const_iterator it = additional_info_silent_energy_.begin();
	it != additional_info_silent_energy_.end(); ++it ){
    out << " " << A( it->width(), it->name() );
  }

  out << "\n";

}

void
DesignSilentStruct::print_scores( std::ostream& out ) const
{

  if( !print_only_additional_ ){
    parent::print_scores( out );
  }

  if( print_additional_ ){
    print_additional_info( out );
  }

}

void
DesignSilentStruct::add_to_additional_silent_energies(
	utility::vector1< core::io::silent::SilentEnergy > const & silent_Es)
{

  for ( utility::vector1< SilentEnergy >::const_iterator it = silent_Es.begin();
	it != silent_Es.end(); ++it ){
		additional_info_silent_energy_.push_back( *it );
	}

}

void
DesignSilentStruct::print_additional_info( std::ostream& out) const
{
  int precision = 2; //number of digits after decimal place

  out << LJ( 50, decoy_tag() );

  for ( utility::vector1< SilentEnergy >::const_iterator it = additional_info_silent_energy_.begin();
	it != additional_info_silent_energy_.end(); ++it ){
    out << " " << F( it->width(), precision, it->value() );
  }

  out << "\n";

}

void
DesignSilentStruct::calculate_additional_info(
  pose::Pose const & pose,
  utility::vector1<Size> const & special_res,
  utility::vector1< std::string > const & score_terms,
  std::map< Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators )
{

	bool separate_out_constraints = false;
	utility::vector1< std::string >::const_iterator cstfind = find( score_terms.begin(), score_terms.end(),"all_cst");
	if( cstfind != score_terms.end() ) separate_out_constraints = true;

  //first write out the relevant score terms for the pose total
  for( utility::vector1< std::string >::const_iterator sco_it = score_terms.begin();
       sco_it != score_terms.end(); ++sco_it ){
    std::string sco_name = *sco_it;
    int width = std::max( 10, (int) sco_name.length() + 3 );

    SilentEnergy new_se;
    if( *sco_it == "all_cst" ) {
      new_se = SilentEnergy ( sco_name, sum_constraint_terms(pose, -1 ), 1, width);
    }
		else if( separate_out_constraints && ( *sco_it == "total_score" ) ){
			core::Real desired_value = pose.energies().total_energies()[ core::scoring::score_type_from_name( *sco_it ) ] - sum_constraint_terms(pose, -1 );
			new_se = SilentEnergy(sco_name, desired_value, 1, width);
		}
    else{
      new_se = SilentEnergy ( sco_name, pose.energies().total_energies()[ core::scoring::score_type_from_name( *sco_it ) ] * pose.energies().weights()[ core::scoring::score_type_from_name( *sco_it ) ], 1 ,width);
    }
    additional_info_silent_energy_.push_back( new_se );

  }

	//pose metric calculators for pose total
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > >::const_iterator totcalc_it = calculators.find( 0 );
	if( totcalc_it != calculators.end() ){

    utility::vector1< std::pair< std::string, std::string > > const & tot_calculators = totcalc_it->second;
		for( utility::vector1< std::pair< std::string, std::string > >::const_iterator calc_it = tot_calculators.begin();
				 calc_it != tot_calculators.end(); ++calc_it){

			std::string calc_name = "tot_" + calc_it->first;
			int width = std::max( 10, (int) calc_name.length() + 3 );

			core::Real calc_value;

			//following lines fairly hacky, but don't know a better solution at the moment
			if( calc_it->first == "hbond_pm" || calc_it->first == "burunsat_pm" || calc_it->first == "NLconts_pm" ){
				basic::MetricValue< core::Size > mval_size;
				pose.metric( calc_it->first, calc_it->second, mval_size );
				calc_value = mval_size.value();
			}
			else{
				basic::MetricValue< core::Real > mval_real;
				pose.metric( calc_it->first, calc_it->second, mval_real );
				calc_value = mval_real.value();
			}

			SilentEnergy new_se( calc_name, calc_value, 1, width);
			additional_info_silent_energy_.push_back( new_se );
		}

	}

	//done with pose totals

  //then write out the relevant scoreterms (and potentially pose metrics) for each of the special residues
	Size spec_res_counter(0);
  for( utility::vector1<Size>::const_iterator res_it = special_res.begin(); res_it != special_res.end(); res_it++ ){

    spec_res_counter++;
    //for convenience, the sequence number of the residue will be written out
    std::stringstream temp;
    temp << spec_res_counter;
    std::string spec_res_name = "SR_" + temp.str();
    //std::cerr << "name for res " << *res_it << " is " << spec_res_name ;
    SilentEnergy res_name(spec_res_name, *res_it, 1, 10);
    additional_info_silent_energy_.push_back( res_name );

    for( utility::vector1< std::string >::const_iterator sco_it = score_terms.begin();
       sco_it != score_terms.end(); ++sco_it ){

      std::string sco_name = spec_res_name + "_" +  *sco_it ;
      int width = std::max( 10, (int) sco_name.length() + 3 );

      SilentEnergy new_se;
      if( *sco_it == "all_cst" ) {
				new_se = SilentEnergy ( sco_name, sum_constraint_terms(pose, *res_it), 1, width);
      }
			else if( separate_out_constraints && ( *sco_it == "total_score" ) ){
				core::Real desired_value = pose.energies().residue_total_energies( *res_it )[ core::scoring::score_type_from_name( *sco_it ) ] - sum_constraint_terms(pose, *res_it );
				new_se = SilentEnergy(sco_name, desired_value, 1, width);
			}
      else{
      new_se = SilentEnergy ( sco_name, pose.energies().residue_total_energies( *res_it )[ core::scoring::score_type_from_name( *sco_it ) ] * pose.energies().weights()[ core::scoring::score_type_from_name( *sco_it ) ], 1 ,width);
      }

      additional_info_silent_energy_.push_back( new_se );
    }//loop over relevant scoreterms


    //if there are calculators that need to be evaluated for this residue, let's do that now
    //note: section still under development, right now only calculators that return reals are supported
    std::map< Size, utility::vector1< std::pair< std::string, std::string > > >::const_iterator res_calc_it = calculators.find( *res_it );
    if( res_calc_it != calculators.end() ){

      utility::vector1< std::pair< std::string, std::string > > calculators_this_res = res_calc_it->second;
      for( utility::vector1< std::pair< std::string, std::string > >::iterator calc_it = calculators_this_res.begin();
					 calc_it != calculators_this_res.end(); calc_it++ ){

				std::string res_calc_name = spec_res_name + "_" + calc_it->first;
				int width = std::max( 10, (int) res_calc_name.length() + 3 );

				core::Real calc_value;

				basic::MetricValue< core::Real > mval_real;
				basic::MetricValue< utility::vector1< core::Size > >mval_sizevec;
				basic::MetricValue< utility::vector1< core::Real > >mval_realvec;

				//following lines fairly hacky, but don't know a better solution at the moment
				if( ( calc_it->first == "hbond_pm") || ( calc_it->first == "burunsat_pm") ){
					pose.metric( calc_it->first, calc_it->second, mval_sizevec );
					calc_value = mval_sizevec.value()[*res_it];
				}
				else if( (calc_it->first == "pstat_pm") || (calc_it->first == "nlpstat_pm" ) ) {
					pose.metric( calc_it->first, calc_it->second, mval_realvec );
					calc_value = mval_realvec.value()[*res_it];
				}
				else {
					pose.metric( calc_it->first, calc_it->second, mval_real );
					calc_value = mval_real.value();
				}
				//std::cerr << " hehe, just executed pose metric for " << calc_it->first << " calculator   ";


				//special case: if this is an interface calculation, we do not want to include constraint terms
				//hacky at the moment, haven't figured out yet how to do this really clean
				if( separate_out_constraints && ( calc_it->first == "interf_E_1_2") ){
					calc_value = calc_value - ( 2 * sum_constraint_terms(pose, *res_it) );
				}
				SilentEnergy new_se( res_calc_name, calc_value, 1, width);
				additional_info_silent_energy_.push_back( new_se );

      } // for calculators this res

    }// if calculators for this res
    //else std::cerr << "resi " << *res_it << " has no calcs." << std::endl;

  }//loop over special residues

}

core::Real
DesignSilentStruct::sum_constraint_terms( pose::Pose const & pose, int which_res){

  using namespace core::scoring;

  EnergyMap all_weights = pose.energies().weights();
  EnergyMap scores;

  if( which_res == -1){ //means we want to know stuff for the whole pose
    scores = pose.energies().total_energies();
  }
  else{ scores = pose.energies().residue_total_energies( which_res ); }

  return scores[ coordinate_constraint ] * all_weights[ coordinate_constraint ] + scores[atom_pair_constraint] * all_weights[ atom_pair_constraint] +
    scores[ angle_constraint ] * all_weights[ angle_constraint ] + scores[ dihedral_constraint ] * all_weights[ dihedral_constraint ];

}



} //namespace enzdes
} //namespace devel
