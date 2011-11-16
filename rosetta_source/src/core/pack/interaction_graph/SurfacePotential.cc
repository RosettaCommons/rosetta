// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/SurfacePotential.cc
/// @brief  Class which keeps reads the residue hydrophobic ASA database file and calculates surface residue energies.
/// @author Ron Jacak

#include <core/pack/interaction_graph/SurfacePotential.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


#include <core/graph/DisjointSets.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pack/interaction_graph/RotamerDots.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>

// AUTO-REMOVED #include <core/scoring/sasa.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <core/types.hh>

// ObjexxFCL Headers

// Numeric Headers

// Utility Headers
#include <utility/io/izstream.hh>

// C++ Headers
#include <sstream>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
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
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
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
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
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
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
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
#include <core/scoring/EnergyMap.hh>
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
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
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
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
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
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/poolfwd.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers




namespace core {
namespace pack {
namespace interaction_graph {

//#define FILE_DEBUG 1

static basic::Tracer TR("core.pack.interaction_graph.SurfacePotential");

const core::Size SurfacePotential::MAX_PATCH_SURFACE_AREA = 1200;
const core::Real SurfacePotential::MAX_SURFACE_ENERGY = 25.0;
const core::Size SurfacePotential::SURFACE_EXPOSED_CUTOFF = 20;
const core::Real SurfacePotential::INTERACTION_RADIUS = 10.0;
const core::Size SurfacePotential::SURFACE_SCORE_BIN_SIZE = 25;
const core::Size SurfacePotential::BURIED_RESIDUE_NO_HSASA_CUTOFF = 24;

const core::Size SurfacePotential::MAX_HPATCH_AREA = 1000;
const core::Real SurfacePotential::MAX_HPATCH_SCORE = 100.0;
const core::Size SurfacePotential::HPATCH_SCORE_BIN_SIZE = 50;

/// @brief set initial value as no instance
SurfacePotential* SurfacePotential::instance_( 0 );


/// @brief static function to get the instance of (pointer to) this singleton class
SurfacePotential* SurfacePotential::get_instance() {
	if ( instance_ == 0 )
		instance_ = new SurfacePotential();
	return instance_;
}

/// @brief private constructor to guarantee the singleton
SurfacePotential::SurfacePotential() {
	read_average_res_hASA_database_file();
	read_hASA_score_database_file();

	read_hpatch_score_database_file();
}


/// @brief Reads in the database file which contains average residue hydrophobic accessible surface areas.
void SurfacePotential::read_average_res_hASA_database_file() {

	utility::io::izstream residue_hASA_ifstream;
	basic::database::open( residue_hASA_ifstream, "average_hASA_by_res_and_neighbor.txt" );

	std::string str_restype;
	Real lte10_asa = 0.0;
	Real lte13_asa = 0.0;
	Real lte16_asa = 0.0;
	Real lte20_asa = 0.0;
	Real lte24_asa = 0.0;

	res_to_average_hASA_.resize( chemical::num_canonical_aas );

	while ( !residue_hASA_ifstream.eof() ) {
		residue_hASA_ifstream >> str_restype >> lte10_asa >> lte13_asa >> lte16_asa >> lte20_asa >> lte24_asa;

		// store the mean hASA values in a vector for faster lookup
		// This way, rather than constructing a key each time we need to get the information on a particular AA type
		// we can just index right into the vector. We waste some memory on storing the same values in the vector
		// but this class is a Singleton and it's only 48 floats and the size of the vectors.
		chemical::AA aa_type = chemical::aa_from_name( str_restype );
		res_to_average_hASA_[ aa_type ].resize( 24, 0.0 ); // we have average hASA values for up to 24 nbs

		for ( Size ii=1; ii <= 10; ++ii ) {  // for nbs 1-10
			res_to_average_hASA_[ aa_type ][ ii ] = lte10_asa;
		}
		for ( Size ii=11; ii <= 13; ++ii ) { // for nbs 11-13
			res_to_average_hASA_[ aa_type ][ ii ] = lte13_asa;
		}
		for ( Size ii=14; ii <= 16; ++ii ) {  // for nbs 14-16
			res_to_average_hASA_[ aa_type ][ ii ] = lte16_asa;
		}
		for ( Size ii=17; ii <= 20; ++ii ) {  // for nbs 17-20
			res_to_average_hASA_[ aa_type ][ ii ] = lte20_asa;
		}
		for ( Size ii=21; ii <= 24; ++ii ) {  // for nbs 21-24
			res_to_average_hASA_[ aa_type ][ ii ] = lte24_asa;
		}

	}

#ifdef FILE_DEBUG
	// quick check to make sure we stored the right things
	for ( Size ii=1; ii <= chemical::num_canonical_aas; ++ii ) {
		TR << chemical::name_from_aa( (chemical::AA)ii ) << ": [ ";
		for ( Size jj=1; jj <= res_to_average_hASA_[ ii ].size(); ++jj ) {
			TR << res_to_average_hASA_[ ii ][ jj ] << ", ";
		}
		TR << "]" << std::endl;
	}
#endif

}


/// @brief Reads in the database file which contains the scores for a distribution of patch sizes.
///
/// @detailed
/// Not assuming any particular length to the database file so that if I want to increase
/// the maximum of the distribution or shrink it, the vector will dynamically resize to what it needs
/// to be.
///
void SurfacePotential::read_hASA_score_database_file() {

	utility::io::izstream hASA_score_ifstream;
	basic::database::open( hASA_score_ifstream, "surface_score.txt" );

	Real amount_hASA = 0.0;
	Real score = 0.0;
	Size num_fields = 20; // keep scores for nb counts of 1, 2, 3 ... x (where x is # fields)

	while ( !hASA_score_ifstream.eof() ) {
		hASA_score_ifstream >> amount_hASA;

		utility::vector1< Real > v( num_fields, 0.0 );
		for ( Size ii=1; ii <= num_fields; ++ii ) {
			hASA_score_ifstream >> score;
			v[ ii ] = score;
		}
		hASA_to_score_.push_back( v );
	}

#ifdef FILE_DEBUG
	TR << "patch energies: [ " << std::endl;
	for ( Size ii=0; ii < hASA_to_score_.size(); ++ii ) {
		TR << "patch size: " << ii * 25 << " [ ";
		for ( Size jj=1; jj <= hASA_to_score_[ ii ].size(); ++jj ) {
			TR << hASA_to_score_[ ii ][ jj ] << ", ";
		}
		TR << "]" << std::endl;
	}
	TR << "]" << std::endl;
#endif
}


/// @brief Reads in the database file for the hpatch score, yet another version of the surface energy.
///
/// @detailed
/// Not assuming any particular length to the database file so that if I want to increase
/// the maximum of the distribution or shrink it, the vector will dynamically resize to what it needs
/// to be.
///
void SurfacePotential::read_hpatch_score_database_file() {

	utility::io::izstream hpatch_score_ifstream;
	basic::database::open( hpatch_score_ifstream, "hpatch_score.txt" );

	Real patch_area = 0.0;
	Real score = 0.0;

	while ( !hpatch_score_ifstream.eof() ) {
		hpatch_score_ifstream >> patch_area >> score;
		patcharea_to_score_.push_back( score );
	}

#ifdef FILE_DEBUG
	TR << "patch area scores: [ " << std::endl;
	for ( Size ii=0; ii < patcharea_to_score_.size(); ++ii ) {
		TR << patcharea_to_score_[ ii ] << ", ";
	}
	TR << "]" << std::endl;
#endif

}


///
/// @begin SurfacePotential::average_residue_hASA
///
/// @brief
/// Returns the average surface energy for the given residue type and number of neighbors.
///
Real SurfacePotential::average_residue_hASA( chemical::AA aa_type, Size num_nbs ) {

	// these checks will only run in debug mode builds
#ifdef FILE_DEBUG
	if ( num_nbs > BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
		std::cout << "Number of neighbors (" << num_nbs << ") outside bounds in SurfacePotential::average_residue_hASA." << std::endl;
	}
	if ( aa_type > chemical::num_canonical_aas ) {
		std::cout << "aatype (" << aa_type << ") outside of canonical 20 in SurfacePotential::average_residue_hASA." << std::endl;
	}
#endif
	assert( num_nbs <= BURIED_RESIDUE_NO_HSASA_CUTOFF );
	if ( aa_type > chemical::num_canonical_aas ) { return 0.0; }
	return res_to_average_hASA_[ aa_type ][ num_nbs ];

}

///
/// @begin SurfacePotential::hASA_patch_energy
///
/// @brief
/// Returns the energy for a given patch size.  The calling function must ensure that an out-of-bounds error will not occur.
///
Real SurfacePotential::hASA_patch_energy( Real patch_area, Size num_nbs ) {

#ifdef FILE_DEBUG
	if ( patch_area > MAX_PATCH_SURFACE_AREA ) {
		std::cout << "patch_area (" << patch_area << ") greater than MAX_PATCH_SURFACE_AREA in SurfacePotential::hASA_patch_energy." << std::endl;
	}
#endif
	assert( patch_area <= MAX_PATCH_SURFACE_AREA );
	return hASA_to_score_[ (Size)(patch_area / SURFACE_SCORE_BIN_SIZE ) ][ num_nbs ];
}


///
/// @begin SurfacePotential::hpatch_score
///
/// @brief
/// Returns the score for a given patch size.  The calling function must ensure that an out-of-bounds error will not occur.
///
Real SurfacePotential::hpatch_score( Real patch_area ) {

#ifdef FILE_DEBUG
	if ( patch_area > MAX_HPATCH_AREA ) {
		std::cout << "patch_area (" << patch_area << ") greater than MAX_HPATCH_AREA in SurfacePotential::hpatch_score." << std::endl;
	}
#endif
	assert( patch_area <= MAX_HPATCH_AREA );
	return patcharea_to_score_[ (Size)(patch_area / HPATCH_SCORE_BIN_SIZE) ];
}


///
/// @begin SurfacePotential::compute_residue_surface_energy
///
/// @brief
/// Calculates the surface energy for a single residue within a Pose object. Used only by the RotamerSet_::compute_one_body_energy_maps
/// function (which, in turn, is only used by the optE protocol).  Nowhere else in mini is this function used.
///
void SurfacePotential::compute_residue_surface_energy( conformation::Residue const & rsd, pose::Pose const & pose, scoring::EnergyMap & emap,
	Size resid, utility::vector1< Size > num_neighbors_ ) {

	// our definition of surface residue is that the residue has fewer than 16 neighbors
	if ( !( pose.energies().residue_neighbors_updated() ) || num_neighbors_[ resid ] > SURFACE_EXPOSED_CUTOFF ) {
		emap[ scoring::surface ] = 0.0;
		return;
	}

	// add this residues hASA to the total
	// this residue must have <= SURFACE_EXPOSED_CUTOFF number of neighbors to get here
	Real total_hASA = 0.0;
	total_hASA += average_residue_hASA( rsd.aa(), num_neighbors_[ resid ] );
#ifdef FILE_DEBUG
	TR << rsd.aa() << resid << " (hASA: " << average_residue_hASA( rsd.aa(), num_neighbors_[ resid ] ) << ", nbs: " << num_neighbors_[ resid ] << ") neighbors... ";
#endif

	// now add the hASA of all se residues within xA (nbr_atom - nbr_atom) distance. Use the nbr_atom coordinates of the passed
	// in (being considered) rotamer, not of the wild type sequence rotamer.  In a small percent of the cases, using the wild type
	// nbr_atom will give a different count than when using the new rotamer nbr_atom position.
	Real distanceBetweenAtoms = 0.0;
	for ( Size res2_position = 1; res2_position < pose.total_residue(); ++res2_position ) {

		if ( resid == res2_position ) { continue; }
		conformation::Residue const & rsd2 = pose.residue( res2_position );

		distanceBetweenAtoms = rsd.xyz( rsd.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );

		// first, have to check again if these two residues are neighbors
		if ( distanceBetweenAtoms <= INTERACTION_RADIUS ) {

			// and finally, check to make sure it's surface exposed
			if ( num_neighbors_[ res2_position ] > BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
				continue; // if it has so many neighbors it probably doesn't add anything to the hSASA (and it would cause out-of-bounds errors below)
			}
			// passed all checks
			total_hASA += average_residue_hASA( rsd2.aa(), num_neighbors_[ res2_position ] );
#ifdef FILE_DEBUG
			TR << rsd2.aa() << res2_position << " " << average_residue_hASA( rsd2.aa(), num_neighbors_[ res2_position ] ) << ", ";
#endif
		}
	}
#ifdef FILE_DEBUG
	TR << std::endl;
#endif

	// now that we know how many surface-exposed neighbors res1 has, get the surface energy for that value
	if ( total_hASA > MAX_PATCH_SURFACE_AREA ) {
		emap[ scoring::surface ] = MAX_SURFACE_ENERGY;
	} else {
		emap[ scoring::surface ] = hASA_patch_energy( total_hASA, num_neighbors_[ resid ] );
	}

#ifdef FILE_DEBUG
	TR << "compute_residue_surface_energy: calculated total_hASA: " << total_hASA << ", surface energy: " << emap[ scoring::surface ] << std::endl;
#endif

}


///
/// @begin compute_pose_surface_energy
///
/// @brief
/// helper method for computing surface score. in the optE protocol we don't care about the total vs. residue
/// level surface scores. so just call that function but discard the values for those variables.
///
void SurfacePotential::compute_pose_surface_energy( pose::Pose const & pose, Real & surface_energy_ ) {

	utility::vector1< Size > num_neighbors_;
	utility::vector1< Real > res_level_energies_;

	compute_pose_surface_energy( pose, surface_energy_, res_level_energies_ );

	return;
}

void SurfacePotential::compute_pose_surface_energy( pose::Pose const & pose, Real & total_surface_energy_,
	utility::vector1< Real > & residue_surface_energy_ ) {

	// the pose has to have been scored at this point for this method to work.  since I can't force a score eval here,
	// return 0.0 for everything if it hasn't been.
	if ( ! pose.energies().residue_neighbors_updated() ) {
		total_surface_energy_ = 0.0;
		residue_surface_energy_.clear();
		return;
	}

	// resize the per-residue surface energy vector
	residue_surface_energy_.clear();
	residue_surface_energy_.resize( pose.n_residue(), 0.0 );

	utility::vector1< Size > num_neighbors_( pose.n_residue(), 0 );

	// first, we need to init the num neighbors array (either using the tenA nb graph or by counting manually)
	for ( core::Size res1_position = 1; res1_position <= pose.n_residue(); ++res1_position ) {
		core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

		// set the number of neighbors vector for later output
		num_neighbors_[ res1_position ] = tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self();
	}

	// need symmetry info to ignore the duplicated units
	core::conformation::symmetry::SymmetryInfoCOP symm_info(NULL);
	if( core::pose::symmetry::is_symmetric( pose ) ) {
		symm_info = (dynamic_cast<const core::conformation::symmetry::SymmetricConformation &>(pose.conformation()).Symmetry_Info());
	}

	// now, we have to loop over all residues and find exposed residues (we can use the num neighbors array to determine
	// which ones are surface exposed)

	for ( core::Size res1_position = 1; res1_position <= pose.n_residue(); ++res1_position ) {

		if( pose.residue( res1_position ).aa() > core::chemical::num_canonical_aas ) continue;
		if( symm_info && !symm_info->bb_is_independent(res1_position) ) continue;

		// reset the counter
		Real total_hASA = 0.0;

		// our definition of surface residue is that the residue has fewer than 16 neighbors. we only assign surface
		// scores to residues with fewer than this number of nbs. but when getting the patch area, neighboring residues
		// can have more nbs than this and still contribute to patch area.
		if ( num_neighbors_[ res1_position ] > SURFACE_EXPOSED_CUTOFF ) {
			continue;
		}

		// passed the surface-exposed check...

		total_hASA += average_residue_hASA( pose.residue( res1_position ).aa(), num_neighbors_[ res1_position ] );
#ifdef FILE_DEBUG
		TR << pose.residue( res1_position ).aa() << res1_position
			<< " (hASA: " << average_residue_hASA( pose.residue( res1_position ).aa(), num_neighbors_[ res1_position ] )
			<< ", nbs: " << num_neighbors_[ res1_position ] << ") neighbors... ";
#endif

		//TR << "Neighbors of residue " << pose.residue( res1_position ).name3() << " " << res1_position << " include " << std::endl;

		// for every Edge in the neighbor graph, figure out if that residue is surface exposed
		//for ( core::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_begin(),
		//	eli_end = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_end(); eli != eli_end; ++eli ) {

			// save the value to simplify code ahead
			//int res2_position = (*eli)->get_other_ind( res1_position );

			// get the other node for this edge, so pass in the res1 node to this method
			//TR << pose.residue( res2_position ).name3() << " " << res2_position << std::endl;

		conformation::Residue const & rsd1 = pose.residue( res1_position );
		Real distanceBetweenAtoms = 0.0;

		for ( Size res2_position = 1; res2_position < pose.total_residue(); ++res2_position ) {
			if( pose.residue( res2_position ).aa() > core::chemical::num_canonical_aas ) continue;
			if( symm_info && !symm_info->bb_is_independent(res2_position) ) continue;

			if ( res2_position == res1_position ) { continue; }
			conformation::Residue const & rsd2 = pose.residue( res2_position );

			distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );

			// first, have to check again if these two residues are neighbors
			if ( distanceBetweenAtoms <= INTERACTION_RADIUS ) {

				// ok, is it surface-exposed, too?
				if ( num_neighbors_[ res2_position ] > BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
					continue; // if it has so many neighbors it probably doesn't add anything to the hSASA (and it would cause out-of-bounds errors below)
				}
				// passed all checks
				total_hASA += average_residue_hASA( rsd2.aa(), num_neighbors_[ res2_position ] );
#ifdef FILE_DEBUG
				TR << rsd2.aa() << res2_position << " " << average_residue_hASA( rsd2.aa(), num_neighbors_[ res2_position ] ) << ", ";
#endif
			}
		}
#ifdef FILE_DEBUG
		TR << std::endl;
#endif

		// now that we know how many surface-exposed neighbors res1 has, get the surface energy for that value
		if ( total_hASA > MAX_PATCH_SURFACE_AREA ) {
			residue_surface_energy_[ res1_position ] = MAX_SURFACE_ENERGY;
		} else {
			residue_surface_energy_[ res1_position ] = hASA_patch_energy( total_hASA, num_neighbors_[ res1_position ] );
		}

#ifdef FILE_DEBUG
		TR << "compute_pose_surface_energy: calculated residue total_hASA: " << total_hASA << ", surface energy: "
			<< residue_surface_energy_[ res1_position ] << std::endl;
#endif

	} // end second res1 loop


	total_surface_energy_ = 0.0;
	for ( Size ii=1; ii < residue_surface_energy_.size(); ++ii ) {
		total_surface_energy_ += residue_surface_energy_[ii];
	}

#ifdef FILE_DEBUG
	TR << "compute_pose_surface_energy: calculated surface energy: " << total_surface_energy_ << ", residue_surfaceE: [ ";
	for ( Size ii=1; ii <= residue_surface_energy_.size(); ++ii ) {
		TR << residue_surface_energy_[ ii ] << ", ";
	}
	TR << "]" << std::endl;
#endif

}


Real
SurfacePotential::compute_pose_hpatch_score(
	pose::Pose const & pose
)
{
	core::Real total_hpatch_score;
	std::map< core::Size, std::pair< core::Real, core::Real > > patch_scores;
	std::map< Size,utility::vector1< id::AtomID > > atoms_in_patches;

	compute_pose_hpatch_score( pose, total_hpatch_score, patch_scores, atoms_in_patches );
	return total_hpatch_score;
}


///
/// @begin compute_pose_hpatch_score
///
/// @brief
/// Uses the src/core/pack/interaction_graph/RotamerDots classes to determine exact SASAs and then uses a graph-based
/// approach to find all the exposed hydrophobic patches on the surface of the pose. Uses the scores in the file
/// hpatch_score.txt to assign a score to each patch and puts the score into the passed in Real.
///
/// Note: If a Pose object has overlapping atoms anywhere, then this function will fail with the following error:
///
/// sin_cos_range ERROR: nan is outside of [-1,+1] sin and cos value legal range
/// sin_cos_range ERROR: nan is outside of [-1,+1] sin and cos value legal range
/// ERROR:: Exit from: src/numeric/trig.functions.hh line: 117
///
///
void SurfacePotential::compute_pose_hpatch_score(
	pose::Pose const & pose,
	core::Real & total_hpatch_score_,
	std::map< core::Size, std::pair< core::Real, core::Real > > & patch_scores_,
	std::map< Size,utility::vector1< id::AtomID > > & atoms_in_patches_
)
{

	using namespace core;
	using namespace core::id;
	using namespace core::pack::interaction_graph;

	// an atomID map is needed for the calc_per_atom_sasa method; it stores the actual calculated sasa for every atom
	core::id::AtomID_Map< core::Real > atom_sasa;
	core::pose::initialize_atomid_map( atom_sasa, pose, (core::Real)0.0 ); // initialize to 0.0 for "not computed"

	utility::vector1< RotamerDotsOP > rdots( pose.total_residue() );
	utility::vector1< InvRotamerDotsOP > invdots( pose.total_residue() );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		rdots[ii] = new RotamerDots( pose.residue(ii).clone(), true /* exclude H's */, true /* use expanded polar atom radii */);
		invdots[ii] = new InvRotamerDots();
	}

	// PointGraph is a one-way graph, which makes it somewhat annoying for iterating over neighbors of a certain
	// position. Only edges to higher-indexed nodes exist. So instead, make a graph which has all the edges at every
	// node to simplify iterating over all neighboring edges.
	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); //create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); //create vertices

	Distance const max_pair_radius = pose::pose_max_nbr_radius( pose );
	Distance const max_ep_radius = rdots[1]->max_atom_radius() + 2 * RotamerDots::probe_radius_;
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, max_pair_radius + max_pair_radius + max_ep_radius /* Angstrom cutoff */ ); //create edges

	// increment the self and residue-residue overlap for each residue
	for ( Size ii = 1; ii <= pose.total_residue(); ++ ii ){
		rdots[ ii ]->increment_self_overlap();
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex(ii ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			Size jj = edge_iter->upper_vertex();
			rdots[ ii ]->increment_both( *rdots[ jj ] );
		}
	}

	// we need to know how many heavy atoms are in the pose before we can construct the disjoint sets object
	Size heavyatom_count = 0;

	for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & rsd = pose.residue( ii );
		heavyatom_count += rsd.nheavyatoms();
	}
	
	graph::DisjointSets ds( heavyatom_count );
	utility::vector1< id::AtomID > ds_index_2_atomid( heavyatom_count );

	// create an AtomID map that will convert an atom in some residue into a DisjointSets index
	id::AtomID_Map< Size > atom_2_ds_index;
	atom_2_ds_index.resize( pose.total_residue() );
	Size ds_index = 1;

	for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & rsd = pose.residue( ii );
		atom_2_ds_index.resize( ii, rsd.nheavyatoms(), 0 );
		for ( Size jj=1; jj <= rsd.nheavyatoms(); ++jj ) {
			id::AtomID const atomid( jj, ii );
			atom_2_ds_index[ atomid ] = ds_index;
			ds_index_2_atomid[ ds_index ] = atomid;
			ds_index++;
		}
	}

	// now iterate over all residues of the pose and find all intra- and inter-residue atom-atom overlaps
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		invdots[ ii ]->setup_from_rotamer_dots( *rdots[ ii ] );
	}

#ifdef FILE_DEBUG
	TR << "compute_pose_hpatch_score(): finding intra- and inter-residue overlaps" << std::endl;
#endif

	std::string carbon_atom = "C", sulfur_atom = "S";
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		// don't try to find overlaps with non-protein residues; this can cause problems.
		if ( ! pose.residue( ii ).is_protein() ) continue;

		// 1. self overlap; iterate over all heavyatoms of residue ii
		for ( Size iia = 1, iia_end = pose.residue_type(ii).nheavyatoms(); iia <= iia_end; ++iia ) {

			// check if iia atom is buried with expanded polar atoms
			if ( rdots[ ii ]->get_atom_sasa( iia ) == 0 )
				continue;

			// check if jj atom is a hydrophobic atom
			std::string const & iia_elem = pose.residue_type(ii).atom_type( iia ).element();
			if ( iia_elem != carbon_atom && iia_elem != sulfur_atom ) // doing char comparison is much faster than string comparison
				continue;

			Real const iia_rad = rdots[ ii ]->get_atom_radius( iia ) + RotamerDots::probe_radius_;
			Vector const & iia_xyz = rdots[ ii ]->rotamer()->xyz( iia );

			// only have to iterate over higher-indexed heavyatoms of iia
			for ( Size iib = iia+1; iib <= iia_end; ++iib ) {

				// check if kk atom is buried with expanded polar atoms
				if ( rdots[ ii ]->get_atom_sasa( iib ) == 0 )
					continue;

				// check if iib atom is a hydrophobic atom
				std::string const & iib_elem = pose.residue_type( ii ).atom_type( iib ).element();
				if ( iib_elem != carbon_atom && iib_elem != sulfur_atom ) // doing char comparison is much faster than string comparison
					continue;

				Real const iib_rad = rdots[ ii ]->get_atom_radius( iib ) + RotamerDots::probe_radius_;
				Vector const & iib_xyz = rdots[ ii ]->rotamer()->xyz( iib );

				if ( iia_xyz.distance_squared( iib_xyz ) < (iia_rad + iib_rad) * (iia_rad + iib_rad) ) {
#ifdef FILE_DEBUG
					/*conformation::Residue const & ii_rsd = pose.residue( ii );
					TR << "compute_pose_hpatch_score(): overlapping atom pair: "
						<< ii_rsd.aa() << " " << ii << "/" << utility::trim( ii_rsd.atom_name( iia ) ) << " - "
						<< ii_rsd.aa() << " " << ii << "/" << utility::trim( ii_rsd.atom_name( iib ) ) << std::endl;*/
#endif

					if ( invdots[ ii ]->atom_overlap_is_exposed( iia, iib ) ) {
						Size iia_dsid = atom_2_ds_index( AtomID( iia, ii ) );
						Size iib_dsid = atom_2_ds_index( AtomID( iib, ii ) );

						// if the two atoms aren't in the same connected component, call union on them
						if ( ds.ds_find( iia_dsid ) != ds.ds_find( iib_dsid ) ) {
							ds.ds_union( iia_dsid, iib_dsid );
						}
					}

				}
			}
		} // end for loop for self overlap

		/// 2. upper neighbors of ii, we'll use jj
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex( ii ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {

			Size jj = edge_iter->upper_vertex();
			// don't try to find overlaps with non-protein residues; this can cause problems.
			if ( ! pose.residue( jj ).is_protein() ) continue;

			for ( Size iia = 1; iia <= pose.residue_type( ii ).nheavyatoms(); ++iia ) {

				if ( rdots[ ii ]->get_atom_sasa( iia ) == 0 )
					continue;

				std::string const & iia_elem = pose.residue_type( ii ).atom_type( iia ).element();
				if ( iia_elem != carbon_atom && iia_elem != sulfur_atom )
					continue;

				Real const iia_rad = rdots[ ii ]->get_atom_radius( iia ) + RotamerDots::probe_radius_;
				Vector const & iia_xyz = rdots[ ii ]->rotamer()->xyz( iia );

				for ( Size jja = 1; jja <= pose.residue_type( jj ).nheavyatoms(); ++jja ) {

					if ( rdots[ jj ]->get_atom_sasa( jja ) == 0 )
						continue;

					std::string const & jja_elem = pose.residue_type( jj ).atom_type( jja ).element();
					if ( jja_elem != carbon_atom && jja_elem != sulfur_atom )
						continue;

					Real const jja_rad = rdots[ jj ]->get_atom_radius( jja ) + RotamerDots::probe_radius_;
					Vector const & jja_xyz = rdots[ jj ]->rotamer()->xyz( jja );

					if ( iia_xyz.distance_squared( jja_xyz ) < (iia_rad+jja_rad) * (iia_rad+jja_rad) ) {
#ifdef FILE_DEBUG
						/*conformation::Residue const & ii_rsd = pose.residue( ii );
						conformation::Residue const & jj_rsd = pose.residue( jj );
						TR << "compute_pose_hpatch_score(): overlapping atom pair: "
							<< ii_rsd.aa() << " " << ii << "/" << utility::trim( ii_rsd.atom_name( iia ) ) << " - "
							<< jj_rsd.aa() << " " << jj << "/" << utility::trim( jj_rsd.atom_name( jja ) ) << std::endl;*/
#endif
						if ( invdots[ ii ]->atom_overlap_is_exposed( iia, *invdots[ jj ], jja ) ) {

							Size iia_dsid = atom_2_ds_index( AtomID( iia, ii ) );
							Size jja_dsid = atom_2_ds_index( AtomID( jja, jj ) );

							// if the two atom aren't in the same connected component, call union on them
							if ( ds.ds_find( iia_dsid ) != ds.ds_find( jja_dsid ) ) {
								ds.ds_union( iia_dsid, jja_dsid );
							}
						}
#ifdef FILE_DEBUG
						/*else {
						 conformation::Residue const & ii_rsd = pose.residue( ii );
						 conformation::Residue const & jj_rsd = pose.residue( jj );
						 TR << "compute_pose_hpatch_score(): overlapping atom pair, but overlap is buried: "
						 << ii_rsd.aa() << " " << ii << "/" << utility::trim( ii_rsd.atom_name( iia ) ) << " - " 
						 << jj_rsd.aa() << " " << jj << "/" << utility::trim( jj_rsd.atom_name( jja ) ) << std::endl;
						 }*/
#endif
						
					}
				}
			}
		} // end for loop over upper neighbors of ii

	} // end for loop over all residues

	// resize the per-patch score map
	patch_scores_.clear();
	total_hpatch_score_ = 0.0; // just in case the value was not init'd

	// tally up the patch area of all the connected components
	std::map< Size, utility::vector1< Size > > sets = ds.sets();
	std::map< Size, utility::vector1< Size > >::iterator it;

	core::Real patch_area = 0.0;
	for ( it = sets.begin() ; it != sets.end(); it++ ) {
		std::ostringstream strstream;
		patch_area = 0.0;

		// only score patches with 4 or more atoms in them. without this filter, the small patches worsen the ability of
		// the total score to discriminate rosetta decoys from natives. this filter can still cause some really small
		// patches to get scores, but there aren't that many of those.
		if ( (*it).second.size() < 4 ) {
			continue;
		}

#ifdef FILE_DEBUG
		TR << "representative: " << (*it).first << " => atoms in CC: [ ";
#endif
		utility::vector1< id::AtomID > atoms( (*it).second.size() );
		for ( Size ii=1; ii <= (*it).second.size(); ++ii ) {
			id::AtomID const & atomid = ds_index_2_atomid[ (*it).second[ii] ];
			atoms[ ii ] = atomid;
#ifdef FILE_DEBUG
			core::conformation::Residue const & rsd = pose.residue( atomid.rsd() );
			if ( pose.pdb_info() ) {
				if ( pose.pdb_info()->chain( atomid.rsd() ) != ' ' ) {
					TR << pose.pdb_info()->chain( atomid.rsd() ) << "/";
				}
				TR << pose.pdb_info()->number( atomid.rsd() ) << "/" << utility::trim( rsd.atom_name( atomid.atomno() ) ) << " + ";
			} else {
				TR << rsd.seqpos() << "/" << utility::trim( rsd.atom_name( atomid.atomno() ) ) << " + ";
			}
#endif
			patch_area += rdots[ atomid.rsd() ]->get_atom_sasa( atomid.atomno() );
		}
		atoms_in_patches_[ (*it).first ] = atoms;
		atoms.clear();

#ifdef FILE_DEBUG
		TR << "], patch_area: " << patch_area;
#endif

		Real score = 0.0;
		if ( patch_area > MAX_HPATCH_AREA ) {
			score = MAX_HPATCH_SCORE;
		} else {
			score = hpatch_score( patch_area );
		}

#ifdef FILE_DEBUG
		TR << ", hpatch_score: " << score << std::endl;
#endif
		total_hpatch_score_ += score;

		patch_scores_[ (*it).first ] = std::make_pair( score, patch_area );
	}

//#ifdef FILE_DEBUG
	// output all of the scores on a single line
	TR << "calculated total hpatch score: " << total_hpatch_score_ << ", individual patch scores: [ ";
	std::map< Size, std::pair< Real, Real > >::iterator scores_iter;
	for ( scores_iter = patch_scores_.begin(); scores_iter != patch_scores_.end(); scores_iter++ ) {
		TR << (*scores_iter).second.first << ", ";
	}
	TR << "]" << std::endl;
//#endif

	// iterate over the connected components again, but this time output only patches greater than or equal to 250A2
	// by using a score cutoff.
	for ( it = sets.begin() ; it != sets.end(); it++ ) {
		std::ostringstream strstream;
		
		std::map< core::Size, std::pair< core::Real, core::Real > >::iterator ps_it = patch_scores_.find( (*it).first );
		if ( ps_it == patch_scores_.end() ) { continue; } // this shouldn't happen though
		Real score = (*ps_it).second.first;
		if ( score < 4.00 ) { continue; }
		
		TR << "large patch, hpatch_score: " << score << ", PyMOL expression: select p" << (*it).first << ", ";
		for ( Size ii=1; ii <= (*it).second.size(); ++ii ) {
			id::AtomID const & atomid = ds_index_2_atomid[ (*it).second[ii] ];
			
			// output PDB numbering if there's a PDBInfo object present
			core::conformation::Residue const & rsd = pose.residue( atomid.rsd() );
			if ( pose.pdb_info() ) {
				if ( pose.pdb_info()->chain( atomid.rsd() ) != ' ' ) {
					TR << pose.pdb_info()->chain( atomid.rsd() ) << "/";
				}
				TR << pose.pdb_info()->number( atomid.rsd() ) << "/" << utility::trim( rsd.atom_name( atomid.atomno() ) ) << " + ";
			} else {
				TR << rsd.seqpos() << "/" << utility::trim( rsd.atom_name( atomid.atomno() ) ) << " + ";
			}
		}
		
		TR << std::endl;
	}
	
	return;

}

} // namespace interaction_graph
} // namespace pack
} // namespace core

