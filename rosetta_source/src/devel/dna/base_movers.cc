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


#include <devel/dna/base_movers.hh>
#include <devel/dna/relax_util.hh>

#include <devel/cartesian_frags/DNA_FragLib.hh>
#include <devel/cartesian_frags/dna_util.hh>

#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/conformation/Residue.hh>

#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>

#include <numeric/random/random.hh>

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
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/find_neighbors.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
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
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
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
#include <core/scoring/dna/BasePartner.fwd.hh>
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
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <devel/cartesian_frags/CartesianFragment.fwd.hh>
#include <devel/cartesian_frags/CartesianFragment.hh>
#include <devel/cartesian_frags/DNA_FragLib.fwd.hh>
#include <devel/cartesian_frags/Direction.hh>
#include <devel/cartesian_frags/SafeID.fwd.hh>
#include <devel/cartesian_frags/SafeID.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
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
#include <utility/tag/Tag.fwd.hh>
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
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzTriple.hh>
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
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
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
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_map.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


namespace devel {
namespace dna {

using namespace core;
using namespace devel::cartesian_frags;
using utility::vector1;

static basic::Tracer tt( "devel.dna.base_movers", basic::t_trace );
static basic::Tracer td( "devel.dna.base_movers", basic::t_debug );
static basic::Tracer ti( "devel.dna.base_movers", basic::t_info );
static basic::Tracer tw( "devel.dna.base_movers", basic::t_warning );

static numeric::random::RandomGenerator RG(53423); // <- Magic number, do not change it!!!

/// @details  Setup a pointgraph for later use in dma calcs
conformation::PointGraphOP
setup_dme_point_graph( pose::Pose const & ref_pose, Real const threshold )
{
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( ref_pose.conformation(), *pg );
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, threshold );
	return pg;
}

/// @details  Calculate the dme using a pointgraph
Real
point_graph_dme( conformation::PointGraph const & pg, pose::Pose const & pose )
{
	Size total(0);
	Real dme(0.0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & i_rsd( pose.residue(i) );
		for ( conformation::PointGraph::UpperEdgeListConstIter
						i_iter     = pg.get_vertex( i ).const_upper_edge_list_begin(),
						i_end_iter = pg.get_vertex( i ).const_upper_edge_list_end();
					i_iter != i_end_iter; ++i_iter ) {
			Size const j = i_iter->upper_vertex();
			Real const reference_distance( std::sqrt( i_iter->data().dsq() ) );
			Real const pose_distance( i_rsd.nbr_atom_xyz().distance( pose.residue(j).nbr_atom_xyz() ) );
			dme += ( reference_distance - pose_distance ) * ( reference_distance - pose_distance );
			++total;
		}
	}
	tt << "dme_nbrs: " << total << std::endl;
	return std::sqrt( dme / total );
}

///////////////////////////////////////////////////////////////////////////////
///
/// @details  Make a base-pair fragment insertion between residue seqpos and seqpos-partner
/// try to patch up the backbone after making the insertion

void
make_base_pair_move(
	Size const seqpos,
	DNA_FragLib const & lib,
	Real const frag_dev_threshold,
	Size const max_tries,
	Real const max_score_increase,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose & pose_inout
)
{
	PROF_START( basic::MAKE_BASE_PAIR_MOVE );

	//// scorefxn is most likely a chainbreak scoring function.
	using namespace conformation;
	using namespace kinematics;
	using namespace scoring::dna;
	using namespace id;

	pose::Pose const start_pose( pose_inout );
	pose::Pose pose( pose_inout ); // our local pose

	// test the use of rsd-rsd distances to gauge the impact of a fragment insertion
	conformation::PointGraphOP pg( setup_dme_point_graph( pose, 10.0 ) );

	Real const start_score( scorefxn( pose ) );

	Size const seqpos_partner( retrieve_base_partner_from_pose( pose )[ seqpos ] );
	chemical::ResidueType const &  rsd_type( pose.residue_type( seqpos         ) );
	chemical::ResidueType const & prsd_type( pose.residue_type( seqpos_partner ) );
	std::string const bp( std::string() + rsd_type.name1() + prsd_type.name1() );
	vector1< CartesianFragment > const & bps( lib.base_pairs(bp) );
	vector1< std::pair< Real, Size > > frag_devs;

	Size top_nn( 0 );
	{ // compute deviations from current for all frags
		Stub const chi  ( torsion_stub( TorsionID( seqpos        , CHI, 1 ), Backward, pose.conformation() ) );
		Stub const chi_p( torsion_stub( TorsionID( seqpos_partner, CHI, 1 ), Backward, pose.conformation() ) );
		RT const current( RT( chi, chi_p ) );
		for ( Size i=1; i<= bps.size(); ++i ) {
			frag_devs.push_back( std::make_pair( current.distance_squared( bps[i].rt(1) ), i ) );
		}
		std::sort( frag_devs.begin(), frag_devs.end() );
		td << "bp_move: frag_devs: " << frag_devs[1].first << ' ' << frag_devs[2].first << std::endl;
		top_nn = std::min( Size(5), frag_devs.size() );
		while ( top_nn < frag_devs.size() && frag_devs[ top_nn ].first < frag_dev_threshold ) ++top_nn;
	}

	// need these for making the fragment insertion
	vector1< Size > offsets;
	offsets.push_back( seqpos );
	offsets.push_back( seqpos_partner );

	Size ntries( 0 );
	Real best_score( 999.9 ), dme(0.0), rmsd(0.0), frag_dev(0.0);
	while ( ntries < max_tries ) { // keep looping until we like the move
		++ntries;

		// make a random basepair fragment insertion
		Size const nn( static_cast< int >( top_nn * RG.uniform() ) + 1 );
		bps[ frag_devs[ nn ].second ].insert( pose.conformation(), offsets );

		//
		dme = point_graph_dme( *pg, pose );
		rmsd = scoring::nbr_atom_rmsd( pose, start_pose );
		frag_dev = frag_devs[ nn ].first;

		tt << "bp_move: choose_frag: nn= " << nn << " top_nn= " <<top_nn << " max_nn= " << bps.size() <<
			" max_frag_dev= " << frag_dev_threshold << " frag_dev= " << frag_dev <<
			" dme= " << dme << " rmsd= " << rmsd << std::endl;

		// by inserting this fragment we have perturbed some of the backbone segments


		// did we perturb the seqpos-1 --> seqpos   connection?
		if ( !rsd_type.is_lower_terminus() && !pose.fold_tree().jump_exists( seqpos-1, seqpos ) ) {
			patch_up_backbone_link( seqpos-1, lib, scorefxn, pose );
		}
		// did we perturb the seqpos   --> seqpos+1 connection?
		if ( !rsd_type.is_upper_terminus() && !pose.fold_tree().jump_exists( seqpos, seqpos+1 ) ) {
			patch_up_backbone_link( seqpos  , lib, scorefxn, pose );
		}

		if ( !prsd_type.is_lower_terminus() && !pose.fold_tree().jump_exists( seqpos_partner-1, seqpos_partner ) ) {
			patch_up_backbone_link( seqpos_partner-1, lib, scorefxn, pose );
		}
		if ( !prsd_type.is_upper_terminus() && !pose.fold_tree().jump_exists( seqpos_partner, seqpos_partner+1 ) ) {
			patch_up_backbone_link( seqpos_partner  , lib, scorefxn, pose );
		}

		Real const final_score( scorefxn( pose ) );
		tt << "bp_move: ntries= " << ntries << " score_increase= " << final_score - start_score <<
			" max_score_increase= " << max_score_increase << " dme,rmsd " << dme << ' ' << rmsd << std::endl;
		scorefxn.show( tt, pose );

		if ( ntries == 1 || final_score < best_score ) {
			best_score = final_score;
			pose_inout = pose;
		}

		if ( final_score - start_score < max_score_increase ) break;

		pose = start_pose; // unmake the move
	}
	td << "make_base_pair_move: score_increase= " << best_score - start_score << " dme= " << dme << " rmsd= " <<
		rmsd << " frag_dev= " << frag_dev << " ntries= " << ntries << std::endl;
	PROF_STOP( basic::MAKE_BASE_PAIR_MOVE );

}


///////////////////////////////////////////////////////////////////////////////
///
/// @details  Make a base-step fragment insertion between residue seqpos and seqpos+1
/// try to patch up the backbone after making the insertion

void
make_base_step_move(
	Size const seqpos,
	DNA_FragLib const & lib,
	Real const frag_dev_threshold,
	Size const max_tries,
	Real const max_score_increase,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose & pose_inout
)
{
	//// scorefxn is most likely a chainbreak scoring function.
	PROF_START( basic::MAKE_BASE_STEP_MOVE );

	using namespace conformation;
	using namespace kinematics;
	using namespace scoring::dna;
	using namespace id;

	pose::Pose const start_pose( pose_inout );
	pose::Pose pose( pose_inout ); // our local pose
	assert( pose.fold_tree().jump_exists( seqpos, seqpos+1 ) );

	// test the use of rsd-rsd distances to gauge the impact of a fragment insertion
	conformation::PointGraphOP pg( setup_dme_point_graph( pose, 10.0 ) );

	Real const start_score( scorefxn( pose ) );

	Size const seqpos_partner( retrieve_base_partner_from_pose( pose )[ seqpos ] );
	chemical::ResidueType const &      rsd_type( pose.residue_type( seqpos   ) );
	chemical::ResidueType const & next_rsd_type( pose.residue_type( seqpos+1 ) );
	//Conformation & conf( pose.conformation() );

	std::string const bs( std::string() + rsd_type.name1() + next_rsd_type.name1() );
	vector1< CartesianFragment > const & bss( lib.base_steps( bs ) );
	vector1< std::pair< Real, Size > > frag_devs;

	Size top_nn( 0 );
	{ // compute deviations from current for all frags
		Stub const chi  ( torsion_stub( TorsionID( seqpos  , CHI, 1 ), Backward, pose.conformation() ) );
		Stub const chi_p( torsion_stub( TorsionID( seqpos+1, CHI, 1 ), Backward, pose.conformation() ) );
		RT const current( RT( chi, chi_p ) );
		for ( Size i=1; i<= bss.size(); ++i ) {
			frag_devs.push_back( std::make_pair( current.distance_squared( bss[i].rt(1) ), i ) );
		}
		std::sort( frag_devs.begin(), frag_devs.end() );
		tt << "bs_move: frag_devs: " << frag_devs[1].first << ' ' << frag_devs[2].first << std::endl;
		top_nn = std::min( Size(5), frag_devs.size() );
		while ( top_nn < frag_devs.size() && frag_devs[ top_nn ].first < frag_dev_threshold ) ++top_nn;
	}

	// need these for making the fragment insertion
	vector1< Size > offsets;
	offsets.push_back( seqpos   );
	offsets.push_back( seqpos+1 );

	Size ntries( 0 );
	Real best_score( 999.9 ), dme(0.0), rmsd(0.0), frag_dev(0.0);

	while ( ntries < max_tries ) { // keep looping until we like the move
		++ntries;

		// make a random basepair fragment insertion
		Size const nn( static_cast< int >( top_nn * RG.uniform() ) + 1 );
		bss[ frag_devs[ nn ].second ].insert( pose.conformation(), offsets );

		//
		dme = point_graph_dme( *pg, pose );
		rmsd = scoring::nbr_atom_rmsd( pose, start_pose );
		frag_dev = frag_devs[ nn ].first;

		tt << "bs_move: choose_frag: nn= " << nn << " top_nn= " <<top_nn << " max_nn= " << bss.size() <<
			" max_frag_dev= " << frag_dev_threshold << " frag_dev= " << frag_dev <<
			" dme= " << dme << " rmsd= " << rmsd << std::endl;

		// by inserting this fragment we have perturbed two of the backbone segments
		patch_up_backbone_link( seqpos, lib, scorefxn, pose );

		patch_up_backbone_link( seqpos_partner-1, lib, scorefxn, pose );

		Real const final_score( scorefxn( pose ) );
		tt << "bs_move: ntries= " << ntries << " score_increase= " << final_score - start_score <<
			" max_score_increase= " << max_score_increase << std::endl;

		scorefxn.show( tt, pose );

		if ( ntries == 1 || final_score < best_score ) {
			best_score = final_score;
			pose_inout = pose;
		}

		if ( final_score - start_score < max_score_increase ) break;

		pose = start_pose; // unmake the move
	}
	td << "make_base_step_move: score_increase= " << best_score - start_score << " dme= " << dme << " rmsd= " <<
		rmsd << " frag_dev= " << frag_dev << " ntries= " << ntries << std::endl;

	PROF_STOP ( basic::MAKE_BASE_STEP_MOVE );
}


/// @details  Make a basepair move using make_base_pair_move
void
BasePairMover::apply( pose::Pose & pose )
{
	make_base_pair_move(
		choose_random_base_pair( pose ), *lib_, frag_dev_threshold_, max_tries_,
		max_score_increase_, *scorefxn_, pose );
}

std::string
BasePairMover::get_name() const {
	return "BasePairMover";
}


/// @details  Make a basestep move using make_base_step_move
void
BaseStepMover::apply( core::pose::Pose & pose )
{
	make_base_step_move(
		choose_random_base_step_jump( pose ), *lib_, frag_dev_threshold_, max_tries_,
		max_score_increase_, *scorefxn_, pose );
}

std::string
BaseStepMover::get_name() const {
	return "BaseStepMover";
}

} // ns dna
} // ns devel
