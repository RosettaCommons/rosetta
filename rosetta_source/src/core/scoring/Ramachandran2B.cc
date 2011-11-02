// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran2B.cc
/// @brief  Neighbor Dependent Ramachandran potential class implementation
/// @author Guoli Wang

// Unit Headers
#include <core/scoring/Ramachandran2B.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
// Auto-header: duplicate removed #include <basic/options/option.hh>

#include <platform/types.hh>
#include <core/types.hh>
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
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
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
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
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
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
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
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
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
#include <ObjexxFCL/FArray2A.hh>
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
#include <execinfo.h>
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
#include <sstream>
#include <string>
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
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/pool/poolfwd.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end



using namespace ObjexxFCL;

namespace core {
namespace scoring {

basic::Tracer T("core.scoring.Ramachandran2B");

Real const Ramachandran2B::binw_( 10.0 );

Ramachandran2B::Ramachandran2B() :
	ram_energ_( n_phi_, n_psi_, n_aa_, 0.0 ),
	ram_entropy_( n_aa_, 0.0 ),
	ram_energ_left_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_left_( n_aa_, n_aa_, 0.0 ),
	ram_energ_right_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_right_( n_aa_, n_aa_, 0.0 ),
	rama_score_limit_( 20 ) // customizable in the future, possibly by command line flags.
{
	read_rama();
}


///////////////////////////////////////////////////////////////////////////////

/// @brief evaluate rama score for each (protein) residue and store that score
/// in the pose.energies() object
void
Ramachandran2B::eval_rama_score_all(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const
{
	if ( scorefxn.has_zero_weight( rama ) ) return; // unnecessary, righ?

	//double rama_sum = 0.0;

	// in pose mode, we use fold_tree.cutpoint info to exclude terminus
	// residues from rama calculation. A cutpoint could be either an artificial
	// cutpoint such as loop cutpoint or a real physical chain break such as
	// multiple-chain complex. For the artificial cutpoint, we may need to
	// calculate rama scores for cutpoint residues, but for the real chain break
	// cutpoint, we don't want to do that. So here we first loop over all the
	// residue in the protein and exclude those ones which are the cutpoints.
	// Then we loop over the cutpoint residues and add rama score for residues
	// at artificial cutpoints, i.e., cut_weight != 0.0, which means that
	// jmp_chainbreak_score is also calculated for this cutpoint. Note that the
	// default value for cut_weight here is dependent on whether
	// jmp_chainbreak_weight is set. This is to ensure that rama score for
	// termini residues are not calculated when jmp_chainbreak_weight is 0.0,
	// e.g normal pose docking.

	int const total_residue = pose.total_residue();

	// retrieve cutpoint info // apl do we actually need this data?
	// if so, Pose must provide it 'cause we're offing all global data
	//
	//kinematics::FoldTree const & fold_tree(
	//		pose.fold_tree() );
	//int const n_cut( fold_tree.num_cutpoint() );

	//FArray1D< Real > cut_weight( n_cut,
	//	scorefxns::jmp_chainbreak_weight == 0.0 ? 0.0 : 1.0 ); // apl need to handle

	//if( cut_weight.size1() == scorefxns::cut_weight.size1() )
	//	cut_weight = scorefxns::cut_weight;

	// exclude chain breaks

	Energies & pose_energies( pose.energies() );

	for ( int ii = 1; ii <= total_residue; ++ii )
	{
		if ( pose.residue(ii).is_protein()  && ! pose.residue(ii).is_terminus()  )
		{
			Real rama_score,dphi,dpsi;
			eval_rama_score_residue(pose.residue(ii), pose.residue(ii-1).aa(), pose.residue(ii+1).aa(), rama_score, dphi, dpsi);
			T << "Rama:eval_all: residue " << ii << " " << pose.residue(ii).name() <<
				" " << ii-1 << " " << pose.residue(ii-1).name() << " " << ii+1 << " " <<
				pose.residue(ii+1).name() << " = " << rama_score << std::endl;
			pose_energies.onebody_energies( ii )[rama] = rama_score;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::write_rama_score_all( Pose const & /*pose*/ ) const
{}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const & rsd,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	//assert( pose.residue(res).is_protein() );
	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 || rsd.is_terminus() ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	eval_rama_score_residue( rsd.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
}

///////////////////////////////////////////////////////////////////////////////
// modified by GL
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const &center,
	chemical::AA const left_aa,
	chemical::AA const right_aa,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	assert( center.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	if( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] ) {
		eval_rama_score_residue( center.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
	} else {
		//if(left.seqpos() == center.seqpos()) {
		//	Ramachandran::RamaE_Upper(center, right_aa, drama_dphi, drama_dpsi);
		//} else if(right.seqpos() == center.seqpos()) {
		//	Ramachandran::RamaE_Lower(center, left_aa, drama_dphi, drama_dpsi);
		//} else {
		Real rama_L(0.0), drama_dphi_L(0.0), drama_dpsi_L(0.0);
		Real rama_R(0.0), drama_dphi_R(0.0), drama_dpsi_R(0.0);
		Real rama_0(0.0), drama_dphi_0(0.0), drama_dpsi_0(0.0);
		rama_L = RamaE_Lower(center, left_aa, drama_dphi_L, drama_dpsi_L);
		rama_R = RamaE_Upper(center, right_aa, drama_dphi_R, drama_dpsi_R);
		rama_0 = RamaE(center, drama_dphi_0, drama_dpsi_0);

		rama = rama_L + rama_R - rama_0;
		drama_dphi = drama_dphi_L + drama_dphi_R - drama_dphi_0;
		drama_dpsi = drama_dpsi_L + drama_dpsi_R - drama_dpsi_0;
		//}
	}
}

void
Ramachandran2B::IdealizeRamaEnergy(
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi,
	Real const entropy,
	FArray2A< Real > const & rama_for_res
) const
{
	using namespace numeric::interpolation::periodic_range::half;
	Real interp_E = bilinearly_interpolated( phi, psi, binw_, n_phi_, rama_for_res, drama_dphi, drama_dpsi );
	// rama = IdealizeRamaEnergy(ram_entropy_(center.aa(), leftIndex, rightIndex), interp_E, drama_dphi, drama_dpsi);
	rama = entropy + interp_E;
	//	std::cout << "Rama::eval_res: " <<  interp_E << " rama " << rama << std::endl;

	if ( ! basic::options::option[basic::options::OptionKeys::corrections::score::rama_not_squared] ) {
		if ( rama > 1.0 ) {
			Real rama_squared = rama * rama;
			if ( rama_squared > rama_score_limit_ ) {
				drama_dphi = 0.0;
				drama_dpsi = 0.0;
				rama = rama_score_limit_;
			} else {
				drama_dphi *= 2.0 * rama;
				drama_dpsi *= 2.0 * rama;
				rama = rama_squared;
			}
		}
	}

	// std::cout << " rama: " << rama << " dphi " << drama_dphi << " dpsi " << drama_dpsi << std::endl;
}

// end modification

///////////////////////////////////////////////////////////////////////////////
// modified by GL according to Andrew's suggestion
Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const
{
	Real drama_dphi, drama_dpsi;
	return RamaE_Lower( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if( ! option[ score::ramaneighbors ] ) {
		return 0.0;
	}

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_left_(1, 1, rsd.aa(), neighbor), zero_index, zero_index );
	Real entropy = ram_entropy_left_(rsd.aa(), neighbor);

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}
Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const
{
	Real drama_dphi, drama_dpsi;
	return RamaE_Upper( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] )
	{
		return 0.0;
	}

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_right_(1, 1, rsd.aa(), neighbor), zero_index, zero_index );
	Real entropy = ram_entropy_right_(rsd.aa(), neighbor);

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd
) const
{
	Real drama_dphi(0.0), drama_dpsi(0.0);
	return RamaE( rsd, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{

	using namespace numeric;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	Real ramaE(0.0);
	eval_rama_score_residue( rsd.aa(), phi, psi, ramaE, drama_dphi, drama_dpsi );
	return ramaE;
}

///////////////////////////////////////////////////////////////////////////////
///
Real
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi
) const
{

	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( res_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	return rama;
}

///////////////////////////////////////////////////////////////////////////////
///
void
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_(1, 1, res_aa ), zero_index, zero_index );
	Real entropy = ram_entropy_(res_aa);

	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
}


// Guoli Wang
void
Ramachandran2B::read_rama()
{
	using namespace basic::options;

	using namespace basic::options::OptionKeys;

	int aa_num( 0 ), aa_num_left( 0 ), aa_num_right( 0 );
	int phi_bin( 0 ), psi_bin( 0 ), ss_type( 0 );
	int tCounts( 0 );
	Real tProb( 0.0 ), tEnergy( 0.0 );

	Size line_count( 0 );


	//utility::io::izstream  iunit;
#ifndef WIN32
#ifndef __CYGWIN__
	clock_t starttime = clock();
#endif
#endif
	std::string energyFileName = basic::options::option[ in::file::rama2b_map ]().name() ; // "wrapHDPprobs36.both";
	T << "Read in ramachandran map: " <<  energyFileName << std::endl;
	utility::io::izstream iRamaEnergy;
	basic::database::open( iRamaEnergy, energyFileName );
	while( ! iRamaEnergy.eof() ) {
		++line_count;
		iRamaEnergy >> aa_num >> aa_num_left >> aa_num_right >> ss_type >> phi_bin >> psi_bin >> tCounts >> tProb >> tEnergy;
		// std::cout << " aa_num " << aa_num << " aa_num_left " << aa_num_left << " aa_num_right " << aa_num_right << " ss_type " << ss_type <<
		// 			" phi_bin " << phi_bin << " psi_bin " << psi_bin << " tProb " << tProb << " tEnergy " << tEnergy << std::endl;
		if(aa_num > n_aa_) continue;

		int phiIndex = phi_bin / 10 + 1;
		int psiIndex = psi_bin / 10 + 1;
		Real entropy = -1.0 * tProb * tEnergy;

		if( aa_num_left == nullaa && aa_num_right == nullaa ) {
			ram_energ_( phiIndex, psiIndex, aa_num ) = tEnergy;
			ram_entropy_( aa_num ) += entropy;
		} else if( aa_num_left != nullaa ) {
			ram_energ_left_( phiIndex, psiIndex, aa_num, aa_num_left ) = tEnergy;
			ram_entropy_left_( aa_num, aa_num_left ) += entropy;
		} else if( aa_num_right != nullaa ) {
			ram_energ_right_( phiIndex, psiIndex, aa_num, aa_num_right ) = tEnergy;
			ram_entropy_right_( aa_num, aa_num_right ) += entropy;
		}
	}

	iRamaEnergy.close();
#ifndef WIN32
#ifndef __CYGWIN__
	clock_t stoptime = clock();
	T << "Reading Rama from database took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << " seconds" << std::endl;
#endif
#endif
}


} // scoring
} // core
