// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairFactory.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

// Package Headers
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPair1B.hh>
#include <core/scoring/etable/count_pair/CountPairIntraRes.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>

#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>
#include <core/scoring/etable/count_pair/CountPairCrossover4.hh>


// Project Headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/PseudoBond.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

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
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
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
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/types.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/etable/etrie/EtableTrie.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
// AUTO-REMOVED #include <utility/file/FileName.fwd.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/PathName.fwd.hh>
// AUTO-REMOVED #include <utility/file/PathName.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.hh>
// AUTO-REMOVED #include <utility/keys/Key.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
// AUTO-REMOVED #include <utility/keys/KeyLess.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.fwd.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.fwd.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.hh>
// AUTO-REMOVED #include <utility/keys/UserKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.hh>
// AUTO-REMOVED #include <utility/options/FileOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileOption.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.hh>
// AUTO-REMOVED #include <utility/options/Option.fwd.hh>
// AUTO-REMOVED #include <utility/options/Option.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.fwd.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.hh>
// AUTO-REMOVED #include <utility/options/PathOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathOption.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.hh>
// AUTO-REMOVED #include <utility/options/RealOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealOption.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.hh>
// AUTO-REMOVED #include <utility/options/StringOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringOption.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.hh>
// AUTO-REMOVED #include <utility/options/mpi_stderr.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/CArray.fwd.hh>
#include <ObjexxFCL/CArrayP.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
// AUTO-REMOVED #include <ObjexxFCL/TypeTraits.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
// AUTO-REMOVED #include <list>
#include <map>
#include <ostream>
// AUTO-REMOVED #include <set>
#include <sstream>
#include <string>
// AUTO-REMOVED #include <utility>
#include <vector>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/option.hh>

//Auto Headers



namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

CountPairFunctionOP
CountPairFactory::create_count_pair_function(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	CPCrossoverBehavior crossover
)
{
	CPResidueConnectionType connection = determine_residue_connection( res1, res2 );
	CountPairFunctionOP cpfxn;

	if ( connection != CP_NO_BONDS )
	{
		switch ( connection ) {
			case CP_NO_BONDS :
				utility_exit();
			break;
			case CP_ONE_BOND: {
				// scope for res1connat, res2connat initializations
				assert( res1.connections_to_residue( res2 ).size() == 1 );
				assert( res2.connections_to_residue( res1 ).size() == 1 );

				Size res1connat = res1.residue_connection( res1.connections_to_residue( res2 )[ 1 ] ).atomno();
				Size res2connat = res2.residue_connection( res2.connections_to_residue( res1 )[ 1 ] ).atomno();
				switch ( crossover ) {
					case CP_CROSSOVER_3 :
						cpfxn = new CountPair1B< CountPairCrossover3 >( res1, res1connat, res2, res2connat );
					break;
					case CP_CROSSOVER_4 :
						cpfxn = new CountPair1B< CountPairCrossover4 >( res1, res1connat, res2, res2connat );
					break;
				}
			}
			break;
			default: {
				CountPairGenericOP gcpfxn = new CountPairGeneric( res1, res2 );
				gcpfxn->set_crossover( crossover == CP_CROSSOVER_3 ? 3 : 4 );
				cpfxn = gcpfxn;
			}
			break;
		}
	} else {
		cpfxn = new CountPairAll;
	}
	return cpfxn;
}


/// @details Instantiate a count-pair function on the *stack* and rely on an
/// invoker to call a function of the newly created count-pair function.  The
/// count-pair function is removed from the stack without new or delete having
/// been called.
void
CountPairFactory::create_count_pair_function_and_invoke(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	CPCrossoverBehavior crossover,
	Invoker & invoker
)
{
	CPResidueConnectionType connection = determine_residue_connection( res1, res2 );

	if ( connection != CP_NO_BONDS )
	{
		switch ( connection ) {
			case CP_NO_BONDS :
				utility_exit();
			break;
			case CP_ONE_BOND: {
				// scope for res1connat, res2connat initializations
				assert( res1.connections_to_residue( res2 ).size() == 1 );
				assert( res2.connections_to_residue( res1 ).size() == 1 );

				Size res1connat = res1.residue_connection( res1.connections_to_residue( res2 )[ 1 ] ).atomno();
				Size res2connat = res2.residue_connection( res2.connections_to_residue( res1 )[ 1 ] ).atomno();
				switch ( crossover ) {
					case CP_CROSSOVER_3 : {
						CountPair1B< CountPairCrossover3 > cpfxn( res1, res1connat, res2, res2connat );
						invoker.invoke( cpfxn );
					}
					break;
					case CP_CROSSOVER_4 : {
						CountPair1B< CountPairCrossover4 > cpfxn( res1, res1connat, res2, res2connat );
						invoker.invoke( cpfxn );
					}
					break;
				}
			}
			break;
			default: {
				CountPairGeneric gcpfxn( res1, res2 );
				gcpfxn.set_crossover( crossover == CP_CROSSOVER_3 ? 3 : 4 );
				invoker.invoke( gcpfxn );
			}
			break;
		}
	} else {
		CountPairAll cpfxn;
		invoker.invoke( cpfxn );
	}
}

CPResidueConnectionType
CountPairFactory::determine_residue_connection(
	conformation::Residue const & res1,
	conformation::Residue const & res2
)
{
	if ( res1.is_pseudo_bonded( res2.seqpos() )) {
		return CP_MULTIPLE_BONDS_OR_PSEUDOBONDS;
	} else if ( res1.is_bonded(res2) ) {
		if ( res1.connections_to_residue( res2 ).size() == 1 ) {
			return CP_ONE_BOND;
		} else {
			return CP_MULTIPLE_BONDS_OR_PSEUDOBONDS;
		}
	}

	return CP_NO_BONDS;
}


CountPairFunctionOP
CountPairFactory::create_intrares_count_pair_function(
	conformation::Residue const & res,
	CPCrossoverBehavior crossover
)
{
	CountPairFunctionOP cpfxn;

	switch ( crossover ) {
		case CP_CROSSOVER_3 :
			cpfxn = new CountPairIntraRes< CountPairCrossover3 >( res );
		break;
		case CP_CROSSOVER_4 :
			cpfxn = new CountPairIntraRes< CountPairCrossover4 >( res );
		break;
	}

	return cpfxn;

}

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core


