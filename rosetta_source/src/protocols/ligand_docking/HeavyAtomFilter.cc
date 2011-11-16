// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/ligand_docking/HeavyAtomFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/HeavyAtomFilter.hh>
#include <protocols/ligand_docking/HeavyAtomFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
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
// AUTO-REMOVED #include <core/chemical/types.hh>
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
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
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
// AUTO-REMOVED #include <core/pose/PDBPoseMap.fwd.hh>
// AUTO-REMOVED #include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/Remarks.fwd.hh>
// AUTO-REMOVED #include <core/pose/Remarks.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/FilterCreator.hh>
#include <protocols/moves/DataMap.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.hh>
// AUTO-REMOVED #include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
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
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <numeric/NumericTraits.hh>
// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
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
#include <set>
#include <sstream>
#include <string>
// AUTO-REMOVED #include <typeinfo>
// AUTO-REMOVED #include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto Headers



namespace protocols {
namespace ligand_docking {

static basic::Tracer heavy_atom_tracer( "protocols.ligand_docking.HeavyAtomFilter" );

bool
HeavyAtomFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(heavy_atom_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	if(	core::pose::num_heavy_atoms(start,end,pose) > heavy_atom_limit_ ){
		heavy_atom_tracer<< "Reached heavy atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
HeavyAtomFilter::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "HeavyAtom" ) {
		heavy_atom_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("heavy_atom_limit") ) ){
		utility_exit_with_message("HeavyAtom filter needs a 'chain' and a 'heavy_atom_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	heavy_atom_limit_ = tag->getOption<core::Size>("heavy_atom_limit");

}

protocols::filters::FilterOP
HeavyAtomFilterCreator::create_filter() const { return new HeavyAtomFilter; }

std::string
HeavyAtomFilterCreator::keyname() const { return "HeavyAtom"; }



} // ligand_docking
} // protocols
