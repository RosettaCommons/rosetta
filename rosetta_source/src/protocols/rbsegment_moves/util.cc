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
/// @author Frank DiMaio
/// @author Srivatsan Raman
#include <protocols/rbsegment_moves/RBSegmentMover.hh>
#include <protocols/rbsegment_moves/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/jumping/Dssp.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Random number generator
#include <numeric/xyzVector.io.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray1D.hh>

//
#include <string>

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
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
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
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
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
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AmbiguousConstraint.fwd.hh>
#include <core/scoring/constraints/BoundConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/StrandPairing.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rbsegment_moves/RBSegment.fwd.hh>
#include <protocols/rbsegment_moves/RBSegment.hh>
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
#include <utility/tag/Tag.fwd.hh>
#include <numeric/IOTraits.hh>
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
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
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
#include <functional>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
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
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


namespace protocols {
namespace rbsegment_moves {

using namespace core;

static numeric::random::RandomGenerator rbseg_RG(186331);
static basic::Tracer TR("protocols::moves::RBSegmentMover");


//
///@brief set up constraints from RB segs
///       currently uses bounded constraints on each CA ... make CST type flag-selectable????
void set_rb_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	utility::vector1< protocols::rbsegment_moves::RBSegment > const & rbsegs ,
	core::id::SequenceMapping const & resmap,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth
                   ) {
	core::scoring::constraints::ConstraintSetOP new_csts =	new core::scoring::constraints::ConstraintSet;

	for (int i =  1; i <= (int)rbsegs.size(); ++i) {
		for (core::Size j=1; j<=rbsegs[i].nContinuousSegments(); ++j) {
			int rbstart = rbsegs[i][j].start(), rbend = rbsegs[i][j].end();

			for (int cst_i=rbstart; cst_i<=rbend; ++cst_i) {

				// make ambiguous cst from ( cst_i-cst_seqwidth , cst_i+cst_seqwidth )
				core::scoring::constraints::ConstraintCOPs CSTs_i;

				//
				for (int cst_ij = std::max(1,cst_i-(int)cst_seqwidth), cst_ij_end=std::min(cst_pose.total_residue(),cst_i+cst_seqwidth);
				     cst_ij<=cst_ij_end; ++cst_ij) {
					// only constrain protein residues
					if (!cst_pose.residue(cst_ij).is_protein()) continue;

					core::scoring::constraints::ConstraintOP newcst_ij = new core::scoring::constraints::CoordinateConstraint(
											 core::id::AtomID( pose.residue(resmap[cst_i]).atom_index("CA"), resmap[cst_i]),
											 core::id::AtomID( pose.residue(pose.total_residue()).atom_index("ORIG"), pose.total_residue()),
											 cst_pose.residue(cst_ij).xyz( "CA" ),
											 new core::scoring::constraints::BoundFunc(0,cst_width,cst_stdev,"xyz") );


					// make sure not randomized <<<< kind of hacky
					core::Real len = (cst_pose.residue(cst_ij).xyz( "CA" )).length();
					if ( len < 500 )
						CSTs_i.push_back( newcst_ij );
				}

				//
				if ( CSTs_i.size() > 0 ) {
					TR << "Adding " << CSTs_i.size() << " ambiguous constraints for res " << resmap[cst_i]
					   << " (= input " << cst_i << ")" << std::endl;
					new_csts->add_constraint( new core::scoring::constraints::AmbiguousConstraint( CSTs_i ) );
				}

			}
		}
	}
	pose.constraint_set( new_csts );
}


///@brief set up constraints on complete pose (not just RB segments)
///       currently uses bounded constraints ... make CST type flag-selectable????
void set_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth
                   ) {
	if (pose.total_residue() != cst_pose.total_residue()) {
		TR.Warning << "set_constraints() error: #res in cst_pose (" << cst_pose.total_residue()
		            << ") != #res in pose (" << pose.total_residue() << ").  Continuing..." << std::endl;
	}
	int nres = (int)std::min( pose.total_residue() , cst_pose.total_residue() );

	core::scoring::constraints::ConstraintSetOP new_csts =	new core::scoring::constraints::ConstraintSet;
	for (int cst_i =  1; cst_i < nres; ++cst_i) {
		// make ambiguous cst from ( cst_i-cst_seqwidth , cst_i+cst_seqwidth )
		core::scoring::constraints::ConstraintCOPs CSTs_i;

		for (int cst_ij=std::max(1,cst_i-(int)cst_seqwidth), cst_ij_end=std::min(nres,cst_i+(int)cst_seqwidth);
				 cst_ij<=cst_ij_end; ++cst_ij) {
			// only constrain protein residues
			if (!cst_pose.residue(cst_ij).is_protein()) continue;

			core::scoring::constraints::ConstraintOP newcst_ij = new core::scoring::constraints::CoordinateConstraint(
				core::id::AtomID( pose.residue(cst_i               ).atom_index("CA"  ), cst_i),
				core::id::AtomID( pose.residue(pose.total_residue()).atom_index("ORIG"), pose.total_residue()),
				cst_pose.residue(cst_ij).xyz( "CA" ),
				new core::scoring::constraints::BoundFunc(0,cst_width,cst_stdev,"xyz")
			);

			// make sure not randomized <<<< kind of hacky
			core::Real len = (cst_pose.residue(cst_ij).xyz( "CA" )).length();
			if ( len < 500 )
				CSTs_i.push_back( newcst_ij );
		}

		//
		if ( CSTs_i.size() > 0 ) {
			TR << "Adding " << CSTs_i.size() << " ambiguous constraints for res " << cst_i << std::endl;
			new_csts->add_constraint( new core::scoring::constraints::AmbiguousConstraint( CSTs_i ) );
		}
	}
	pose.constraint_set( new_csts );
}


///
///@brief Helper function to set up a pose; unlike alt version keep loops (use cutpoint variants)
///   unlike version in loops_main, this uses RBSegment structure to build multi-level topology
/// returns jump residues
utility::vector1< core::Size > setup_pose_rbsegs_keep_loops(
              core::pose::Pose &pose,
              utility::vector1< protocols::rbsegment_moves::RBSegment > const &rbsegs,
              loops::Loops const &loops,
              core::kinematics::MoveMapOP mm) {
 	using namespace core::kinematics;

 	core::Size nres( pose.total_residue()-1 );
 	core::Size vrtid = nres+1;
 	core::Size nrbsegs( rbsegs.size() );

 	// "star" topology fold tree
	utility::vector1< core::Size > cuts, jump_res;
	utility::vector1< std::pair<core::Size,core::Size> > jumps;
	for( loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
		if (it->start() == 1 || it->stop() == nres) continue;
		cuts.push_back( it->cut() );
	}
	cuts.push_back( nres );
	//cuts.push_back( nres+1 ); //?
	for (int i=1; i <= (int)nrbsegs; ++i ) {
		// jump from vrt to midpt of 1st seg
		core::Size midpt = (rbsegs[i][1].start()+rbsegs[i][1].end()) / 2;
		jumps.push_back (std::pair<core::Size,core::Size>( vrtid, midpt ) );
		jump_res.push_back ( midpt );
	}
	for (int i=1; i <= (int)nrbsegs; ++i ) {
		for (int j=2; j<= (int)rbsegs[i].nContinuousSegments() ; ++j ) {
			// jump from midpt of 1st seg here
			core::Size midpt1 = (rbsegs[i][1].start()+rbsegs[i][1].end()) / 2;
			core::Size midpt2 = (rbsegs[i][j].start()+rbsegs[i][j].end()) / 2;
			jumps.push_back (std::pair<core::Size,core::Size>( midpt1, midpt2 ) );
			jump_res.push_back ( midpt2 );
		}
	}

	ObjexxFCL::FArray2D_int fjumps( 2, jumps.size() );
	ObjexxFCL::FArray1D_int fcuts ( cuts.size() );
	for ( Size i=1; i<=jumps.size(); ++i ) {
		fjumps(1,i) = std::min( jumps[i].first , jumps[i].second );
		fjumps(2,i) = std::max( jumps[i].first , jumps[i].second );
	}
	for ( Size i=1; i<=cuts.size(); ++i ) {
		fcuts(i) = cuts[i];
	}
	kinematics::FoldTree f;
	bool valid_tree = f.tree_from_jumps_and_cuts( nres+1, jumps.size(), fjumps, fcuts );
	runtime_assert( valid_tree );
	f.reorder( vrtid );
	TR << "New fold tree: " << f << std::endl;

	pose.fold_tree( f );

	// movemap
	mm->clear();
	for( loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it )
		for ( int j=(int)(it->start()) ; j<=(int)(it->stop()); ++j )
			mm->set_bb(j,true);
	for (int i=1; i <= (int)nrbsegs; ++i )
		mm->set_jump(i,true);

	// cb variants
	for( protocols::loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
	 	Size const loop_cut(it->cut());
	 	if ( !pose.residue(loop_cut).is_upper_terminus() ) {
	 		if ( ! pose.residue(loop_cut).has_variant_type(core::chemical::CUTPOINT_LOWER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, loop_cut );
	 		if ( ! pose.residue(loop_cut+1).has_variant_type(core::chemical::CUTPOINT_UPPER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, loop_cut+1 );
	 	}
	}

	return jump_res;
}


///
///@brief Helper function to set up a loop-removed pose
void setup_pose_from_rbsegs(
             utility::vector1< protocols::rbsegment_moves::RBSegment > const &rbsegs ,
             core::pose::Pose const &pose_in ,
             core::pose::Pose &pose_out ,
             core::id::SequenceMapping &resmap,
             core::kinematics::MoveMap &mm,
             bool fix_ligand  ) {
	using namespace core::kinematics;

	core::Size nres( pose_in.total_residue() );
	core::Size nres_rb = 0;

	// each ligand is its own rb-segment (for the purposes of fold-tree generation)
	utility::vector1< protocols::rbsegment_moves::RBSegment > rbsegs_with_ligands = rbsegs;
	core::Size last_peptide_res = nres;
	while ( !pose_in.residue( last_peptide_res ).is_protein() )
		last_peptide_res--;
	for (core::Size i=last_peptide_res+1; i<=nres; ++i) {
		if ( pose_in.residue( i ).aa() != core::chemical::aa_vrt ) {
			rbsegs_with_ligands.push_back( protocols::rbsegment_moves::RBSegment( i, i, 'X' ) );
			TR << "setup_pose_from_rbsegs: Ligand at " << i << std::endl;
		}
	}

	// count rb reses
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i )
		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j)
			nres_rb += rbsegs_with_ligands[i][j].end() - rbsegs_with_ligands[i][j].start() + 1;

	resmap.resize( nres, nres_rb );

	pose_out.clear();
	int rb_ctr = 0;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {

		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {
			core::Size rb_start  = rbsegs_with_ligands[i][j].start() ,
			           rb_end    = rbsegs_with_ligands[i][j].end()   ;
			int nsegment_reses = rb_end - rb_start + 1;
			rb_ctr++;

			char secstruct = rbsegs_with_ligands[i][j].char_type();

			if (nsegment_reses > 0) {
				core::conformation::Residue new_res_0 = pose_in.residue( rb_start );
				pose_out.append_residue_by_jump( new_res_0 , pose_out.total_residue(), "" , "" , true );
				resmap[rb_start] = pose_out.total_residue();

				for (int k=1; k<nsegment_reses; ++k) {
					core::conformation::Residue new_res_k = pose_in.residue( rb_start+k );
					pose_out.append_residue_by_bond( new_res_k );
					resmap[rb_start + k] = pose_out.total_residue();

					// set secstruct
					if ( (secstruct == 'H' || secstruct == 'E') && k!=nsegment_reses-1 )
						pose_out.set_secstruct( pose_out.total_residue(), secstruct );
				}
			}
		}
	}

	// virtual res as root
	core::chemical::ResidueTypeSetCAP const &residue_set(
									core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )  );
	core::chemical::ResidueTypeCAPs const & rsd_type_list( residue_set->name3_map("VRT") );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );
	pose_out.append_residue_by_jump( *new_res , 1 );

	// "star" topology fold tree
	core::kinematics::FoldTree newF;
	core::Size nstart = 1, njump = 1;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {
		core::Size in_start  = rbsegs_with_ligands[i][1].start() ,
				       in_end    = rbsegs_with_ligands[i][1].end();
		core::Size out_start  = nstart,
				       out_end    = nstart + (in_end - in_start);
		core::Size out_midpt  = (out_end + out_start)/2;

		newF.add_edge( nres_rb+1, out_midpt, njump );
		newF.add_edge( out_midpt, out_start,    Edge::PEPTIDE );
		newF.add_edge( out_midpt, out_end  , Edge::PEPTIDE );
		nstart = out_end+1;

		// let all jumps except those to ligands move
		if ( !fix_ligand || i <= rbsegs.size() )
			mm.set_jump( njump , true );

		njump++;

		// in strands let bb move
		if (rbsegs_with_ligands[i].isSheet()) {
			for ( core::Size j=out_start; j<=out_end; ++j) {
				mm.set_bb( j, true );
			}
		}

		core::Size in_midpt = out_midpt;
		for (core::Size j=2; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {
			in_start  = rbsegs_with_ligands[i][j].start();
			in_end    = rbsegs_with_ligands[i][j].end();
			out_start  = nstart;
			out_end    = nstart + (in_end - in_start);
			out_midpt  = (out_end + out_start)/2;

			newF.add_edge( in_midpt, out_midpt, njump );
			newF.add_edge( out_midpt, out_start, Edge::PEPTIDE );
			newF.add_edge( out_midpt, out_end  , Edge::PEPTIDE );
			nstart = out_end+1;
			mm.set_jump( njump , false );
			njump++;
		}
	}

	newF.reorder( nres_rb+1 );
	pose_out.fold_tree( newF );

	TR << "New fold tree: " << newF << std::endl;
}



///@brief Helper function to restore a fully-connected pose
void restore_pose_from_rbsegs(
             utility::vector1< protocols::rbsegment_moves::RBSegment > const &rbsegs ,
             core::pose::Pose const &pose_in ,
             core::pose::Pose &pose_out /* input/output */ )
{
	using namespace core::kinematics;

	// each ligand is its own rb-segment (for the purposes of fold-tree generation)
	utility::vector1< protocols::rbsegment_moves::RBSegment > rbsegs_with_ligands = rbsegs;
	/*
	int nres( pose_out.total_residue() );
	// keep ligand confs from starting pose
	int last_peptide_res = nres;
	while ( !pose_out.residue( last_peptide_res ).is_protein() )
		last_peptide_res--;
	for (core::Size i=last_peptide_res+1; i<=nres; ++i) {
		if ( pose_out.residue( i ).aa() != core::chemical::aa_vrt ) {
			rbsegs_with_ligands.push_back( protocols::rbsegment_moves::RBSegment( i, i, 'X' ) );
			TZ << "restore_pose_from_rbsegs: Ligand at " << i << std::endl;
		}
	}
	*/

	int res_rb_counter = 1;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {

		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {

			core::Size rb_start  = rbsegs_with_ligands[i][j].start() ,
								 rb_end    = rbsegs_with_ligands[i][j].end()   ;
			int nsegment_reses = rb_end - rb_start + 1;

			// xyz copy??
			if (nsegment_reses > 0) {
					pose_out.copy_segment( nsegment_reses, pose_in,  rb_start, res_rb_counter  );
			}
			// copy secstruct
			for (int k=(int)rb_start; k<=(int)rb_end; ++k)
				pose_out.set_secstruct( k, pose_in.secstruct( res_rb_counter+k-rb_start ) );
			res_rb_counter += nsegment_reses;
		}

	}

	// virtual res as root
	core::chemical::ResidueTypeSetCAP const &residue_set(
									core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )  );
	core::chemical::ResidueTypeCAPs const & rsd_type_list( residue_set->name3_map("VRT") );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );
	pose_out.append_residue_by_jump( *new_res , pose_out.total_residue()/2 );

	// make the virt atom the root
	core::kinematics::FoldTree newF(pose_out.fold_tree());
	newF.reorder( pose_out.total_residue() );
	pose_out.fold_tree( newF );

	TR << "New fold tree: " << newF << std::endl;
}

////
//// set up foldtree, variants, movemap, etc.
void guess_rbsegs_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< RBSegment > & rigid_segs,
	utility::vector1< RBSegment > & rb_chunks,
	protocols::loops::Loops & loops
) {
	core::Size nres = pose.total_residue()-1; // terminal VRT

	rigid_segs.clear();
	rb_chunks.clear();
	loops.clear();

	// dssp parse
	protocols::jumping::Dssp secstruct( pose );
	ObjexxFCL::FArray1D< char > dssp_pose( nres );
	secstruct.dssp_reduced (dssp_pose);
	//secstruct.insert_ss_into_pose( pose );

	// find all helices > 5 residues
	//          strands > 3 residues
	utility::vector1< RBSegment > simple_segments;
	bool in_helix=false, in_strand=false;
	int ss_start = -1;
	for (int i=1; i<=(int)nres; ++i) {
		// strand end
		if (dssp_pose(i) != 'E' && in_strand) {
			in_strand = false;
			if (i-ss_start >= 3)
				simple_segments.push_back( RBSegment( ss_start, i-1, 'E' ) );
		}
		// helix end
		if (dssp_pose(i) != 'H' && in_helix) {
			in_helix = false;
			if (i-ss_start >= 5)
				simple_segments.push_back( RBSegment( ss_start, i-1, 'H' ) );
		}
		// strand start
		if (dssp_pose(i) == 'E' && !in_strand)  {
			in_strand = true;
			ss_start = i;
		}
		// helix start
		if (dssp_pose(i) == 'H' && !in_helix)  {
			in_helix = true;
			ss_start = i;
		}
	}

	// put at least 2 "loop residues" inbetween each RB segment.
	// always eat into the helix (even if this leaves a helix < 3 reses)
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (simple_segments[i+1][1].start() - simple_segments[i][1].end() <= 2) {
			if (simple_segments[i][1].char_type() == 'H') {
				simple_segments[i][1].set_end( simple_segments[i+1][1].start()-3 );
			} else if (simple_segments[i+1][1].char_type() == 'H') {
				simple_segments[i+1][1].set_start( simple_segments[i][1].end()+3 );
			} else {
				// eat into longer strand (will only need to chomp 1 res)
				if (simple_segments[i+1][1].length() > simple_segments[i][1].length() )
					simple_segments[i+1][1].set_start( simple_segments[i][1].end()+3 );
				else
					simple_segments[i][1].set_end( simple_segments[i+1][1].start()-3 );
			}
		}
	}

	// check for "1-residue" loops at termini; extend RB segments if necessary
	if (simple_segments[1][1].start() == 2) simple_segments[1][1].set_start(1);
	if (simple_segments[simple_segments.size()][1].end() == nres-1) simple_segments[simple_segments.size()][1].set_end(nres);


	// auto-gen loops
	if ( simple_segments.size() > 1 ) {
		std::sort( simple_segments.begin(), simple_segments.end(), RB_lt());
		int start_res=1, end_res=simple_segments[1][1].start()-1;
		int cutpt = (start_res+end_res)/2;
		int nsegs = simple_segments.size();

		if (end_res >= start_res)
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
		for (int i=1; i<nsegs; ++i) {
			start_res = simple_segments[i][1].end()+1;
			end_res   = simple_segments[i+1][1].start()-1;
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
		}
		start_res = simple_segments[nsegs][1].end()+1;
		end_res   = nres;
		cutpt = (start_res+end_res)/2;
		if (end_res >= start_res)
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );

		// TODO: split loops on cutpoints from original pose
	}


	// now combine paired strands into a compound segment
	//   look for NH--O distance < 2.6A (?)
	utility::vector1< utility::vector1< int > > compound;
	utility::vector1< core::Size > parent_seg( simple_segments.size(), 0 );
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (simple_segments[i][1].char_type() != 'E') continue;

		utility::vector1< int > this_seg( 1, i );
		for (int j=i+1; j<=(int)simple_segments.size(); ++j) {
			if (simple_segments[j][1].char_type() != 'E') continue;

			// foreach res in i,j
			bool found = false;
			for (int ii=(int)simple_segments[i][1].start(); ii<=(int)simple_segments[i][1].end() && !found; ++ii)
			for (int jj=(int)simple_segments[j][1].start(); jj<=(int)simple_segments[j][1].end() && !found; ++jj) {
				core::Real d2=10;

				if (pose.residue(ii).aa() != core::chemical::aa_pro && pose.residue(jj).aa() != core::chemical::aa_pro)
					d2 = std::min(
						(pose.residue(ii).atom("H").xyz() - pose.residue(jj).atom("O").xyz()).length_squared() ,
						(pose.residue(ii).atom("O").xyz() - pose.residue(jj).atom("H").xyz()).length_squared() );
				else if (pose.residue(jj).aa() != core::chemical::aa_pro)
					d2 = (pose.residue(ii).atom("O").xyz() - pose.residue(jj).atom("H").xyz()).length_squared();
				else if (pose.residue(ii).aa() != core::chemical::aa_pro)
					d2 = (pose.residue(jj).atom("O").xyz() - pose.residue(ii).atom("H").xyz()).length_squared();

				if (d2 < 2.6*2.6) {
					this_seg.push_back(j);
					if (parent_seg[i] == 0)
						parent_seg[i] = i;
					if (parent_seg[j] == 0)
						parent_seg[j] = parent_seg[i];
					else {
						// tricky case ... j is already mapped
						// in this case map everything mapped to i to parent_seg[j]
						for (int k=1; k<j; ++k)
							if ((int)parent_seg[k] == i) parent_seg[k] = parent_seg[j];
					}

					TR << "Merging " << j << " (" << simple_segments[j][1].start() << "," << simple_segments[j][1].end() << ") ";
					TR << " to " << i << " (" <<  simple_segments[i][1].start() << "," << simple_segments[i][1].end() << ") " << std::endl;
					found = true;
				}
			}
		}
	}

	// make the compound segments
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (((int)parent_seg[i]) != i) continue; // not compound or already added
		utility::vector1< RBSegment > thisLockSeg;
		for (core::Size j=i; j<=simple_segments.size(); ++j) {
			if ( ((int)parent_seg[j]) == i )
				thisLockSeg.push_back( simple_segments[ j ] );
		}
		rigid_segs.push_back( RBSegment( thisLockSeg ) );
	}
	// add in all other simple segs
	for (int i=1; i<=(int)simple_segments.size(); ++i)
		if (parent_seg[i] == 0)
			rigid_segs.push_back( simple_segments[i] );

	// sort loops & rbsegs, choose cutpoints
	std::sort( rigid_segs.begin(), rigid_segs.end(), RB_lt());
	std::sort( loops.v_begin(), loops.v_end(), protocols::loops::Loop_lt());
	loops.auto_choose_cutpoints( pose );

	// define chunks
	rb_chunks = rigid_segs;
	for (int i=1; i<=(int)rigid_segs.size(); ++i) {
		for (int j=1; j<=(int)rigid_segs[i].nContinuousSegments() ; ++j ) {
			core::Size c_low=1, c_high=nres;
			for (int k=1; k<=(int)loops.size(); ++k) {
				if (loops[k].cut() < rigid_segs[i][j].start() && loops[k].cut() > c_low )
					c_low = loops[k].cut();
				if (loops[k].cut()+1 > rigid_segs[i][j].end() && loops[k].cut()+1 < c_high )
					c_high = loops[k].cut()+1;
			}
			rb_chunks[i][j].set_start(c_low);
			rb_chunks[i][j].set_end(c_high);
		}
	}
}

//// remap all segs
void remap_rb_segments(
                 utility::vector1< RBSegment > const &rbsegs,
                 utility::vector1< RBSegment > &rbsegs_remap,
                 core::id::SequenceMapping const &resmap ) {
	for ( RBConsIt it_seg = rbsegs.begin(), it_seg_end = rbsegs.end(); it_seg != it_seg_end; ++it_seg ) {
		RBSegment it_remap = it_seg->remap( resmap );
		rbsegs_remap.push_back( it_remap );
	}
}



}
}
