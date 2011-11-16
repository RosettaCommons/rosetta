// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/DomainAssembly.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/DomainAssembly.hh>
#include <protocols/protein_interface_design/movers/DomainAssemblyCreator.hh>

// Package headers

// Project headers
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/moves/ReturnSidechainMover.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/moves/SwitchResidueTypeSetMover.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
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
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/types.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
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
#include <core/scoring/ScoreFunction.hh>
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
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
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
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/mpistream.hh>
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
#include <utility/tag/Tag.fwd.hh>
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
#include <ios>
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
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers

// Auto-header: duplicate removed #include <core/scoring/ScoreFunctionFactory.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.DomainAssembly" );

std::string
DomainAssemblyCreator::keyname() const
{
	return DomainAssemblyCreator::mover_name();
}

protocols::moves::MoverOP
DomainAssemblyCreator::create_mover() const {
	return new DomainAssembly;
}

std::string
DomainAssemblyCreator::mover_name()
{
	return "DomainAssembly";
}

DomainAssembly::DomainAssembly() :
	protocols::moves::Mover( DomainAssemblyCreator::mover_name() )
{}

DomainAssembly::DomainAssembly(
	core::Size const linker_start,
	core::Size const linker_end,
	FragSetOP fragset_large,
	FragSetOP fragset_small
) :
	protocols::moves::Mover( DomainAssemblyCreator::mover_name() ),
	linker_start_( linker_start ),
	linker_end_( linker_end ),
	fragset_large_( fragset_large ),
	fragset_small_( fragset_small ),
	fragments_set_( true )
{}

DomainAssembly::~DomainAssembly() {}

protocols::moves::MoverOP
DomainAssembly::clone() const {
	return( protocols::moves::MoverOP( new DomainAssembly( *this ) ) );
}

protocols::moves::MoverOP
DomainAssembly::fresh_instance() const
{
	return protocols::moves::MoverOP( new DomainAssembly );
}

void
DomainAssembly::apply( core::pose::Pose & pose )
{
	protocols::moves::ReturnSidechainMover recover_sidechains( pose );

	if( !fragments_set_ ) {
		TR<<"*******WARNING WARNING********: fragments not set, skipping domain assembly"<<std::endl;
		runtime_assert( fragments_set_ );
		return;
	}

	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap();
	mm->set_bb( false );
	for( core::Size i=linker_start_; i<=linker_end_; ++i ) mm->set_bb( i, true );
/* //The following code only makes sense if we do 'blind' prediction stuff. In all cases where we have a starting structure,
//it's probably a good idea to start near that.
	for(i = 1; i <= flexible_regions_.size(); i++ ){
		ir = std::max(0, std::min( (int)pose_full_centroid.total_residue() ,  (int) flexible_regions_[i] ) );

		std::cout << "Linker: " << ir << std::endl;

		Real const init_phi  ( -150.0 );
		Real const init_psi  (  150.0 );
		Real const init_omega(  180.0 );
		pose_full_centroid.set_phi( ir ,  init_phi  );
		pose_full_centroid.set_psi( ir ,  init_psi  );
		pose_full_centroid.set_omega( ir, init_omega);
		movemap->set_bb(ir, true);
	}
*/
	core::scoring::ScoreFunctionCOP scorefxn( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	protocols::moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
	protocols::moves::SwitchResidueTypeSetMover to_fullatom( core::chemical::FA_STANDARD );
	to_centroid.apply( pose );
	protocols::abinitio::AbrelaxApplication abrelax_app;
	do {
		protocols::abinitio::ClassicAbinitio abinit( fragset_large_, fragset_small_, mm );
		abinit.init( pose );
		abinit.apply( pose );
	} while( !abrelax_app.check_filters( pose ) );

	//recover sidechains from starting structures
	to_fullatom.apply( pose );
	recover_sidechains.apply( pose );
//	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
//	(*scorefxn)( pose );
//	scorefxn->accumulate_residue_total_energies( pose );

//Repack regions around the linker
	core::pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	task->restrict_to_repacking();

	for ( core::Size i = linker_end_+1; i <= pose.total_residue(); ++i) {
		if ( !pose.residue(i).is_protein() ) continue;
		if( pose.residue(i).type().name() == "CYD" ) {
			task->nonconst_residue_task( i ).prevent_repacking();
			continue;
		}

		core::conformation::Residue const resi( pose.residue( i ) );
		core::Size j;
		for( j = 1; j<=linker_end_; ++j ) {
			core::conformation::Residue const resj( pose.residue( j ) );

			core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if( distance <= 8.0 ) break;
		}
		if( j>linker_end_ ) task->nonconst_residue_task( i ).prevent_repacking();
	}
	for ( core::Size i = 1; i <= linker_start_ - 1; ++i) {
		if ( !pose.residue(i).is_protein() ) continue;
		if( pose.residue(i).type().name() == "CYD" ) {
			task->nonconst_residue_task( i ).prevent_repacking();
			continue;
		}

		core::conformation::Residue const resi( pose.residue( i ) );
		core::Size j;
		for( j = linker_start_; j<=pose.total_residue(); ++j ) {
			core::conformation::Residue const resj( pose.residue( j ) );

			core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if( distance <= 8.0 ) break;
		}
		if( j>pose.total_residue() ) task->nonconst_residue_task( i ).prevent_repacking();
	}
	//in case there is a resfile, information in this resfile overrides the computed task
	if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() )
		core::pack::task::parse_resfile(pose, *task);

	pack::pack_rotamers( pose, *scorefxn, task);
	(*scorefxn)( pose );
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
}

std::string
DomainAssembly::get_name() const {
	return DomainAssemblyCreator::mover_name();
}

void
DomainAssembly::parse_my_tag( TagPtr const tag, DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	linker_start_ = protocols::rosetta_scripts::get_resnum( tag, pose, "linker_start_" );
	linker_end_   = protocols::rosetta_scripts::get_resnum( tag, pose, "linker_end_" );
	runtime_assert( linker_end_ > linker_start_ );
	runtime_assert( linker_start_ > 0 );
	runtime_assert( linker_end_ < pose.total_residue() );

	std::string const frag_large_fname( tag->getOption< std::string >( "frag9", "frag9" ) );
	std::string const frag_small_fname( tag->getOption< std::string >( "frag3", "frag3" ) );

	using namespace core::fragment;
	using namespace basic::options;
	fragset_large_ = FragmentIO(option[ OptionKeys::abinitio::number_9mer_frags ] ).read_data( frag_large_fname );
	fragset_small_ = FragmentIO(option[ OptionKeys::abinitio::number_3mer_frags ] ).read_data( frag_small_fname );
	fragments_set_ = true;
}

} //movers
} //protein_interface_design
} //protocols
