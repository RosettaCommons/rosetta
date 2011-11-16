// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NoesyModule.cc
/// @brief main hook-up for the automatic NOESY assignment module
/// @detailed
///	  handling of input-output options
///   class NoesyModule:
///       read input files
///       write assignments, constraints
///       run assignment
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/NoesyModule.hh>
//o a;kdfj h hd dd
// Package Headers
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/PeakFileFormat_Sparky.hh>
//#include <devel/noesy_assign/DistanceScoreMover.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh> // REQUIRED FOR WINDOWS

// for switching residue type set to centroid
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.fwd.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/prof.hh>

//// C++ headers
// AUTO-REMOVED #include <math.h> //for isnan



OPT_2GRP_KEY( File, noesy, in, resonances )
OPT_2GRP_KEY( FileVector, noesy, in, peaks )
OPT_2GRP_KEY( FileVector, noesy, in, peak_resonance_pairs )
OPT_2GRP_KEY( Boolean, noesy, in, use_assignments )
OPT_2GRP_KEY( File, noesy, in, decoys )

OPT_2GRP_KEY( File, noesy, out, resonances )
OPT_2GRP_KEY( File, noesy, out, peaks )
OPT_2GRP_KEY( File, noesy, out, talos )
OPT_2GRP_KEY( Boolean, noesy, out, names )
OPT_2GRP_KEY( Boolean, noesy, out, separate_peak_files )
OPT_2GRP_KEY( Boolean, noesy, out, unambiguous )
OPT_2GRP_KEY( String, noesy, out, format )
OPT_2GRP_KEY( Real, noesy, out, minVC )

OPT_1GRP_KEY( Boolean, noesy, no_decoys )

// OPT_1GRP_KEY( Boolean, noesy, no_symm )
// OPT_1GRP_KEY( Boolean, noesy, no_cs )
// OPT_1GRP_KEY( Boolean, noesy, no_upper )
// OPT_1GRP_KEY( Boolean, noesy, no_remove_diagonal )
// OPT_1GRP_KEY( Boolean, noesy, no_calibrate )

bool protocols::noesy_assign::NoesyModule::options_registered_( false );


bool protocols::noesy_assign::NoesyModule::cmdline_options_activated() {
	return basic::options::option[ basic::options::OptionKeys::noesy::in::resonances ].user()
		&& basic::options::option[ basic::options::OptionKeys::noesy::in::peaks ].user();
}

void protocols::noesy_assign::NoesyModule::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;

  PeakAssignmentParameters::register_options();
	PeakFileFormat_xeasy::register_options();
  if ( options_registered_ ) return;
  options_registered_ = true;

  //....
  NEW_OPT( noesy::in::resonances, "file with assigned chemical shifts", "" );

  NEW_OPT( noesy::in::peaks, "file with noesy peaks", "" );
	NEW_OPT( noesy::in::peak_resonance_pairs, "pairs of files that belong together: cc.peaks cc.prot ilv.peaks ilv.prot", "" );
  NEW_OPT( noesy::in::use_assignments, "when reading peaks the already existing assignments are not ignored", false );
  NEW_OPT( noesy::in::decoys, "silent file with decoys used for 3D structural compatibility test", "" );

  NEW_OPT( noesy::out::resonances, "the parsed resonances file with translated atom names etc.", "cs_out.dat" );
  NEW_OPT( noesy::out::peaks, "the parsed peaks file with assignments", "NOE_out.dat" );
	NEW_OPT( noesy::out::talos, "write the resonances also as talos file", "cs_prot.tab" );
	NEW_OPT( noesy::out::names, "write atom-names rather then resonance ID for assignments", true );
	NEW_OPT( noesy::out::separate_peak_files, "write peaks to independent files with out::peaks as prefix", false );
	NEW_OPT( noesy::out::unambiguous, "write only the assignment with hightest VC", false );
	NEW_OPT( noesy::out::minVC, "write only assignments that contribute more than X to the peak-volume", 0.0 );
	NEW_OPT( noesy::out::format, "write as xeasy or sparky file", "xeasy" );
  NEW_OPT( noesy::no_decoys, "check comp. with decoys", false );

//   NEW_OPT( noesy::no_symm, "check comp. with decoys", false );
//   NEW_OPT( noesy::no_cs, "check comp. with decoys", false );
//   NEW_OPT( noesy::no_upper, "check upper", false );
//   NEW_OPT( noesy::no_remove_diagonal, "", false );
//   NEW_OPT( noesy::no_calibrate, "don't calibrate the NOE distance bound", false );


}


static basic::Tracer tr("protocols.noesy_assign.NoesyModule");

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
//using namespace OptionKeys;
///// templates

#include <protocols/noesy_assign/CrossPeakList.impl.hh>
#include <protocols/noesy_assign/NoesyModule.impl.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/XYZ_Func.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
// AUTO-REMOVED #include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.fwd.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/PairingsList.fwd.hh>
#include <protocols/jumping/PairingsList.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/noesy_assign/CrossPeak.fwd.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakInfo.fwd.hh>
#include <protocols/noesy_assign/CrossPeakInfo.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/DistanceScoreMover.fwd.hh>
#include <protocols/noesy_assign/DistanceScoreMover.hh>
#include <protocols/noesy_assign/NoesyModule.fwd.hh>
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>
// Auto-header: duplicate removed #include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.hh>
#include <protocols/noesy_assign/PeakFileFormat.fwd.hh>
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
#include <protocols/noesy_assign/StructureDependentPeakCalibrator.hh>
#include <protocols/noesy_assign/StructureIndependentPeakCalibrator.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
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
#include <numeric/random/random.fwd.hh>
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
#include <bitset>
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
#include <time.h>
// AUTO-REMOVED #include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers




namespace protocols {
namespace noesy_assign {


///Constructor   - read input files / requires options to be initialized
NoesyModule::NoesyModule( std::string const& fasta_sequence ) :
  crosspeaks_( NULL ),
  main_resonances_( new ResonanceList( fasta_sequence ) )
{
	read_input_files();
	runtime_assert( options_registered_ );
	// moved this option into PeakAssignmentParameters.cc
	//	skip_network_analysis_ = basic::options::option[  options::OptionKeys::noesy::no_network ]();
}


///delete all data and read input files again...  aka fresh_instance()
void NoesyModule::reset() {
	crosspeaks_ = NULL;
	main_resonances_ = new ResonanceList( main_resonances_->sequence() );
	read_input_files();
}

///read peak- and resonance files
void NoesyModule::read_input_files() {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_READ_INPUT );

	//read resonances
	if ( option[ options::OptionKeys::noesy::in::resonances ].user() ) {
		utility::io::izstream input_file( option[ options::OptionKeys::noesy::in::resonances ]() );
		utility::io::ozstream output_file( option[ options::OptionKeys::noesy::out::resonances ]() );
		if ( input_file.good() ) {
			main_resonances_->read_from_stream( input_file );
			main_resonances_->write_to_stream( output_file );
			if ( option[ options::OptionKeys::noesy::out::talos ].user() ) {
				utility::io::ozstream talos_file( option[ options::OptionKeys::noesy::out::talos ]() );
				main_resonances_->write_talos_format( talos_file, true /*backbone only*/ );
			}
		} else {
			tr.Error << "cannot read " << input_file << std::endl;
			throw utility::excn::EXCN_FileNotFound( option[ options::OptionKeys::noesy::in::resonances ]() );
		}
	}

  { //scope
		//read peak lists
		crosspeaks_ = new CrossPeakList();
		Size nfiles( option[ options::OptionKeys::noesy::in::peaks ]().size() );

		//loop-all files: read one file at a time
		for ( core::Size ifile = 1; ifile <= nfiles; ++ifile ) {  //as many as there are files
			std::string file( option[ options::OptionKeys::noesy::in::peaks ]()[ ifile ] );
			utility::io::izstream input_file( file );
			if ( input_file.good() ) {
				if ( main_resonances_->size() < 1 ) {
					throw utility::excn::EXCN_BadInput( "attempt to read Peak-Files with option -noesy:in:peaks without a global resonance file: -noesy:in:resonances" );
				}
				PeakFileFormat_xeasy format;
				format.set_filename( option[ options::OptionKeys::noesy::in::peaks ]()[ ifile ].base() );
				format.set_ignore_assignments( !option[ options::OptionKeys::noesy::in::use_assignments ]() );
				tr.Info << "reading " << file << "... " << std::endl;
				crosspeaks_->read_from_stream( input_file, format, main_resonances_ );
			} else {
				tr.Error << "cannot read " << file << std::endl;
			}
		}
	} // scope end

	//specific resonances for a given peak-file
	if ( option[ options::OptionKeys::noesy::in::peak_resonance_pairs ].user() ) { //scope
		Size n_pair_files(  option[ options::OptionKeys::noesy::in::peak_resonance_pairs ]().size() );
		if ( n_pair_files % 2 != 0 ) {
			throw utility::excn::EXCN_BadInput( "odd number of entries in option -noesy:in:peak_resonance_pairs, always provide pairs of files  <*.peaks> <*.prot>" );
		}
		for ( core::Size ifile = 1; ifile <= n_pair_files; ifile += 2 ) {
			ResonanceListOP resonances = new ResonanceList( main_resonances_->sequence() );
			utility::io::izstream res_input_file( option[ options::OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile+1 ] );
			if ( res_input_file.good() ) {
				resonances->read_from_stream( res_input_file );
			} else {
				tr.Error << "cannot read " << res_input_file << std::endl;
				throw utility::excn::EXCN_FileNotFound( option[ options::OptionKeys::noesy::in::resonances ]() );
			}
			std::string file( option[ options::OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile ] );
			utility::io::izstream input_file( file );
			if ( input_file.good() ) {
				PeakFileFormat_xeasy format;
				format.set_filename( option[ options::OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile ].base() );
				format.set_ignore_assignments( !option[ options::OptionKeys::noesy::in::use_assignments ]() );
				tr.Info << "reading " << file << "... " << std::endl;
				crosspeaks_->read_from_stream( input_file, format, resonances );
			} else {
				tr.Error << "cannot read " << file << std::endl;
			}
		}
	} // scope end
}

///@brief write peak assignments into peak-file (sparky, cyana)
void NoesyModule::write_assignments( std::string file_name ) {
	ProfileThis doit( NOESY_ASSIGN_WRITE_ASSIGNMENTS );
	if ( file_name == "use_cmd_line" ) {
		file_name = option[ options::OptionKeys::noesy::out::peaks ]();
	}
	PeakFileFormatOP format;
	std::string const format_str( option[ options::OptionKeys::noesy::out::format ]() );
	if ( format_str == "xeasy" ) {
		format = new PeakFileFormat_xeasy( main_resonances_ );
	} else if ( format_str == "sparky" ) {
		format = new PeakFileFormat_Sparky( main_resonances_ );
	} else utility_exit_with_message( "NOE_data output format "+format_str+" is not known! ");
	format->set_write_atom_names( option[ options::OptionKeys::noesy::out::names ]() );
	format->set_write_only_highest_VC( option[ options::OptionKeys::noesy::out::unambiguous ]() );
	format->set_min_VC_to_write( option[ options::OptionKeys::noesy::out::minVC ]() );
	if ( option[ options::OptionKeys::noesy::out::separate_peak_files ]() ) {
		crosspeaks_->write_peak_files( file_name, *format );
	} else {
		utility::io::ozstream output_file( file_name );
		crosspeaks_->write_to_stream( output_file, *format );
	}
}

///@brief assign peaks ( no explicit decoys - wrapper )
void NoesyModule::assign( Size cycle ) {
  using namespace options;
  using namespace options::OptionKeys::noesy;
  core::io::silent::SilentFileData sfd;
  if ( !option[ no_decoys ]() ) {
		if ( option[ in::decoys ].user() ) sfd.read_file( option[ in::decoys ]() );
		if ( sfd.size() == 0 && option[ OptionKeys::in::file::silent ].user() ) sfd.read_file( option[ OptionKeys::in::file::silent ]()[ 1 ] );
  }
  assign( sfd.begin(), sfd.end(), cycle );
}

///@brief generate constraint files from assignments
void NoesyModule::generate_constraint_files(
	 core::pose::Pose const& pose,
	 std::string const& cst_fa_file,
	 std::string const& cst_centroid_file,
   core::Size min_seq_separation
) const {

	PROF_START( NOESY_ASSIGN_GEN_CST );

	using namespace core::scoring::constraints;
  core::pose::Pose centroid_pose = pose;
	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID );

	ConstraintSetOP cstset = new ConstraintSet;
	ConstraintSetOP centroid_cstset = new ConstraintSet;
	tr.Info << "generate constraints..." << std::endl;
	crosspeaks_->generate_fa_and_cen_constraints( cstset, centroid_cstset, pose, centroid_pose, min_seq_separation );

	PROF_STOP( NOESY_ASSIGN_GEN_CST );

	tr.Info << "write constraints..." << std::endl;
	PROF_START( NOESY_ASSIGN_WRITE_CST );

	ConstraintIO::write_constraints( cst_fa_file, *cstset, pose );
  ConstraintIO::write_constraints( cst_centroid_file, *centroid_cstset, centroid_pose );

	PROF_STOP( NOESY_ASSIGN_WRITE_CST );
}

}
}


