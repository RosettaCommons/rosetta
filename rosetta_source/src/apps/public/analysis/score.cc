// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/moves/ScoreMover.hh>
#include <protocols/moves/TailsScoreMover.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <protocols/jobdist/standard_mains.hh>

#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/NullMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/ProlineFixMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/ConstraintSetMover.hh>
// AUTO-REMOVED #include <protocols/electron_density/util.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>

// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <core/scoring/ScoreFunctionFactory.hh> // getScoreFunction
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/SuperimposeMover.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/krassk.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/Frame.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.fwd.hh>
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
// AUTO-REMOVED #include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
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
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
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
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/ozstream.fwd.hh>
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
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
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
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
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
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
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
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

//Auto Headers


//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end



using namespace core;
using namespace basic::options;

static basic::Tracer TR("protocols.moves.ScoreMover");

namespace score_app { BooleanOptionKey linmin( "score_app:linmin" );
											BooleanOptionKey superimpose_to_native( "score_app:superimpose_to_native" );	}

int
main( int argc, char * argv [] )
{
	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;

	ScoreMover::register_options();
	protocols::jobdist::register_options_universal_main();
  option.add_relevant( in::file::fullatom );
  option.add_relevant( relax::fast );
	option.add( score_app::linmin, "Do a quick round of linmin before reporting the score" );
  option.add_relevant( score_app::linmin );
  option.add_relevant( out::output                  );
  option.add_relevant( out::nooutput                );
  option.add_relevant( in::file::fullatom           );
  option.add_relevant( rescore::verbose             );
  option.add_relevant( in::file::repair_sidechains  );
	option.add( score_app::superimpose_to_native, "superimpose structure to native" );


// scoring should by default not produce output files - that's so annoying
// unless of coures the user insists.

	// initialize core
	devel::init(argc, argv);

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << " Rosetta Tool:   score   -  rescores PDBs and silent files, extracts PDBs from silent files, assembles PDBs into silent files. " << std::endl;
	std::cout << " Usage:                                                                  " << std::endl;
	std::cout << "   PDB input:      -in:file:s *.pdb   or  " << std::endl;
	std::cout << "                   -in:file:l  list_of_pdbs  " << std::endl;
	std::cout << "                   -no_optH                                    Dont change positions of Hydrogen atoms! (default true, specify false if you want optH)" << std::endl;
	std::cout << "   Silent input:   -in:file:silent silent.out                  silent input filesname " << std::endl;
	std::cout << "                   -in:file:tags                               specify specific tags to be extracted, if left out all will be taken " << std::endl;
	std::cout << "                   -in:file:fullatom                           for full atom structures " << std::endl;
	std::cout << "                   -in:file:binary_silentfile                  for non-ideal structures (such as from looprelax) " << std::endl;
	std::cout << "                   -in:file:silent_optH                        Call optH when reading silent files (useful for HisD/HisE determination)" << std::endl;
	std::cout << "                   -score_app:linmin                           Run a quick linmin before scoring" << std::endl;
	std::cout << "   Native:         -in:file:native                             native PDB (rms, maxsub and gdtm scores will be calculated)" << std::endl;
	std::cout << "   Scorefunction:  -score:weights  weights                     weight set or weights file " << std::endl;
	std::cout << "                   -score:patch  patch                         patch set " << std::endl;
	std::cout << "                   -score:optH_weights                         Weights file for optH (default standard.wts w/ sc12 patch)" << std::endl;
	std::cout << "                   -score:optH_patch                           Weights patch file for optH" << std::endl;
	std::cout << "                   -rescore:verbose                            display score breakdown " << std::endl;
	std::cout << "   Output:         -out:nooutput                               don't print PDB structures (default now) " << std::endl;
	std::cout << "                   -out:output                                 force printing of PDB structures " << std::endl;
	std::cout << "                   -out:file:silent                            write silent-out file " << std::endl;
	std::cout << "                   -out:file:scorefile name                    write scorefile (default default.sc)" << std::endl;
	std::cout << "                   -out:prefix  myprefix                       prefix the output structures with a string " << std::endl;
	std::cout << "  Examples: " << std::endl;
	std::cout << "   score  -database ~/minirosetta_database -in:file:silent silent.out -in::file::binary_silentfile -in::file::fullatom -native 1a19.pdb " << std::endl;
	std::cout << "   Will rescore all structures in silent.out, in full atom mode and accounting for nonideal structure if present. Additionally " << std::endl;
	std::cout << " 	 it will print a PDB for every structure with -out:output flag " << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	//The following lines are to ensure one can rescore the pcs energy term (that uses TopologyClaimer)
	if ( option[ broker::setup ].user() ) {
		protocols::topology_broker::TopologyBrokerOP top_bro_OP = new  topology_broker::TopologyBroker();
		try {
			add_cmdline_claims(*top_bro_OP, false /* do_I_need_fragments */);
		}
		catch ( utility::excn::EXCN_Exception &excn )  {
			excn.show( TR.Error );
			utility_exit();
		}
	}

	// do not output pdb by default, unless with -out:output flag
	if( !option[ out::output ].user() ){
		option[ out::nooutput ].value( true );
	}

	// get scorefxn and add constraints if defined
	core::scoring::ScoreFunctionOP sfxn = core::scoring::getScoreFunction();
	if( option[ in::file::fullatom ]() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sfxn );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *sfxn );
	}

	// now add density scores from cmd line
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sfxn );
	}

	// create a ScoreMover

	ScoreMover* scoretmp;
	if( option[ krassk::tail_mode])
	   {
			scoretmp = new moves::TailsScoreMover(sfxn);
	   }
	   else
	   {
			scoretmp = new moves::ScoreMover(sfxn);
	   }

	if(  option[ rescore::verbose ] )	{
		scoretmp->set_verbose( true );
	} else {
		scoretmp->set_verbose( false );
	}

	// save it to a mover that will be passed to job_distributor
	MoverOP mover = scoretmp;

	// do sth more than just scoring
	if ( option[ score_app::linmin ]() || option[ in::file::repair_sidechains ]() ) {
		assert( sfxn );
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		if ( option[ in::file::repair_sidechains ]() ) {
			protocols::moves::ProlineFixMoverOP pfm = new ProlineFixMover;
			seqmov->add_mover( pfm );
		}
		if ( option[ score_app::linmin ]() ) {
			core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
			movemap->set_bb( true ); movemap->set_chi( true );
			protocols::moves::MinMoverOP minmover = new protocols::moves::MinMover(
				movemap, sfxn, "linmin", 1e-4,
				true /*use_nblist*/, false /*deriv_check*/, false /*verbose driv check*/ );
			seqmov->add_mover( minmover );
		}

		seqmov->add_mover( mover );
		mover = seqmov;
	} // if infile remediation necessary

	// add constraints from cmd line
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() || option[ OptionKeys::constraints::cst_file ].user() ) {
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
			protocols::moves::ConstraintSetMoverOP loadCsts( new protocols::moves::ConstraintSetMover );
			if( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			} else {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
			}
			seqmov->add_mover( loadCsts );
			seqmov->add_mover( mover );
			mover = seqmov;
	}

	// set pose for density scoring if a map was input
	//   + (potentially) dock map into density
	if ( option[ edensity::mapfile ].user() ) {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::electron_density::SetupForDensityScoringMover );
		seqmov->add_mover( mover );
		mover = seqmov;
	}

	// set pose for symmetry
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::moves::symmetry::SetupForSymmetryMover );
		seqmov->add_mover( mover );
		mover = seqmov;
	}

	if ( option[ score_app::superimpose_to_native ]() ) {
		if ( !option[ in::file::native ].user() ) {
				TR << "No native specified. Cannot align to native..." << '\n';
		} else {
			// read native structure
			core::pose::Pose native;
			core::import_pose::pose_from_pdb( native, option[ basic::options::OptionKeys::in::file::native ] );
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
			seqmov->add_mover( new protocols::moves::SuperimposeMover( native ) );
			seqmov->add_mover( mover );
			mover = seqmov;
		}
	}

	// operate this mover and output pdbs/scorefile
	protocols::jobdist::universal_main( *mover );

	return 0;
}

