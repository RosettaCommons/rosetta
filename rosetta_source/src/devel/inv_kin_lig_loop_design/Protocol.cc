// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Protocol.cc
///
/// @brief
/// @author
#include <devel/inv_kin_lig_loop_design/Protocol.hh>
#include <devel/inv_kin_lig_loop_design/Mover.hh>
#include <devel/inv_kin_lig_loop_design/JumpManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/ccd_closure.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/pack/pack_rotamers.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
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
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
#include <core/coarse/Translator.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/DOF_ID_Mask.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
//#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/types.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
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
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <devel/inv_kin_lig_loop_design/Cloner.hh>
#include <devel/inv_kin_lig_loop_design/Fragment.hh>
#include <devel/inv_kin_lig_loop_design/Loop.hh>
#include <devel/inv_kin_lig_loop_design/ResID.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/string_util.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
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
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
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
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/CArray.fwd.hh>
#include <ObjexxFCL/CArrayP.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1.all.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2.all.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3P.fwd.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4P.fwd.hh>
#include <ObjexxFCL/FArray4.all.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArray5.fwd.hh>
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5P.fwd.hh>
#include <ObjexxFCL/FArray5.all.fwd.hh>
// Auto-header: duplicate removed #include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/shared_ptr.hpp>



#define REPORTCYCLES1(c0,m0)                    " Cycles:\t" << c0 << "/" << m0
#define REPORTCYCLES2(c0,m0,c1,m1)              " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1
#define REPORTCYCLES3(c0,m0,c1,m1,c2,m2)        " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1 << "\t" << c2 << "/" << m2
#define REPORTCYCLES4(c0,m0,c1,m1,c2,m2,c3,m3)  " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1 << "\t" << c2 << "/" << m2 << "\t" << c3 << "/" << m3


namespace devel {

  namespace inv_kin_lig_loop_design {


    namespace {

      bool contains(vector<Loop> const& loops, int seqpos ) {
				for( Size k = 0; k < loops.size(); ++k ) {
					if( loops[k].lo <= seqpos && seqpos <= loops[k].hi ) {
						return true;
					}
				}
				return false;
      }

      bool is_anchor(vector<Loop> const& loops, int seqpos ) {
				for( Size k = 0; k < loops.size(); ++k ) {
					if( loops[k].to == seqpos ) {
						return true;
					}
				}
				return false;
      }

      core::kinematics::MoveMap get_move_map(int num_res, vector<Loop> const& loops) {
				core::kinematics::MoveMap rval;
				for( int i = 1; i <= num_res; ++i ) {
					rval.set_bb( i, contains(loops,i) );
					rval.set_chi( i, is_anchor(loops,i) );

					//cout << "get_move_map: " << i << " " << contains(loops,i) << " " << is_anchor(loops,i) << endl;

				}
				rval.set_jump( false );
				return rval;
      }

    }

    //     Protocol::Protocol() {
    //     }

    // const double Protocol::DEFAULT_MIN_ACCEPTANCE_RATE = 5e-2;
    //     const double Protocol::DEFAULT_MIN_WEIGHT_REP      = 1e-2;

    void Protocol::setParams(TagPtr tag0, vector<Loop> const& loops_) {

      loops = loops_;

      max_cycles               = tag0->getTag("param")->getOption<int>("max_cycles",10);
      max_cycles_start         = tag0->getTag("param")->getOption<int>("max_cycles_start",10);
      max_cycles_start_3mer    = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer",1000);
      max_cycles_start_3mer_anchor = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer_anchor",100);
      max_cycles_start_3mer_1mer   = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer_1mer",100);
      max_cycles_start_1mer    = tag0->getTag("param")->getOption<int>("max_cycles_start_1mer",100);
      max_cycles_ramp          = tag0->getTag("param")->getOption<int>("max_cycles_ramp",10);
      max_cycles_ramp_sm       = tag0->getTag("param")->getOption<int>("max_cycles_ramp_sm",1);
      max_cycles_ramp_sm_min   = tag0->getTag("param")->getOption<int>("max_cycles_ramp_sm_min",10);
      max_cycles_design        = tag0->getTag("param")->getOption<int>("max_cycles_design",10);
      max_cycles_design_sm     = tag0->getTag("param")->getOption<int>("max_cycles_design_sm",10);
      max_cycles_design_sm_min = tag0->getTag("param")->getOption<int>("max_cycles_design_sm_min",10);

    } // Protocol::setParams

    void Protocol::phase_lores() {

      //setLoRes();

      core::optimization::AtomTreeMinimizer atm;

      double temp = 1;

      core::pose::PoseOP pose_init( new core::pose::Pose );
      (*pose_init) = (*pose);

      protocols::moves::MonteCarlo min_start(*pose, *score_fxn_lores, temp);

      // StageI / start - find an acceptable starting configuration for subsequent ljrep optimization

      JumpManager jm;
      Mover mover(pose.get());

      for( int cycles_start = 0; cycles_start < max_cycles_start; ++cycles_start ) {

				cout << ">>> StageI / start - " << REPORTCYCLES1(cycles_start,max_cycles_start) << endl;

				(*pose) = (*pose_init);

				protocols::moves::MonteCarlo min_start_3mer(*pose,*score_fxn_lores,temp);

				cout << ">>> StageI / start / 3mer - inserting 3mers " << endl;

				for( int cycles_start_3mer = 0; cycles_start_3mer < max_cycles_start_3mer /* && get_acceptance_rate(min_start_3mer,min_acceptance_rate_start_3mer)*/ ; ++cycles_start_3mer ) {

					//	    if( getVerbose() ) {
					cout << "start_3mer: " << REPORTCYCLES2(cycles_start,max_cycles_start,cycles_start_3mer,max_cycles_start_3mer) << endl;
					//	    }

					(*pose) = (*pose_init);

					vector<Loop> loops_shuffled = loops; // need a list of all loops not just anchored loops
					random_shuffle(loops_shuffled.begin(),loops_shuffled.end());

					//g.setUseNb(false);

					// start with a completely random configuration
					for( Size k = 0; k < loops_shuffled.size(); ++k ) { // !!! first generate a completely random
						const Loop& loop = loops_shuffled[k];

						mover.randomFragments( loop.lo, loop.hi, 1 );
						mover.randomFragments( loop.lo, loop.hi, 3 );

						if( loop.to != 0 ) {

							mover.randomRotamer( loop.to );

							if( loop.tag->hasTag("hbond") ) {
								jm.set_random_hbond_jump( *pose, loop );
							}

						} // cycles_start_3mer_anchor

					}

					protocols::moves::MonteCarlo min_start_3mer_inner(*pose,*score_fxn_lores,temp);

					min_start_3mer_inner.boltzmann( *pose );
					min_start_3mer_inner.recover_low( *pose );

					for( Size k = 0; k < loops_shuffled.size(); ++k ) {

						const Loop& loop = loops_shuffled[k];
						int const a = loop.lo;
						int const b = loop.hi;

						if( loop.to != 0 ) {

							for( int cycles_start_3mer_anchor = 0; cycles_start_3mer_anchor < max_cycles_start_3mer_anchor; ++cycles_start_3mer_anchor ) {

								mover.randomRotamer( loop.to );
								min_start_3mer_inner.boltzmann( *pose );

							}

						}

						min_start_3mer_inner.recover_low( *pose );

						for( int cycles_start_3mer_1mer = 0; cycles_start_3mer_1mer < max_cycles_start_3mer_1mer; ++cycles_start_3mer_1mer ) {

							cout << "start_3mer_1mer: " << REPORTCYCLES3(cycles_start,max_cycles_start,cycles_start_3mer,max_cycles_start_3mer,cycles_start_3mer_1mer,max_cycles_start_3mer_1mer) << endl;

							mover.insertFragment(a,b,1);
							atm.run(*pose,move_map,*score_fxn_lores,core::optimization::MinimizerOptions("dfpmin",1,true,false,false) );

							min_start_3mer_inner.boltzmann(*pose);

							pose->energies().show_total_headers( cout ); cout << endl;
							pose->energies().show_totals( cout ); cout << endl;

						} // cycles_start_3mer_anchor

						min_start_3mer_inner.recover_low(*pose); // ???

					} // k

					min_start_3mer.boltzmann(*pose);

				} // cycles_start_3mer

				min_start_3mer.recover_low(*pose);

				// Stage I / start / 1mer - optimize the given conformation by 1mer & ccd moves
				cout << ">>> StageI / start / 1mer - inserting 1mers" << endl;

				protocols::moves::MonteCarlo min_start_1mer(*pose,*score_fxn_lores,temp);

				for( int cycles_start_1mer = 0; cycles_start_1mer < max_cycles_start_1mer /*&& get_acceptance_rate(min_start_1mer,DEFAULT_MIN_ACCEPTANCE_RATE)*/ ; ++cycles_start_1mer ) {

					cout << "Stage I / start / 1mer " << REPORTCYCLES2(cycles_start,max_cycles_start,cycles_start_1mer,max_cycles_start_1mer) << endl;

					pose->energies().show_total_headers( cout ); cout << endl;
					pose->energies().show_totals( cout ); cout << endl;

					const int k_rand_loop = rand() % loops.size(); // this is biased towards shorter loops

					const Loop& loop = loops[k_rand_loop];

					mover.smallMoves(loop.lo,loop.hi);
					//mover.insertFragment(a,b,1);

					atm.run(*pose,move_map,*score_fxn_lores,core::optimization::MinimizerOptions("dfpmin",0.01,true,false,false) );

					/*
						int
						fast_ccd_loop_closure(
						core::pose::Pose & pose,
						core::kinematics::MoveMap const & mm,
						int const loop_begin,
						int const loop_end,
						int const cutpoint,
						int const ii_cycles,
						core::Real const tolerance,
						bool const rama_check,
						core::Real const max_rama_score_increase,
						core::Real const max_total_delta_helix,
						core::Real const max_total_delta_strand,
						core::Real const max_total_delta_loop,
						core::Real & forward_deviation, // output
						core::Real & backward_deviation, // output
						core::Real & torsion_delta,
						core::Real & rama_delta
						);
					*/

					if( loop.to ) {
						core::Real fdev, rdev, tdelta, rdelta;
						protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.lo-1, loop.to-1, loop.lo, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
						//cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;
						protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.to+1, loop.hi+1, loop.hi, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
						//cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;

						(*score_fxn_lores)( *pose );
						pose->energies().show_total_headers( cout ); cout << endl;
						pose->energies().show_totals( cout ); cout << endl;

					}
					else {
						assert( false );
					}

					atm.run(*pose,move_map,*score_fxn_lores,core::optimization::MinimizerOptions("dfpmin",0.01,true,false,false) );

					pose->energies().show_total_headers( cout ); cout << endl;
					pose->energies().show_totals( cout ); cout << endl;

					min_start_1mer.boltzmann(*pose);
					//min_start_1mer.recover_low(*pose);

				} // cycles_start_1mer

				min_start_1mer.recover_low(*pose);
				min_start.boltzmann(*pose );

      } // cycles_start

      min_start.recover_low(*pose);

    }

    void Protocol::phase_hires() {

      // let's do this in a random order
      vector<Loop> loops_shuffled = loops;
      random_shuffle(loops_shuffled.begin(),loops_shuffled.end());

      Mover mover(pose.get());
      core::optimization::AtomTreeMinimizer atm;

      double temp = 1;

      const double MIN_WEIGHT_REP = 0.001;
      const double MAX_WEIGHT_REP = score_fxn_hires->get_weight( core::scoring::fa_rep );

      // LJRep Ramp up

      for( int cycles_ramp = 0; cycles_ramp <= max_cycles_ramp; ++cycles_ramp ) {

				cout << ">>> Stage II / start / ramp - " << REPORTCYCLES1(cycles_ramp,max_cycles_ramp) << endl;

				double f;
				if( cycles_ramp == 0 ) {
					f = MIN_WEIGHT_REP;
				}
				else {
					f = MAX_WEIGHT_REP*(static_cast<double>(cycles_ramp + 1) / max_cycles_ramp);
				}

				score_fxn_hires->set_weight( core::scoring::fa_rep, f );

				protocols::moves::MonteCarlo min_ramp(*pose, *score_fxn_hires, temp);

				min_ramp.boltzmann( *pose );

				for( Size k = 0; k < loops.size(); ++k ) {
					const Loop& loop = loops_shuffled[k];

					int const a = loop.lo;
					int const b = loop.hi;

					cout << "hi: " << a << " " << b << endl;

					for( int cycles_ramp_sm = 0; cycles_ramp_sm < max_cycles_ramp_sm; ++cycles_ramp_sm ) {

						mover.smallMoves(a,b);

						for( int cycles_ramp_sm_min = 0; cycles_ramp_sm_min < max_cycles_ramp_sm_min ; ++ cycles_ramp_sm_min ) {

							if( loop.to ) {
								core::Real fdev, rdev, tdelta, rdelta;
								protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.lo-1, loop.to-1, loop.lo, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
								cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;
								protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.to+1, loop.hi+1, loop.hi, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
								cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;

							}
							else {
								assert( false );
							}

							atm.run(*pose,move_map,*score_fxn_lores,core::optimization::MinimizerOptions("dfpmin",1,true,false,false) );

							min_ramp.boltzmann(*pose); // !!! don't recover here

						} // cycles_ramp_min

						min_ramp.boltzmann(*pose);

					} // cycles_ramp_small_moves

				} // cycles_ramp_loop

				min_ramp.recover_low(*pose);

			} // cycles_rep

      // stage 3 - design

      utility::vector1<bool> residues_allowed_to_be_packed( pose->n_residue(), false);

      core::pack::task::PackerTaskOP design_task( core::pack::task::TaskFactory::create_packer_task( *pose ) );
      //design_task->initialize_from_command_line(); // .read_resfile().or_include_current( true );

      for( Size k = 0; k < loops.size(); ++k ) {
				for( int i = loops[k].lo; i <= loops[k].hi ; ++i ) {
					if( i != loops[k].to ) {
						residues_allowed_to_be_packed[i] = true;
					}
				}
      }

      design_task->restrict_to_residues( residues_allowed_to_be_packed );

      for( int cycles_design = 0; cycles_design < max_cycles_design; ++cycles_design ) {

				cout << "Protocol::execute - DESIGN - starting cycle #" << cycles_design << " / " << max_cycles_design << endl;

				cout << "calling pack_rotamers..." << endl;
				core::pack::pack_rotamers( *pose, *score_fxn_hires, design_task );
				cout << "done calling pack_rotamers" << endl;

				protocols::moves::MonteCarlo min_design_sm(*pose, *score_fxn_hires, temp);

				for( Size k = 0; k < loops_shuffled.size(); ++k ) {

					cout << "Protocol::execute - DESIGN - small moves + min after design cycle# " << cycles_design << " small moves + insertions into loop#" << k << endl;

					Loop& loop = loops_shuffled[k];

					int const a = loop.lo;
					int const b = loop.hi;

					protocols::moves::MonteCarlo min_design_sm(*pose,*score_fxn_hires,temp);

					for( int cycles_design_sm = 0; cycles_design_sm < max_cycles_design_sm; ++cycles_design_sm ) {

						mover.smallMoves(a,b);

						for( int cycles_design_sm_min = 0; cycles_design_sm_min < max_cycles_design_sm_min; ++cycles_design_sm_min ) {

							if( loop.to ) {
								//cout << "calling ccd on: " << a << " " << loop.to << " " << b << endl;

								core::Real fdev, rdev, tdelta, rdelta;
								protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.lo-1, loop.to-1, loop.lo, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
								//cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;
								protocols::loops::fast_ccd_loop_closure(*pose,move_map, loop.to+1, loop.hi+1, loop.hi, 100, 0.1, true, 5.0, 5.0, 5.0, 5.0, fdev, rdev, tdelta, rdelta );
								//cout << "fdev=" << fdev << " rdev=" << rdev << " tdelta=" << tdelta << " rdelta=" << rdelta << endl;

							}
							else {
								assert( false );
							}

							atm.run(*pose,move_map,*score_fxn_lores,core::optimization::MinimizerOptions("dfpmin",1,true,false,false) );

							min_design_sm.boltzmann(*pose);

						} // cycles_design_sm_min

					} // cycles_design_sm

				} // k

				min_design_sm.recover_low(*pose);

      } // cycles_design

    }


    void Protocol::apply(core::pose::PoseOP pose_) {
      pose = pose_;
      move_map = get_move_map(pose->n_residue(),loops);

      const double CHAINBREAK_WEIGHT = 10;
      const double RAMA_WEIGHT       = 10;
      const double DUNBRACK_WEIGHT   = 1;
      const double LORES_REP_WEIGHT  = 0.001;

      score_fxn_lores = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
      score_fxn_lores->set_weight( core::scoring::chainbreak, CHAINBREAK_WEIGHT );
      score_fxn_lores->set_weight( core::scoring::rama,       RAMA_WEIGHT       );
      score_fxn_lores->set_weight( core::scoring::fa_dun,     DUNBRACK_WEIGHT   );
      score_fxn_lores->set_weight( core::scoring::fa_rep,     LORES_REP_WEIGHT  );

      //score_fxn_lores = core::scoring::ScoreFunctionFactory::create_score_function( "murphp_lores" );
      //score_fxn_lores = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::SOFT_REP_WTS );

      score_fxn_hires = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS );
      score_fxn_hires->set_weight( core::scoring::chainbreak, CHAINBREAK_WEIGHT );
      score_fxn_hires->set_weight( core::scoring::rama,       RAMA_WEIGHT       );

      //for( int i = 0; i <

      cout << "Protocol::apply: running lo_res phase" << endl;

      phase_lores();

      cout << "Protocol::apply: running hi_res_phase" << endl;

      phase_hires();


      // Score new structure.  Cached energies (including *residue* energies)
      // must be up-to-date in order to get sensible output.  If you remove these
      // lines, you *must* insert equivalent logic at the end of all apply() methods
      // (or at least for all movers that might be passed to this function).
      (*score_fxn_hires)( *pose );
      //score_fxn_hires->accumulate_residue_total_energies( *pose );


      /*

      core::scoring::Energies const& energies = pose->energies();

      cout << "Protocol::apply::getting score_types:" << endl;

      core::scoring::ScoreTypes score_types = score_fxn_hires->get_nonzero_weighted_scoretypes();


      cout << "Protocol::apply::ScoreTypes: " << score_types.size() << ": ";
      for( Size i = 1; i <= score_types.size(); ++i ) {
			cout << score_types[i] << ",";
      }
      cout << endl;


      cout << "Protocol::apply::Energies:" << endl;

      for( Size i = 1; i <= pose->n_residue(); ++i ) {

			//core::scoring::EnergyMap emap = energies();

			core::scoring::EMapVector const& emap_1b = energies.onebody_energies(i);

			cout << "i " << i << " emap_1b:\t";
			for( Size k = 1; k <= score_types.size(); ++k ) {
			if( score_types[k] == core::scoring::chainbreak ) {
	    cout << "cb: ";
	    cout << emap_1b[ score_types[k] ] << "\t";
			}
			}
			cout << " ";



			core::scoring::EMapVector const& emap_2b = energies.residue_total_energies(i);
			cout << " emap_2b:\t";
			for( Size k = 1; k <= score_types.size(); ++k ) {
			if( score_types[k] == core::scoring::chainbreak ) {
	    cout << emap_2b[ score_types[k] ] << "\t";
			}
			}

			cout << endl;

      }


      for( Size k = 1; k <= score_types.size(); ++k ) {
			cout << score_types[k] << " " << energies.total_energies()[ score_types[k] ] << endl;
      }

      //       for( int i = 1; i <= pose->n_residues(); ++i ) {
      //	for( int j = 1; j <= pose->n_residues(); ++j ) {
      //	}
      //       }

      //cout << "done with pose->energies" << endl;

      */

    }


  }

}
