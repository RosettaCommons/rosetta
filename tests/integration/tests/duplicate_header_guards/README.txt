This test checks for duplicate header guards, which are usually the result of copy/pasting.
For example, the following duplicates existed at the time this test was written.

$ grep '#define INCLUDED' -r src | awk -F: '{print $2}' | sort | uniq -c | sort -nk1 | awk 'int($1)>=2'
      2 #define INCLUDED_basic_gpu_GPU_cc
      2 #define INCLUDED_core_pack_annealer_FixbbSimAnnealer_fwd_hh
      2 #define INCLUDED_core_pack_annealer_FixbbSimAnnealer_hh
      2 #define INCLUDED_core_pack_interaction_graph_MultiplexedAnnealableGraph_hh
      2 #define INCLUDED_core_pack_min_pack_hh
      2 #define INCLUDED_core_pose_Pose_hh
      2 #define INCLUDED_core_scoring_hbonds_HBondEnergy_fwd_hh
      2 #define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc
      2 #define INCLUDED_core_scoring_methods_CenPairMotifEnergy_hh
      2 #define INCLUDED_core_select_residue_selector_GlycanResidueSelector_fwd_hh
      2 #define INCLUDED_core_sequence_SequenceProfile_hh
      2 #define INCLUDED_protocols_cryst_spacegroup_hh
      2 #define INCLUDED_protocols_forge_constraints_InverseRotamersCstGeneratorCreator_hh
      2 #define INCLUDED_protocols_forge_constraints_InvrotTreeCstGeneratorCreator_hh
      2 #define INCLUDED_protocols_moves_ResidueVicinityCstGeneratorCreator_hh
      2 #define INCLUDED_protocols_STEPWISE_Clusterer_FWD_HH
      2 #define INCLUDED_utility_io_irstream_hh
      3 #define INCLUDED_protocols_simple_filters_ReportFilter_hh





Breaking up the command into digestible fragments:

Find all header guards:
$ grep '#define INCLUDED' -r src | head
src/devel/replica_docking/ModulatedMoverCreator.hh:#define INCLUDED_devel_replica_docking_ModulatedMoverCreator_hh
src/devel/replica_docking/TempInterpolatorFactory.hh:#define INCLUDED_devel_replica_docking_TempInterpolatorFactory_hh
src/devel/replica_docking/UnbiasedRigidBodyMover.fwd.hh:#define INCLUDED_devel_replica_docking_UnbiasedRigidBodyMover_fwd_hh
src/devel/replica_docking/IrmsdFilter.hh:#define INCLUDED_devel_replica_docking_IrmsdFilter_hh
src/devel/replica_docking/CaIrmsdFilter.hh:#define INCLUDED_devel_replica_docking_CaIrmsdFilter_hh
src/devel/replica_docking/CaIrmsdFilterCreator.hh:#define INCLUDED_devel_replica_docking_CaIrmsdFilterCreator_hh
src/devel/replica_docking/UnbiasedRigidBodyMoverCreator.hh:#define INCLUDED_devel_replica_docking_UnbiasedRigidBodyMoverCreator_hh
src/devel/replica_docking/TempWeightedMetropolisHastingsMoverCreator.hh:#define INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMoverCreator_hh
src/devel/replica_docking/InteractionScoreFilter.hh:#define INCLUDED_devel_replica_docking_InteractionScoreFilter_hh
src/devel/replica_docking/AddEncounterConstraintMover.hh:#define INCLUDED_devel_replica_docking_AddEncounterConstraintMover_hh

Remove filenames from lines:
$ grep '#define INCLUDED' -r src | awk -F: '{print $2}' | head
#define INCLUDED_devel_replica_docking_ModulatedMoverCreator_hh
#define INCLUDED_devel_replica_docking_TempInterpolatorFactory_hh
#define INCLUDED_devel_replica_docking_UnbiasedRigidBodyMover_fwd_hh
#define INCLUDED_devel_replica_docking_IrmsdFilter_hh
#define INCLUDED_devel_replica_docking_CaIrmsdFilter_hh
#define INCLUDED_devel_replica_docking_CaIrmsdFilterCreator_hh
#define INCLUDED_devel_replica_docking_UnbiasedRigidBodyMoverCreator_hh
#define INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMoverCreator_hh
#define INCLUDED_devel_replica_docking_InteractionScoreFilter_hh
#define INCLUDED_devel_replica_docking_AddEncounterConstraintMover_hh

Count the number of each instance:
(requires sorting to ensure that identical items are binned together)
$ grep '#define INCLUDED' -r src | awk -F: '{print $2}' | sort | uniq -c | head
      1 #define INCLUDED_apps_benchmark_FastRelax_bench_hh
      1 #define INCLUDED_apps_benchmark_init_util_hh
      1 #define INCLUDED_apps_benchmark_InteractionGraph_bench_hh
      1 #define INCLUDED_apps_benchmark_performance_benchmark_hh
      1 #define INCLUDED_apps_benchmark_ScoreEach_bench_hh
      1 #define INCLUDED_apps_pilot_blivens_disulfides_hh
      1 #define INCLUDED_apps_pilot_christophe_PCS_main_hh
      1 #define INCLUDED_apps_pilot_dgront_JumpSpecificAbrelax_hh
      1 #define INCLUDED_apps_pilot_hpark_sampling_movers_hh
      1 #define INCLUDED_apps_pilot_phil_phil_HH

