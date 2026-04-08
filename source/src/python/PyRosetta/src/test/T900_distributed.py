# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases, run_distributed_cluster_test_cases


run_test_cases(
    "pyrosetta.tests.bindings.init.test_init_files",
    "pyrosetta.tests.bindings.core.test_pose",
    "pyrosetta.tests.distributed.test_concurrency",
    "pyrosetta.tests.distributed.test_dask",
    "pyrosetta.tests.distributed.test_gil",
    "pyrosetta.tests.distributed.test_smoke",
    "pyrosetta.tests.distributed.test_viewer",
    "pyrosetta.tests.numeric.test_alignment",
)
run_distributed_cluster_test_cases(
    "test_smoke.SmokeTest.test_smoke",
    "test_smoke.SmokeTest.test_invalid_tasks",
    "test_smoke.SmokeTest.test_ignore_errors",
    "test_smoke_multi.SmokeTestMulti.test_smoke_multi",
    "test_smoke_multi.SmokeTestMulti.test_smoke_multi_from_instance",
    "test_save_all.SaveAllTest.test_save_all",
    "test_save_all.SaveAllTest.test_save_all_dry_run",
    "test_serialization.SerializationTest.test_serialization",
    "test_scores.ScoresTest.test_detached_scores",
    "test_scores.ScoresTest.test_detached_scores_in_protocol",
    "test_scores.ScoresTest.test_detached_scores_with_reserve_scores",
    "test_scores.ScoresTest.test_secure_packages_billiard",
    "test_clients.MultipleClientsTest.test_clients",
    "test_resources.ResourcesTest.test_resources",
    "test_resources.ResourcesTest.test_resources_clients",
    "test_generator.GeneratorTest.test_generate_builtin_clients",
    "test_generator.GeneratorTest.test_generate_multi_user_clients",
    "test_generator.GeneratorTest.test_generate_partition_clients",
    "test_generator.GeneratorTest.test_generate_user_client",
    "test_io.IOTest.test_io",
    "test_init_file_signer.TestInitFileSigner.test_init_file_signer",
    "test_priorities.PrioritiesTest.test_priorities",
    "test_retries.RetriesTest.test_retries_persistent_errors",
    "test_retries.RetriesTest.test_retries_succeed_on_last_retry",
    "test_retries.RetriesTest.test_no_retries",
    "test_retries.RetriesTest.test_retries_api",
    "test_task_keys.TaskKeysTest.test_task_keys",
    "test_worker_preemption.WorkerPreemptionTest.test_disk_task_registry",
    "test_worker_preemption.WorkerPreemptionTest.test_memory_task_registry",
    "test_worker_preemption.WorkerPreemptionTest.test_disk_max_task_replicas_all",
    "test_worker_preemption.WorkerPreemptionTest.test_disk_max_task_replicas_int",
    "test_logging.LoggingTest.test_logging",
    "test_reproducibility.TestReproducibility.test_reproducibility_minimizer_nstruct",
    "test_reproducibility.TestReproducibility.test_reproducibility_minimizer_nstruct_filter_results",
    "test_reproducibility.TestReproducibility.test_reproducibility_packer_nstruct",
    "test_reproducibility.TestReproducibility.test_reproducibility_packer_nstruct_filter_results",
    "test_reproducibility.TestReproducibility.test_reproducibility_packer_separate",
    "test_reproducibility.TestReproducibility.test_reproducibility_packer_separate_filter_results",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_filter_results",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_decoy_ids",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_decoy_ids_filter_results",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_from_reproduce",
    "test_reproducibility_multi.TestReproducibilityMulti.test_reproducibility_from_reproduce_filter_results",
    "test_reproducibility_pose_dataframe.TestReproducibilityPoseDataFrame.test_reproducibility_from_reproduce",
    "test_reproducibility_task_updates.TestReproducibilityTaskUpdates.test_reproduce_task_updates",
    "test_reproducibility_task_updates.TestReproducibilityTaskUpdates.test_reproduce_task_updates_norm_task_options",
    "test_reproducibility_task_updates.TestReproducibilityTaskUpdates.test_reproduce_task_updates_with_init_file",
    "test_reproducibility_task_updates.TestReproducibilityTaskUpdates.test_reproduce_task_updates_norm_task_options_with_init_file",
    "test_reproducibility_remodel_task_updates.TestReproducibilityRemodelTaskUpdates.test_reproduce_remodel_task_updates",
    "test_reproducibility_remodel_task_updates.TestReproducibilityRemodelTaskUpdates.test_reproduce_remodel_task_updates_norm_task_options",
    "test_reproducibility_remodel_task_updates.TestReproducibilityRemodelTaskUpdates.test_reproduce_remodel_task_updates_with_init_file",
    "test_reproducibility_remodel_task_updates.TestReproducibilityRemodelTaskUpdates.test_reproduce_remodel_task_updates_norm_task_options_with_init_file",
    "test_runtime.RuntimeTest.test_timing_multi_instance",
    "test_runtime.RuntimeTest.test_timing_single_instance",
)
