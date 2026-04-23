# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from utils.distributed import parse_args, run_distributed_cluster_test_cases


def main(wait: bool, streaming: bool, timeout: int) -> None:
    run_distributed_cluster_test_cases(
        "test_runtime.RuntimeTest.test_timing_multi_instance",
        "test_runtime.RuntimeTest.test_timing_single_instance",
        wait=wait,
        streaming=streaming,
        timeout=timeout,
    )

if __name__ == "__main__":
    args = parse_args()
    main(args.wait, args.streaming, args.timeout)
