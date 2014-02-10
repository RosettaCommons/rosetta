from __jd2_all_at_once_ import *

import rosetta.utility.tag

from types import MethodType

def parse_rosettascript(self, protocol_string, starting_pose = None):
    """Parses a rosetta script string into a wrapper executor mover."""
    tag = rosetta.utility.tag.Tag.create(protocol_string)

    if starting_pose is None:
        protocol = self.parse_protocol_tag(tag)
    else:
        protocol = self.parse_protocol_tag(starting_pose, tag)

    return RosettaScriptWrapper(protocol)

# temporary commenting out due to change in upstream
#DockDesignParser.parse_rosettascript = MethodType(parse_rosettascript, None, DockDesignParser)

class RosettaScriptWrapper(object):
    def __init__(self, parsed_protocol):
        self.parsed_protocol = parsed_protocol

    def apply(self, pose, include_partial_result = True):
        """Runs the stored protocol and returns protocol results stored in jd2.

        If include_partial_result is True return partial application results even
        if non MS_SUCCESS mover status is returned.

        Destructively interacts with jd2."""

        #Reset jd2 to clear job data.
        JobDistributor.get_instance().restart()

        self.parsed_protocol.apply(pose)
        status = self.parsed_protocol.get_last_move_status()

        if not include_partial_result and status != rosetta.protocols.moves.MoverStatus.MS_SUCCESS:
            return None
        else:
            return dict(
                JobDistributor.get_instance().current_job().get_string_real_pairs().items() +
                JobDistributor.get_instance().current_job().get_string_string_pairs().items())

    def apply_to_failure(self, pose, limit=10000):
        """Interates the stored protocol until FAIL_DO_NOT_RETRY is returned."""
        run_count = 0

        while (limit is None) or (run_count < limit):
            run_count += 1

            results = self.apply(pose)
            status = self.parsed_protocol.get_last_move_status()

            if status == rosetta.protocols.moves.MoverStatus.MS_SUCCESS:
                yield (results, pose)
            elif status == rosetta.protocols.moves.MoverStatus.FAIL_RETRY:
                continue
            else:
                break
