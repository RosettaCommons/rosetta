import threading

import traitlets

import pyrosetta
import pyrosetta.rosetta.basic.options
import pyrosetta.rosetta.protocols.rosetta_scripts as rosetta_scripts
import pyrosetta.rosetta.protocols.moves as moves

import pyrosetta.distributed
import pyrosetta.distributed.tasks.taskbase as taskbase
import pyrosetta.distributed.packed_pose as packed_pose


def validate(protocol_xml):
    """Perform schema and parse validation for the given protocol xml."""
    try:
        test_task = BaseRosettaScriptsTask(protocol_xml)
        test_task.maybe_setup()
    except RuntimeError as error:
        raise error


class BaseRosettaScriptsTask(taskbase.TaskBase):
    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def parser(self):
        if not getattr(self, "_parser", None):
            BaseRosettaScriptsTask._parser = \
                    rosetta_scripts.RosettaScriptsParser()

        return self._parser

    protocol_xml = traitlets.CUnicode()

    def __init__(self, protocol_xml):
        super().__init__(protocol_xml=protocol_xml)

    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def setup(self):
        self.default_options = pyrosetta.rosetta.basic.options.process()
        self.tag = self.parser.create_tag_from_xml_string(
                self.protocol_xml, self.default_options)

        # Validate by parsing
        self.parser.parse_protocol_tag(self.tag, self.default_options)
        self.protocol_lock = threading.Lock()

    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def parsed_protocol(self):
        return self.parser.parse_protocol_tag(self.tag, self.default_options)

    def execute(self, pack_or_pose):
        return packed_pose.to_packed(self.apply(pack_or_pose))


class MultioutputRosettaScriptsTask(BaseRosettaScriptsTask):
    @pyrosetta.distributed.requires_init
    def apply(self, pack_or_pose):
        """Apply task generating pose objects."""
        protocol = self.parsed_protocol

        wpose = packed_pose.to_pose(pack_or_pose)

        with self.protocol_lock:
            protocol.apply(wpose)

            if protocol.get_last_move_status() != moves.MoverStatus.MS_SUCCESS:
                return

            while wpose:
                yield wpose
                wpose = protocol.get_additional_output()


class SingleoutputRosettaScriptsTask(BaseRosettaScriptsTask):
    @pyrosetta.distributed.requires_init
    def apply(self, pack_or_pose):
        """Apply task returning a pose object."""
        protocol = self.parsed_protocol

        wpose = packed_pose.to_pose(pack_or_pose)

        with self.protocol_lock:
            protocol.apply(wpose)

            if protocol.get_last_move_status() != moves.MoverStatus.MS_SUCCESS:
                return
            else:
                return wpose
