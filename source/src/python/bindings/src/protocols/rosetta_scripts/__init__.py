from __rosetta_scripts_all_at_once_ import *

from types import MethodType

def parse_rosettascript(self, protocol_string, starting_pose = None):
    """Parses a rosetta script string into a wrapper executor mover."""
    tag = rosetta.utility.tag.Tag.create(protocol_string)

    if starting_pose is None:
        protocol = self.parse_protocol_tag(tag)
    else:
        protocol = self.parse_protocol_tag(starting_pose, tag)

    return RosettaScriptWrapper(protocol)

RosettaScriptsParser.parse_rosettascript = MethodType(parse_rosettascript, None, RosettaScriptsParser)
