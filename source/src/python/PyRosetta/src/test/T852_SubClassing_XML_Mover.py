#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Hope Woods
## @author Sergey Lyskov

## @brief  Demo for using PyRosetta Movers in XML scripts
##         Based on Hope Wood work for PyRosetta BootCamp

from __future__ import print_function
from __future__ import annotations


import sys

import pyrosetta
import pyrosetta.rosetta as rosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')



# Create Custom Mover ----------------------------------------------------------------------------------------------------------
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta import core, protocols, basic, utility
from pyrosetta.rosetta.protocols.moves import Mover, MoverCreator, MoverFactory
from pyrosetta.rosetta.utility.tag import XMLSchemaDefinition, xsct_positive_integer, XMLSchemaComplexTypeGenerator, XMLSchemaAttribute, XMLSchemaType
from pyrosetta.rosetta.std import list_utility_tag_XMLSchemaAttribute_t
from pyrosetta.rosetta.core.scoring import ScoreFunction, ScoreFunctionFactory, attributes_for_parse_score_function_w_description, parse_score_function
from pyrosetta.rosetta.protocols import moves as protocols_moves
from typing import Optional, Union


class BootCampMoverPy(Mover):
    """Python version of BootCampMover with getters/setters, parse_my_tag, and schema."""

    clones_ = list()

    def __init__(self, num_iterations: int = 500):

        Mover.__init__(self)

        self.num_iterations_: int = int(num_iterations)

    # --------- identity / cloning ---------
    @staticmethod
    def mover_name() -> str:
        return "BootCampMoverPy"

    def get_name(self) -> str:
        name = BootCampMoverPy.mover_name()
        return name

    def clone(self):
        copy = BootCampMoverPy(self.num_iterations_)
        BootCampMoverPy.clones_.append(copy)
        return copy

    def fresh_instance(self):
        return BootCampMoverPy()

    # --------- getters / setters ----------
    def get_num_iterations(self):
        return self.num_iterations_

    def set_num_iterations(self, iters: int):
        iters = int(iters)
        assert iters > 0
        self.num_iterations_ = iters

    # --------- RosettaScripts parse--------
    def parse_my_tag(self, tag, datamap):
        if tag.hasOption("num_iterations"):
            iters = tag.get_option_int("num_iterations", 1)
            self.set_num_iterations(iters)

    @staticmethod
    def provide_xml_schema(xsd):
        attrs = list_utility_tag_XMLSchemaAttribute_t()

        num_iterations_attr = XMLSchemaAttribute.attribute_w_default("num_iterations",
                                               XMLSchemaType(xsct_positive_integer),
                                               "Number of iterations for loop", "10")
        attrs.append(num_iterations_attr)
        attributes_for_parse_score_function_w_description(attrs, "ScoreFunction to use for sampling")

        protocols.moves.xsd_type_definition_w_attributes(
            xsd, "BootCampMoverPy",
            "BootCampMover: for an example",
            attrs
        )

    # --------------- apply ----------------
    def apply(self, pose):

        for i in range(self.num_iterations_):
            print("Hello World from PyRosetta XML Mover!")


# Register Mover ---------------------------------------------------------------------------------------------------------------
from pyrosetta.rosetta import protocols

#global var
_py_mover_creators_ = []

class BootCampMoverCreatorPy(protocols.moves.MoverCreator):
    instances_ = list()

    def __init__(self):
        protocols.moves.MoverCreator.__init__(self)

    def create_mover(self):
        mover = BootCampMoverPy()
        self.instances_.append(mover)
        return mover

    def keyname(self):
        return BootCampMoverPy.mover_name()

    def provide_xml_schema(self, xsd):
        print("creator provide_xml_schema is called")
        BootCampMoverPy.provide_xml_schema(xsd)

def register():
    factory = protocols.moves.MoverFactory.get_instance()
    creator = BootCampMoverCreatorPy()
    factory.factory_register(creator)

    _py_mover_creators_.append(creator)
    print("registration is happening")


# Main -------------------------------------------------------------------------------------------------------------------------
from pyrosetta import init, pose_from_sequence
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta import protocols, utility
from pyrosetta.rosetta.utility import vector1_std_string
from pyrosetta.rosetta.utility.options import OptionCollection
from pyrosetta.rosetta.basic.datacache import DataMap

# Embedded RosettaScripts XML:
EMBEDDED_XML = """<ROSETTASCRIPTS>
  <SCOREFXNS>
  </SCOREFXNS>
  <MOVERS>
    <!-- The tag name MUST match your MoverCreator.keyname() -->
    <BootCampMoverPy name="bcm" num_iterations="4"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="bcm"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
"""



# Register your mover BEFORE parsing XML so schema/type are known
register()
init()

pose = pose_from_sequence("ARNDCEQGHILKMFPSTWYV")


# Parse XML string and build the overall protocol mover
xmlobj = XmlObjects.create_from_string(EMBEDDED_XML)
protocol = xmlobj.get_mover("ParsedProtocol")  # the top-level protocol
protocol.apply(pose)
