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

# Create Custom Metric ----------------------------------------------------------------------------------------------------------
from pyrosetta.rosetta.core import simple_metrics
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.std import list_utility_tag_XMLSchemaAttribute_t
from pyrosetta.rosetta.std import map_unsigned_long_double

class PerResidueBootCampMetric(simple_metrics.PerResidueRealMetric):
    """Python PerResidueBfactorBootCampMetric simple metric"""
    clones_ = list()

    def __init__(self, num: float = 10.5):
        simple_metrics.PerResidueRealMetric.__init__(self)
        self.num_: float = num

    # Name of the class
    def name(self) -> str:
        return self.class_name()

    # Name of the metric
    def metric(self) -> str:
        return "bootcamp"

    @staticmethod
    def class_name() -> str:
        return "PerResidueBootCampMetric"

    def parse_my_tag(self, tag, datamap) -> None:
        self.num_ = tag.get_option_string("num", "10.5")

    @classmethod
    def provide_xml_schema(cls, xsd) -> None:
        from pyrosetta.rosetta.utility.tag import XMLSchemaAttribute, XMLSchemaType
        from pyrosetta.rosetta.utility.tag import xsct_real

        attrlist = list_utility_tag_XMLSchemaAttribute_t()

        attrlist.append(XMLSchemaAttribute.required_attribute(
            "num",
            XMLSchemaType(xsct_real),
            "The number to report as your metric"))

        description = '''
            PerResidue simple metric that reports 10.5 for every residue
            '''
        simple_metrics.xsd_simple_metric_type_definition_w_attributes(
                xsd,
                cls.class_name(),
                description, attrlist)

    def calculate(self, pose) -> map_unsigned_long_double:
        residue_selector = self.get_selector()
        subset = residue_selector.apply(pose)
        selection = select.get_residues_from_subset(subset)

        metric_map = map_unsigned_long_double()

        for resi in selection:
            rt = pose.residue_type(resi)

            metric_map[ int(resi) ] = float(self.num_)
        return metric_map


# Register Metric ---------------------------------------------------------------------------------------------------------------
from pyrosetta.rosetta.core import simple_metrics

#global var
_py_metric_creators_ = []

class PerResidueBootCampMetricCreator(simple_metrics.SimpleMetricCreator):
    instances_ = list()

    def __init__(self):
        simple_metrics.SimpleMetricCreator.__init__(self)

    def create_simple_metric(self):
        metric = PerResidueBootCampMetric()
        self.instances_.append(metric)
        return metric

    def keyname(self):
        return PerResidueBootCampMetric.class_name()

    def provide_xml_schema(self, xsd):
        PerResidueBootCampMetric.provide_xml_schema(xsd)

def register():

    factory = simple_metrics.SimpleMetricFactory.get_instance()
    creator = PerResidueBootCampMetricCreator()
    factory.factory_register(creator)

    _py_metric_creators_.append(creator)


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
  <SIMPLE_METRICS>
    <PerResidueBootCampMetric name="bootcamp" num="13.7"/>
  </SIMPLE_METRICS>
  <MOVERS>
    <RunSimpleMetrics name="metrics" metrics="bootcamp"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="metrics"/>
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
