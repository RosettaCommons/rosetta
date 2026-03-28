# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import pyrosetta
import unittest

from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector,
    ChainSelector,
    NotResidueSelector,
    OrResidueSelector,
    ResidueIndexSelector,
    TrueResidueSelector,
)
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects


pyrosetta.init(extra_options="-constant_seed", set_logging_handler="logging")

class TestResidueSelectors(unittest.TestCase):
    def test_residue_selectors(self):
        pose = pyrosetta.pose_from_sequence("TESTING/MANY/CHAINS")

        # Define default residue selectors
        chA = ChainSelector("A")
        chB = ChainSelector("B")
        chC = ChainSelector("C")

        # Test AndResidueSelector
        chA_and_chB = AndResidueSelector(chA, chB)
        chA_and_chB_test = chA & chB
        self.assertListEqual(chA_and_chB.get_residues(pose), chA_and_chB_test.get_residues(pose))
        self.assertListEqual([], chA_and_chB_test.get_residues(pose))

        chA_or_chB = OrResidueSelector(chA, chB)
        chB_or_chC = OrResidueSelector(chB, chC)
        chB = AndResidueSelector(chA_or_chB, chB_or_chC)
        chB_test = chA_or_chB & chB_or_chC
        self.assertListEqual(chB.get_residues(pose), chB_test.get_residues(pose))
        self.assertNotEqual([], chB_test.get_residues(pose))

        # Test OrResidueSelector
        chA_or_chC = OrResidueSelector(chA, chC)
        chA_or_chC_test = chA | chC
        self.assertListEqual(chA_or_chC.get_residues(pose), chA_or_chC_test.get_residues(pose))

        chA_or_chB_or_chC = OrResidueSelector()
        chA_or_chB_or_chC.add_residue_selector(chA)
        chA_or_chB_or_chC.add_residue_selector(chB)
        chA_or_chB_or_chC.add_residue_selector(chC)
        chA_or_chB_or_chC_test = chA | chB | chC
        self.assertListEqual(chA_or_chB_or_chC.get_residues(pose), chA_or_chB_or_chC_test.get_residues(pose))
        self.assertListEqual(TrueResidueSelector().get_residues(pose), chA_or_chB_or_chC_test.get_residues(pose))

        # Test NotResidueSelector
        not_chA = NotResidueSelector(chA)
        not_chA_test = ~chA
        self.assertListEqual(not_chA.get_residues(pose), not_chA_test.get_residues(pose))

        # Test exclusive OrResidueSelector
        chAB = ChainSelector("A,B")
        chBC = ChainSelector("B,C")
        chAB_xor_chBC = AndResidueSelector(OrResidueSelector(chAB, chBC), NotResidueSelector(AndResidueSelector(chAB, chBC)))
        chA_or_chC = OrResidueSelector(chA, chC)
        chAB_xor_chBC_test = chAB ^ chBC
        self.assertListEqual(chAB_xor_chBC.get_residues(pose), chA_or_chC.get_residues(pose))
        self.assertListEqual(chAB_xor_chBC.get_residues(pose), chAB_xor_chBC_test.get_residues(pose))

        selector_1 = ResidueIndexSelector("2,4,6,8,10")
        selector_2 = ResidueIndexSelector("1,4,6,8,9")
        xor_selector_test = selector_1 ^ selector_2
        self.assertListEqual([1, 2, 9, 10], xor_selector_test.get_residues(pose))

        # Test selection logic using python operators
        selector_1 = ResidueIndexSelector("2,4,6,8,10")
        selector_2 = ResidueIndexSelector("1,4,6,8,9")
        selector_3 = ResidueIndexSelector("1,3,5,7,9")
        selector_4 = ResidueIndexSelector("3,8,15")
        logic_selector = OrResidueSelector(
            AndResidueSelector(selector_1, NotResidueSelector(selector_2)),
            AndResidueSelector(selector_3, NotResidueSelector(OrResidueSelector(selector_2, selector_4))),
        )
        logic_selector_test = (selector_1 & ~selector_2) | (selector_3 & ~(selector_2 | selector_4))
        self.assertListEqual(logic_selector.get_residues(pose), logic_selector_test.get_residues(pose))

        # Test compatibility with LogicResidueSelector in XML
        xml_str = """
        <RESIDUE_SELECTORS>
            <Index name="selector_1" resnums="2,4,6,8,10"/>
            <Index name="selector_2" resnums="1,4,6,8,9"/>
            <Index name="selector_3" resnums="1,3,5,7,9"/>
            <Index name="selector_4" resnums="3,8,15"/>
            <Logic name="logic_selector" selector="(selector_1 AND (NOT selector_2)) OR (selector_3 AND (NOT (selector_2 OR selector_4)))"/>
        </RESIDUE_SELECTORS>
        """
        xml_obj = XmlObjects.create_from_string(xml_str)
        selector_1 = xml_obj.get_residue_selector("selector_1")
        selector_2 = xml_obj.get_residue_selector("selector_2")
        selector_3 = xml_obj.get_residue_selector("selector_3")
        selector_4 = xml_obj.get_residue_selector("selector_4")
        logic_selector = xml_obj.get_residue_selector("logic_selector")
        logic_selector_test = (selector_1 & ~selector_2) | (selector_3 & ~(selector_2 | selector_4))
        self.assertListEqual(logic_selector.get_residues(pose), logic_selector_test.get_residues(pose))

        # Test in-place AndResidueSelector
        chA_or_chB = OrResidueSelector(chA, chB)
        chB_or_chC = OrResidueSelector(chB, chC)
        chB = AndResidueSelector(chA_or_chB, chB_or_chC)
        chB_test = chA_or_chB
        chB_test &= chB_or_chC
        self.assertListEqual(chB.get_residues(pose), chB_test.get_residues(pose))

        # Test in-place OrResidueSelector
        chA_or_chC = OrResidueSelector(chA, chC)
        chA_or_chC_test = chA
        chA_or_chC_test |= chC
        self.assertListEqual(chA_or_chC.get_residues(pose), chA_or_chC_test.get_residues(pose))

        # Test in-place exclusive OrResidueSelector
        selector_1 = ResidueIndexSelector("2,4,6,8,10")
        selector_2 = ResidueIndexSelector("1,4,6,8,9")
        xor_selector_test = selector_1
        xor_selector_test ^= selector_2
        self.assertListEqual([1, 2, 9, 10], xor_selector_test.get_residues(pose))

    def test_residue_selector_ctors(self):
        pose = pyrosetta.pose_from_sequence("TESTING/MANY/CHAINS")
        #Not
        primary_selector = ResidueIndexSelector("2,4,6")
        not_primary_selector = NotResidueSelector(primary_selector)
        # Case 1: Pass NotSelector to NotSelector
        selector_1 = NotResidueSelector(not_primary_selector)
        # Case 2: Pass AndSelector to NotSelector
        selector_2 = NotResidueSelector(AndResidueSelector(not_primary_selector,TrueResidueSelector()))
        self.assertListEqual(selector_1.get_residues(pose), selector_2.get_residues(pose))

if __name__ == "__main__":
    unittest.main()
