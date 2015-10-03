// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/gasteiger/GasteigerAtomTyper.cxxtest.hh
/// @brief  test suite for core::chemical::gasteiger::GasteigerAtomTyper.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/residue_io.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// Project headers
#include <test/core/init_util.hh>

using namespace core;
using namespace core::chemical;
using namespace core::chemical::gasteiger;

static basic::Tracer TR("core.chemical.gasteiger.GasteigerAtomTyper.cxxtest");

template< class Graph >
void dump_resgraph( Graph const & graph ) {
	typedef typename Graph::vertex_iterator Viter;
	typedef typename Graph::out_edge_iterator OEiter;
	Viter aitr,aitr_end;
	for ( boost::tie( aitr, aitr_end ) = boost::vertices( graph ); aitr != aitr_end; ++aitr ) {
		TR << " Atom " << graph[ *aitr ].name() << " -- ";
		OEiter eitr, eitr_end;
		for ( boost::tie( eitr, eitr_end ) = boost::out_edges( *aitr, graph ); eitr !=  eitr_end; ++eitr_end ) {
			TR << " " << graph[ boost::target( *eitr, graph ) ].name();
		}
		TR << std::endl;
	}
}

void dump_gast_types( core::chemical::ResidueType const & restype ) {
	for ( core::Size ii(1); ii <= restype.natoms(); ++ii ) {
		TR << ii << " Atom " << restype.atom_name(ii) << " " ;
		if ( restype.atom(ii).gasteiger_atom_type() ) {
			TR << restype.atom(ii).gasteiger_atom_type()->get_name() << std::endl;
		} else {
			TR << " unassigned" << std::endl;
		}
	}
}

// --------------- Test Class --------------- //

class GasteigerAtomTyperTests : public CxxTest::TestSuite {

public:

	GasteigerAtomTypeSetCOP atom_type_set_;

	// --------------- Suite-level Fixture --------------- //

	GasteigerAtomTyperTests() {
	}

	virtual ~GasteigerAtomTyperTests() {}


	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();

		atom_type_set_ = ChemicalManager::get_instance()->gasteiger_atom_type_set();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_safety() {
		// The Gasteiger typing should not crash when presented with any of the standard residue types.
		core::chemical::ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set("fa_standard") );
		//utility::vector1< ResidueTypeCOP > const & residues( residue_set->residue_types_DO_NOT_USE() );
		utility::vector1< ResidueTypeCOP > const & residues( residue_set->base_residue_types() );
		for ( core::Size ii(1); ii <= residues.size(); ++ii ) {
			core::chemical::ResidueTypeOP restype( new core::chemical::ResidueType( *(residues[ii]) ) );
			if ( restype->name() != core::chemical::residue_type_base_name( *restype ) ) {
				continue; //Ignore patched residue types
			}
			// AMW: ignore mineral surface types. Can't handle these (for now?)
			if ( restype->name() == "PHO" || restype->name() == "HYD" ||restype->name() == "CAL" ||restype->name() == "CO3" ) {
				continue;
			}
			TR << "Gasteiger Typing " << restype->name() << std::endl;
			core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );
			//No tests, just make sure it doesn't crash.
		}
	}

	void test_cannonicals() {
		core::chemical::ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set("fa_standard") );
		core::chemical::ResidueTypeOP restype;

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "ASP" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("N").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" ); //Special Cased
		TS_ASSERT_EQUALS( restype->atom("CA").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("C").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("O").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		TS_ASSERT_EQUALS( restype->atom("CB").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("CG").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("OD1").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		TS_ASSERT_EQUALS( restype->atom("OD2").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		TS_ASSERT_EQUALS( restype->atom("H").gasteiger_atom_type()->get_name(), "H_S" );
		TS_ASSERT_EQUALS( restype->atom("HA").gasteiger_atom_type()->get_name(), "H_S" );
		TS_ASSERT_EQUALS( restype->atom("1HB").gasteiger_atom_type()->get_name(), "H_S" );
		TS_ASSERT_EQUALS( restype->atom("2HB").gasteiger_atom_type()->get_name(), "H_S" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "GLY:NtermProteinFull" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("N").gasteiger_atom_type()->get_name(), "N_TeTeTeTe" ); //Not special cased, protonated amine
		TS_ASSERT_EQUALS( restype->atom("CA").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("C").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("O").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "GLY:CtermProteinFull" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("N").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" ); // Still special cased
		TS_ASSERT_EQUALS( restype->atom("CA").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("C").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("O").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" ); //Still trigonal

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "GLN" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CD").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("OE1").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		TS_ASSERT_EQUALS( restype->atom("NE2").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "TYR" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CG").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CD1").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CD2").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CE1").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CE2").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CZ").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		// We don't conjugate the phenolic into the ring
		TS_ASSERT_EQUALS( restype->atom("OH").gasteiger_atom_type()->get_name(), "O_Te2Te2TeTe" );
		TS_ASSERT_EQUALS( restype->atom("HH").gasteiger_atom_type()->get_name(), "H_S" );
		TS_ASSERT_EQUALS( restype->atom("HD2").gasteiger_atom_type()->get_name(), "H_S" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "HIS" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CG").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("ND1").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CD2").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("CE1").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("NE2").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("HE2").gasteiger_atom_type()->get_name(), "H_S" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "HIS_D" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("ND1").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("NE2").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("HD1").gasteiger_atom_type()->get_name(), "H_S" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "TRP" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CD1").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("NE1").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("CD2").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "ARG" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("NE").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("CZ").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("NH1").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("NH2").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "LYS" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("NZ").gasteiger_atom_type()->get_name(), "N_TeTeTeTe" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "MET" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CG").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("SD").gasteiger_atom_type()->get_name(), "S_Te2Te2TeTe" );
		TS_ASSERT_EQUALS( restype->atom("CE").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "CYS" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("CB").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( restype->atom("SG").gasteiger_atom_type()->get_name(), "S_Te2Te2TeTe" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "GUA" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("P").gasteiger_atom_type()->get_name(), "P_TeTeTeTePi" );
		TS_ASSERT_EQUALS( restype->atom("OP1").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		TS_ASSERT_EQUALS( restype->atom("OP2").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );
		//TODO: Current typing makes this a O_Tr2TrTrPi2 -- is there delocalization into the phosphorus?
		//TS_ASSERT_EQUALS( restype->atom("O5'").gasteiger_atom_type()->get_name(), "O_Te2Te2TeTe" );
		TS_ASSERT_EQUALS( restype->atom("C2").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("C4").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("C5").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("C6").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("C8").gasteiger_atom_type()->get_name(), "C_TrTrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("N1").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		// The N2 exocyclic nitrogen should be conjugated into the ring.
		TS_ASSERT_EQUALS( restype->atom("N2").gasteiger_atom_type()->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( restype->atom("N3").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" );
		TS_ASSERT_EQUALS( restype->atom("N7").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" );
	}

	void test_virtuals() {
		core::chemical::ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set("fa_standard") );
		core::chemical::ResidueTypeOP restype;

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "NA" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("NA").gasteiger_atom_type()->get_name(), "Na_" );
		TS_ASSERT_EQUALS( restype->atom("V1").gasteiger_atom_type()->get_name(), "FAKE" );
		TS_ASSERT_EQUALS( restype->atom("V2").gasteiger_atom_type()->get_name(), "FAKE" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "FE" ) ) ); // Formal charge of 3, transition metal
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("FE").gasteiger_atom_type()->get_name(), "Fe_" );

		restype = core::chemical::ResidueTypeOP( new core::chemical::ResidueType( residue_set->name_map( "PRO" ) ) );
		core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( restype->atom("NV").gasteiger_atom_type()->get_name(), "FAKE" );
		TS_ASSERT_EQUALS( restype->atom("CD").gasteiger_atom_type()->get_name(), "C_TeTeTeTe" );
	}

	void test_charge_states() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		///////////// Thiolate

		core::chemical::ResidueTypeOP thiolate_test = read_topology_file("core/chemical/gasteiger/THI.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *thiolate_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( thiolate_test->atom("S1").gasteiger_atom_type()->get_name(), "S_Te2Te2Te2Te" ); // Negatively charged singly bonded sulfur

	}

	// Test to make sure we do decent typing when presented with a residue with missing hydrogens.
	// The items are know corner cases.
	void test_implicit_hydrogens() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		///////////// Sulfurs

		core::chemical::ResidueTypeOP sulfur_test = read_topology_file("core/chemical/gasteiger/GST.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *sulfur_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( sulfur_test->atom("S1").gasteiger_atom_type()->get_name(), "S_TeTeTeTePiPi" ); // Sulfone, 1 implicit hydrogen
		TS_ASSERT_EQUALS( sulfur_test->atom("S2").gasteiger_atom_type()->get_name(), "S_Tr2Tr2TrPi" ); // Thioketone
		TS_ASSERT_EQUALS( sulfur_test->atom("S3").gasteiger_atom_type()->get_name(), "S_Te2Te2TeTe" ); // Thiol, 1 implicit hydrogen
		//TS_ASSERT_EQUALS( sulfur_test->atom("S4").gasteiger_atom_type()->get_name(), "S_Tr2TrTrPi" ); // Sulfoxide, 1 implicit hydrogen (double bonded)
		//RM: Rosetta agrees with BCL, but I think this should probably be S_Te2TeTeTePi for double bonded, or S_Te2TeTeTe if S-O charge separated form.

		//TS_ASSERT_EQUALS( sulfur_test->atom("S5").gasteiger_atom_type()->get_name(), "S_Tr2Tr2TrPi" );
		// S5 was N=S=O, for a 2 implicit hydrogen sulfone, but both Rosetta and BCL choked on it - I had to remove the oxygen

		//TS_ASSERT_EQUALS( sulfur_test->atom("S6").gasteiger_atom_type()->get_name(), "S_Tr2Tr2TrPi" ); // N=S
		// S6 was intended to be a 2 implicit hydrogen sulfoxide (see above), but in retrospect it can't be distinquished from a thioketone-type compound.

		///////////// Sulfoxides

		core::chemical::ResidueTypeOP sulfoxide_test = read_topology_file("core/chemical/gasteiger/SXX.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *sulfoxide_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( sulfoxide_test->atom("S1").gasteiger_atom_type()->get_name(), "S_Te2TeTeTePi" ); // Sulfoxide, double bonded
		TS_ASSERT_EQUALS( sulfoxide_test->atom("S2").gasteiger_atom_type()->get_name(), "S_Te2TeTeTe" );   // Sulfoxide, charge separated
		TS_ASSERT_EQUALS( sulfoxide_test->atom("O1").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );  // Sulfoxide, double bonded
		TS_ASSERT_EQUALS( sulfoxide_test->atom("O2").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); // Sulfoxide, charge separated

		///////////// Phosphates/Phosphines

		core::chemical::ResidueTypeOP phos_test = read_topology_file("core/chemical/gasteiger/PXX.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *phos_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( phos_test->atom("P3").gasteiger_atom_type()->get_name(), "P_Te2TeTeTe" ); // Phosphine, 2 implicit hydrogen
		TS_ASSERT_EQUALS( phos_test->atom("P4").gasteiger_atom_type()->get_name(), "P_Te2TeTeTe" ); // Phosphine, 1 implicit hydrogen
		TS_ASSERT_EQUALS( phos_test->atom("P5").gasteiger_atom_type()->get_name(), "P_Te2TeTeTe" ); // Phosphine, 0 implicit hydrogen

		//TS_ASSERT_EQUALS( phos_test->atom("P6").gasteiger_atom_type()->get_name(), "P_Tr2TrTrPi" );  // Phosphate, 2 implicit hydrogen
		// P6 can't be distinguished from a phosphine-like compound lacking hydrogens

		//TS_ASSERT_EQUALS( phos_test->atom("P1").gasteiger_atom_type()->get_name(), "P_TrTrTrPi" );   // Phosphate, 1 implicit hydrogen
		//RM: Rosetta agrees with BCL, but I think this should probably be P_TeTeTeTePi with an implicit hydrogen instead of a "carboxylate" like hybridization

		TS_ASSERT_EQUALS( phos_test->atom("P2").gasteiger_atom_type()->get_name(), "P_TeTeTeTePi" ); // Phosphate, 0 implicit hydrogen

		//////////// Nitrogens

		core::chemical::ResidueTypeOP nitr_test = read_topology_file("core/chemical/gasteiger/NXX.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *nitr_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( nitr_test->atom("N1").gasteiger_atom_type()->get_name(), "N_Te2TeTeTe" ); // (uncharged) Amine, 2 implicit hydrogen
		TS_ASSERT_EQUALS( nitr_test->atom("N2").gasteiger_atom_type()->get_name(), "N_Te2TeTeTe" ); // (uncharged) Amine, 1 implicit hydrogen
		TS_ASSERT_EQUALS( nitr_test->atom("N3").gasteiger_atom_type()->get_name(), "N_Te2TeTeTe" ); // (uncharged) Amine, 0 implicit hydrogen

		TS_ASSERT_EQUALS( nitr_test->atom("N4").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" ); // Nitroso, no impilict hydrogens

		TS_ASSERT_EQUALS( nitr_test->atom("N5").gasteiger_atom_type()->get_name(), "N_TrTrTrPiPi" ); // Nitro, two double bonded oxygens
		TS_ASSERT_EQUALS( nitr_test->atom("N6").gasteiger_atom_type()->get_name(), "N_TrTrTrPi" );   // Nitro, charge separated

		//////////// Azides

		core::chemical::ResidueTypeOP azo_test = read_topology_file("core/chemical/gasteiger/AZO.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *azo_test, atom_type_set_, /*keep_existing=*/ false );

		TS_ASSERT_EQUALS( azo_test->atom("N1").gasteiger_atom_type()->get_name(), "N_Tr2TrTrPi" ); // -*N*=[N+]=[N-] (two doubles)
		TS_ASSERT_EQUALS( azo_test->atom("N4").gasteiger_atom_type()->get_name(), "N_DiDiPiPi" ); // -N=*[N+]*=[N-]
		TS_ASSERT_EQUALS( azo_test->atom("N3").gasteiger_atom_type()->get_name(), "N_Di2DiPi2Pi" ); // -N=[N+]=*[N-]*
		// RM: I would have though N_Tr2Tr2TrPi, but the sp(Di) hybridization probably works better with resonance forms

		TS_ASSERT_EQUALS( azo_test->atom("N2").gasteiger_atom_type()->get_name(), "N_Te2Te2TeTe" ); // -*[N-]*-[N+]#N  (single-triple)
		TS_ASSERT_EQUALS( azo_test->atom("N5").gasteiger_atom_type()->get_name(), "N_DiDiPiPi" ); // -[N-]-*[N+]*#N
		TS_ASSERT_EQUALS( azo_test->atom("N6").gasteiger_atom_type()->get_name(), "N_Di2DiPiPi" ); // -[N-]-[N+]#*N*

	}


	// Test various oxides, particularly in the case of formal charge separation
	void test_oxides() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		///////////// Sulfates and Sulfites

		core::chemical::ResidueTypeOP sulfur_test = read_topology_file("core/chemical/gasteiger/SOx.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *sulfur_test, atom_type_set_, /*keep_existing=*/ false, /*allow_unknown=*/ true );

		dump_gast_types( *sulfur_test );

		// Regular sulfate
		TS_ASSERT_EQUALS( sulfur_test->atom("S1").gasteiger_atom_type()->get_name(), "S_TeTeTeTePiPi" );
		TS_ASSERT_EQUALS( sulfur_test->atom("O1").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O5").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );  //S=O
		TS_ASSERT_EQUALS( sulfur_test->atom("O6").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );  //S=O
		// Sulfate, 1 oxygen charge separated
		TS_ASSERT_EQUALS( sulfur_test->atom("S2").gasteiger_atom_type()->get_name(), "S_TeTeTeTePi" );
		TS_ASSERT_EQUALS( sulfur_test->atom("O7").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );  //S=O
		TS_ASSERT_EQUALS( sulfur_test->atom("O8").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O9").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		// Sulfate, 2 oxygen charge separated
		TS_ASSERT_EQUALS( sulfur_test->atom("S3").gasteiger_atom_type()->get_name(), "S_TeTeTeTe" );
		TS_ASSERT_EQUALS( sulfur_test->atom("O2").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" );  //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O10").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O11").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		// Regular sulfite
		TS_ASSERT_EQUALS( sulfur_test->atom("S4").gasteiger_atom_type()->get_name(), "S_Te2TeTeTePi" );
		TS_ASSERT_EQUALS( sulfur_test->atom("O3").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O12").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" ); //S=O
		// Sulfite, charge separated
		TS_ASSERT_EQUALS( sulfur_test->atom("S5").gasteiger_atom_type()->get_name(), "S_Te2TeTeTe" );
		TS_ASSERT_EQUALS( sulfur_test->atom("O4").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" );  //S-O(-)
		TS_ASSERT_EQUALS( sulfur_test->atom("O13").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //S-O(-)

		core::chemical::ResidueTypeOP phosphorus_test = read_topology_file("core/chemical/gasteiger/POx.params",
			atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

		core::chemical::gasteiger::assign_gasteiger_atom_types( *phosphorus_test, atom_type_set_, /*keep_existing=*/ false, /*allow_unknown=*/ true );

		/// Regular phosphate
		TS_ASSERT_EQUALS( phosphorus_test->atom("P1").gasteiger_atom_type()->get_name(), "P_TeTeTeTePi" );
		TS_ASSERT_EQUALS( phosphorus_test->atom("O1").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		TS_ASSERT_EQUALS( phosphorus_test->atom("O3").gasteiger_atom_type()->get_name(), "O_Tr2Tr2TrPi" );  //P=O
		TS_ASSERT_EQUALS( phosphorus_test->atom("O4").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		// Phosphate, charge separated
		TS_ASSERT_EQUALS( phosphorus_test->atom("P2").gasteiger_atom_type()->get_name(), "P_TeTeTeTe" );
		TS_ASSERT_EQUALS( phosphorus_test->atom("O5").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		TS_ASSERT_EQUALS( phosphorus_test->atom("O6").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		TS_ASSERT_EQUALS( phosphorus_test->atom("O7").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		// Phosphonite
		TS_ASSERT_EQUALS( phosphorus_test->atom("P3").gasteiger_atom_type()->get_name(), "P_Te2TeTeTe" );
		TS_ASSERT_EQUALS( phosphorus_test->atom("O8").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)
		TS_ASSERT_EQUALS( phosphorus_test->atom("O9").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" );  //P-O(-)
		// Phosphonite, protonated
		TS_ASSERT_EQUALS( phosphorus_test->atom("P4").gasteiger_atom_type()->get_name(), "P_TeTeTeTe" );
		TS_ASSERT_EQUALS( phosphorus_test->atom("O2").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" );  //P-O(-)
		TS_ASSERT_EQUALS( phosphorus_test->atom("O10").gasteiger_atom_type()->get_name(), "O_Te2Te2Te2Te" ); //P-O(-)

	}

};


