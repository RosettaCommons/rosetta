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
	for( boost::tie( aitr, aitr_end ) = boost::vertices( graph ); aitr != aitr_end; ++aitr ) {
		TR << " Atom " << graph[ *aitr ].name() << " -- ";
		OEiter eitr, eitr_end;
		for( boost::tie( eitr, eitr_end ) = boost::out_edges( *aitr, graph ); eitr !=  eitr_end; ++eitr_end ) {
			TR << " " << graph[ boost::target( *eitr, graph ) ].name();
		}
		TR << std::endl;
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
		utility::vector1< ResidueTypeCOP > const & residues( residue_set->residue_types_DO_NOT_USE() );
		for( core::Size ii(1); ii <= residues.size(); ++ii ) {
			core::chemical::ResidueTypeOP restype( new core::chemical::ResidueType( *(residues[ii]) ) );
			if( restype->name() != core::chemical::residue_type_base_name( *restype ) ) {
				continue; //Ignore patched residue types
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

		// This is being typed as a free amine, rather than as an amide (because we don't have connection info)
		TS_ASSERT_EQUALS( restype->atom("N").gasteiger_atom_type()->get_name(), "N_Te2TeTeTe" );
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

};


