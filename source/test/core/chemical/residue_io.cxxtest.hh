// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/residue_io.cxxtest.hh
/// @brief unit tests for the residue_io file
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/residue_io.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>

#include <core/types.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

#include <test/UTracer.hh>

static basic::Tracer TR("core.chemical.restype_io.cxxtest");

class residue_io_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_graph_out() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		ResidueTypeCOP rsd_ref( rsd_types->name_mapOP("TYR") );

		test::UTracer UT("core/chemical/TYR.dot.u");

		write_graphviz( *rsd_ref, UT );
	}

	void test_param_from_istream() {


		//create a 'params file' for some silly ligand as a string, hardcoded

		// NAME LG1
		// IO_STRING LG1 Z
		// TYPE LIGAND
		// AA UNK
		// ATOM  C3  CH1   X   0.14
		// ATOM  C2  CH1   X   0.14
		// ATOM  C1  CH1   X   0.14
		// ATOM  O1  OOC   X   -0.53
		// ATOM  P1  Phos  X   1.73
		// ATOM  O4  OOC   X   -0.53
		// ATOM  O2  OOC   X   -0.53
		// ATOM  O3  OOC   X   -0.53
		// BOND_TYPE  C1   C2  1
		// BOND_TYPE  C3   C2  1
		// BOND_TYPE  C2   O1  1
		// BOND_TYPE  C3   P1  1
		// BOND_TYPE  O4   P1  2
		// BOND_TYPE  P1   O2  1
		// BOND_TYPE  P1   O3  1
		// CHI 1  P1   C3   C2   C1
		// CHI 2  C2   C3   P1   O4
		// NBR_ATOM  C3
		// NBR_RADIUS 3.328921
		// ICOOR_INTERNAL    C3     0.000000    0.000000    0.000000   C3    C2    C1
		// ICOOR_INTERNAL    C2     0.000000  180.000000    1.511813   C3    C2    C1
		// ICOOR_INTERNAL    C1     0.000000   69.809952    1.510875   C2    C3    C1
		// ICOOR_INTERNAL    O1   120.854877   70.146622    1.407682   C2    C3    C1
		// ICOOR_INTERNAL    P1   175.284710   68.793359    1.817238   C3    C2    C1
		// ICOOR_INTERNAL    O4  -142.931649   70.965715    1.510131   P1    C3    C2
		// ICOOR_INTERNAL    O2   119.487828   70.765017    1.510934   P1    C3    O4
		// ICOOR_INTERNAL    O3   120.276972   70.005096    1.511683   P1    C3    O2

		std::string const LG1("NAME LG1\nIO_STRING LG1 Z\nTYPE LIGAND\nAA UNK\nATOM  C3  CH1   X   0.14\nATOM  C2  CH1   X   0.14\nATOM  C1  CH1   X   0.14\nATOM  O1  OOC   X   -0.53\nATOM  P1  Phos  X   1.73\nATOM  O4  OOC   X   -0.53\nATOM  O2  OOC   X   -0.53\nATOM  O3  OOC   X   -0.53\nBOND_TYPE  C1   C2  1   \nBOND_TYPE  C3   C2  1   \nBOND_TYPE  C2   O1  1   \nBOND_TYPE  C3   P1  1   \nBOND_TYPE  O4   P1  2   \nBOND_TYPE  P1   O2  1   \nBOND_TYPE  P1   O3  1   \nCHI 1  P1   C3   C2   C1 \nCHI 2  C2   C3   P1   O4 \nNBR_ATOM  C3 \nNBR_RADIUS 3.328921\nICOOR_INTERNAL    C3     0.000000    0.000000    0.000000   C3    C2    C1 \nICOOR_INTERNAL    C2     0.000000  180.000000    1.511813   C3    C2    C1 \nICOOR_INTERNAL    C1     0.000000   69.809952    1.510875   C2    C3    C1 \nICOOR_INTERNAL    O1   120.854877   70.146622    1.407682   C2    C3    C1 \nICOOR_INTERNAL    P1   175.284710   68.793359    1.817238   C3    C2    C1 \nICOOR_INTERNAL    O4  -142.931649   70.965715    1.510131   P1    C3    C2 \nICOOR_INTERNAL    O2   119.487828   70.765017    1.510934   P1    C3    O4 \nICOOR_INTERNAL    O3   120.276972   70.005096    1.511683   P1    C3    O2 \n");

		//fill an istream with this params file
		std::istringstream ligand_stream;
		ligand_stream.str(LG1);

		//try reading this params file
		//TS_ASSERT a bunch of stuff about the resulting ResidueTypeOP

		core::chemical::ResidueTypeSetCOP rtsCOP(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ));

		core::chemical::ResidueTypeOP LG1_RTOP(core::chemical::read_topology_file( ligand_stream, "dummy_filename", rtsCOP ));

		//This is a collection of tests based narrowly on "functions of ResidueType that are easy to use and verify by glancing at the params string", improvements welcome
		TS_ASSERT(LG1_RTOP);
		TS_ASSERT_EQUALS(LG1_RTOP->nbr_atom(), LG1_RTOP->atom_index("C3"));
		TS_ASSERT_EQUALS(LG1_RTOP->nbr_radius(), 3.328921);
		TS_ASSERT_EQUALS(LG1_RTOP->natoms(), 8);
		TS_ASSERT_EQUALS(LG1_RTOP->nbonds(), 7);
		TS_ASSERT(  LG1_RTOP->atoms_are_bonded(LG1_RTOP->atom_index("C3"), LG1_RTOP->atom_index("P1")));
		TS_ASSERT(! LG1_RTOP->atoms_are_bonded(LG1_RTOP->atom_index("C3"), LG1_RTOP->atom_index("O3")));
		TS_ASSERT(! LG1_RTOP->atom_is_polar_hydrogen(LG1_RTOP->atom_index("C3")));

		TS_ASSERT(LG1_RTOP->has("C3"));
		TS_ASSERT(LG1_RTOP->has("C2"));
		TS_ASSERT(LG1_RTOP->has("C1"));
		TS_ASSERT(LG1_RTOP->has("O1"));
		TS_ASSERT(LG1_RTOP->has("P1"));
		TS_ASSERT(LG1_RTOP->has("O4"));
		TS_ASSERT(LG1_RTOP->has("O2"));
		TS_ASSERT(LG1_RTOP->has("O3"));

		TS_ASSERT_EQUALS(LG1_RTOP->n_polymeric_residue_connections(), 0);
		TS_ASSERT(LG1_RTOP->finalized());

		TS_ASSERT(! LG1_RTOP->is_polymer());
		TS_ASSERT(! LG1_RTOP->is_sidechain_thiol());
		TS_ASSERT(! LG1_RTOP->is_disulfide_bonded());
		TS_ASSERT(! LG1_RTOP->is_sidechain_amine());
		TS_ASSERT(! LG1_RTOP->is_protein());
		TS_ASSERT(! LG1_RTOP->is_alpha_aa());
		TS_ASSERT(! LG1_RTOP->is_beta_aa());
		TS_ASSERT(! LG1_RTOP->is_gamma_aa());
		TS_ASSERT(! LG1_RTOP->is_sri());
		TS_ASSERT(! LG1_RTOP->is_triazolemer());
		TS_ASSERT(! LG1_RTOP->is_d_aa());
		TS_ASSERT(! LG1_RTOP->is_l_aa());
		TS_ASSERT(! LG1_RTOP->is_achiral_backbone());
		TS_ASSERT(! LG1_RTOP->is_DNA());
		TS_ASSERT(! LG1_RTOP->is_RNA());
		TS_ASSERT(! LG1_RTOP->is_coarse());
		TS_ASSERT(! LG1_RTOP->is_NA());
		TS_ASSERT(! LG1_RTOP->is_peptoid());
		TS_ASSERT(! LG1_RTOP->is_carbohydrate());
		TS_ASSERT(  LG1_RTOP->is_ligand());
		TS_ASSERT(! LG1_RTOP->is_lipid());
		TS_ASSERT(! LG1_RTOP->is_metal());
		TS_ASSERT(! LG1_RTOP->is_metalbinding());
		TS_ASSERT(! LG1_RTOP->is_membrane());
		TS_ASSERT(! LG1_RTOP->is_surface());
		TS_ASSERT(! LG1_RTOP->has_sc_orbitals());
		TS_ASSERT(! LG1_RTOP->is_polar());
		TS_ASSERT(! LG1_RTOP->is_charged());
		TS_ASSERT(! LG1_RTOP->is_aromatic());
		TS_ASSERT(! LG1_RTOP->is_cyclic());
		TS_ASSERT(! LG1_RTOP->is_terminus());
		TS_ASSERT(! LG1_RTOP->is_lower_terminus());
		TS_ASSERT(! LG1_RTOP->is_upper_terminus());
		TS_ASSERT(! LG1_RTOP->is_branch_point());
		TS_ASSERT(! LG1_RTOP->is_acetylated_nterminus());
		TS_ASSERT(! LG1_RTOP->is_methylated_cterminus());
		TS_ASSERT(! LG1_RTOP->is_virtual_residue());
		TS_ASSERT(! LG1_RTOP->is_adduct());

		TS_ASSERT_EQUALS( LG1_RTOP->name(), "LG1");
		TS_ASSERT_EQUALS( LG1_RTOP->name3(), "LG1");
		TS_ASSERT_EQUALS( LG1_RTOP->name1(), 'Z');
		TS_ASSERT_EQUALS( LG1_RTOP->aa(), core::chemical::aa_unk);

	}

};
