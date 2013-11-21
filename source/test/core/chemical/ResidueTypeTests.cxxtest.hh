// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeTests.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Matthew O'Meara


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>

// AUTO-REMOVED #include <core/chemical/ElementSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalTypeSet.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>

//Auto Headers


using std::endl;
using std::string;
using basic::Tracer;
using core::chemical::AtomTypeSetCAP;
using core::chemical::AtomType;
using core::chemical::ChemicalManager;
using core::chemical::ElementSetCAP;
using core::chemical::MMAtomTypeSetCAP;
using core::chemical::orbitals::OrbitalTypeSetCAP;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeSet;
using core::chemical::FA_STANDARD;
using utility::vector1;

static Tracer TR("core.chemical.ResidueTypeTests.cxxtest");

void add_atom(
		ResidueType & rsd,
		AtomTypeSetCAP & atom_types,
		string const& name,
		string const& type,
		string const& mm_type,
		core::Real const& charge
){
	core::Size natoms = rsd.natoms();
	core::Size nheavyatoms = rsd.nheavyatoms();
	core::Size nhbond_acceptors = rsd.n_hbond_acceptors();
	core::Size nhbond_donors = rsd.n_hbond_donors();
	core::Size Haro_size = rsd.Haro_index().size();
	core::Size Hpol_size = rsd.Hpol_index().size();
	core::Size Hpos_polar = rsd.Hpos_polar().size();
	core::Size Hpos_apolar = rsd.Hpos_apolar().size();
	core::Size accpt_pos = rsd.accpt_pos().size();

	TS_ASSERT(atom_types->has_atom(type));
	AtomType const& at = (*atom_types)[ atom_types->atom_type_index(type)];
	++natoms;
	if(at.is_heavyatom()) ++nheavyatoms;
	if(at.is_acceptor()){
		++nhbond_acceptors;
		++accpt_pos;
	}
	if(at.is_donor()) ++nhbond_donors;
	if(at.is_haro()) ++Haro_size;
	if(at.is_polar_hydrogen()){
		++Hpol_size;
		if(std::abs(charge > 1.0e-3)) ++Hpos_polar;
	}else if(std::abs(charge > 1.0e-3)) ++Hpos_apolar;

	rsd.add_atom( name, type, mm_type, charge);

	TS_ASSERT_EQUALS(rsd.natoms(), natoms);
	//TS_ASSERT_EQUALS(rsd.nheavyatoms(), nheavyatoms); // Only updated with finalize
	TS_ASSERT_EQUALS(rsd.n_hbond_acceptors(), nhbond_acceptors);
	TS_ASSERT_EQUALS(rsd.n_hbond_donors(), nhbond_donors);
	//TS_ASSERT_EQUALS(rsd.Haro_index().size(), Haro_size); // Only updated with finalize
	//TS_ASSERT_EQUALS(rsd.Hpol_index().size(), Hpol_size); // Only updated with finalize
	//TS_ASSERT_EQUALS(rsd.Hpos_polar().size(), Hpos_polar); // Only updated with finalize
	//TS_ASSERT_EQUALS(rsd.Hpos_apolar().size(), Hpos_apolar); // Only updated with finalize
	//TS_ASSERT_EQUALS(rsd.accpt_pos().size(), accpt_pos); // Only updated with finalize
}

void add_bond(ResidueType & rsd, std::string const& a1, std::string const& a2){
	rsd.add_bond(a1,a2);
	//TS_ASSERT_EQUALS( rsd.path_distance(rsd.atom_index(a1), rsd.atom_index(a2)), 1); // Only updated with finalize
}

void do_retyping(ResidueType const & ref, bool build_emap = true) {
	core::chemical::ElementMap emap;
	if( build_emap ) {
		for( core::Size ii(1); ii <= ref.natoms(); ++ii ) {
			emap[ ref.vd_from_index( ii ) ] = ref.atom_type(ii).element();
		}
	}
	ResidueType rsd( ref );
	rsd.retype_atoms(emap);
	for( core::Size ii(1); ii <= ref.natoms(); ++ii ) {
		TS_ASSERT_EQUALS( rsd.atom_type(ii).name(), ref.atom_type(ii).name());
	}
}

class ResidueTypeTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_residue_type_initialization() {
		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCAP atom_types = cm->atom_type_set(tag);
		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCAP element_types = cm->element_set(tag);
		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCAP mm_atom_types = cm->mm_atom_type_set(tag);
		// minirosetta_database/chemical/orbital_type_sets/<tag>
		OrbitalTypeSetCAP orbital_types = cm->orbital_type_set(tag);

		ResidueType rsd( atom_types, element_types, mm_atom_types, orbital_types);

		TS_ASSERT_EQUALS( &(*atom_types), &rsd.atom_type_set());
		TS_ASSERT_EQUALS(atom_types, rsd.atom_type_set_ptr());
		TS_ASSERT_EQUALS(&(*element_types), &rsd.element_set());

		TS_ASSERT_EQUALS( rsd.natoms(), (core::Size) 0);
		TS_ASSERT_EQUALS( rsd.nheavyatoms(), (core::Size) 0);

		rsd.add_property("PROTEIN");
		TS_ASSERT(rsd.has_property("PROTEIN"));

		rsd.add_numeric_property("foo",1.5);
		TS_ASSERT_EQUALS(rsd.get_numeric_property("foo"),1.5);

/// Build an Alanine
		add_atom( rsd, atom_types," N  ", "Nbb", "NH1", -0.47);
		add_atom( rsd, atom_types," CA ", "CAbb", "CT1", 0.07);
		add_atom( rsd, atom_types," C  ", "CObb", "C", 0.51);
		add_atom( rsd, atom_types," O  ", "OCbb", "O", -0.51);
		add_atom( rsd, atom_types," CB ", "CH3", "CT3", -0.27);
		add_atom( rsd, atom_types," H  ", "HNbb", "H", 0.31);
		add_atom( rsd, atom_types," HA ", "Hapo", "HB", 0.09);
		add_atom( rsd, atom_types,"1HB ", "Hapo", "HA", 0.09);
		add_atom( rsd, atom_types,"2HB ", "Hapo", "HA", 0.09);
		add_atom( rsd, atom_types,"3HB ", "Hapo", "HA", 0.09);

		rsd.set_lower_connect_atom("N");
		rsd.set_upper_connect_atom("C");

		TS_ASSERT_EQUALS( rsd.lower_connect_atom(), rsd.atom_index("N"));
		TS_ASSERT_EQUALS( rsd.upper_connect_atom(), rsd.atom_index("C"));

		add_bond( rsd, "N", "CA");
		add_bond( rsd, "N", "H");
		add_bond( rsd, "CA", "C");
		add_bond( rsd, "CA", "CB");
		add_bond( rsd, "CA", "HA");
		add_bond( rsd, "C", "O");
		add_bond( rsd, "CB", "1HB");
		add_bond( rsd, "CB", "2HB");
		add_bond( rsd, "CB", "3HB");


		rsd.nbr_atom("CB");
		TS_ASSERT_EQUALS(rsd.nbr_atom(), rsd.atom_index("CB"));
		rsd.nbr_radius(3.4473);
		TS_ASSERT_EQUALS(rsd.nbr_radius(), 3.4473);
		rsd.add_actcoord_atom("CB");
		TS_ASSERT_EQUALS(rsd.actcoord_atoms().size(), (core::Size) 1);

		rsd.set_icoor("N", 0.000000, 0.000000, 0.000000, "N",  "CA", "C");
		rsd.set_icoor("CA", 0.000000, 180.000000, 1.458001, "N",  "CA", "C");
		rsd.set_icoor("C", 0.000000, 68.800003, 1.523258, "CA",  "N", "C");
		rsd.set_icoor("UPPER",  149.999985,  63.800007, 1.328685, "C",  "CA", "N");
		rsd.set_icoor("O", -180.000000, 59.200005, 1.231015, "C",  "CA", "UPPER");
		rsd.set_icoor("CB", -123.028732, 69.625412, 1.521736, "CA",  "N", "C");
		rsd.set_icoor("1HB",-179.994110, 70.562309, 1.090040, "CB",  "CA", "N");
		rsd.set_icoor("2HB", 119.960403, 70.475021, 1.090069, "CB",  "CA", "1HB");
		rsd.set_icoor("3HB", 120.089012, 70.493805, 1.088803, "CB",  "CA", "2HB");
		rsd.set_icoor("HA", -119.013138, 71.295197, 1.090078, "CA",  "N", "CB");
		rsd.set_icoor("LOWER", -150.000000,  58.300003, 1.328685, "N",  "CA", "C");
		rsd.set_icoor("H", -180.000000,  60.849998,   1.010000, "N",  "CA", "LOWER");

		rsd.finalize();

		//TS_ASSERT_EQUALS( rsd.chi_atoms( rsd.atom_index("CA")).size(), (core::Size) 3); // we have to "add_chi" first and ALA has no chi listed

		TS_ASSERT_EQUALS( rsd.path_distance(rsd.atom_index("N"), rsd.atom_index("C")), 2);
		TS_ASSERT_EQUALS( rsd.path_distance(rsd.atom_index("N"), rsd.atom_index("O")), 3);

		TS_ASSERT_EQUALS( rsd.attached_H_begin( rsd.atom_index("N")), (core::Size) 6);
		TS_ASSERT_EQUALS( rsd.attached_H_end( rsd.atom_index("N")), (core::Size) 6);
		TS_ASSERT_EQUALS( rsd.attached_H_begin( rsd.atom_index("CB")), (core::Size) 8);
		TS_ASSERT_EQUALS( rsd.attached_H_end( rsd.atom_index("CB")), (core::Size) 10);
		TS_ASSERT_EQUALS( rsd.number_bonded_hydrogens( rsd.atom_index("N")), (core::Size) 1);
		TS_ASSERT_EQUALS( rsd.number_bonded_hydrogens( rsd.atom_index("CB")), (core::Size) 3);
		TS_ASSERT_EQUALS( rsd.attached_H_end().size(), (core::Size) 5);
		TS_ASSERT_EQUALS( rsd.attached_H_begin().size(), (core::Size) 5);
        utility::vector1<core::Size> H_begin(rsd.attached_H_begin());
        TS_ASSERT_EQUALS(rsd.atom_name(H_begin[1]), " H  ");
        TS_ASSERT_EQUALS(rsd.atom_name(H_begin[4]), "1HB ");
        utility::vector1<core::Size> H_end(rsd.attached_H_end());
        TS_ASSERT_EQUALS(rsd.atom_name(H_end[1]), " H  ");
        TS_ASSERT_EQUALS(rsd.atom_name(H_end[2]), " HA ");
		TS_ASSERT_EQUALS( rsd.bonded_neighbor( rsd.atom_index("N")).size(), (core::Size) 2);
		TS_ASSERT_EQUALS( rsd.bonded_neighbor( rsd.atom_index("CA")).size(), (core::Size) 4);
		TS_ASSERT_EQUALS( rsd.bonded_neighbor( rsd.atom_index("O")).size(), (core::Size) 1);
		TS_ASSERT_EQUALS( rsd.bonded_neighbor( rsd.atom_index("CA")), rsd.nbrs( rsd.atom_index("CA")));
		TS_ASSERT( ! rsd.heavyatom_has_polar_hydrogens( rsd.atom_index("CB"))); // should be > 0
		TS_ASSERT( rsd.heavyatom_is_an_acceptor( rsd.atom_index("O")));
		TS_ASSERT( ! rsd.atom_is_polar_hydrogen( rsd.atom_index("O")));
		TS_ASSERT_EQUALS(rsd.mainchain_atoms().size(), (core::Size) 0); // Number of MC atoms should be 4? // Is this even used?
		//TS_ASSERT_EQUALS(rsd.mainchain_atom( rsd.atom_index("N")), (core::Size) 1); // Index of MC atom should be 1
		TS_ASSERT( ! rsd.has( "BURRITO" ));


		TS_ASSERT_EQUALS(rsd.all_bb_atoms().size(), (core::Size) 0); // Why are all atoms being counted as side chain atoms?
		TS_ASSERT_EQUALS(rsd.all_sc_atoms().size(), (core::Size) 10); // Why are all atoms being called side chain atoms?
		TS_ASSERT_EQUALS(rsd.Hpos_polar_sc().size(), (core::Size) 1);
		TS_ASSERT_EQUALS(rsd.accpt_pos_sc().size(), (core::Size) 1);


		TS_ASSERT_EQUALS(rsd.atom_type(1).name(), "Nbb");
		TS_ASSERT_EQUALS(rsd.atom_type(2).name(), "CAbb");
		TS_ASSERT_EQUALS(rsd.atom_type(5).name(), "CH3");

		TS_ASSERT_EQUALS(rsd.mm_atom_type(9).name(), "HA");
		TS_ASSERT_EQUALS(rsd.mm_atom_type(3).name(), "C");

		TS_ASSERT_EQUALS(rsd.atom_base(2), (core::Size) 1);
		TS_ASSERT_EQUALS(rsd.atom_base(3), (core::Size) 2);
		TS_ASSERT_EQUALS(rsd.atom_base(6), (core::Size) 1);
		TS_ASSERT_EQUALS(rsd.abase2(3), (core::Size) 0);
		TS_ASSERT_EQUALS(rsd.abase2(5), (core::Size) 0);
		rsd.set_atom_base("N", "CA");
		rsd.finalize();
		TS_ASSERT_EQUALS(rsd.atom_base(rsd.atom_index("N")), (core::Size) 2);
		TS_ASSERT_EQUALS(rsd.atom_base(rsd.atom_index("CA")), (core::Size) 1);


    	core::Size center=0;
    	core::Size nbr1=0;
    	core::Size nbr=0;
    	rsd.select_orient_atoms(center, nbr1, nbr);

    	TS_ASSERT_EQUALS(center, (core::Size) 5);
    	TS_ASSERT_EQUALS(nbr1, (core::Size) 2);
    	TS_ASSERT_EQUALS(nbr, (core::Size) 10);
        rsd.add_cut_bond(rsd.atom_name(1), rsd.atom_name(2));
        rsd.add_cut_bond(rsd.atom_name(2), rsd.atom_name(3));
        rsd.add_cut_bond(rsd.atom_name(2), rsd.atom_name(4));
        rsd.finalize();

        utility::vector1<core::Size> neigh(rsd.cut_bond_neighbor(2));
        TS_ASSERT_EQUALS(rsd.atom_name(neigh[1]), " N  ");
        TS_ASSERT_EQUALS(rsd.atom_name(neigh[2]), " C  ");
        TS_ASSERT_EQUALS(rsd.atom_name(neigh[3]), " O  ");


		rsd.delete_atom("3HB");
		rsd.delete_atom("2HB");
		rsd.delete_atom("1HB");
		rsd.delete_atom("HA");
		rsd.delete_atom("H");
		rsd.delete_atom("CB");

		TS_ASSERT_EQUALS(rsd.natoms(), 4);

		rsd.finalize();

		rsd.delete_atom("O"); // can't delete these
		rsd.delete_atom("C");
		rsd.delete_atom("CA");
		rsd.delete_atom("N");

		TS_ASSERT_EQUALS(rsd.natoms(), 0);
		//rsd.finalize(); // residue type doesn't support removing the backbone...

	}

	void test_retyping() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		AtomTypeSetCAP atom_types = cm->atom_type_set(tag);
		ElementSetCAP element_types = cm->element_set(tag);
		MMAtomTypeSetCAP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCAP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSet rsd_types;

		utility::io::izstream paramslist("core/chemical/params/retype_list.txt");
		std::string filename;
		paramslist >> filename;
		while( paramslist ) {
			TR << "Retyping " << filename << std::endl;
			core::chemical::ResidueTypeOP rsd = read_topology_file("core/chemical/"+filename,
					atom_types, element_types, mm_atom_types, orbital_types, &rsd_types);
			do_retyping(*rsd);
			paramslist >> filename;
		}

		// Now test autodiscovery of elements.
		filename = "params/1aq1.mol2.params";
		ResidueTypeOP rsd = read_topology_file("core/chemical/"+filename,
				atom_types, element_types, mm_atom_types, orbital_types, &rsd_types);
		do_retyping(*rsd,false);

		// TODO: Test partial reassignment.
	}
};
