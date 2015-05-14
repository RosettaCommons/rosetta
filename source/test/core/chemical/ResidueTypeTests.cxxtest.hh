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
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/util.hh>

#include <core/chemical/sdf/mol_writer.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/Stub.hh>

#include <utility/graph/RingDetection.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/bond_support.hh>
#include <utility/graph/BFS_prune.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.io.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

// C++ Headers
#include <string>
#include <ostream>
#include <boost/graph/adjacency_list.hpp>

using std::endl;
using std::string;
using basic::Tracer;
using core::chemical::AtomTypeSetCOP;
using core::chemical::AtomType;
using core::chemical::ChemicalManager;
using core::chemical::ElementSetCOP;
using core::chemical::MMAtomTypeSetCOP;
using core::chemical::orbitals::OrbitalTypeSetCOP;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeSet;
using core::chemical::FA_STANDARD;
using utility::vector1;

static Tracer TR("core.chemical.ResidueTypeTests.cxxtest");

void add_atom(
		ResidueType & rsd,
		AtomTypeSetCOP & atom_types,
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
		if(std::abs(charge) > 1.0e-3) ++Hpos_polar;
	}else if(std::abs(charge) > 1.0e-3) ++Hpos_apolar;

	rsd.add_atom( name, type, mm_type, charge);

	TS_ASSERT_EQUALS(rsd.natoms(), natoms);
	TS_ASSERT_EQUALS(rsd.n_hbond_acceptors(), nhbond_acceptors);
	TS_ASSERT_EQUALS(rsd.n_hbond_donors(), nhbond_donors);
}

void add_bond(ResidueType & rsd, std::string const& a1, std::string const& a2){
	rsd.add_bond(a1,a2);
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
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

		ResidueType rsd( atom_types, element_types, mm_atom_types, orbital_types);

		TS_ASSERT_EQUALS( &(*atom_types), &rsd.atom_type_set());
		TS_ASSERT_EQUALS(atom_types, rsd.atom_type_set_ptr());
		TS_ASSERT_EQUALS(&(*element_types), &rsd.element_set());

		TS_ASSERT_EQUALS( rsd.natoms(), (core::Size) 0);
		TS_ASSERT_EQUALS( rsd.nheavyatoms(), (core::Size) 0);

		rsd.add_property( "PROTEIN" );
		TS_ASSERT( rsd.has_property( "PROTEIN" ) );

		rsd.add_numeric_property( "foo", 1.5 );
		TS_ASSERT_EQUALS( rsd.get_numeric_property( "foo" ), 1.5 );

		// Build an "Alanine"
		add_atom( rsd, atom_types," N ", "Nbb", "NH1", -0.47);
		add_atom( rsd, atom_types," CA ", "CAbb", "CT1", 0.07);
		add_atom( rsd, atom_types," C ", "CObb", "C", 0.51);
		add_atom( rsd, atom_types," O ", "OCbb", "O", -0.51);
		add_atom( rsd, atom_types," CB ", "CH3", "CT3", -0.27);
		add_atom( rsd, atom_types," H ", "HNbb", "H", 0.31);
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

		rsd.nbr_atom("CB");
		rsd.nbr_radius(3.4473);
		rsd.add_actcoord_atom("CB");


		rsd.finalize();


		TS_ASSERT_EQUALS(rsd.nbr_atom(), rsd.atom_index("CB"));
		TS_ASSERT_EQUALS(rsd.nbr_radius(), 3.4473);
		TS_ASSERT_EQUALS(rsd.actcoord_atoms().size(), (core::Size) 1);

		//TS_ASSERT_EQUALS( rsd.chi_atoms( rsd.atom_index("CA")).size(), (core::Size) 3); // we have to "add_chi" first and ALA has no chi listed

		// For the purposes of a unit test, we'll just add_chi (and add_nu, too) anyway, even though both are meaning-
		// less for alanine. ~Labonte
		rsd.add_chi(1, "N", "CA", "CB", "1HB");
		rsd.add_chi("N", "CA", "CB", "2HB"); // alternate designations, just for kicks
		rsd.add_chi("N", "CA", "CB", "3HB");
		rsd.add_nu(1, "C", "CA", "CB", "1HB"); // not really a nu angle, but it doesn't matter for testing purposes

		rsd.finalize();

		TS_ASSERT_EQUALS(rsd.chi_atoms().size(), 3); // We specified three chi angles above.
		TS_ASSERT_EQUALS(rsd.chi_atoms(3).size(), 4); // There should be four atom indices for any torsion angle.
		TS_ASSERT_EQUALS(rsd.nu_atoms().size(), 1); // We specified one nu angle above.
		TS_ASSERT_EQUALS(rsd.nu_atoms(1).at(1), rsd.atom_index("C")); // 1st atom index in the list should be that of C.

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
		TS_ASSERT_EQUALS(rsd.atom_name(H_begin[1]), " H ");
		TS_ASSERT_EQUALS(rsd.atom_name(H_begin[4]), "1HB ");
		utility::vector1<core::Size> H_end(rsd.attached_H_end());
		TS_ASSERT_EQUALS(rsd.atom_name(H_end[1]), " H ");
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
		TS_ASSERT_EQUALS(rsd.atom_name(neigh[1]), " N ");
		TS_ASSERT_EQUALS(rsd.atom_name(neigh[2]), " C ");
		TS_ASSERT_EQUALS(rsd.atom_name(neigh[3]), " O ");

		rsd.delete_atom("3HB");
		rsd.delete_atom("2HB");
		rsd.delete_atom("1HB");
		rsd.delete_atom("HA");
		rsd.delete_atom("H");
		rsd.delete_atom("CB");

		TS_ASSERT_EQUALS(rsd.natoms(), 4);

		//rsd.finalize(); // can't remove the above atoms if they are part of a CHI definition and then finalize

		rsd.delete_atom("O"); // can't delete these
		rsd.delete_atom("C");
		rsd.delete_atom("CA");
		rsd.delete_atom("N");

		TS_ASSERT_EQUALS(rsd.natoms(), 0);
		//rsd.finalize(); // residue type doesn't support removing the backbone...
	}

	void test_vd_has() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::chemical::ResidueType rsd(pose.residue_type(1));

		for( core::Size ii(1); ii <= rsd.natoms(); ++ii ) {
			core::chemical::VD vd( rsd.atom_vertex( ii ) );
			TR << "Atom: " << rsd.atom_name(ii) << " VD: " << vd << std::endl;
			TS_ASSERT( rsd.has( vd ) );
		}

		//To make sure we are not testing the same residue type.
		TS_ASSERT( pose.residue_type(1).name() != pose.residue_type(2).name() );
		TS_ASSERT( ! rsd.has( pose.residue_type(2).atom_vertex( 1 ) ) );
	}

	void test_atom_add() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

		core::chemical::ResidueType rsd(atom_types, element_types, mm_atom_types, orbital_types);
		rsd.add_atom("ONE");
		rsd.add_atom("MULT", "aroC", "VIRT", -1.5);
		TS_ASSERT( rsd.has("ONE") );
		TS_ASSERT( rsd.has("MULT") );
		Atom multi( rsd.atom( "MULT" ) );
		TS_ASSERT_EQUALS( multi.name(), "MULT" );
		TS_ASSERT_EQUALS( multi.charge(), -1.5 );
		TS_ASSERT_EQUALS( multi.element_type()->get_atomic_number(), 6 );
		TS_ASSERT_EQUALS( (*atom_types)[multi.atom_type_index()].name(), "aroC" );
		TS_ASSERT_EQUALS( (*mm_atom_types)[multi.mm_atom_type_index()].name(), "VIRT" );
	}

	void test_chi_assignment() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		utility::io::izstream paramslist("core/chemical/params/retype_list.txt");
		std::string filename("params/type_test1.params");
		paramslist >> filename;
		while( paramslist ) {
			TR << "------- Redoing chis for " << filename << std::endl;
			core::chemical::ResidueTypeOP rsd = read_topology_file("core/chemical/"+filename,
					atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));
			rsd->finalize();
			if( TR.Debug.visible() ) {
				print_chis(TR.Debug, *rsd);
			}
			core::chemical::ResidueTypeOP copy( new core::chemical::ResidueType( *rsd) );
			core::chemical::find_bonds_in_rings( *copy ); //Ring bond annotations needed
			copy->autodetermine_chi_bonds();
			copy->finalize();
			if( TR.Debug.visible() ) {
				TR.Debug << " Post chi redo: " << std::endl;
				print_chis(TR.Debug, *copy);
			}
			TS_ASSERT_EQUALS(rsd->nchi(), copy->nchi());
			// All the chis should be present, but they need not be in the same order.
			for( core::Size ii(1); ii <= rsd->nchi(); ++ii ) {
				AtomIndices const & rsd_chi( rsd->chi_atoms(ii) );
				bool chi_found=false;
				for( core::Size jj(1); jj <= copy->nchi(); ++jj ) {
					AtomIndices const & copy_chi( copy->chi_atoms(jj) );
					// Atom indicies should match
					if( rsd_chi[1] == copy_chi[1] && rsd_chi[2] == copy_chi[2] &&
							rsd_chi[3] == copy_chi[3] && rsd_chi[4] == copy_chi[4]) {
						TS_ASSERT_EQUALS(rsd->is_proton_chi(ii), copy->is_proton_chi(jj));
						chi_found=true;
						break;
					}
				}
				TS_ASSERT(chi_found);
				if( ! chi_found ) {
					TR << "---- Chis before" << std::endl;
					print_chis(TR, *rsd);
					TR << "---- Chis after" << std::endl;
					print_chis(TR, *copy);
				}
			}
			paramslist >> filename;
		} // while( paramslist )
		paramslist.close();
	} //  test_chi_assignment

	void test_icoor_reassignment() {
		// As matching up the icoor tree exactly is not necessarily required,
		// we instead check to make sure that the produced icoor table will produce
		// the same xyz coordinates that were fed in to produce it.
		using namespace core::chemical;
		core::Real delta2 = 0.0005*0.0005; // Half the resolution of the PDB format, squared

		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		utility::io::izstream paramslist("core/chemical/params/retype_list.txt");
		std::string filename;
		paramslist >> filename;
		while( paramslist ) {
			TR << "Rebuilding Icoord from xyz " << filename << std::endl;
			core::chemical::ResidueTypeOP restype = read_topology_file("core/chemical/"+filename,
					atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));
			core::chemical::ResidueTypeCOP restype_ref( core::chemical::ResidueTypeOP( new core::chemical::ResidueType(*restype) ) );
			restype->name( restype_ref->name() + "_IcoorRedo"); // For debugging purposes.

			restype->assign_internal_coordinates();


			restype->fill_ideal_xyz_from_icoor();

			//As we may not have picked the identical root (or atom tree ordering)
			//make sure we're oriented in the same reference frame.
			//Code based on Residue::orient_onto_residue()
			core::Size center, nbr1, nbr2;
			restype->select_orient_atoms( center, nbr1, nbr2 );

			core::Size  src_center = restype_ref->atom_index( restype_ref->atom_name( center ));
			core::Size  src_nbr1 = restype_ref->atom_index( restype_ref->atom_name( nbr1 ));
			core::Size  src_nbr2 = restype_ref->atom_index( restype_ref->atom_name( nbr2 ));
			using core::kinematics::Stub;
			core::Vector const
			rot_midpoint ( 0.5 * (     restype->atom(     nbr1 ).ideal_xyz() +     restype->atom(     nbr2 ).ideal_xyz() ) ),
			src_midpoint ( 0.5 * ( restype_ref->atom( src_nbr1 ).ideal_xyz() + restype_ref->atom( src_nbr2 ).ideal_xyz() ) );

			core::kinematics::Stub rot_stub( restype->atom( center ).ideal_xyz(),
					rot_midpoint,
					restype->atom( nbr1 ).ideal_xyz() );

			core::kinematics::Stub src_stub( restype_ref->atom( src_center ).ideal_xyz(),
					src_midpoint,
					restype_ref->atom( src_nbr1 ).ideal_xyz() );
			bool mistake=false;
			TS_ASSERT_EQUALS( restype_ref->natoms(), restype->natoms() );
			for ( core::Size i=1; i<= restype->natoms(); ++i ) {
				std::string atom_name(restype->atom_name(i));
				TS_ASSERT( restype_ref->has( atom_name ) );
				core::Vector const old_xyz( restype->atom(i).ideal_xyz() );
				core::Vector const new_xyz( src_stub.local2global( rot_stub.global2local( old_xyz ) ) );
				core::Vector const original_xyz( restype_ref->atom(atom_name).ideal_xyz() );
				TS_ASSERT_DELTA( original_xyz.distance_squared( new_xyz ), 0, delta2 );
				if( original_xyz.distance_squared( new_xyz ) > delta2 ) {
					TR << "Atom " << atom_name << " inappropriately placed " << new_xyz << std::endl;
					TR << "     " << "    "    << "              should be " << original_xyz << std::endl;
					TR << "     " << "    "    << "  (without translation) " << old_xyz << std::endl;
					TR << "Distance shift of " << original_xyz.distance( new_xyz ) << std::endl;
					mistake=true;
				}
				//restype2->atom(i).ideal_xyz( new_xyz );
			}
			if( mistake ) {
				core::chemical::sdf::MolWriter writer;
				writer.output_residue("core/chemical/"+ filename + "_orig.mol", core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( *restype_ref, true) ) ) );
				writer.output_residue("core/chemical/"+ filename + "_after.mol", restype );

				TR << "### For Restype Before " << std::endl;
				TR << "\n";
				// Here we assume that the neighbor atom is the ICOOR root.
				core::chemical::pretty_print_atomicoor(TR, restype_ref->icoor( 1 ), *restype_ref);
				TR << std::endl;
				TR << "### For Restype After" << std::endl;
				TR << "\n";
				core::chemical::pretty_print_atomicoor(TR, restype->icoor( restype->nbr_atom() ), *restype);
				TR << std::endl;
			}
			paramslist >> filename;
		}
	}

	void test_variant_matching()
	{
		using namespace basic::options;
		using namespace core::chemical;

		// Set up test ResidueTypes.
		ChemicalManager * manager( ChemicalManager::get_instance() );
		AtomTypeSetCOP atom_types = manager->atom_type_set( FA_STANDARD );
		ElementSetCOP element_types = manager->element_set( "default" );
		MMAtomTypeSetCOP mm_atom_types = manager->mm_atom_type_set( FA_STANDARD );
		OrbitalTypeSetCOP orbital_types = manager->orbital_type_set( FA_STANDARD );

		ResidueType res1( atom_types, element_types, mm_atom_types, orbital_types );
		res1.add_variant_type( UPPER_TERMINUS_VARIANT );
		res1.add_variant_type( LOWER_TERMINUS_VARIANT );

		ResidueType res2( res1 );

		ResidueType res3( res1 );
		res3.add_variant_type( PROTONATED );

		ResidueType res4( res1 );
		res4.add_variant_type( ADDUCT_VARIANT );

		ResidueType res5( res1 );
		res5.add_variant_type( REPLONLY );
		res5.add_variant_type( SPECIAL_ROT );

		// Check all three varieties of comparisons.
		TS_TRACE( "Testing variant comparisons using VariantType enums." );

		TS_ASSERT( variants_match( res1, res2 ) );
		TS_ASSERT( ! variants_match( res1, res3 ) );
		option[ OptionKeys::pH::pH_mode ]( true );
		TS_ASSERT( variants_match( res1, res3 ) );  // These should match in pH mode.
		TS_ASSERT( ! variants_match( res1, res4 ) );
		TS_ASSERT( ! variants_match( res1, res5 ) );

		TS_ASSERT( nonadduct_variants_match( res1, res2 ) );
		TS_ASSERT( ! nonadduct_variants_match( res1, res3 ) );
		TS_ASSERT( nonadduct_variants_match( res1, res4 ) );
		TS_ASSERT( ! nonadduct_variants_match( res1, res5 ) );

		vector1< VariantType > exceptions;
		exceptions.push_back( REPLONLY );
		exceptions.push_back( SPECIAL_ROT );
		TS_ASSERT( variants_match_with_exceptions( res1, res2, exceptions ) );
		TS_ASSERT( ! variants_match_with_exceptions( res1, res3, exceptions ) );
		TS_ASSERT( ! variants_match_with_exceptions( res1, res4, exceptions ) );
		TS_ASSERT( variants_match_with_exceptions( res1, res5, exceptions ) );

		// Check for comparisons when "custom" VariantTypes are involved.
		TS_TRACE( "Testing variant comparisons using custom VariantTypes." );

		res2.enable_custom_variant_types();
		TS_ASSERT( variants_match( res1, res2 ) );

		res2.add_variant_type( "FUNKY" );
		TS_ASSERT( ! variants_match( res1, res2 ) );

		res3.enable_custom_variant_types();
		TS_ASSERT( ! variants_match( res2, res3 ) );

		res3.add_variant_type( "TUBULAR" );
		TS_ASSERT( ! variants_match( res2, res3 ) );

		res2.add_variant_type( "TUBULAR" );
		res3.add_variant_type( "FUNKY" );
		TS_ASSERT( variants_match( res2, res3 ) );
		TS_ASSERT( ! nonadduct_variants_match( res2, res3 ) );
	}

	//test functions in bond support and ring detection. If this fails, those functions were altered
	void test_odd_proline_ring_detection(){
		ChemicalManager * cm(ChemicalManager::get_instance());
		string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

		ResidueType rsd( atom_types, element_types, mm_atom_types, orbital_types);

		// Build a  weird "Proline" no hydrogens
		add_atom( rsd, atom_types," N  ", "Npro", "NH1", -0.47);
		add_atom( rsd, atom_types," CA ", "CAbb", "CT1", 0.07);
		add_atom( rsd, atom_types," C  ", "CObb", "C", 0.51);
		add_atom( rsd, atom_types," O  ", "OCbb", "O", -0.51);
		add_atom( rsd, atom_types," CB ", "CH2", "CP2", -0.27);
		add_atom( rsd, atom_types," CG ", "CH2", "CP2", -0.27);
		add_atom( rsd, atom_types," CD ", "CH2", "CP3", -0.27);
		add_atom(rsd, atom_types, " CN ", "CH2", "CP3", -0.27);

		add_bond( rsd, "N", "CA");
		add_bond( rsd, "N", "CD");
		add_bond( rsd, "N", "CN");
		add_bond( rsd, "CN", "CD");
		add_bond( rsd, "CA", "C");
		add_bond( rsd, "CA", "CB");
		add_bond( rsd, "C", "O");
		add_bond( rsd, "CB", "CG");
		add_bond( rsd, "CG", "CD");


		//Do we need icoors? Not for this test!
		rsd.finalize(); //gotta finalize to work on it
		core::chemical::find_bonds_in_rings(rsd);
		core::chemical::ResidueGraph const full_residue_graph = rsd.graph();


		core::chemical::VD  source = rsd.atom_vertex( rsd.atom_index(" CA ") );
		core::chemical::VD  target = rsd.atom_vertex( rsd.atom_index("C") );

		core::chemical::ED edge_ca_c = core::chemical::get_bond(rsd, source, target);
		core::chemical::VD target_cb = rsd.atom_vertex( rsd.atom_index(" CB ") );

		core::chemical::ED edge_ca_cb = core::chemical::get_bond(rsd, source, target_cb);

		core::chemical::Bond const & bond_ca_c = rsd.bond(edge_ca_c);
		core::chemical::Bond const & bond_ca_cb = rsd.bond(edge_ca_cb);


		TS_ASSERT_EQUALS(bond_ca_c.ringness(), core::chemical::BondNotInRing);
		TS_ASSERT_EQUALS(bond_ca_cb.ringness(), core::chemical::BondInRing);


		utility::vector1<core::chemical::VD> check_that_get_bond_and_get_connecting_atoms_works = core::chemical::get_connecting_atoms(rsd, edge_ca_c);

		TS_ASSERT_EQUALS(check_that_get_bond_and_get_connecting_atoms_works[1], source);
		TS_ASSERT_EQUALS(check_that_get_bond_and_get_connecting_atoms_works[2], target);


	}

};
