// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RDMolToRestypeTests.cxxtest.hh
/// @brief unit tests for RDMolToRestype functionality
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/chemical/rdkit/RDMolToRestype.hh>

// Project Headers

#include <core/chemical/ResidueGraphTypes.hh>

#include <core/chemical/sdf/mol_writer.hh>

#include <basic/Tracer.hh>

// Platform Headers

// C++ Headers

// External library headers

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>

#include <core/chemical/Element.hh> // AUTO IWYU For Element

static basic::Tracer TR("core.chemical.rdkit.RDMolToRestypeTests.cxxtest");

class RDMolToRestypeTests : public CxxTest::TestSuite {

public:


	void setUp() {
		core_init();

	}

	void tearDown() {}

	void test_convert() {
		TR << "In test_convert." << std::endl;

		using namespace core::chemical;

		std::string smiles("CC(=O)Oc1ccccc1C(=O)O" ); //Aspirin

		// We want to take ownership of the returned molecule
		::RDKit::RWMolOP rdmol(  ::RDKit::SmilesToMol( smiles ) );
		TS_ASSERT( rdmol.get() != 0 );
		::RDKit::MolOps::addHs( *rdmol );
		TS_ASSERT( rdmol.get() != 0 );
		::RDKit::DGeomHelpers::EmbedMolecule( *rdmol ); // Generate crude 3D coordinates

		core::chemical::rdkit::RDMolToRestype converter( *rdmol );

		core::chemical::MutableResidueTypeOP restype( converter.generate_restype() );

		VDIndexMapping const & mapping( converter.index_to_vd().reverse() );

		TS_ASSERT_EQUALS( restype->natoms(), rdmol->getNumAtoms() );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), rdmol->getNumHeavyAtoms() );

		TR << "Atoms  " <<  restype->natoms() << " " << rdmol->getNumAtoms() << std::endl;

		for ( VIterPair itr( restype->atom_iterators() ); itr.first != itr.second; ++itr.first ) {
			if ( mapping[ *itr.first ] == mapping.invalid_entry() ) { continue; }
			TS_ASSERT_EQUALS( restype->atom( *itr.first ).element_type()->get_atomic_number(),
				(core::Size) rdmol->getAtomWithIdx( mapping[*itr.first]  )->getAtomicNum() );
			TR << "Atom " << restype->atom( *itr.first ).element_type()->get_atomic_number() << " " << (core::Size) rdmol->getAtomWithIdx( mapping[*itr.first]  )->getAtomicNum() << std::endl;
		}

		for ( EIterPair itr( restype->bond_iterators() ); itr.first != itr.second; ++itr.first ) {
			VD source( boost::source( *itr.first, restype->graph() ) );
			VD target( boost::target( *itr.first, restype->graph() ) );
			if ( mapping[ source ] == mapping.invalid_entry() ) { continue; }
			if ( mapping[ target ] == mapping.invalid_entry() ) { continue; }
			::RDKit::Bond const * bond( rdmol->getBondBetweenAtoms( mapping[source], mapping[target] ) );
			TS_ASSERT( bond != 0 );
			double bondtype( bond->getBondTypeAsDouble() );
			if ( 1.4 <= bondtype && bondtype <= 1.6 ) {
				TS_ASSERT_EQUALS( restype->bond( *itr.first ).bond_name(), AromaticBond );
				TR << "Bond ARO" << std::endl;
			} else {
				TS_ASSERT_EQUALS( double( restype->bond( *itr.first ).bond_name() ), bondtype );
				TR <<  "Bond " << restype->bond( *itr.first ).bond_name() << " " << bondtype << std::endl;
			}
		}

		// Dump to tracer for debugging purposes.
		core::chemical::sdf::MolWriter writer;
		TR.Debug << "\n$$$$$\n";
		writer.output_residue( TR.Debug, restype );
		TR.Debug << std::endl;

	}

};
