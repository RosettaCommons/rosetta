// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RestypeToRDMolTests.cxxtest.hh
/// @brief unit tests for RestypeToRDMol functionality
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers

#include <core/chemical/rdkit/RestypeToRDMol.hh>

// Project Headers

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

#include <basic/Tracer.hh>

// C++ Headers

// External library headers

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

static basic::Tracer TR("core.chemical.rdkit.RestypeToRDMolTests");

class RestypeToRDMolTests : public CxxTest::TestSuite {

	core::chemical::MutableResidueTypeOP restype_;

public:


	void setUp() {
		core_init();

		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP rsd_types( cm->residue_type_set(FA_STANDARD) );

		restype_ = read_topology_file( "core/chemical/params/U01.params", rsd_types );
	}

	void tearDown() {}

	void test_convert() {
		TR << "In test_convert." << std::endl;

		using namespace core::chemical;

		core::chemical::rdkit::RestypeToRDMol converter( *restype_ );

		::RDKit::RWMOL_SPTR rdmol( converter.Mol() );
		VDIndexMapping const & mapping( converter.vd_to_index() );

		// By default, the converter removes hydrogens.
		TS_ASSERT_EQUALS( restype_->nheavyatoms(), rdmol->getNumAtoms() );
		TS_ASSERT_EQUALS( restype_->nheavyatoms(), rdmol->getNumHeavyAtoms() );

		TR << "Atoms  " <<  restype_->natoms() << " " << rdmol->getNumAtoms() << std::endl;

		for ( VIterPair itr( restype_->atom_iterators() ); itr.first != itr.second; ++itr.first ) {
			if ( mapping[ *itr.first ] == mapping.invalid_entry() ) { continue; }
			TS_ASSERT_EQUALS( restype_->atom( *itr.first ).element_type()->get_atomic_number(),
				(core::Size) rdmol->getAtomWithIdx( mapping[*itr.first]  )->getAtomicNum() );
			TR << "Atom " << restype_->atom( *itr.first ).element_type()->get_atomic_number() << " " << (core::Size) rdmol->getAtomWithIdx( mapping[*itr.first]  )->getAtomicNum() << std::endl;
		}

		for ( EIterPair itr( restype_->bond_iterators() ); itr.first != itr.second; ++itr.first ) {
			VD source( boost::source( *itr.first, restype_->graph() ) );
			VD target( boost::target( *itr.first, restype_->graph() ) );
			if ( mapping[ source ] == mapping.invalid_entry() ) { continue; }
			if ( mapping[ target ] == mapping.invalid_entry() ) { continue; }
			::RDKit::Bond const * bond( rdmol->getBondBetweenAtoms( mapping[source], mapping[target] ) );
			TS_ASSERT( bond != 0 );
			double bondtype( bond->getBondTypeAsDouble() );
			if ( 1.4 <= bondtype && bondtype <= 1.6 ) {
				// RDKit detects aromaticity, though test structure is kekulized
				//TS_ASSERT_EQUALS( restype_->bond( *itr.first ).bond_name(), AromaticBond );
				TR << "Bond ARO" << std::endl;
			} else {
				TS_ASSERT_EQUALS( double( restype_->bond( *itr.first ).bond_name() ), bondtype );
				TR <<  "Bond " << restype_->bond( *itr.first ).bond_name() << " " << bondtype << std::endl;
			}
		}
	}

	// Test if we can adequately deal with charged species and the like
	// Mainly, this is testing that we don't crash on these.
	void test_charge_states() {
		using namespace core::chemical;

		utility::vector1< std::string > filenames;
		filenames.push_back( "purine.sdf" );
		filenames.push_back( "pyridineNoxide.sdf" );
		filenames.push_back( "nmethylpyridine.sdf" );
		filenames.push_back( "benzoate.sdf" );
		filenames.push_back( "ATP.sdf" );
		filenames.push_back( "nitro.sdf" );

		for ( core::Size ii(1); ii <= filenames.size(); ++ii ) {
			TR << "Testing file: " << filenames[ii] << std::endl;
			utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("core/chemical/rdkit/"+filenames[ii]) );
			TS_ASSERT( restypes.size() > 0 );
			core::chemical::MutableResidueTypeOP restype( restypes[1] );

			core::chemical::rdkit::RestypeToRDMol converter( *restype, /* neutralize= */ false );
			::RDKit::RWMOL_SPTR rdmol( converter.Mol() );
			TS_ASSERT( rdmol );
			TS_ASSERT_THROWS_NOTHING( ::RDKit::MolOps::sanitizeMol(*rdmol) );
			// That's it - just make sure that things convert properly and gives a sanitizable molecule.

			core::chemical::rdkit::RestypeToRDMol converter_nochg( *restype, /* neutralize= */ true );
			::RDKit::RWMOL_SPTR rdmol_nochg( converter_nochg.Mol() );
			TS_ASSERT( rdmol_nochg );
			TS_ASSERT_THROWS_NOTHING( ::RDKit::MolOps::sanitizeMol(*rdmol_nochg) );
			// That's it - just make sure that things convert properly and gives a sanitizable molecule.
		}
	}

	// Test if we deal with normalize charges appropriately.
	void test_charged_reprotonation() {

		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP rts( cm->residue_type_set(FA_STANDARD));

		{
			TR << "Testing benzoate reprotonation." << std::endl;
			utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("core/chemical/rdkit/benzoate.sdf") );
			TS_ASSERT( restypes.size() > 0 );
			core::chemical::MutableResidueTypeOP benzoate( restypes[1] );

			core::chemical::rdkit::RestypeToRDMol convert_benz( *benzoate, /* neutralize= */ true );
			::RDKit::RWMOL_SPTR benz( convert_benz.Mol() );
			TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *benz ), "O=C(O)c1ccccc1" );
		}

		// // Glutamate doesn't work, as the "aromatic" bonds on the sidechain aren't understood by RDKit.
		//core::chemical::rdkit::RestypeToRDMol convert_glu( rts->name_map("GLU"), /* normalize= */ true );
		//::RDKit::RWMOL_SPTR glu( convert_glu.Mol() );
		//TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *glu ), "" );

		TR << "Testing lysine reprotonation." << std::endl;
		core::chemical::MutableResidueTypeOP lys_rt( utility::pointer::make_shared< core::chemical::MutableResidueType >(rts->name_map("LYS")) );
		core::chemical::rdkit::RestypeToRDMol convert_lys( *lys_rt, /* neutralize= */ true );
		::RDKit::RWMOL_SPTR lys( convert_lys.Mol() );
		TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *lys, /*doIsomericSmiles*/ false ), "NCCCCC(N)C=O" );

		{
			TR << "Testing ATP reprotonation." << std::endl;
			utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("core/chemical/rdkit/ATP.sdf") );
			TS_ASSERT( restypes.size() > 0 );
			core::chemical::MutableResidueTypeOP ATP( restypes[1] );
			core::chemical::rdkit::RestypeToRDMol convert_atp( *ATP, /* neutralize= */ true );

			::RDKit::RWMOL_SPTR atp( convert_atp.Mol() );
			TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *atp, /*doIsomericSmiles*/ false ), "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O" );
		}

		{
			TR << "Testing Nitro reprotonation." << std::endl;
			utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("core/chemical/rdkit/nitro.sdf") );
			TS_ASSERT( restypes.size() > 0 );
			core::chemical::MutableResidueTypeOP nitro( restypes[1] );
			core::chemical::rdkit::RestypeToRDMol convert_nitro( *nitro, /* neutralize= */ true );

			::RDKit::RWMOL_SPTR nitro_mol( convert_nitro.Mol() );
			TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *nitro_mol ), "Nc1cc(C(=O)O)cc([N+](=O)[O-])c1" );
		}

		// // Non-Ring atom marked aromatic
		// core::chemical::rdkit::RestypeToRDMol convert_gua( rts->name_map("GUA"), /* normalize= */ true );
		// ::RDKit::RWMOL_SPTR gua( convert_gua.Mol() );
		// TS_ASSERT_EQUALS( ::RDKit::MolToSmiles( *gua ), "" );
	}
};
