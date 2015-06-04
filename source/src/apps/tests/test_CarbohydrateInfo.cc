// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    apps/tests/test_CarbohydrateInfo.cc
/// @brief   Application source code for testing CarbohydrateInfo.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/RingConformerSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>

// C++ headers
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace pose;
using namespace import_pose;
using namespace chemical;
using namespace gasteiger;
using namespace conformation;


string const PATH = "input/";


void
test_sugar( Pose const & sugar )
{
	cout << endl << sugar << endl;

	cout << "Sequences:" << endl;
	Size const n_chains( sugar.conformation().num_chains() );
	for ( core::uint i = 1; i <= n_chains; ++i ) {
		cout << " Chain " << i << ": ";
		cout << sugar.chain_sequence( i ) << endl;
	}
	file::FileName filename( sugar.pdb_info()->name() );
	filename.path( "output/" );
	cout << "Writing file: " << filename.ext( "gws" ) << endl;
	io::carbohydrates::dump_gws( sugar, filename.ext( "gws" ) );

	cout << endl << "Residue Info:" << endl;
	Size const n_res( sugar.total_residue() );
	for ( core::uint i = 1; i <= n_res; ++i ) {
		Residue const & res( sugar.residue( i ) );
		cout << "PDB ID: " << sugar.pdb_info()->pose2pdb( i ) << ": ";
		// BEGIN_HACK
		//res.show( cout, true );  // Show verbose output, including atomic details.

		// I would like to make a flag, -assign_gasteiger_atom_types, that does this on creation of the ResidueTypeSet,
		// but in the meantime, make a copy of the ResidueType and assign the GasteigerAtomTypes for testing.
		res.show( cout );
		ResidueType copy_of_type( res.type() );
		GasteigerAtomTypeSetCOP gasteiger_set( ChemicalManager::get_instance()->gasteiger_atom_type_set() );
		assign_gasteiger_atom_types( copy_of_type, gasteiger_set, false );
		cout << " Atomic Details:" << endl;
		Size n_atoms = copy_of_type.natoms();
		for ( core::uint j = 1; j <= n_atoms; ++j ) {
			cout << "  Atom " << j << ": ";
			copy_of_type.atom( j ).show( cout );
		}
		// END_HACK

		cout << endl << endl;
	}
}


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose maltotriose, isomaltose, lactose, amylopectin, glycopeptide, glucosamine, N_linked_14_mer, psicose,
				neuraminate, Lex, GalCer, Me_glycoside, maltobiose;
		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );


		cout << "Importing maltotriose:" << endl;

		pose_from_pdb( maltotriose, PATH + "maltotriose.pdb" );

		test_sugar( maltotriose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing isomaltose:" << endl;

		pose_from_pdb( isomaltose, PATH + "isomaltose.pdb" );

		test_sugar( isomaltose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing lactose:" << endl;

		pose_from_pdb( lactose, PATH + "lactose.pdb" );

		test_sugar( lactose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating maltotriose from sequence: alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp:" << endl;

		make_pose_from_saccharide_sequence(
				maltotriose, "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp", *residue_set );

		// TODO: Figure out what is wrong from sugar poses made from sequences.
		// Note to self: It's a PDBInfo issue.
		//test_sugar(maltotriose);
		cout << endl << maltotriose << endl;

		cout << "Sequences:" << endl;
		cout << " Chain 1: ";
		cout << maltotriose.chain_sequence( 1 ) << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing branched amylopectin fragment:" << endl;

		pose_from_pdb( amylopectin, PATH + "amylopectin_fragment.pdb" );

		test_sugar( amylopectin );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing N-glycosylated sample:" << endl;

		pose_from_pdb( glycopeptide, PATH + "glycosylated_peptide.pdb" );

		test_sugar( glycopeptide );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing glucosamine:" << endl;

		pose_from_pdb( glucosamine, PATH + "GlcN.pdb" );

		test_sugar( glucosamine );

		cout << *glucosamine.residue( 1 ).type().ring_conformer_set();


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating glucosamine from sequence: GlcpN:" << endl;

		make_pose_from_saccharide_sequence( glucosamine, "GlcpN", *residue_set );

		cout << endl << glucosamine << endl;

		cout << "Sequences:" << endl;
		cout << " Chain 1: ";
		cout << glucosamine.chain_sequence( 1 ) << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing N-glycosylated 14-mer:" << endl;

		pose_from_pdb( N_linked_14_mer, PATH + "N-linked_14-mer_glycan.pdb" );

		test_sugar( N_linked_14_mer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing beta-D-psicopyranose:" << endl;

		pose_from_pdb(psicose, PATH + "beta-psicose.pdb");

		test_sugar( psicose );

		cout << *psicose.residue( 1 ).type().ring_conformer_set();


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing alpha-neuraminic acid:" << endl;

		pose_from_pdb( neuraminate, PATH + "Neu.pdb" );

		test_sugar( neuraminate );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing Lewisx:" << endl;

		pose_from_pdb( Lex, PATH + "Lex.pdb" );

		test_sugar( Lex );
		
		
		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating Lewisx from sequence: beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpN:" << endl;

		//make_pose_from_saccharide_sequence( Lex, "beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpN", *residue_set );
		// TEMP
		vector1< ResidueTypeCOP > residue_types(
				residue_types_from_saccharide_sequence( "beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpN", *residue_set ) );

		for ( core::uint i( 1 ); i <= residue_types.size();  ++i ) {
			cout << residue_types[ i ]->name() << endl;
		}
		// END TEMP

		cout << endl << Lex << endl;

		cout << "Sequences:" << endl;
		cout << " Chain 1: ";
		cout << Lex.chain_sequence( 1 ) << endl;
		cout << " Chain 2: ";
		cout << Lex.chain_sequence( 2 ) << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing Lewisx with alternate GlcNAc atom names:" << endl;

		pose_from_pdb( Lex, PATH + "Lex_alternate_names.pdb" );

		test_sugar( Lex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing maltobiose using GLYCAM residue names:" << endl;

		pose_from_pdb( maltobiose, PATH + "GLYCAM_maltobiose.pdb" );

		test_sugar( maltobiose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing GalCer:" << endl;

		pose_from_pdb( GalCer, PATH + "GalCer.pdb" );

		test_sugar( GalCer );

		/*// FIXME
		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing a sample methyl glycoside:" << endl;

		pose_from_pdb( Me_glycoside, PATH + "Me_glycoside.pdb" );

		test_sugar( Me_glycoside );
		*/

	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
