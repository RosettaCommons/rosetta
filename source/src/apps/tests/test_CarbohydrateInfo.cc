// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/tests/test_CarbohydrateInfo.cc
/// @brief   Application source code for testing CarbohydrateInfo.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace pose;
using namespace import_pose;
using namespace chemical;
using namespace conformation;


// Constants
int const SUCCESS( 0 );
int const FAILURE( -1 );

string const PATH( "input/" );


void
test_sugar( Pose const & sugar )
{
	cout << endl << sugar << endl;

	cout << "Sequences:" << endl;
	Size const n_chains( sugar.conformation().num_chains() );
	for ( core::uint i( 1 ); i <= n_chains; ++i ) {
		cout << " Chain " << i << ": ";
		cout << sugar.chain_sequence( i ) << endl;
	}
	file::FileName filename( sugar.pdb_info()->name() );
	filename.path( "output/" );
	cout << "Writing file: " << filename.ext( "gws" ) << endl;
	core::io::carbohydrates::dump_gws( sugar, filename.ext( "gws" ) );
	cout << "Writing file: " << filename.ext( "pdb" ) << endl;
	sugar.dump_pdb( filename.ext( "pdb" ) );

	cout << endl << "Residue Info:" << endl;
	Size const n_res( sugar.size() );
	for ( core::uint i( 1 ); i <= n_res; ++i ) {
		Residue const & res( sugar.residue( i ) );
		cout << "PDB ID: " << sugar.pdb_info()->pose2pdb( i ) << ": ";
		res.show( cout, true );  // Show verbose output, including atomic details.
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
		Pose maltotriose, isomaltose, lactose, amylopectin, glycopeptide, glucosamine, N_linked_14_mer, fluoro_sugar, free_14_mer,
			O_linked, psicose, glucuronic_acid, neuraminate, bacillosamine, Murp, Rhof, Lex, SLex, GalCer, UDP_D_Glc,
			target57, maltobiose, Me_glycoside, Me_glycoside_sequence, Me_glycoside_3mer, C_linked, Ac_sugar,
			ketopentofuranose, ketohexofuranose, Kdo, Kdn, whacky_sugar, pdb_code_pdb, bad_pdb, lactyl_sugar, ligand_sugar, phosphorylated_sugar;

		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );


		cout << "Importing maltotriose:" << endl;

		pose_from_file( maltotriose, PATH + "maltotriose.pdb", PDB_file );

		test_sugar( maltotriose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing isomaltose:" << endl;

		pose_from_file( isomaltose, PATH + "isomaltose.pdb", PDB_file );

		test_sugar( isomaltose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing lactose (with a wonky 4H5 ring):" << endl;

		pose_from_file( lactose, PATH + "lactose.pdb", PDB_file );

		test_sugar( lactose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating maltotriose from sequence: alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp:" << endl;

		make_pose_from_saccharide_sequence(
			maltotriose, "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp", *residue_set );
		maltotriose.pdb_info()->name( "maltotriose" );

		test_sugar( maltotriose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing branched amylopectin fragment:" << endl;

		pose_from_file( amylopectin, PATH + "amylopectin_fragment.pdb", PDB_file );

		test_sugar( amylopectin );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing glucosamine:" << endl;

		pose_from_file( glucosamine, PATH + "GlcN.pdb", PDB_file );

		test_sugar( glucosamine );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating glucosamine from sequence: GlcpN:" << endl;

		make_pose_from_saccharide_sequence( glucosamine, "GlcpN", *residue_set );
		glucosamine.pdb_info()->name( "glucosamine" );

		test_sugar( glucosamine );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing N-glycosylated 14-mer:" << endl;

		pose_from_file( N_linked_14_mer, PATH + "N-linked_14-mer_glycan.pdb", PDB_file );

		test_sugar( N_linked_14_mer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating free 14-mer from sequence: " <<
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
			"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-" << endl;

		make_pose_from_saccharide_sequence( free_14_mer,
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
			"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-", *residue_set );

		test_sugar( free_14_mer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating O-linked glycan from sequences: ASA and Glcp" << endl;

		make_pose_from_sequence( O_linked, "ASA", *residue_set );
		pose::carbohydrates::glycosylate_pose( O_linked, 2, "Glcp" );

		test_sugar( O_linked );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating core 5 O-linked glycan by glycosylating sequence: ASA" << endl;

		make_pose_from_sequence( O_linked, "ASA", *residue_set );
		pose::carbohydrates::glycosylate_pose_by_file( O_linked, 2, "core_5_O-glycan" );

		test_sugar( O_linked );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating N-linked, branched glycan from sequences:" << endl;

		make_pose_from_sequence( N_linked_14_mer, "ANASA", *residue_set );
		pose::carbohydrates::glycosylate_pose( N_linked_14_mer, 2,
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
			"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-" );

		test_sugar( N_linked_14_mer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing beta-D-psicopyranose:" << endl;

		pose_from_file(psicose, PATH + "beta-psicose.pdb", PDB_file );

		test_sugar( psicose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing beta-D-glucopyranuronic acid:" << endl;

		pose_from_file( glucuronic_acid, PATH + "glucuronic_acid.pdb", PDB_file );

		test_sugar( glucuronic_acid );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing alpha-neuraminic acid:" << endl;

		pose_from_file( neuraminate, PATH + "Neu.pdb", PDB_file );

		test_sugar( neuraminate );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating 2,4-diacetylbacillosamine:" << endl;

		make_pose_from_saccharide_sequence( bacillosamine, "->3)-beta-Bacp2,4Ac", *residue_set );

		test_sugar( bacillosamine );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing Lewisx:" << endl;

		pose_from_file( Lex, PATH + "Lex.pdb", PDB_file );

		test_sugar( Lex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating Lewisx from sequence: beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc:" << endl;

		make_pose_from_saccharide_sequence( Lex, "beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc", *residue_set );

		test_sugar( Lex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating muramic acid from sequence: alpha-Murp:" << endl;

		make_pose_from_saccharide_sequence( Murp, "alpha-Murp", *residue_set );

		test_sugar( Murp );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating rhodinofuranose from sequence: alpha-D-Rhof:" << endl;

		make_pose_from_saccharide_sequence( Rhof, "->5)-alpha-D-Rhof", *residue_set );

		test_sugar( Rhof );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing Lewisx with alternate GlcNAc atom names:" << endl;

		pose_from_file( Lex, PATH + "Lex_alternate_names.pdb", PDB_file );

		test_sugar( Lex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing sialyl-Lewisx:" << endl;

		pose_from_file( SLex, PATH + "SLex.pdb", PDB_file );

		test_sugar( SLex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating sialyl-Lewisx from sequence: " <<
			"->8)-alpha-Neup5Ac-(2->6)-beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc:" << endl;

		make_pose_from_saccharide_sequence(
			SLex, "->8)-alpha-Neup5Ac-(2->6)-beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc", *residue_set );

		test_sugar( SLex );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing maltobiose using GLYCAM residue names:" << endl;

		pose_from_file( maltobiose, PATH + "GLYCAM_maltobiose.pdb", PDB_file );

		test_sugar( maltobiose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing GalCer:" << endl;

		pose_from_file( GalCer, PATH + "GalCer.pdb", PDB_file );

		test_sugar( GalCer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing UDP-D-Glc:" << endl;

		pose_from_file( UDP_D_Glc, PATH + "UDP-alpha-D-Glcp.pdb", PDB_file );

		test_sugar( UDP_D_Glc );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating UDP-D-Glc through glycosylation:" << endl;

		make_pose_from_sequence( UDP_D_Glc, "Z[UDP:non-conjugated]", *residue_set );
		pose::carbohydrates::glycosylate_pose( UDP_D_Glc, 1, "a-D-Glcp-" );

		test_sugar( UDP_D_Glc );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing CAPRI Round 27 Target 57, a crazy, heparin-like hexamer with uronic acids, sulfates, and "
			"sulfonamidos:" << endl;

		pose_from_file( target57, PATH + "target57.pdb", PDB_file );

		test_sugar( target57 );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing a sample methyl glycoside:" << endl;

		pose_from_file( Me_glycoside, PATH + "Me_glycoside.pdb", PDB_file );

		test_sugar( Me_glycoside );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating a single methyl glycoside from sequence alpha-D-ManpOMe:" << endl;

		make_pose_from_saccharide_sequence( Me_glycoside_sequence, "alpha-D-ManpOMe", *residue_set );

		test_sugar( Me_glycoside_sequence );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating a 3-mer branched chain with methyl glycoside from sequence" <<
			" alpha-D-Galp-(1->2)-[alpha-D-Abep-(1->3)]-alpha-D-ManpOMe" <<
			" (based on PDB 1MFA):" << endl;

		make_pose_from_saccharide_sequence( Me_glycoside_3mer,
			"alpha-D-Galp-(1->2)-[alpha-D-Abep-(1->3)]-alpha-D-ManpOMe", *residue_set );

		test_sugar( Me_glycoside_3mer );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating C-linked glycan from sequences: AWAAWA and alpha-D-Manp" << endl;

		make_pose_from_sequence( C_linked, "AWAAWA", *residue_set );
		pose::carbohydrates::glycosylate_pose( C_linked, 2, "a-D-Manp-" );

		test_sugar( C_linked );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating aldohexofuranoses from sequences" << endl;

		// TODO: Fix defaults so that people don't have to keep specifying the non-reducing-end linkage. ~JWL
		make_pose_from_saccharide_sequence( ketohexofuranose, "->3)-a-D-Glcf-(1->3)-b-D-Idof", *residue_set );

		test_sugar( ketohexofuranose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating an O-acetylated sugar from sequence:" << endl;

		make_pose_from_saccharide_sequence( Ac_sugar, "alpha-d-Glcp6Ac", *residue_set );

		test_sugar( Ac_sugar );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating aldopentofuranose dimer from sequence:" << endl;

		make_pose_from_saccharide_sequence( ketopentofuranose, "->3)-a-D-Xulf-(2->1)-b-D-Rulf", *residue_set );

		test_sugar( ketopentofuranose );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating ketodeoxyoctanoate from sequence:" << endl;

		make_pose_from_saccharide_sequence( Kdo, "->8)-a-Kdop", *residue_set );

		test_sugar( Kdo );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating ketodeoxynonanoate from sequence:" << endl;

		make_pose_from_saccharide_sequence( Kdn, "->9)-a-Kdnp", *residue_set );

		test_sugar( Kdn );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating whacky sugar from sequence to really stretch the system:" << endl;

		make_pose_from_saccharide_sequence( whacky_sugar,
			"b-D-Fruf-(2->8)-a-Neup5Ac-(2->4)-b-D-GlcpNS6S-(1->4)-[a-D-Xylp-(1->3)]-b-L-GulpA-(1->5)-b-D-Psip",
			*residue_set );

		test_sugar( whacky_sugar );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing a .pdb file using PDB 3-letter codes, "
			"including one that cannot have position 3 as the default main-chain connection:" << endl;

		pose_from_file( pdb_code_pdb, PATH + "pdb_w_pdb_3_letter_codes.pdb", PDB_file );

		test_sugar( pdb_code_pdb );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating lactyl sugar from sequence:" << endl;

		make_pose_from_saccharide_sequence( lactyl_sugar, "->4)-alpha-d-Glcp3Lac", *residue_set );
		lactyl_sugar.pdb_info()->name( "Glcp3Lac" );

		test_sugar( lactyl_sugar );


		cout << "---------------------------------------------------------------------------------------------" << endl;

		cout << "Creating fluoro sugar from sequence:" << endl;

		make_pose_from_saccharide_sequence( fluoro_sugar, "->3)-alpha-d-Glcp2F", *residue_set );
		fluoro_sugar.pdb_info()->name( "Glcp2F" );

		test_sugar( fluoro_sugar );


		cout << "Creating phosphorylated sugar from sequence:" << endl;

		make_pose_from_saccharide_sequence( phosphorylated_sugar, "->4)-alpha-d-Glcp6P", *residue_set );
		phosphorylated_sugar.pdb_info()->name( "Glcp6P" );

		test_sugar( phosphorylated_sugar );


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing a .pdb file with bad LINK records:" << endl;

		pose_from_file( bad_pdb, PATH + "pdb_w_bad_links.pdb", PDB_file );

		cout << ".pdb file with bad LINK records imported successfully." << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating a non-polymer linear monosaccharide:" << endl;

		make_pose_from_sequence( ligand_sugar, "X[D-Gly]", *residue_set );
		ligand_sugar.pdb_info( pointer::make_shared< PDBInfo >( ligand_sugar ) );
		ligand_sugar.pdb_info()->name( "D-Gly" );

		test_sugar( ligand_sugar );

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return FAILURE;
	}
	return SUCCESS;
}
