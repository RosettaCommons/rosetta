// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmcif/util.cc
/// @brief Functions for MMCIF writing.
/// @author Andy Watkins (andy.watkins2@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Unit Header
#include <core/io/mmcif/cif_writer.hh>

// Package headers
#include <core/io/rcsb/ExperimentalTechnique.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// When you move PDBReader and PoseUnbuilder, take these.
#include <core/pose/Pose.fwd.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>


#include <core/io/Remarks.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/gemmi_util.hh>
#include <gemmi/to_cif.hpp>

// Numeric headers

// External headers
#include <gemmi/cif.hpp>

#include <core/io/AtomInformation.hh> // AUTO IWYU For AtomInformation

// C++ headers

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.io.mmcif.cif_writer.hh" );

namespace core {
namespace io {
namespace mmcif {

using utility::gemmi_add_table;

using utility::to_string;

/// @brief Dump an mmCIF from a pose to a file.
///  Return success or failure.
bool
dump_cif(
	core::pose::Pose const & pose,
	std::string const & file_name,
	StructFileRepOptionsCOP options )

{

	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );

	converter.init_from_pose( pose );
	StructFileRepOP cif_sfr =  converter.sfr();

	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "cif_writer:dump_cif Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_cif( file, cif_sfr, *options);
	file.close();
	return true;
}

/// @brief Dump an mmCIF from a pose to an ostream.
StructFileRepOP
dump_cif(
	core::pose::Pose const & pose,
	std::ostream & out,
	StructFileRepOptionsCOP options)

{

	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );

	converter.init_from_pose( pose );
	StructFileRepOP cif_sfr =  converter.sfr();


	dump_cif( out, cif_sfr, *options);
	return cif_sfr;
}

///@brief Return an mmCIF-format string from a pose.
std::string
dump_cif( core::pose::Pose const & pose )
{
	std::ostringstream out;
	dump_cif( pose, out );
	return out.str();
}

/// @brief Dump an mmCIF from a pose, optionally extracting extra info.
StructFileRepOP
dump_cif(
	core::pose::Pose const & pose,
	std::string const & jd2_job_data,
	std::ostream & out)

{
	StructFileRepOptionsOP options =  utility::pointer::make_shared< StructFileRepOptions >();
	core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter();

	converter.init_from_pose( pose );
	StructFileRepOP cif_sfr =  converter.sfr();
	cif_sfr->append_to_additional_string_output( jd2_job_data );

	dump_cif( out, cif_sfr, *options);
	return cif_sfr;

}


bool
dump_cif(
	std::string const & file_name,
	StructFileRepOP sfr,
	StructFileRepOptions const & options
) {

	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "cif_writer:dump_cif Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}

	dump_cif(file, sfr, options);
	file.close();
	return true;
}

/// @brief Main dump_cif function.
///  Return a string with contents of a CIF file, extracted from a Pose to a StructFileRep via PoseToStructFileRepConverter
std::string
dump_cif(
	StructFileRepOP sfr,
	StructFileRepOptions const & options )

{
	std::ostringstream out;
	dump_cif( out, sfr, options );
	return out.str();
}


void
dump_cif(
	std::ostream & out,
	StructFileRepOP sfr,
	StructFileRepOptions const & options
) {
	////// NOTES:
	// It's important to wrap all strings being passed with gemmi:cif::quote()
	// This only will wrap things in quotes if it needs it.
	// The utility function gemmi_add_row() will do this for you automatically.
	using utility::gemmi_add_row;

	gemmi::cif::Document cifdoc;
	gemmi::cif::Block & block = cifdoc.add_new_block("Rosetta");

	// "header information", i.e., is from the Title Section of the PDB file.
	if ( options.preserve_header() ) {
		gemmi::cif::Loop & citation = gemmi_add_table(block, "citation", {"title"});
		gemmi_add_row( citation, { sfr->header()->title() } );

		gemmi::cif::Loop & entry = gemmi_add_table(block, "entry", {"id"});
		gemmi_add_row( entry, { sfr->header()->idCode() } );

		gemmi::cif::Loop & entity = gemmi_add_table(block, "entity", {"pdbx_description"});
		for ( auto e: sfr->header()->compounds() ) {
			gemmi_add_row( entity, {e.second} );
		}

		gemmi::cif::Loop & struct_keywords = gemmi_add_table(block, "struct_keywords", {"pdbx_keywords","text"});
		gemmi_add_row( struct_keywords, {sfr->header()->classification(), utility::join(sfr->header()->keywords(), ", ")} );

		gemmi::cif::Loop & database_PDB_rev = gemmi_add_table(block, "database_PDB_rev", {"date_original"});
		gemmi_add_row( database_PDB_rev, {sfr->header()->deposition_date()} );

		utility::vector1< std::string > expt_tech;
		for ( auto e: sfr->header()->experimental_techniques() ) {
			expt_tech.push_back( rcsb::experimental_technique_to_string(e) );
		}
		gemmi::cif::Loop & exptl = gemmi_add_table(block, "exptl", {"method"} );
		gemmi_add_row( exptl, {utility::join(expt_tech,", ")} );
	}

	// HETNAM
	gemmi::cif::Loop & chem_comp = gemmi_add_table(block, "chem_comp", {"id","name"} );
	for ( auto const & elem : sfr->heterogen_names() ) {
		gemmi_add_row( chem_comp, {elem.first, elem.second} );
	}

	// LINK
	gemmi::cif::Loop & struct_conn = gemmi_add_table(block, "struct_conn", {
			"ptnr1_label_atom_id",
			"ptnr1_label_comp_id",
			"ptnr1_label_asym_id",
			"ptnr1_label_seq_id",
			"pdbx_ptnr1_PDB_ins_code",
			"ptnr2_label_atom_id",
			"ptnr2_label_comp_id",
			"ptnr2_label_asym_id",
			"ptnr2_label_seq_id",
			"pdbx_ptnr2_PDB_ins_code",
			"pdbx_dist_value",
			"conn_type_id"
	} );

	for ( auto const & elem : sfr->ssbond_map() ) {
		for ( auto const & iter2 : elem.second ) {
			std::vector< std::string > vec;

			vec.push_back( iter2.name1 );
			vec.push_back( iter2.resName1 );
			vec.emplace_back(1,iter2.chainID1 );
			std::stringstream ss;
			ss << iter2.resSeq1;
			vec.push_back( ss.str() );
			vec.emplace_back(1,iter2.iCode1 );
			vec.push_back( iter2.name2 );
			vec.push_back( iter2.resName2 );
			vec.emplace_back(1,iter2.chainID2 );
			ss.str( std::string() );
			ss << iter2.resSeq2;
			vec.push_back( ss.str() );
			vec.emplace_back(1,iter2.iCode2 );
			ss.str( std::string() );
			ss << iter2.length;
			vec.push_back( ss.str() );
			vec.emplace_back("disulf" );

			gemmi_add_row( struct_conn, vec );
		}
	}

	for ( auto const & elem : sfr->link_map() ) {
		for ( auto const & iter2 : elem.second ) {
			std::vector< std::string > vec;
			vec.push_back( iter2.name1 );
			vec.push_back( iter2.resName1 );
			vec.emplace_back(1,iter2.chainID1 );
			std::stringstream ss;
			ss << iter2.resSeq1;
			vec.push_back( ss.str() );
			vec.emplace_back(1,iter2.iCode1 );
			vec.push_back( iter2.name2 );
			vec.push_back( iter2.resName2 );
			vec.emplace_back(1,iter2.chainID2 );
			ss.str( std::string() );
			ss << iter2.resSeq2;
			vec.push_back( ss.str() );
			vec.emplace_back(1,iter2.iCode2 );
			ss.str( std::string() );
			ss << iter2.length;
			vec.push_back( ss.str() );
			vec.emplace_back("covale" );

			gemmi_add_row( struct_conn, vec );
		}
	}

	// CRYST1
	gemmi::cif::Loop & cell = gemmi_add_table(block, "cell", {
		"length_a",
		"length_b",
		"length_c",
		"angle_alpha",
		"angle_beta",
		"angle_gamma"
	} );
	std::vector< std::string > realvec;
	std::stringstream ss;
	ss << sfr->crystinfo().A();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	ss << sfr->crystinfo().B();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	ss << sfr->crystinfo().C();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	ss << sfr->crystinfo().alpha();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	ss << sfr->crystinfo().beta();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	ss << sfr->crystinfo().gamma();
	realvec.push_back( ss.str() );
	ss.str( std::string() );
	gemmi_add_row( cell, realvec );

	gemmi::cif::Loop & symmetry = gemmi_add_table(block, "symmetry", {"space_group_name_H-M"} );
	gemmi_add_row( symmetry, {sfr->crystinfo().spacegroup()} );

	// We cannot support REMARKs yet because these are stored in DIVERSE places.
	// There isn't a coherent "REMARKs" object. AMW TODO

	gemmi::cif::Loop & atom_site = gemmi_add_table(block, "atom_site", {
			"group_PDB",
			"id",
			"label_atom_id",
			"label_alt_id",
			"label_comp_id",
			"label_asym_id",
			"label_seq_id",
			"pdbx_PDB_ins_code",
			"Cartn_x",
			"Cartn_y",
			"Cartn_z",
			"occupancy",
			"B_iso_or_equiv",
			"type_symbol",
			"pdbx_PDB_model_num",
			"auth_asym_id" // https://mmcif.wwpdb.org says it's manditory
			// "label_entity_id" // https://mmcif.wwpdb.org says it's manditory, but we don't track that info
	} );

	// ATOM/HETATM
	for ( Size i = 0; i < sfr->chains().size(); ++i ) {
		for ( Size j = 0; j < sfr->chains()[i].size(); ++j ) {
			AtomInformation const & ai = sfr->chains()[ i ][j];

			std::vector< std::string > vec;
			vec.emplace_back(ai.isHet ? "HETATM" : "ATOM" ); // group_PDB
			std::stringstream ss2;
			ss2 << ai.serial;
			vec.push_back( ss2.str() ); // id
			vec.push_back( utility::strip(ai.name) ); // label_atom_id
			vec.emplace_back(1, ai.altLoc == ' ' ? '?' : ai.altLoc ); // label_alt_id
			vec.push_back( ai.resName ); // label_comp_id
			vec.emplace_back(1, ai.chainID == ' ' ? '?' : ai.chainID ); // label_asym_id
			ss2.str(std::string());
			ss2 << ai.resSeq;
			vec.push_back( ss2.str() ); // label_seq_id
			vec.emplace_back(1, ai.iCode == ' ' ? '?' : ai.iCode ); // pdbx_PDB_ins_code

			// NOTE: we will never be writing out "     nan" here.
			ss2.str(std::string());
			ss2 << ai.x;
			vec.push_back( ss2.str() ); // Cartn_x

			ss2.str(std::string());
			ss2 << ai.y;
			vec.push_back( ss2.str() ); // Cartn_y

			ss2.str(std::string());
			ss2 << ai.z;
			vec.push_back( ss2.str() ); // Cartn_z

			ss2.str(std::string());
			ss2 << ai.occupancy;
			vec.push_back( ss2.str() ); // occupancy

			ss2.str(std::string());
			ss2 << ai.temperature;
			vec.push_back( ss2.str() ); // B_iso_or_equiv

			vec.push_back( utility::strip(ai.element) ); // type_symbol
			if ( sfr->modeltag().empty() ) {
				vec.push_back( "?" ); //  pdbx_PDB_model_num
			} else {
				vec.push_back( sfr->modeltag() ); //pdbx_PDB_model_num
			}

			vec.emplace_back(1, ai.chainID == ' ' ? '?' : ai.chainID ); // auth_asym_id

			gemmi_add_row( atom_site, vec );
		}
	}

	// Pose Energies Table
	if ( sfr->score_table_labels().size() > 0 && options.output_pose_energies_table() ) {
		std::vector<std::string> columns;
		columns.push_back( "label" );
		for ( core::Size i =1; i <= sfr->score_table_labels().size(); ++i ) {
			columns.push_back( sfr->score_table_labels()[ i ] );
		}

		gemmi::cif::Loop & pose_energies = gemmi_add_table(block, "pose_energies", columns );

		//Add the Score Weights as a Row
		std::vector< std::string > weights;
		weights.emplace_back("weights");
		for ( core::Size i = 1; i <= sfr->score_table_weights().size(); ++i ) {
			weights.push_back( to_string( sfr->score_table_weights()[ i ] ) );
		}
		weights.emplace_back("NA");
		gemmi_add_row( pose_energies, weights );
		//Add the rest of the Rows
		for ( core::Size i = 1; i <= sfr->score_table_lines().size(); ++i ) {
			gemmi_add_row( pose_energies, sfr->score_table_lines()[ i ] );
		}
	}

	// Pose Arbitrary String and Float Data.
	if ( (sfr->pose_cache_string_data().size() > 0 || sfr->pose_cache_real_data().size() > 0) && options.output_pose_cache() ) {
		gemmi::cif::Loop & pose_cache = gemmi_add_table(block, "pose_cache", {"key","value"} );

		if ( sfr->pose_cache_string_data().size() > 0 ) {
			for ( auto & it : sfr->pose_cache_string_data() ) {
				std::vector< std::string > row(2, "");
				row[0] = it.first;
				row[1] = it.second;
				gemmi_add_row( pose_cache, row );
			}
		}

		if ( sfr->pose_cache_real_data().size() > 0 ) {
			for ( auto & it : sfr->pose_cache_real_data() ) {
				std::vector< std::string> row(2, "");
				row[0] = it.first;
				row[1] = utility::to_string( it.second ); //PDB Writing of this data had no rounding of decimal places, so I'm not doing it here either (JAB).
				gemmi_add_row( pose_cache, row );
			}
		}
	}

	//Even more pose info.
	if ( sfr->pdb_comments().size() > 0 ) {
		gemmi::cif::Loop & pose_comments = gemmi_add_table(block, "pose_comments", {"type","comment"} );

		using namespace std;
		map< string, string > const comments = sfr->pdb_comments();
		for ( auto const & comment : comments ) {
			std::vector< std::string > row(2, "");
			row[0] = comment.first;
			row[1] = comment.second;
			gemmi_add_row( pose_comments, row );
		}
	}

	//Remarks are now split into specific types in the mmCIF format.  Here we just place remarks as Rosetta remarks.
	if ( sfr->remarks()->size() > 0 ) {
		gemmi::cif::Loop & rosetta_remarks = gemmi_add_table(block, "rosetta_remarks", {"num","remark"} );

		for ( Size i=0; i<sfr->remarks()->size(); ++i ) {
			RemarkInfo const & ri( sfr->remarks()->at(i) );
			std::vector< std::string > row(2, "");
			row[0] = utility::pad_left( ri.num, 3 ); //("%3d", ri.num);
			row[1] = ri.value;
			gemmi_add_row(rosetta_remarks, row );
		}
	}

	//Finally, the single string output.  If you have a better way to do this, please refactor this.
	if ( !sfr->additional_string_output().empty() || ! sfr->foldtree_string().empty() ) {
		gemmi::cif::Loop & rosetta_additional = gemmi_add_table(block, "rosetta_additional", {"type","output"} );

		if ( ! sfr->foldtree_string().empty() ) {
			std::vector< std::string > out_vec;

			out_vec.emplace_back("fold_tree");
			out_vec.push_back(sfr->foldtree_string());
			gemmi_add_row( rosetta_additional, out_vec);
		}
		if ( ! sfr->additional_string_output().empty() ) {
			std::vector< std::string > out_vec;

			out_vec.emplace_back("etc" );
			out_vec.push_back( sfr->additional_string_output() );
			gemmi_add_row( rosetta_additional, out_vec );
		}
	}

	gemmi::cif::WriteOptions gemmi_options;
	gemmi_options.misuse_hash = true; // Use hash separation, like the wwPDB does.
	gemmi_options.prefer_pairs = true; // Write single row tables as pairs, rather than loops
	gemmi_options.align_loops = 16; // Space pad things to align columns (up to 16 columns wide)

	gemmi::cif::write_cif_to_stream(out,cifdoc,gemmi_options);
}

} //core
} //io
} //mmcif
