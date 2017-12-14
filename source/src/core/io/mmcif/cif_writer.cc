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


#include <core/io/mmcif/cif_writer.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Package headers
#include <core/io/pdb/pdb_reader.hh>  // TODO: Pull out pseudo-duplicated code and move to sfr_storage.cc.

// When you move PDBReader and PoseUnbuilder, take these.
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueConnection.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>

#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

#include <core/io/Remarks.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <cifparse/CifFile.h>
using CifFileOP = utility::pointer::shared_ptr<CifFile>;

// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>


#include <basic/Tracer.hh>

static basic::Tracer TR( "core.io.mmcif.cif_writer.hh" );

namespace core {
namespace io {
namespace mmcif {

using basic::Error;
using basic::Warning;
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
	utility::io::ozstream & out)

{
	StructFileRepOptionsOP options =  StructFileRepOptionsOP( new StructFileRepOptions );
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
	// ALERT: do not delete pointers. block.WriteTable takes ownership.
	CifFile cifFile;
	cifFile.AddBlock( "Rosetta" );
	Block& block = cifFile.GetBlock( "Rosetta" );

	// "header information", i.e., is from the Title Section of the PDB file.
	if ( options.preserve_header() ) {
		ISTable* citation = new ISTable( "citation" );
		citation->AddColumn( "title" );
		citation->AddRow( std::vector< std::string >( 1, sfr->header()->title() ) );
		block.WriteTable( citation );

		ISTable* entry = new ISTable( "entry" );
		entry->AddColumn( "id" );
		entry->AddRow( std::vector< std::string >( 1, sfr->header()->idCode() ) );
		block.WriteTable( entry );

		ISTable* entity = new ISTable( "entity" );
		entity->AddColumn( "pdbx_description" );
		for ( Size i = 0; i < sfr->header()->compounds().size(); ++i ) {
			entity->AddRow( std::vector< std::string >( 1, sfr->header()->compounds()[ i ].second ) );
		}
		block.WriteTable( entity );

		ISTable* struct_keywords = new ISTable( "struct_keywords" );
		struct_keywords->AddColumn( "pdbx_keywords" );
		struct_keywords->AddColumn( "text" );
		std::vector< std::string > struct_keywords_vec;
		struct_keywords_vec.push_back( sfr->header()->classification() );
		std::string keyword_str = "";
		// keywords is a std::set, so <
		std::list< std::string >::const_iterator iter, end;
		for ( iter = sfr->header()->keywords().begin(),
				end = --sfr->header()->keywords().end(); iter != end; ++iter ) {
			keyword_str += *iter + ", ";
		}
		keyword_str += *++iter;

		struct_keywords_vec.push_back( keyword_str );
		struct_keywords->AddRow( struct_keywords_vec );
		block.WriteTable( struct_keywords );

		ISTable* database_PDB_rev = new ISTable( "database_PDB_rev" );
		database_PDB_rev->AddColumn( "date_original" );
		database_PDB_rev->AddRow( std::vector< std::string >( 1, sfr->header()->deposition_date() ) );
		block.WriteTable( database_PDB_rev );

		ISTable* exptl = new ISTable( "exptl" );
		exptl->AddColumn( "method" );
		std::string tech_str = "";
		// keywords is a std::list, so <
		std::list< HeaderInformation::ExperimentalTechnique >::const_iterator iter3, end3;
		for ( iter3 = sfr->header()->experimental_techniques().begin(),
				end3 = --sfr->header()->experimental_techniques().end(); iter3 != end3; ++iter3 ) {
			tech_str += HeaderInformation::experimental_technique_to_string(
				*iter3 ) + ", ";
		}
		tech_str += HeaderInformation::experimental_technique_to_string( *++iter3 );
		exptl->AddRow( std::vector< std::string >( 1, tech_str ) );
		block.WriteTable( exptl );
	}

	// HETNAM
	ISTable* chem_comp = new ISTable( "chem_comp" );
	chem_comp->AddColumn( "id" );
	chem_comp->AddColumn( "name" );
	for ( std::map< std::string, std::string >::const_iterator iter = sfr->heterogen_names().begin(),
			end = sfr->heterogen_names().end(); iter != end; ++iter ) {
		std::vector< std::string > vec;
		vec.push_back( iter->first );
		vec.push_back( iter->second );
		chem_comp->AddRow( vec );
	}
	block.WriteTable( chem_comp );

	// LINK
	ISTable* struct_conn = new ISTable( "struct_conn" );
	struct_conn->AddColumn( "ptnr1_label_atom_id" );
	struct_conn->AddColumn( "ptnr1_label_comp_id" );
	struct_conn->AddColumn( "ptnr1_label_asym_id" );
	struct_conn->AddColumn( "ptnr1_label_seq_id" );
	struct_conn->AddColumn( "pdbx_ptnr1_PDB_ins_code" );
	struct_conn->AddColumn( "ptnr2_label_atom_id" );
	struct_conn->AddColumn( "ptnr2_label_comp_id" );
	struct_conn->AddColumn( "ptnr2_label_asym_id" );
	struct_conn->AddColumn( "ptnr2_label_seq_id" );
	struct_conn->AddColumn( "pdbx_ptnr2_PDB_ins_code" );
	struct_conn->AddColumn( "pdbx_dist_value" );
	struct_conn->AddColumn( "conn_type_id" );

	for ( std::map< std::string, utility::vector1< SSBondInformation > >::const_iterator iter = sfr->ssbond_map().begin(), end = sfr->ssbond_map().end(); iter != end; ++iter ) {
		for ( auto const & iter2 : iter->second ) {
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

			struct_conn->AddRow( vec );
		}
	}

	for ( std::map< std::string, utility::vector1< LinkInformation > >::const_iterator iter = sfr->link_map().begin(), end = sfr->link_map().end(); iter != end; ++iter ) {
		for ( auto const & iter2 : iter->second ) {
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

			struct_conn->AddRow( vec );
		}
	}
	block.WriteTable( struct_conn );

	// CRYST1
	ISTable* cell = new ISTable( "cell" );
	cell->AddColumn( "length_a" );
	cell->AddColumn( "length_b" );
	cell->AddColumn( "length_c" );
	cell->AddColumn( "angle_alpha" );
	cell->AddColumn( "angle_beta" );
	cell->AddColumn( "angle_gamma" );
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
	cell->AddRow( realvec );
	block.WriteTable( cell );

	ISTable* symmetry = new ISTable( "symmetry" );
	symmetry->AddColumn( "space_group_name_H-M" );
	symmetry->AddRow( std::vector< std::string >( 1, sfr->crystinfo().spacegroup() ) );
	block.WriteTable( symmetry );

	// We cannot support REMARKs yet because these are stored in DIVERSE places.
	// There isn't a coherent "REMARKs" object. AMW TODO

	ISTable* atom_site = new ISTable( "atom_site" );
	atom_site->AddColumn( "group_PDB" );
	atom_site->AddColumn( "id" );
	atom_site->AddColumn( "auth_atom_id" );
	atom_site->AddColumn( "label_alt_id" );
	atom_site->AddColumn( "auth_comp_id" );
	atom_site->AddColumn( "auth_asym_id" );
	atom_site->AddColumn( "auth_seq_id" );
	atom_site->AddColumn( "pdbx_PDB_ins_code" );
	atom_site->AddColumn( "Cartn_x" );
	atom_site->AddColumn( "Cartn_y" );
	atom_site->AddColumn( "Cartn_z" );
	atom_site->AddColumn( "occupancy" );
	atom_site->AddColumn( "B_iso_or_equiv_esd" );
	atom_site->AddColumn( "type_symbol" );
	atom_site->AddColumn( "pdbx_PDB_model_num" );

	// ATOM/HETATM
	for ( Size i = 0; i < sfr->chains().size(); ++i ) {
		for ( Size j = 0; j < sfr->chains()[i].size(); ++j ) {
			AtomInformation const & ai = sfr->chains()[ i ][j];

			std::vector< std::string > vec;
			vec.emplace_back(ai.isHet ? "HETATM" : "ATOM" );
			std::stringstream ss;
			ss << ai.serial;
			vec.push_back( ss.str() );
			vec.push_back( utility::strip(ai.name) );
			vec.emplace_back(1, ai.altLoc == ' ' ? '?' : ai.altLoc );
			vec.push_back( ai.resName );
			vec.emplace_back(1, ai.chainID == ' ' ? '?' : ai.chainID );
			ss.str(std::string());
			ss << ai.resSeq;
			vec.push_back( ss.str() );
			vec.emplace_back(1, ai.iCode == ' ' ? '?' : ai.iCode );

			// NOTE: we will never be writing out "     nan" here.
			ss.str(std::string());
			ss << ai.x;
			vec.push_back( ss.str() );

			ss.str(std::string());
			ss << ai.y;
			vec.push_back( ss.str() );

			ss.str(std::string());
			ss << ai.z;
			vec.push_back( ss.str() );

			ss.str(std::string());
			ss << ai.occupancy;
			vec.push_back( ss.str() );

			ss.str(std::string());
			ss << ai.temperature;
			vec.push_back( ss.str() );

			vec.push_back( utility::strip(ai.element) );
			vec.push_back( sfr->modeltag() );

			atom_site->AddRow( vec );
		}
	}
	block.WriteTable( atom_site );

	// Pose Energies Table
	if ( sfr->score_table_labels().size() > 0 && options.output_pose_energies_table() ) {
		ISTable* pose_energies = new ISTable( "pose_energies" );

		//Add the Columns
		pose_energies->AddColumn( "label" );
		for ( core::Size i =1; i <= sfr->score_table_labels().size(); ++i ) {
			pose_energies->AddColumn( sfr->score_table_labels()[ i ] );
		}

		//Add the Score Weights as a Row
		std::vector< std::string > weights;
		weights.emplace_back("weights");
		for ( core::Size i = 1; i <= sfr->score_table_weights().size(); ++i ) {
			weights.push_back( to_string( sfr->score_table_weights()[ i ] ) );
		}
		weights.emplace_back("NA");
		pose_energies->AddRow(weights);
		//Add the rest of the Rows
		for ( core::Size i = 1; i <= sfr->score_table_lines().size(); ++i ) {
			pose_energies->AddRow( sfr->score_table_lines()[ i ] );
		}
		block.WriteTable( pose_energies );
	}

	// Pose Arbitrary String and Float Data.
	if ( (sfr->pose_cache_string_data().size() > 0 || sfr->pose_cache_float_data().size() > 0) && options.output_pose_cache() ) {
		ISTable* pose_cache = new ISTable( "pose_cache_data" );

		if ( sfr->pose_cache_string_data().size() > 0 ) {
			for ( auto & it : sfr->pose_cache_string_data() ) {
				std::vector< std::string > row(2, "");
				row[0] = it.first;
				row[1] = it.second;
				pose_cache->AddRow( row );
			}
		}

		if ( sfr->pose_cache_float_data().size() > 0 ) {
			for ( auto & it : sfr->pose_cache_float_data() ) {
				std::vector< std::string> row(2, "");
				row[0] = it.first;
				row[1] = utility::to_string( it.second ); //PDB Writing of this data had no rounding of decimal places, so I'm not doing it here either (JAB).
				pose_cache->AddRow( row );
			}
		}
		block.WriteTable( pose_cache );

	}

	//Even more pose info.
	if ( sfr->pdb_comments().size() > 0 ) {
		ISTable* pose_comments = new ISTable( "pose_comments" );
		pose_comments->AddColumn( "type" );
		pose_comments->AddColumn( "comment" );


		using namespace std;
		map< string, string > const comments = sfr->pdb_comments();
		for ( auto const & comment : comments ) {
			std::vector< std::string > row(2, "");
			row[0] = comment.first;
			row[1] = comment.second;
			pose_comments->AddRow( row );
		}

		block.WriteTable( pose_comments );
	}

	//Remarks are now split into specific types in the mmCIF format.  Here we just place remarks as Rosetta remarks.
	if ( sfr->remarks()->size() > 0 ) {
		ISTable* rosetta_remarks = new ISTable( "rosetta_remarks" );
		rosetta_remarks->AddColumn( "num"); //Feel free to change the name of this.
		rosetta_remarks->AddColumn( "remark");

		for ( Size i=0; i<sfr->remarks()->size(); ++i ) {
			RemarkInfo const & ri( sfr->remarks()->at(i) );
			std::vector< std::string > row(2, "");
			row[0] = utility::pad_left( ri.num, 3 ); //("%3d", ri.num);
			row[1] = ri.value;
			rosetta_remarks->AddRow( row );
		}

		block.WriteTable( rosetta_remarks );
	}

	//Finally, the single string output.  If you have a better way to do this, please refactor this.
	if ( !sfr->additional_string_output().empty() || ! sfr->foldtree_string().empty() ) {
		ISTable* rosetta_additional = new ISTable( "rosetta_additional" );
		rosetta_additional->AddColumn( "type" );
		rosetta_additional->AddColumn( "output");

		if ( ! sfr->foldtree_string().empty() ) {
			std::vector< std::string > out;

			out.emplace_back("fold_tree");
			out.push_back(sfr->foldtree_string());
			rosetta_additional->AddRow(out);
		}
		if ( ! sfr->additional_string_output().empty() ) {
			std::vector< std::string > out;

			out.emplace_back("etc" );
			out.push_back( sfr->additional_string_output() );
			rosetta_additional->AddRow( std::vector< std::string >( 1, sfr->additional_string_output() ) );
		}
		block.WriteTable( rosetta_additional );
	}

	cifFile.Write( out );
}

} //core
} //io
} //mmcif
