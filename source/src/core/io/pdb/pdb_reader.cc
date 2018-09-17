// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/io/pdb/pdb_reader.cc
/// @brief   Function definitions for reading of .pdb files.
/// @author  Sergey Lyskov (Sergey.Lyskov@jhu.edu)
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    We may perhaps decide in the future to wrap this functionality in a class.  ~Labonte


// Unit headers
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/RecordCollection.hh>

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileReaderOptions.hh>  // TODO: Rename after refactoring is complete.
#include <core/io/HeaderInformation.hh>
#include <core/io/Remarks.hh>
#include <core/io/NomenclatureManager.hh>

// Utility headers
#include <utility/string_constants.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/Tracer.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ header
#include <algorithm>
#include <sstream>
#include <map>


static basic::Tracer TR( "core.io.pdb.pdb_reader" );


namespace core {
namespace io {
namespace pdb {

// Convert a .pdb file line into a Record data structure.
Record
create_record_from_pdb_line( std::string const & line )
{
	std::string resized_line( line ) ;
	if ( resized_line.size() < 80 ) {  // allow for longer Rosetta-specific lines
		resized_line.resize( 80, ' ' );  // standard .pdb line width
	}

	// All PDB record names are the first 6 characters of a line.
	std::string const & record_type( resized_line.substr( 0, 6 ) );

	Record record( RecordCollection::record_from_record_type( record_type ) );
	for ( auto & field : record ) {
		field.second.set_value_from_pdb_line( resized_line );
	}

	return record;
}


// Create a list of .pdb format records from the lines from a .pdb file.
utility::vector1< Record >
create_records_from_pdb_lines( utility::vector1< std::string > const & lines )
{
	runtime_assert( ! lines.empty() );  // We're wasting time if there's no data here....

	utility::vector1< Record > records( lines.size() );
	std::transform( lines.begin(), lines.end(), records.begin(), create_record_from_pdb_line );
	return records;
}

// Create a list of .pdb format records from the entire contents of a .pdb file.
utility::vector1< Record >
create_records_from_pdb_file_contents( std::string const & pdb_contents )
{
	utility::vector1< std::string > lines( utility::split_by_newlines( pdb_contents ) );
	return create_records_from_pdb_lines( lines );
}


/// @remarks  The bulk of SFR generation occurs within this function.
// Create a representation of structural file data from a list of .pdb format records with options.
StructFileRep
create_sfr_from_pdb_records( utility::vector1< Record > & records, StructFileReaderOptions const & options )
{
	using namespace std;
	StructFileRep sfr;

	map< char, ChainAtoms > chain_atoms_map;  // A map of chain ID to every atom in the chain.
	Size ter_record_count = 0;
	utility::vector1< char > chain_list;
	std::map< char, Size > chain_to_idx;

	// Prepare for multi-model .pdbs.
	map< std::pair< Size, Size >, char > modelchain_to_chain;  // TODO: What is this?
	string const & chain_letters( utility::UPPERCASE_ALPHANUMERICS );
	for ( Size i = 0; i < chain_letters.size(); ++i ) {
		modelchain_to_chain[ std::pair< Size, Size >( 0, i ) ] = chain_letters[ i ];
		modelchain_to_chain[ std::pair< Size, Size >( 1, i ) ] = chain_letters[ i ];
	}
	Size modelidx = 1;
	bool modeltags_present = false;


	// Loop over all PDB records.
	Size const n_records( records.size() );
	utility::vector1< Record > unknow_records;  // place to store unrecognized records for 2nd pass for non-RCSB
	bool stop_reading_coordinate_section( false );  // might be set to true when ENDMDL is reached
	for ( core::uint i = 1; i <= n_records; ++i ) {
		std::string const record_type( records[ i ][ "type" ].value );
		TR.Trace << "Record type: " << record_type << endl;

		// Title Section //////////////////////////////////////////////////////
		// Record contains "header information", i.e., is from the Title Section of the PDB file.
		if ( (record_type == "HEADER" || record_type == "KEYWDS" ||
				record_type == "TITLE " || record_type == "COMPND" ||
				record_type == "EXPDTA" ) && (! options.read_only_ATOM_entries() ) ) {
			// TODO: Add rest of Title Section records.
			sfr.header()->store_record( records[ i ] );

			// Record contains a remark from the Title Section of the PDB file.
		} else if ( record_type == "REMARK" && ! options.read_only_ATOM_entries() )  {
			RemarkInfo ri;
			ri.num = atoi( records[ i ][ "remarkNum" ].value.c_str() );
			ri.value = records[ i ][ "value" ].value;

			//Added by DANIEL to skip reading the PDBinfo-LABEL, which comes "nfo-LABEL:"
			//Those are read in a different way, using: core.import_pose().read_additional_pdb_data()
			if ( ( ri.value.size() >= 10 ) && ( ri.value.substr( 0, 10 ) == "nfo-LABEL:" ) ) {
				continue;
			}
			sfr.remarks()->push_back( ri );


			// Primary Structure Section //////////////////////////////////////////
			// Record contains cross-references from PDB sequence fragments to a corresponding database sequence.
		} else if ( ( record_type == "DBREF " || record_type == "DBREF1" || record_type == "DBREF2" ||
				record_type == "SEQADV" ) && ( ! options.read_only_ATOM_entries() ) ) {
			//sfr.primary_struct_info()->store_sequence_database_refs( records[ i ] );  // TODO
			continue;  // TEMP

			// Record contains a linear (or cyclic) primary sequence declaration.
		} else if ( record_type == "SEQRES" && ! options.read_only_ATOM_entries() ) {
			store_chain_sequence_record_in_sfr( records[ i ], sfr );

			// Record that a residue is modified and how.
		} else if ( record_type == "MODRES" && ! options.read_only_ATOM_entries() ) {
			store_mod_res_record_in_sfr( records[ i ], sfr );


			// Heterogen Section //////////////////////////////////////////////////
			// Record contains heterogen nomenclature information.
		} else if ( record_type == "HETNAM" && ! options.read_only_ATOM_entries() ) {
			store_heterogen_name_record_in_sfr( records[ i ], sfr );

			// Record contains heterogen synonym information.
		} else if ( record_type == "HETSYN" && ! options.read_only_ATOM_entries() ) {
			store_heterogen_synonym_record_in_sfr( records[ i ], sfr );

			// Record contains formula information.
		} else if ( record_type == "FORMUL" && ! options.read_only_ATOM_entries() ) {
			store_formula_record_in_sfr( records[ i ], sfr );


			// Secondary Structure Section ////////////////////////////////////////
			// Record contains helix definitions.
		} else if ( record_type == "HELIX " && ! options.read_only_ATOM_entries() ) {
			// TODO: Store HELIX record types here.
			continue;

			// Record contains sheet definitions.
		} else if ( record_type == "SHEET " && ! options.read_only_ATOM_entries() ) {
			// TODO: Store SHEET record types here.
			continue;


			// Connectivity Annotation Section ////////////////////////////////////
			// Record contains disulfide linkage information.
		} else if ( record_type == "SSBOND" && ! options.read_only_ATOM_entries() ) {
			store_ssbond_record_in_sfr( records[ i ], sfr );

			// Record contains nonstandard polymer linkage information.
		} else if ( record_type == "LINK  " && ! options.read_only_ATOM_entries() ) {
			store_link_record_in_sfr( records[ i ], sfr );

		} else if ( record_type == "CISPEP" && ! options.read_only_ATOM_entries() ) {
			store_cis_peptide_record_in_sfr( records[ i ], sfr );


			// Miscellaneous Features Section /////////////////////////////////////
		} else if ( record_type == "SITE  " && ! options.read_only_ATOM_entries() ) {
			// TODO: Store SITE record types here.
			continue;

			// Crystallographic and Coordinate Transformation Section /////////////
			// Record contains crystal information.
		} else if ( records[i]["type"].value == "CRYST1" && ! options.read_only_ATOM_entries() )  {
			store_crystallographic_parameter_record_in_sfr( records[ i ], sfr );


			// Coordinate Section /////////////////////////////////////////////////
			// Record contains multimodel PDBs.
		} else if ( record_type == "MODEL " ) {
			if ( stop_reading_coordinate_section ) continue;

			// store the serial number as the filename, which will become the PDBInfo name of the pose
			std::string temp_model = ObjexxFCL::strip_whitespace( records[i]["serial"].value ) ;
			sfr.modeltag() = temp_model.c_str();
			if ( options.new_chain_order() ) {
				if ( modeltags_present ) {
					// second model... all chains should be present...
					for ( Size model_idx=2; model_idx*chain_to_idx.size()<chain_letters.size(); ++model_idx ) {
						for ( Size chain_idx=1; chain_idx <= chain_to_idx.size(); ++chain_idx ) {
							TR << "REARRANGE CHAINS " << model_idx << " " << chain_idx << " ";
							TR << (model_idx-1)*chain_to_idx.size()+chain_idx << std::endl;
							modelchain_to_chain[std::pair<Size, Size>(model_idx, chain_idx)] =
								chain_letters[(model_idx-1)*chain_to_idx.size() + chain_idx - 1];
						}
					}
					++modelidx;
					if ( modelidx > 8 ) utility_exit_with_message("quitting: too many MODELs");
				} else {
					modeltags_present = true;
				}
			}

			// Record contains atom information.
		} else if ( record_type == "ATOM  " || record_type == "HETATM" ) {
			if ( stop_reading_coordinate_section ) continue;
			if ( record_type == "HETATM" && options.read_only_ATOM_entries() ) continue;

			// TODO: Refactor?
			Record & R( records[ i ] );

			AtomInformation ai;
			ai.isHet = ( R[ "type" ].value == "HETATM" );  // probably not worth storing here
			ai.serial = atoi( R[ "serial" ].value.c_str() );
			ai.name = R[ "name" ].value;
			ai.altLoc = 0;
			if ( R[ "altLoc" ].value.size() > 0 ) {
				ai.altLoc = R[ "altLoc" ].value[ 0 ];
			}

			ai.resName = R[ "resName" ].value;
			ai.chainID = 0;
			if ( R[ "chainID" ].value.size() > 0 ) {
				ai.chainID = R[ "chainID" ].value[ 0 ];
			}
			if ( options.new_chain_order() ) {
				if ( R["chainID"].value.size() > 0 ) {
					char chainid = R["chainID"].value[0];
					if ( chain_to_idx.find(chainid) == chain_to_idx.end() ) {
						chain_to_idx[chainid] = chain_to_idx.size()+1;
						TR << "found new chain " << chainid << " " << chain_to_idx.size() << std::endl;
					}
					ai.chainID = modelchain_to_chain[std::pair<Size, Size>(modelidx, chain_to_idx[chainid])];
				}
			}

			ai.resSeq = atoi( R["resSeq"].value.c_str() );
			ai.iCode = 0;
			if ( R["iCode"].value.size() > 0 ) ai.iCode = R["iCode"].value[0];

			// how can you check properly if something will successfully convert to a number !?!?!?
			bool force_no_occupancy = false;
			if ( R["x"].value == "     nan" ) {
				ai.x =0.0;
				force_no_occupancy=true;
			} else {
				ai.x = atof( R["x"].value.c_str() );
			}
			if ( R["y"].value == "     nan" ) {
				ai.y =0.0;
				force_no_occupancy=true;
			} else {
				ai.y = atof( R["y"].value.c_str() );
			}
			if ( R["z"].value == "     nan" ) {
				ai.z =0.0;
				force_no_occupancy=true;
			} else {
				ai.z = atof( R["z"].value.c_str() );
			}

			// check that the occupancy column actually exists. If it doesn't, assume full occupancy.
			// otherwise read it.
			if ( R["occupancy"].value == "      " ) {
				ai.occupancy = 1.0;
			} else {
				ai.occupancy = atof( R["occupancy"].value.c_str() );
			}
			if ( force_no_occupancy ) {
				ai.occupancy = -1.0;
			}

			ai.temperature = atof( R["tempFactor"].value.c_str() );
			// if ( ai.resName == "CYS" ) {
			//  TR << "In reader " << R["tempFactor"] << endl;
			// }
			ai.segmentID = R["segmentID"].value;
			ai.element = R["element"].value;
			ai.terCount = ter_record_count;

			chain_atoms_map[ ai.chainID ].push_back( ai );
			if ( find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}

		} else if ( ! stop_reading_coordinate_section && ( record_type == "TER   " || record_type == "END   " ) ) {
			++ter_record_count;

		} else if ( record_type == "ENDMDL" )  {
			if ( options.obey_ENDMDL() )  {
				TR.Warning << "Hit ENDMDL; not reading further coordinate section records." << endl;
				stop_reading_coordinate_section = true;
			}
		} else /*UNKNOW record*/ {
			unknow_records.push_back( records[ i ] );
		}
	}

	sfr.header()->finalize_parse();

	for ( Size i = 1; i <= chain_list.size(); ++i ) {
		sfr.chains().push_back( chain_atoms_map.find( chain_list[ i ] )->second );
	}

	// Now check through the list of unknown record lines to see if there are any Rosetta-specific things that we care
	// about loading.
	if ( options.pdb_comments() ) {
		store_unknown_records_in_sfr( unknow_records, sfr );
	}

	return sfr;
}

// Create a representation of structural file data from a list of .pdb format records.
StructFileRep
create_sfr_from_pdb_records( utility::vector1< Record > & records )
{
	StructFileReaderOptions options;
	return create_sfr_from_pdb_records( records, options );
}


// Create a representation of structural file data from .pdb file contents with options.
StructFileRep
create_sfr_from_pdb_file_contents( std::string const & pdb_contents, StructFileReaderOptions const & options )
{
	utility::vector1< Record > records( create_records_from_pdb_file_contents( pdb_contents ) );
	return create_sfr_from_pdb_records( records, options );
}

// Create a representation of structural file data from .pdb file contents.
StructFileRep
create_sfr_from_pdb_file_contents( std::string const & pdb_contents )
{
	StructFileReaderOptions options;
	return create_sfr_from_pdb_file_contents( pdb_contents, options );
}


// The 2 functions below are covered by import_pose.cc.
// Create a representation of structural file data from a .pdb file by file.
//StructFileRep
//create_sfr_from_pdb_file( utility::io::izstream const & file )
//{
// TODO: Slurp file
//return ;
//}

// Create a representation of structural file data from a .pdb file by filename.
//StructFileRep
//create_sfr_from_pdb_file( std::string const & filename )
//{
// return create_sfr_from_pdb_file( utility::io::izstream( filename ) );
//}


// .pdb Record Storage Functions //////////////////////////////////////////////
// Convert .pdb SEQRES record into SFR data.
void
store_chain_sequence_record_in_sfr( Record seqres_record, StructFileRep & sfr )
{
	sfr.chain_sequences()[ seqres_record[ "chainID" ].value[ 0 ] ].push_back( seqres_record[ "resName1" ].value );
	for ( uint i( 2 ); i <= 13; ++i ) {  // There are 13 total resName fields per SEQRES record.
		std::ostringstream field_name;
		field_name << "resName" << i;
		if ( seqres_record[ field_name.str() ].value == "   " ) {  // an empty field
			break;
		} else {
			sfr.chain_sequences()[ seqres_record[ "chainID" ].value[ 0 ] ].push_back(
				seqres_record[ field_name.str() ].value );
		}
	}
	TR.Debug << "SEQRES record information stored successfully." << std::endl;
}

// Convert .pdb MODRES record into SFR data.
void
store_mod_res_record_in_sfr( Record modres_record, StructFileRep & sfr )
{
	using namespace std;

	ModifiedResidueInformation modres_info;
	modres_info.resName = modres_record[ "resName" ].value;
	modres_info.chainID = modres_record[ "chainID" ].value[ 0 ];
	modres_info.seqNum = atof( modres_record[ "seqNum" ].value.c_str() );
	modres_info.iCode = modres_record[ "iCode" ].value[ 0 ];
	modres_info.stdRes = modres_record[ "stdRes" ].value;
	string const comment_w_spaces( modres_record[ "comment" ].value );
	modres_info.comment = utility::trim( comment_w_spaces );

	string const resID(
		modres_record[ "seqNum" ].value + modres_record[ "iCode" ].value + modres_record[ "chainID" ].value );
	sfr.modres_map()[ resID ] = modres_info;

	TR.Debug << "MODRES record information stored successfully." << endl;
}


// Parse .pdb HETNAM text field to extract full resID and convert into SFR data.
/// @remarks Called by store_heterogen_name_record_in_sfr().
void
store_base_residue_type_name_in_sfr( std::string const & hetID, std::string const & text_field, StructFileRep & sfr )
{
	using namespace std;

	if ( text_field.size() < 8 ) { return; }  // must contain the 6-char. resID, a space, and at least 1 more char.

	string const chainID( string( text_field, 0, 1 ) );  // 1 character for chainID
	string const resSeq( string( text_field, 1, 4 ) );  // 4 characters for resSeq
	string const iCode( string( text_field, 5, 1 ) );  // 1 character for iCode
	string const key( resSeq + iCode + chainID );  // a resID, as defined elsewhere in StructFileRep

	// The name starts after 7th character; any word after a comma is ignored.
	uint const comma_location( text_field.find( ", " ) );  // if not found, will be set as string::npos
	string const needed_residue_type_base_name( string( text_field, 7, comma_location - 7 ) );  // npos-7 >> 80!

	sfr.residue_type_base_names()[ key ] = make_pair( hetID, needed_residue_type_base_name );
}

// Convert .pdb HETNAM record into SFR data.
/// @remarks  Heterogen "names" for carbohydrates (from "Rosetta-ready" PDB
/// files) additionally have the name field parsed to extract the base (non-
/// variant) ResidueType needed for a particular residue.
void
store_heterogen_name_record_in_sfr( Record hetnam_record, StructFileRep & sfr )
{
	using namespace std;

	string const & hetID( hetnam_record[ "hetID" ].value );
	string text( hetnam_record[ "text" ].value );
	utility::trim( text );

	if ( hetID.empty() ) {
		TR.Warning << "PDB HETNAM record is missing an heterogen ID field." << endl;
		return;
	}
	if ( text.empty() ) {
		TR.Warning << "PDB HETNAM chemical name field is an empty string." << endl;
		return;
	}

	if ( sfr.heterogen_names().count( hetID ) == 0 ) {
		sfr.heterogen_names()[ hetID ] = text;
	} else {
		// Check whether we are continuing a hyphenated word or adding a new one.
		if ( sfr.heterogen_names()[ hetID ][ sfr.heterogen_names()[ hetID ].size() - 1 ] == '-' ) {
			// The last character is a hyphen.
			// (If we were using C++11, I would just call back()....)
			sfr.heterogen_names()[ hetID ] += text;
		} else {
			sfr.heterogen_names()[ hetID ] += " " + text;
		}
	}

	store_base_residue_type_name_in_sfr( hetID, text, sfr );
}

// Convert .pdb HETSYN record into SFR data.
void
store_heterogen_synonym_record_in_sfr( Record hetsyn_record, StructFileRep & sfr )
{
	using namespace std;
	using namespace utility;

	string const synonyms_field( hetsyn_record[ "hetSynonyms" ].value );
	vector1< string > synonyms( string_split_simple( synonyms_field, ';' ) );
	Size const n_synonyms( synonyms.size() );
	for ( uint i( 1 ); i <= n_synonyms; ++i ) {
		trim( synonyms[ i ] );
	}

	sfr.heterogen_synonyms()[ hetsyn_record[ "hetID" ].value ] = synonyms;

	TR.Debug << "HETSYN record information stored successfully." << endl;
}

// Convert .pdb FORMUL record into SFR data.
void
store_formula_record_in_sfr( Record formul_record, StructFileRep & sfr )
{
	using namespace std;

	string formula( formul_record[ "asterisk" ].value + formul_record[ "text" ].value  );
	formula.erase( formula.find_last_not_of( " " ) + 1 );  // R-trim spaces.

	sfr.heterogen_formulae()[ formul_record[ "hetID" ].value ] = formula;

	TR.Debug << "FORMUL record information stored successfully." << endl;
}


// Convert .pdb SSBOND record into SFR data.
/// @author Watkins
void
store_ssbond_record_in_sfr( Record ssbond_record, StructFileRep & sfr )
{
	using namespace std;
	using namespace utility;
	using ObjexxFCL::stripped;

	SSBondInformation ssbond;
	vector1< SSBondInformation > ssbonds;

	// If an SSBOND defines a connection across crystal cells, there is no way that we can handle this.
	// Ignore such lines.
	// (1555 is equivalent to a translation matrix of [0,0,0].)
	if ( ( ( ssbond_record[ "sym1" ].value != "  1555" ) && ( ssbond_record[ "sym1" ].value != "      " ) ) ||
			( ( ssbond_record[ "sym2" ].value != "  1555" ) && ( ssbond_record[ "sym2" ].value != "      " ) ) ) {
		TR.Debug << "Throwing out SSBOND across crystal cells." << endl;
		return;
	}

	// Extract values from SSBOND ssbond_record fields.
	ssbond.resName1 = ssbond_record[ "resName1" ].value;
	ssbond.chainID1 = ssbond_record[ "chainID1" ].value[ 0 ];
	ssbond.resSeq1 = atof( ssbond_record[ "resSeq1" ].value.c_str() );
	ssbond.iCode1 = ssbond_record[ "iCode1" ].value[ 0 ];

	//ssbond.resID1 = stripped( ssbond_record[ "resSeq1" ].value + ssbond_record[ "iCode1" ].value + ssbond_record[ "chainID1" ].value );
	ssbond.resID1 = ssbond_record[ "resSeq1" ].value + ssbond_record[ "iCode1" ].value + ssbond_record[ "chainID1" ].value;

	ssbond.resName2 = ssbond_record[ "resName2" ].value;
	ssbond.chainID2 = ssbond_record[ "chainID2" ].value[ 0 ];
	ssbond.resSeq2 = atof( ssbond_record[ "resSeq2" ].value.c_str() );
	ssbond.iCode2 = ssbond_record[ "iCode2" ].value[ 0 ];

	//ssbond.resID2 = stripped( ssbond_record[ "resSeq2" ].value + ssbond_record[ "iCode2" ].value + ssbond_record[ "chainID2" ].value );
	ssbond.resID2 = ssbond_record[ "resSeq2" ].value + ssbond_record[ "iCode2" ].value + ssbond_record[ "chainID2" ].value;

	// An old PDB standard would put two symmetry operations here;
	// we do not support this (yet?)

	ssbond.length = atof( ssbond_record[ "length" ].value.c_str() );  // bond length

	// If key is found in the links map, add this new linkage information to the links already keyed to this residue.
	if ( sfr.ssbond_map().count( ssbond.resID1 ) ) {
		ssbonds = sfr.ssbond_map()[ ssbond.resID1 ];
	}
	ssbonds.push_back( ssbond );

	sfr.ssbond_map()[ ssbond.resID1 ] = ssbonds;

	TR.Debug << "SSBOND record information stored successfully." << endl;
}


// Convert .pdb LINK record into SFR data.
void
store_link_record_in_sfr( Record link_record, StructFileRep & sfr )
{
	using namespace std;
	using namespace utility;
	using ObjexxFCL::stripped;

	LinkInformation link;
	vector1< LinkInformation > links;

	// If a LINK defines a connection across crystal cells, there is no way that we can handle this.
	// Ignore such lines.
	// (1555 is equivalent to a translation matrix of [0,0,0].)
	if ( ( ( link_record[ "sym1" ].value != "  1555" ) && ( link_record[ "sym1" ].value != "      " ) ) ||
			( ( link_record[ "sym2" ].value != "  1555" ) && ( link_record[ "sym2" ].value != "      " ) ) ) {
		TR.Debug << "Throwing out LINK across crystal cells." << endl;
		return;
	}

	// Extract values from LINK link_record fields.
	link.name1 = link_record[ "name1" ].value;  // 1st atom name
	link.resName1 = link_record[ "resName1" ].value;
	link.chainID1 = link_record[ "chainID1" ].value[ 0 ];
	link.resSeq1 = atof( link_record[ "resSeq1" ].value.c_str() );
	link.iCode1 = link_record[ "iCode1" ].value[ 0 ];

	//link.resID1 = stripped( link_record[ "resSeq1" ].value + link_record[ "iCode1" ].value + link_record[ "chainID1" ].value );
	link.resID1 = ( link_record[ "resSeq1" ].value + link_record[ "iCode1" ].value + link_record[ "chainID1" ].value );

	link.name2 = link_record[ "name2" ].value;  // 2nd atom name
	link.resName2 = link_record[ "resName2" ].value;
	link.chainID2 = link_record[ "chainID2" ].value[ 0 ];
	link.resSeq2 = atof( link_record[ "resSeq2" ].value.c_str() );
	link.iCode2 = link_record[ "iCode2" ].value[ 0 ];

	//link.resID2 = stripped( link_record[ "resSeq2" ].value + link_record[ "iCode2" ].value + link_record[ "chainID2" ].value );
	link.resID2 = ( link_record[ "resSeq2" ].value + link_record[ "iCode2" ].value + link_record[ "chainID2" ].value );

	link.length = atof( link_record[ "length" ].value.c_str() );  // bond length

	// If key is found in the links map, add this new linkage information to the links already keyed to this residue.
	if ( sfr.link_map().count( link.resID1 ) ) {
		links = sfr.link_map()[ link.resID1 ];
	}
	links.push_back( link );

	if ( links.size() > 1 ) {
		// The links need to be sorted such that higher-numbered residues come later.
		auto sort_func = []( LinkInformation const & lhs, LinkInformation const & rhs ) {
			return ( lhs.chainID2 < rhs.chainID2 ) || ( lhs.chainID2 == rhs.chainID2 && lhs.resSeq2 < rhs.resSeq2 );
		};
		sort( links.begin(), links.end(), sort_func );
	}

	sfr.link_map()[ link.resID1 ] = links;

	TR.Debug << "LINK record information stored successfully." << endl;
}

// Convert .pdb CISPEP record into SFR data.
void
store_cis_peptide_record_in_sfr( Record cispep_record, StructFileRep & sfr )
{
	using namespace std;

	CisPeptideInformation cis_pep;

	cis_pep.pep1 = cispep_record[ "pep1" ].value;
	cis_pep.chainID1 = cispep_record[ "chainID1" ].value[ 0 ];
	cis_pep.seqNum1 = atof( cispep_record[ "seqNum1" ].value.c_str() );
	cis_pep.icode1 = cispep_record[ "icode1" ].value[ 0 ];

	cis_pep.pep2 = cispep_record[ "pep2" ].value;
	cis_pep.chainID2 = cispep_record[ "chainID2" ].value[ 0 ];
	cis_pep.seqNum2 = atof( cispep_record[ "seqNum2" ].value.c_str() );
	cis_pep.icode2 = cispep_record[ "icode2" ].value[ 0 ];

	cis_pep.measure = atof( cispep_record[ "measure" ].value.c_str() );

	string const resID(
		cispep_record[ "seqNum1" ].value + cispep_record[ "icode1" ].value + cispep_record[ "chainID1" ].value );

	sfr.cispep_map()[ resID ] = cis_pep;
	TR.Debug << "CISPEP record information stored successfully." << endl;
}


// Convert .pdb CRYST1 record into SFR data.
void
store_crystallographic_parameter_record_in_sfr( Record crystal_record, StructFileRep & sfr )
{
	CrystInfo ci;
	ci.A( atof( crystal_record[ "a" ].value.c_str() ) );
	ci.B( atof( crystal_record[ "b" ].value.c_str() ) );
	ci.C( atof( crystal_record[ "c" ].value.c_str() ) );
	ci.alpha( atof( crystal_record[ "alpha" ].value.c_str() ) );
	ci.beta( atof( crystal_record[ "beta" ].value.c_str() ) );
	ci.gamma( atof( crystal_record[ "gamma" ].value.c_str() ) );
	ci.spacegroup( crystal_record[ "spacegroup" ].value );
	sfr.crystinfo() = ci;

	TR.Debug << "CRYST1 record information stored successfully." << std::endl;
}

// Parse and store unknown record types into SFR data.
/// @details  When a .pdb file is converted into a vector1 of records, unrecognized record types are labeled and stored
/// as UNKNOW records.  This function parses out Rosetta-specific information and stores it in the SFR.
void
store_unknown_records_in_sfr( utility::vector1< Record > unknown_records, StructFileRep & sfr )
{
	using namespace std;
	using namespace utility;

	// This is lame and hacky, but better than re-reading the whole file, which was happening before....  ~Labonte
	string line;
	bool reading_comments( false );
	Size const n_unknown_records( unknown_records.size() );
	TR << "Parsing " << n_unknown_records <<
		" .pdb records with unknown format to search for Rosetta-specific comments." << endl;
	for ( uint i( 1 ); i <= n_unknown_records; ++i ) {
		TR.Debug << unknown_records[ i ] << endl;
		line = create_pdb_line_from_record( unknown_records[ i ] );
		TR.Debug << line << endl;
		if ( startswith( line, "##Begin comments##" ) ) {
			TR.Debug << "Comments found; reading..." << endl;
			reading_comments = true;
		} else if ( startswith( line, "##End comments##" ) ) {
			TR.Debug << "Finished reading comments" << endl;
			reading_comments = false;
		} else if ( reading_comments ) {
			uint const space_location( line.find_first_of( " " ) );
			if ( space_location == string::npos ) { continue; }
			string const & key( line.substr( 0, space_location ) );
			string const & value( trim( line.substr( space_location + 1 ) ) );
			sfr.pdb_comments()[ key ] = value ;
		}
	}
}

}  // namespace pdb
}  // namespace io
}  // namespace core
