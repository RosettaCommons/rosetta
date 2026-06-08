// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmcif/StructFileReader.hh
/// @author Andy Watkins (andy.watkins2@gmail.com)


// Unit headers
#include <core/io/mmcif/cif_reader.hh>
#include <core/io/StructFileReaderOptions.hh>

// Package headers
#include <core/io/pdb/pdb_reader.hh>  // TODO: Pull out pseudo-duplicated code and move to sfr_storage.cc.

// When you move PDBReader and PoseUnbuilder, take these.
#include <utility/string_util.hh>

#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>



// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/gemmi_util.hh>

// Numeric headers

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp>

#include <core/io/AtomInformation.hh> // AUTO IWYU For AtomInformation


// C++ headers
#include <cstdlib>
#include <algorithm>



namespace core {
namespace io {
namespace mmcif {

static basic::Tracer TR( "core.io.mmcif.cif_reader" );

using utility::find_gemmi_column;

StructFileRepOP create_sfr_from_cif_file( gemmi::cif::Document & cifdoc, StructFileReaderOptions const & options ) {
	// Note that pretty much every entry should be unwrapped with the appropriate accessor:
	using gemmi::cif::as_string; // Takes care of unquoting, use even if it's a simple string (e.g. atom names can have odd characters)
	using utility::as_char; // More robust verson
	using gemmi::cif::as_number; // Real
	using gemmi::cif::as_int; // Size

	// NO TER OR END
	// INSERTION CODE DEFAULTS TO '?' SO TURN IT INTO A SPACE!.

	StructFileRepOP sfr( new StructFileRep );

	std::map<std::string, ChainAtoms> atom_chain_map;
	std::vector< std::string > chain_list; // preserve order
	std::map<std::string, Size> chain_to_idx;

	std::map<std::pair<Size, Size>, char> modelchain_to_chain;
	std::string const & chain_letters( utility::UPPERCASE_ALPHANUMERICS );
	for ( Size i = 0; i < chain_letters.size(); ++i ) {
		modelchain_to_chain[std::pair<Size, Size>(0, i)] = chain_letters[i];
		modelchain_to_chain[std::pair<Size, Size>(1, i)] = chain_letters[i];
	}
	Size modelidx = 1;
	bool modeltags_present = false;


	sfr->header() = utility::pointer::make_shared< HeaderInformation >();

	bool read_pdb_header = options.read_pdb_header();

	if ( cifdoc.blocks.empty() ) {
		TR.Error << "Cannot interpret CIF file without any block information: `" << cifdoc.source << "`" << std::endl;
		return sfr; // Empty contents
	}

	gemmi::cif::Block & block = cifdoc.blocks[0];
	TR.Debug << "Reading structure from CIF block " << block.name << std::endl;

	// "header information", i.e., is from the Title Section of the PDB file.
	if ( read_pdb_header && ! options.read_only_ATOM_entries() ) {
		gemmi::cif::Table citation = block.find("_citation.", {"title"});
		if ( citation.size() > 0 ) {
			sfr->header()->store_title( as_string(citation[0][0]) );
		}

		gemmi::cif::Table entry = block.find("_entry.", {"id"});
		if ( entry.size() > 0 ) {
			sfr->header()->store_idCode( as_string(entry[0][0]) );
		}

		gemmi::cif::Table entity = block.find("_entity.", {"pdbx_description"});
		for ( Size ii = 0; ii < entity.size(); ++ii ) {
			sfr->header()->store_compound( as_string(entity[ii][0]) );
		}

		gemmi::cif::Table struct_keywords_pdbx_keywords = block.find("_struct_keywords.", {"pdbx_keywords"});
		if ( struct_keywords_pdbx_keywords.size() > 0 ) {
			sfr->header()->store_classification( as_string(struct_keywords_pdbx_keywords[0][0]) );
		}

		gemmi::cif::Table struct_keywords_text = block.find("_struct_keywords.", {"text"});
		if ( struct_keywords_pdbx_keywords.size() > 0 ) {
			sfr->header()->store_keywords( as_string(struct_keywords_text[0][0]) );
		}

		gemmi::cif::Table database_PDB_rev = block.find("_database_PDB_rev.", {"date_original"});
		if ( database_PDB_rev.size() > 0 ) {
			sfr->header()->store_deposition_date( as_string(database_PDB_rev[0][0]) );
		}

		gemmi::cif::Table exptl = block.find("_exptl.", {"method"});
		if ( exptl.size() > 0 ) {
			sfr->header()->store_experimental_techniques( as_string(exptl[0][0]) );
		}

		gemmi::cif::Table author = block.find("_audit_author.", {"name"});
		if ( author.size() ) {
			sfr->header()->store_authors( as_string(author[0][0]) );
		}

		sfr->header()->finalize_parse();
	}


	// We cannot support REMARKs yet because these are stored in DIVERSE places.
	// There isn't a coherent "REMARKs" object. AMW TODO

	// HETNAM
	if ( ! options.read_only_ATOM_entries() ) {
		gemmi::cif::Table chem_comp = block.find( "_chem_comp.", {"name","id"} );
		for ( Size ii = 0; ii < chem_comp.size(); ++ii ) {
			std::string name = as_string(chem_comp[ii][0]);
			std::string hetID = as_string(chem_comp[ii][1]);
			utility::trim( name );

			sfr->heterogen_names()[ hetID ] = name;

			pdb::store_base_residue_type_name_in_sfr( hetID, name, *sfr );
		}
	}

	// LINK
	if ( ! options.read_only_ATOM_entries() && block.has_mmcif_category("_struct_conn") ) {
		gemmi::cif::Table struct_conn = block.find_mmcif_category("_struct_conn");
		int conn_type_id = find_gemmi_column(struct_conn, "conn_type_id");
		int ptnr1_label_atom_id = find_gemmi_column(struct_conn, "ptnr1_label_atom_id");
		int ptnr1_auth_comp_id = find_gemmi_column(struct_conn, "ptnr1_auth_comp_id" );
		int ptnr1_auth_asym_id = find_gemmi_column(struct_conn, "ptnr1_auth_asym_id" );
		int ptnr1_auth_seq_id = find_gemmi_column(struct_conn, "ptnr1_auth_seq_id");
		int ptnr1_label_comp_id = find_gemmi_column(struct_conn, "ptnr1_label_comp_id");
		int ptnr1_label_asym_id = find_gemmi_column(struct_conn, "ptnr1_label_asym_id" );
		int ptnr1_label_seq_id = find_gemmi_column(struct_conn, "ptnr1_label_seq_id");
		int pdbx_ptnr1_PDB_ins_code = find_gemmi_column(struct_conn, "pdbx_ptnr1_PDB_ins_code");
		int ptnr2_auth_comp_id = find_gemmi_column(struct_conn, "ptnr2_auth_comp_id" );
		int ptnr2_auth_asym_id = find_gemmi_column(struct_conn, "ptnr2_auth_asym_id" );
		int ptnr2_auth_seq_id = find_gemmi_column(struct_conn, "ptnr2_auth_seq_id");
		int ptnr2_label_comp_id = find_gemmi_column(struct_conn, "ptnr2_label_comp_id");
		int ptnr2_label_asym_id = find_gemmi_column(struct_conn, "ptnr2_label_asym_id" );
		int ptnr2_label_seq_id = find_gemmi_column(struct_conn, "ptnr2_label_seq_id");
		int pdbx_ptnr2_PDB_ins_code = find_gemmi_column(struct_conn, "pdbx_ptnr2_PDB_ins_code");
		int pdbx_dist_value = find_gemmi_column(struct_conn, "pdbx_dist_value");

		for ( Size ii = 0; ii < struct_conn.size(); ++ii ) {
			gemmi::cif::Table::Row row = struct_conn[ii];
			if ( conn_type_id >= 0 && row[conn_type_id] == "disulf" ) {
				SSBondInformation ssbond;
				utility::vector1< SSBondInformation > ssbonds;

				// Others: conn_type_id  for "covale" etc--alert JWL!

				// Prefer 'author' annotations if available.
				if ( ptnr1_auth_comp_id >= 0 ) {
					ssbond.resName1 = as_string(row[ptnr1_auth_comp_id]);
				} else if ( ptnr1_label_comp_id >= 0 ) {
					ssbond.resName1 = as_string(row[ptnr1_label_comp_id]);
				} else {
					TR.Warning << "Can't find ptnr1_auth_comp_id or ptnr1_label_comp_id in disulfide annotation" << std::endl;
					continue;
				}
				if ( ptnr1_auth_asym_id >= 0 ) {
					ssbond.chainID1 = as_char(row[ptnr1_auth_asym_id], ' ');
				} else if ( ptnr1_label_asym_id >= 0 ) {
					ssbond.chainID1 = as_char(row[ptnr1_label_asym_id], ' ');
				} else {
					TR.Warning << "Can't find ptnr1_auth_asym_id or ptnr1_label_asym_id in disulfide annotation" << std::endl;
					continue;
				}
				if ( ptnr1_auth_seq_id >= 0 ) {
					ssbond.resSeq1 = as_int(row[ptnr1_auth_seq_id]);
				} else if ( ptnr1_label_seq_id >= 0 ) {
					ssbond.resSeq1 = as_int(row[ptnr1_label_seq_id]);
				} else {
					TR.Warning << "Can't find ptnr1_auth_seq_id or ptnr1_label_seq_id in disulfide annotation" << std::endl;
					continue;
				}

				if ( pdbx_ptnr1_PDB_ins_code >= 0 ) {
					ssbond.iCode1 = as_char(row[pdbx_ptnr1_PDB_ins_code], ' ');
				} else {
					ssbond.iCode1 = ' ';
				}

				ssbond.resID1 = ResID( ssbond.resSeq1, ssbond.iCode1, ssbond.chainID1 );

				if ( ptnr2_auth_comp_id >= 0 ) {
					ssbond.resName2 = as_string(row[ptnr2_auth_comp_id]);
				} else if ( ptnr2_label_comp_id >= 0 ) {
					ssbond.resName2 = as_string(row[ptnr2_label_comp_id]);
				} else {
					TR.Warning << "Can't find ptnr2_auth_comp_id or ptnr2_label_comp_id in disulfide annotation" << std::endl;
					continue;
				}
				if ( ptnr2_auth_asym_id >= 0 ) {
					ssbond.chainID2 = as_char(row[ptnr2_auth_asym_id], ' ');
				} else if ( ptnr2_label_asym_id >= 0 ) {
					ssbond.chainID2 = as_char(row[ptnr2_label_asym_id], ' ');
				} else {
					TR.Warning << "Can't find ptnr2_auth_asym_id or ptnr2_label_asym_id in disulfide annotation" << std::endl;
					continue;
				}
				if ( ptnr2_auth_seq_id >= 0 ) {
					ssbond.resSeq2 = as_int( row[ptnr2_auth_seq_id] );
				} else if ( ptnr2_label_seq_id >= 0 ) {
					ssbond.resSeq2 = as_int( row[ptnr2_label_seq_id] );
				} else {
					TR.Warning << "Can't find ptnr2_auth_seq_id or ptnr2_label_seq_id in disulfide annotation" << std::endl;
					continue;
				}

				if ( pdbx_ptnr2_PDB_ins_code >= 0 ) {
					ssbond.iCode2 = as_char(row[pdbx_ptnr2_PDB_ins_code], ' ');
				} else {
					ssbond.iCode2 = ' ';
				}

				ssbond.resID2 = ResID( ssbond.resSeq2, ssbond.iCode2, ssbond.chainID2 );

				if ( pdbx_dist_value >= 0 ) {
					ssbond.length = as_number( row[pdbx_dist_value] );
				} else {
					TR.Warning << "Can't find pdbx_dist_value in disulfide annotation" << std::endl;
					continue;
				}

				// If key is found in the links map, add this new linkage information to the links already keyed to this residue.
				if ( sfr->ssbond_map().count( ssbond.resID1 ) ) {
					ssbonds = sfr->ssbond_map()[ ssbond.resID1 ];
				}
				ssbonds.push_back( ssbond );

				if ( ssbonds.size() > 1 ) {
					// The ssbonds need to be sorted such that higher-numbered residues come later.
					auto sort_func = []( SSBondInformation const & lhs, SSBondInformation const & rhs ) {
						return ( lhs.chainID2 < rhs.chainID2 ) || ( lhs.chainID2 == rhs.chainID2 && lhs.resSeq2 < rhs.resSeq2 );
					};
					sort( ssbonds.begin(), ssbonds.end(), sort_func );
				}
				sfr->ssbond_map()[ ssbond.resID1 ] = ssbonds;

				if ( TR.Debug.visible() ) {
					TR.Debug << "SSBOND record information stored successfully." << std::endl;
				}

			} else { //
				LinkInformation link;
				utility::vector1< LinkInformation > links;
				if ( conn_type_id >= 0 && (
						// hydrogen bonds should never be represented as LINKs
						row[conn_type_id] == "hydrog" ||
						row[conn_type_id] == "saltbr" ||
						// We treat metal coordination separately (thanks, -auto_setup_metals!)
						row[conn_type_id] == "metalc" ||
						// Mismatched base pairs aren't treated in the input at all, but
						// conceivably RNA code might be interested
						row[conn_type_id] == "mismat" ||
						// Hmm... try skipping.
						row[conn_type_id] == "modres" )
						) {
					continue;
				}

				if ( ptnr1_label_atom_id >= 0 ) {
					link.name1 = as_string(row[ptnr1_label_atom_id]);
				} else {
					// ???
				}

				// Prefer 'author' annotations if available.
				if ( ptnr1_auth_comp_id >= 0 ) {
					link.resName1 = as_string(row[ptnr1_auth_comp_id]);
				} else if ( ptnr1_label_comp_id >= 0 ) {
					link.resName1 = as_string(row[ptnr1_label_comp_id]);
				} else {
					TR.Warning << "Can't find ptnr1_auth_comp_id or ptnr1_label_comp_id in LINK annotation" << std::endl;
					continue;
				}
				if ( ptnr1_auth_asym_id >= 0 ) {
					link.chainID1 = as_char(row[ptnr1_auth_asym_id], ' ');
				} else if ( ptnr1_label_asym_id >= 0 ) {
					link.chainID1 = as_char(row[ptnr1_label_asym_id], ' ');
				} else {
					TR.Warning << "Can't find ptnr1_auth_asym_id or ptnr1_label_asym_id in LINK annotation" << std::endl;
					continue;
				}
				if ( ptnr1_auth_seq_id >= 0 ) {
					link.resSeq1 = as_int(row[ptnr1_auth_seq_id]);
				} else if ( ptnr1_label_seq_id >= 0 ) {
					link.resSeq1 = as_int(row[ptnr1_label_seq_id]);
				} else {
					TR.Warning << "Can't find ptnr1_auth_seq_id or ptnr1_label_seq_id in LINK annotation" << std::endl;
					continue;
				}

				if ( pdbx_ptnr1_PDB_ins_code >= 0 ) {
					link.iCode1 = as_char(row[pdbx_ptnr1_PDB_ins_code], ' ');
				} else {
					link.iCode1 = ' ';
				}

				link.resID1 = ResID( link.resSeq1, link.iCode1, link.chainID1);

				if ( ptnr2_auth_comp_id >= 0 ) {
					link.resName2 = as_string(row[ptnr2_auth_comp_id]);
				} else if ( ptnr2_label_comp_id >= 0 ) {
					link.resName2 = as_string(row[ptnr2_label_comp_id]);
				} else {
					TR.Warning << "Can't find ptnr2_auth_comp_id or ptnr2_label_comp_id in LINK annotation" << std::endl;
					continue;
				}
				if ( ptnr2_auth_asym_id >= 0 ) {
					link.chainID2 = as_char(row[ptnr2_auth_asym_id], ' ');
				} else if ( ptnr2_label_asym_id >= 0 ) {
					link.chainID2 = as_char(row[ptnr2_label_asym_id], ' ');
				} else {
					TR.Warning << "Can't find ptnr2_auth_asym_id or ptnr2_label_asym_id in LINK annotation" << std::endl;
					continue;
				}
				if ( ptnr2_auth_seq_id >= 0 ) {
					link.resSeq2 = as_int(row[ptnr2_auth_seq_id]);
				} else if ( ptnr2_label_seq_id >= 0 ) {
					link.resSeq2 = as_int(row[ptnr2_label_seq_id]);
				} else {
					TR.Warning << "Can't find ptnr2_auth_seq_id or ptnr2_label_seq_id in LINK annotation" << std::endl;
					continue;
				}

				if ( pdbx_ptnr2_PDB_ins_code >= 0 ) {
					link.iCode2 = as_char(row[pdbx_ptnr2_PDB_ins_code], ' ');
				} else {
					link.iCode2 = ' ';
				}

				link.resID2 = ResID( link.resSeq2, link.iCode2, link.chainID2 );

				if ( pdbx_dist_value >= 0 ) {
					link.length = as_number( row[pdbx_dist_value] );
				} else {
					link.length = 0;
				}

				// If key is found in the links map, add this new linkage information to the links already keyed to this residue.
				if ( sfr->link_map().count( link.resID1 ) ) {
					links = sfr->link_map()[ link.resID1 ];
				}
				links.push_back( link );

				if ( links.size() > 1 ) {
					// The links need to be sorted such that higher-numbered residues come later.
					auto sort_func = []( LinkInformation const & lhs, LinkInformation const & rhs ) {
						return ( lhs.chainID2 < rhs.chainID2 ) || ( lhs.chainID2 == rhs.chainID2 && lhs.resSeq2 < rhs.resSeq2 );
					};
					sort( links.begin(), links.end(), sort_func );
				}

				sfr->link_map()[ link.resID1 ] = links;

				if ( TR.Debug.visible() ) {
					TR.Debug << "LINK record information stored successfully." << std::endl;
				}
			}
		}
	}

	// CRYST1
	if ( ! options.read_only_ATOM_entries() ) {
		gemmi::cif::Table cell = block.find("_cell.", {"length_a","length_b","length_c","angle_alpha","angle_beta","angle_gamma"});
		gemmi::cif::Table symmetry = block.find("_symmetry.", {"space_group_name_H-M"});
		if ( cell.size() > 0 && symmetry.size() > 0 ) {
			CrystInfo ci;
			ci.A( as_number( cell[0][0] ) );
			ci.B( as_number( cell[0][1] ) );
			ci.C( as_number( cell[0][2] ) );
			ci.alpha( as_number( cell[0][3] ) );
			ci.beta(  as_number( cell[0][4] ) );
			ci.gamma( as_number( cell[0][5] ) );
			ci.spacegroup( as_string(symmetry[0][0]) );
			sfr->crystinfo() = ci;
		}
	}

	// ATOM/HETATM

	std::string last_model = "";
	if ( block.has_mmcif_category("_atom_site") ) {
		gemmi::cif::Table atom_site = block.find_mmcif_category("_atom_site");

		int pdbx_PDB_model_num  = find_gemmi_column(atom_site, "pdbx_PDB_model_num");
		int group_PDB = find_gemmi_column(atom_site, "group_PDB" );
		int id = find_gemmi_column(atom_site, "id" );
		int auth_atom_id = find_gemmi_column(atom_site, "auth_atom_id" );
		int label_atom_id = find_gemmi_column(atom_site, "label_atom_id" );
		int label_alt_id = find_gemmi_column(atom_site, "label_alt_id" );
		int auth_comp_id = find_gemmi_column(atom_site, "auth_comp_id");
		int label_comp_id = find_gemmi_column(atom_site, "label_comp_id");
		int auth_asym_id  = find_gemmi_column(atom_site, "auth_asym_id");
		int label_asym_id = find_gemmi_column(atom_site, "label_asym_id");
		int auth_seq_id = find_gemmi_column(atom_site, "auth_seq_id" );
		int label_seq_id = find_gemmi_column(atom_site, "label_seq_id" );
		int pdbx_PDB_ins_code = find_gemmi_column(atom_site, "pdbx_PDB_ins_code" );
		int Cartn_x = find_gemmi_column(atom_site, "Cartn_x" );
		int Cartn_y = find_gemmi_column(atom_site, "Cartn_y" );
		int Cartn_z = find_gemmi_column(atom_site, "Cartn_z" );
		int occupancy = find_gemmi_column(atom_site, "occupancy" );
		int B_iso_or_equiv = find_gemmi_column(atom_site, "B_iso_or_equiv" );
		int type_symbol = find_gemmi_column(atom_site, "type_symbol" );

		for ( Size ii = 0; ii < atom_site.size(); ++ii ) {
			gemmi::cif::Table::Row row = atom_site[ii];

			AtomInformation ai;
			std::string temp_model = last_model;

			if ( pdbx_PDB_model_num >= 0 ) {
				temp_model = as_string( row[pdbx_PDB_model_num] );
				temp_model = ObjexxFCL::strip_whitespace( temp_model );
			}

			// PUT OBEY ENDMDL RECOGNITION HERE!!!!!!!
			if ( last_model != "" && last_model != temp_model ) {
				break;
			}
			last_model = temp_model;

			// store the serial number as the filename, which will become the PDBInfo name of the pose
			sfr->modeltag() = temp_model;
			if ( options.new_chain_order() ) {
				if ( modeltags_present ) {
					if ( options.obey_ENDMDL() ) {
						// That's enough records for now.
						break;
					}
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

			bool is_het = (group_PDB >= 0 && as_string(row[group_PDB]) == "HETATM" );
			if ( is_het && options.read_only_ATOM_entries() ) continue;

			ai.isHet = is_het;
			ai.serial = as_int( row[id], 0 ); // Mandatory
			if ( auth_atom_id >= 0 ) {
				ai.name = as_string( row[auth_atom_id] );
			} else if ( label_atom_id >= 0 ) {
				ai.name = as_string( row[label_atom_id] );
			}
			ai.altLoc = 0;
			if ( label_alt_id >=0 ) {
				ai.altLoc = as_char(row[label_alt_id], 0);
			}

			if ( auth_comp_id >= 0 ) {
				ai.resName = as_string(row[auth_comp_id]);
			} else {
				ai.resName = as_string(row[label_comp_id]); // Mandatory
			}

			ai.chainID = " ";
			if ( auth_asym_id >= 0 ) {
				ai.chainID = as_string(row[auth_asym_id]);
			} else if ( label_asym_id >= 0 ) {
				ai.chainID = as_string(row[label_asym_id]); // Mandatory
			}
			if ( options.new_chain_order() ) {
				std::string chainid = ai.chainID;
				if ( chain_to_idx.find(chainid) == chain_to_idx.end() ) {
					chain_to_idx[chainid] = chain_to_idx.size()+1;
					TR << "found new chain " << chainid << " " << chain_to_idx.size() << std::endl;
				}
				ai.chainID = modelchain_to_chain[std::pair<Size, Size>(modelidx, chain_to_idx[chainid])];
			}

			if ( auth_seq_id >= 0 ) {
				ai.resSeq = as_int( row[auth_seq_id], 0 );
			} else {
				ai.resSeq = as_int( row[label_seq_id], 0 ); // Mandatory
			}
			ai.iCode = ' ';
			if ( pdbx_PDB_ins_code >= 0 ) {
				ai.iCode = as_char( row[pdbx_PDB_ins_code], ' ');
			}

			bool force_no_occupancy = false;
			ai.x = as_number(row[Cartn_x]);
			if ( std::isnan(ai.x) ) {
				ai.x =0.0;
				force_no_occupancy=true;
			}
			ai.y = as_number(row[Cartn_y]);
			if ( std::isnan(ai.y) ) {
				ai.y =0.0;
				force_no_occupancy=true;
			}
			ai.z = as_number(row[Cartn_z]);
			if ( std::isnan(ai.z) ) {
				ai.z =0.0;
				force_no_occupancy=true;
			}

			// check that the occupancy column actually exists. If it doesn't, assume full occupancy.
			// otherwise read it.
			if ( occupancy < 0 ) {
				ai.occupancy = 1.0;
			} else {
				ai.occupancy = as_number( row[occupancy] );
				// On error
			}
			if ( force_no_occupancy || std::isnan(ai.occupancy) ) ai.occupancy = -1.0;

			if ( B_iso_or_equiv >= 0 ) {
				ai.temperature = as_number( row[B_iso_or_equiv] );
			}
			ai.segmentID = "    ";
			ai.element = as_string( row[type_symbol] ); // Mandatory
			ai.terCount = 0;

			atom_chain_map[ai.chainID].push_back(ai);
			if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}
		}
	} else {
		TR.Warning << "atom_site table apparently not present in mmCIF file - Structure will not have coordinates!" << std::endl;
		//  std::vector<std::string> table_names;
		//  block.GetTableNames(table_names);
		//  TR.Warning << " Tables present:";
		//  for ( std::string const & table: table_names ) {
		//   TR.Warning << " " << table;
		//  }
		//  TR.Warning << std::endl;
		//  TR.Warning << " Note that the final table in the cif file may not be recognized - adding a dummy entry (like `_citation.title  \"\"`) to the end of the file may help." << std::endl;
	}

	for ( std::string const & i : chain_list ) { // std::vector
		sfr->chains().push_back( atom_chain_map.find( i )->second );
	}

	return sfr;
}

} // namespace mmcif
} // namespace io
} // namespace core
