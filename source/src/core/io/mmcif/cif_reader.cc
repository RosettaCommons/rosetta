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
#include <core/io/pdb/build_pose_as_is.hh>
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
#include <utility/string_constants.hh>
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



using basic::Error;
using basic::Warning;
namespace core {
namespace io {
namespace mmcif {

static basic::Tracer TR( "core.io.mmcif.cif_reader" );


StructFileRepOP create_sfr_from_cif_file_op( CifFileOP cifFile, StructFileReaderOptions const & options ) {

	// NO TER OR END
	// INSERTION CODE DEFAULTS TO '?' SO TURN IT INTO A SPACE!.

	StructFileRepOP sfr( new StructFileRep );

	std::map<char, ChainAtoms> atom_chain_map;
	std::vector< char > chain_list; // preserve order
	std::map<char, Size> chain_to_idx;

	std::map<std::pair<Size, Size>, char> modelchain_to_chain;
	std::string const & chain_letters( utility::UPPERCASE_ALPHANUMERICS );
	for ( Size i = 0; i < chain_letters.size(); ++i ) {
		modelchain_to_chain[std::pair<Size, Size>(0, i)] = chain_letters[i];
		modelchain_to_chain[std::pair<Size, Size>(1, i)] = chain_letters[i];
	}
	Size modelidx = 1;
	bool modeltags_present = false;


	sfr->header() = HeaderInformationOP( new HeaderInformation() );

	bool read_pdb_header = options.read_pdb_header();

	TR.Debug << "Reading structure from CIF block " << cifFile->GetFirstBlockName() << std::endl;
	Block& block = cifFile->GetBlock( cifFile->GetFirstBlockName() );

	// "header information", i.e., is from the Title Section of the PDB file.
	if ( read_pdb_header && ! options.read_only_ATOM_entries() ) {
		if ( block.IsTablePresent( "citation" ) ) {
			ISTable& citation  = block.GetTable( "citation" );

			sfr->header()->store_title( citation( 0, "title" ) );
		}

		if ( block.IsTablePresent( "entry" ) ) {
			ISTable& entry = block.GetTable( "entry" );
			sfr->header()->store_idCode( entry( 0, "id" ) );
		}

		if ( block.IsTablePresent( "entity" ) ) {
			ISTable& entity = block.GetTable("entity");
			for ( Size i = 0; i <= entity.GetLastRowIndex(); ++i ) {
				sfr->header()->store_compound( entity( i, "pdbx_description" ) );
			}
		}

		if ( block.IsTablePresent( "keywords" ) ) {
			ISTable& keywords = block.GetTable( "struct_keywords" );
			sfr->header()->store_classification( keywords( 0, "pdbx_keywords" ) );
			sfr->header()->store_keywords( keywords( 0, "text" ) );
		}

		if ( block.IsTablePresent( "database_PDB_rev" ) ) {
			ISTable& database_PDB_rev = block.GetTable( "database_PDB_rev" );
			sfr->header()->store_deposition_date( database_PDB_rev( 0, "date_original" ) );
		}

		if ( block.IsTablePresent( "exptl" ) ) {
			ISTable& exptl = block.GetTable( "exptl" );
			sfr->header()->store_experimental_techniques( exptl( 0, "method" ) );
		}

		sfr->header()->finalize_parse();
	}


	// We cannot support REMARKs yet because these are stored in DIVERSE places.
	// There isn't a coherent "REMARKs" object. AMW TODO

	// HETNAM
	if ( block.IsTablePresent( "chem_comp" ) && ! options.read_only_ATOM_entries() ) {
		ISTable& chem_comp = block.GetTable("chem_comp");
		for ( Size i = 0; i <= chem_comp.GetLastRowIndex(); ++i ) {
			std::string name = chem_comp( i, "name" );
			string const & hetID( chem_comp( i, "id" ) );
			utility::trim( name );

			sfr->heterogen_names()[ hetID ] = name;

			pdb::store_base_residue_type_name_in_sfr( hetID, name, *sfr );
		}
	}

	// LINK
	if ( block.IsTablePresent( "struct_conn" ) && ! options.read_only_ATOM_entries() ) {
		ISTable& struct_conn = block.GetTable("struct_conn");
		for ( Size i = 0; i <= struct_conn.GetLastRowIndex(); ++i ) {
			if ( struct_conn( i, "conn_type_id" ) == "disulf" ) {
				SSBondInformation ssbond;
				utility::vector1< SSBondInformation > ssbonds;

				// Others: conn_type_id  for "covale" etc--alert JWL!

				// Extract values from record fields.
				//ssbond.name1 = struct_conn( i, "ptnr1_label_atom_id" );
				// Prefer 'author' annotations if available.
				if ( struct_conn.IsColumnPresent( "ptnr1_auth_comp_id" ) ) {
					ssbond.resName1 = struct_conn( i, "ptnr1_auth_comp_id" );
				} else {
					ssbond.resName1 = struct_conn( i, "ptnr1_label_comp_id" );
				}
				if ( struct_conn.IsColumnPresent( "ptnr1_auth_asym_id" ) ) {
					ssbond.chainID1 = struct_conn( i, "ptnr1_auth_asym_id" )[0];
				} else {
					ssbond.chainID1 = struct_conn( i, "ptnr1_label_asym_id" )[0];
				}
				if ( struct_conn.IsColumnPresent( "ptnr1_auth_seq_id" ) ) {
					ssbond.resSeq1 = atof( struct_conn( i, "ptnr1_auth_seq_id" ).c_str() );
				} else {
					ssbond.resSeq1 = atof( struct_conn( i, "ptnr1_label_seq_id" ).c_str() );
				}
				ssbond.iCode1 = struct_conn( i, "pdbx_ptnr1_PDB_ins_code" )[0] == '?' ? ' ' : struct_conn( i, "pdbx_ptnr1_PDB_ins_code" )[0];

				std::stringstream strstr1;
				strstr1 << std::setw( 4 ) << std::right << ssbond.resSeq1 << ssbond.iCode1 << ssbond.chainID1;
				ssbond.resID1 = strstr1.str();

				//ssbond.name2 = struct_conn( i, "ptnr2_label_atom_id" );
				if ( struct_conn.IsColumnPresent( "ptnr2_auth_comp_id" ) ) {
					ssbond.resName2 = struct_conn( i, "ptnr2_auth_comp_id" );
				} else {
					ssbond.resName2 = struct_conn( i, "ptnr2_label_comp_id" );
				}
				if ( struct_conn.IsColumnPresent( "ptnr2_auth_asym_id" ) ) {
					ssbond.chainID2 = struct_conn( i, "ptnr2_auth_asym_id" )[0];
				} else {
					ssbond.chainID2 = struct_conn( i, "ptnr2_label_asym_id" )[0];
				}
				if ( struct_conn.IsColumnPresent( "ptnr2_auth_seq_id" ) ) {
					ssbond.resSeq2 = atof( struct_conn( i, "ptnr2_auth_seq_id" ).c_str() );
				} else {
					ssbond.resSeq2 = atof( struct_conn( i, "ptnr2_label_seq_id" ).c_str() );
				}
				ssbond.iCode2 = struct_conn( i, "pdbx_ptnr2_PDB_ins_code" )[0] == '?' ? ' ' : struct_conn( i, "pdbx_ptnr2_PDB_ins_code" )[0];

				std::stringstream strstr2;
				strstr2 << std::setw( 4 ) << std::right << ssbond.resSeq2 << ssbond.iCode2 << ssbond.chainID2;
				ssbond.resID2 = strstr2.str();

				ssbond.length = atof( struct_conn( i, "pdbx_dist_value" ).c_str() ); // bond length

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

			} else {
				LinkInformation link;
				utility::vector1< LinkInformation > links;

				// hydrogen bonds should never be represented as LINKs
				if ( struct_conn( i, "conn_type_id" ) == "hydrog" ) continue;
				if ( struct_conn( i, "conn_type_id" ) == "saltbr" ) continue;
				// We treat metal coordination separately (thanks, -auto_setup_metals!)
				if ( struct_conn( i, "conn_type_id" ) == "metalc" ) continue;
				// Mismatched base pairs aren't treated in the input at all, but
				// conceivably RNA code might be interested
				if ( struct_conn( i, "conn_type_id" ) == "mismat" ) continue;
				// Hmm... try skipping.
				if ( struct_conn( i, "conn_type_id" ) == "modres" ) continue;

				// Extract values from record fields.
				link.name1 = struct_conn( i, "ptnr1_label_atom_id" );

				// Prefer 'author' annotations if available.
				if ( struct_conn.IsColumnPresent( "ptnr1_auth_comp_id" ) ) {
					link.resName1 = struct_conn( i, "ptnr1_auth_comp_id" );
				} else {
					link.resName1 = struct_conn( i, "ptnr1_label_comp_id" );
				}

				if ( struct_conn.IsColumnPresent( "ptnr1_auth_asym_id" ) ) {
					link.chainID1 = struct_conn( i, "ptnr1_auth_asym_id" )[0];
				} else {
					link.chainID1 = struct_conn( i, "ptnr1_label_asym_id" )[0];
				}

				if ( struct_conn.IsColumnPresent( "ptnr1_auth_seq_id" ) ) {
					link.resSeq1 = atof( struct_conn( i, "ptnr1_auth_seq_id" ).c_str() );
				} else {
					link.resSeq1 = atof( struct_conn( i, "ptnr1_label_seq_id" ).c_str() );
				}

				link.iCode1 = struct_conn( i, "pdbx_ptnr1_PDB_ins_code" )[0] == '?' ? ' ' : struct_conn( i, "pdbx_ptnr1_PDB_ins_code" )[0];

				//  LINK     type    1   6   name1           13  16   altLoc1       17  17   resName1      18  20   chainID1     22  22   resSeq1      23  26   iCode1      27  27   name2        43  46   altLoc2      47  47   resName2    48  50   chainID2    52  52   resSeq2     53  56   iCode2      57  57   sym1       60  65   sym2        67  72   length      74  78

				std::stringstream strstr1;
				strstr1 << std::setw( 4 ) << std::right << link.resSeq1 << link.iCode1 << link.chainID1;
				link.resID1 = strstr1.str();
				//struct_conn( i, "ptnr1_label_seq_id" )
				//+ ( struct_conn( i, "pdbx_ptnr1_PDB_ins_code" ) == "?" ) ? " " : struct_conn( i, "pdbx_ptnr1_PDB_ins_code" )
				//+ struct_conn( i, "ptnr1_label_asym_id" );

				link.name2 = struct_conn( i, "ptnr2_label_atom_id" );
				link.resName2 = struct_conn( i, "ptnr2_auth_comp_id" );
				link.chainID2 = struct_conn( i, "ptnr2_auth_asym_id" )[0];
				link.resSeq2 = atof( struct_conn( i, "ptnr2_auth_seq_id" ).c_str() );
				link.iCode2 = struct_conn( i, "pdbx_ptnr2_PDB_ins_code" )[0] == '?' ? ' ' : struct_conn( i, "pdbx_ptnr2_PDB_ins_code" )[0];

				std::stringstream strstr2;
				strstr2 << std::setw(4) << std::right << link.resSeq2 << link.iCode2 << link.chainID2;
				link.resID2 = strstr2.str();
				//struct_conn( i, "ptnr2_label_seq_id" )
				//+ ( struct_conn( i, "pdbx_ptnr2_PDB_ins_code" ) == "?" ) ? " " : struct_conn( i, "pdbx_ptnr2_PDB_ins_code" )
				//+ struct_conn( i, "ptnr2_label_asym_id" );

				link.length = struct_conn( i, "pdbx_dist_value" )[0] == '?' ? 0 : atof( struct_conn( i, "pdbx_dist_value" ).c_str() );  // bond length

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
	if ( block.IsTablePresent( "cell" ) && ! options.read_only_ATOM_entries() ) {
		ISTable& cell = block.GetTable("cell");
		CrystInfo ci;
		ci.A( atof( cell(0, "length_a" ).c_str() ) );
		ci.B( atof( cell(0, "length_b" ).c_str() ) );
		ci.C( atof( cell(0, "length_c" ).c_str() ) );
		ci.alpha( atof( cell(0, "angle_alpha" ).c_str() ) );
		ci.beta( atof( cell(0, "angle_beta" ).c_str() ) );
		ci.gamma( atof( cell(0, "angle_gamma" ).c_str() ) );
		ISTable& symmetry = block.GetTable("symmetry");
		ci.spacegroup( symmetry( 0, "space_group_name_H-M" ) );
		sfr->crystinfo() = ci;
	}

	// ATOM/HETATM

	std::string last_model = "";
	if ( block.IsTablePresent( "atom_site" ) ) {
		ISTable& atom_site = block.GetTable("atom_site");
		for ( Size i = 0; i <= atom_site.GetLastRowIndex(); ++i ) {
			AtomInformation ai;

			std::string temp_model = last_model;
			if ( atom_site.IsColumnPresent("pdbx_PDB_model_num") ) {
				temp_model = atom_site( i, "pdbx_PDB_model_num" );
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

			bool is_het = (atom_site.IsColumnPresent("group_PDB") && atom_site( i, "group_PDB" ) == "HETATM" );
			if ( is_het && options.read_only_ATOM_entries() ) continue;

			ai.isHet = is_het;
			ai.serial = atoi( atom_site( i, "id" ).c_str() ); // Mandatory
			if ( atom_site.IsColumnPresent("auth_atom_id") ) {
				ai.name = atom_site( i, "auth_atom_id" );
			} else if ( atom_site.IsColumnPresent("label_atom_id") ) {
				ai.name = atom_site( i, "label_atom_id");
			}
			ai.altLoc = 0;
			if ( atom_site( i, "label_alt_id" ).size() > 0 ) { // Mandatory
				ai.altLoc = atom_site( i, "label_alt_id" )[ 0 ];
			}

			if ( atom_site.IsColumnPresent("auth_comp_id") ) {
				ai.resName = atom_site( i, "auth_comp_id" );
			} else {
				ai.resName = atom_site( i, "label_comp_id"); // Mandatory
			}

			ai.chainID = ' ';
			if ( atom_site.IsColumnPresent("auth_asym_id") && atom_site( i, "auth_asym_id" ).size() > 0 ) {
				ai.chainID = atom_site( i, "auth_asym_id" )[0];
			} else if ( atom_site( i, "label_asym_id" ).size() > 0 ) {
				ai.chainID = atom_site( i, "label_asym_id" )[0]; // Mandatory
			}
			if ( options.new_chain_order() ) {
				char chainid = ai.chainID;
				if ( chain_to_idx.find(chainid) == chain_to_idx.end() ) {
					chain_to_idx[chainid] = chain_to_idx.size()+1;
					TR << "found new chain " << chainid << " " << chain_to_idx.size() << std::endl;
				}
				ai.chainID = modelchain_to_chain[std::pair<Size, Size>(modelidx, chain_to_idx[chainid])];
			}

			if ( atom_site.IsColumnPresent("auth_seq_id") ) {
				ai.resSeq = atoi( atom_site( i, "auth_seq_id" ).c_str() );
			} else {
				ai.resSeq = atoi( atom_site( i, "label_seq_id" ).c_str() ); // Mandatory
			}
			ai.iCode = ' ';
			if ( atom_site.IsColumnPresent("pdbx_PDB_ins_code") && atom_site( i, "pdbx_PDB_ins_code" ).size() > 0 && atom_site( i, "pdbx_PDB_ins_code" )[0] != '?' ) {
				ai.iCode = atom_site( i, "pdbx_PDB_ins_code" )[0];
			}

			// how can you check properly if something will successfully convert to a number !?!?!?
			bool force_no_occupancy = false;
			if ( atom_site( i, "Cartn_x" ) == "     nan" ) {
				ai.x =0.0;
				force_no_occupancy=true;
			} else {
				ai.x = atof( atom_site( i, "Cartn_x" ).c_str() );
			}
			if ( atom_site( i, "Cartn_y" ) == "     nan" ) {
				ai.y =0.0;
				force_no_occupancy=true;
			} else {
				ai.y = atof( atom_site( i, "Cartn_y" ).c_str() );
			}
			if ( atom_site( i, "Cartn_z" ) == "     nan" ) {
				ai.z =0.0;
				force_no_occupancy=true;
			} else {
				ai.z = atof( atom_site( i, "Cartn_z" ).c_str() );
			}

			// check that the occupancy column actually exists. If it doesn't, assume full occupancy.
			// otherwise read it.
			if ( ! atom_site.IsColumnPresent("occupancy") || atom_site( i, "occupancy" ) == "      " ) {
				ai.occupancy = 1.0;
			} else {
				ai.occupancy = atof( atom_site( i, "occupancy" ).c_str() );
			}
			if ( force_no_occupancy ) ai.occupancy = -1.0;

			if ( atom_site.IsColumnPresent("B_iso_or_equiv") ) {
				ai.temperature = atof( atom_site( i, "B_iso_or_equiv" ).c_str() );
			}
			ai.segmentID = "    ";
			ai.element = atom_site( i, "type_symbol" ); // Mandatory
			ai.terCount = 0;

			atom_chain_map[ai.chainID].push_back(ai);
			if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}
		}
	} else {
		TR.Warning << "atom_site table apparently not present in mmCIF file - Structure will not have coordinates!" << std::endl;
		std::vector<std::string> table_names;
		block.GetTableNames(table_names);
		TR.Warning << " Tables present:";
		for ( std::string const & table: table_names ) {
			TR.Warning << " " << table;
		}
		TR.Warning << std::endl;
		TR.Warning << " Note that the final table in the cif file may not be recognized - adding a dummy entry (like `_citation.title  \"\"`) to the end of the file may help." << std::endl;
	}

	for ( char i : chain_list ) { // std::vector
		sfr->chains().push_back( atom_chain_map.find( i )->second );
	}

	return sfr;
}

} // namespace mmcif
} // namespace io
} // namespace core
