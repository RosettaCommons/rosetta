// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// comments in header file!
/// @file   src/core/chemical/sdf/mmCIFParser.cc
/// @author Steven Combs


#include <core/chemical/mmCIF/mmCIFParser.hh>
#include <core/types.hh>
#include <string>
#include <core/chemical/sdf/MolFileIOData.hh>

//Utility functions
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>
#include <utility/gemmi_util.hh>

//external CIF includes
#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp> // for as_number

namespace core {
namespace chemical {
namespace mmCIF {

//Load up the tracer for this class
static basic::Tracer TR( "core.io.mmCIF.mmCIFParser" );

using utility::find_gemmi_column;

mmCIFParser::mmCIFParser() {
	bond_string_to_sdf_size_[ "SING"] = 1;
	bond_string_to_sdf_size_[ "DOUB"] = 2;
	bond_string_to_sdf_size_[ "TRIP"] = 3;
}

sdf::MolFileIOMoleculeOP
mmCIFParser::parse(
	std::string const & lines,
	std::string const & pdb_id
) {
	sdf::MolFileIOMoleculeOP molecule;

	try {
		gemmi::cif::Document cifdoc = gemmi::cif::read_memory(lines.c_str(), lines.size(), "CIFFILE" );

		for ( gemmi::cif::Block & block: cifdoc.blocks ) {
			//You have to add the "#" string to the pdb_id. This is because of how
			//the cifFile parser interpets the lines. For whatever reason, it adds
			//a # to the data_??? block, where ??? is the 3 letter code for the
			//file
			if ( block.name == pdb_id + "#" ) {
				molecule = get_molfile_molecule( block );
				molecule->name( pdb_id );
				return molecule;
			}
		}
	} catch (std::runtime_error const & e) { // gemmi errors
		TR.Error << "Error reading Residue Type definition CIF file: " << e.what() << std::endl;
		return molecule;
	}
	return molecule;
}

utility::vector1< sdf::MolFileIOMoleculeOP>
mmCIFParser::parse(std::string const &filename){
	utility::vector1< sdf::MolFileIOMoleculeOP> molecules;

	try {
		gemmi::cif::Document cifdoc = gemmi::cif::read_file(filename);
		for ( gemmi::cif::Block & block: cifdoc.blocks ) {
			molecules.push_back( get_molfile_molecule( block ) );
		}
	} catch (std::runtime_error const & e) { // gemmi errors
		TR.Error << "Error reading Residue Type definition CIF file `" << filename << "`: " << e.what() << std::endl;
	}

	return molecules;
}

sdf::MolFileIOMoleculeOP
mmCIFParser::get_molfile_molecule( gemmi::cif::Block & block ) {
	// Note that pretty much every entry should be unwrapped with the appropriate accessor:
	using gemmi::cif::as_string; // Takes care of unquoting, use even if it's a simple string (e.g. atom names can have odd characters)
	using utility::as_char; // More robust version
	using gemmi::cif::as_number; // Real
	using gemmi::cif::as_int; // Size

	sdf::MolFileIOMoleculeOP molecule( new sdf::MolFileIOMolecule() );

	molecule->name( block.name );
	//only proceed if the tables for bonds and atoms are present
	if ( !block.has_mmcif_category("_chem_comp_atom") ) {
		TR.Error << "Cannot parse CIF file. No atom block (chem_comp_atom) found for " << block.name << std::endl;
		return molecule;
	}
	gemmi::cif::Table atom_comp = block.find_mmcif_category("_chem_comp_atom");
	if ( atom_comp.size() == 0 ) {
		TR.Error << "Cannot parse CIF file. Empty atom block (chem_comp_atom) found for " << block.name << std::endl;
		return molecule;
	}

	int type_symbol = find_gemmi_column(atom_comp,"type_symbol");
	if ( type_symbol < 0 ) {
		TR.Error << "Cannot parse CIF file. Missing element information in atom block (chem_comp_atom) for " << block.name << std::endl;
		return molecule;
	}

	int pdbx_model_Cartn_x_ideal = find_gemmi_column(atom_comp,"pdbx_model_Cartn_x_ideal");
	int pdbx_model_Cartn_y_ideal = find_gemmi_column(atom_comp,"pdbx_model_Cartn_y_ideal");
	int pdbx_model_Cartn_z_ideal = find_gemmi_column(atom_comp,"pdbx_model_Cartn_z_ideal");
	int model_Cartn_x = find_gemmi_column(atom_comp,"model_Cartn_x");
	int model_Cartn_y = find_gemmi_column(atom_comp,"model_Cartn_y");
	int model_Cartn_z = find_gemmi_column(atom_comp,"model_Cartn_z");
	int charge = find_gemmi_column(atom_comp,"charge");

	int atom_name_id = find_gemmi_column(atom_comp,"atom_id");
	if ( atom_name_id < 0 ) {
		atom_name_id = find_gemmi_column(atom_comp,"pdbx_component_atom_id");
	}
	if ( atom_name_id < 0 ) {
		TR.Error << "Can't find atom id column (atom_id/pdbx_component_atom_id) in chem_comp_atom table for " << block.name << std::endl;
		return molecule;
	}

	// A map of atom names to their chemical symbols -- includes ignored atoms
	std::map< std::string, std::string > name_to_element_map;

	for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
		gemmi::cif::Table::Row row = atom_comp[ii];

		//set atom name
		std::string atom_name = as_string( row[atom_name_id] );
		name_to_element_map[ atom_name ] = as_string( row[type_symbol] );
	}

	// Pick up if we need/want to treat this as a polymeric type
	bool is_peptide_linking = false;
	bool is_nucleic_linking = false;
	gemmi::cif::Table chem_comp = block.find( "_chem_comp.", {"type"} );
	if ( chem_comp.size() > 0 ) {
		std::string type = utility::uppercased( as_string(chem_comp[0][0]) );
		// TODO: Need better handling of variants (e.g. 1ZN, B5I)
		if ( type == "L-PEPTIDE LINKING" ) {
			TR.Debug << "Found L-peptide RT" << std::endl;// named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "PROTEIN POLYMER L_AA" );
			is_peptide_linking = true;
			is_nucleic_linking = false;
		} else if ( type == "D-PEPTIDE LINKING" ) {
			TR.Debug << "Found D-peptide RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "PROTEIN POLYMER D_AA" );
			is_peptide_linking = true;
			is_nucleic_linking = false;
		} else if ( type == "RNA LINKING" ) {
			TR.Debug << "Found D-RNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "RNA POLYMER D_RNA" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else if ( type == "L-RNA LINKING" ) {
			TR.Debug << "Found L-RNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "RNA POLYMER L_RNA" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		}  else if ( type == "DNA LINKING" ) {
			TR.Debug << "Found D-DNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "DNA POLYMER" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else if ( type == "L-DNA LINKING" ) {
			TR.Debug << "Found L-DNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "DNA POLYMER" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else {
			is_peptide_linking = false;
			is_nucleic_linking = false;
		}
	}

	//this is to map the atom names to a core::Size value (id). This is added when adding
	//bonds to the MolIO object.
	std::map< std::string, core::Size > atom_name_to_id;

	//prefer the ideal coordinates, but if not found, use cartesian coordinates
	int x_id = pdbx_model_Cartn_x_ideal, y_id = pdbx_model_Cartn_y_ideal, z_id = pdbx_model_Cartn_z_ideal;
	if ( x_id < 0 || y_id < 0 || z_id < 0 || gemmi::cif::is_null( atom_comp[0][x_id] ) ) {
		x_id = model_Cartn_x;
		y_id = model_Cartn_y;
		z_id = model_Cartn_z;
	}
	if ( x_id < 0 || y_id < 0 || z_id < 0 || gemmi::cif::is_null( atom_comp[0][z_id] ) ) {
		utility_exit_with_message( "No usable coordinates for mmCIF file for " + block.name );
	}

	std::set< std::string > atoms_to_ignore = get_atoms_to_ignore(block, atom_name_id, name_to_element_map, is_peptide_linking, is_nucleic_linking );

	//start atom block
	Size index = 1;
	for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
		sdf::MolFileIOAtomOP atom( new sdf::MolFileIOAtom());
		//atom id is whatever the number we are on +1
		// AMW: meaning... that we don't need the below, nor would we actually want it
		// because needing to skip (due to non-patched polymers) screws it up.
		// ditto our index, ii... so use another.
		atom->index( index );

		//set atom name
		std::string atom_name = as_string( atom_comp[ ii][ atom_name_id ] );
		TR.Trace << "Examining atom entry " << atom_name << std::endl;
		if ( atoms_to_ignore.count( atom_name ) ) {
			TR.Trace << "Ignoring." << std::endl;
			continue;
		}

		atom->name( atom_name );

		//set map to index
		atom_name_to_id[ atom_name ] = index;
		//set element name
		atom->element( as_string( atom_comp[ii][type_symbol] ) );

		//get the xyz cordinates
		core::Real x = as_number( atom_comp[ii][ x_id ] );
		core::Real y = as_number( atom_comp[ii][ y_id ] );
		core::Real z = as_number( atom_comp[ii][ z_id ] );
		//set xyz coordinates
		atom->position( core::Vector( x, y, z ) );

		if ( charge >= 0 ) {
			atom->formal_charge( as_int( atom_comp[ii][charge], 0 ) ); // Default zero if present and null
		} else {
			atom->formal_charge( 0 );
		}

		molecule->add_atom( atom );
		// only increment if we actually get here.
		index++;
	}

	gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2","value_order"} );
	if ( bond_comp.size() == 0 && atom_comp.size() > 1 ) {
		TR.Error << "Cannot parse CIF file. No bond block (chem_comp_bond) found for multi-atom entry " << block.name << std::endl;
	} // else one atom entry without bond block

	for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
		std::string source = as_string( bond_comp[ii][0] ); //atom 1
		std::string target = as_string( bond_comp[ii][1] ); //atom 2
		core::Size bond_type( bond_string_to_sdf_size_[ as_string(bond_comp[ii][2]) ] );

		sdf::MolFileIOBondOP bond( new sdf::MolFileIOBond() );

		TR.Trace << "Examining bond entry " << source << " " << target << std::endl;

		if ( atoms_to_ignore.count(source) > 0 || atoms_to_ignore.count(target) > 0 ) {
			continue;
		}

		TR.Trace << "Keeping bond entry " << source << " " << target << std::endl;

		bond->atom1( atom_name_to_id[ source ] );
		bond->atom2( atom_name_to_id[ target ] );

		bond->sdf_type( bond_type); //bond order

		molecule->add_bond( bond);
	}

	annotate_polymeric_connections(*molecule, block, atom_name_id, name_to_element_map, is_peptide_linking, is_nucleic_linking);
	TR.Debug << "LOWER: `" << molecule->get_lower_atom() << "` UPPER: `" << molecule->get_upper_atom() << "`" << std::endl;

	return molecule;
}

/// Utilities for get_atoms_to_ignore() and annotate_polymeric_connections()

template< class C >
utility::vector1< std::string >
find_elements( C const & container, std::string const & elem, std::map< std::string, std::string > const & name_to_element_map ) {
	utility::vector1< std::string > found;
	for ( std::string const & atm: container ) {
		if ( name_to_element_map.at(atm) == elem ) {
			found.push_back( atm );
		}
	}
	return found;
}

template< class C >
utility::vector1< std::string >
find_heavy( C const & container, std::map< std::string, std::string > const & name_to_element_map ) {
	utility::vector1< std::string > found;
	for ( std::string const & atm: container ) {
		if ( name_to_element_map.at(atm) != "H" ) {
			found.push_back( atm );
		}
	}
	return found;
}

utility::vector1< std::string >
get_attached_atoms( std::string const & atm, gemmi::cif::Block& block ) {
	utility::vector1< std::string > found;
	gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
	for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
		std::string source = gemmi::cif::as_string( bond_comp[ii][0] ); //atom 1
		std::string target = gemmi::cif::as_string( bond_comp[ii][1] ); //atom 2
		if ( source == atm ) {
			found.push_back( target );
		}
		if ( target == atm ) {
			found.push_back( source );
		}
	}
	return found;
}

std::set< std::string >
mmCIFParser::get_atoms_to_ignore(
	gemmi::cif::Block& block,
	int atom_name_id,
	std::map< std::string, std::string > const & name_to_element_map,
	bool is_peptide_linking,
	bool is_nucleic_linking
) {
	using gemmi::cif::as_string; // Takes care of unquoting, use even if it's a simple string (e.g. atom names can have odd characters)

	std::set< std::string > atoms_to_ignore;
	if ( !is_peptide_linking && !is_nucleic_linking) {
		// Ligands keep all the atoms
		return atoms_to_ignore;
	}

	gemmi::cif::Table atom_comp = block.find_mmcif_category("_chem_comp_atom");
	int leaving_flag = find_gemmi_column(atom_comp,"pdbx_leaving_atom_flag");
	if ( leaving_flag >= 0 ) {
		for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
			gemmi::cif::Table::Row row = atom_comp[ii];

			if ( as_string(row[leaving_flag]) == "Y" ) {
				std::string atom_name = as_string(row[atom_name_id]);
				TR.Trace << "Ignoring " << atom_name << " as it's flagged as a leaving atom." << std::endl;
				atoms_to_ignore.insert( atom_name );
			}
		}
		// If the annotation is valid, trust it
		// Though there are some edge cases (e.g. BVP) where hydrogens attached to leaving heavyatoms aren't annotated as leaving -- catch this
		for ( std::string const & to_ignore: find_heavy(atoms_to_ignore, name_to_element_map) ) {
			for ( std::string const & hydro: find_elements(get_attached_atoms(to_ignore, block), "H", name_to_element_map) ) {
				if ( atoms_to_ignore.count(hydro) == 0 ) {
					TR.Trace << "Ignoring " << hydro << " as a hydrogen connected to a heavy atom being ignored." << std::endl;
					atoms_to_ignore.insert( hydro ); // Safe, as find_heavy() above means we're not iterating over this continer directly
				}
			}
		}
		return atoms_to_ignore;
	}

	////////////////////////////////////////
	// Now we're into fall-back heuristics

	gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );

	if ( is_peptide_linking ) {
		TR << "Using fall-back heuristics for ignored peptide atoms" << std::endl;
		std::set< std::string > no_hydro = {"OXT", "N" };
		core::Size n_heavy_to_OXT = 0;

		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			if ( no_hydro.count(source) > 0  ) {
				if ( name_to_element_map.at(target) == "H" ) {
					TR.Trace << "Ignoring hydrogen " << target << " due to its bond to " << source << std::endl;
					atoms_to_ignore.insert( target );
				}
			}
			if ( no_hydro.count(target) > 0  ) {
				if ( name_to_element_map.at(source) == "H" ) {
					TR.Trace << "Ignoring hydrogen " << source << " due to its bond to " << target << std::endl;
					atoms_to_ignore.insert( source );
				}
			}

			if ( source == "OXT" ) {
				if ( name_to_element_map.at(target) != "H" ) {
					n_heavy_to_OXT += 1;
				}
			} else if ( target == "OP3" || target == "O3P" ) {
				if ( name_to_element_map.at(source) != "H" ) {
					n_heavy_to_OXT += 1;
				}
			}
		}

		if ( n_heavy_to_OXT < 2	) {
			atoms_to_ignore.insert( "OXT" );
		} else {
			TR.Trace << "Not ignoring OXT, as it contains a pendant group." << std::endl;
		}
	}

	if ( is_nucleic_linking ) {
		TR << "Using fall-back heuristics for ignored nucleotide atoms" << std::endl;
		std::set< std::string > no_hydro = {"P", "O3'", "OP3", "O3P", "OP2", "O2P", "OP1", "O1P"  };

		core::Size n_heavy_to_O3P = 0;

		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			if ( no_hydro.count(source) > 0 ) {
				if ( name_to_element_map.at(target) == "H" ) {
					TR.Trace << "Ignoring hydrogen " << target << " due to its bond to " << source << std::endl;
					atoms_to_ignore.insert( target );
				}
			}
			if ( no_hydro.count(target) > 0 ) {
				if ( name_to_element_map.at(source) == "H" ) {
					TR.Trace << "Ignoring hydrogen " << target << " due to its bond to " << target << std::endl;
					atoms_to_ignore.insert( source );
				}
			}

			if ( source == "OP3" || source == "O3P" ) {
				if ( name_to_element_map.at(target) != "H" ) {
					n_heavy_to_O3P += 1;
				}
			} else if ( target == "OP3" || target == "O3P" ) {
				if ( name_to_element_map.at(source) != "H" ) {
					n_heavy_to_O3P += 1;
				}
			}
		}

		if ( n_heavy_to_O3P < 2	) {
			atoms_to_ignore.insert( "OP3" );
			atoms_to_ignore.insert( "O3P" );
		} else {
			TR.Trace << "Not ignoring OP3/O3P, as it contains a pendant group." << std::endl;
		}
	}

	return atoms_to_ignore;
}

void
mmCIFParser::annotate_polymeric_connections(
	sdf::MolFileIOMolecule & molecule,
	gemmi::cif::Block& block,
	int atom_name_id,
	std::map< std::string, std::string > const & name_to_element_map,
	bool is_peptide_linking,
	bool is_nucleic_linking
) {
	if ( !is_peptide_linking && !is_nucleic_linking ) {
		return; // non-polymeric, return early
	}

	using gemmi::cif::as_string; // Takes care of unquoting, use even if it's a simple string (e.g. atom names can have odd characters)

	std::set< std::string > leaving, backbone, n_term, c_term;
	gemmi::cif::Table atom_comp = block.find_mmcif_category("_chem_comp_atom");
	int leaving_flag = find_gemmi_column(atom_comp,"pdbx_leaving_atom_flag");
	int backbone_flag = find_gemmi_column(atom_comp,"pdbx_backbone_atom_flag");
	int n_term_flag = find_gemmi_column(atom_comp,"pdbx_n_terminal_atom_flag");
	int c_term_flag = find_gemmi_column(atom_comp,"pdbx_c_terminal_atom_flag");

	for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
		gemmi::cif::Table::Row row = atom_comp[ii];

		if ( leaving_flag >= 0 && as_string(row[leaving_flag]) == "Y" ) {
			leaving.insert( as_string(row[atom_name_id]) );
		}
		if ( backbone_flag >= 0 && as_string(row[backbone_flag]) == "Y" ) {
			backbone.insert( as_string(row[atom_name_id]) );
		}
		if ( n_term_flag >= 0 && as_string(row[n_term_flag]) == "Y" ) {
			n_term.insert( as_string(row[atom_name_id]) );
		}
		if ( c_term_flag >= 0 && as_string(row[c_term_flag]) == "Y" ) {
			c_term.insert( as_string(row[atom_name_id]) );
		}
	}

	if ( is_peptide_linking ) {
		std::string n_term_atm;

		/////////////////// PEPTIDE LOWER

		// TODO: Do we always want to fall back, or should that be conditional based on whether particular columns are or are not found?
		if ( !n_term.empty() ) {
			// Ideally we'd like a nitrogen -- if there's multiple, just pick the first.
			utility::vector1< std::string > const & found_N = find_elements( n_term, "N", name_to_element_map );
			TR.Trace << "Found annotated N-terminal nitrogens: " << found_N << std::endl;
			if ( found_N.size() >= 1 ) {
				n_term_atm = found_N[1];
				TR.Debug << "Picking N-terminal connection point as (first) N-term annotated nitrogen: " << n_term_atm << std::endl;
			}
			if ( n_term_atm.empty() ) {
				// Failing that, pick the first heavy atom.
				utility::vector1< std::string > const & found_heavy = find_heavy( n_term, name_to_element_map );
				TR.Trace << "Found annotated N-terminal heavy: " << found_heavy << std::endl;
				if ( found_heavy.size() >= 1 ) {
					n_term_atm = found_heavy[1];
					TR.Debug << "Picking N-terminal connection point as (first) N-term annotated heavy atom: " << n_term_atm << std::endl;
				}
			}
		}
		if ( n_term_atm.empty() && !leaving.empty() ) {
			// Next, try to find a (unique) nitrogen attached to a leaving hydrogen.
			utility::vector1< std::string > found_N;
			for ( std::string const & atm: leaving ) {
				if ( name_to_element_map.at(atm) == "H" ) {
					found_N.append( find_elements( get_attached_atoms( atm, block ), "N", name_to_element_map ) );
				}
			}
			std::set< std::string > found_N_set( found_N.begin(), found_N.end() ); // Need to deduplicate
			TR.Trace << "Found nitrogens with attached leaving Hs: " << found_N << std::endl;
			if ( found_N_set.size() == 1 ) {
				n_term_atm = *(found_N.begin());
				TR.Debug << "Picking N-terminal connection point as unique nitrogen attached to leaving hydrogen: " << n_term_atm << std::endl;
			}
		}
		if ( n_term_atm.empty() && !backbone.empty() ) {
			// Next, try to find a (unique) nitrogen in the annotated backbone.
			utility::vector1< std::string > const & found_N = find_elements( backbone, "N", name_to_element_map );
			TR.Trace << "Found Backbone nitrogens: " << found_N << std::endl;
			if ( found_N.size() == 1 ) {
				n_term_atm = found_N[1];
				TR.Debug << "Picking N-terminal connection point as unique backbone nitrogen: " << n_term_atm << std::endl;
			}
		}
		if ( n_term_atm.empty() ) {
			// Fallback -- look for an atom named 'N'
			if ( name_to_element_map.count("N") ) {
				n_term_atm = "N";
				TR.Debug << "Picking N-terminal connection point based on name: " << n_term_atm << std::endl;
			}
		}

		if ( !n_term_atm.empty() ) {
			TR.Debug << "Setting LOWER atom to " << n_term_atm << std::endl;
			molecule.set_lower_atom( n_term_atm );
			// The attached hydrogens to N are generally for the free molecule, and not in the correct geometry for the LOWER connect.
		} else {
			TR.Warning << "Could not find N-terminal atom in nominally peptide residue." << std::endl;
		}

		/////////////////// PEPTIDE UPPER

		std::string c_term_atm, c_term_connect;

		if ( !c_term.empty() && !leaving.empty() ) {
			// Look for a unique heavy atom which is in the leaving set
			utility::vector1< std::string > leaving_heavy;

			for ( std::string const & atm: c_term ) {
				if ( name_to_element_map.at(atm) != "H" && leaving.count(atm) > 0 ) {
					leaving_heavy.push_back(atm);
					c_term_connect = atm;
				}
			}
			TR.Trace << "Found C-term annotated leaving heavy atoms " << leaving_heavy << std::endl;
			if ( leaving_heavy.size() == 1 ) {
				c_term_connect = leaving_heavy[1];
				TR.Debug << "Picking UPPER location as unique C-term-annotated leaving heavy atom: " << c_term_connect << std::endl;
			} else if ( leaving_heavy.size() > 1 ) {
				// Is the leaving oxygen unique?
				utility::vector1< std::string > leaving_oxy = find_elements( leaving_heavy, "O", name_to_element_map );
				TR.Trace << "Found C-term annotated leaving oxygens " << leaving_oxy << std::endl;
				if ( leaving_oxy.size() == 1 ) {
					c_term_connect = leaving_oxy[1];
					TR.Debug << "Picking UPPER location as unique C-term-annotated leaving oxygen: " << c_term_connect << std::endl;
				}
			}
		}
		if ( c_term_connect.empty() ) {
			// Fallback -- look for an atom named 'OXT'
			if ( name_to_element_map.count("OXT") ) {
				c_term_connect = "OXT";
				TR.Debug << "Picking C-terminal leaving atom based on name: " << c_term_connect << std::endl;
			}
		}

		if ( c_term_connect.empty() ) {
			// Fallback -- look for atom named 'C'
			if ( name_to_element_map.count("C") ) {
				c_term_atm = "C";
				TR.Debug << "Picking C-terminal atom based on name: " << c_term_atm << std::endl;
			}
		}
		if ( c_term_connect.empty() && c_term_atm.empty() ) {
			// Fallback -- look for atom named 'P' -- we may have a phosphonate
			if ( name_to_element_map.count("P") ) {
				c_term_atm = "P";
				TR.Debug << "Picking C-terminal atom based on name: " << c_term_atm << std::endl;
			}
		}

		// Find the c_term atom from the c_term connect atom.
		if ( c_term_atm.empty() && !c_term_connect.empty() ) {
			// Find connected heavy atom
			utility::vector1< std::string > connected_heavy;
			for ( std::string const & atm: get_attached_atoms( c_term_connect, block ) ) {
				if ( name_to_element_map.at(atm) != "H" ) {
					connected_heavy.push_back( atm );
				}
			}
			TR.Trace << "Atoms found connected with the c-term leaving group: " << connected_heavy << std::endl;
			if ( connected_heavy.size() == 1 ) {
				c_term_atm = connected_heavy[1];
			}
		}

		if ( !c_term_atm.empty() ) {
			TR.Debug << "Setting UPPER atom to " << n_term_atm << std::endl;
			molecule.set_upper_atom( c_term_atm );
		} else {
			TR.Warning << "Could not find C-terminal atom in nominally peptide residue." << std::endl;
		}

		/////////////////// ATTACHED H

		//ICOOR_INTERNAL    H   -180.000000   60.849998    1.010000   N     CA  LOWER
		// Since we (may have) deleted ALL hydrogens attached to N, we should add back in an ideal
		// one. We can't add its coordinates yet, though, since only the residue type
		// knows where LOWER is -- The to-ResidueType conversion should update the ICOOR appropriately.

		if ( !n_term_atm.empty() ) {
			utility::vector1< std::string > N_attached_atoms = get_attached_atoms( n_term_atm, block );
			utility::vector1< std::string > N_attached_heavy = find_heavy( N_attached_atoms, name_to_element_map );
			utility::vector1< std::string > N_attached_H = find_elements( N_attached_atoms, "H", name_to_element_map );

			// Add an H if there isn't one already, and we're not in a proline-like situation
			if ( N_attached_H.empty() && N_attached_heavy.size() == 1 ) {
				TR.Debug << "Adding missing H atom to " << n_term_atm << std::endl;
				sdf::MolFileIOAtomOP H_atom( new sdf::MolFileIOAtom());
				auto index = molecule.get_free_index();
				H_atom->index( index );
				H_atom->name( "H" );
				H_atom->element( "H" );
				// faked xyz cordinates
				H_atom->position( core::Vector( 0, 0, 0 ) );
				H_atom->formal_charge( 0 );
				molecule.add_atom( H_atom );

				sdf::MolFileIOBondOP H_bond( new sdf::MolFileIOBond() );
				H_bond->atom1( index );
				H_bond->atom2( molecule.atom( n_term_atm )->index() );
				H_bond->sdf_type( 1 ); //bond order
				molecule.add_bond( H_bond );
			}
		}

	} else if ( is_nucleic_linking ) {

	} else {

	}

}


}
}
}
