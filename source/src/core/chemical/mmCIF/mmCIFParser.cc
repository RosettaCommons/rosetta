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

	//only proceed if the tables for bonds and atoms are present
	if ( !block.has_mmcif_category("_chem_comp_atom") ) {
		TR.Error << "Cannot parse CIF file. No atom block (chem_comp_atom) found for " << block.name << std::endl;
		return molecule;
	}

	/////////////////// Standard Residue-level data
	molecule->name( block.name );

	std::string const * name3 = block.find_value("_chem_comp.three_letter_code"); // Non-owning raw pointer, null if not found.
	if ( name3 ) {
		molecule->name3( as_string(name3) );
	}
	std::string const * name1 = block.find_value("_chem_comp.one_letter_code"); // Non-owning raw pointer, null if not found.
	if ( name1 && name1 ) {
		molecule->name1( as_string(name1) );
	}

	//////////////////// Standard Atom-level data

	// There's another possible issue. to pre-pick about. We absolutely NEED N,
	// because we need to be very specific about adding and deleting atoms.
	// also... residue types without N are very likely to be a poor representative
	// of "L-PEPTIDE LINKING"
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
	int partial_charge = find_gemmi_column(atom_comp,"partial_charge");


	int atom_name_id = find_gemmi_column(atom_comp,"atom_id");
	if ( atom_name_id < 0 ) {
		atom_name_id = find_gemmi_column(atom_comp,"pdbx_component_atom_id");
	}
	if ( atom_name_id < 0 ) {
		TR.Error << "Can't find atom id column (atom_id/pdbx_component_atom_id) in chem_comp_atom table for " << block.name << std::endl;
		return molecule;
	}

	// A map of atom names to their chemical symbols
	std::map< std::string, std::string > name_to_element_map;

	bool N_found = false;
	bool P_found = false;
	bool is_peptide_linking = true;
	bool is_nucleic_linking = true;
	for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
		gemmi::cif::Table::Row row = atom_comp[ii];

		//set atom name
		std::string atom_name = as_string( row[atom_name_id] );
		name_to_element_map[ atom_name ] = as_string( row[type_symbol] );

		if ( atom_name == "N" ) N_found = true;
		if ( atom_name == "P" ) P_found = true;
	}
	if ( ! N_found ) is_peptide_linking = false;
	if ( ! P_found ) is_nucleic_linking = false;

	// It's not nucleic linking if we can't establish that the phosphate P is bonded to
	// O5'. There are some 'reversed' types -- that plausibly we might choose later to
	// read in (AMW TODO) as patched base types -- like T3P that we should USUALLY treat
	// as ligands instead.
	if ( is_nucleic_linking ) {
		bool P_O5P_bond_found = false;
		gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			// Could imagine getting 'all Hs' by finding, instead, the
			// names that match H[number] -- but why not wait, for now.
			if ( ( source == "P" && target == "O5'" ) || ( source == "O5'" && target == "P" ) ) {
				P_O5P_bond_found = true;
			}
		}
		if ( !P_O5P_bond_found ) {
			is_nucleic_linking = false;
		}
	}


	// It's possible (NA8 is one example) that a peptide linking residue will have an
	// atom bonded to OXT that is not named C. Catch this.
	std::string rename_to_C = "C";
	if ( is_peptide_linking ) {
		gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			// If we already have a "C", then don't rename something else to it (even if our 'C' isn't bonded to OXT)
			if ( source == "C" || target == "C" ) {
				rename_to_C = "C";
				break;
			}

			// We only want to rename carbons. Ignore any hydrogens attached, or any non-carbon atoms
			if ( source == "OXT" && name_to_element_map[ target ] == "C" ) {
				rename_to_C = target;
			} else if ( target == "OXT" && name_to_element_map[ source ] == "C" ) {
				rename_to_C = source;
			}
		}
	}


	// Get the chem_comp table first, because this will help us
	// look out for extraneous atoms common in CIF entries -- extra nitrogen H
	// and OH terminus on C
	gemmi::cif::Table chem_comp_type = block.find( "_chem_comp.", {"type"} );
	if ( chem_comp_type.size() > 0 ) {
		std::string type = as_string(chem_comp_type[0][0]);
		if ( type == "L-PEPTIDE LINKING" && is_peptide_linking ) {
			TR.Debug << "Found L-peptide RT" << std::endl;// named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "PROTEIN POLYMER L_AA" );
			is_peptide_linking = true;
			is_nucleic_linking = false;
		} else if ( type == "D-PEPTIDE LINKING" && is_peptide_linking ) {
			TR.Debug << "Found D-peptide RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "PROTEIN POLYMER D_AA" );
			is_peptide_linking = true;
			is_nucleic_linking = false;
		} else if ( type == "RNA LINKING" && is_nucleic_linking ) {
			TR.Debug << "Found D-RNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "RNA POLYMER D_RNA" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else if ( type == "L-RNA LINKING" && is_nucleic_linking ) {
			TR.Debug << "Found L-RNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "RNA POLYMER L_RNA" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		}  else if ( type == "DNA LINKING" && is_nucleic_linking ) {
			TR.Debug << "Found D-DNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "DNA POLYMER" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else if ( type == "L-DNA LINKING" && is_nucleic_linking ) {
			TR.Debug << "Found L-DNA RT" << std::endl;//named " << molecule->name() << std::endl;
			molecule->add_str_str_data( "Rosetta Properties", "DNA POLYMER" );
			is_peptide_linking = false;
			is_nucleic_linking = true;
		} else {
			is_peptide_linking = false;
			is_nucleic_linking = false;
		}
	}



	// Sometimes OP3/O3P is used for non-term-deletable phosphate oxygens. (THX)
	//bool interesting_upper_behavior = false;

	utility::vector1< std::string > O3P_connected;

	// Before we actually LOOK at the atoms (or bonds) for real, we need to know
	// what atoms are bonded to the polymeric termini or to to-be-deleted
	// atoms -- so we can ignore them too as appropriate
	utility::vector1< std::string > possible_atoms_to_skip;
	if ( is_peptide_linking ) {
		gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			if ( source == rename_to_C ) source = "C";
			if ( target == rename_to_C ) target = "C";

			// Could imagine getting 'all Hs' by finding, instead, the
			// names that match H[number] -- but why not wait, for now.
			if ( source == "OXT" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to OXT " << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "OXT" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to OXT " << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( source == "N" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to N " << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "N" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to N " << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
		}
	} else if ( is_nucleic_linking ) {
		gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			if ( source == "P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  P" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  P" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}

			if ( source == "O3'" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  O3'" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "O3'" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  O3'" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}

			if ( source == "OP3" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  OP3" << std::endl;
				possible_atoms_to_skip.push_back( target );
				O3P_connected.push_back( target );
			}
			if ( target == "OP3" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  OP3" << std::endl;
				possible_atoms_to_skip.push_back( source );
				O3P_connected.push_back( source );
			}
			if ( source == "O3P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  O3P" << std::endl;
				possible_atoms_to_skip.push_back( target );
				O3P_connected.push_back( target );
			}
			if ( target == "O3P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  O3P" << std::endl;
				possible_atoms_to_skip.push_back( source );
				O3P_connected.push_back( source );
			}
			if ( source == "OP2" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  OP2" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "OP2" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  OP2" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( target == "O2P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  O2P" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( source == "O2P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  O2P" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( source == "OP1" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  OP1" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "OP1" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  OP1" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( source == "O1P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to  O1P" << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "O1P" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to  O1P" << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
		}
	} else {
		gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2"} );
		for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
			std::string source = as_string( bond_comp[ii][0] ); //atom 1
			std::string target = as_string( bond_comp[ii][1] ); //atom 2

			if ( source == "O2A" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to O2A " << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "O2A" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to O2A " << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( source == "O2B" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to O2B " << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "O2B" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to O2B " << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
			if ( source == "O3B" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << target << " due to its bond to O3B " << std::endl;
				possible_atoms_to_skip.push_back( target );
			}
			if ( target == "O3B" ) {
				TR.Trace << "It may be appropriate to skip the maybe-hydrogen " << source << " due to its bond to O3B " << std::endl;
				possible_atoms_to_skip.push_back( source );
			}
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

	// Loop over atom block and check to see if any heavyatoms are bound to OP3/O3P
	bool interesting_pendant = false;
	for ( Size ii = 0; ii < atom_comp.size(); ++ii ) {
		std::string atom_name = as_string( atom_comp[ ii ][ atom_name_id ] );
		std::string element = as_string( atom_comp[ii][ type_symbol ] );
		if ( O3P_connected.contains( atom_name ) ) {
			if ( element != "H" && atom_name != "P" ) {
				TR.Trace << "There is an OP3-bonded heavyatom: " << atom_name << std::endl;
				interesting_pendant = true;
				break;
			}
		}
	}

	utility::vector1< std::string > actual_atoms_to_skip;
	//start atom block
	Size index = 1;
	TR.Trace << "possible_atoms_to_skip: " << possible_atoms_to_skip << std::endl;
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
		if ( is_peptide_linking && atom_name == "OXT" ) continue;
		if ( is_nucleic_linking && !interesting_pendant && atom_name == "OP3" ) continue;
		if ( is_nucleic_linking && !interesting_pendant && atom_name == "O3P" ) continue;

		if ( atom_name == rename_to_C ) atom_name = "C";

		atom->name( atom_name );

		//set map to index
		atom_name_to_id[ atom_name ] = index;
		//set element name
		atom->element( as_string( atom_comp[ii][type_symbol] ) );

		TR.Trace << "Type symbol for atom  " << atom_name << " is " <<  as_string(atom_comp[ii][type_symbol]) << std::endl;
		if ( possible_atoms_to_skip.contains( atom_name ) && atom->element() == "H" ) {
			actual_atoms_to_skip.push_back( atom_name );
			continue;
		}

		TR.Trace << "Keeping atom entry " << atom_name << std::endl;

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
		if ( partial_charge >= 0 ) {
			atom->partial_charge( as_number( atom_comp[ii][partial_charge], 0 ) ); // Default zero if present and null
		}

		molecule->add_atom( atom );
		// only increment if we actually get here.
		index++;
	}

	TR.Trace << "actual_atoms_to_skip: " << actual_atoms_to_skip << std::endl;

	// pdb_BVP messes with this logic because it names H3' "H1" and HO3' H3'.
	// So we need to place an additional requirement for H skipping.
	// This is actually a very hard problem.

	gemmi::cif::Table bond_comp = block.find( "_chem_comp_bond.", {"atom_id_1","atom_id_2","value_order"} );
	for ( core::Size ii(0); ii < bond_comp.size(); ++ii ) {
		std::string source = as_string( bond_comp[ii][0] ); //atom 1
		std::string target = as_string( bond_comp[ii][1] ); //atom 2
		core::Size bond_type( bond_string_to_sdf_size_[ as_string(bond_comp[ii][2]) ] );

		sdf::MolFileIOBondOP bond( new sdf::MolFileIOBond() );

		if ( source == rename_to_C ) source = "C";
		if ( target == rename_to_C ) target = "C";

		TR.Trace << "Examining bond entry " << source << " " << target << std::endl;

		if ( is_peptide_linking && source == "OXT" ) continue;
		if ( is_peptide_linking && target == "OXT" ) continue;
		if ( is_nucleic_linking && !interesting_pendant && ( source == "OP3" || source == "O3P" ) ) continue;
		if ( is_nucleic_linking && !interesting_pendant && ( target == "OP3" || target == "O3P" ) ) continue;
		if ( actual_atoms_to_skip.contains( source ) ) continue;
		if ( actual_atoms_to_skip.contains( target ) ) continue;

		TR.Trace << "Keeping bond entry " << source << " " << target << std::endl;

		bond->atom1( atom_name_to_id[ source ] );
		bond->atom2( atom_name_to_id[ target ] );

		bond->sdf_type( bond_type); //bond order

		molecule->add_bond( bond);
	}
	if ( bond_comp.size() == 0 && atom_comp.size() > 1 ) {
		TR.Error << "Cannot parse CIF file. No bond block (chem_comp_bond) found for multi-atom entry " << block.name << std::endl;
	} // else one atom entry without bond block

	// Note, for polymer generalization we may want to find the block here that
	// describes polymeric connections, polymer type, et cetera -- and use it
	// to add UPPER and LOWER entries. Otherwise, we may need/want to do so
	// using molfile comments records... let's see, when we look at the parser,
	// how possible that will be.

	//ICOOR_INTERNAL    H   -180.000000   60.849998    1.010000   N     CA  LOWER
	// Since we deleted ALL hydrogens attached to N, we should add back in an ideal
	// one. We can't add its coordinates yet, though, since only the residue type
	// knows where LOWER is. Hey, actually, LOWER needs to be assigned in the
	// list too!
	if ( is_peptide_linking ) {
		sdf::MolFileIOAtomOP H_atom( new sdf::MolFileIOAtom());
		//atom id is whatever the number we are on +1, because we dont do 0 based index
		H_atom->index( index );
		H_atom->name( "H" );
		atom_name_to_id[ "H" ] = index;
		H_atom->element( "H" );
		// faked xyz cordinates
		H_atom->position( core::Vector( 0, 0, 0 ) );
		H_atom->formal_charge( 0 );
		molecule->add_atom( H_atom );

		sdf::MolFileIOBondOP H_bond( new sdf::MolFileIOBond() );
		H_bond->atom1( atom_name_to_id[ "H" ] );
		H_bond->atom2( atom_name_to_id[ "N" ] );
		H_bond->sdf_type( 1 ); //bond order
		molecule->add_bond( H_bond );
	}

	return molecule;
}


}
}
}
