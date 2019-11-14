// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/mol_writer.cc
/// @author Sam DeLuca
/// @details this class outputs a residue in the form of a V3000 molfile, for details of the spec, see http://www.symyx.com/downloads/public/ctfile/ctfile.pdf

#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/sdf/mol_util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <utility/string_util.hh>
#include <iomanip>
#include <iostream>
#include <algorithm>

// Boost Headers
#include <boost/format.hpp>

#include <utility/vector1.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer tr( "core.chemical.sdf.mol_writer" );


using std::setfill;
using std::setw;
using std::dec;
MolWriter::MolWriter() : line_header_("M  V30 "),ctab_mode_(V2000)
{}

MolWriter::MolWriter(std::string const & ctab_mode) : line_header_("M  V30 ")
{
	if ( ctab_mode == "V2000" ) {
		ctab_mode_ = V2000;
	} else {
		ctab_mode_ = V3000;
	}
}

void
MolWriter::output_residue(std::ostream & output_stream, core::conformation::ResidueCOP residue) {
	output_residue( output_stream, *residue);
}

void
MolWriter::output_residue(std::ostream & output_stream, core::chemical::ResidueTypeCOP residue_type) {
	output_residue( output_stream, *residue_type);
}

void
MolWriter::output_residue(std::ostream & output_stream, core::chemical::MutableResidueTypeCOP residue_type) {
	output_residue( output_stream, *residue_type);
}

void MolWriter::output_residue(std::ostream & output_stream, core::conformation::Residue const & residue)
{
	std::map< std::string, core::Vector > coords;
	for ( core::Size ii(1); ii <= residue.natoms(); ++ ii ) {
		coords[ residue.atom_name(ii) ] = residue.xyz(ii);
	}

	core::chemical::MutableResidueType mut_type( residue.type() );
	output_residue_impl(output_stream, mut_type, coords );
}

void MolWriter::output_residue(std::ostream & output_stream, core::chemical::ResidueType const & residue_type)
{
	core::chemical::MutableResidueType mut_type( residue_type );
	output_residue_impl(output_stream, mut_type );
}

void
MolWriter::output_residue(std::ostream & output_stream, core::chemical::MutableResidueType const & residue_type) {
	output_residue_impl( output_stream, residue_type );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MolWriter::output_residue_impl(std::ostream & output_stream, core::chemical::MutableResidueType const & residue, std::map< std::string, core::Vector > const & coords )
{
	std::list<std::string> prepared_lines;
	std::list<std::string> metadata = this->compose_metadata(residue);
	std::list<std::string> ctab = this->compose_ctab(residue,coords);
	std::list<std::string> job_data = this->compose_job_info();

	prepared_lines.insert(prepared_lines.end(),metadata.begin(),metadata.end());
	prepared_lines.insert(prepared_lines.end(),ctab.begin(),ctab.end());
	if ( ! basic::options::option[ basic::options::OptionKeys::out::file::no_extra_sdf_data ] ) {
		// These entries should be kept up-to-date with those processed in for MolFileIOData.cc
		std::list<std::string> naming = this->compose_naming(residue);
		std::list<std::string> atomnames = this->compose_atomnames(residue);
		std::list<std::string> typeinfo = this->compose_typeinfo(residue);
		std::list<std::string> nbr_atom = this->compose_nbr_atom(residue);
		std::list<std::string> properties = this->compose_rosetta_properties(residue);

		prepared_lines.insert(prepared_lines.end(),naming.begin(),naming.end());
		prepared_lines.insert(prepared_lines.end(),atomnames.begin(),atomnames.end());
		prepared_lines.insert(prepared_lines.end(),typeinfo.begin(),typeinfo.end());
		prepared_lines.insert(prepared_lines.end(),nbr_atom.begin(),nbr_atom.end());
		prepared_lines.insert(prepared_lines.end(),properties.begin(),properties.end());
	}
	prepared_lines.insert(prepared_lines.end(),job_data.begin(),job_data.end());
	if ( ctab_mode_ == V2000 ) {
		prepared_lines.emplace_back("$$$$\n");
	}


	for ( std::string const & line : prepared_lines ) {
		output_stream << line;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::list<std::string> MolWriter::compose_metadata(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;

	std::string const name_line = residue.name()+"\n";
	//for specification see page 12 and 35 of the sybl molfile specification document
	std::string const info_line = "  Rosetta           3D                              \n";
	std::string const name3_line = residue.name3()+"\n";
	std::string counts_line;
	if ( ctab_mode_ == V2000 ) {

		counts_line = str(boost::format("%|1|%|2|%|3|%|4|%|5|%|6|%|7|%|8|%|9|%|10|%|11|%|12|\n") %
			boost::io::group(setfill(' '),dec,setw(3),residue.natoms()) % //atom count
			boost::io::group(setfill(' '),dec,setw(3),residue.nbonds()) % //bond count
			boost::io::group(setfill(' '),dec,setw(3),0) % //number of atom lists
			boost::io::group(setfill(' '),dec,setw(3),0) % //obselete field
			boost::io::group(setfill(' '),dec,setw(3),0) %  // chiral flag
			boost::io::group(setfill(' '),dec,setw(3),0) % // number of stexts
			boost::io::group(setfill(' '),dec,setw(3),0) % // obselete field
			boost::io::group(setfill(' '),dec,setw(3),0) %  // obselete field
			boost::io::group(setfill(' '),dec,setw(3),0) %  // obselete field
			boost::io::group(setfill(' '),dec,setw(3),0) %  // obselete field
			boost::io::group(setfill(' '),dec,setw(3),999) %  //extra properties
			boost::io::group(setfill(' '),setw(6),"V2000")); //version string;


	} else {
		counts_line = "  0  0  0     0  0            999 V3000\n";
	}


	lines.push_back(name_line);
	lines.push_back(info_line);
	lines.push_back(name3_line);
	lines.push_back(counts_line);
	return lines;
}

std::list<std::string> MolWriter::compose_ctab(core::chemical::MutableResidueType const & residue, std::map< std::string, core::Vector > const & coords )
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN CTAB\n";
	std::string end_header = line_header_+"END CTAB\n";

	core::Size n_atoms = residue.natoms();
	core::Size n_bonds = residue.nbonds();

	std::string counts = line_header_+"COUNTS "+ utility::to_string<core::Size>(n_atoms)+" "+
		utility::to_string<core::Size>(n_bonds)+" " + "0" + " " + "0" + " "+ "0"+"\n";

	//compose bonds and atoms, append all this to the ctab

	std::list<std::string> atom_lines = this->compose_atoms(residue, coords);
	std::list<std::string> bond_lines = this->compose_bonds(residue);
	std::list<std::string> prop_lines = this->compose_properties(residue);

	if ( ctab_mode_ == V3000 ) {
		lines.push_back(begin_header);
		lines.push_back(counts);
	}

	lines.insert(lines.end(),atom_lines.begin(),atom_lines.end());
	lines.insert(lines.end(),bond_lines.begin(),bond_lines.end());
	lines.insert(lines.end(),prop_lines.begin(),prop_lines.end());

	if ( ctab_mode_ == V3000 ) {
		lines.push_back(end_header);
	}
	lines.emplace_back("M  END\n");

	return lines;
}

std::list<std::string> MolWriter::compose_atoms(core::chemical::MutableResidueType const & residue, std::map< std::string, core::Vector > const & coords )
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN ATOM"+"\n";
	std::string end_header = line_header_+"END ATOM"+"\n";

	if ( ctab_mode_ == V3000 ) {
		lines.push_back(begin_header);
	}
	for ( VD vd: residue.all_atoms() ) {
		core::Vector xyz_coords( residue.atom(vd).ideal_xyz() );
		if ( coords.count( residue.atom_name(vd) ) ) {
			xyz_coords = coords.at( residue.atom_name(vd) );
		}
		core::chemical::AtomType const & atom_type = residue.atom_type(vd);
		std::string element = atom_type.element();
		//Rosetta stores elements as allcaps for whatever reason, this will turn CL -> Cl
		if ( element.size() == 2 ) {
			element[1] = tolower(element[1]);
		}
		core::Real charge = residue.atom(vd).charge();
		std::string atom_string;
		core::Size hydrogen_count = 0;
		core::Size all_bound_count = 0;
		if ( ! residue.atom(vd).is_hydrogen() ) {
			hydrogen_count = residue.bonded_hydrogens(vd).size();
			all_bound_count = residue.bonded_neighbors(vd).size();
		}

		if ( ctab_mode_ == V2000 ) {
			atom_string = str(boost::format("%|1|%|2|%|3| %|4|%|5|%|6|%|7|%|8|%|9|%|10|%|11|%|12|%|13|%|14|%|15|%|16|\n") %
				boost::io::group(setfill(' '),dec,setw(10),std::fixed,std::setprecision(4),xyz_coords.x()) % // x coord
				boost::io::group(setfill(' '),dec,setw(10),std::fixed,std::setprecision(4),xyz_coords.y()) % // y coord
				boost::io::group(setfill(' '),dec,setw(10),std::fixed,std::setprecision(4),xyz_coords.z()) % // z coord
				boost::io::group(setfill(' '),setw(3),std::left,element) % // atom symbol
				boost::io::group(setfill(' '),dec,setw(2),0) % // mass difference
				boost::io::group(setfill(' '),dec,setw(3),0) % // charge
				boost::io::group(setfill(' '),dec,setw(3),0) % // atom stereo parity
				boost::io::group(setfill(' '),dec,setw(3),hydrogen_count+1) % // hydrogen count + 1
				boost::io::group(setfill(' '),dec,setw(3),0) % //stereo care box
				boost::io::group(setfill(' '),dec,setw(3),all_bound_count) % //valance
				boost::io::group(setfill(' '),dec,setw(3),0) % //h0 designator
				boost::io::group(setfill(' '),dec,setw(3),0) % //unused field
				boost::io::group(setfill(' '),dec,setw(3),0) % //unused field
				boost::io::group(setfill(' '),dec,setw(3),0) % // atom atom mapping
				boost::io::group(setfill(' '),dec,setw(3),0) % // inversion/retention
				boost::io::group(setfill(' '),dec,setw(3),0)); //exact change flag
		} else {
			atom_string = line_header_ + " " + std::to_string(residue.atom_index(vd))+" "+
				element+" "+ utility::to_string<core::Real>(xyz_coords.x())+ " "+
				utility::to_string<core::Real>(xyz_coords.y())+ " " +
				utility::to_string<core::Real>(xyz_coords.z())+ " " +
				"0"+" "+"CHG="+utility::to_string<core::Real>(charge)+"\n";
		}

		lines.push_back(atom_string);
	}
	if ( ctab_mode_ == V3000 ) {
		lines.push_back(end_header);
	}
	return lines;
}

std::list<std::string> MolWriter::compose_bonds(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN BOND"+"\n";
	std::string end_header = line_header_+"END BOND"+"\n";

	if ( ctab_mode_ == V3000 ) {
		lines.push_back(begin_header);
	}
	//bonds are non-directional so we want to make sure they only get counted once.
	//(bond from atom 1 to atom 5 is the same as the one from 5 to 1)
	//stick all the bonds in a set and then iterate through the set to output
	//the smaller atom index is always designated lower in the struct


	utility::vector1<BondData> bond_data_set;
	for ( VD vd: residue.all_atoms() ) {
		utility::vector1< VD > const bonded_neighbors = residue.bonded_neighbors(vd);

		for ( core::Size neighbor_index = 1; neighbor_index <= bonded_neighbors.size(); ++neighbor_index ) {
			VD neighbor = bonded_neighbors[neighbor_index];
			core::Size type = residue.bond( vd, neighbor ).bond_name();
			BondData bond(residue.atom_index(vd),residue.atom_index(neighbor),type);

			if ( std::find(bond_data_set.begin(),bond_data_set.end(),bond)!= bond_data_set.end() ) {
				continue;
			}
			//std::cout << bond.lower << " " <<  bond.upper <<" " <<bond.bondType << std::endl;
			bond_data_set.push_back(bond);
		}
	}

	core::Size bond_index = 1;
	std::string bond_line;
	for ( BondData const & current_bond : bond_data_set ) {
		if ( ctab_mode_ == V2000 ) {
			bond_line = str(boost::format("%|1|%|2|%|3|%|4|%|5|%|6|%|7|\n") %
				boost::io::group(setfill(' '),dec,setw(3),current_bond.lower) % // lower bond
				boost::io::group(setfill(' '),dec,setw(3),current_bond.upper) % // upper bond
				boost::io::group(setfill(' '),dec,setw(3),current_bond.bondType) % // bond type
				boost::io::group(setfill(' '),dec,setw(3),0) % // bond stereo
				boost::io::group(setfill(' '),dec,setw(3),0) % // unused
				boost::io::group(setfill(' '),dec,setw(3),0) % // bond topology
				boost::io::group(setfill(' '),dec,setw(3),0)); //reaction center status
		} else {
			bond_line = line_header_ + utility::to_string<core::Size>(bond_index) + " " +
				utility::to_string<core::Size>(current_bond.bondType) + " " +
				utility::to_string<core::Size>(current_bond.lower) + " " +
				utility::to_string<core::Size>(current_bond.upper) + " " + "\n";
		}

		lines.push_back(bond_line);
		++bond_index;
	}

	if ( ctab_mode_==V3000 ) {
		lines.push_back(end_header);
	}
	return lines;
}

std::list<std::string> MolWriter::compose_properties(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;
	if ( ctab_mode_ == V3000 ) {
		// Properties block is V2000 specific.
		return lines;
	}

	////////////
	// Charges:
	utility::vector1< VD > charged_atoms;
	for ( VD atm: residue.all_atoms() ) {
		if ( residue.atom( atm ).formal_charge() != 0 ) {
			charged_atoms.push_back( atm );
		}
	}
	// V2000 CHG lines can only have 8 atoms max per line
	for ( core::Size b(1), e(9); b <= charged_atoms.size(); b += 8, e+=8 ) {
		core::Size const end( std::min<core::Size>( e, charged_atoms.size()+1 ) );
		core::Size const nentries = end - b;
		runtime_assert( nentries >=1 && nentries <= 8 );
		std::string line( boost::str( boost::format("M  CHG%3d") % nentries ) );
		for ( core::Size n(b); n < end; ++n ) {
			line.append( boost::str( boost::format(" %3d %3d")
				% residue.atom_index( charged_atoms[n] )
				% residue.atom( charged_atoms[n] ).formal_charge() ));
		}
		line.append( "\n" );
		lines.push_back( line );
	}

	////////////
	// Other Properties?

	return lines;
}

std::list<std::string> MolWriter::compose_atomnames(core::chemical::MutableResidueType const & residue_type)
{
	std::list<std::string> lines;

	std::string header = "> <Atom Names>\n";
	std::string type_data = "";
	utility::vector1< VD > const & all_atoms( residue_type.all_atoms() ); // Hopefully this is consistent enough for our purposes.
	for ( core::Size index =1; index <= all_atoms.size(); ++index ) {
		std::string const & atom_name = residue_type.atom_name( all_atoms[index] );
		std::string data_string = "("+utility::to_string<core::Size>(index)+","+atom_name+") ";
		type_data.append(data_string);
	}
	type_data.append("\n");

	lines.push_back(header);
	lines.push_back(type_data);
	lines.push_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_typeinfo(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta AtomTypes>\n";
	std::string type_data = "";
	for ( VD atm: residue.all_atoms() ) {
		std::string const & atom_type_name = residue.atom_type(atm).name();
		std::string data_string = "("+utility::to_string<core::Size>(residue.atom_index(atm))+","+atom_type_name+") ";
		type_data.append(data_string);
	}
	type_data.append("\n");

	lines.push_back(header);
	lines.push_back(type_data);
	lines.emplace_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_nbr_atom(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta nbr_atom>\n";
	std::string nbr_atom = std::to_string( residue.atom_index( residue.nbr_vertex() ) ) + "\n";

	lines.push_back(header);
	lines.push_back(nbr_atom);
	lines.emplace_back("\n");

	std::string header2 = "> <Rosetta nbr_radius>\n";
	std::string nbr_radius = utility::to_string<core::Real>(residue.nbr_radius()) + "\n";

	lines.push_back(header2);
	lines.push_back(nbr_radius);
	lines.emplace_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_naming(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta Name>\n";
	std::string name = residue.name() + "\n";

	lines.push_back(header);
	lines.push_back(name);
	lines.emplace_back("\n");

	if ( residue.name3() != residue.name().substr(0,3) || residue.name1() != 'Z' ) {
		std::string header2 = "> <Rosetta IO_string>\n";
		std::string io_string = residue.name3() + " " + residue.name1() + "\n";

		lines.push_back(header2);
		lines.push_back(io_string);
		lines.emplace_back("\n");
	}
	if ( residue.interchangeability_group() != residue.name3() ) {
		std::string header3 = "> <Rosetta Interchangeability Group>\n";
		std::string group = residue.interchangeability_group() + "\n";

		lines.push_back(header3);
		lines.push_back(group);
		lines.emplace_back("\n");
	}
	if ( residue.aa() != core::chemical::aa_unk ) {
		std::string header4 = "> <Rosetta AA>\n";
		std::string aa = utility::to_string<core::chemical::AA>( residue.aa() ) + "\n";

		lines.push_back(header4);
		lines.push_back(aa);
		lines.emplace_back("\n");
	}

	return lines;
}

std::list<std::string> MolWriter::compose_rosetta_properties(core::chemical::MutableResidueType const & residue)
{
	std::list<std::string> lines;

	utility::vector1< std::string > const & properties( residue.properties().get_list_of_properties() );

	if ( properties.size() ) {
		std::string header = "> <Rosetta Properties>\n";
		lines.push_back(header);

		for ( core::Size ii(1); ii <= properties.size(); ++ii ) {
			std::string property_line = properties[ ii ] + "\n";
			lines.push_back( property_line );
		}

		lines.emplace_back("\n" );
	}

	return lines;
}

std::list<std::string> MolWriter::compose_job_info()
{
	std::list<std::string> lines;
	for ( std::map<std::string,std::string>::const_iterator data_it = job_data_.begin(); data_it != job_data_.end(); ++data_it ) {
		std::string header_name(data_it->first);
		std::string data(data_it->second+"\n");
		std::string header("> <"+header_name+">\n");
		lines.push_back(header);
		lines.push_back(data);
		lines.emplace_back("\n");
	}
	return lines;
}


}
}
}
