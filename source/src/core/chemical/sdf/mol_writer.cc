// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/mol_writer.cc
/// @author Sam DeLuca
/// @details this class outputs a residue in the form of a V3000 molfile, for details of the spec, see http://www.symyx.com/downloads/public/ctfile/ctfile.pdf

#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/sdf/mol_util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <set>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

namespace core {
namespace chemical {
namespace sdf {

static thread_local basic::Tracer tr( "core.chemical.sdf.mol_writer" );


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

void MolWriter::output_residue(utility::io::ozstream & output_stream, core::conformation::ResidueCOP residue)
{
	std::list<std::string> prepared_lines;
	std::list<std::string> metadata = this->compose_metadata(residue);
	std::list<std::string> ctab = this->compose_ctab(residue);
	std::list<std::string> job_data = this->compose_job_info();

	prepared_lines.insert(prepared_lines.end(),metadata.begin(),metadata.end());
	prepared_lines.insert(prepared_lines.end(),ctab.begin(),ctab.end());
	if ( ! basic::options::option[ basic::options::OptionKeys::out::file::no_extra_sdf_data ] ) {
		// These entries should be kept up-to-date with those processed in for MolFileIOData.cc
		std::list<std::string> naming = this->compose_naming(residue);
		std::list<std::string> typeinfo = this->compose_typeinfo(residue);
		std::list<std::string> nbr_atom = this->compose_nbr_atom(residue);
		std::list<std::string> properties = this->compose_rosetta_properties(residue);

		prepared_lines.insert(prepared_lines.end(),naming.begin(),naming.end());
		prepared_lines.insert(prepared_lines.end(),typeinfo.begin(),typeinfo.end());
		prepared_lines.insert(prepared_lines.end(),nbr_atom.begin(),nbr_atom.end());
		prepared_lines.insert(prepared_lines.end(),properties.begin(),properties.end());
	}
	prepared_lines.insert(prepared_lines.end(),job_data.begin(),job_data.end());
	if ( ctab_mode_ == V2000 ) {
		prepared_lines.push_back("$$$$\n");
	}


	BOOST_FOREACH ( std::string line, prepared_lines ) {
		output_stream << line;
	}

}


void MolWriter::output_residue(utility::io::ozstream & output_stream, core::chemical::ResidueTypeOP residue_type)
{
	using namespace core::conformation;
	ResidueCOP residue_ptr( ResidueOP( new Residue(*residue_type,false) ) );
	//core::conformation::ResidueCOP residue_ptr= &residue;
	//std::cout <<residue_ptr->name3() <<std::endl;
	output_residue(output_stream,residue_ptr);
}


void MolWriter::output_residue(std::string file_name,core::conformation::ResidueCOP residue)
{
	utility::io::ozstream outfile;

	outfile.open(file_name.c_str(),std::ios::out | std::ios::binary);
	if ( !outfile ) {
		throw utility::excn::EXCN_FileNotFound("Cannot open file"+file_name);
	}
	output_residue(outfile,residue);
	outfile.close();
}


void MolWriter::output_residue(std::string file_name, core::chemical::ResidueTypeOP residue_type)
{
	utility::io::ozstream outfile;
	outfile.open(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !outfile ) {
		throw utility::excn::EXCN_FileNotFound("Cannot open file"+file_name);
	}
	output_residue(outfile,residue_type);
	outfile.close();
}

std::list<std::string> MolWriter::compose_metadata(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;

	std::string const name_line = residue->name()+"\n";
	//for specification see page 12 and 35 of the sybl molfile specification document
	std::string const info_line = "  Rosetta           3D                              \n";
	std::string const name3_line = residue->name3()+"\n";
	std::string counts_line;
	if ( ctab_mode_ == V2000 ) {

		counts_line = str(boost::format("%|1|%|2|%|3|%|4|%|5|%|6|%|7|%|8|%|9|%|10|%|11|%|12|\n") %
			boost::io::group(setfill(' '),dec,setw(3),residue->natoms()) % //atom count
			boost::io::group(setfill(' '),dec,setw(3),residue->type().nbonds()) % //bond count
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

std::list<std::string> MolWriter::compose_ctab(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN CTAB\n";
	std::string end_header = line_header_+"END CTAB\n";

	core::Size n_atoms = residue->natoms();
	core::Size n_bonds = residue->type().nbonds();

	std::string counts = line_header_+"COUNTS "+ utility::to_string<core::Size>(n_atoms)+" "+
		utility::to_string<core::Size>(n_bonds)+" " + "0" + " " + "0" + " "+ "0"+"\n";

	//compose bonds and atoms, append all this to the ctab

	std::list<std::string> atom_lines = this->compose_atoms(residue);
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
	lines.push_back("M  END\n");

	return lines;
}

std::list<std::string> MolWriter::compose_atoms(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN ATOM"+"\n";
	std::string end_header = line_header_+"END ATOM"+"\n";

	if ( ctab_mode_ == V3000 ) {
		lines.push_back(begin_header);
	}
	for ( core::Size index = 1; index <= residue->natoms(); ++index ) {
		core::Vector xyz_coords(residue->xyz(index));
		core::chemical::AtomType const atom_type = residue->atom_type(index);
		core::chemical::ResidueType const & residue_type = residue->type();
		std::string element = atom_type.element();
		//Rosetta stores elements as allcaps for whatever reason, this will turn CL -> Cl
		if ( element.size() == 2 ) {
			element[1] = tolower(element[1]);
		}
		core::Real charge = residue_type.atom(index).charge();
		std::string atom_string;
		core::Size hydrogen_count = 0;
		core::Size heavy_bond_count = 0;
		if ( index <= residue_type.nheavyatoms() ) {
			hydrogen_count = residue_type.number_bonded_hydrogens(index);
			heavy_bond_count = residue_type.number_bonded_heavyatoms(index);
		} else {
			hydrogen_count = 0;
			heavy_bond_count = 0;
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
				boost::io::group(setfill(' '),dec,setw(3),hydrogen_count+heavy_bond_count) % //valance
				boost::io::group(setfill(' '),dec,setw(3),0) % //h0 designator
				boost::io::group(setfill(' '),dec,setw(3),0) % //unused field
				boost::io::group(setfill(' '),dec,setw(3),0) % //unused field
				boost::io::group(setfill(' '),dec,setw(3),0) % // atom atom mapping
				boost::io::group(setfill(' '),dec,setw(3),0) % // inversion/retention
				boost::io::group(setfill(' '),dec,setw(3),0)); //exact change flag
		} else {
			atom_string = line_header_ + " " + utility::to_string<core::Size>(index)+" "+
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

std::list<std::string> MolWriter::compose_bonds(core::conformation::ResidueCOP residue)
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
	for ( core::Size index = 1; index <= residue->natoms(); ++index ) {
		core::chemical::AtomIndices const bonded_neighbors = residue->bonded_neighbor(index);
		utility::vector1<core::chemical::BondName> const bonded_neighbor_types = residue->type().bonded_neighbor_types(index);
		debug_assert(bonded_neighbors.size()== bonded_neighbor_types.size());

		for ( core::Size neighbor_index = 1; neighbor_index <= bonded_neighbors.size(); ++neighbor_index ) {
			core::Size type = bonded_neighbor_types[neighbor_index];
			core::Size neighbor = bonded_neighbors[neighbor_index];
			BondData bond(index,neighbor,type);

			if ( std::find(bond_data_set.begin(),bond_data_set.end(),bond)!= bond_data_set.end() ) {
				continue;
			}
			//std::cout << bond.lower << " " <<  bond.upper <<" " <<bond.bondType << std::endl;
			bond_data_set.push_back(bond);
		}
	}

	core::Size bond_index = 1;
	std::string bond_line;
	BOOST_FOREACH ( BondData current_bond, bond_data_set ) {
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

std::list<std::string> MolWriter::compose_properties(core::conformation::ResidueCOP residue)
{
	runtime_assert( residue != 0 );
	std::list<std::string> lines;
	if ( ctab_mode_ == V3000 ) {
		// Properties block is V2000 specific.
		return lines;
	}

	////////////
	// Charges:
	utility::vector1< core::Size > charged_atoms;
	for ( core::Size ii(1); ii <= residue->natoms(); ++ii ) {
		if ( residue->type().atom(ii).formal_charge() != 0 ) {
			charged_atoms.push_back( ii );
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
				% charged_atoms[n]
				% residue->type().atom(charged_atoms[n]).formal_charge() ));
		}
		line.append( "\n" );
		lines.push_back( line );
	}

	////////////
	// Other Properties?

	return lines;
}

std::list<std::string> MolWriter::compose_typeinfo(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta AtomTypes>\n";
	std::string type_data = "";
	core::chemical::ResidueType const & residue_type = residue->type();
	for ( core::Size index =1; index <= residue->natoms(); ++index ) {
		std::string atom_type_name = residue_type.atom_type(index).name();
		std::string data_string = "("+utility::to_string<core::Size>(index)+","+atom_type_name+") ";
		type_data.append(data_string);
	}
	type_data.append("\n");

	lines.push_back(header);
	lines.push_back(type_data);
	lines.push_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_nbr_atom(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta nbr_atom>\n";
	std::string nbr_atom = utility::to_string<core::Size>(residue->nbr_atom()) + "\n";

	lines.push_back(header);
	lines.push_back(nbr_atom);
	lines.push_back("\n");

	std::string header2 = "> <Rosetta nbr_radius>\n";
	std::string nbr_radius = utility::to_string<core::Real>(residue->nbr_radius()) + "\n";

	lines.push_back(header2);
	lines.push_back(nbr_radius);
	lines.push_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_naming(core::conformation::ResidueCOP residue)
{
	core::chemical::ResidueType const & restype( residue->type() );
	std::list<std::string> lines;

	std::string header = "> <Rosetta Name>\n";
	std::string name = residue->name() + "\n";

	lines.push_back(header);
	lines.push_back(name);
	lines.push_back("\n");

	if ( residue->name3() != residue->name().substr(0,3) || residue->name1() != 'Z' ) {
		std::string header2 = "> <Rosetta IO_string>\n";
		std::string io_string = residue->name3() + ' ' + residue->name1() + "\n";

		lines.push_back(header2);
		lines.push_back(io_string);
		lines.push_back("\n");
	}
	if ( restype.interchangeability_group() != restype.name3() ) {
		std::string header3 = "> <Rosetta Interchangeability Group>\n";
		std::string group = restype.interchangeability_group() + "\n";

		lines.push_back(header3);
		lines.push_back(group);
		lines.push_back("\n");
	}
	if ( restype.aa() != core::chemical::aa_unk ) {
		std::string header4 = "> <Rosetta AA>\n";
		std::string aa = utility::to_string<core::chemical::AA>( restype.aa() ) + "\n";

		lines.push_back(header4);
		lines.push_back(aa);
		lines.push_back("\n");
	}

	return lines;
}

std::list<std::string> MolWriter::compose_rosetta_properties(core::conformation::ResidueCOP residue)
{
	core::chemical::ResidueType const & restype( residue->type() );
	std::list<std::string> lines;

	utility::vector1< std::string > const & properties( restype.properties().get_list_of_properties() );

	if ( properties.size() ) {
		std::string header = "> <Rosetta Properties>\n";
		lines.push_back(header);

		for ( core::Size ii(1); ii <= properties.size(); ++ii ) {
			std::string property_line = properties[ ii ] + "\n";
			lines.push_back( property_line );
		}

		lines.push_back( "\n" );
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
		lines.push_back("\n");
	}
	return lines;
}


}
}
}
