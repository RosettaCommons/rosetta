// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsInputFile.cc
///
/// @brief Read all input from a .npc input file, and hold the data in the class
/// One file per lanthanide data
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsInputFile.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static basic::Tracer TR_PcsInputFile( "protocols.scoring.methods.pcs.PcsInputFile" );

PcsInputFile::PcsInputFile():
	filename_(""), weight_(0)
{
	utility_exit_with_message( "You shouldn't call the empty constructor for PcsInputFile class" );
}

PcsInputFile::~PcsInputFile(){
}

PcsInputFile::PcsInputFile(PcsInputFile const & other):
	filename_(other.filename_), weight_(other.weight_)
{
	PcsInputLine_all_ = other.PcsInputLine_all_;
}

PcsInputFile &
PcsInputFile::operator=( PcsInputFile const & other ){
	if ( this != &other ) {
		PcsInputLine_all_ = other.PcsInputLine_all_;
	}
	return *this;
}

std::string
PcsInputFile::get_filename() const{
	return(filename_);
}

core::Real
PcsInputFile::get_weight() const{
	return(weight_);
}

PcsInputFile::PcsInputFile(std::string const & filename, core::Real const my_weight):
	filename_(std::string(filename)), weight_(my_weight)
{
	read_PCS_file();
}

void
PcsInputFile::read_PCS_file(){
	core::Size residue_num;
	std::string atom_name;
	core::Real PCS_experimental;
	core::Real PCS_tolerance;
	std::ifstream myfile;
	std::string line;
	core::Size line_number(0);

	TR_PcsInputFile << "Opening file '" << get_filename().c_str() << "'" << std::endl;
	myfile.open (get_filename().c_str(), std::ios::in);
	if ( !myfile.is_open () ) {
		std::cerr << "Unable to open the file '" << get_filename().c_str()  <<"'" << std::endl;
		utility_exit();
	}

	while ( getline( myfile, line ) ) {
		line_number++;
		std::istringstream line_stream( line ,std::istringstream::in);
		if ( (line_stream >> residue_num >> atom_name >> PCS_experimental >> PCS_tolerance).fail() ) {
			TR_PcsInputFile << "Ignoring line " <<line_number << ": `" << line << "` from file " << get_filename() <<std::endl;
			continue;
		}

		PcsInputLine_all_.push_back( PcsInputLine( residue_num, atom_name, PCS_experimental, PCS_tolerance ) );
	}

	myfile.close();
}

utility::vector1<PcsInputLine> &
PcsInputFile::get_PcsInputLine_all(){
	return (PcsInputLine_all_);
}

std::ostream &
operator<<(std::ostream& out, const PcsInputFile &me){
	utility::vector1<PcsInputLine>::iterator it;
	utility::vector1<PcsInputLine> pcs_i_l_a;
	pcs_i_l_a = me.PcsInputLine_all_;
	core::Size i, n;

	i = 1;
	n = pcs_i_l_a.size();
	out<< "Found the following " << n << " PCS value in the file " << me.get_filename() << std::endl;
	for ( it = pcs_i_l_a.begin(); it != pcs_i_l_a.end(); ++it ) {

		if ( (i == 1) || (i == 2) ) {
			out << *it;
			i++;
			continue;
		}
		if ( i == (n-1) ) {
			out << std::endl;
			out << *it;
			i++;
			continue;
		}
		if ( i == (n) ) {
			out << *it;
			i++;
			continue;
		}
		out << ".";
		i++;
	}

	return out;
}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
