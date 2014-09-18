// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/methods/pcs2/PcsLineData.cc
 ///
 /// @brief  Class that hold a line data of the input file
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified February 2010
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsInputLine.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers
//#include <iostream>
#include <iomanip>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

static thread_local basic::Tracer TR_PcsInputLine( "protocols.scoring.methods.pcs.PcsInputLine" );

PcsInputLine::PcsInputLine() :
	residue_num_(0),
	atom_name_(""),
	PCS_experimental_(0),
	PCS_tolerance_(0)
{
	utility_exit_with_message( "You shouldn't call the empty constructor for PcsInputLine class" );
}

PcsInputLine::~PcsInputLine(){
}

PcsInputLine::PcsInputLine(PcsInputLine const & other):
	residue_num_(other.residue_num_),
	atom_name_(other.atom_name_),
	PCS_experimental_(other.PCS_experimental_),
	PCS_tolerance_(other.PCS_tolerance_)
{
}

PcsInputLine &
PcsInputLine::operator=( PcsInputLine const & other )
{
	if ( this != &other ) {
	//All data member are const, nothing to copy
	}
	return *this;
}

core::Size
PcsInputLine::get_residue_num() const{
	return residue_num_;
}

std::string
PcsInputLine::get_atom_name() const{
	return atom_name_;
}

core::Real
PcsInputLine::get_PCS_experimental() const{
	return PCS_experimental_;
}


core::Real
PcsInputLine::get_PCS_tolerance() const{
	 	return PCS_tolerance_;
}

PcsInputLine::PcsInputLine(core::Size residue_num,
													 std::string atom_name,
													 core::Real PCS_experimental,
													 core::Real PCS_tolerance
													 ) :
	residue_num_( residue_num ),
	atom_name_(atom_name),
	PCS_experimental_(PCS_experimental),
	PCS_tolerance_(PCS_tolerance)
{
}

std::ostream &
operator<<(std::ostream& out, const PcsInputLine &me){
		out << "Residue: " << std::setw(4) << me.get_residue_num();
		out << "   Atom: " << std::setw(4) << me.get_atom_name();
		out << "   PCS: " << std::setw(7) << me.get_PCS_experimental();
		out << "   Tolerance: " << std::setw(7) << me.get_PCS_tolerance()<< std::endl;
		return out;
}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
