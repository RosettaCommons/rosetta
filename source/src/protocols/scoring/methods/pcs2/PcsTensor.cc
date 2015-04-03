// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 ///
 /// @file protocols/scoring/methods/pcs2/PcsTensor.cc
 ///
 /// @brief Hold chi-tensor information of PCS
 ///
 /// @details
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references JBNMR 2008  41:179-189 schmitz et all will explains the tensor convention used
 ///
 /// @authorv Christophe Schmitz
 ///
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsTensor.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <iostream>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

PcsTensor::PcsTensor(){
	utility_exit_with_message( "You shouldn't call the empty constructor for PcsTensor" );
}

PcsTensor::~PcsTensor(){
}

PcsTensor::PcsTensor(PcsTensor const & other){
	a_ = other.a_;
	b_ = other.b_;
	g_ = other.g_;
	ax_ = other.ax_;
	rh_ = other.rh_;
	chi_xx_ = other.chi_xx_;
	chi_yy_ = other.chi_yy_;
	chi_xy_ = other.chi_xy_;
	chi_xz_ = other.chi_xz_;
	chi_yz_ = other.chi_yz_;
	label_ = other.label_;
}

PcsTensor &
PcsTensor::operator=(PcsTensor const & other){

	if ( this != &other ) {
		a_ = other.a_;
		b_ = other.b_;
		g_ = other.g_;
		ax_ = other.ax_;
		rh_ = other.rh_;
		chi_xx_ = other.chi_xx_;
		chi_yy_ = other.chi_yy_;
		chi_xy_ = other.chi_xy_;
		chi_xz_ = other.chi_xz_;
		chi_yz_ = other.chi_yz_;
		label_ = other.label_;
	}
	return *this;
}

///////////////////////////////////////////////
/// @brief The constructeur use the chi matrix parameters (not the alpha beta gamma Ax and Rh component...)
///////////////////////////////////////////////
PcsTensor::PcsTensor(core::Real const chi_xx,
										 core::Real const chi_xy,
										 core::Real const chi_xz,
										 core::Real const chi_yy,
										 core::Real const chi_yz,
										 std::string const label)
{
	chi_xx_ = chi_xx;
	chi_yy_ = chi_yy;
	chi_xy_ = chi_xy;
	chi_xz_ = chi_xz;
	chi_yz_ = chi_yz;
	label_ = label;

	set_abgar();
}

void
PcsTensor::reset_tensor(core::Real const chi_xx,
												core::Real const chi_xy,
												core::Real const chi_xz,
												core::Real const chi_yy,
												core::Real const chi_yz){
	chi_xx_ = chi_xx;
	chi_yy_ = chi_yy;
	chi_xy_ = chi_xy;
	chi_xz_ = chi_xz;
	chi_yz_ = chi_yz;

	set_abgar();
}

std::ostream &
operator<<(std::ostream& out, const PcsTensor &me ){

	out << "The tensor parameters of " << "'" <<me.label_ << "' are:" << std::endl;
	out << "Delta Chi xx: " << me.chi_xx_ << std::endl;
	out << "Delta Chi yy: " << me.chi_yy_ << std::endl;
	out << "              => Delta Chi zz: " << me.get_chi_zz() << std::endl;
	out << "Delta Chi xy: " << me.chi_xy_ << std::endl;
	out << "Delta Chi xz: " << me.chi_xz_ << std::endl;
	out << "Delta Chi yz: " << me.chi_yz_ << std::endl;
	/*
	out << "In the UTR representation: " << std::endl;
	out << "Alpha:   "  << me.a_ << std::endl;
	out << "Beta:    "  << me.b_ << std::endl;
	out << "Gamma:   "  << me.g_ << std::endl;
	out << "Axial:   "  << me.ax_ << std::endl;
	out << "Rhombic: "  << me.rh_ << std::endl;
	out << "        => Xxx:   "  << me.get_delta_X_xx() << std::endl;
	out << "        => Xyy:   "  << me.get_delta_X_yy() << std::endl;
	out << "        => Xzz:   "  << me.get_delta_X_zz() << std::endl;
	*/

	return out;
}

void
PcsTensor::reset_from_ref(PcsTensor & other){
	a_ = other.a_;
  b_ = other.b_;
  g_ = other.g_;
  ax_ = other.ax_;
  rh_ = other.rh_;

  chi_xx_ = other.chi_xx_;
  chi_yy_ = other.chi_yy_;
  chi_xy_ = other.chi_xy_;
  chi_xz_ = other.chi_xz_;
  chi_yz_ = other.chi_yz_;

	label_ = other.label_;
}

void
PcsTensor::set_abgar(){
	//This has to be updated once I get the translation right
	//from the chi parameters, i should get the a b g ax rh (abgar) parameters
	//TODO
	a_ = 0;
	b_ = 0;
	g_ = 0;
	ax_ = 0;
	rh_ = 0;
}

std::string const &
PcsTensor::get_label() const{
	return label_;
}

core::Real
PcsTensor::get_delta_X_xx() const{
	return(rh_ / 2.0 - ax_ / 3.0);
}

core::Real
PcsTensor::get_delta_X_yy() const{
	return(-rh_ / 2.0 - ax_ / 2.0);
}

core::Real
PcsTensor::get_delta_X_zz() const{
	return(2.0/3.0 * ax_);
}

core::Real
PcsTensor::get_chi_zz() const{
	return(-chi_xx_ -chi_yy_);
}


core::Real
PcsTensor::get_chi_xx() const{
	return(chi_xx_);
}

core::Real
PcsTensor::get_chi_yy() const{
	return(chi_yy_);
}

core::Real
PcsTensor::get_chi_xy() const{
	return(chi_xy_);
}

core::Real
PcsTensor::get_chi_xz() const{
	return(chi_xz_);
}

core::Real
PcsTensor::get_chi_yz() const{
	return(chi_yz_);
}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
