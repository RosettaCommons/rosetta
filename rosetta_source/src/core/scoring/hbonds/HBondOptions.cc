// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondOptions.cc
/// @brief  HBondOptions class, holds the options for the hbond energy function
/// @author Matthew O'Meara

/// @detail
/// To add an additional option for hydrogen bonds do the following:
///
/// In HBondOptions.hh:
/// 1) add it to the default constructor
/// 2) add it to the copy constructor
/// 3) add a getter and a setter
/// 4) add it to operator==
/// 5) add it to the private data
/// 6) add it to HBondOptions::show

// Unit Headers
#include <core/scoring/hbonds/HBondOptions.hh>

// Package Headers
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh> //pba
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh> //pba

// C++ Headers
#include <string>
#include <iostream>


namespace core {
namespace scoring {
namespace hbonds {

HBondOptions::HBondOptions( std::string params_db_tag ):
	exclude_DNA_DNA_( true ),
	use_hb_env_dep_ ( true ),
	use_hb_env_dep_DNA_( true ),
	smooth_hb_env_dep_( true ),
	decompose_bb_hb_into_pair_energies_( false ),
	params_database_tag_(params_db_tag),
	use_incorrect_deriv_( true ),
	use_sp2_chi_penalty_( false ),
	Mbhbond_( false ) //pba
{
	using namespace basic::options;
	if (option.has(OptionKeys::membrane::Mhbond_depth) &&
		option[OptionKeys::membrane::Mhbond_depth].user()){
		Mbhbond_ = option[OptionKeys::membrane::Mhbond_depth];//pba
	}
	use_incorrect_deriv_ = option[OptionKeys::corrections::score::use_incorrect_hbond_deriv];
	use_sp2_chi_penalty_ = option[OptionKeys::corrections::score::hb_sp2_chipen ];
}



HBondOptions::HBondOptions():
	exclude_DNA_DNA_( true ),
	use_hb_env_dep_ ( true ),
	use_hb_env_dep_DNA_( true ),
	smooth_hb_env_dep_( true ),
	decompose_bb_hb_into_pair_energies_( false ),
	params_database_tag_("standard_params"),
	use_incorrect_deriv_( false ),
	use_sp2_chi_penalty_( false ),
	Mbhbond_( false ) //pba
{
	using namespace basic::options;
	if (option.has(OptionKeys::score::hbond_params) &&
		option[OptionKeys::score::hbond_params].user()){
		params_database_tag_ = option[OptionKeys::score::hbond_params];
	}
	if (option.has(OptionKeys::membrane::Mhbond_depth) &&
		option[OptionKeys::membrane::Mhbond_depth].user()){
		Mbhbond_ = option[OptionKeys::membrane::Mhbond_depth];//pba
	}
	use_incorrect_deriv_ = option[OptionKeys::corrections::score::use_incorrect_hbond_deriv];
	use_sp2_chi_penalty_ = option[OptionKeys::corrections::score::hb_sp2_chipen ];
}

/// copy constructor
HBondOptions::HBondOptions( HBondOptions const & src ):
	ReferenceCount( src )
{
	*this = src;
}

HBondOptions::~HBondOptions(){}

/// copy operator
HBondOptions const &
HBondOptions::operator=( HBondOptions const & src )
{
	exclude_DNA_DNA_ = src.exclude_DNA_DNA_;
	use_hb_env_dep_ = src.use_hb_env_dep_;
	use_hb_env_dep_DNA_ = src.use_hb_env_dep_DNA_;
	smooth_hb_env_dep_ = src.smooth_hb_env_dep_;
	decompose_bb_hb_into_pair_energies_ = src.decompose_bb_hb_into_pair_energies_;
	params_database_tag_ = src.params_database_tag_;
	use_incorrect_deriv_ = src.use_incorrect_deriv_;
	use_sp2_chi_penalty_ = src.use_sp2_chi_penalty_;
	Mbhbond_ = src.Mbhbond_; //pba
	return *this;
}

///
bool
HBondOptions::exclude_DNA_DNA() const
{
	return exclude_DNA_DNA_;
}

///
void
HBondOptions::exclude_DNA_DNA( bool const setting )
{
	exclude_DNA_DNA_ = setting;
}

///
bool
HBondOptions::use_hb_env_dep_DNA() const
{
	return use_hb_env_dep_DNA_;
}

///
void
HBondOptions::use_hb_env_dep_DNA( bool const setting )
{
	use_hb_env_dep_DNA_ = setting;
}

///
bool
HBondOptions::use_hb_env_dep() const
{
	return use_hb_env_dep_;
}

///
void
HBondOptions::use_hb_env_dep( bool const setting )
{
	use_hb_env_dep_ = setting;
}

///
bool
HBondOptions::smooth_hb_env_dep() const
{
	return smooth_hb_env_dep_;
}

///
void
HBondOptions::smooth_hb_env_dep( bool const setting )
{
	smooth_hb_env_dep_ = setting;
}

///
bool
HBondOptions::decompose_bb_hb_into_pair_energies() const
{
	return decompose_bb_hb_into_pair_energies_;
}

///
void
HBondOptions::decompose_bb_hb_into_pair_energies( bool const setting )
{
	decompose_bb_hb_into_pair_energies_ = setting;
}

///
std::string const &
HBondOptions::params_database_tag() const
{
	return params_database_tag_;
}

///
void
HBondOptions::params_database_tag( std::string const & setting )
{
	params_database_tag_ = setting;
}

///pba
bool
HBondOptions::Mbhbond() const
{
	return Mbhbond_;
}

///
void
HBondOptions::Mbhbond( bool const setting )
{
	Mbhbond_ = setting;
}

///
bool
HBondOptions::use_incorrect_deriv() const
{
	return use_incorrect_deriv_;
}

///
void
HBondOptions::use_incorrect_deriv( bool const setting )
{
	use_incorrect_deriv_ = setting;
}

bool HBondOptions::use_sp2_chi_penalty() const
{
	return use_sp2_chi_penalty_;
}

void HBondOptions::use_sp2_chi_penalty( bool setting )
{
	use_sp2_chi_penalty_ = setting;
}


bool
operator==( HBondOptions const & a, HBondOptions const & b )
{
	return ( ( a.exclude_DNA_DNA_ == b.exclude_DNA_DNA_ ) &&
		( a.use_hb_env_dep_ == b.use_hb_env_dep_ ) &&
		( a.use_hb_env_dep_DNA_ == b.use_hb_env_dep_DNA_ ) &&
		( a.smooth_hb_env_dep_ == b.smooth_hb_env_dep_ ) &&
		( a.decompose_bb_hb_into_pair_energies_ == b.decompose_bb_hb_into_pair_energies_ ) &&
		( a.params_database_tag_ == b.params_database_tag_ ) &&
		( a.use_incorrect_deriv_ == b.use_incorrect_deriv_ ) &&
		( a.use_sp2_chi_penalty_ == b.use_sp2_chi_penalty_ ) &&
		( a.Mbhbond_ == b.Mbhbond_ ) ); //pba
}

bool
operator!=( HBondOptions const & a, HBondOptions const & b )
{
	return !( a == b );
}

std::ostream &
operator<< ( std::ostream & out, const HBondOptions & options ){
	options.show( out );
	return out;
}

///
void
HBondOptions::show( std::ostream & out ) const
{
	out <<"HBondOptions::show: exclude_DNA_DNA: "
		<<( exclude_DNA_DNA_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: use_hb_env_dep: "
		<<( use_hb_env_dep_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: use_hb_env_dep_DNA: "
		<<( use_hb_env_dep_DNA_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: smooth_hb_env_dep: "
		<<( smooth_hb_env_dep_ ? "true" : "false " ) << std::endl;
	out <<"HBondOptions::show: decompose_bb_hb_into_pair_energies: "
		<<( decompose_bb_hb_into_pair_energies_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: params_database_tag_: "
		<< params_database_tag_ << std::endl;
	out <<"HBondOptions::show: use_incorrect_deriv_: "
		<<( use_incorrect_deriv_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: use_sp2_chi_penalty_: "
		<<( use_sp2_chi_penalty_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: Mbhbond: "
		<<( Mbhbond_ ? "true" : "false " ) << std::endl; //pba
}

}
}
}
