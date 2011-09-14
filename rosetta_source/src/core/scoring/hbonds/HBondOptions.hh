// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondOptions.hh
/// @brief  HBondOptions class, holds the options for the hbond energy function
/// @author Matthew O'Meara

/// @detail
/// To add an additional option for hydrogen bonds do the following:
/// 1) add it to the default constructor
/// 2) add it to the copy constructor
/// 3) add a getter and a setter
/// 4) add it to operator==
/// 5) add it to the private data
/// 6) add it to HBondOptions::show


#ifndef INCLUDE_core_scoring_hbonds_HBondOptions_HH
#define INCLUDE_core_scoring_hbonds_HBondOptions_HH

// Unit headers
#include <core/scoring/hbonds/HBondOptions.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace core {
namespace scoring {
namespace hbonds {

class HBondOptions : public utility::pointer::ReferenceCount {
public:

	HBondOptions();

	HBondOptions( std::string params_db_tag );

	~HBondOptions();

	/// copy constructor
	HBondOptions( HBondOptions const & src );

	/// copy operator
	HBondOptions const &
	operator=( HBondOptions const & src );

	/// @brief Double counted hbonds include:
  /// @brief  - Hydrogen bonds to self
	/// @brief  - Backbone - sidechain hydrogen bonds where the backbone partner is forming a backbone - backbone hydrogen bond.
	/// @brief Turning off this exclusion rule is useful for collecting statistics on hydrogen bond site satisfaction
	bool
	exclude_self_hbonds() const;

	///
	void
	exclude_self_hbonds( bool const setting );

	///
	bool
	exclude_DNA_DNA() const;

	///
	void
	exclude_DNA_DNA( bool const setting );

	///
	bool
	use_hb_env_dep_DNA() const;

	///
	void
	use_hb_env_dep_DNA( bool const setting );

	///
	bool
	use_hb_env_dep() const;

	///
	void
	use_hb_env_dep( bool const setting );

	///
	bool
	smooth_hb_env_dep() const;

	///
	void
	smooth_hb_env_dep( bool const setting );

	///
	bool
	decompose_bb_hb_into_pair_energies() const;

	///
	void
	decompose_bb_hb_into_pair_energies( bool const setting );

	///
	//std::string const &
	//HBEval_fname() const;

	///
	void
	params_database_tag( std::string const & setting );

	///
	std::string const &
	params_database_tag() const;

	///pba
	bool
	Mbhbond() const;

	///pba
	void
	Mbhbond( bool const setting );

	///
	bool
	use_incorrect_deriv() const;

	///
	void
	use_incorrect_deriv( bool const setting );

	bool use_sp2_chi_penalty() const;

	void use_sp2_chi_penalty( bool setting );


	friend
	bool
	operator==( HBondOptions const & a, HBondOptions const & b );

	friend
	bool
	operator!=( HBondOptions const & a, HBondOptions const & b );

	friend
	std::ostream &
	operator<< ( std::ostream & out, const HBondOptions & options );

	///
	void
	show( std::ostream & out ) const;

private:

	bool exclude_DNA_DNA_;
	bool exclude_self_hbonds_;
	bool use_hb_env_dep_;
	bool use_hb_env_dep_DNA_;
	bool smooth_hb_env_dep_;
	bool decompose_bb_hb_into_pair_energies_;
	std::string params_database_tag_;
	bool use_incorrect_deriv_;
	bool use_sp2_chi_penalty_;
	bool Mbhbond_; //pba
};


} // hbonds
} // scoring
} // core

#endif // INCLUDE_core_scoring_hbonds_HBondOptions_HH
