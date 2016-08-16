// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBondTypeManager.hh
/// @brief  HBond enumeration type manager
/// @author Matthew O'Meara


#ifndef INCLUDED_core_scoring_hbonds_HBondTypeManager_hh
#define INCLUDED_core_scoring_hbonds_HBondTypeManager_hh

// Package headers
#include <core/scoring/hbonds/types.hh>
#include <core/chemical/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace hbonds {

class HBondTypeManager {
public:
	////////////////////
	/// Bond Weight Type
	////////////////////
	static
	HBondWeightType
	weight_type_from_name( std::string const & name );

	static
	std::string
	name_from_weight_type( HBondWeightType score_type );

	static
	bool
	is_weight_type( std::string const & name );


	///////////////////
	/// Derivative Type
	///////////////////
	static
	HBDerivType
	deriv_type_from_name( std::string const & name );

	static
	std::string
	name_from_deriv_type( HBDerivType score_type );

	static
	bool
	is_deriv_type( std::string const & name );

	///////////////////////
	/// Donor Chemical Type
	///////////////////////
	static
	HBDonChemType
	don_chem_type_from_name( std::string const & name );

	static
	std::string
	name_from_don_chem_type( HBDonChemType score_type );

	static
	bool
	is_don_chem_type( std::string const & name );

	//////////////////////////
	/// Acceptor Chemical Type
	//////////////////////////
	static
	HBAccChemType
	acc_chem_type_from_name( std::string const & name );

	static
	std::string
	name_from_acc_chem_type( HBAccChemType score_type );

	static
	bool
	is_acc_chem_type( std::string const & name );

	////////////////////////////
	/// Sequence Separation Type
	////////////////////////////
	static
	HBSeqSep
	seq_sep_type_from_name( std::string const & name );

	static
	std::string
	name_from_seq_sep_type( HBSeqSep score_type );

	static
	bool
	is_seq_sep_type( std::string const & name );

	////////////////////////
	//Acceptor Hybridization
	////////////////////////
	static
	chemical::Hybridization
	hybridization_type_from_name( std::string const & name );

	static
	std::string
	name_from_hybridization_type( chemical::Hybridization );

	static
	bool
	is_hybridization_type( std::string const & name);

	///////////////////////////
	///Geometric Dimension Type
	///////////////////////////
	static
	HBGeoDimType
	geo_dim_type_from_name( std::string const & name );

	static
	std::string
	name_from_geo_dim_type( HBGeoDimType score_type );

	static
	bool
	is_geo_dim_type( std::string const & name );


private:
	static void setup_type_names();

private:
	static bool initialized_;

	/// lookup map from string name to enum type
	static std::map< std::string, HBondWeightType > name2weight_type_;
	static utility::vector1< std::string > weight_type2name_;

	static std::map< std::string, HBDerivType > name2deriv_type_;
	static utility::vector1< std::string > deriv_type2name_;

	static std::map< std::string, HBDonChemType > name2don_chem_type_;
	static utility::vector1< std::string > don_chem_type2name_;

	static std::map< std::string, HBAccChemType > name2acc_chem_type_;
	static utility::vector1< std::string > acc_chem_type2name_;

	static std::map< std::string, HBSeqSep > name2seq_sep_type_;
	static utility::vector1< std::string > seq_sep_type2name_;

	static std::map< std::string, chemical::Hybridization > name2hybridization_type_;
	static utility::vector1< std::string > hybridization_type2name_;

	static std::map< std::string, HBGeoDimType > name2geo_dim_type_;
	static utility::vector1< std::string > geo_dim_type2name_;

}; // HBondTypeManager

} // hbonds
} // scoring
} // core

#endif
