// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/hbonds/HBondTypeManager.hh
/// @brief HBond enumeration type manager
/// @author Matthew O'Meara

// Unit headers
#include <core/scoring/hbonds/HBondTypeManager.hh>

// Project headers
#include <core/scoring/hbonds/types.hh>
#include <core/chemical/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1_bool.hh>

// C++ headers
#include <map>
#include <string>
#include <iostream>

#include <utility/vector1.hh>




namespace core {
namespace scoring {
namespace hbonds{

using namespace std;
using utility::vector1;
using namespace chemical;

bool HBondTypeManager::initialized_( false );

map< string, HBondWeightType > HBondTypeManager::name2weight_type_;
vector1< string > HBondTypeManager::weight_type2name_;

map< string, HBDerivType > HBondTypeManager::name2deriv_type_;
vector1< string > HBondTypeManager::deriv_type2name_;

map< string, HBDonChemType > HBondTypeManager::name2don_chem_type_;
vector1< string > HBondTypeManager::don_chem_type2name_;

map< string, HBAccChemType > HBondTypeManager::name2acc_chem_type_;
vector1< string > HBondTypeManager::acc_chem_type2name_;

map< string, HBSeqSep > HBondTypeManager::name2seq_sep_type_;
vector1< string > HBondTypeManager::seq_sep_type2name_;

map< string, Hybridization > HBondTypeManager::name2hybridization_type_;
vector1< string > HBondTypeManager::hybridization_type2name_;

map< string, HBGeoDimType > HBondTypeManager::name2geo_dim_type_;
vector1< string > HBondTypeManager::geo_dim_type2name_;


/// @brief initialize the ScoreType name vector and map
///
/// @details initialize all the SCORETYPE string name into the vector then set up
/// the look-up map from string name to enum type
void
HBondTypeManager::setup_type_names()
{
	if ( initialized_ ) return;
	initialized_ = true;


	name2weight_type_["hbw_NONE"] = hbw_NONE;
	name2weight_type_["hbw_SR_BB"] = hbw_SR_BB;
	name2weight_type_["hbw_LR_BB"] = hbw_LR_BB;
	name2weight_type_["hbw_SR_BB_SC"] = hbw_SR_BB_SC;
	name2weight_type_["hbw_LR_BB_SC"] = hbw_LR_BB_SC;
	name2weight_type_["hbw_SC"] = hbw_SC;

	name2deriv_type_["hbderiv_NONE"] = hbderiv_NONE;
	name2deriv_type_["hbderiv_ABE_GO"] = hbderiv_ABE_GO;
	name2deriv_type_["hbderiv_ABE_GO_GEOMSOL_OCC_ACC"] = hbderiv_ABE_GO_GEOMSOL_OCC_ACC;
	name2deriv_type_["hbderiv_ABE_GO_GEOMSOL_OCC_DON"] = hbderiv_ABE_GO_GEOMSOL_OCC_DON;

	name2don_chem_type_["hbdon_NONE"] = hbdon_NONE;
	name2don_chem_type_["hbdon_PBA"] = hbdon_PBA;
	name2don_chem_type_["hbdon_CXA"] = hbdon_CXA;
	name2don_chem_type_["hbdon_IMD"] = hbdon_IMD;
	name2don_chem_type_["hbdon_IME"] = hbdon_IME;
	name2don_chem_type_["hbdon_IND"] = hbdon_IND;
	name2don_chem_type_["hbdon_AMO"] = hbdon_AMO;
	name2don_chem_type_["hbdon_GDE"] = hbdon_GDE;
	name2don_chem_type_["hbdon_GDH"] = hbdon_GDH;
	name2don_chem_type_["hbdon_AHX"] = hbdon_AHX;
	name2don_chem_type_["hbdon_HXL"] = hbdon_HXL;
	name2don_chem_type_["hbdon_H2O"] = hbdon_H2O;
	name2don_chem_type_["hbdon_GENERIC_BB"] = hbdon_GENERIC_BB;
	name2don_chem_type_["hbdon_GENERIC_SC"] = hbdon_GENERIC_SC;

	name2acc_chem_type_["hbacc_NONE"] = hbacc_NONE;
	name2acc_chem_type_["hbacc_PBA"] = hbacc_PBA;
	name2acc_chem_type_["hbacc_CXA"] = hbacc_CXA;
	name2acc_chem_type_["hbacc_CXL"] = hbacc_CXL;
	name2acc_chem_type_["hbacc_IMD"] = hbacc_IMD;
	name2acc_chem_type_["hbacc_IME"] = hbacc_IME;
	name2acc_chem_type_["hbacc_AHX"] = hbacc_AHX;
	name2acc_chem_type_["hbacc_HXL"] = hbacc_HXL;
	name2acc_chem_type_["hbacc_PCA_DNA"] = hbacc_PCA_DNA;
	name2acc_chem_type_["hbacc_PES_DNA"] = hbacc_PES_DNA;
	name2acc_chem_type_["hbacc_RRI_DNA"] = hbacc_RRI_DNA;
	name2acc_chem_type_["hbacc_PCA_RNA"] = hbacc_PCA_RNA;
	name2acc_chem_type_["hbacc_PES_RNA"] = hbacc_PES_RNA;
	name2acc_chem_type_["hbacc_RRI_RNA"] = hbacc_RRI_RNA;
	name2acc_chem_type_["hbacc_H2O"] = hbacc_H2O;
	name2acc_chem_type_["hbacc_GENERIC_SP2BB"] = hbacc_GENERIC_SP2BB;
	name2acc_chem_type_["hbacc_GENERIC_SP2SC"] = hbacc_GENERIC_SP2SC;
	name2acc_chem_type_["hbacc_GENERIC_SP3BB"] = hbacc_GENERIC_SP3BB;
	name2acc_chem_type_["hbacc_GENERIC_SP3SC"] = hbacc_GENERIC_SP3SC;
	name2acc_chem_type_["hbacc_GENERIC_RINGBB"] = hbacc_GENERIC_RINGBB;
	name2acc_chem_type_["hbacc_GENERIC_RINGSC"] = hbacc_GENERIC_RINGSC;

	name2seq_sep_type_["seq_sep_other"] = seq_sep_other;
	name2seq_sep_type_["seq_sep_M4"] = seq_sep_M4;
	name2seq_sep_type_["seq_sep_M3"] = seq_sep_M3;
	name2seq_sep_type_["seq_sep_M2"] = seq_sep_M2;
	name2seq_sep_type_["seq_sep_PM1"] = seq_sep_PM1;
	name2seq_sep_type_["seq_sep_P2"] = seq_sep_P2;
	name2seq_sep_type_["seq_sep_P3"] = seq_sep_P3;
	name2seq_sep_type_["seq_sep_P4"] = seq_sep_P4;

	name2hybridization_type_["SP2_HYBRID"] = chemical::SP2_HYBRID;
	name2hybridization_type_["SP3_HYBRID"] = chemical::SP3_HYBRID;
	name2hybridization_type_["RING_HYBRID"] = chemical::RING_HYBRID;
	name2hybridization_type_["UNKNOWN_HYBRID"] = chemical::UNKNOWN_HYBRID;

	name2geo_dim_type_["hbgd_NONE"] = hbgd_NONE;
	name2geo_dim_type_["hbgd_AHdist"] = hbgd_AHdist;
	name2geo_dim_type_["hbgd_cosBAH"] = hbgd_cosBAH;
	name2geo_dim_type_["hbgd_cosAHD"] = hbgd_cosAHD;
	name2geo_dim_type_["hbgd_AHD"] = hbgd_AHD;
	name2geo_dim_type_["hbgd_chi"] = hbgd_chi;

	assert( name2weight_type_.size() == hbw_MAX );
	assert( name2deriv_type_.size() == hbderiv_MAX );
	assert( name2don_chem_type_.size() == hbdon_MAX );
	assert( name2acc_chem_type_.size() == hbacc_MAX );
	assert( name2seq_sep_type_.size() == seq_sep_MAX );
	assert( name2geo_dim_type_.size() == hbgd_MAX );
	assert( name2hybridization_type_.size() == chemical::HYBRID_MAX );

	weight_type2name_.resize( hbw_MAX );
	for ( map< string, HBondWeightType >::const_iterator iter = name2weight_type_.begin(),
		iter_end = name2weight_type_.end(); iter != iter_end; ++iter ) {
		weight_type2name_[ iter->second ] = iter->first;
	}

	deriv_type2name_.resize( hbderiv_MAX );
	for ( map< string, HBDerivType >::const_iterator iter = name2deriv_type_.begin(),
		iter_end = name2deriv_type_.end(); iter != iter_end; ++iter ) {
		deriv_type2name_[ iter->second ] = iter->first;
	}

	don_chem_type2name_.resize(hbdon_MAX );
	for ( map< string, HBDonChemType >::const_iterator iter = name2don_chem_type_.begin(),
		iter_end = name2don_chem_type_.end(); iter != iter_end; ++iter ) {
		don_chem_type2name_[ iter->second ] = iter->first;
	}

	acc_chem_type2name_.resize(hbacc_MAX );
	for ( map< string, HBAccChemType >::const_iterator iter = name2acc_chem_type_.begin(),
		iter_end = name2acc_chem_type_.end(); iter != iter_end; ++iter ) {
		acc_chem_type2name_[ iter->second ] = iter->first;
	}

	seq_sep_type2name_.resize(seq_sep_MAX );
	for ( map< string, HBSeqSep >::const_iterator iter = name2seq_sep_type_.begin(),
		iter_end = name2seq_sep_type_.end(); iter != iter_end; ++iter ) {
		seq_sep_type2name_[ iter->second ] = iter->first;
	}

	hybridization_type2name_.resize(chemical::HYBRID_MAX);
	for ( map< string, Hybridization >::const_iterator
					iter = name2hybridization_type_.begin(),
					iter_end = name2hybridization_type_.end();
				iter != iter_end; ++iter ) {
		hybridization_type2name_[iter->second ] = iter->first;
	}

	geo_dim_type2name_.resize(hbgd_MAX );
	for ( map< string, HBGeoDimType >::const_iterator iter = name2geo_dim_type_.begin(),
		iter_end = name2geo_dim_type_.end(); iter != iter_end; ++iter ) {
		geo_dim_type2name_[ iter->second ] = iter->first;
	}
}

//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBondWeightType string name and return its enum type
HBondWeightType
HBondTypeManager::weight_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBondWeightType >::const_iterator iter( name2weight_type_.find( name ) );
	if ( iter == name2weight_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond weight type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_weight_type( HBondWeightType type )
{
	setup_type_names();
	return weight_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_weight_type( string const & name )
{
	setup_type_names();
	map< string, HBondWeightType >::const_iterator iter( name2weight_type_.find( name ) );
	return iter != name2weight_type_.end();
}



//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBDerivType string name and return its enum type
HBDerivType
HBondTypeManager::deriv_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBDerivType >::const_iterator iter( name2deriv_type_.find( name ) );
	if ( iter == name2deriv_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond deriv type "+name);
	}
	return iter->second;
}

string
HBondTypeManager::name_from_deriv_type( HBDerivType type )
{
	setup_type_names();
	return deriv_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_deriv_type( string const & name )
{
	setup_type_names();
	map< string, HBDerivType >::const_iterator iter( name2deriv_type_.find( name ) );
	return iter != name2deriv_type_.end();
}




//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBDonChemType string name and return its enum type
HBDonChemType
HBondTypeManager::don_chem_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBDonChemType >::const_iterator iter( name2don_chem_type_.find( name ) );
	if ( iter == name2don_chem_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond donor chemical type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_don_chem_type( HBDonChemType type )
{
	setup_type_names();
	return don_chem_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_don_chem_type( string const & name )
{
	setup_type_names();
	map< string, HBDonChemType >::const_iterator iter( name2don_chem_type_.find( name ) );
	return iter != name2don_chem_type_.end();
}


//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBAccChemType string name and return its enum type
HBAccChemType
HBondTypeManager::acc_chem_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBAccChemType >::const_iterator iter( name2acc_chem_type_.find( name ) );
	if ( iter == name2acc_chem_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond acceptor chemical type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_acc_chem_type( HBAccChemType type )
{
	setup_type_names();
	return acc_chem_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_acc_chem_type( string const & name )
{
	setup_type_names();
	map< string, HBAccChemType >::const_iterator iter( name2acc_chem_type_.find( name ) );
	return iter != name2acc_chem_type_.end();
}

//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBSeqSep string name and return its enum type
HBSeqSep
HBondTypeManager::seq_sep_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBSeqSep >::const_iterator iter( name2seq_sep_type_.find( name ) );
	if ( iter == name2seq_sep_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond sequence separation type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_seq_sep_type( HBSeqSep type )
{
	setup_type_names();
	return seq_sep_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_seq_sep_type( string const & name )
{
	setup_type_names();
	map< string, HBSeqSep >::const_iterator iter( name2seq_sep_type_.find( name ) );
	return iter != name2seq_sep_type_.end();
}

//////////////////////////////////////////////////////////////////////////////
/// @brief given a chemical::Hybridization string name return its enum type
/// In the perfect world this would live in the chemical/ChemicalTypeManager.hh,
/// but the chemical types are currently not managed!
chemical::Hybridization
HBondTypeManager::hybridization_type_from_name( string const & name )
{
	setup_type_names();
	map< string, Hybridization >::const_iterator iter( name2hybridization_type_.find( name ) );
	if ( iter == name2hybridization_type_.end() ) {
		utility_exit_with_message("unrecognized chemical hybridization type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_hybridization_type( Hybridization type )
{
	setup_type_names();
	return hybridization_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_hybridization_type( string const & name )
{
	setup_type_names();
	map< string, Hybridization >::const_iterator iter( name2hybridization_type_.find( name ) );
	return iter != name2hybridization_type_.end();
}

//////////////////////////////////////////////////////////////////////////////
/// @brief give a HBGeoDimType string name and return its enum type
HBGeoDimType
HBondTypeManager::geo_dim_type_from_name( string const & name )
{
	setup_type_names();
	map< string, HBGeoDimType >::const_iterator iter( name2geo_dim_type_.find( name ) );
	if ( iter == name2geo_dim_type_.end() ) {
		utility_exit_with_message("unrecognized hydrogen bond geometric dimension type '"+name+"'");
	}
	return iter->second;
}

string
HBondTypeManager::name_from_geo_dim_type( HBGeoDimType type )
{
	setup_type_names();
	return geo_dim_type2name_[ type ];
}

///@brief
bool
HBondTypeManager::is_geo_dim_type( string const & name )
{
	setup_type_names();
	map< string, HBGeoDimType >::const_iterator iter( name2geo_dim_type_.find( name ) );
	return iter != name2geo_dim_type_.end();
}



} // hbonds
} // scoring
} // core
