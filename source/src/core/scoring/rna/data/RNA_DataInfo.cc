// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_DataInfo.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/data/RNA_DataInfo.hh>

// Package headers

// Project headers

// Utility headers

#include <utility/vector1.hh>


// C++

///////////////////////////////////////////////////////
// Keep track of information from, e.g., chemical
// accessibility experiments -- useful for scoring.
///////////////////////////////////////////////////////

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// ObjexxFCL serialization headers
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace rna {
namespace data {


/// @details Copy constructors must copy all data, not just some...
RNA_DataInfo::RNA_DataInfo( RNA_DataInfo const & src ):
	utility::pointer::ReferenceCount( src )
{
	*this = src;
}

RNA_DataInfo &
RNA_DataInfo::operator = ( RNA_DataInfo const & src ){

	rna_data_         = src.rna_data_;
	rna_reactivities_ = src.rna_reactivities_;
	backbone_burial_  = src.backbone_burial_;
	backbone_exposed_ = src.backbone_exposed_;
	return  *this;

}

////////////////////////////////////////////////////////
void
RNA_DataInfo::zero()
{
	rna_data_.clear();
	rna_reactivities_.clear();
}


} //data
} //rna
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::data::RNA_DataInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( rna_data_ ) ); // RNA_Data
	arc( CEREAL_NVP( backbone_burial_ ) ); // ObjexxFCL::FArray1D<_Bool>
	arc( CEREAL_NVP( backbone_exposed_ ) ); // ObjexxFCL::FArray1D<_Bool>
	arc( CEREAL_NVP( rna_reactivities_ ) ); // RNA_Reactivities
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::data::RNA_DataInfo::load( Archive & arc ) {
	arc( rna_data_ ); // RNA_Data
	arc( backbone_burial_ ); // ObjexxFCL::FArray1D<_Bool>
	arc( backbone_exposed_ ); // ObjexxFCL::FArray1D<_Bool>
	arc( rna_reactivities_ ); // RNA_Reactivities
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::data::RNA_DataInfo );
CEREAL_REGISTER_TYPE( core::scoring::rna::data::RNA_DataInfo )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::rna::data::RNA_Datum::RNA_Datum() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::data::RNA_Datum::save( Archive & arc ) const {
	arc( CEREAL_NVP( position_ ) ); // Size
	arc( CEREAL_NVP( edge_ ) ); // Size
	arc( CEREAL_NVP( weight_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::data::RNA_Datum::load( Archive & arc ) {
	arc( position_ ); // Size
	arc( edge_ ); // Size
	arc( weight_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::data::RNA_Datum );

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::rna::data::RNA_Reactivity::RNA_Reactivity() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::data::RNA_Reactivity::save( Archive & arc ) const {
	arc( CEREAL_NVP( position_ ) ); // Size
	arc( CEREAL_NVP( type_ ) ); // enum core::scoring::rna::data::RNA_ReactivityType
	arc( CEREAL_NVP( value_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::data::RNA_Reactivity::load( Archive & arc ) {
	arc( position_ ); // Size
	arc( type_ ); // enum core::scoring::rna::data::RNA_ReactivityType
	arc( value_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::data::RNA_Reactivity );
CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_rna_data_RNA_DataInfo )
#endif // SERIALIZATION
