// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_DataPotential.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_DataInfo_hh
#define INCLUDED_core_scoring_rna_RNA_DataInfo_hh

#include <core/types.hh>
#include <core/scoring/rna/data/RNA_DataInfo.fwd.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structure that holds experimental information from RNA structure mapping experiments.
//
// Originally held information in 'processed' form -- not the reactivities themselves, but
//  a guess as to whether residues (or edges) were buried or not.
//
// In 2014, decided to start including actual experimental values and implementing score functions
//  that compute log-odds scores for structural features that correlate with the data.
//
// In the future, should probably be able to deprecate backbone_burial_, backbone_exposed_ (FArrays! Dang that's old)
//  and perhaps RNA_Data, which requires user to dial in weight for exposing base edges -- those were
//  totally experimental (but nicely generic).


namespace core {
namespace scoring {
namespace rna {
namespace data {

class RNA_Datum {
public:

	RNA_Datum( Size const position, Size const edge, Real const weight ){
		position_ = position;
		edge_ = edge;
		weight_ = weight;
	}

	Size const & position() const { return position_; };
	Size const & edge() const { return edge_; };
	Real const & weight() const { return weight_; };

	void set_position( Size const setting ) { position_ = setting; };
	void set_edge( Size const setting ) { edge_ = setting; };
	void set_weight( Real const setting ) { weight_ = setting; };

private:
	Size position_;
	Size edge_;
	Real weight_;
};

typedef utility::vector1< RNA_Datum > RNA_Data;

class RNA_Reactivity {
public:

	RNA_Reactivity( Size const position, RNA_ReactivityType type, Real const value ){
		position_ = position;
		type_ = type;
		value_ = value;
	}

	Size position() const { return position_; };
	RNA_ReactivityType type() const { return type_; };
	Real value() const { return value_; };

	void set_position( 	Size const setting ) { position_ = setting; }
	void set_type( RNA_ReactivityType const setting ) { type_ = setting; }
	void set_value( Real const setting ) { value_ = setting; }

private:
	Size position_;
	RNA_ReactivityType type_;
	Real value_;
};

typedef utility::vector1< RNA_Reactivity > RNA_Reactivities;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA chemical mapping data inside the pose.
class RNA_DataInfo : public utility::pointer::ReferenceCount {

public:

	RNA_DataInfo(){}

  RNA_DataInfo( RNA_DataInfo const & src );

	~RNA_DataInfo(){}

  RNA_DataInfoOP
  clone() const
  {
    return RNA_DataInfoOP( new RNA_DataInfo( *this ) );
  }


	RNA_DataInfo &
	operator = ( RNA_DataInfo const & src );

  Size
  size() const {
    return rna_data_.size();
  }

  void
  zero();

	RNA_Data const & rna_data() const { return rna_data_; }

	void
	add_datum( RNA_Datum const & rna_datum ){ rna_data_.push_back( rna_datum ); }

	RNA_Reactivities const & rna_reactivities() const { return rna_reactivities_; }

	void
	add_reactivity( RNA_Reactivity const & rna_reactivity ){ rna_reactivities_.push_back( rna_reactivity ); }

	ObjexxFCL::FArray1D < bool > const & backbone_burial() const { return backbone_burial_; }

	void set_backbone_burial( ObjexxFCL::FArray1D < bool > const & backbone_burial ){ backbone_burial_ = backbone_burial; }

	ObjexxFCL::FArray1D < bool > const & backbone_exposed() const { return backbone_exposed_; }

	void set_backbone_exposed( ObjexxFCL::FArray1D < bool > const & backbone_exposed ){ backbone_exposed_ = backbone_exposed; }

private:

	RNA_Data rna_data_;
	ObjexxFCL::FArray1D < bool > backbone_burial_;
	ObjexxFCL::FArray1D < bool > backbone_exposed_;

	RNA_Reactivities rna_reactivities_;

};

} //data
} //rna
} //scoring
} //core

#endif
