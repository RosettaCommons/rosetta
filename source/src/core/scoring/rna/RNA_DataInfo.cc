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
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_DataInfo.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_DataInfo.fwd.hh>

// Package headers
// AUTO-REMOVED #include <core/chemical/rna/RNA_Util.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility headers
// AUTO-REMOVED #include <numeric/xyzMatrix.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>


// C++

///////////////////////////////////////////////////////
// Keep track of information from, e.g., chemical
// accessibility experiments -- useful for scoring.
///////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {


/// @details Copy constructors must copy all data, not just some...
RNA_DataInfo::RNA_DataInfo( RNA_DataInfo const & src ) :
	CacheableData()
{
	*this = src;
}

RNA_DataInfo &
RNA_DataInfo::operator=( RNA_DataInfo const & src ){

	rna_data_ = src.rna_data_;
	backbone_burial_ = src.backbone_burial_;
	backbone_exposed_ = src.backbone_exposed_;
	return  *this;

}

////////////////////////////////////////////////////////
void
RNA_DataInfo::zero()
{
	rna_data_.clear();
}


}
}
}
