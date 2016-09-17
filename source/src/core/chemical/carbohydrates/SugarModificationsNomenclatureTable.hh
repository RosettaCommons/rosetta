// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/SugarModificationsNomenclatureTable.hh
/// @brief   Definitions for SugarModificationsNomenclatureTable.
/// @author  Labonte <JWLabonte@jhu.edu


#ifndef INCLUDED_core_chemical_carbohydrates_SugarModificationsNomenclatureTable_HH
#define INCLUDED_core_chemical_carbohydrates_SugarModificationsNomenclatureTable_HH

// Project header
#include <core/types.hh>

// C++ headers
#include <string>
#include <map>


namespace core {
namespace chemical {
namespace carbohydrates {

/// @brief  A structure for storing information related to the nomenclature of modified sugars.
struct SugarModificationsNomenclatureTableRow {
	std::string substituent_full_name;  // e.g., "acetylamino"; used for saccharide IUPAC names
	std::string implies;  // usually, "deoxy"
	std::string short_affix;  // e.g., "NAc"; used for saccharide abbreviations
	std::string patch_name;  // e.g., "AcNH"; Rosetta's patch name for this modification
	core::uint default_position;  // e.g., "2"; the position assumed if absent from the short affix
	bool has_inherent_position;  //  e.g., "0" = false, "1" = true; the position is defined by the modification
	//std::string reducing_end_suffix;  // e.g., "ose"; used for saccharide IUPAC names
	//std::string glycoside_suffix;  // e.g., "oside"; used for saccharide IUPAC names
	//std::string following_word_or_phrase;  // e.g., "sulfate"; used for saccharide IUPAC names
};  // struct SugarModificationsNomenclatureTableRow


typedef std::map< std::string, SugarModificationsNomenclatureTableRow > SugarModificationsNomenclatureTable;

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_SugarModificationsNomenclatureTable_HH
