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
 /// @file protocols/scoring/methods/pcs2/PcsInputCenter.hh
 ///
 /// @brief This class hold the information for each paramagnetic center
 /// Multiple data set with different lanthanide can be measured for each paramagnetic center
 ///
 /// @details
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorv Christophe Schmitz
 ///
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsInputCenter_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsInputCenter_hh

// Package headers
#include  <protocols/scoring/methods/pcs2/PcsInputFile.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <map>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

//////////////////////////////////////////////////////////////
/// @brief PcsInputCenter contain all the input information for one paramagnetic center.
/// It can contain multiple data set
class PcsInputCenter : public utility::pointer::ReferenceCount{
private:
	std::map< std::string, PcsInputFile > PcsInputFile_all_;

public:
	PcsInputCenter(); //Construct

	virtual ~PcsInputCenter(); //Destruct

	PcsInputCenter(PcsInputCenter const & other); //Copy

	PcsInputCenter &
	operator=( PcsInputCenter const & other ); //=

	PcsInputCenter(utility::vector1<std::string> const & filenames,  utility::vector1<core::Real> const & weight);

	/// @brief Give me all PcsInputFile
	std::map< std::string, PcsInputFile > &
	get_PcsInputFile_all();

	/// @brief Print me
	friend std::ostream &
	operator<<(std::ostream& out, const PcsInputCenter &me);

};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
