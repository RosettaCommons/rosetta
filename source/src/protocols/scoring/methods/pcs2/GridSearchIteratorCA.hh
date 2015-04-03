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
 /// @file protocols/scoring/methods/pcs2/GridSearchIteratorCA.hh
 ///
 /// @brief
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


#ifndef INCLUDED_protocols_scoring_methods_pcs2_GridSearchIteratorCA_hh
#define INCLUDED_protocols_scoring_methods_pcs2_GridSearchIteratorCA_hh


// Unit headers

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

//Auto Headers
// Numeric headers
// Objexx headers

// C++ headers


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

class GridSearchIteratorCA{
public:

	GridSearchIteratorCA(); //Construct

	~GridSearchIteratorCA(); //Destruct

	GridSearchIteratorCA(GridSearchIteratorCA const & other); //copy

	GridSearchIteratorCA & // =
	operator=(GridSearchIteratorCA const & other);

	GridSearchIteratorCA(core::pose::Pose const & pose);

	/// @brief give me the next x-y-z coordinate to visit
	/// bool return FALSE if everything has been visited
	bool
	next_center(core::Real &x,
							core::Real &y,
							core::Real &z);

	void
	reset();


private:
	void
	set_vec(utility::vector1<core::Real> & x_vec,
					utility::vector1<core::Real> & y_vec,
					utility::vector1<core::Real> & z_vec,
					core::Size index,
					core::pose::Pose const & pose);

	utility::vector1<core::Real>  x_vec_;
	utility::vector1<core::Real>  y_vec_;
	utility::vector1<core::Real>  z_vec_;
	core::Size res_num_cur_;
	core::Size res_num_total_;
};

} //namespace pcs2
} //namespace methods
} //namespace scoring
} //namespace methods

#endif
