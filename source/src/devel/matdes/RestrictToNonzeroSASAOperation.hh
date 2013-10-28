// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/RestrictToNonzeroSASAOperation.hh
/// @brief  Restrict design to only residues that have SASA>0 in the monomeric state
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_devel_matdes_RestrictToNonzeroSASAOperation_hh
#define INCLUDED_devel_matdes_RestrictToNonzeroSASAOperation_hh

// Unit Headers
#include <devel/matdes/RestrictToNonzeroSASAOperation.fwd.hh>
#include <devel/matdes/RestrictToNonzeroSASAOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/types.hh>

// Utility Headers

// C++ Headers

namespace devel { 
namespace matdes {

class RestrictToNonzeroSASAOperation : public core::pack::task::operation::TaskOperation {
public:
	RestrictToNonzeroSASAOperation( core::Real probe_radius = 2.2 , core::Size ncomp = 1 );

	virtual ~RestrictToNonzeroSASAOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	void parse_tag( TagCOP tag , DataMap & );
	void parse_def( utility::lua::LuaObject const & def );

private:

 	core::Real probe_radius_; 
	core::Size ncomp_;

};

} //namespace matdes 
} //namespace devel

#endif 
