// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file --path--/--class--
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

#include <--path--/--class--.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>


--namespace--

///@brief --brief--
class --class--: public core::pack::task::operation::TaskOperation {
public:

	--class--();

	--class--(--class-- const & src);

	virtual ~--class--();

	core::pack::task::operation::TaskOperationOP
	clone() const;

	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );


	//////////////////////

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;



private:


};



--end_namespace--


#endif //INCLUDED_--path--_--class--_fwd_hh

