// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file --path--/--class--Creator.hh
/// @brief --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--Creator_hh
#define INCLUDED_--path_underscore--_--class--Creator_hh

#include <protocols/filters/FilterCreator.hh>

--namespace--

class --class--Creator : public protocols::filters::FilterCreator {
public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
	static std::string filter_name();
};

--end_namespace--

#endif //INCLUDED_--path_underscore--_--class--_fwd_hh



