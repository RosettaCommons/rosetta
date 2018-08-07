// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

#include <--path--/--class--.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

--namespace--

/// @brief --brief--
class --class-- : public utility::pointer::ReferenceCount {

public:

	--class--();
	--class--(--class-- const & src);

	virtual ~--class--();

	--class--OP
	clone() const;

private:

};


--end_namespace--



#endif //INCLUDED_--path_underscore--_--class--_hh





