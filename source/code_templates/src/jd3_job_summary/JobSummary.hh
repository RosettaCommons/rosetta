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

#ifndef INCLUDED_--path_underscore--_--class--_HH
#define INCLUDED_--path_underscore--_--class--_HH

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <protocols/jd3/JobSummary.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


--namespace--

///@brief --brief--
class --class--: public protocols::jd3::JobSummary{

public:

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

--end_namespace--

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( --namespace_underscore--_--class-- )
#endif // SERIALIZATION

#endif //--path_underscore--_--class--_HH

