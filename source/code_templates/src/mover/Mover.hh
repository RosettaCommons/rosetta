// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


--namespace--

///@brief --brief--
class --class-- : public protocols::moves::Mover {

public:

	--class--();

	// copy constructor
	--class--( --class-- const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~--class--();

	virtual void
	apply( core::pose::Pose & pose );


public:
	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );



	//--class-- & operator=( --class-- const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;


private:



};

std::ostream &operator<< (std::ostream &os, --class-- const &mover);


--end_namespace--


#endif //--path--_--class--_hh







