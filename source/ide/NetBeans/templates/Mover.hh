// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file .../%<%NAME%>%.hh
/// @brief
/// @author %<%USER%>% ()


#ifndef INCLUDED_ _%<%NAME%>%_hh
#define INCLUDED_ _%<%NAME%>%_hh

#include < /%<%NAME%>%.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


// Forward
class %<%NAME%>% : public protocols::moves::Mover {
	
public :
	
	%<%NAME%>%();
	%<%NAME%>%(%<%NAME%>% const & src);
	
	virtual~%<%NAME%>%();
	
	virtual void
	apply(core::pose::Pose & pose);
	
	
public:
	
	std::string
	get_name() const;
	
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose
	);	
	
	protocols::moves::MoverOP
	clone() const;
	
	%<%NAME%>% & operator=( %<%NAME%>% const & src);
	
	virtual moves::MoverOP fresh_instance() const;
	
private:
	

	
};





#endif	//INCLUDED_ %<%CLASSNAME%>%.hh







