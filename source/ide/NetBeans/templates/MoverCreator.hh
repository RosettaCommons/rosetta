// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file /%<%NAME%>%.hh
/// @brief 
/// @author %<%USER%>% ()

#ifndef INCLUDED_ _%<%NAME%>%_hh
#define	INCLUDED_ _%<%NAME%>%_hh

#include <protocols/moves/MoverCreator.hh>
		
class %<%NAME%>% : public protocols::moves::MoverCreator {
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
	
	
};	
