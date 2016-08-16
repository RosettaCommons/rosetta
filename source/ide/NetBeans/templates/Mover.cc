// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file /%<%NAME%>%.cc
/// @brief 
/// @author %<%USER%>% ()

#include < /%<%NAME%>%.hh>
#include < /%<%NAME%>%Creator.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR("%<%NAME%>%");

%<%NAME%>%::%<%NAME%>%():
		protocols::moves::Mover("%<%NAME%>%")
{
	
}

%<%NAME%>%::~%<%NAME%>%(){}
		
%<%NAME%>%::%<%NAME%>%(%<%NAME%>% const & src):
		protocols::moves::Mover(src)
{
	
}

void
%<%NAME%>%::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& ,
	const Filters_map& ,
	const Movers_map& ,
	const Pose& )
{

}

protocols::moves::MoverOP
%<%NAME%>%::clone() const{
	return protocols::moves::MoverOP( new %<%NAME%>%(*this) );
}

%<%NAME%>% & operator=( %<%NAME%>% const & src){
	return %<%NAME%>%(src);
}

moves::MoverOP
%<%NAME%>%::fresh_instance() const
{
	return protocols::moves::MoverOP( new %<%NAME%>% );
}

std::string
%<%NAME%>%::get_name() const {
	return "%<%NAME%>%";
}





void
%<%NAME%>%::apply(core::pose::Pose& pose){
	
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
%<%NAME%>%Creator::create_mover() const {
	return new %<%NAME%>%;
}

std::string
%<%NAME%>%Creator::keyname() const {
	return %<%NAME%>%Creator::mover_name();
}

std::string
%<%NAME%>%Creator::mover_name(){
	return "%<%NAME%>%";
}