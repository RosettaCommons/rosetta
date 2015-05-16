// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/AntibodyInfoRMLoader.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//unit headers
#include <protocols/antibody/AntibodyInfoRMLoader.hh>
#include <protocols/antibody/AntibodyInfoRMLoaderCreator.hh>

//package headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyInfoRMOptions.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// numeric headers
#include <numeric/xyzVector.hh>

//C++ headers
#include <istream>

//external headers
#include <boost/lexical_cast.hpp>

namespace protocols {
namespace antibody {

AntibodyInfoRMLoader::AntibodyInfoRMLoader() {}
AntibodyInfoRMLoader::~AntibodyInfoRMLoader() {}

/// @details Takes a locator id and istream, ensures that the correct number of points are present and performs error
/// checking that points are of the correct type before using the points to construct SurfaceParameters
/// @throws EXCN_Msg_Exception
utility::pointer::ReferenceCountOP
AntibodyInfoRMLoader::create_resource(
	basic::resource_manager::ResourceOptions const &,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream & 
) const
{
	using boost::lexical_cast;
	using core::Real;
	using core::Size;
	using numeric::xyzVector;
	using std::getline;
	using std::string;
	using utility::excn::EXCN_Msg_Exception;
	using utility::split;
	using utility::vector1;
	
	Size const number_of_points(3);
	Size const number_of_dimensions(3);
	
	Size lines_read(0);
	vector1< xyzVector<Real> > surf_coords(number_of_points);

	string line;
	while (getline( istream, line ))
	{
		++lines_read;
		if (lines_read > number_of_points)
		{
			break;
		}
		
		vector1< string > point_coords ( split( line ) );
		if (point_coords.size() != number_of_dimensions)
		{
			std::ostringstream err;
			err << "The input coordinates must specify a point in three-dimensional space, but the point on line ";
			err << lines_read << " of " << locator_id << " has " << point_coords.size() << " coordinates listed.";
			throw EXCN_Msg_Exception(err.str());
		}
		surf_coords[lines_read] = xyzVector<Real>(lexical_cast<Real>(point_coords[1]),lexical_cast<Real>(point_coords[2]),
			lexical_cast<Real>(point_coords[3]));
	}
	
	if (lines_read > number_of_points)
	{
		throw EXCN_Msg_Exception( "AntibodyInfoRMLoader expected to be given exactly three points to " \
			"define the periodicity of the surface, but more than three points were provided in "  + locator_id + ".");
	}
	
	else if (lines_read < number_of_points)
	{
		throw EXCN_Msg_Exception( "AntibodyInfoRMLoader expected to be given exactly three points to " \
			"define the periodicity of the surface, but less than three points were provided in "  + locator_id + ".");
	}
	
	return utility::pointer::ReferenceCountOP( new SurfaceParameters(surf_coords[1], surf_coords[2], surf_coords[3]) );
}

basic::resource_manager::ResourceOptionsOP
AntibodyInfoRMLoader::default_options() const
{
	return basic::resource_manager::ResourceOptionsOP( new SurfaceVectorOptions() );
}

basic::resource_manager::ResourceLoaderOP AntibodyInfoRMLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new AntibodyInfoRMLoader() );
}

std::string AntibodyInfoRMLoaderCreator::loader_type() const
{
	return "SurfaceVector";
}

} // namespace surface_docking
} // namespace protocols
