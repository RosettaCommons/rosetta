// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file   core/io/silent/SilentFileLoader.hh
/// @brief  Load a silent file from an input stream and a SilentFileRMOptions object
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/io/silent/SilentFileLoader.hh>
#include <core/io/silent/SilentFileLoaderCreator.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Project Headers
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileRMOptions.hh>
#include <core/pose/Pose.hh>

// Platform Headers
#include <core/types.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

static basic::Tracer TR( "core.io.silent.SilentFileLoader" );

namespace core {
namespace io {
namespace silent {

using pose::Pose;
using pose::PoseOP;
using basic::resource_manager::ResourceOP;
using basic::resource_manager::ResourceLoaderOP;
using basic::resource_manager::ResourceOptions;
using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::LocatorID;
using std::istream;
using std::string;

///// SilentFileLoaderCreator /////
SilentFileLoaderCreator::SilentFileLoaderCreator() {}

SilentFileLoaderCreator::~SilentFileLoaderCreator() {}

ResourceLoaderOP
SilentFileLoaderCreator::create_resource_loader() const {
	return ResourceLoaderOP( new SilentFileLoader );
}

string
SilentFileLoaderCreator::loader_type() const {
	return "SilentFile";
}

//// SilentFileLoader /////
SilentFileLoader::SilentFileLoader() {}

SilentFileLoader::~SilentFileLoader() {}

SilentFileLoader::SilentFileLoader( SilentFileLoader const &) : ResourceLoader() {}

ResourceOP
SilentFileLoader::create_resource(
	ResourceOptions const & options,
	LocatorID const & locator_id,
	istream & istream
) const {

	if ( ! dynamic_cast< SilentFileRMOptions const * >( &options ) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, 
			"SilentFileLoader expected to get a SilentFileRMOptions object, "
			"but was given a ResourceOptions of type '" + options.type() + "', "
			"which has the name '" + options.name() + "'." );
	}
	SilentFileRMOptions const & resource_options(
		static_cast< SilentFileRMOptions const & >( options ));

	utility::vector1< std::string > lines;
	std::string line;
	while ( getline(istream, line) ) {
		lines.push_back(line);
	}

	SilentStructOP ss( SilentStructFactory::get_instance()->get_silent_struct_in( resource_options.opts() ));

	SilentFileData container( resource_options.opts() );

	if ( !(ss->init_from_lines(lines, container)) ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "SilentFileLoader failed to load silent file with locator_id '" + locator_id + "'." );
	}

	PoseOP pose( new Pose() );
	ss->fill_pose(*pose);

	if ( pose->size() == 0 ) {
		TR.Warning
			<< "Loading Pose with SilentFileLoader for the "
			<< "locator_id '" << locator_id << "', "
			<< "but the resulting pose does not have any residues." << std::endl;
		TR.Warning
			<< "Also note that the input stream has '" << lines.size() << "'" << std::endl;
	}

	return pose;
}

ResourceOptionsOP
SilentFileLoader::default_options(
) const {
	return ResourceOptionsOP( new SilentFileRMOptions() );
}

} // namespace
} // namespace
} // namespace

