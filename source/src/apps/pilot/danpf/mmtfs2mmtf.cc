// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/danpf/mmtfs2mmtf.cc
/// @brief list of mmtfs to a single mmtf (using models)
/// @author Danny Farrell (danpf@uw.edu)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/io/StructFileReaderOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/io/mmtf/mmtf_writer.hh>
#include <core/io/mmtf/mmtf_reader.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>

static basic::Tracer TR("mmtfs2mmtf");


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		//Enables quick one of running as same behavior as rosetta_scripts application.
		std::string const out_mmtf_file(option[ out::file::o ]());
		utility::vector1< utility::file::FileName > mmtfs(option[ in::file::s ]());
		utility::vector1<core::io::StructFileRepOP> sfrs;
		for ( auto const & filename : mmtfs ) {
			TR << "getting: " << filename << std::endl;
			core::io::StructFileReaderOptionsCOP opts( utility::pointer::make_shared<core::io::StructFileReaderOptions const >() );;
			core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( filename, *opts ));
			sfrs.push_back(mmtfsfr);
		}
		core::io::StructFileRepOptions options;
		core::io::mmtf::set_mmtf_default_options(options);
		mmtf::StructureData sd(core::io::mmtf::sfrs_to_sd(sfrs, options));
		mmtf::encodeToFile(sd, out_mmtf_file);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
