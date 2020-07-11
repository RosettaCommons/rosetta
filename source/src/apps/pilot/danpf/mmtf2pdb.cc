// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/danpf/mmtf2pdb.cc
/// @brief mmtf to pdb file using rosetta
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
#include <utility/io/ozstream.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>

static basic::Tracer TR("mmtf2pdb");


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		//Enables quick one of running as same behavior as rosetta_scripts application.
		utility::vector1< utility::file::FileName > mmtfs(option[ in::file::s ]());
		for ( auto const & filename : mmtfs ) {
			core::io::StructFileReaderOptionsCOP opts( utility::pointer::make_shared<core::io::StructFileReaderOptions const >() );;
			core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( filename, *opts ));
			std::string data( core::io::pdb::create_pdb_contents_from_sfr(*mmtfsfr, opts) );
			std::string pdb_fn(filename.base() + ".pdb");
			utility::io::ozstream out(pdb_fn.c_str(), std::ios::out | std::ios::binary);
			out.write( data.c_str(), data.size() );
			out.close();
		}
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
