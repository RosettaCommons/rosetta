// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/util.cc
/// @brief  Utility functions for helical bundle construction.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// C++ headers
#include <stdio.h>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/conversions.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

//static numeric::random::RandomGenerator RG(192923);  // <- Magic number, do not change it!

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

	static basic::Tracer TR("protocols.helical_bundle.util");

	void write_minor_helix_params (
		std::string const &filename,
		utility::vector1 < core::Real > const &r1,
		core::Real const &omega1,
		core::Real const &z1,
		utility::vector1 < core::Real > const &delta_omega1,
		utility::vector1 < core::Real > const &delta_z1
	) {
		using namespace utility::io;

		ozstream outfile;
		outfile.open( filename );

		runtime_assert_string_msg( outfile.good(), "In protocols::helical_bundle::write_minor_helix_params: Unable to open file for write!" );

		char linebuffer[256];
		std::string line_out;

		for(core::Size i=1, imax=r1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "r1\t%.12f\n", r1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		sprintf(linebuffer, "omega1\t%.12f\nz1\t%.12f\n", omega1, z1);
		line_out=linebuffer;
		outfile.write(line_out, line_out.length() );

		for(core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "delta_omega1\t%.12f\n", delta_omega1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		for(core::Size i=1, imax=delta_z1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "delta_z1\t%.12f\n", delta_z1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		outfile.close();

		return;
	}

	void read_minor_helix_params (
		std::string const &filename,
		utility::vector1 < core::Real > &r1,
		core::Real &omega1,
		core::Real &z1,
		utility::vector1 < core::Real > &delta_omega1,
		utility::vector1 < core::Real > &delta_z1
	) {
		using namespace utility::io;

		r1.clear();
		delta_omega1.clear();
		delta_z1.clear();
		omega1 = 0.0;
		z1 = 0.0;

		std::string filename_formatted = filename;
		if(utility::file::file_extension(filename_formatted)!="crick_params") filename_formatted+= ".crick_params";

		izstream infile;
		infile.open( filename_formatted );
		if(!infile.good()) {
			filename_formatted = "protocol_data/crick_parameters/" + utility::file::file_basename(filename_formatted) + ".crick_params";
			basic::database::open( infile, filename_formatted );
			runtime_assert_string_msg( infile.good(), "In protocols::helical_bundle::read_minor_helix_params: Unable to open .crick_params file for read!" );
		}

		if(TR.Debug.visible()) TR.Debug << "Reading " << filename_formatted << std::endl;

		while(true) {
			std::string current_line = "";
			infile.getline(current_line);//Get the current line and store it in current_line.
			if(infile.eof()) break;

			if(TR.Debug.visible()) TR.Debug << current_line << std::endl;

			if(current_line[0] == '#') continue; //Ignore lines that start with a pound sign.

			std::stringstream ss(current_line);
			char buffer [25];
			ss.getline(buffer, 25, '\t');
			std::string strbuffer=buffer;
			//TR.Debug << strbuffer << "." << std::endl; //DELETE ME
			if(strbuffer=="r1") {
				ss.getline(buffer, 25);
				r1.push_back( atof(buffer) );
			} else if (strbuffer=="delta_omega1") {
				ss.getline(buffer, 25);
				delta_omega1.push_back( atof(buffer) );
			} else if (strbuffer=="delta_z1") {
				ss.getline(buffer, 25);
				delta_z1.push_back( atof(buffer) );
			} else if (strbuffer=="omega1") {
				ss.getline(buffer, 25);
				omega1=atof(buffer);
			} else if (strbuffer=="z1") {
				ss.getline(buffer, 25);
				z1=atof(buffer);
			} else continue;

		}

		infile.close();

		if(TR.Debug.visible()) {
			TR.Debug << "Finished reading " << filename_formatted << std::endl;
			for(core::Size i=1, imax=r1.size(); i<=imax; ++i) {
				TR.Debug << "r1[" << i << "]\t" << r1[i] << std::endl;
			}
			TR.Debug << "omega1\t" << omega1 << std::endl;
			TR.Debug << "z1\t" << z1 << std::endl;
			for(core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i) {
				TR.Debug << "delta_omega1[" << i << "]\t" << delta_omega1[i] << std::endl;
			}
			for(core::Size i=1, imax=delta_z1.size(); i<=imax; ++i) {
				TR.Debug << "delta_z1[" << i << "]\t" << delta_z1[i] << std::endl;
			}
		}

		runtime_assert_string_msg( (r1.size() == delta_omega1.size()) && ( r1.size() == delta_z1.size() ), "In protocols::helical_bundle::read_minor_helix_params: the r1, delta_omega1, and delta_z1 lists are of different lengths.  Check the .crick_params file, since something is clearly wonky." );

		return;
	}


} //helical_bundle
} //namespace protocols
