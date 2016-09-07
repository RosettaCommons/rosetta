// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/AnchoredDesign/Anchor.cc
/// @brief Anchor methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/anchored_design/Anchor.hh>

// Project Headers

#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// C++ Headers
#include <string>
#include <iostream>

// option key includes
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.AnchoredDesign.Anchor" );

namespace protocols {
namespace anchored_design {

std::string const bad_anchorfile("YOU_FORGOT_TO_SPECIFY_AN_ANCHOR_FILE");

// default ctor sets all to 0.  You'll need to run the setter and run read_anchorfile yourself.
protocols::anchored_design::Anchor::Anchor():
	start_(0),
	end_(0),
	anchorfile_(bad_anchorfile)
{}

// input constructor for arbitrary anchor, if you know resids
protocols::anchored_design::Anchor::Anchor( core::Size const start, core::Size const end):
	start_(start),
	end_(end),
	anchorfile_("NO_ANCHORFILE_NEEDED")
{}

// input constructor from pose and anchorfile
protocols::anchored_design::Anchor::Anchor( core::pose::Pose const & pose, std::string const & filename):
	start_(0),
	end_(0),
	anchorfile_(filename)
{
	read_anchorfile(pose, filename);
}

// input constructor from pose, assumes anchorfile
protocols::anchored_design::Anchor::Anchor( core::pose::Pose const & pose):
	start_(0),
	end_(0),
	anchorfile_(bad_anchorfile)
{
	read_options();
	read_anchorfile(pose);
}

/// @details virtual destructors make c++ happy
protocols::anchored_design::Anchor::~Anchor()= default;

/// @brief copy ctor
Anchor::Anchor( Anchor const & rhs ) :
	utility::pointer::ReferenceCount(rhs)
{
	*this = rhs;
}

/// @brief assignment operator
Anchor & Anchor::operator=( Anchor const & rhs ){

	//abort self-assignment
	if ( this == &rhs ) return *this;

	start_ = rhs.start();
	end_ = rhs.end();
	anchorfile_ = rhs.get_filename();
	return *this;
}

/// @details read_anchorfile from internally stored anchorfile (from the option system)
void protocols::anchored_design::Anchor::read_anchorfile(
	core::pose::Pose const & pose)
{
	read_anchorfile(pose, anchorfile_);
}

/// @details read_anchorfile reads in an anchor file (formatted as whitespace-delineated PDB-relevant start end chain information).  It queries the pose('s members) for the mapping from PDB numbers to pose resids.
void protocols::anchored_design::Anchor::read_anchorfile(
	core::pose::Pose const & pose,
	std::string const & filename )
{

	//find the anchor file, open it, handle error
	utility::io::izstream anchorfile( filename );
	if ( !anchorfile ) {
		Error() << "Can't open anchorfile file: " << filename << std::endl;
		Error() << "Use -AnchoredDesign::anchor to specify" << std::endl;
		utility_exit();
	}

	core::Size PDBstart;
	core::Size PDBend;
	char chain;

	anchorfile >> chain >> PDBstart >> PDBend;

	//if the stream fails, bad formatting
	if ( anchorfile.fail() ) {
		Error() << "Can't parse anchor file.  Using the PDB numbering, format is:\n    chain start end \n  "
			<< "Only the first line is read." << std::endl;
		utility_exit();
	}
	//use info to modify Anchor's data members
	using core::pose::PDBInfo;
	using core::pose::PDBPoseMap;
	if ( !pose.pdb_info() ) { //if that's not a NULL pointer
		utility_exit_with_message("Cannot read anchor file because pose has no PDBInfo");
	}

	PDBPoseMap const & pose_map = pose.pdb_info()->pdb2pose();
	start_ = pose_map.find(chain, PDBstart);
	end_ = pose_map.find(chain, PDBend);

}


/// @brief setter for filename for anchorfile
void protocols::anchored_design::Anchor::set_filename(std::string const & filename) {anchorfile_ = filename;}

/// @brief getter for filename for anchorfile
std::string const & protocols::anchored_design::Anchor::get_filename() const {return anchorfile_;}

/// @brief read option system into internal data
void protocols::anchored_design::Anchor::read_options() {
	anchorfile_ = basic::options::option[ basic::options::OptionKeys::AnchoredDesign::anchor ].value();
}

}//namespace protocols
}//namespace anchored_design
