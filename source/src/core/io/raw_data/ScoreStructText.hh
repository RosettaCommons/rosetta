// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/raw_data/ScoreStructText.hh
///
/// @brief Write out only the scores (equivalent to a scorefile)
/// @author James Thompson, Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_ScoreStructText_hh
#define INCLUDED_core_io_raw_data_ScoreStructText_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/io/raw_data/RawStruct.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace io {
namespace raw_data {

/////////////////////////////////////////////////////////////////////////////
// holds all the data for a single entry in a silent file
class ScoreStructText : public RawStruct {

public:
	// constructor
	ScoreStructText();

	ScoreStructText(
		core::pose::Pose, // pose,
		std::string tag = "empty_tag"
	);

	// destructor
	~ScoreStructText() {}

	//ScoreStructText & operator= (ScoreStructText const & src);


	/// @brief Fill a Pose with the conformation information in this RawStruct and the FA_STANDARD
	/// ResidueTypeSet. This is a virtual method which must be implemented by classes derived from RawStruct.
	void fill_pose(
		core::pose::Pose & //pose
	);

	/// @brief Fill a Pose with the conformation information in this RawStruct and the ResidueTypeSet
	/// provided by the caller. This is a virtual method which must be implemented by classes derived from RawStruct.
	virtual void fill_pose(
		core::pose::Pose &, //pose,
		core::chemical::ResidueTypeSet const& residue_set
	);

	virtual Real get_debug_rmsd();
	virtual void print_conformation( std::ostream& out ) const;

}; // class ScoreStructText

} // namespace silent
} // namespace io
} // namespace core

#endif
