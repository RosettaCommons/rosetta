// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/ScoreStructJSON.hh
///
/// @brief Write out only the scores in JSON format
/// @author Luki Goldschmidt

#ifndef INCLUDED_core_io_raw_data_ScoreStructJSON_hh
#define INCLUDED_core_io_raw_data_ScoreStructJSON_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/io/raw_data/RawStruct.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace io {
namespace raw_data {

class ScoreStructJSON : public RawStruct {

public:
	ScoreStructJSON();
	ScoreStructJSON( core::pose::Pose, std::string tag );

	// Not implemented for this format  
	void fill_pose( core::pose::Pose & ) {}
	void fill_pose( core::pose::Pose &, core::chemical::ResidueTypeSet const & ) {}
	void print_conformation( std::ostream & ) const {}
	Real get_debug_rmsd() { return 0; }

	void print_header(
		std::ostream& /* out */,
		std::map < std::string, core::Real > const & /* score_map */,
		std::map < std::string, std::string > const & /* string_map */,
		bool /* print_sequence */ ) const
	{
		// No header in JSON
	}
		
	void print_scores(
		std::ostream& out,
		std::map < std::string, core::Real > const & score_map,
		std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () ) ) const;

}; // class ScoreStructJSON

} // namespace silent
} // namespace io
} // namespace core

#endif
