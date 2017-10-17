// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/PoseOutputSpecification.hh
/// @brief  Definition of the %PoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_standard_PoseOutputSpecification_HH
#define INCLUDED_protocols_jd3_standard_PoseOutputSpecification_HH

// Unit headers
#include <protocols/jd3/standard/PoseOutputSpecification.fwd.hh>

// Package headers
#include <protocols/jd3/JobResultID.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace standard {

/// @brief The %PoseOutputSpecification
class PoseOutputSpecification : utility::pointer::ReferenceCount
{
public:

	PoseOutputSpecification();
	PoseOutputSpecification( JobResultID const & result_id );
	virtual ~PoseOutputSpecification();

	std::string const & outputter_type() const;
	void outputter_type( std::string const & setting );

private:
	std::string outputter_type_;

};

} // namespace standard
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_standard_PoseOutputSpecification_HH
