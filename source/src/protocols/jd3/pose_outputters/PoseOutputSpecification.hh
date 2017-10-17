// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PoseOutputSpecification.hh
/// @brief  Definition of the %PoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_PoseOutputSpecification_HH
#define INCLUDED_protocols_jd3_pose_outputters_PoseOutputSpecification_HH

// Unit headers
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.fwd.hh>

// Package headers
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/JobOutputIndex.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief The %PoseOutputSpecification
class PoseOutputSpecification : public output::OutputSpecification
{
public:

	PoseOutputSpecification();
	PoseOutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index );
	virtual ~PoseOutputSpecification();

	std::string const & outputter_type() const;
	void outputter_type( std::string const & setting );

private:
	std::string outputter_type_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_pose_outputters_PoseOutputSpecification )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_pose_outputters_PoseOutputSpecification_HH
