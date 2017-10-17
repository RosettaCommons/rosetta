// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.hh
/// @brief  Definition of the %SilentFilePoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_SilentFilePoseOutputSpecification_HH
#define INCLUDED_protocols_jd3_pose_outputters_SilentFilePoseOutputSpecification_HH

// Unit headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.fwd.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>

// Package headers
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/JobOutputIndex.hh>

// Project headers
#include <core/io/silent/SilentFileOptions.fwd.hh>

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

/// @brief The %SilentFilePoseOutputSpecification giving all of the details of how
/// to write a Pose out to disk as a silent file
class SilentFilePoseOutputSpecification : public PoseOutputSpecification
{
public:

	SilentFilePoseOutputSpecification();
	SilentFilePoseOutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index );
	virtual ~SilentFilePoseOutputSpecification();

	core::io::silent::SilentFileOptions const & sf_opts() const;
	void sf_opts( core::io::silent::SilentFileOptions const & setting );

	std::string out_fname() const;
	void out_fname( std::string const & setting );

	std::string const & pose_tag() const;
	void pose_tag( std::string const & setting );

	core::Size buffer_limit() const;
	void buffer_limit( core::Size setting );

private:

	core::io::silent::SilentFileOptionsOP sf_opts_;
	std::string out_fname_;
	std::string pose_tag_;
	core::Size buffer_limit_;
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
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_pose_outputters_SilentFilePoseOutputSpecification )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_pose_outputters_SilentFilePoseOutputSpecification_HH
