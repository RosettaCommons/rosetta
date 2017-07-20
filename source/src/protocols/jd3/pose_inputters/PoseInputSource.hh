// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PoseInputSource.hh
/// @brief  Declaration of the %PoseInputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PoseInputSource_hh
#define INCLUDED_protocols_jd3_pose_inputters_PoseInputSource_hh

//unit headers
#include <protocols/jd3/pose_inputters/PoseInputSource.fwd.hh>
#include <protocols/jd3/InputSource.hh>

//project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <string>

// Utility headers
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace pose_inputters {

/// @brief The %PoseInputSource is a small class for holding data about
/// the starting Pose for a Job and where it comes from (i.e. which
/// of the PoseInputters claims responsibility for creating a Pose for this
/// instance). The "input_tag" is a string description of the input source and will
/// be used as the "job_tag" to control output -- the input tag should not
/// include the file extension.  It is perfectly reasonable for complex
/// PoseInputters to subclass from PoseInputSource to tuck more complex
/// data in the PoseInputSource, though, the string-string map ought to
/// provide considerable flexibility in storing data without deriving
/// new subclasses.
class PoseInputSource : public InputSource
{
public:
	typedef std::map< std::string, std::string > StringStringMap;

	PoseInputSource() :
		InputSource()
	{}
	PoseInputSource( std::string origin ); // move-constructed

public:
	bool operator == ( PoseInputSource const & rhs ) const;
	bool operator != ( PoseInputSource const & rhs ) const;
	bool operator <  ( PoseInputSource const & rhs ) const;

	StringStringMap const & string_string_map() const;
	void store_string_pair( std::string const & key, std::string const & value );

private:
	std::string origin_;
	std::string input_tag_;
	StringStringMap string_string_map_;
	core::Size pose_id_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_pose_inputters_PoseInputSource )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_pose_inputters_PoseInputSource_HH
