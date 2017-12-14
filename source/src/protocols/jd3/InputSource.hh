// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/InputSource.hh
/// @brief  Declaration of the %InputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_InputSource_hh
#define INCLUDED_protocols_jd3_InputSource_hh

//unit headers
#include <protocols/jd3/InputSource.fwd.hh>

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

/// @brief The %InputSource is a small class for holding data about
/// the starting Pose for a Job and where it comes from (i.e. which
/// of the Inputters claims responsibility for creating a Pose for this
/// instance). The "input_tag" is a string description of the input source and will
/// be used as the "job_tag" to control output -- the input tag should not
/// include the file extension.  It is perfectly reasonable for complex
/// Inputters to subclass from InputSource to tuck more complex
/// data in the InputSource, though, the string-string map ought to
/// provide considerable flexibility in storing data without deriving
/// new subclasses.
class InputSource : public utility::pointer::ReferenceCount
{
public:
	InputSource();
	InputSource( std::string const & origin );
	~InputSource() override;

	bool operator == ( InputSource const & rhs ) const;
	bool operator != ( InputSource const & rhs ) const;
	bool operator <  ( InputSource const & rhs ) const;

	std::string const & input_tag() const;
	//StringStringMap const & string_string_map() const = 0;
	std::string const & origin() const;
	core::Size pose_id() const;

	void input_tag( std::string const & setting );
	//virtual void store_string_pair( std::string const & key, std::string const & value ) = 0;
	void origin( std::string const & setting );
	void pose_id( core::Size setting );

private:
	std::string origin_;
	std::string input_tag_;
	core::Size pose_id_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_InputSource )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_InputSource_HH
