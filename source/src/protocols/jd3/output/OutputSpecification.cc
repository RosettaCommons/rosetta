// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/OutputSpecification.cc
/// @brief  Method definitions of the %OutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/output/OutputSpecification.hh>
#include <utility>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace output {

OutputSpecification::OutputSpecification() = default;

OutputSpecification::OutputSpecification(
	JobResultID const & result_id,
	JobOutputIndex const & output_index
) :
	result_id_( result_id ),
	output_index_( output_index )
{}

OutputSpecification::~OutputSpecification() = default;

JobResultID
OutputSpecification::result_id() const
{
	return result_id_;
}

void OutputSpecification::result_id( JobResultID const & setting )
{
	result_id_ = setting;
}

JobOutputIndex
OutputSpecification::output_index() const
{
	return output_index_;
}

void
OutputSpecification::output_index( JobOutputIndex const & setting )
{
	output_index_ = setting;
}

std::string const &
OutputSpecification::jd_output_suffix() const
{
	return jd_output_suffix_;
}

void
OutputSpecification::jd_output_suffix( std::string const & setting )
{
	jd_output_suffix_ = setting;
}

std::string
OutputSpecification::suffix_from_jd_with_sep() const
{
	return jd_output_suffix().empty() ? "" : ( "." + jd_output_suffix() );
}

} // namespace output
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::output::OutputSpecification::save( Archive & arc ) const {
	arc( CEREAL_NVP( result_id_ ) ); // JobResultID
	arc( CEREAL_NVP( output_index_ ) ); // struct protocols::jd3::JobOutputIndex
	arc( CEREAL_NVP( jd_output_suffix_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::output::OutputSpecification::load( Archive & arc ) {
	arc( result_id_ ); // JobResultID
	arc( output_index_ ); // struct protocols::jd3::JobOutputIndex
	arc( jd_output_suffix_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::output::OutputSpecification );
CEREAL_REGISTER_TYPE( protocols::jd3::output::OutputSpecification )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_output_OutputSpecification )
#endif // SERIALIZATION
