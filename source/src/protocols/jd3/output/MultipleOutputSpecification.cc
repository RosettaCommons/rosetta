// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/MultipleOutputSpecification.cc
/// @brief  Method Definitions for the %MultipleOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/output/MultipleOutputSpecification.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace output {

MultipleOutputSpecification::MultipleOutputSpecification() : OutputSpecification() {}

MultipleOutputSpecification::MultipleOutputSpecification(
	JobResultID const & result_id,
	JobOutputIndex const & output_index
) :
	OutputSpecification( result_id, output_index )
{}

MultipleOutputSpecification::~MultipleOutputSpecification() = default;

void
MultipleOutputSpecification::append_specification( OutputSpecificationOP spec )
{
	output_specifications_.push_back( spec );
}

utility::vector1< OutputSpecificationOP > const &
MultipleOutputSpecification::output_specifications() const
{
	return output_specifications_;
}

void MultipleOutputSpecification::result_id( JobResultID const & setting )
{
	OutputSpecification::result_id( setting );
	for ( auto spec : output_specifications_ ) {
		spec->result_id( setting );
	}
}

/// @brief Set the output_index for myself and for the specifications I hold
void
MultipleOutputSpecification::output_index( JobOutputIndex const & setting )
{
	OutputSpecification::output_index( setting );
	for ( auto spec : output_specifications_ ) {
		spec->output_index( setting );
	}
}

/// @brief Set the jd_output_suffix for myself and for the specifications I hold
void
MultipleOutputSpecification::jd_output_suffix( std::string const & setting )
{
	OutputSpecification::jd_output_suffix( setting );
	for ( auto spec : output_specifications_ ) {
		spec->jd_output_suffix( setting );
	}
}

} // namespace output
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::output::MultipleOutputSpecification::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::output::OutputSpecification >( this ) );
	arc( CEREAL_NVP( output_specifications_ ) ); // utility::vector1<OutputSpecificationOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::output::MultipleOutputSpecification::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::output::OutputSpecification >( this ) );
	arc( output_specifications_ ); // utility::vector1<OutputSpecificationOP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::output::MultipleOutputSpecification );
CEREAL_REGISTER_TYPE( protocols::jd3::output::MultipleOutputSpecification )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_output_MultipleOutputSpecification )
#endif // SERIALIZATION
