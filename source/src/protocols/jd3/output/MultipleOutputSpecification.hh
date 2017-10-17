// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/MultipleOutputSpecification.hh
/// @brief  Definition of the %MultipleOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_output_MultipleOutputSpecification_HH
#define INCLUDED_protocols_jd3_output_MultipleOutputSpecification_HH

// Unit headers
#include <protocols/jd3/output/MultipleOutputSpecification.fwd.hh>

// Package headers
#include <protocols/jd3/output/OutputSpecification.hh>

// Utility headers
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace output {

/// @brief The %MultipleOutputSpecification
class MultipleOutputSpecification : public OutputSpecification
{
public:

	MultipleOutputSpecification();
	MultipleOutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index );
	virtual ~MultipleOutputSpecification();

	void append_specification( OutputSpecificationOP spec );

	utility::vector1< OutputSpecificationOP > const &
	output_specifications() const;

	/// @brief Set the result_id for myself and for the specifications I hold
	void result_id( JobResultID const & setting ) override;

	/// @brief Set the output_index for myself and for the specifications I hold
	void output_index( JobOutputIndex const & setting ) override;

	/// @brief Set the jd_output_suffix for myself and for the specifications I hold
	void jd_output_suffix( std::string const & setting ) override;

private:
	utility::vector1< OutputSpecificationOP > output_specifications_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace output
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_output_MultipleOutputSpecification )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_output_MultipleOutputSpecification_HH
