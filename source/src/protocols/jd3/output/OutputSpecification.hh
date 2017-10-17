// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/OutputSpecification.hh
/// @brief  Definition of the %OutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_output_OutputSpecification_HH
#define INCLUDED_protocols_jd3_output_OutputSpecification_HH

// Unit headers
#include <protocols/jd3/output/OutputSpecification.fwd.hh>

// Package headers
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/CompletedJobOutput.fwd.hh>

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
namespace output {

/// @brief The %OutputSpecification
class OutputSpecification : utility::pointer::ReferenceCount
{
public:

	OutputSpecification();
	OutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index );
	virtual ~OutputSpecification();

	JobResultID result_id() const;
	virtual void result_id( JobResultID const & setting );

	JobOutputIndex output_index() const;
	virtual void output_index( JobOutputIndex const & setting );

	std::string const & jd_output_suffix() const;
	virtual void jd_output_suffix( std::string const & setting );

	// a dot separator character for an output filename if the JD has provided a suffix
	std::string suffix_from_jd_with_sep() const;

private:
	JobResultID    result_id_;
	JobOutputIndex output_index_;
	std::string jd_output_suffix_;
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
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_output_OutputSpecification )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_output_OutputSpecification_HH
