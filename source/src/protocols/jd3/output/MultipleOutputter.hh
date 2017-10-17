// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/MultipleOutputter.hh
/// @brief  Definition of the %MultipleOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_output_MultipleOutputter_HH
#define INCLUDED_protocols_jd3_output_MultipleOutputter_HH

// Unit headers
#include <protocols/jd3/output/MultipleOutputter.fwd.hh>

// Package headers
#include <protocols/jd3/output/ResultOutputter.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace output {

/// @brief The %MultipleOutputter class is a vector of ResultOutputters. It should be
/// used alongside the MultipleOutputSpecification class, which is a vector of ResultSpecifications
/// and will hand each ResultSpecification in that vector to the corresponding ResultOutputter
/// in its own vector.
class MultipleOutputter : public ResultOutputter
{
public:

	MultipleOutputter();
	virtual ~MultipleOutputter();

	/// @brief Invoke write_output on all of the ResultOutputters this %MultipleOutputter contains.
	/// This class expects the OutputSpecification to be of type MultipleOutputSpecification and will
	/// hand the ResultSpecifications that this class contains to the corresponding ResultOutputter.
	void write_output(
		OutputSpecification const & specification,
		JobResult const & result
	) override;

	/// @brief Invoke flush on all of the ResultOutputters this %MultipleOutputter contains
	void flush() override;

	void append_outputter( ResultOutputterOP outputter );
	utility::vector1< ResultOutputterOP > const & outputters() const;

private:
	utility::vector1< ResultOutputterOP > outputters_;
};

} // namespace output
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_output_MultipleOutputter_HH
