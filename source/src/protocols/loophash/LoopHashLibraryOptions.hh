// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/LoopHashLibraryOptions.hh
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoopHashLibraryOptions_hh
#define INCLUDED_protocols_loophash_LoopHashLibraryOptions_hh

// unit headers
#include <protocols/loophash/LoopHashLibraryOptions.fwd.hh>

// core headers
#include <core/types.hh>

// project headers
#include <basic/resource_manager/ResourceOptions.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ headers

namespace protocols {
namespace loophash {

/// @brief %LoopHashLibraryOptions encapsulates the options associated with LoopHashLibrary %resources.
/// @details These options are read in from a resource definition file and each LoopHashLibrary resource has a
/// corresponding %LoopHashLibraryOptions instance.
class LoopHashLibraryOptions : public basic::resource_manager::ResourceOptions
{
public:
	/// @brief Construct the %LoopHashLibraryOptions.
	LoopHashLibraryOptions();

	/// @brief Destructor.
	virtual ~LoopHashLibraryOptions();

	/// @brief Read the configuration of the LoopHashLibrary %resource from the tag generated from the resource definition
	/// file.
	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr tag
	);

	/// @brief Return the name of this class (LoopHashLibraryOptions).
	virtual
	std::string
	type() const;

	/// @brief Return a %vector of loop sizes that will be used to construct a LoopHashLibrary.
	utility::vector1<core::Size>
	loop_sizes() const;

	/// @brief Set the loop sizes that will be used to construct a LoopHashLibrary to a %vector of loop lengths.
	void
	loop_sizes(
		utility::vector1<core::Size> loop_sizes
	);

private:
	/// @brief %Vector of %Sizes that stores the loop lengths that will be used to construct a LoopHashLibrary.
	utility::vector1<core::Size> loop_sizes_;
};

} // namespace loophash
} // namespace protocols

#endif //INCLUDED_protocols_loops_LoopHashLibraryOptions_hh
