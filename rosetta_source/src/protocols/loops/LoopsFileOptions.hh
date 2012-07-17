// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopsFileOptions.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_loops_LoopsFileOptions_hh
#define INCLUDED_protocols_loops_LoopsFileOptions_hh

//unit headers
#include <protocols/loops/LoopsFileOptions.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers

namespace protocols {
namespace loops {

class LoopsFileOptions : public basic::resource_manager::ResourceOptions
{
public:
	LoopsFileOptions();
	virtual ~LoopsFileOptions();

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr tag
	);

	virtual
	std::string
	type() const;

	bool prohibit_single_residue_loops() const;
	void prohibit_single_residue_loops( bool setting );

private:
	bool prohibit_single_residue_loops_;

};

} // namespace loops
} // namespace protocols



#endif //INCLUDED_protocols_loops_LoopsFileOptions_hh
