// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ExhaustiveAssemblyMover.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_devel_sewing_sampling_ExhaustiveAssemblyMover_HH
#define INCLUDED_devel_sewing_sampling_ExhaustiveAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/ExhaustiveAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/AssemblyMover.hh>

#include <protocols/sewing/conformation/Assembly.hh> // for accessing segments_ and all_segments_


//Protocol headers
#include <core/pose/Pose.hh>

#include <protocols/sewing/conformation/Model.hh>

namespace protocols {
namespace sewing  {

class ExhaustiveAssemblyMover : public AssemblyMover {

private:

	typedef AssemblyMover parent;

public:

	ExhaustiveAssemblyMover();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	virtual
	AssemblyOP
	generate_assembly();

	void
	parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	);

private:
	core::Size bogus_var_for_constructor_;

	std::map< int, Model > models;

};

} // sewing
} // protocols

#endif
