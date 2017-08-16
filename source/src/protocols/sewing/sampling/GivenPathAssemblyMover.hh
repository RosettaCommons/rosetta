// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file GivenPathAssemblyMover.hh
///
/// @brief This Mover generates assemblies by taking the best edge possible for each addition. This process is repeated
/// cycles times and the best overall Assembly that was found is returned
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_sampling_GivenPathAssemblyMover_HH
#define INCLUDED_protocols_sewing_sampling_GivenPathAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/GivenPathAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/AssemblyMover.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace sewing  {

class GivenPathAssemblyMover : public AssemblyMover {

private:

	typedef AssemblyMover parent;

public:

	GivenPathAssemblyMover();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	AssemblyOP
	generate_assembly() override;

	// void
	// follow_random_edge_from_node(
	//  AssemblyOP assembly,
	//  ModelNode const * node
	// );

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:

	utility::vector1<int> model_ids_;
	utility::vector1<ModelNode const *> nodes_;

};

} //sewing
} //protocols

#endif
