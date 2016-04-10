// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RehashAssemblyMover.hh
///
/// @brief A Mover that uses more sewing to close gaps from discontinous assemblies
/// @author Tim Jacobs

#ifdef NOT_IN_SCONS_DEPRECATED

#ifndef INCLUDED_protocols_sewing_sampling_RehashAssemblyMover_HH
#define INCLUDED_protocols_sewing_sampling_RehashAssemblyMover_HH

//Unit headers
#include <protocols/sewing/sampling/RehashAssemblyMover.fwd.hh>

//Package headers
#include <protocols/sewing/sampling/AssemblyMover.hh>
#include <protocols/sewing/conformation/Assembly.hh>

//Protocol headers
#include <core/pose/Pose.hh>

#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>

#include <protocols/loops/Loops.hh>


namespace protocols {
namespace sewing  {
	
class RehashAssemblyMover : public AssemblyMover {
	
public:

	RehashAssemblyMover();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;
	
	std::string
	get_name() const;

	virtual
	bool
	generate_assembly(
		AssemblyOP & assembly
	);

	//bool
	//rearrange_assembly(
	//	AssemblyOP & assembly
	//) const;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

private:

	core::Real max_loop_distance_;
	std::map< int, Model > bridge_models_;

};
	
} //sewing
} //protocols

#endif

#endif
