// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh
/// @brief  TaskOperation class that finds an interface based on InterfaceVectorDefinition
/// and leaves it mobile in the PackerTask.  Serves mostly to wrap InterfaceVectorDefinition
/// into a TaskOperation. see src/core/select/util/interface_vector_calculate.hh
/// @author Ben Stranges (stranges@unc.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToInterfaceVectorOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToInterfaceVectorOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <utility/vector1.hh>


// Utility Headers
//#include <core/types.hh>

// C++ Headers
//#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictToInterfaceVectorOperation : public InterfaceTaskOperation
{
public:
	typedef InterfaceTaskOperation parent;

	//empty contstructor for parser, uses jumps
	RestrictToInterfaceVectorOperation();

	RestrictToInterfaceVectorOperation( core::Size const lower_chain_id, core::Size const upper_chain_id );

	//full constructor
	RestrictToInterfaceVectorOperation( core::Size const lower_chain, core::Size const upper_chain,
		core::Real CB_dist_cutoff,
		core::Real nearby_atom_cutoff,
		core::Real vector_angle_cutoff,
		core::Real vector_dist_cutoff);

	//basic jump constructor
	RestrictToInterfaceVectorOperation( utility::vector1_int const movable_jumps );

	//full constructor for jumps
	RestrictToInterfaceVectorOperation( utility::vector1_int const movable_jumps ,
		core::Real CB_dist_cutoff,
		core::Real nearby_atom_cutoff,
		core::Real vector_angle_cutoff,
		core::Real vector_dist_cutoff);


	// //if you want to use chain characters this is probably the best way, define the calculator separately
	// RestrictToInterfaceVectorOperation( std::string const & calculator );

	virtual ~RestrictToInterfaceVectorOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToInterfaceVector"; }

	//@brief, setters for the calculator.
	void upper_chain( core::Size upper_chain);
	void lower_chain( core::Size lower_chain);

	void upper_chain( utility::vector1<core::Size> upper_chain);
	void lower_chain( utility::vector1<core::Size> lower_chain);

	/// Commenting out to fix PyRosetta build  void jump_num( int jump_num);
	void CB_dist_cutoff( core::Real CB_dist_cutoff);
	void nearby_atom_cutoff(core::Real nearby_atom_cutoff);
	void vector_angle_cutoff(core::Real vector_angle_cutoff);
	void vector_dist_cutoff(core::Real vector_dist_cutoff);
	/// @brief parse_tag function for rosetta scripts
	void parse_tag( TagCOP tag , DataMap & );

	/*
	// Used to make the eventual inheritance frodm protocols::toolbox::task_operations::InterfaceTaskOperation easier.
	void
	setup_interface_chains_from_jumps( core::pose::Pose const & pose );
	*/

private:

	/// @brief private data used to pass to the definition function
	//chain ids of the interface lower=chain1 upper=chain2 for most purposes.
	utility::vector1<core::Size> lower_chains_;
	utility::vector1<core::Size> upper_chains_;
	bool jump_active_; //is the jump deffinition being used

	//cutoffs for various restrictions
	core::Real CB_dist_cutoff_; //distance for big CB cutoffs
	core::Real nearby_atom_cutoff_; // used for finding atoms that are close
	core::Real vector_angle_cutoff_; // used for cutoff for res1 CB to res2 CB angle cutoff
	core::Real vector_dist_cutoff_; // used for distance between CBs for vector
	//char upper_chain_char_, lower_chain_char;
	//int jump_vector_; //what jump is the interface across
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictToInterfaceVectorOperation_HH
