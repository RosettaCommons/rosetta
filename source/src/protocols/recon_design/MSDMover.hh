// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/MSDMover.hh
/// @brief Multistate design mover used for RECON multistate design
/// Takes in multiple poses, applies residue linking constraints based on
/// sequence of all input poses and runs a design submover that has been specified in the tag
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_recon_design_MSDMover_hh
#define INCLUDED_protocols_recon_design_MSDMover_hh

// Unit headers
#include <protocols/recon_design/MSDMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>

#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace recon_design {

/// @brief Multistate design mover used for RECON multistate design
/// Takes in multiple poses, applies residue linking constraints based on
/// sequence of all input poses and runs a design submover that has been specified in the tag
/// Only accessible through recon application.
class MSDMover : public moves::VectorPoseMover {

public:
	/// @brief empty constructor fills values with the values
	/// read in from the commandline
	MSDMover();

	/// @brief Constructor takes in the design mover and a list of
	/// resfiles that specifies the residues to constrain to one another
	MSDMover( protocols::moves::MoverOP mover,
		utility::vector1< std::string > const & resfiles );

	~MSDMover() override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	/// @brief Read options from RosettaScripts
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &
	) override;

	/// @brief Specify XML schema
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Movers derived from VectorPoseMover must define two apply methods:
	/// apply and apply_mpi, depending on whether the protocol is run in MPI or not
	/// Running in MPI will distribute each pose to a separate node, so the mover will
	/// operate on one pose
	void apply( Pose & pose ) override;
	void apply_mpi( Pose & pose ) override;

	/// @brief Based on the sequence of the other poses, apply a residue type constraint to
	/// encourage poses to adopt same sequence. Returns a COP of the constraints that were
	/// applied so they can be removed later
	utility::vector1< core::scoring::constraints::ConstraintCOP >
	apply_linked_constraints( Pose & pose, utility::vector1< utility::vector1< std::string > > other_pose_sequences,
		utility::vector1< core::Size > my_designable_residues );

	/// @brief Initialize mover by checking that input poses were passed correctly,
	/// a design mover was specified, and finding the pose given to apply() in the
	/// poses_ vector
	void setup_mover ( Pose & pose );

	/// Getters and setters
	moves::MoverOP design_mover();
	void design_mover( moves::MoverOP design_mover );

	utility::vector1< std::string > resfiles();
	void resfiles ( utility::vector1< std::string > resfiles );

	void weight( core::Real weight );
	core::Real weight();

	utility::vector1< utility::vector1< core::Size > > designable_residues();

private:
	/// @brief Populates designable_residues_ with a list of designable residues corresponding
	/// to the poses in poses_, corresponding element-wise (designables_residues[0] matches
	/// to poses_[0], etc)
	void parse_resfiles();

	/// @brief The design mover can't be given its task operations on a global
	/// level in the RosettaScript, bc each of the input poses can have
	/// different designable residues. This function assigns the design mover
	/// its tasks based on the designable residues of the current pose
	void update_packer_task();

	/// @brief Assign the index of the pose that will be operated on by
	/// the current mover
	void set_current_pose ( core::Size current_pose );

	/// @brief Get the resfile corresponding to the pose at index. If
	/// only one resfile is present then it will be returned regardless
	/// of the value of index.
	std::string resfile_at ( core::Size index );

	/// @brief Runs the post design mover, if present
	void run_post_design_mover ( core::pose::Pose & pose );

	/// Private class members

	/// Mover used for design. Typically a PackRotamersMover but not necessar
	protocols::moves::MoverOP design_mover_;

	/// Mover applied after design
	protocols::moves::MoverOP post_mover_;

	/// List of resfile file names, must be either 1 resfile, or same length as poses_
	utility::vector1< std::string > resfiles_;

	/// Weight to give the residue linking constraints.
	core::Real weight_;

	/// Index of the pose to be designed
	core::Size current_pose_;

	/// Vector of designable residues for each pose in poses_.
	/// Must be same length as poses_.
	/// Also each pose must have same number of designable residues
	/// (i.e. designable_residues[1].size() == designable_residues[2].size())
	utility::vector1< utility::vector1< core::Size > > designable_residues_;

	/// Output extra debug messages
	bool debug_;
};


} // recon_design
} // protocols

#endif //INCLUDED_protocols_recon_design_MSDMover_HH
