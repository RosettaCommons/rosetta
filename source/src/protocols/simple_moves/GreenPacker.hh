// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/GreenPacker.hh
/// @brief  packing mover that makes extensive reuse of rotamer pair energies class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_GreenPacker_hh
#define INCLUDED_protocols_simple_moves_GreenPacker_hh

/// Unit headers
#include <protocols/simple_moves/GreenPacker.fwd.hh>

/// Package headers
#include <protocols/moves/Mover.hh>

/// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

/// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace simple_moves {

class MinimalRotamer : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MinimalRotamer();
	typedef core::chemical::AA          AA;
	typedef core::chemical::ResidueType ResidueType;
	typedef core::conformation::Residue Residue;
	typedef core::Size            Size;
	typedef core::Vector          Vector;
	typedef core::Real				Real;

public:
	MinimalRotamer( Residue const & );

	bool
	same( MinimalRotamer const & ) const;

	bool
	same_residue_type( MinimalRotamer const & ) const;

	AA
	aa() const;

private:

	bool
	has_ideal_geometry(
		core::conformation::Residue const & res
	) const;

	bool
	atom_is_ideal(
		core::conformation::Residue const & res,
		Size const atom_id
	) const;

	void
	record_chi( core::conformation::Residue const & res );

	void
	record_internal_geometry( core::conformation::Residue const & res );

	void
	record_internal_geometry(
		core::conformation::Residue const & res,
		Size const atom_id
	);

	bool
	chi_matches_coords( core::conformation::Residue const & res, Size chi_index ) const;

	bool
	same_chi( MinimalRotamer const & other ) const;

	bool
	same_nonideal_geometry( MinimalRotamer const & other ) const;


private:
	ResidueType const & residue_type_;
	bool ideal_geometry_;
	utility::vector1< Real > chi_;

	// Record of internal geometry using the tree-geometry defined in the
	// ResidueType.
	// steal xyzVector to record internal geometry:
	// x = 0 = d, y = 1 = theta, z = 2 = phi
	utility::vector1< Vector > internal_geometry_;

	static Size const d = 0;
	static Size const theta = 1;
	static Size const phi = 2;

	/// No default constructor or assignment operator
	/// copy c-tor is fine
	MinimalRotamer();
	MinimalRotamer const & operator = ( MinimalRotamer const & );
};

/// @brief Interface class used to break a pose down into a set of component "groups"
/// where intra-group rotamer-pair energies are preserved between calls to the
/// GreenPacker.  E.g. in rigid-body docking between two proteins, chains 1 and 2 define
/// groups 1 and 2.  In rigid-body docking between two domains of the same chain, those
/// residues upstream of jump 1 define group 1, and those downstream of jump 1 define group 2.
/// In loop modelling, the static background is group 1, and the loop itself is group 0, since
/// loop residues will have their bb dofs change regularly between repackings.
class GroupDiscriminator : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Size       Size;

public:
	virtual ~GroupDiscriminator();

	virtual
	protocols::simple_moves::GroupDiscriminatorOP clone() const = 0;

	virtual
	Size
	group_id( Pose const & pose, Size seqpos ) const = 0;

};

class UserDefinedGroupDiscriminator : public GroupDiscriminator
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Size       Size;

public:
	virtual ~UserDefinedGroupDiscriminator();

	virtual
	protocols::simple_moves::GroupDiscriminatorOP clone() const;

	virtual
	Size
	group_id( Pose const & pose, Size seqpos ) const;

	void
	set_group_ids( utility::vector1<Size > const & group_ids_input );

private:

	utility::vector1< Size > group_ids_;

};

class ChainGroupDiscriminator : public GroupDiscriminator
{

public:
	virtual ~ChainGroupDiscriminator();

	virtual
	protocols::simple_moves::GroupDiscriminatorOP clone() const;

	virtual
	Size
	group_id( Pose const & pose, Size seqpos ) const;
};

class GreenPacker : public protocols::moves::Mover
{
public:
	/// Types
	typedef core::graph::Graph                             Graph;
	typedef core::graph::GraphOP                           GraphOP;
	typedef core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph PrecomputedPairEnergiesInteractionGraph;
	typedef core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP PrecomputedPairEnergiesInteractionGraphOP;
	typedef core::pack::rotamer_set::RotamerSets           RotamerSets;
	typedef core::pack::rotamer_set::RotamerSetsOP         RotamerSetsOP;
	typedef core::pack::task::TaskFactory                  TaskFactory;
	typedef core::pack::task::TaskFactoryOP                TaskFactoryOP;
	typedef core::pack::task::PackerTaskOP                 PackerTaskOP;
	typedef core::pose::Pose                               Pose;
	typedef core::scoring::EnergyMap                       EnergyMap;
	typedef core::scoring::ScoreFunction                   ScoreFunction;
	typedef core::scoring::ScoreFunctionOP                 ScoreFunctionOP;
	typedef core::scoring::ScoreTypes                      ScoreTypes;
	typedef core::scoring::methods::LongRangeTwoBodyEnergy LongRangeTwoBodyEnergy;
	typedef core::Size                                     Size;
	typedef core::Real                                     Real;
	typedef core::Vector                                   Vector;

public:
	GreenPacker();
	virtual ~GreenPacker();

	virtual
	void
	apply( Pose & );
	virtual std::string get_name() const;

	// Undefined, commentin out to make PyRosetta compile
	// void reset();

	void
	set_scorefunction( ScoreFunction const & );

	void
	set_group_discriminator( protocols::simple_moves::GroupDiscriminatorOP );

	void
	set_task_factory( TaskFactoryOP );

	void
	set_reference_round_task_factory( TaskFactoryOP );

private:
	/// Private methods
	void
	set_weights_for_sfxn(
		ScoreFunction & sfxn,
		ScoreTypes const & scoretypes,
		EnergyMap const & weights
	) const;


	void setup_reference_data( Pose & );
	void repack( Pose & pose );

	void split_pose_into_groups( Pose & pose );
	void create_reference_packer_task( Pose & pose );
	void create_reference_packer_neighbor_graph( Pose & pose );
	void create_reference_rotamers( Pose & pose );
	void compute_reference_intragroup_rpes( Pose & pose );

	void create_fresh_task( Pose & pose );
	void create_fresh_packer_neighbor_graph( Pose & pose );
	void create_fresh_rotamers( Pose & pose );
	void find_reference_and_current_rotamer_correspondence( Pose & pose );
	void initialize_internal_correspondence_data( Pose & pose );
	void compute_energies( Pose & pose );
	void run_sa( Pose & pose );
	void cleanup();

	void add_precomputed_energies( Pose & pose, PrecomputedPairEnergiesInteractionGraphOP pig );
	void compute_absent_energies( Pose & pose, PrecomputedPairEnergiesInteractionGraphOP pig );

	void
	compute_absent_srci_energies_for_residue_pair(
		Pose & pose,
		PrecomputedPairEnergiesInteractionGraphOP pig,
		Size lower_res,
		Size upper_res
	);

	void
	compute_absent_lrci_energies_for_residue_pair(
		Pose & pose,
		LongRangeTwoBodyEnergy const & lre,
		PrecomputedPairEnergiesInteractionGraphOP pig,
		Size lower_res,
		Size upper_res
	);

	void drop_inter_group_edges( Pose & pose, GraphOP packer_neighbor_graph ) const;
	void drop_intra_group_edges( Pose & pose, GraphOP packer_neighbor_graph ) const;

	void store_reference_pose_geometry( Pose & pose );
	void compare_input_pose_geometry_to_reference( Pose & pose );


private:
	ScoreFunctionOP full_sfxn_;
	ScoreFunctionOP ci_sfxn_;
	ScoreFunctionOP cd_sfxn_;

	//@brief true at construction, or after a call to reset, before a call to apply.
	bool create_reference_data_;

	protocols::simple_moves::GroupDiscriminatorOP group_discriminator_;
	utility::vector1< Size > group_ids_;
	std::vector< utility::vector1< Size > > group_members_; // index by zero to represent group 0

	TaskFactoryOP reference_task_factory_;
	TaskFactoryOP task_factory_;

	////////////////////////////////////////////
	/// Data captured from the reference apply()

	/// @brief task used in construction of the reference data
	PackerTaskOP reference_task_;
	GraphOP reference_packer_neighbor_graph_;

	/// @brief rotamer sets created by the reference task
	RotamerSetsOP reference_rotamer_sets_;
	utility::vector1< Size > reference_resid_2_moltenres_;
	utility::vector1< Size > reference_moltenres_2_resid_;

	/// @brief the stored intra-group RPEs from the context independent components of the
	/// score function
	PrecomputedPairEnergiesInteractionGraphOP ci_rpes_;

	/// @brief the internal geometry of the rotamers created in the first round
	utility::vector1< utility::vector1< protocols::simple_moves::MinimalRotamerOP > > original_rotamers_;

	/// @brief the internal geometry of the backbone used in the first round
	utility::vector1< utility::vector1< Real > > orig_bb_tors_; // compare input tors against

	/// @brief for debugging; a set of coordinates of each residue for backbone atom #1 (N in proteins, Phos for NA's)
	/// so that if the backbone DOFs do change, and the structure of the input pose cannot be superimposed
	/// back onto the original coordinates, then a debug runtime_assert statementmet will get caught.
	utility::vector1< Vector > original_bb_rep_coords_;

	//////////////////////////////////////////////////////////////////////////////////////////////
	/// Data for the "current" packing -- valid between the beginning and end of a call to "apply"

	PackerTaskOP current_task_;
	RotamerSetsOP current_rotamer_sets_;
	PrecomputedPairEnergiesInteractionGraphOP current_ig_;

	/// @brief the internal geometry of the rotamers created in the first round
	utility::vector1< utility::vector1< protocols::simple_moves::MinimalRotamerOP > > current_rotamers_;

	GraphOP current_packer_neighbor_graph_;
	GraphOP current_inter_group_packer_neighbor_graph_;
	GraphOP current_intra_group_packer_neighbor_graph_;

	/// correspondence between current rotamer set and original rotamer set.
	utility::vector1< utility::vector1< Size > > orig_rot_2_curr_rot_;
	utility::vector1< utility::vector1< Size > > curr_rot_2_orig_rot_;
	utility::vector1< utility::vector1< Size > > curr_rotamers_with_correspondence_;
	utility::vector1< utility::vector1< Size > > curr_rotamers_without_correspondence_;


};

} // namespace moves
} // namespace protocols

#endif // INCLUDED_protocols_simple_moves_GreenPacker_HH
