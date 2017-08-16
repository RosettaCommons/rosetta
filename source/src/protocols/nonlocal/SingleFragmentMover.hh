// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/SingleFragmentMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_SINGLEFRAGMENTMOVER_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_SINGLEFRAGMENTMOVER_HH

// Unit header
#include <protocols/nonlocal/SingleFragmentMover.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/nonlocal/Chunk.hh>
#include <protocols/nonlocal/Policy.hh>

namespace protocols {
namespace nonlocal {

class SingleFragmentMover : public protocols::moves::Mover {
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::fragment::FragSetOP FragSetOP;
	typedef core::fragment::Frame Frame;
	typedef core::kinematics::FoldTree FoldTree;
	typedef core::kinematics::FoldTreeOP FoldTreeOP;
	typedef core::kinematics::MoveMapOP MoveMapOP;
	typedef core::pose::Pose Pose;

	typedef boost::unordered_map<Size, Frame> FrameMap;
	typedef utility::vector1<Chunk> Chunks;

public:
	typedef protocols::moves::Mover Mover;
	typedef protocols::moves::MoverOP MoverOP;
	/// @brief No-argument constructor required by RosettaScripts. The caller is
	/// responsible for initializing the instance
	SingleFragmentMover();

	/// @brief Creates a new SingleFragmentMover that selects uniformly from the
	/// available fragments at the selected insertion position.
	SingleFragmentMover(const FragSetOP& fragments,
		const MoveMapOP& movable);

	/// @brief Creates a new SingleFragmentMover that selects fragments at the
	/// selected insertion position using the given policy.
	SingleFragmentMover(const FragSetOP& fragments,
		const MoveMapOP& movable,
		const PolicyOP& policy);

	~SingleFragmentMover() override = default;

	/// @brief Performs a single fragment insertion on <pose>, drawn from the set
	/// of fragments specified in the constructor. Respects the underlying
	/// kinematics of the system, as determined by the Pose's FoldTree and the
	/// user-specified MoveMap. Fragment insertions will only occur in allowable
	/// regions of the pose. No moves will span jumps in the FoldTree.
	///
	/// Assumptions:
	///   - <pose> has been instantiated (i.e. constructed from sequence or an
	///     alternate source) elsewhere
	///   - <pose>'s FoldTree is valid
	///   - <pose> is centroid-level (warning if full-atom)
	///   - The combination of FoldTree and MoveMap provide an unambiguous
	///     definition of what's movable and what's not
	void apply(Pose& pose) override;

	/// @brief Returns the name of this mover

	/// @brief Creates a new instance using the copy constructor
	MoverOP clone() const override;

	/// @brief Creates a new instance using the default constructor
	MoverOP fresh_instance() const override;

	/// @brief Mover-specific parsing required by RosettaScripts
	void parse_my_tag(utility::tag::TagCOP tag,
		basic::datacache::DataMap& data,
		const protocols::filters::Filters_map& filters,
		const protocols::moves::Movers_map& movers,
		const Pose& pose) override;

	/// @brief Returns true if this instance is in a usable state, false otherwise
	bool valid() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	void initialize(const FragSetOP& fragments,
		const MoveMapOP& movable,
		const PolicyOP& policy);

	/// @brief Creates a position-indexable list of Frames
	void initialize_library();

	/// @brief Creates a set of chunks by examining the stored MoveMap and FragSet
	/// in conjunction with <tree>
	void initialize_chunks(const FoldTree& tree);

	/// @brief Returns a randomly chosen chunk with uniform probability
	const Chunk* random_chunk() const;

	/// @brief If <pose> is fullatom, converts it to centroid and returns true.
	/// Otherwise, takes no action and returns false.
	bool to_centroid(Pose* pose) const;


	/// @brief The set of fragments to apply to the pose
	FragSetOP fragments_;

	/// @brief Defines restrictions on which degrees of freedom in the system can
	/// be modified.
	MoveMapOP movable_;

	/// @brief Selects the fragment to be inserted at <insertion_pos> given
	/// knowledge of the fragment library and the current state of the pose.
	PolicyOP policy_;

	/// @brief FoldTree used to initialize <chunks_> in a previous call to apply()
	FoldTreeOP previous_tree_;

	/// @brief Provides index-based access to the data contained in the FragSet
	FrameMap library_;

	/// @brief Regions of sequence on which to perform fragment insertion
	Chunks chunks_;

	/// @brief Probability of selecting chunk c_i. Proportional to chunk length.
	utility::vector1<Real> probs_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_SINGLEFRAGMENTMOVER_HH_
