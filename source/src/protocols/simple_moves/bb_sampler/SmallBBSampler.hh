// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/bb_sampler/SmallBBSampler.hh
/// @brief A bb sampler that samples within a range of a starting angle.  Similar to small mover.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_hh
#define INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_hh

#include <protocols/simple_moves/bb_sampler/SmallBBSampler.fwd.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>

namespace protocols {
namespace simple_moves {
namespace bb_sampler {

///@brief A bb sampler that samples within a range of a starting angle.  Similar to small mover.
/// Samples an angle +/- the current angle.  If angle max is not set, will sample any angle.
///
/// new_angle = old_angle +/- ( angle_max/2 )
///
class SmallBBSampler : public bb_sampler::BBDihedralSampler {

public:

	SmallBBSampler();

	SmallBBSampler( core::id::MainchainTorsionType torsion_type);
	SmallBBSampler( core::id::MainchainTorsionType torsion_type, core::Real all_ss_angle_max);

public:

	core::Real
	get_torsion( core::pose::Pose const & pose, core::Size resnum ) const;

	void
	set_torsion_to_pose( core::pose::Pose & pose, core::Size resnum ) const;


public:

	/// @brief Sets the maximum angle of perturbation, independent of
	/// secondary structure. new_angle = old_angle +/- ( angle_max/2 )
	///
	/// Example:
	///     bbmover.angle_max(25)
	/// See also:
	///     ShearMover
	///     SmallMover
	void
	set_angle_max( core::Angle const angle );

	/// @brief Sets the max angle of perturbation for residues with <type>
	/// secondary structure.  (<type> must be 'H', 'E', or 'L'.) new_angle = old_angle +/- ( angle_max/2 )
	///
	/// Example:
	///     bbmover.angle_max('H', 25)
	///
	/// See also:
	///     ShearMover
	///     SmallMover
	void
	set_angle_max( char const type, core::Angle const angle );

	// Note: Pass in by value for one-direction assignment.
	/// @brief Sets the max angle of perturbation, for secondary structures
	/// 'H', 'E', and 'L'.
	/// new_angle = old_angle +/- ( angle_max/2 )
	///
	void
	set_angle_max( std::map< char, core::Angle > angle_max_in );

	/// @brief Gets the max angle of perturbation for residues with <type>
	/// secondary structure.  (<type> must be 'H', 'E', or 'L'.)
	///
	/// Example:
	///     bbmover.angle_max('H')
	///
	/// See also:
	///     ShearMover
	///     SmallMover
	core::Real
	get_angle_max(char const type) const ;


	SmallBBSampler(SmallBBSampler const & src);

	virtual ~SmallBBSampler();

	SmallBBSamplerOP
	clone() const;

private:

	/// max allowed angle-change as a function of ss type
	std::map< char, core::Angle > angle_max_;

};


} //protocols
} //simple_moves
} //bb_sampler



#endif //INCLUDED_protocols_simple_moves_bb_sampler_SmallBBSampler_hh





