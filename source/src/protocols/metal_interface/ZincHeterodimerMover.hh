// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/metal_interface/ZincHeterodimerMover.hh
/// @brief ZincHeterodimerMover protocol main mover
/// @author Steven Lewis

#ifndef INCLUDED_protocols_metal_interface_ZincHeterodimerMover_hh
#define INCLUDED_protocols_metal_interface_ZincHeterodimerMover_hh

#include <protocols/metal_interface/ZincHeterodimerMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/Edge.hh> //composition = .hh
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh> //needed?
#include <protocols/moves/Mover.hh> //inheritance = .hh
#include <utility/vector1.hh>


namespace protocols {
namespace metal_interface {

/// @details ZincHeterodimerMover is a mover for the ZincHeterodimerDesign protocol
class ZincHeterodimerMover : public protocols::moves::Mover {

public:

	/// @brief constructor
	/// @param[in] metal_site this is the residues directly involved in metal binding (plus the metal itself) - these should not be repacked.
	/// @param[in] fixed_to_metal this is a core::kinematic::Edge which represents a bond from the metal center to a liganding residue on the "fixed" partner (the one without rigid-body freedom)
	/// @param[in] metal_to_mobile this is a core::kinematic::Edge which represents a bond from the metal center to a liganding residue on the mobile partner with rigid body freedom.  The bond will be rotated by RotateJumpAxisMover
	ZincHeterodimerMover(
		utility::vector1< core::Size > const & metal_site,
		core::kinematics::Edge const & fixed_to_metal,
		core::kinematics::Edge const & metal_to_mobile
	);

	virtual ~ZincHeterodimerMover();

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:
	/// @brief constructor helper function - generates scorefunctions
	void generate_scorefunctions();
	/// @brief constructor helper function - generates TaskFactory
	void generate_factory();
	/// @brief apply helper function (ugly complex interface, I know)
	void copy_to_centroid(core::pose::Pose const & pose,
												core::pose::Pose & centroid,
												core::kinematics::FoldTree const & centroid_tree,
												core::Size const metal_res);

	//generated data
	core::scoring::ScoreFunctionOP centroid_scorefunction_;
	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryCOP factory_;

	/// @brief fixed_to_metal this is a core::kinematic::Edge which represents a bond from the metal center to a liganding residue on the "fixed" partner (the one without rigid-body freedom)
	core::kinematics::Edge const fixed_to_metal_;
	/// @brief metal_to_mobile this is a core::kinematic::Edge which represents a bond from the metal center to a liganding residue on the mobile partner with rigid body freedom.  The bond will be rotated by RotateJumpAxisMover
	core::kinematics::Edge const metal_to_mobile_;
	/// @brief the residues involved in the metal site (generally 4 residues plus the metal for 5 total)
	utility::vector1< core::Size > const metal_site_;

};//end ZincHeterodimerMover

}//metal_interface
}//protocols

#endif // INCLUDED_protocols_metal_interface_ZincHeterodimerMover_HH
