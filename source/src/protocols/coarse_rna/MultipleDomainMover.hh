// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_MultipleDomainMover_HH
#define INCLUDED_protocols_rna_MultipleDomainMover_HH

#include <protocols/coarse_rna/MultipleDomainMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <protocols/coarse_rna/CoarseRNA_LoopCloser.fwd.hh>
#include <protocols/coarse_rna/CoarseRNA_LoopCloser.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.fwd.hh>

#include <core/types.hh>

//// C++ headers
#include <string>


using core::Size;
using core::Real;

namespace protocols {
namespace farna {

/// @brief The RNA de novo structure modeling protocol
class MultipleDomainMover: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object
	MultipleDomainMover( core::pose::Pose const & pose, protocols::coarse_rna::CoarseRNA_LoopCloserOP rna_loop_closer  );

	/// @brief Clone this object
	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new MultipleDomainMover(*this) );
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	Size
	apply_and_return_jump( core::pose::Pose & pose );

	Size
	apply_at_domain( core::pose::Pose & pose, Size const & n );

	void
	randomize_pose_rigid_bodies( core::pose::Pose & pose );

	void
	slide_back_to_origin( core::pose::Pose & pose );

	Size
	num_domains(){ return num_domains_;}

	void
	update_rot_trans_mag( Real const & rot_mag, Real const & trans_mag );


private:
	void
	initialize( core::pose::Pose const & pose, protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map );

	void
	setup_jump_numbers_and_partner( core::pose::Pose const & pose );

	void
	setup_ok_for_centroid_calculation( protocols::toolbox::AtomLevelDomainMapOP & atom_level_domain_map );

	void
	randomize_orientations( core::pose::Pose & pose );

	void
	try_to_slide_into_contact( core::pose::Pose & pose );

	void
	close_all_loops( core::pose::Pose & pose );

	void
	initialize_rigid_body_movers();

	numeric::xyzVector< core::Real >
	get_centroid( core::pose::Pose const & pose );

private:

	bool verbose_;
	Real rot_mag_;
	Real trans_mag_;
	Size num_domains_;
	protocols::coarse_rna::CoarseRNA_LoopCloserOP rna_loop_closer_; //Later can make this a "general" loop closer.
	utility::vector1< int > jump_numbers_;
	utility::vector1< protocols::rigid::Partner > partner_;
	utility::vector1< bool > ok_for_centroid_calculation_;
	utility::vector1< protocols::rigid::RigidBodyPerturbMoverOP > rb_movers_;

}; // class MultipleDomainMover


} //farna
} //protocols

#endif
