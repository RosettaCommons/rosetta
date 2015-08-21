// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgHydrater.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgHydrater_HH
#define INCLUDED_protocols_magnesium_MgHydrater_HH

#include <protocols/moves/Mover.hh>
#include <protocols/magnesium/MgHydrater.fwd.hh>
#include <protocols/magnesium/MgWaterHydrogenPacker.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/UniformRotationSampler.fwd.hh>

namespace protocols {
namespace magnesium {

class MgHydrater: public moves::Mover {

public:

	//constructor
	MgHydrater();

	//constructor
	MgHydrater( utility::vector1< Size > const & mg_res_list );

	//destructor
	~MgHydrater();

	void set_use_fast_frame_heuristic( bool const & setting ){ use_fast_frame_heuristic_ = setting; }
	bool use_fast_frame_heuristic() const { return use_fast_frame_heuristic_; }

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "MgHydrater"; }

	void set_excise_mini_pose( bool const & setting ){ excise_mini_pose_ = setting; }
	bool excise_mini_pose() const { return excise_mini_pose_; }

	void set_verbose( bool const & setting ){ verbose_ = setting; }
	bool verbose() const { return verbose_; }

private:

	void
	hydrate_magnesium( core::pose::Pose & pose, core::Size const i );

	void
	hydrate_magnesium_in_pose( core::pose::Pose & pose, core::Size const i,
		bool force_full_shell = true );

	numeric::xyzMatrix< core::Real >
	set_frame( core::Vector const & orig, core::Vector const & xyz1, core::Vector const & xyz2 ) const;

	bool
	hydrate_magnesium_with_orbital_frame( core::pose::Pose & pose,
		core::Size const i,
		utility::vector1< core::id::AtomID > const & nbr_atom_ids,
		numeric::xyzMatrix< core::Real > const & R,
		bool force_full_shell = true ) const;
private:

	utility::vector1< core::Size > mg_res_list_;
	bool excise_mini_pose_;
	bool use_fast_frame_heuristic_;
	numeric::UniformRotationSamplerCOP urs_;
	MgWaterHydrogenPackerOP mg_water_hydrogen_packer_;
	bool verbose_;
};

} //magnesium
} //protocols

#endif
