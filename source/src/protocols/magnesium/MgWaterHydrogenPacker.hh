// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgWaterHydrogenPacker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgWaterHydrogenPacker_HH
#define INCLUDED_protocols_magnesium_MgWaterHydrogenPacker_HH

#include <protocols/moves/Mover.hh>
#include <protocols/magnesium/MgWaterHydrogenPacker.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/UniformRotationSampler.fwd.hh>

namespace protocols {
namespace magnesium {

class MgWaterHydrogenPacker: public moves::Mover {

public:

	//constructor
	MgWaterHydrogenPacker();

	//constructor
	MgWaterHydrogenPacker( utility::vector1< Size > const & mg_res_list );

	//destructor
	~MgWaterHydrogenPacker();

public:

	using Mover::apply;

	virtual void apply( core::pose::Pose & pose );

	void
	apply( core::pose::Pose & pose,
		std::pair< core::Size, core::Size > const & mg_water );


	virtual std::string get_name() const{ return "MgWaterHydrogenPacker"; }

	void remove_waters_except_mg_bound( core::pose::Pose & pose ) const;

	void set_excise_mini_pose( bool const & setting ){ excise_mini_pose_ = setting; }
	bool excise_mini_pose() const { return excise_mini_pose_; }

	void set_use_fast_heuristic( bool const & setting ){ use_fast_heuristic_ = setting; }
	bool use_fast_heuristic() const { return use_fast_heuristic_; }

private:

	void
	pack_mg_water_hydrogens_in_pose( core::pose::Pose & pose,
		std::pair< core::Size, core::Size > const & mg_water_pair );

	core::Real
	get_heuristic_water_hydrogen_score( utility::vector1< core::Vector > const & acc_vecs,
		utility::vector1< core::Vector > const & rep_vecs,
		core::Vector const & mg_vec,
		core::Vector const & OH1c,
		core::Vector const & OH2c,
		numeric::xyzMatrix< core::Real > const & R ) const;

	void
	find_water_neighbor_vecs( core::pose::Pose const & pose,
		core::Size const water_res,
		core::Size const mg_res,
		utility::vector1< core::Vector > & acc_vecs,
		utility::vector1< core::Vector > & rep_vecs,
		core::Vector & mg_vec ) const;

	bool
	rotate_water_away_from_magnesium( core::pose::Pose & pose,
		core::Size const seqpos,
		core::Vector const & O,
		core::Vector const & OH1c,
		core::Vector const & OH2c,
		core::Vector const & MG,
		numeric::xyzMatrix< core::Real > const & R ) const;
private:

	utility::vector1< core::Size > mg_res_list_;
	utility::vector1< std::pair< core::Size, core::Size > > mg_water_pairs_;
	bool excise_mini_pose_;
	bool use_fast_heuristic_;
	numeric::UniformRotationSamplerCOP urs_;

};

} //magnesium
} //protocols

#endif
