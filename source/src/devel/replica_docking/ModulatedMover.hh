// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_ModulatedMover_hh
#define INCLUDED_devel_replica_docking_ModulatedMover_hh

#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/HamiltonianExchange.hh>
#include <protocols/canonical_sampling/HamiltonianExchange.fwd.hh>

#include <devel/replica_docking/TempInterpolator.hh>
#include <devel/replica_docking/TempInterpolator.fwd.hh>

#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/ReferenceCount.hh>


namespace devel {
namespace replica_docking {

class ModulatedMover : public protocols::canonical_sampling::ThermodynamicMover {
	typedef ThermodynamicMover Parent;
	//  typedef std::map< std::string, devel::replica_docking::TempInterpolatorOP > Interpolators;
	typedef std::map< std::string, devel::replica_docking::TempInterpolatorBaseOP > Interpolators;
	typedef utility::vector1< protocols::canonical_sampling::ThermodynamicMoverOP > MoverOPs;
	typedef utility::vector1< core::Size > GridCoord;

public:
	ModulatedMover();
	ModulatedMover( ModulatedMover const & );

	~ModulatedMover() override;
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	bool preserve_detailed_balance() const override { return true; }

	void set_preserve_detailed_balance( bool ) override {};

	utility::vector1< core::id::TorsionID_Range > torsion_id_ranges( core::pose::Pose & ) override {
		return utility::vector1< core::id::TorsionID_Range>();
	}

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	utility::tag::TagCOP generate_mover_tag(
		core::Size temp_level,
		std::string const& prefix,
		std::map< std::string, std::string > const& common_options
	) const;

	void
	initialize_simulation(
		core::pose::Pose& pose,
		protocols::canonical_sampling::MetropolisHastingsMover const& mhm,
		core::Size cycle
	) override;

	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & mhm
	) override;

	//   virtual void
	//   observe_after_metropolis(
	//      protocols::canonical_sampling::MetropolisHastingsMover const & mhm
	//   );

private:

	//  protocols::canonical_sampling::ThermodynamicMoverOP mover_;
	//  protocols::canonical_sampling::TemperatureControllerOP tempering_; // inteprete tempering as TemperatureController does not work for n_temp_levels, current_temp and so on, so here I directly interprete it as HamiltonianExchange
	protocols::canonical_sampling::HamiltonianExchangeOP tempering_;
	//  core::Size n_temp_levels_ ;
	MoverOPs movers_;
	//  Interpolators interpolators_; // OP or like this?
	Interpolators interps_1_;
	Interpolators interps_2_; // for the 2nd dimension
	std::string mover_name_; //mover_creator, mover_type_name, mover_creator_key
};

} // namespace replica_docking
}

#endif
