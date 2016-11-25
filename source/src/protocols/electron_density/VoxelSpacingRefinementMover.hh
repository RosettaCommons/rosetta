// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file

#ifndef INCLUDED_protocols_electron_density_VoxelSpacingRefinementMover_hh
#define INCLUDED_protocols_electron_density_VoxelSpacingRefinementMover_hh

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace electron_density {


/// Bfactor multifunc
class VoxelSpacingMultifunc : public core::optimization::Multifunc {
public:
	VoxelSpacingMultifunc( core::pose::Pose const &pose );

	~VoxelSpacingMultifunc() override = default;

	core::Real
	operator ()( core::optimization::Multivec const & vars ) const override;

	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const override;

	void
	dump( core::optimization::Multivec const & x1, core::optimization::Multivec const & x2 ) const override;

	void
	getMapSpacingAndOrigin( core::optimization::Multivec & x1, bool aniso );

	void
	foldInChanges( core::pose::Pose &pose, core::optimization::Multivec & y );

private:
	ObjexxFCL::FArray3D< double > rhoC_; // f calc
};


/// mover to fit B factors
class VoxelSpacingRefinementMover : public moves::Mover {
public:
	VoxelSpacingRefinementMover();
	~VoxelSpacingRefinementMover() override = default;

	void init();

	void apply( core::pose::Pose & ) override;

	// XRW TEMP  std::string get_name() const override { return "VoxelSpacingRefinement"; }

	moves::MoverOP clone() const override { return moves::MoverOP( new VoxelSpacingRefinementMover( *this ) ); }
	moves::MoverOP fresh_instance() const override { return moves::MoverOP( new VoxelSpacingRefinementMover ); }

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string minimizer_, mapout_;
	bool aniso_;
	core::Size max_iter_;
};

} // moves
} // protocols

#endif
