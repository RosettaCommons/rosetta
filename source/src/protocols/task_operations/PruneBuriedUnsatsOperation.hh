// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/PruneBuriedUnsatsOperation.hh
/// @brief  Eliminate rotamers that contain a sidechain buried unsat and make no other sidechain h-bonds
/// @detail This is intended to speed up packing by eliminating useless polar rotamers. This can
///         also help to reduce the number of buried unsats in designs because Rosetta can't pack them.
/// @author Longxing Cao -- original idea
/// @author Brian Coventry ( bcov@uw.edu ) -- code implementation


#ifndef INCLUDED_protocols_task_operations_PruneBuriedUnsatsOperation_hh
#define INCLUDED_protocols_task_operations_PruneBuriedUnsatsOperation_hh


// Unit Headers
#include <protocols/task_operations/PruneBuriedUnsatsOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace task_operations {

class PruneBuriedUnsats_RotamerSetsOperation : public core::pack::rotamer_set::RotamerSetsOperation
{

public:

	PruneBuriedUnsats_RotamerSetsOperation(
		bool allow_even_trades,
		core::Real atomic_depth_probe_radius,
		core::Real atomic_depth_resolution,
		core::Real atomic_depth_cutoff,
		core::Real minimum_hbond_energy,
		core::scoring::ScoreFunctionOP scorefxn_sc,
		core::scoring::ScoreFunctionOP scorefxn_bb
	);

	~PruneBuriedUnsats_RotamerSetsOperation() override;

	core::pack::rotamer_set::RotamerSetsOperationOP
	clone() const override;

	void
	alter_rotamer_sets(
		core::pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::pack::task::PackerTask const & ptask,
		utility::graph::GraphCOP,
		core::pack::rotamer_set::RotamerSets & rotamer_sets
	) override;


private: // data

	/// Allow residues that satisfy an unsat and create a new unsatisfiable one
	bool allow_even_trades_;

	core::Real atomic_depth_probe_radius_;
	core::Real atomic_depth_resolution_;
	core::Real atomic_depth_cutoff_;

	core::Real minimum_hbond_energy_;

	core::scoring::ScoreFunctionOP scorefxn_sc_;
	core::scoring::ScoreFunctionOP scorefxn_bb_;


};


//////////////////////////////////////////////////////////////////////////////////////////////
class PruneBuriedUnsatsOperation : public core::pack::task::operation::TaskOperation {
public:

	/// @brief default constructor
	PruneBuriedUnsatsOperation();

	/// @brief copy constructor
	PruneBuriedUnsatsOperation( PruneBuriedUnsatsOperation const & other );

	/// @brief destructor
	~PruneBuriedUnsatsOperation() override;

	PruneBuriedUnsatsOperation &
	operator=( PruneBuriedUnsatsOperation const & other );

	/// @brief make clone
	core::pack::task::operation::TaskOperationOP clone() const override;


	/// @brief apply
	void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const override;


	void parse_tag( TagCOP tag , DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "PruneBuriedUnsats"; }


public:

	bool allow_even_trades() const { return allow_even_trades_; }
	core::Real atomic_depth_probe_radius() const { return atomic_depth_probe_radius_; }
	core::Real atomic_depth_resolution() const { return atomic_depth_resolution_; }
	core::Real atomic_depth_cutoff() const { return atomic_depth_cutoff_; }
	core::Real minimum_hbond_energy() const { return minimum_hbond_energy_; }

	void allow_even_trades( bool allow_even_trades ) { allow_even_trades_ = allow_even_trades; }
	void atomic_depth_probe_radius( core::Real atomic_depth_probe_radius ) { atomic_depth_probe_radius_ = atomic_depth_probe_radius; }
	void atomic_depth_resolution( core::Real atomic_depth_resolution ) { atomic_depth_resolution_ = atomic_depth_resolution; }
	void atomic_depth_cutoff( core::Real atomic_depth_cutoff ) { atomic_depth_cutoff_ = atomic_depth_cutoff; }
	void minimum_hbond_energy( core::Real minimum_hbond_energy ) { minimum_hbond_energy_ = minimum_hbond_energy; }


private: // data

	/// Allow residues that satisfy an unsat and create a new unsatisfiable one
	bool allow_even_trades_;

	core::Real atomic_depth_probe_radius_;
	core::Real atomic_depth_resolution_;
	core::Real atomic_depth_cutoff_;

	core::Real minimum_hbond_energy_;

	core::scoring::ScoreFunctionOP scorefxn_sc_;
	core::scoring::ScoreFunctionOP scorefxn_bb_;

};


} // TaskOperations
} // protocols


#endif
