// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/constraints/InverseRotamersRCG.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2010


#ifndef INCLUDED_protocols_forge_constraints_InverseRotamersRCG_hh
#define INCLUDED_protocols_forge_constraints_InverseRotamersRCG_hh

//unit headers
#include <protocols/forge/constraints/InverseRotamersRCG.fwd.hh>

//package headers
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/build/Interval.fwd.hh>

//project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>

#include <list>

#ifdef PYROSETTA
#include <protocols/forge/build/Interval.hh>
#endif

namespace protocols {
namespace forge {
namespace constraints {


/// @brief a RemodelConstraintGenerator that creates AmbiguousConstraints for all positions
/// in a remodeled region towards a list of inverse rotamers. For every remodel position/inverse rotamer
/// pair, there will be one MultiConstraint consisting of three CoordinateConstraints. the three
/// coordinate constraints will be between:
/// 1) remodel res N  - invrot N coords
/// 2) remodel res Ca - invrot Ca coords
/// 3) remodel res Cb - invrot Cb coords
/// All of these MultiConstraints are combined to form one AmbiguousConstraint.
/// In effect, this RCG should bias the remodel trajectory such that
/// one remodel residue backbone overlays with one inverse rotamer backbone
class InverseRotamersRCG : public remodel::RemodelConstraintGenerator
{

public:
	// constructors and virtual functions
	InverseRotamersRCG();

	/// @brief copy constructor
	InverseRotamersRCG( InverseRotamersRCG const & rval );

	/// @brief value constructor which initializes and then calls init()
	InverseRotamersRCG(
		core::Size const lstart,
		core::Size const lstop,
		std::list< core::conformation::ResidueCOP > const & inverse_rotamers
	);

	~InverseRotamersRCG();

	void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	// XRW TEMP  virtual std::string
	// XRW TEMP  get_name() const;

	protocols::moves::MoverOP
	fresh_instance() const override;

	protocols::moves::MoverOP
	clone() const override;

	void
	generate_remodel_constraints( core::pose::Pose const & pose ) override;

public:
	// public member functions
	void
	set_constraint_func(
		core::scoring::func::FuncOP constraint_func );

	void
	clear_inverse_rotamers();

	using remodel::RemodelConstraintGenerator::init;

	/// @brief initialize with given lstart, lstop, and inverse rotamer residue list
	void init( core::Size const lstart,
		core::Size const lstop,
		std::list< core::conformation::ResidueCOP > const & inverse_rotamers );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1< forge::build::Interval > intervals_;

	std::list< core::conformation::ResidueCOP > inverse_rotamers_;

	core::scoring::func::FuncOP constraint_func_;
	core::Real func_sd_;
}; //class InverseRotamersRCG


} //namespace remodel
} //namespace forge
} //namespace protocols


#endif // INCLUDED_protocols_forge_remodel_InverseRotamersRCG_HH
