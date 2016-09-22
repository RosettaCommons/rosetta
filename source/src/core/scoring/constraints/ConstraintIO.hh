// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IO-functionality for Constraints
/// @brief
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_core_scoring_constraints_ConstraintIO_hh
#define INCLUDED_core_scoring_constraints_ConstraintIO_hh

// Unit headers
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

// Package headers
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/FuncFactory.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

class ConstraintIO : public utility::SingletonBase< ConstraintIO >
{
public:
	friend class utility::SingletonBase< ConstraintIO >;

public:
	static ConstraintSetOP read_constraints(
		std::string const & filename,
		ConstraintSetOP cst_set,
		pose::Pose const& pose );
	static ConstraintSetOP read_constraints(
		std::istream & data,
		ConstraintSetOP cst_set,
		pose::Pose const& pose );

	static void write_constraints( std::ostream&, ConstraintSet const& cst_set, core::pose::Pose const& );
	static void write_constraints( std::string const& filename, ConstraintSet const& cst_set, core::pose::Pose const& );

	static func::FuncFactory & get_func_factory();
	static ConstraintFactory & get_cst_factory();

	static ConstraintOP parse_atom_pair_constraint(
		std::istream & data,
		core::pose::Pose pose
	);

	static ConstraintSetOP read_constraints_new(
		std::string const & fname,
		ConstraintSetOP cset,
		pose::Pose const & pose
	);

	static ConstraintSetOP read_constraints_new(
		std::istream & data,
		ConstraintSetOP cset,
		pose::Pose const & pose
	);

	/// @brief read one individual constraint defined.
	static ConstraintOP read_individual_constraint_new(
		std::istream & data,
		core::pose::Pose const& pose,
		func::FuncFactory const & func_factory
	);

	/// @brief read one individual constraint defined.
	static ConstraintOP read_individual_constraint_new(
		std::istream & data,
		core::pose::Pose const& pose,
		func::FuncFactory const & func_factory,
		std::string type /*cst -type*/
	);

	static ConstraintOP parse_coordinate_constraint(
		std::istream & data,
		core::pose::Pose pose
	);

	static void parse_residue(
		pose::Pose const& pose,
		std::string const & residue_string,
		Size & residue_num
	);

	static Size
	parse_residue( pose::Pose const& pose, int const resnum, char const chain = 'A' );

	// gkt - tmp hack for BoundFunc, should be private
protected:
	static void read_cst_bindingsites( std::istream &data, std::string& next_section,  ConstraintSet&, pose::Pose const&  );
	static void read_cst_atom_pairs( std::istream &data, std::string& next_section,  ConstraintSet&, pose::Pose const&  );
	static void read_cst_coordinates( std::istream &data, std::string& next_section,  ConstraintSet&, pose::Pose const&  );
	static void read_cst_angles( std::istream &data, std::string& next_section,  ConstraintSet&, pose::Pose const&  );

private:
	ConstraintIO () {};

	// should the func_factory_ just be a singleton instead of being held inside this class?
	func::FuncFactory func_factory_;
	//static ConstraintFactory cst_factory_;

};

} //constraints
} //scoring
} //core

#endif
