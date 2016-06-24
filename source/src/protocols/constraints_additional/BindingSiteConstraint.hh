// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file source/src/protocols/constraints_additional/BindingSiteConstraint.cc
/// @author Frank DiMaio (?)
/// @brief This constrains some set of three or more atoms to maintain the same geometry as the starting pose. Constrainted binding sites are read in using the ConstraintIO class (maybe only with the defunct section-based constraints?)  Binding-site constraints use the same weight as atom-pair constraints, and (should) work with any protocol where atom-pair constraints are respected. Constraints on sidechain atoms are automatically converted to constraints on centroids when performing centroid-level manipulation.

///File format example:
/*
[ bindingsites ]
   CG 37  CD 38 CA 39 OG 40  CD1 41
   CA 69  CA 70 CA 71 CA 72  CA  73
[ atompairs ]
   CA  3 CA 49 HARMONIC 4.44 2.0
   CA  5 CA 71 HARMONIC 4.68 2.0
*/

#ifndef INCLUDED_protocols_constraints_additional_BindingSiteConstraint_hh
#define INCLUDED_protocols_constraints_additional_BindingSiteConstraint_hh

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

class BindingSiteConstraint : public core::scoring::constraints::Constraint {
public:

	/// null constructor
	BindingSiteConstraint( ) :
		core::scoring::constraints::Constraint( core::scoring::atom_pair_constraint )  /// ? TO DO -- give own scoretype
	{ }

	/// ctor from atom list + input pose
	BindingSiteConstraint(
		utility::vector1< AtomID > const & atms,
		core::pose::Pose const &start_pose,
		core::scoring::ScoreType scoretype = core::scoring::atom_pair_constraint   /// ? TO DO -- give own scoretype
	);

	/// ctor from a vector of atom positions (in lieu of a pose)
	BindingSiteConstraint(
		utility::vector1< AtomID > const & atms,
		ObjexxFCL::FArray2D< core::Real >  tgt_pos,
		ObjexxFCL::FArray2D< core::Real >  tgt_pos_centroid,
		core::scoring::ScoreType scoretype = core::scoring::atom_pair_constraint   /// ? TO DO -- give own scoretype
	);

	virtual core::scoring::constraints::ConstraintOP clone() const {
		return core::scoring::constraints::ConstraintOP( new BindingSiteConstraint( atms_, tgt_pos_, tgt_pos_centroid_ ) );
	}

	bool operator == ( core::scoring::constraints::Constraint const & other ) const;

	bool same_type_as_me( core::scoring::constraints::Constraint const & other ) const;

	void
	score( core::scoring::func::XYZ_Func const & xyz, core::scoring::EnergyMap const &, core::scoring::EnergyMap & emap ) const;

	// do some pre-scoring calculations
	void setup_for_scoring( core::scoring::func::XYZ_Func const & xyz, core::scoring::ScoreFunction const &scfxn ) const;

	// align the atoms
	//   ... placing a vector  -- from each atom to the the rotated >target< atoms -- in the database
	void pre_align( utility::vector1< numeric::xyzVector< core::Real > > const & templ_atms,
		utility::vector1< bool > const & ) const;

	// call the setup_for_derivatives for each constraint
	void setup_for_derivatives(  core::scoring::func::XYZ_Func const & xyz, core::scoring::ScoreFunction const &scfxn ) const;

	// atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const;

	std::string type() const;


	Size
	natoms() const;

	virtual
	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;


	AtomID const &
	atom( Size const n ) const;

	void show( std::ostream& out ) const;

	void show_def( std::ostream& out, core::pose::Pose const & pose ) const;
	void read_def( std::istream& in, core::pose::Pose const & pose, core::scoring::func::FuncFactory const & func_factory );

	Size show_violations( std::ostream & out, core::pose::Pose const & pose, core::Size verbose_level, core::Real threshold = 1.0 ) const;

protected:
	void init( core::pose::Pose const & start_pose );

private:
	// data
	utility::vector1< AtomID > atms_;
	ObjexxFCL::FArray2D< core::Real >  tgt_pos_;
	ObjexxFCL::FArray2D< core::Real >  tgt_pos_centroid_;

	// map of pos->tgt in rotated struct
	// This is an inappropriate use of global data
	static std::map< AtomID , numeric::xyzVector< core::Real > > rot_db;

	// database mapping constraints to RB transformations
	// static std::map< AtomID , ?? > transformDB;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraints_additional_BindingSiteConstraint )
#endif // SERIALIZATION


#endif
