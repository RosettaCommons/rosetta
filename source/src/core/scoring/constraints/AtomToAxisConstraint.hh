// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   AtomToAxisConstraint.fwd.hh
/// @author Jack Maguire

#ifndef INCLUDED_core_scoring_constraints_AtomToAxisConstraint_hh
#define INCLUDED_core_scoring_constraints_AtomToAxisConstraint_hh

#include <core/scoring/constraints/AtomToAxisConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
//#include <cstdlib>
//#include <iosfwd>
//#include <map>
//#include <utility>

#ifdef    SERIALIZATION
// Cereal headers
#include <utility/vector1.srlz.hh>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

using XYZ = numeric::xyzVector< core::Real >;

struct XYZs {
	XYZ atom;
	XYZ axis1;
	XYZ axis2;
};


class AtomToAxisConstraint : public Constraint {
public:

	// default c-tor
	AtomToAxisConstraint();

	///c-tor
	AtomToAxisConstraint(
		AtomID const & a,
		utility::vector1< AtomID > const & axis1,
		utility::vector1< AtomID > const & axis2,
		core::scoring::func::FuncOP const & func,
		ScoreType scoretype = atom_pair_constraint
	);

	AtomToAxisConstraint( AtomToAxisConstraint const & ) = default;

	/// @brief Create a deep copy of this %AtomToAxisConstraint, cloning its Func
	ConstraintOP clone() const override;

	/// @brief Create a deep copy of this %AtomToAxisConstraint except that the copy should
	/// use the input Func instead of its existing one.
	ConstraintOP clone( func::FuncOP func ) const override;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. nullptr = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = nullptr
	) const override;

	/// @brief Compare a) the class types (w/ same_type_as_me), b) the atoms being constrained,
	/// c) the score_type being used, and d) the Func objects (the FuncOPs may point at different
	/// objects, but as long as those objects are equal, that counts)
	/// @details Beware of noise in our floating point comparison checks
	bool operator == ( Constraint const & other ) const override;
	bool same_type_as_me( Constraint const & other ) const override;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	Real
	score(
		Vector const & xyz1,
		Vector const & xyz2
	) const;

	void
	score( core::scoring::func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const override;

	Real score( pose::Pose const & pose ) const override;

	XYZs
	calc_xyzs( pose::Pose const & pose ) const;

	Real dist( pose::Pose const & pose ) const override;

	XYZs
	calc_xyzs( core::scoring::func::XYZ_Func const & xyz ) const;

	Real dist( core::scoring::func::XYZ_Func const & xyz ) const override;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const override;

	std::string type() const override {
		return "AtomPair";
	}


	Size
	natoms() const override
	{
		return 1 + axis1_.size() + axis2_.size();
	}

	AtomID const &
	atom( core::Size const ii ) const override {
		if ( ii == 1 ) return atom_;

		core::Size const jj = ii - 1;
		if ( jj <= axis1_.size() ) {
			return axis1_[ jj ];
		} else {
			return axis2_[ jj - axis1_.size() ];
		}
	}

	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const override;


	AtomID const & get_atom() const { return atom_; }

	void show( std::ostream& out ) const override;
	void show_def( std::ostream& out, pose::Pose const& pose ) const override;

	void read_def( std::istream& in, pose::Pose const& pose, func::FuncFactory const & func_factory ) override;
	// //@brief set constraint such that the pose doesn't violate it.
	// virtual void steal( pose::Pose& );

	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const override;

	func::Func const & get_func() const override;

	void setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const override;

private:
	// functions
	Real
	func( Real const theta ) const;

	// deriv
	Real
	dfunc( Real const theta ) const;

protected:
	/// @brief Setter for the derived class
	void set_func( func::FuncOP setting );

	void atom( AtomID const & newid ) {
		atom_ = newid;
	}

	void add_atom_to_axis1( AtomID const & newid ) {
		axis1_.push_back( newid );
	}

	void add_atom_to_axis2( AtomID const & newid ) {
		axis2_.push_back( newid );
	}

private:
	AtomID atom_;
	utility::vector1< AtomID > axis1_;
	utility::vector1< AtomID > axis2_;

	core::scoring::func::FuncOP func_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class AtomToAxisConstraint

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_AtomToAxisConstraint )
#endif // SERIALIZATION


#endif
