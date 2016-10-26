// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/Constraint.hh
/// @brief Base class definition for Constraint class hierarchy.

#ifndef INCLUDED_core_scoring_constraints_Constraint_hh
#define INCLUDED_core_scoring_constraints_Constraint_hh

// Unit header
#include <core/scoring/constraints/Constraint.fwd.hh>

// package headers
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>

// project headers
#include <core/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// #include <core/id/SequenceMapping.hh> Removing the SequenceMapping-as-default-NULL parameter to remove this #include

// utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// numeric headers
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>

// C++ Headers
#include <sstream>

#ifdef WIN32
#include <core/id/SequenceMapping.hh>
#endif

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief Actually a *restraint*, like a virtual rubber band between a pair of atoms.
///
/// @details All Constraints are expected to be immutable once created,
/// meaning their internal data (state) should not change over their lifetime.
/// This allows Constraints to be shared between copies of Poses (e.g. in Monte Carlo),
/// and is important for both speed (with thousands of constraints) and correctness.
///
/// To "change" a constraint, remove the old one and add a new and different one.
/// The steal() methods have been removed because it is
/// incompatible with the idea of immutable constraints.

class Constraint : public utility::pointer::ReferenceCount {

public:
	typedef id::AtomID AtomID;

public:

	/// @brief Constructor for Constraint class.
	Constraint( ScoreType const & t );

	/// @brief Virtual destructor.
	virtual ~Constraint();

	/// @brief Copies the data from this %Constraint into a new object and returns
	/// an OP to the new object. Intended to be implemented by derived classes and
	/// used by pose.add_constraint.  This function must return a *deep copy* of
	/// itself -- meaning that if this %Constraint holds pointers to other %Constraints
	/// that it must invoke clone on those %Constraints as well.  If the %Constraint
	/// holds a FuncOP, then the Func should also be cloned.
	virtual ConstraintOP clone() const = 0;

	/// @brief Clone the constraint, but where a new Func object is to be used instead.
	virtual ConstraintOP clone( func::FuncOP ) const;

	/// @brief Copies the data from this Constraint into a new object and returns
	/// an OP atoms are mapped to atoms with the same name in dest pose ( e.g.
	/// for switch from centroid to fullatom ) if a sequence_mapping is present
	/// it is used to map residue numbers .. NULL = identity mapping to the new
	/// object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const& /*src*/,
		pose::Pose const& /*dest*/,
		id::SequenceMappingCOP map ) const;

	/// @brief Returns the number of atoms involved in defining this constraint.
	/// If the constraint doesn't depend on particular atoms (e.g. a residue type constraint)
	/// this function can return zero
	/// @details Note that this function isn't actually used by the constraint scoring machenery.
	/// If you're calling it on a generic Constraint (as opposed to specifically on a derived class)
	/// you're probably doing something wrong.
	virtual
	Size
	natoms() const = 0;

	/// @brief Returns the AtomID referred to by index.
	/// @details Note that this function isn't actually used by the constraint scoring machenery.
	/// If you're calling it on a generic Constraint (as opposed to specifically on a derived class)
	/// you're probably doing something wrong.
	virtual
	AtomID const &
	atom( Size const index ) const = 0;

	/// @brief Returns the pose numbers of the residues involved in this constraint, in no particular order.
	/// @details Used in determining one-body/two-body/multi-body status.
	/// For historical reasons, the default uses a simple protocol based on natoms()/atom() -
	/// feel free to reimplement more efficiently.
	virtual
	utility::vector1< core::Size >
	residues() const;

	/// @brief This method is totally redundant with read_def  YAY
	/// DON'T USE THIS ONE.. Most Constraint classes have not overloaded this one, but read_def ! OL
	virtual
	void read_constraint( std::istream & /*in*/, core::pose::Pose const & /*pose*/ );

	/// @brief Returns the ScoreType that this Constraint object will use.
	ScoreType const &
	score_type() const
	{
		return score_type_;
	}

	/// @brief initialize this Constraint from the given std::istream. It's amazing
	/// that there are three functions for doing this inside of Constraint.hh.
	/// SO WHAT IS THIS SUPPOSED TO DO ? not overloaded by e.g.,  AtomPairConstraint or CoordinateConstraint,
	// -- use read_def() if in doubt.
	virtual void read_data( std::istream & );

	/// @brief apply a resid remapping to this constraint, returns the remapped
	/// constraint Does this return an owning pointer to this constraint or a
	/// copy? Documentation would be nice.
	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const & /*seqmap*/ ) const;

	/// @brief Calculates a score for this constraint using XYZ_Func, and puts
	/// the UNWEIGHTED score into emap. Although the current set of weights
	/// currently is provided, Constraint objects should put unweighted scores
	/// into emap because the ScoreFunction will do the weighting itself.
	virtual
	void
	score( core::scoring::func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const = 0;

	/// @brief Returns the unweighted score of this constraint computed over the given pose.
	virtual
	Real
	score( pose::Pose const& pose ) const;

	/// @brief Returns the weighted score of this constraint computed over the given pose.
	virtual
	Real
	score( pose::Pose const& pose,  EnergyMap const & weights ) const;

	/// @brief return the raw "distance" before that distance is handed to the FUNC object
	/// @details - If such a distance doesn't make sense for this constraint, just return 0
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & /*xyz*/ ) const = 0;

	/// @brief return the raw "distance" before that distance is handed to the FUNC object
	/// @details - If such a distance doesn't make sense for this constraint, just return 0
	virtual
	Real
	dist( core::pose::Pose const & /*pose*/ ) const;

	/// @brief Returns a unique string identified for this constraint. Used in several
	/// places, including the ConstraintIO class.
	virtual
	std::string type() const;

	// do some pre-scoring calculations -- does nothing by default
	virtual void setup_for_scoring( core::scoring::func::XYZ_Func const &, ScoreFunction const & ) const;

	// call the setup_for_derivatives for each constraint -- does nothing by default
	virtual void setup_for_derivatives( core::scoring::func::XYZ_Func const &, ScoreFunction const & ) const;

	/// @brief Fill the f1 and f2 vectors, necessary for considering the
	/// derivative this constraint during minimization. (someone please reference
	/// Bill Wedermeyer's paper here, as I'm in an airport and can't fill it in
	/// myself!)
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz_func,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const = 0;

	/// @brief This method is intended to show the value of the Constraint function
	/// evaluated over some reasonable range of values. For example, a constraint
	/// between pairs of atoms might show the values of the Constraint function
	/// between 4 and 12 angstroms.
	virtual void show( std::ostream & /*out*/ ) const;

	/// @brief Prints the definition of a Constraint to the given std::ostream,
	/// using the given Pose, and the given func::FuncFactory. This method is intended
	/// to be overridden by derived classes if they'd like to use the
	/// ConstraintIO machinery. It's also not clear why this method takes a Pose,
	/// other than to be symmetric with read_def.
	virtual void show_def( std::ostream & /*out*/, pose::Pose const & ) const;

	/// @brief Reads the definition of a Constraint from the given std::istream,
	/// using the given Pose, and the given func::FuncFactory. This method is intended
	/// to be overridden by derived classes if they'd like to use the
	/// ConstraintIO machinery.
	virtual void read_def( std::istream &, pose::Pose const &, core::scoring::func::FuncFactory const & );

	/// @brief take coordinates, distances, angles, etc from given pose
	virtual void steal_def( pose::Pose const& );

	/// @brief Convenience function, returns the results of show() as a string.
	/// Not to be overriden by derived classes.
	std::string to_string() const;

	/// @brief Prints the violations of this constraint to the given
	/// std::ostream. What are violations? It's not defined, and it depends on
	/// the constraint and the function!  also - wtf is threshold? it was defined
	/// as a Size in CoordinateConstraint, I don't know which definition is the
	/// right one. Documentation would be nice ...
	virtual Size show_violations(
		std::ostream & out,
		pose::Pose const &,
		Size,
		Real threshold = 1
	) const;

	/// @brief Returns the func::Func object associated with this Constraint object.
	virtual
	core::scoring::func::Func const & get_func() const;


	/// @brief Equality operator to test whether two constraints are functionally
	/// identical.
	///
	/// @details This operator should use floating point comparison and should not
	/// decide that two floats are identical if they are within some epsilon > 0.
	/// This method allows developes to remove specific constraints from Poses, even
	/// if the constraints have been cloned.  Remapped constraints should not be
	/// considered identical -- i.e., if cst1 is between residues i and j and
	/// cst2 is between residues i+1 and j+1.  All subclasses of Constraint must
	/// implement this method.
	virtual
	bool operator == ( Constraint const & /*other*/ ) const = 0;

	/// @brief Inequality operator to test whether two constraints are not functionally
	/// identical.
	bool operator != ( Constraint const & other ) const;

	/// @brief Determine if the calling class has the same type as the
	/// input class, using dynamic casts.  This is important in ensuring
	/// that two constraints are equal: both this and other must check
	/// that the other is the same as it.  This is not an optional method
	/// and every class should implement it, regaurdless of whether a
	/// parent class implements it.
	virtual
	bool
	same_type_as_me( Constraint const & other ) const = 0;

	virtual
	core::Size choose_effective_sequence_separation(
		core::kinematics::ShortestPathInFoldTree const& sp,
		numeric::random::RandomGenerator&
	);


	virtual
	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& ) const;


private:
	ScoreType const score_type_;

	/// @brief Utility method for producing useful error messages and exiting
	/// from program. Declared const which is funny, because exiting the program
	/// certainly changes the state of this object! This might be replaced with
	/// exception handling if we ever start using those inside of mini.
	void unimplemented_method_error( std::string const & method_name ) const;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Constraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class Constraint


std::ostream& operator<< ( std::ostream & out, Constraint const & cst );


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_Constraint )
#endif // SERIALIZATION


#endif
