// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/Constraint.hh
/// @brief Base class definition for Constraint class hierarchy.

#ifndef INCLUDED_core_scoring_constraints_Constraint_hh
#define INCLUDED_core_scoring_constraints_Constraint_hh

// Unit header
#include <core/scoring/constraints/Constraint.fwd.hh>

// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// AUTO-REMOVED #include <core/id/AtomID.hh>
#include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>

//Utility Headers
#include <numeric/xyzVector.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ Headers
// AUTO-REMOVED #include <ostream>

#include <core/id/AtomID.fwd.hh>
// #include <core/id/SequenceMapping.fwd.hh> Commented by Brian Weitzner to fix compilation
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <utility/vector1.hh>
#include <sstream>

#ifdef WIN32
	#include <core/id/SequenceMapping.hh>
#endif


namespace core {
namespace scoring {
namespace constraints {

/// @brief Actually a *restraint*, like a virtual rubber band between a pair of atoms.
/// @details All Constraints are expected to be immutable once created,
/// meaning their internal data (state) should not change over their lifetime.
/// This allows Constraints to be shared between copies of Poses (e.g. in Monte Carlo),
/// and is important for both speed (with thousands of contraints) and correctness.
///
/// To "change" a constraint, remove the old one and add a new and different one.
/// The steal() methods have been removed because it is
/// incompatible with the idea of immutable constraints.

class Constraint : public utility::pointer::ReferenceCount {

public:
	typedef id::AtomID AtomID;

public:
	/// @brief Constructor for Constraint class.
	Constraint( ScoreType const & t ): score_type_(t) {}

	/// @brief Virtual destructor.
	virtual
	~Constraint();

	/// @brief Copies the data from this Constraint into a new object and returns
	/// an OP to the new object. Intended to be implemented by derived classes and
	/// used by pose.add_constraint
	virtual ConstraintOP clone() const = 0;

	virtual ConstraintOP clone( FuncOP ) const {
		unimplemented_method_error( std::string("clone" ) );
		return NULL;
	}


	/// @brief Copies the data from this Constraint into a new object and returns
	/// an OP atoms are mapped to atoms with the same name in dest pose ( e.g.
	/// for switch from centroid to fullatom ) if a sequence_mapping is present
	/// it is used to map residue numbers .. NULL = identity mapping to the new
	/// object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const& /*src*/,
		pose::Pose const& /*dest*/,
		id::SequenceMappingCOP map=NULL ) const {
		unimplemented_method_error( std::string("remapped_clone" ) );
		if ( !map ) return NULL; // to make compile happy
		return NULL;
	}

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
	// DON'T USE THIS ONE.. Most Constraint classes have not overloaded this one, but read_def ! OL
	virtual
	void read_constraint( std::istream & /*in*/, core::pose::Pose const & /*pose*/) {
		unimplemented_method_error( std::string( "read_constraint" ) );
	}

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
	virtual void read_data( std::istream & ) {}

	/// @brief apply a resid remapping to this constraint, returns the remapped
	/// constraint Does this return an owning pointer to this constraint or a
	/// copy? Documentation would be nice.
	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &/*seqmap*/ ) const
	{
		unimplemented_method_error( std::string( "remap_resid" ) );
		return NULL;
	}


	/// @brief return the "raw" distance before handed to the FUNC object
	virtual
	core::Real
	dist( core::pose::Pose const & /*pose*/ ) const {
		unimplemented_method_error( std::string( "dist" ) );
		return -1.0;
	}

	virtual
	core::Real
	dist( XYZ_Func const & /*xyz*/ ) const {
		unimplemented_method_error( std::string( "dist" ) );
		return -1.0;
	};

	/// @brief Calculates a score for this constraint using XYZ_Func, and puts
	/// the UNWEIGHTED score into emap. Although the current set of weights
	/// currently is provided, Constraint objects should put unweighted scores
	/// into emap because the ScoreFunction will do the weighting itself.
	virtual
	void
	score( XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const = 0;

	/// @brief Returns a unique string identified for this constraint. Used in several
	/// places, including the ConstraintIO class.
	virtual
	std::string type() const {
		return "UNKNOWN_TYPE";
	}

	// do some pre-scoring calculations -- does nothing by default
 	virtual void setup_for_scoring( XYZ_Func const &, ScoreFunction const & ) const {}

	// call the setup_for_derivatives for each constraint -- does nothing by default
 	virtual void setup_for_derivatives( XYZ_Func const &, ScoreFunction const & ) const {}

	/// @brief Returns the score of this constraint computed over the given conformation.
	/// Not necessarily implemented in all derived classes, as it's redundant with
	/// the score( XYZ_Func, EnergyMap, EnergyMap )  method defined above. Returns
	/// 0.0 if not implemented.
	virtual
	Real
	score( conformation::Conformation const &	) const { return 0.0; }

	/// @brief Fill the f1 and f2 vectors, necessary for considering the
	/// derivative this constraint during minimization. (someone please reference
	/// Bill Wedermeyer's paper here, as I'm in an airport and can't fill it in
	/// myself!)
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		XYZ_Func const & xyz_func,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const = 0;

	/// @brief This method is intended to show the value of the Constraint function
	/// evaluated over some reasonable range of values. For example, a constraint
	/// between pairs of atoms might show the values of the Constraint function
	/// between 4 and 12 angstroms.
	virtual void show( std::ostream & /*out*/ ) const {
		unimplemented_method_error( std::string( "show" ) );
	}

	/// @brief Prints the definition of a Constraint to the given std::ostream,
	/// using the given Pose, and the given FuncFactory. This method is intended
	/// to be overridden by derived classes if they'd like to use the
	/// ConstraintIO machinery. It's also not clear why this method takes a Pose,
	/// other than to be symmetric with read_def.
	virtual void show_def( std::ostream & /*out*/, pose::Pose const & ) const {
		unimplemented_method_error( std::string( "show_def" ) );
	}

	// @brief Reads the definition of a Constraint from the given std::istream,
	// using the given Pose, and the given FuncFactory. This method is intended
	// to be overridden by derived classes if they'd like to use the
	// ConstraintIO machinery.
	virtual void read_def( std::istream &, pose::Pose const &, FuncFactory const & ) {
		unimplemented_method_error( std::string( "read_def" ) );
	}

	// @brief take coordinates, distances, angles, etc from given pose
	///
	virtual void steal_def( pose::Pose const& ) {
		unimplemented_method_error( std::string( "steal_def" ) );
	}

	/// @brief Convenience function, returns the results of show() as a string.
	/// Not to be overriden by derived classes.
	std::string to_string() const {
		std::ostringstream out;
		show(out);
		return out.str();
	}

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

	/// @brief Returns the Func object associated with this Constraint object.
	virtual
	Func const & get_func() const {
		unimplemented_method_error( std::string( "get_func" ) );
		static HarmonicFunc dummy_func( 0.0, 0.0);
		return dummy_func; // satisfy compiler
	}


	/// @brief possibility to do object comparison instead
	/// of pointer comparison
	virtual
	bool operator == ( Constraint const & /*other*/ ) const{
		unimplemented_method_error( std::string( "== operator" ) );
		return false;
	}

	/// @brief possibility to do object comparison instead
	/// of pointer comparison
	bool operator != ( Constraint const & other ) const{
		return !(*this == other);
	}

	virtual
	core::Size choose_effective_sequence_separation(
			core::kinematics::ShortestPathInFoldTree const& sp,
			numeric::random::RandomGenerator&
	) {
		return effective_sequence_separation( sp );
	}


	virtual
	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& ) const {
		return 0;
	}

private:
	ScoreType const score_type_;

	/// @brief Utility method for producing useful error messages and exiting
	/// from program. Declared const which is funny, because exiting the program
	/// certainly changes the state of this object! This might be replaced with
	/// exception handling if we ever start using those inside of mini.
	void unimplemented_method_error( std::string const & method_name ) const {
		utility_exit_with_message(
			"Called Constraint::" + method_name + " method from derived class " +
			type() + "," + "ended up in Constraint::" + method_name + "\n"
		);
	}
}; // class Constraint


inline std::ostream& operator<< ( std::ostream & out, Constraint & cst ) {
	cst.show( out );
	return out;
}


} // constraints
} // scoring
} // core

#endif
