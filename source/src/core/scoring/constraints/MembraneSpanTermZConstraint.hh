// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueConstraint.hh
///
/// @brief constraint the Z axis distance between TM span edges
/// @details in hetero fold and dock (w/w- design) rosetta commonly shifts the chains along the Z axis
/// if you prefer them not shofting, this constraint penalizes the pose by the distance on the Z axis between the span edges
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)


#ifndef INCLUDED_core_scoring_constraints_MembraneSpanTermZConstraint_hh
#define INCLUDED_core_scoring_constraints_MembraneSpanTermZConstraint_hh

#include <core/scoring/constraints/MembraneSpanTermZConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.hh>


#include <utility/vector1.hh>

// Membrane headers
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief This class favors a particular residue identity at a particular position by reducing its res_type energy.
class MembraneSpanTermZConstraint : public Constraint
{
public:
	typedef core::Real Real;
public:
	MembraneSpanTermZConstraint();

	MembraneSpanTermZConstraint(
		pose::Pose const & pose
	);

	virtual ~MembraneSpanTermZConstraint();

	virtual
	Size
	natoms() const { return 0; }

	virtual
	AtomID const &
	atom( Size const ) const {
		utility_exit_with_message("MembraneSpanTermZConstraint is not atom-based!.");
		return core::id::GLOBAL_BOGUS_ATOM_ID;  // required for compilation on Windows
	}

	virtual
	utility::vector1< core::Size >
	residues() const;

	void
	show( std::ostream & out ) const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const & ) const;

	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual
	void
	score( func::XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const;

	core::Real
	calc_score() const;

	core::Real
	score_delta( core::Real a, core::Real b ) const;

	/// @details Return 1.0 if constraint will get a bonus, 0.0 if not
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	virtual ConstraintOP
	clone() const;
private:
	utility::vector1< core::conformation::ResidueCOP > in_res_;
	utility::vector1< core::conformation::ResidueCOP > out_res_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // RotamerConstraint

typedef utility::pointer::shared_ptr< MembraneSpanTermZConstraint > MembraneSpanTermZConstraintOP;
typedef utility::pointer::shared_ptr< MembraneSpanTermZConstraint const > MembraneSpanTermZConstraintCOP;

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_MembraneSpanTermZConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_constraints_MembraneSpanTermZConstraint_HH
