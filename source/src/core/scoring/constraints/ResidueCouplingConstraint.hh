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
/// @brief implements favor_native_residue constraints, assigning an energy bonus
/// to a residue if it matches the identity within the constraint object
/// @author Moritz Ertelt


#ifndef INCLUDED_core_scoring_constraints_ResidueCouplingConstraint_hh
#define INCLUDED_core_scoring_constraints_ResidueCouplingConstraint_hh

#include <core/scoring/constraints/ResidueCouplingConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/id/AtomID.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <numeric/MathNTensor.hh>
#include <numeric/MathNTensor_io.hh>

namespace core {
namespace scoring {
namespace constraints {

typedef numeric::MathNTensor< double, 4 > CouplingTensor;

typedef utility::pointer::shared_ptr < CouplingTensor > CouplingTensorOP;
typedef utility::pointer::shared_ptr < const CouplingTensor > CouplingTensorCOP;

typedef utility::pointer::shared_ptr < std::vector<int> > indexListOP;
typedef utility::pointer::shared_ptr < const std::vector<int> > indexListCOP;

/// @brief This class favors a particular residue identity at two particular positions by adding the scaled tensor value to the energy map.
class ResidueCouplingConstraint : public Constraint
{
public:
	typedef core::Real Real;
	ResidueCouplingConstraint();

	ResidueCouplingConstraint(
		pose::Pose const & pose,
		Size const& seqpos1,
		Size const& seqpos2,
		Size const& tensor_index1,
		Size const& tensor_index2,
		CouplingTensorCOP const& tensor_const,
		Real const& strength,
		std::string const& alphabet
	);

	~ResidueCouplingConstraint() override;

	Size
	natoms() const override;

	AtomID const &
	atom( Size const ) const override;

	utility::vector1< core::Size >
	residues() const override;

	void
	set_alphabet(std::string const& alphabet);

	std::string
	get_alphabet() const;

	void
	show( std::ostream & out ) const override;

	/// and not just pointers
	bool operator == ( Constraint const & other ) const override;
	bool same_type_as_me( Constraint const & ) const override;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const override;

	/// @details Return 1.0 if constraint will get a bonus/penalty, 0.0 if not
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const override;

	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const override;

	ConstraintOP
	clone() const override;

private:
	Size seqpos1_;
	Size seqpos2_;
	Size tensor_index1_;
	Size tensor_index2_;
	CouplingTensorCOP tensor_const_;
	Real strength_;
	std::string alphabet_;
	std::map< char, long unsigned int > amino_acids_;
}; // RotamerConstraint


} // namespace constraints
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_constraints_ResidueCouplingConstraint_HH
