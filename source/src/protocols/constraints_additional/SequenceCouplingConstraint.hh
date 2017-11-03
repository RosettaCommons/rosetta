// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraints_additional/SequenceCouplingConstraint.hh
/// @brief  This is a constraint that refers to a core::sequence::SequenceCoupling? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design).
/// @author HetuKamisetty

#ifndef INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_hh
#define INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_hh

#include <protocols/constraints_additional/SequenceCouplingConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/sequence/SequenceCoupling.fwd.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

/// @brief
class SequenceCouplingConstraint : public core::scoring::constraints::Constraint {
public:
	typedef core::sequence::SequenceCoupling SequenceCoupling;
	typedef core::sequence::SequenceCouplingOP SequenceCouplingOP;
	typedef core::sequence::SequenceCouplingCOP SequenceCouplingCOP;
	typedef core::id::SequenceMapping SequenceMapping;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Conformation Conformation;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::func::FuncFactory FuncFactory;
	typedef core::scoring::func::XYZ_Func XYZ_Func;
	typedef core::scoring::constraints::ConstraintOP ConstraintOP;

public:

	SequenceCouplingConstraint();

	SequenceCouplingConstraint(
		Pose const &,
		core::Size,
		core::Size,
		SequenceCouplingOP profile = nullptr
	);

	SequenceCouplingConstraint(
		core::Size,
		core::Size,
		SequenceCouplingOP profile = nullptr
	);

	~SequenceCouplingConstraint() override;

	ConstraintOP clone() const override;

	bool operator == ( core::scoring::constraints::Constraint const & other ) const override;

	bool same_type_as_me( core::scoring::constraints::Constraint const & other ) const override;

	std::string type() const override { return "SequenceCoupling"; }

	/// @brief used by ConstraintIO and ConstraintFactory to construct this constraint from a input file stream (constraint file)
	void
	read_def(
		std::istream &,
		Pose const &,
		FuncFactory const &
	) override;

	void show_def( std::ostream &, Pose const & ) const override;

	void show( std::ostream & out ) const override;

	core::Size seqpos1() const { return seqpos1_; }
	core::Size seqpos2() const { return seqpos2_; }

	void set_sequence_coupling( SequenceCouplingOP );
	SequenceCouplingOP sequence_coupling();
	SequenceCouplingCOP sequence_coupling() const;

	core::Size natoms() const override { return 0; };
	AtomID const & atom( Size const ) const override {
		utility_exit_with_message("SequenceCouplingConstraint is not atom-based!.");
		return core::id::GLOBAL_BOGUS_ATOM_ID;  // required for compilation on Windows
	}

	utility::vector1< core::Size > residues() const override;

	//virtual ConstraintOP remap_resid( SequenceMapping const & ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual void
	score(
		XYZ_Func const &,
		EnergyMap const &,
		EnergyMap &
	) const override;

	/// @details "Distance" for SequenceCouplingConstraint isn't all that simple
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & ) const override { return 0; }

	virtual void
	fill_f1_f2(
		AtomID const &,
		XYZ_Func const &,
		core::Vector &,
		core::Vector &,
		EnergyMap const &
	) const override;

private:
	core::Size seqpos1_;
	core::Size seqpos2_;
	SequenceCouplingOP sequence_coupling_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace constraints_additional
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraints_additional_SequenceCouplingConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_HH
