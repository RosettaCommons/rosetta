// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraints_additional/SequenceProfileConstraint.hh
/// @brief  This is a constraint that refers to a core::sequence::SequenceProfile? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design).
/// @author ashworth

#ifndef INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_hh
#define INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_hh

#include <protocols/constraints_additional/SequenceProfileConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

// AUTO-REMOVED #include <core/chemical/AA.hh>
#include <core/sequence/SequenceProfile.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

//Auto Headers
#include <utility/vector1_bool.hh>


namespace protocols {
namespace constraints_additional {

///@brief
class SequenceProfileConstraint : public core::scoring::constraints::Constraint {
public:
	typedef core::sequence::SequenceProfile SequenceProfile;
	typedef core::sequence::SequenceProfileOP SequenceProfileOP;
	typedef core::sequence::SequenceProfileCOP SequenceProfileCOP;
	typedef core::id::SequenceMapping SequenceMapping;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Conformation Conformation;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::constraints::FuncFactory FuncFactory;
	typedef core::scoring::constraints::XYZ_Func XYZ_Func;
	typedef core::scoring::constraints::ConstraintOP ConstraintOP;

public:

	SequenceProfileConstraint();

	SequenceProfileConstraint(
		Pose const &,
		core::Size,
		SequenceProfileOP profile = NULL
	);

	SequenceProfileConstraint(
		core::Size,
		utility::vector1< AtomID > const &,
		SequenceProfileOP profile = NULL
	);

	virtual ~SequenceProfileConstraint();

	virtual ConstraintOP clone() const;

	virtual std::string type() const { return "SequenceProfile"; }

	///@brief used by ConstraintIO and ConstraintFactory to construct this constraint from a input file stream (constraint file)
	virtual void
	read_def(
		std::istream &,
		Pose const &,
		FuncFactory const &
	);

	virtual void show_def( std::ostream &, Pose const & ) const;

	virtual void show( std::ostream & out ) const;

	core::Size seqpos() const { return seqpos_; }

	void set_sequence_profile( SequenceProfileOP );
	SequenceProfileOP sequence_profile();
	SequenceProfileCOP sequence_profile() const;

	virtual core::Size natoms() const;
	virtual ConstraintOP remap_resid( SequenceMapping const & ) const;
	virtual AtomID const & atom( core::Size const ) const;
	utility::vector1< AtomID > const & atom_ids() const;

	virtual void
	score(
		XYZ_Func const &,
		EnergyMap const &,
		EnergyMap &
	) const;

	virtual void
	fill_f1_f2(
		AtomID const &,
		XYZ_Func const &,
		core::Vector &,
		core::Vector &,
		EnergyMap const &
	) const;

private:
	core::Size seqpos_;
	SequenceProfileOP sequence_profile_;
	utility::vector1< AtomID > atom_ids_;
};


} // namespace constraints_additional
} // namespace protocols

#endif // INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_HH
