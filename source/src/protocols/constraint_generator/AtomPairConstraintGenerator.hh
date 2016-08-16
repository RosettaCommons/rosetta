// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/AtomPairConstraintGenerator.hh
/// @brief Generates atom pair constraints for a set of residues from the current or reference pose
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/AtomPairConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace constraint_generator {

///@brief Generates atom pair constraints for a set of residues from the current or reference pose
class AtomPairConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	struct MappedAtom {
	public:
		MappedAtom( core::Size const pose_atomno, core::Size const ref_atomno ):
			pose_atom( pose_atomno ), ref_atom( ref_atomno ) {}
		core::Size pose_atom;
		core::Size ref_atom;
	private:
		MappedAtom();
	};
	typedef utility::vector1< MappedAtom > MappedAtoms;

public:
	AtomPairConstraintGenerator();

	virtual ~AtomPairConstraintGenerator();

	static std::string
	class_name() { return "AtomPairConstraintGenerator"; }

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	void
	set_residue_selector( core::select::residue_selector::ResidueSelector const & selector );

	void
	set_sd( core::Real const sd );

	void
	set_ca_only( bool const ca_only );

	void
	set_max_distance( core::Real const max_dist );

	void
	set_min_seq_sep( core::Size const min_seq_sep );

	void
	set_weight( core::Real const weight );

	void
	set_reference_pose( core::pose::PoseCOP ref_pose );

private:
	core::scoring::constraints::ConstraintCOPs
	generate_constraints( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_atom_pair_constraints(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & subset ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_atom_pair_constraints(
		core::pose::Pose const & pose,
		core::pose::Pose const & ref_pose,
		core::select::residue_selector::ResidueSubset const & subset ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_atom_pair_constraints(
		core::pose::Pose const & pose,
		core::pose::Pose const & ref_pose,
		core::select::residue_selector::ResidueSubset const & subset,
		core::id::SequenceMapping const & seqmap ) const;

	core::id::SequenceMapping
	create_sequence_mapping( core::pose::Pose const & pose, core::pose::Pose const & ref_pose ) const;

	MappedAtoms
	atoms_to_constrain(
		core::conformation::Residue const & pose_rsd,
		core::conformation::Residue const & ref_rsd ) const;

	void
	add_constraints(
		core::scoring::constraints::ConstraintCOPs & csts,
		core::Size const pose_resid1,
		core::Size const pose_resid2,
		core::conformation::Residue const & ref_ires,
		core::conformation::Residue const & ref_jres,
		MappedAtoms const & iatoms,
		MappedAtoms const & jatoms ) const;

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::pose::PoseCOP reference_pose_;

	core::Real sd_;
	core::Real weight_;

	bool ca_only_;

	core::Real max_distance_;
	core::Size min_seq_sep_;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_AtomPairConstraintGenerator_hh

