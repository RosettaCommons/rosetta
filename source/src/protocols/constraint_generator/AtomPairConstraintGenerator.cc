// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/AtomPairConstraintGenerator.cc
/// @brief Generates atom pair constraints for a set of residues from the current or reference pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/AtomPairConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.AtomPairConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
AtomPairConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new AtomPairConstraintGenerator );
}

std::string
AtomPairConstraintGeneratorCreator::keyname() const
{
	return AtomPairConstraintGenerator::class_name();
}

AtomPairConstraintGenerator::AtomPairConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( AtomPairConstraintGenerator::class_name() ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	reference_pose_(),
	sd_( 0.5 ),
	weight_( 1.0 ),
	ca_only_( true ),
	max_distance_( 12.0 ),
	min_seq_sep_( 8 )
{
}

AtomPairConstraintGenerator::~AtomPairConstraintGenerator()
{}

protocols::constraint_generator::ConstraintGeneratorOP
AtomPairConstraintGenerator::clone() const
{
	return AtomPairConstraintGeneratorOP( new AtomPairConstraintGenerator( *this ) );
}

void
AtomPairConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	bool const use_native = tag->getOption< bool >( "native", false );
	if ( use_native ) {
		set_reference_pose( get_native_pose() );
		if ( ! reference_pose_ ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "'native' option for AtomPairConstraintGenerator specified, but no native pose is availible." );
		}
	}

	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( selector );

	set_sd( tag->getOption< core::Real >( "sd", sd_ ) );
	set_weight( tag->getOption< core::Real >( "weight", weight_ ) );
	set_ca_only( tag->getOption< bool >( "ca_only", ca_only_ ) );
	set_max_distance( tag->getOption< core::Real >( "max_distance", max_distance_ ) );
	set_min_seq_sep( tag->getOption< core::Size >( "min_seq_sep", min_seq_sep_ ) );

	if ( !selector_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "AtomPairConstraintGenerator requires a residue selector, but one is not set.\n" );
	}
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	debug_assert( selector_ );
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	return generate_atom_pair_constraints( pose, subset );
}

void
AtomPairConstraintGenerator::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

void
AtomPairConstraintGenerator::set_sd( core::Real const sd )
{
	sd_ = sd;
}

void
AtomPairConstraintGenerator::set_ca_only( bool const ca_only )
{
	ca_only_ = ca_only;
}

void
AtomPairConstraintGenerator::set_max_distance( core::Real const max_dist )
{
	max_distance_ = max_dist;
}

void
AtomPairConstraintGenerator::set_min_seq_sep( core::Size const min_seq_sep )
{
	min_seq_sep_ = min_seq_sep;
}

void
AtomPairConstraintGenerator::set_weight( core::Real const weight )
{
	weight_ = weight;
}

void
AtomPairConstraintGenerator::set_reference_pose( core::pose::PoseCOP ref_pose )
{
	reference_pose_ = ref_pose;
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	if ( reference_pose_ ) {
		return generate_atom_pair_constraints( pose, *reference_pose_, subset );
	} else {
		return generate_atom_pair_constraints( pose, pose, subset );
	}
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	core::id::SequenceMapping const seqmap = create_sequence_mapping( pose, ref_pose );
	return generate_atom_pair_constraints( pose, ref_pose, subset, seqmap );
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::id::SequenceMapping const & seqmap ) const
{
	core::scoring::constraints::ConstraintCOPs csts;

	core::Size const nres = compute_nres_in_asymmetric_unit( pose );
	for ( Size ires=1; ires<=nres; ++ires ) {
		if ( !subset[ ires ] ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
		core::Size const ref_ires = seqmap[ ires ];
		if ( ref_ires == 0 ) {
			TR.Debug << "Residue " << ires << " not found in reference pose, skipping" << std::endl;
			continue;
		}
		MappedAtoms const iatoms = atoms_to_constrain( pose.residue( ires ), ref_pose.residue( ref_ires ) );

		for ( Size jres=ires+min_seq_sep_; jres<=pose.total_residue(); ++jres ) {
			if ( !subset[ jres ] ) continue;
			if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
			core::Size const ref_jres = seqmap[ jres ];
			if ( ref_jres == 0 ) {
				TR.Debug << "Residue " << jres << " not found in reference pose, skipping" << std::endl;
				continue;
			}

			MappedAtoms const jatoms = atoms_to_constrain( pose.residue( jres ), ref_pose.residue( ref_jres ) );
			add_constraints( csts, ires, jres, ref_pose.residue( ref_ires ), ref_pose.residue( ref_jres ), iatoms, jatoms );
		} // jres loop
	} // ires loop
	return csts;
}

void
AtomPairConstraintGenerator::add_constraints(
	core::scoring::constraints::ConstraintCOPs & csts,
	core::Size const pose_resid1,
	core::Size const pose_resid2,
	core::conformation::Residue const & ref_ires,
	core::conformation::Residue const & ref_jres,
	MappedAtoms const & iatoms,
	MappedAtoms const & jatoms ) const
{
	for ( MappedAtoms::const_iterator iatom=iatoms.begin(); iatom!=iatoms.end(); ++iatom ) {
		for ( MappedAtoms::const_iterator jatom=jatoms.begin(); jatom!=jatoms.end(); ++jatom ) {
			core::Real const dist = ref_ires.xyz( iatom->ref_atom ).distance( ref_jres.xyz( jatom->ref_atom ) );
			if ( dist > max_distance_ ) continue;

			core::id::AtomID const pose_atom1( iatom->pose_atom, pose_resid1 );
			core::id::AtomID const pose_atom2( jatom->pose_atom, pose_resid2 );
			core::scoring::func::FuncOP sog_func( new core::scoring::func::SOGFunc( dist, sd_ ) );
			core::scoring::func::FuncOP weighted_func = scalar_weighted( sog_func, weight_ );
			core::scoring::constraints::ConstraintCOP newcst(
				new core::scoring::constraints::AtomPairConstraint( pose_atom1, pose_atom2, weighted_func ) );
			csts.push_back( newcst );
			TR.Debug << "atom_pair_constraint generated for residue " << pose_atom1 << " and " << pose_atom2 << ", distance=" << dist << " with weight " << weight_ << std::endl;
		}
	}
}

core::id::SequenceMapping
AtomPairConstraintGenerator::create_sequence_mapping( core::pose::Pose const & pose, core::pose::Pose const & ref_pose ) const
{
	bool const same_length = ( pose.total_residue() == ref_pose.total_residue() );
	bool const same_sequence = ( pose.sequence() == ref_pose.sequence() );

	if ( same_length && same_sequence ) {
		return core::id::SequenceMapping::identity( pose.total_residue() );
	} else { // !same_sequence || !same_length
		TR << "Input structure and native differ in ";
		if ( !same_length ) TR << "length and sequence ";
		else if ( !same_sequence ) TR << "sequence ";
		TR << "- aligning on PDB identity or sequence." << std::endl;
		return core::pose::sequence_map_from_pdbinfo( pose, ref_pose );
	}
}

AtomPairConstraintGenerator::MappedAtoms
AtomPairConstraintGenerator::atoms_to_constrain(
	core::conformation::Residue const & pose_rsd,
	core::conformation::Residue const & ref_rsd ) const
{
	utility::vector1< MappedAtom > atoms;
	if ( pose_rsd.has( "CA" ) && ref_rsd.has( "CA" ) ) {
		atoms.push_back( MappedAtom( pose_rsd.type().atom_index( "CA" ), ref_rsd.type().atom_index( "CA" ) ) );
	}
	if ( !ca_only_ ) {
		if ( ( pose_rsd.nbr_atom() != atoms[1].pose_atom ) && ( ref_rsd.nbr_atom() != atoms[1].ref_atom ) ) {
			atoms.push_back( MappedAtom( pose_rsd.nbr_atom(), ref_rsd.nbr_atom() ) );
		}
	}
	return atoms;
}

} //protocols
} //constraint_generator

