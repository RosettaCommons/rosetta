// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/CoordinateConstraintGenerator.cc
/// @brief Generates coodinate constraints
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/constraint_generator/CoordinateConstraintGenerator.hh>
#include <protocols/constraint_generator/CoordinateConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/sequence/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.CoordinateConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
CoordinateConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new CoordinateConstraintGenerator );
}

std::string
CoordinateConstraintGeneratorCreator::keyname() const
{
	return CoordinateConstraintGenerator::class_name();
}

CoordinateConstraintGenerator::CoordinateConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( CoordinateConstraintGenerator::class_name() ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	refpose_(),
	bounded_( false ),
	bounded_width_( 0.0 ),
	sd_( 0.5 ),
	ca_only_( false ),
	sidechain_( false ),
	ambiguous_hnq_( false ),
	align_reference_( false )
{
}

CoordinateConstraintGenerator::~CoordinateConstraintGenerator() = default;

ConstraintGeneratorOP
CoordinateConstraintGenerator::clone() const
{
	return ConstraintGeneratorOP( new CoordinateConstraintGenerator( *this ) );
}

void
CoordinateConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	// ref pose
	bool const use_native = tag->getOption< bool >( "native", false );
	if ( use_native ) {
		set_reference_pose( get_native_pose() );
		if ( ! refpose_ ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "'native' option for CoordinateConstraintGenerator specified, but no native pose is availible." );
		}
	}

	set_ca_only( tag->getOption< bool >( "ca_only", ca_only_ ) );
	set_sidechain( tag->getOption< bool >( "sidechain", sidechain_ ) );
	set_sd( tag->getOption< core::Real >( "sd", sd_ ) );
	set_ambiguous_hnq( tag->getOption< bool >( "ambiguous_hnq", ambiguous_hnq_ ) );
	set_bounded( tag->getOption< bool >( "bounded", bounded_ ) );
	set_bounded_width( tag->getOption< core::Real >( "bounded_width", bounded_width_ ) );
	align_reference_ = tag->getOption< bool >( "align_reference", align_reference_ );

	core::select::residue_selector::ResidueSelectorCOP selector =
		core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( selector );

	if ( !selector_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "CoordinateConstraintGenerator::parse_tag(): Error obtaining ResidueSelector from tag\n" );
	}
}

core::scoring::constraints::ConstraintCOPs
CoordinateConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring::constraints;

	// get sequence map
	core::id::SequenceMapping const seq_map = generate_seqmap( pose );
	TR.Debug << "Generated sequence map: " << seq_map << std::endl;

	// get reference pose
	core::pose::PoseOP reference_pose = prepare_constraint_target_pose( pose, seq_map );

	// get residue subset
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );

	// root atom id
	core::id::AtomID const root_atomid( 1, pose.fold_tree().root() );

	core::scoring::constraints::ConstraintCOPs csts;
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( !subset[ resid ] ) continue;

		// don't constrain virtual root residue to itself
		if ( ( resid == root_atomid.rsd() ) && ( pose.residue( resid ).aa() == core::chemical::aa_vrt ) ) {
			continue;
		}

		core::Size const ref_resid = seq_map[ resid ];
		if ( ref_resid == 0 ) {
			TR.Debug << "apply(): skipping residue " << resid << " because it was not found in the sequence map." << std::endl;
			continue;
		}
		debug_assert( ref_resid <= reference_pose->size() );

		create_residue_constraints( csts, root_atomid, pose.residue( resid ), reference_pose->residue( ref_resid ) );
	}
	return csts;
}

void
CoordinateConstraintGenerator::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

void
CoordinateConstraintGenerator::set_reference_pose( core::pose::PoseCOP refpose )
{
	refpose_ = refpose;
}

void
CoordinateConstraintGenerator::set_sidechain( bool const sidechain )
{
	sidechain_ = sidechain;
}

void
CoordinateConstraintGenerator::set_ca_only( bool const ca_only )
{
	ca_only_ = ca_only;
}

void
CoordinateConstraintGenerator::set_ambiguous_hnq( bool const ambiguous_hnq )
{
	ambiguous_hnq_ = ambiguous_hnq;
}

void
CoordinateConstraintGenerator::set_sd( core::Real const sd )
{
	sd_ = sd;
}

void
CoordinateConstraintGenerator::set_bounded( bool const bounded )
{
	bounded_ = bounded;
}

void
CoordinateConstraintGenerator::set_bounded_width( core::Real const bound_width )
{
	bounded_width_ = bound_width;
}

std::pair< std::string, std::string >
compute_hnq_atoms( core::conformation::Residue const & pose_rsd, std::string const & atom_name )
{
	if ( pose_rsd.aa() == core::chemical::aa_asn ) {
		return std::make_pair( " OD1", " ND2" );
	} else if ( pose_rsd.aa() == core::chemical::aa_gln ) {
		return std::make_pair( " OE1", " NE2" );
	} else if ( pose_rsd.aa() == core::chemical::aa_his &&
			(atom_name == " ND1" || atom_name == " CD2") ) {
		return std::make_pair( " ND1", " CD2" );
	} else if ( pose_rsd.aa() == core::chemical::aa_his &&
			(atom_name == " NE2" || atom_name == " CE1") ) {
		return std::make_pair( " NE2", " CE1" );
	} else {
		utility_exit_with_message("Logic error in AtomCoordinateConstraints");
		return std::make_pair( "", "" );
	}
}

core::scoring::func::FuncOP
CoordinateConstraintGenerator::create_function() const
{
	if ( bounded_ ) {
		return core::scoring::func::FuncOP( new core::scoring::constraints::BoundFunc( 0, bounded_width_, sd_, "xyz" ) );
	} else {
		return core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( 0.0, sd_ ) );
	}
}

core::scoring::constraints::ConstraintOP
CoordinateConstraintGenerator::create_coordinate_constraint(
	core::id::AtomID const & pose_atomid,
	core::id::AtomID const & ref_atom,
	core::id::AtomID const & root_atom,
	core::conformation::Residue const & ref_rsd ) const
{
	return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint( pose_atomid, root_atom, ref_rsd.xyz( ref_atom.atomno() ), create_function() ) );
}

core::scoring::constraints::ConstraintOP
CoordinateConstraintGenerator::create_ambiguous_constraint(
	core::id::AtomID const & pose_atomid,
	core::conformation::Residue const & ref_rsd,
	core::id::AtomID const & root_atomid,
	std::pair< std::string, std::string > const & hnq_atoms ) const
{
	core::scoring::constraints::AmbiguousConstraintOP amb_constr( new core::scoring::constraints::AmbiguousConstraint );
	core::id::AtomID const atomid1( ref_rsd.atom_index( hnq_atoms.first ), ref_rsd.seqpos() );
	amb_constr->add_individual_constraint( create_coordinate_constraint( pose_atomid, atomid1, root_atomid, ref_rsd ) );
	core::id::AtomID const atomid2( ref_rsd.atom_index( hnq_atoms.second ), ref_rsd.seqpos() );
	amb_constr->add_individual_constraint( create_coordinate_constraint( pose_atomid, atomid2, root_atomid, ref_rsd ) );
	return amb_constr;
}

core::Size
compute_ref_atom(
	std::string const & atomname,
	core::conformation::Residue const & pose_rsd,
	core::conformation::Residue const & ref_rsd )
{
	if ( ! ref_rsd.has( atomname ) ) {
		TR.Debug << "Skip generating coordinate constraints for atom " << atomname << " of residue " << pose_rsd.seqpos() << " (" << pose_rsd.name() <<
			") - not found in residue " << ref_rsd.seqpos() << " (" << ref_rsd.name() << ") of reference structure." << std::endl;
		return 0;
	}
	return ref_rsd.atom_index( atomname );
}

void
CoordinateConstraintGenerator::create_residue_constraints(
	core::scoring::constraints::ConstraintCOPs & csts,
	core::id::AtomID const & root_atomid,
	core::conformation::Residue const & pose_rsd,
	core::conformation::Residue const & ref_rsd ) const
{
	core::Size last_atom = pose_rsd.last_backbone_atom();

	if ( sidechain_ ) {
		last_atom = pose_rsd.nheavyatoms();
	}

	for ( core::Size atom=1; atom<=last_atom; ++atom ) {
		std::string const & atom_name = pose_rsd.atom_name( atom );
		TR.Debug << "Considering constraining atom " << atom_name << " residue " << pose_rsd.seqpos() << std::endl;

		// skip ahead if we only want CA
		if ( ca_only_ && ( atom_name.find( "CA" ) == std::string::npos ) ) continue;
		TR.Debug << "Going to constrain atom " << atom_name << " residue " << pose_rsd.seqpos() << std::endl;

		core::Size const ref_atom = compute_ref_atom( atom_name, pose_rsd, ref_rsd );
		if ( ref_atom == 0 ) continue;

		core::id::AtomID const ref_atomid( ref_atom, ref_rsd.seqpos() );
		core::id::AtomID const pose_atomid( atom, pose_rsd.seqpos() );

		// Rely on shortcutting evaluation to speed things up - get to else clause as soon as possible.
		if ( ambiguous_hnq_ && sidechain_ &&
				( (pose_rsd.aa() == core::chemical::aa_asn && ref_rsd.aa() == core::chemical::aa_asn &&
				( atom_name == " OD1" || atom_name == " ND2" ) ) ||
				(pose_rsd.aa() == core::chemical::aa_gln && ref_rsd.aa() == core::chemical::aa_gln &&
				( atom_name == " OE1" || atom_name == " NE2" )) ||
				(pose_rsd.aa() == core::chemical::aa_his && ref_rsd.aa() == core::chemical::aa_his &&
				(atom_name == " ND1" || atom_name == " NE2" ||
				atom_name == " CD2" || atom_name == " CE1") ) ) ) {
			csts.push_back( create_ambiguous_constraint( pose_atomid, ref_rsd, root_atomid, compute_hnq_atoms( pose_rsd, atom_name ) ) );
		} else {
			csts.push_back( create_coordinate_constraint( pose_atomid, ref_atomid, root_atomid, ref_rsd ) );
		}
	}
}

core::pose::Pose const &
CoordinateConstraintGenerator::select_constraint_target_pose( core::pose::Pose const & pose ) const
{
	if ( refpose_ ) {
		return *refpose_;
	} else { // !refpose
		return pose;
	}
}

core::id::SequenceMapping
CoordinateConstraintGenerator::generate_seqmap( core::pose::Pose const & pose ) const
{
	return generate_seqmap_from_poses( pose, select_constraint_target_pose( pose ) );
}

core::pose::PoseOP
CoordinateConstraintGenerator::prepare_constraint_target_pose(
	core::pose::Pose const & input,
	core::id::SequenceMapping const & seq_map ) const
{
	core::pose::PoseOP constraint_target_pose = select_constraint_target_pose( input ).clone();
	debug_assert( constraint_target_pose );

	bool const has_virtual_root = ( input.residue( input.fold_tree().root() ).aa() == core::chemical::aa_vrt );

	// If a ref pose is given, and pose doens't have a virtual root,
	// align the native pose to the input pose to avoid rotation/translation based
	//  errors.
	if ( refpose_ && ( align_reference_ || !has_virtual_root ) ) {
		align_reference_pose( *constraint_target_pose, input, seq_map );
	}

	// Warn about not having a virtual root (but go ahead with constraints).
	if ( !has_virtual_root ) {
		TR.Warning << "WARNING: Adding coordinate constraints to a pose without a virtual root - results may not be as expected." << std::endl;
	}

	return constraint_target_pose;
}

void
CoordinateConstraintGenerator::align_reference_pose(
	core::pose::Pose & reference,
	core::pose::Pose const & input_pose,
	core::id::SequenceMapping const & seq_map ) const
{
	core::id::SequenceMapping rev_seq_map( seq_map ); // constraint_target_pose -> pose mapping
	rev_seq_map.reverse();
	core::sequence::calpha_superimpose_with_mapping( reference, input_pose, rev_seq_map );
}

void
CoordinateConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CoordinateConstraintGenerator::provide_xml_schema( xsd );
}

void
CoordinateConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist

		+ XMLSchemaAttribute::attribute_w_default( "native", xsct_rosetta_bool, "Constrain to native coordinates", "false" )
		+ XMLSchemaAttribute( "ca_only", xsct_rosetta_bool, "Only apply constraints to alpha carbons" )
		+ XMLSchemaAttribute( "sidechain", xsct_rosetta_bool, "Apply constraints to side chain atoms" )
		+ XMLSchemaAttribute( "sd", xsct_real, "Standard deviation for constraints" )
		+ XMLSchemaAttribute( "ambiguous_hnq", xsct_rosetta_bool, "If true, these residues could be flipped without affecting the score" )
		+ XMLSchemaAttribute( "bounded", xsct_rosetta_bool, "Should we use a bounded function for this constraint?" )
		+ XMLSchemaAttribute( "bounded_width", xsct_real, "Width of bounded function if using" )
		+ XMLSchemaAttribute( "align_reference", xsct_rosetta_bool, "Align reference pose to constraint target pose?" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Selector specifying residues to constrain" );
	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Adds coordinate constraints to the specified residues",
		attlist );
}

} //protocols
} //constraint_generator






