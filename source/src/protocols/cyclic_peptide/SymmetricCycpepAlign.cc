// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/SymmetricCycpepAlign.cc
/// @brief Given a quasi-symmetric cyclic peptide, this mover aligns the peptide so that the cyclic symmetry axis lies along the Z-axis and the centre of mass is at the origin.
/// It then optionally removes all but one symmetry repeat, so that true symmetry may be set up with the SetupForSymmetry mover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/SymmetricCycpepAlign.hh>
#include <protocols/cyclic_peptide/SymmetricCycpepAlignCreator.hh>

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>

// Protocols headers
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>

// Numeric headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.SymmetricCycpepAlign" );

namespace protocols {
namespace cyclic_peptide {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SymmetricCycpepAlign::SymmetricCycpepAlign():
	protocols::moves::Mover( SymmetricCycpepAlign::mover_name() ),
	symmetry_repeats_(2),
	mirror_symmetry_(false),
	auto_detect_symmetry_(false),
	angle_threshold_(10.0),
	trim_to_single_repeat_(false),
	repeat_to_preserve_(1),
	last_symmetry_repeats_(0),
	last_symmetry_mirror_(false),
	invert_(false)
	//TODO initialize vars here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SymmetricCycpepAlign::SymmetricCycpepAlign( SymmetricCycpepAlign const & src ):
	protocols::moves::Mover( src ),
	symmetry_repeats_(src.symmetry_repeats_),
	mirror_symmetry_(src.mirror_symmetry_),
	auto_detect_symmetry_(src.auto_detect_symmetry_),
	angle_threshold_(src.angle_threshold_),
	trim_to_single_repeat_(src.trim_to_single_repeat_),
	repeat_to_preserve_(src.repeat_to_preserve_),
	last_symmetry_repeats_(src.last_symmetry_repeats_),
	last_symmetry_mirror_(src.last_symmetry_mirror_),
	invert_(src.invert_)
	//TODO copy vars here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SymmetricCycpepAlign::~SymmetricCycpepAlign(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
SymmetricCycpepAlign::apply(
	core::pose::Pose &pose )
{
	core::Size symmrepeats( symmetry_repeats() );
	bool mirrorsymm( mirror_symmetry() );
	if ( auto_detect_symmetry() ) {
		if ( !do_auto_detection_of_symmetry( pose, symmrepeats, mirrorsymm ) ) {
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			return;
		}
		set_last_symmetry( symmrepeats, mirrorsymm );
	} else {
		if ( !do_symmetry_checks( pose, symmrepeats, mirrorsymm ) ) {
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			return;
		}
	}

	align_to_origin( pose );
	align_to_zaxis( pose, symmrepeats, mirrorsymm );
	if ( trim_to_single_repeat() ) {
		do_trim_to_single_repeat( pose, symmrepeats, repeat_to_preserve() );
	}

	set_last_move_status( protocols::moves::MS_SUCCESS );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
SymmetricCycpepAlign::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SymmetricCycpepAlign::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	bool auto_detect( tag->getOption<bool>("auto_detect_symmetry", "false") );
	runtime_assert_string_msg(
		!auto_detect || ( !tag->hasOption("symmetry_repeats") && !tag->hasOption("mirror_symmetry") ),
		"Error in protocols::cyclic_peptide::SymmetricCycpepAlign::parse_my_tag(): If the auto_detect_symmetry=\"true\" option is used, then the \"mirror_symmetry\" and \"symmetry_repeats\" options cannot be used."
	);

	set_auto_detect_symmetry(auto_detect);
	set_symmetry(
		tag->getOption<core::Size>( "symmetry_repeats", symmetry_repeats() ),
		tag->getOption<core::Size>( "mirror_symmetry", mirror_symmetry() )
	);
	set_angle_threshold( tag->getOption<core::Real>("angle_threshold", angle_threshold() ) );
	set_trim_info(
		tag->getOption<bool>( "trim_to_single_repeat", trim_to_single_repeat() ),
		tag->getOption<core::Size>( "repeat_to_preserve", repeat_to_preserve() )
	);
	set_invert( tag->getOption<bool>( "invert", invert() ) );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
SymmetricCycpepAlign::fresh_instance() const
{
	return protocols::moves::MoverOP( new SymmetricCycpepAlign );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SymmetricCycpepAlign::clone() const
{
	return protocols::moves::MoverOP( new SymmetricCycpepAlign( *this ) );
}

std::string
SymmetricCycpepAlign::get_name() const {
	return mover_name();
}

std::string
SymmetricCycpepAlign::mover_name() {
	return "SymmetricCycpepAlign";
}

void
SymmetricCycpepAlign::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "auto_detect_symmetry", xsct_rosetta_bool, "If true, the symmetry of this peptide will be detected automatically.  False by default.  (Note that if this is set to \"true\", the \"symmetry_repeats\" and \"mirror_symmetry\" options cannot be used.)" )
		+ XMLSchemaAttribute( "symmetry_repeats", xsct_non_negative_integer, "The number of symmetry repeats.  For example, to specify c6 or c6/m symmetry, set this to \"6\".  Defaults to \"2\" (for c2 symmetry)." )
		+ XMLSchemaAttribute( "mirror_symmetry", xsct_rosetta_bool, "If true, this indicates that the pose possesses mirror symmetry.  For example, to specify c6/m symmetry, set \"symmetry_repeats\" to \"6\" and \"mirror_symmetry\" to \"true\".  If set to \"false\", this would indicate c6 symmetry.  Defaults to \"false\" (for c2 symmetry with no mirroring)." )
		+ XMLSchemaAttribute( "angle_threshold", xsct_real, "The cutoff, in degrees, for considering two dihedral values to be equal.  Defaults to 10 degrees.  This is used when confirming that a quasi-symmetric peptide has the indicated symmetry, or for detecting the symmetry of the peptide." )
		+ XMLSchemaAttribute( "trim_to_single_repeat", xsct_rosetta_bool, "If true, all geometry in the peptide (including all crosslinkers) will be deleted, except for the polypeptide part of a single symmetry repeat.  Defaults to \"false\".  (This is useful for setting up a peptide for the SetupForSymmetry mover)." )
		+ XMLSchemaAttribute( "repeat_to_preserve", xsct_non_negative_integer, "If \"trim_to_single_repeat\" is true, then this is the symmetry repeat that will NOT be deleted (i.e. the one that will be preserved).  Defaults to \"1\"." )
		+ XMLSchemaAttribute( "invert", xsct_rosetta_bool, "If true, the peptide normal is aligned to the negative Z-axis; if false, it is aligned to the positive.  Default false." )
		;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Given a quasi-symmetric cyclic peptide, this mover aligns the peptide so that the cyclic symmetry axis lies along the Z-axis and the centre of mass is at the origin.  It then optionally removes all but one symmetry repeat, so that true symmetry may be set up with the SetupForSymmetry mover.", attlist );
}

/// @brief Set the number of symmetry repeats and whether we're using mirror symmetry.  For example, for c4
/// symmetry, inputs are "4", "false".  For c4/m, they'd be "4", "true".
void
SymmetricCycpepAlign::set_symmetry(
	core::Size const repeats_in,
	bool const mirror_in
) {
	runtime_assert_string_msg( repeats_in > 1, "Error in protocols::cyclic_peptide::SymmetricCycpepAlign::set_symmetry(): The number of symmetry repeats must be greater than 1." );
	symmetry_repeats_ = repeats_in;
	mirror_symmetry_ = mirror_in;
}

/// @brief Set whether we're auto-detecting symmetry.
///
void
SymmetricCycpepAlign::set_auto_detect_symmetry(
	bool const setting
) {
	auto_detect_symmetry_ = setting;
}

/// @brief Set the angle threshold for counting poses as symmetric.
/// @details Two dihedral values from different residues must fall within this cutoff, in degrees, for them to
/// be considered the "same".  Default 10 degrees.
void
SymmetricCycpepAlign::set_angle_threshold(
	core::Real const &setting
) {
	runtime_assert_string_msg( setting > 0.0, "Error in protocols::cyclic_peptide::SymmetricCycpepAlign::set_angle_threshold(): The angle threshold must be greater than zero." );
	angle_threshold_ = setting;
}

/// @brief Set whether we're going to delete all geometry that isn't the protein backbone of a single
/// symmetry repeat, and which repeat to preserve.
/// @param[in] do_trim If set to true, the peptide will be trimmed down to a single repeat.  False by default (no trimming).
/// @param[in] repeat_to_preserve If do_trim is set to true, this is the repeat that should be preserved.  1 by default.
void
SymmetricCycpepAlign::set_trim_info(
	bool const do_trim,
	core::Size const repeat_to_preserve/*=1*/
) {
	trim_to_single_repeat_ = do_trim;
	runtime_assert_string_msg( repeat_to_preserve > 0, "Error in protocols::cyclic_peptide::SymmetricCycpepAlign::set_trim_info(): The index of the symmetry repeat to preserve must be greater than zero." );
	repeat_to_preserve_ = repeat_to_preserve;
}

/////////////// Private methods ///////

/// @brief Set the symmetry of the last auto-detected peptide symmetry.
///
void
SymmetricCycpepAlign::set_last_symmetry(
	core::Size const last_symm_in,
	bool const last_symm_mirror_in
) {
	last_symmetry_repeats_ = last_symm_in;
	last_symmetry_mirror_ = last_symm_mirror_in;
}

/// @brief Given a quasi-symmetric pose, figure out its symmetry.
/// @details Starts with the maximum possible, and tries every possible symmetry type, favouring mirror symmetry over non-mirror symmetry.
/// @param[in] pose The quasi-symmetric pose.
/// @param[out] symmrepeats The number of symmetry repeats.
/// @param[out] mirrorsymm Is this a mirror symmetry type (cN/m) or not (cN)?
bool
SymmetricCycpepAlign::do_auto_detection_of_symmetry(
	core::pose::Pose const &pose,
	core::Size &symmrepeats,
	bool &mirrorsymm
) const {
	core::Size const protrescount( count_protein_residues( pose ) );
	core::select::residue_selector::ResidueSelectorCOP selector( select_protein_residues( pose ) );

	CycpepSymmetryFilter symfilt;
	symfilt.set_angle_threshold( angle_threshold() );
	symfilt.set_selector( selector );

	bool found(false);
	for ( core::Size i( protrescount ); i>=2; --i ) { //Count backwards from the number of residues in the peptide.
		if ( protrescount % i != 0 ) continue;
		symfilt.set_symm_repeats(i);
		if ( i % 2 == 0 ) {
			symfilt.set_mirror_symm( true );
			if ( symfilt.apply( pose ) ) {
				found=true;
				break;
			}
		}
		symfilt.set_mirror_symm( false );
		if ( symfilt.apply( pose ) ) {
			found=true;
			break;
		}
	}

	if ( !found ) {
		TR.Warning << "Warning: The SymmetricCycpepAlign mover attempted to auto-detect the symmetry of the peptide, but the peptide was not at all symmetric!" << std::endl;
	} else {
		symmrepeats = symfilt.symm_repeats();
		mirrorsymm = symfilt.mirror_symm();
		TR << "The peptide was found to have c" << symmrepeats << (mirrorsymm ? "/m " : " ") << " symmetry." << std::endl;
	}

	return found;
}

/// @brief Given a quasi-symmetric pose, confirm that it has the specified symmetry.
///
bool
SymmetricCycpepAlign::do_symmetry_checks(
	core::pose::Pose const &pose,
	core::Size const symmrepeats,
	bool const mirrorsymm
) const {
	std::stringstream symmtype;
	symmtype << "c" << symmrepeats;
	if ( mirrorsymm ) symmtype << "/m";


	core::Size const protrescount( count_protein_residues( pose ) );
	core::select::residue_selector::ResidueSelectorCOP selector( select_protein_residues( pose ) );

	if ( protrescount % symmrepeats != 0 ) {
		TR.Warning << "Warning: The number of residues in the peptide is incompatible with " + symmtype.str() + " symmetry." << std::endl;
		return false;
	}
	CycpepSymmetryFilter symfilt;
	symfilt.set_angle_threshold( angle_threshold() );
	symfilt.set_selector( selector );
	symfilt.set_symm_repeats( symmrepeats );
	symfilt.set_mirror_symm( mirrorsymm );
	if ( !symfilt.apply(pose) ) {
		TR.Warning << "Warning: The peptide does not have " + symmtype.str() + " symmetry." << std::endl;
		return false;
	}

	TR << "Confirmed that the peptide has " << symmtype.str() << " symmetry." << std::endl;
	return true;
}

/// @brief Given a pose, count the protein residues in it.
///
core::Size
SymmetricCycpepAlign::count_protein_residues(
	core::pose::Pose const &pose
) const {
	core::Size const nres( pose.total_residue() );
	if ( nres == 0 ) return 0;

	core::Size counter(0);
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		if ( pose.residue_type(ir).is_protein() ) ++counter;
	}
	return counter;
}

/// @brief Given a pose, return a ResidueSelector that selects all the protein residues in the pose.
/// @details Uses Jared's ReturnResidueSubetSelector, which was a good idea.  Thanks, Jared!
core::select::residue_selector::ResidueSelectorCOP
SymmetricCycpepAlign::select_protein_residues(
	core::pose::Pose const &pose
) const {
	using namespace core::select::residue_selector;

	core::Size const nres( pose.total_residue() );

	ResidueSubset subset( nres );
	for ( core::Size i=1; i<=nres; ++i ) {
		subset[i] = pose.residue_type(i).is_protein();
	}
	ReturnResidueSubsetSelectorOP selector( new ReturnResidueSubsetSelector );
	selector->set_residue_subset( subset );

	return ResidueSelectorCOP( selector );
}

/// @brief Given a pose, center it on the origin.
///
void
SymmetricCycpepAlign::align_to_origin( core::pose::Pose &pose ) const {
	core::Size atcounter(0);
	numeric::xyzVector<core::Real> centvect;
	centvect.x(0); centvect.y(0); centvect.z(0);
	for ( core::Size ir=1, irmax=pose.total_residue(); ir<=irmax; ++ir ) {
		if ( pose.residue_type(ir).is_protein() ) {
			for ( core::Size ia=1, iamax=pose.residue_type(ir).first_sidechain_atom(); ia<iamax; ++ia ) {
				if ( pose.residue_type(ir).atom_is_hydrogen(ia) ) continue;
				centvect += pose.residue(ir).xyz(ia);
				++atcounter;
			}
		}
	}

	centvect /= static_cast<core::Real>( atcounter );
	centvect *= -1.0;

	numeric::xyzMatrix<core::Real> rotmatrix( numeric::xyzMatrix<core::Real>::identity() );

	pose.apply_transform_Rx_plus_v( rotmatrix, centvect );
}

/// @brief Given a pose centered on the origin, align it to the z-axis based on its symmetry.
///
void
SymmetricCycpepAlign::align_to_zaxis(
	core::pose::Pose &pose,
	core::Size const symmrepeats,
	bool const mirrorsymm
) const {
	if ( symmrepeats == 2 && mirrorsymm ) return; //Special case: any axis can be considered the symmetry axis in the c2/m symmetric case.

	core::Size const protein_residues( count_protein_residues(pose) ); //Get the number of protein residues in the peptide.

	core::Size counter(0), first_repeat_res(0), last_repeat_res(0), first_repeat_res_equiv(0), last_repeat_res_equiv(0);
	core::Size const last_res_expected(protein_residues / symmrepeats);
	core::Size const first_res_equiv_expected( mirrorsymm ? 2*last_res_expected + 1 : last_res_expected + 1 );
	core::Size const last_res_equiv_expected( mirrorsymm ? 3*last_res_expected : 2*last_res_expected );

	for ( core::Size ir=1, irmax=pose.total_residue(); ir<=irmax; ++ir ) {
		if ( !pose.residue_type(ir).is_protein() ) continue;
		++counter;
		if ( counter == 1 ) first_repeat_res = ir;
		if ( counter == last_res_expected ) last_repeat_res = ir;
		if ( counter == first_res_equiv_expected ) first_repeat_res_equiv = ir;
		if ( counter == last_res_equiv_expected ) last_repeat_res_equiv = ir;
	}
	debug_assert( first_repeat_res != 0 && last_repeat_res != 0 && first_repeat_res_equiv != 0 && last_repeat_res_equiv != 0 );

	numeric::xyzVector< core::Real > v1( pose.residue(first_repeat_res).xyz(1) - pose.residue(first_repeat_res_equiv).xyz(1) );
	numeric::xyzVector< core::Real > v2( pose.residue(last_repeat_res).xyz(1) - pose.residue(last_repeat_res_equiv).xyz(1) );
	v1.normalize(); v2.normalize();
	core::Real const dotpdt( v1.dot(v2) );
	if ( dotpdt > 0.999 && dotpdt < 1.001 ) { //Vectors too close to being parallel; choose new v2.
		v2 = pose.residue(last_repeat_res).xyz(2) - pose.residue(last_repeat_res_equiv).xyz(2);
		v2.normalize();
	}

	numeric::xyzVector< core::Real > const peptide_normal( v1.cross(v2).normalize() );
	numeric::xyzVector< core::Real > const zaxis( 0.0, 0.0, invert() ? -1.0 : 1.0 );

	core::Real const costheta( peptide_normal.dot(zaxis) );
	if ( costheta > 0.999999 && costheta < 1.000001 ) return; //Already aligned in this case.

	//The cross product of the peptide normal and the Z-axis is the vector about which we want to rotate.
	numeric::xyzVector< core::Real > const crosspdt( peptide_normal.cross(zaxis) );
	core::Real const sinsqtheta( crosspdt.length_squared() );

	//This matrix helps us to calculate the rotation matrix:
	numeric::xyzMatrix< core::Real > const vmatrix(
		numeric::xyzMatrix< core::Real>::rows(
		0.0, -1.0*crosspdt.z(), crosspdt.y(),
		crosspdt.z(), 0.0, -1.0*crosspdt.x(),
		-1.0*crosspdt.y(), crosspdt.x(), 0.0
		)
	);

	//Calculate the rotation matrix:
	numeric::xyzMatrix< core::Real > const rotmatrix(
		numeric::xyzMatrix<core::Real>::identity() +
		vmatrix +
		((1-costheta)/sinsqtheta)*(vmatrix*vmatrix)
	);

	//We won't be translating, since our pose is already aligned to the origin:
	numeric::xyzVector< core::Real > const transvect( 0.0, 0.0, 0.0);

	//Apply the transformation:
	pose.apply_transform_Rx_plus_v( rotmatrix, transvect );
}

/// @brief Given a pose that has been properly centered on the origin and aligned to the z-axis, delete all but
/// a single symmetry repeat.
void
SymmetricCycpepAlign::do_trim_to_single_repeat(
	core::pose::Pose &pose,
	core::Size const symmrepeats,
	core::Size const repeat_to_preserve
) const {
	debug_assert( repeat_to_preserve > 0 ); //Should be guaranteed true.
	core::pose::Pose newpose;

	core::Size const pose_rescount( count_protein_residues( pose ) );
	debug_assert( pose_rescount % symmrepeats == 0 ); //Should be guaranteed true.
	core::Size const residues_per_repeat( pose_rescount / symmrepeats );

	core::Size const start_res( (repeat_to_preserve - 1)*residues_per_repeat + 1 );
	core::Size const end_res( repeat_to_preserve*residues_per_repeat );

	core::Size counter(0);
	for ( core::Size ir=1, irmax=pose.total_residue(); ir<=irmax; ++ir ) {
		if ( !pose.residue_type(ir).is_protein() ) continue;
		++counter;

		if ( counter == start_res ) {
			newpose.append_residue_by_jump( pose.residue(ir), 1 );
		} else if ( counter > start_res && counter <= end_res ) {
			newpose.append_residue_by_bond( pose.residue(ir) ); //NOTE!  Assumes normal connectivity!
		}
	}

	std::stringstream new_fold_tree;
	new_fold_tree << "FOLD_TREE EDGE " << residues_per_repeat / 2 << " 1 -1 EDGE " << residues_per_repeat / 2 << " " << residues_per_repeat << " -1";

	core::kinematics::FoldTree foldtree;
	new_fold_tree >> foldtree;
	newpose.fold_tree(foldtree);

	pose=newpose;
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
SymmetricCycpepAlignCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SymmetricCycpepAlign );
}

std::string
SymmetricCycpepAlignCreator::keyname() const
{
	return SymmetricCycpepAlign::mover_name();
}

void SymmetricCycpepAlignCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymmetricCycpepAlign::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, SymmetricCycpepAlign const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //cyclic_peptide
