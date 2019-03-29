// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPalette.cc
/// @brief  PackerPalette: a class for storing the set of ResidueTypes
/// and VariantTypes that the packer uses by default, in the absence of any
/// TaskOperations that limit the set actually used.
/// @details The PackerPalette says, "Here are the types that you're
/// allowed to use, and which are on in the absence of TaskOperations."
/// TaskOperations then prune this, turning OFF types that have been
/// enabled.  This allows users to turn on noncanonicals for design, and
/// then use TaskOperations with the same commutativity rules (turning OFF
/// types only) that are used for canonicals, making mixed design with
/// canonicals and noncanonicals much easier.\nThis was implemented as
/// part of the 2016 Chemical XRW (eXtreme Rosetta Workshop).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

// Unit Headers
#include <core/pack/palette/PackerPalette.hh>

//Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
//Serialization functions for particular classes:
#include <core/chemical/ResidueTypeSet.srlz.hh>
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.PackerPalette" );

/// @brief Default constructor
PackerPalette::PackerPalette() :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< PackerPalette >(),
	restore_pre_talaris_behaviour_(false),
	icoor_05_2009_(false),
	pH_mode_(false),
	base_residue_type_names_(),
	residue_type_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ),
	variant_type_names_(),
	special_behaviours_(),
	terminal_types_( core::chemical::get_terminal_varianttypes() ),
	non_terminal_types_()
	//TODO -- initialize all private member vars here.
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	restore_pre_talaris_behaviour_ = option[ mistakes::restore_pre_talaris_2013_behavior ]();
	icoor_05_2009_ = option[corrections::chemical::icoor_05_2009]();
	pH_mode_ = option[pH::pH_mode]();
	initialize_special_behaviours(); //Sets this to be a map of each behaviour to FALSE.
	initialize_non_terminal_types();
}

/// @brief Constructor with ResidueTypeSet specifier.
/// @details Needed for packing with anything but the default (fa_standard) ResidueTypeSet.
PackerPalette::PackerPalette(
	core::chemical::ResidueTypeSetCOP residue_type_set_in
) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< PackerPalette >(),
	restore_pre_talaris_behaviour_(false),
	icoor_05_2009_(false),
	pH_mode_(false),
	base_residue_type_names_(),
	residue_type_set_( residue_type_set_in ),
	variant_type_names_(),
	special_behaviours_(),
	terminal_types_( core::chemical::get_terminal_varianttypes() ),
	non_terminal_types_()
	//TODO -- initialize all private member vars here.
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	restore_pre_talaris_behaviour_ = option[ mistakes::restore_pre_talaris_2013_behavior ]();
	icoor_05_2009_ = option[corrections::chemical::icoor_05_2009]();
	pH_mode_ = option[pH::pH_mode]();
	initialize_special_behaviours(); //Sets this to be a map of each behaviour to FALSE.
	initialize_non_terminal_types();
}


/// @brief Copy constructor
PackerPalette::PackerPalette(
	PackerPalette const &src
) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< PackerPalette >(),
	restore_pre_talaris_behaviour_(src.restore_pre_talaris_behaviour_),
	icoor_05_2009_(src.icoor_05_2009_),
	pH_mode_(src.pH_mode_),
	base_residue_type_names_( src.base_residue_type_names_ ),
	residue_type_set_( src.residue_type_set_ ),
	variant_type_names_( src.variant_type_names_ ),
	special_behaviours_( src.special_behaviours_ ),
	terminal_types_( src.terminal_types_ ),
	non_terminal_types_( src.non_terminal_types_ )
	//TODO -- copy all private member vars here.
{
	debug_assert( special_behaviours_.size() == (static_cast<core::Size>(END_OF_BEHAVIOUR_LIST) - 1) );
}

/// @brief Destructor
PackerPalette::~PackerPalette() {}

/// @brief The initialize_residue_level_task ("apply") function -- called during rotamer setup for the packer to get the
/// list of all ResidueTypes that are allowed.
/// @details Derived classes don't get to implement this to override the default behaviour (which is just to set up all
/// allowed ResidueTypes that result from the combination of the base types listed and the VariantTypes listed), but
/// specific variations on the default behaviour, provided by enum, are permitted.
/// Note that the PackerPalette is not meant to be used for position-specific setup, despite having access to the pose
/// residue.  Use TaskOperations for that.
/// @param [in] existing_residue The existing residue, for reference (though this should be largely unneeded).
/// It's largely only used for variant type matching.
/// @param [out] residue_type_list A std::list of ResidueTypeCOPs, which is cleared and populated by this function.
void
PackerPalette::initialize_residue_level_task(
	core::conformation::Residue const & existing_residue,
	std::list< chemical::ResidueTypeCOP > & residue_type_list
) const {
	residue_type_list.clear();

	// Bool to determine whether decide_what_to_do_with_base_type() function ends up checking the existing type or not.
	bool existing_type_processed( false );
	debug_assert( special_behaviours_.count( FORCE_EXISTING_BASE_TYPE ) );
	if ( !special_behaviours_.at( FORCE_EXISTING_BASE_TYPE ) ) {
		for ( core::Size i = 1, imax = base_residue_type_names_.size(); i <= imax; ++i ) {
			decide_what_to_do_with_base_type(
				existing_residue, residue_type_list, base_residue_type_names_[ i ].second, existing_type_processed );
		}
		if ( TR.Debug.visible() && residue_type_list.empty() ) {
			TR.Debug << "Empty list of residues for designing ";
			TR.Debug << existing_residue.name3() << existing_residue.seqpos() << std::endl;
		}
	}
	if ( ! existing_type_processed ) {
		//If we didn't already do the existing type, do it now.
		TR.Trace << existing_residue.name3() << existing_residue.seqpos() << "'s type not yet processed." << std::endl;
		decide_what_to_do_with_existing_type( existing_residue, residue_type_list );
	}
}


/// @brief Function to parse XML tags, implemented by derived classes.
/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
/// program execution.
void
PackerPalette::parse_my_tag(
	utility::tag::TagCOP const &tag,
	basic::datacache::DataMap const &/*datamap*/
) {
	if ( TR.Warning.visible() ) {
		TR.Warning << "Warning!  The RosettaScripts parse_my_tag() function was called, but this PackerPalette ";
		if ( tag.get() != nullptr ) TR.Warning << "( of type " << tag->getName() << ") ";
		TR.Warning << "does not implement this function!" << std::endl;
		TR.Warning << "Program execution continues, but the PackerPalette may not have been properly set up." << std::endl;
	}
}

void
PackerPalette::set_residue_type_set (
	core::chemical::ResidueTypeSetCOP new_type_set
) {
	residue_type_set_ = new_type_set;
}

/// @brief Add a base residue type name to the list of base residue type names.
void
PackerPalette::add_base_residue_type( std::string const &name ) {
	for ( core::Size i=1, imax=base_residue_type_names_.size(); i<=imax; ++i ) {
		runtime_assert_string_msg( base_residue_type_names_[i].first != name,
			"Error in core::pack::palette::PackerPalette::add_base_residue_type(): The type " +
			name + " has already been added." );
	}

	core::chemical::ResidueTypeCOP type_to_add( core::chemical::ResidueTypeFinder( *residue_type_set_ ).residue_base_name( name ).get_representative_type() );

	runtime_assert_string_msg( type_to_add != nullptr, "Error in core::pack::palette::PackerPalette::add_base_residue_type(): No base type corresponding to \"" + name + "\" could be found." );

	base_residue_type_names_.push_back( std::make_pair( name, type_to_add ) );
	if ( TR.Debug.visible() ) TR.Debug << "Adding base type " << name << " to PackerPalette." << std::endl;
}

/// @brief Clear the list of base residue types.
void
PackerPalette::clear_base_residue_type_list() {
	base_residue_type_names_.clear();
}

/// @brief   Set up the 'expanded' NCBB default base types
/// @author Andy Watkins (amw579@stanford.edu)
void
PackerPalette::set_up_expanded_base_types() {
	set_up_default_base_types();
	// Design aramids
	// Currently, it turns out, the only aramids added are the ones that should be defaults.
	// Really, we are doing 12 additions here, one for each backbone.
	add_base_residue_types_by_properties( utility::tools::make_vector1< chemical::ResidueProperty >( chemical::ARAMID ) );
	// Ditto BAAs
	add_base_residue_types_by_properties( utility::tools::make_vector1< chemical::ResidueProperty >( chemical::BETA_AA ) );
	// Oligoureas
	add_base_residue_types_by_properties( utility::tools::make_vector1< chemical::ResidueProperty >( chemical::OLIGOUREA ) );

	// A limited number of peptoids -- there isn't quite so clean a set of 'canonical' peptoids.
	/*
	* This is actually a living nightmare b/c we don't version these libraries and bad things happen to good people
	add_base_residue_type( "001" );
	//add_base_residue_type( "017" );
	add_base_residue_type( "3101" );
	//add_base_residue_type( "201" ); AMW: sarcosine is dead, long live sarcosine
	add_base_residue_type( "202" );
	add_base_residue_type( "203" );
	add_base_residue_type( "204" );
	add_base_residue_type( "205" );
	add_base_residue_type( "208" );
	add_base_residue_type( "211" );
	add_base_residue_type( "303" );
	add_base_residue_type( "305" );
	add_base_residue_type( "313" );
	add_base_residue_type( "314" );
	//add_base_residue_type( "317" );
	//add_base_residue_type( "318" );
	*/

	// AMW: once there are "canonical aramids" versus ridiculous "noncanonical aramids" (and so forth),
	// you should ensure that we only add the "canonical" ones here.

}

/// @brief   Set up the default base types:
/// @details Clears the list of base types first.
/// @author  Labonte <JWLabonte@jhu.edu>
void
PackerPalette::set_up_default_base_types() {
	using namespace core::chemical;

	clear_base_residue_type_list();

	if ( !residue_type_set_ ) return;

	// Default L-alpha amino acids
	// Temporarily, we are hard-coding this list. The reason
	// is fundamentally terrible: many unit and integration tests
	// are intimately dependent on the *order of base types
	// in the default configuration of the PackerTask*. Let's
	// defer doing things the right way until PackerPalettes are
	// clearly established -- this will also help us find legitimate
	// changes more easily.

	add_base_residue_type( "ALA" );
	add_base_residue_type( "CYS" );
	add_base_residue_type( "ASP" );
	add_base_residue_type( "GLU" );
	add_base_residue_type( "PHE" );
	add_base_residue_type( "GLY" );
	add_base_residue_type( "HIS" );
	if ( residue_type_set_->mode() != core::chemical::CENTROID_ROT_t ) add_base_residue_type( "HIS_D" );
	add_base_residue_type( "ILE" );
	add_base_residue_type( "LYS" );
	add_base_residue_type( "LEU" );
	add_base_residue_type( "MET" );
	add_base_residue_type( "ASN" );
	add_base_residue_type( "PRO" );
	add_base_residue_type( "GLN" );
	add_base_residue_type( "ARG" );
	add_base_residue_type( "SER" );
	add_base_residue_type( "THR" );
	add_base_residue_type( "VAL" );
	add_base_residue_type( "TRP" );
	add_base_residue_type( "TYR" );
	//add_base_residue_types_by_properties( { CANONICAL_AA } );

	// Protonation variants in pH mode.
	if ( pH_mode() ) {
		add_base_residue_type( "ASP_P1" );
		add_base_residue_type( "ASP_P2" );
		add_base_residue_type( "GLU_P1" );
		add_base_residue_type( "GLU_P2" );
		add_base_residue_type( "HIS_P" );
		add_base_residue_type( "LYS_D" );
		add_base_residue_type( "TYR_D" );
	}

	// Default nucleic acids (DNA)
	add_base_residue_types_by_properties( { CANONICAL_NUCLEIC, DNA } );

	// Default carbohydrates (all non-modified sugars)
	if ( basic::options::option[ basic::options::OptionKeys::in::include_sugars ] ) {
		add_base_residue_types_by_properties( { CARBOHYDRATE } );
	}
}

/// @brief  Add a group of base ResidueTypes and names to the PackerPalette by properties.
/// @author Labonte <JWLabonte@jhu.edu>
void
PackerPalette::add_base_residue_types_by_properties(
	utility::vector1< core::chemical::ResidueProperty > const & properties )
{
	utility::vector1< std::string > base_names;
	base_names.append( get_base_names_by_properties( { properties } ) );
	for ( auto const & base_name : base_names ) {
		add_base_residue_type( base_name );
	}
}


/// @brief Does the PackerPalette have a base residue type?
bool
PackerPalette::has_base_residue_type(
	std::string const &name
) {
	for ( auto const & entry : base_residue_type_names_ ) {
		if ( entry.first == name ) return true;
	}
	return false;
}

/// @brief Set the defaults for special behaviours.
void
PackerPalette::set_up_default_special_behaviours() {

	initialize_special_behaviours();
	special_behaviours_[ ONLY_DESIGN_POLYMER_RESIDUES ] = true;
	special_behaviours_[ ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ] = true;
	special_behaviours_[ KEEP_EXISTING_BASE_TYPE ] = true;
	special_behaviours_[ FORCE_EXISTING_BASE_TYPE ] = false;
	special_behaviours_[ ALL_DNA_TYPES_ON ] = true;
	special_behaviours_[ ONLY_ALPHA_AA_AT_ALPHA_POSITIONS ] = true;
	special_behaviours_[ ONLY_BETA_AA_AT_BETA_POSITIONS ] = true;
	special_behaviours_[ ONLY_GAMMA_AA_AT_GAMMA_POSITIONS ] = true;
	special_behaviours_[ ONLY_DNA_AT_DNA_POSITIONS ] = true;
	special_behaviours_[ ONLY_OLIGOUREA_AT_OLIGOUREA_POSITIONS ] = true;
	special_behaviours_[ ONLY_ARAMID_AT_ARAMID_POSITIONS ] = true;
	special_behaviours_[ ONLY_SACCHARIDES_AT_SACCHARIDE_POSITIONS ] = true;
	special_behaviours_[ ONLY_LIGAND_AT_LIGAND_POSITIONS ] = true;
	special_behaviours_[ ONLY_MATCHING_LIGAND_NAMES ] = true;
	special_behaviours_[ EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ] = true;
	special_behaviours_[ STRIP_VIRTUAL_SIDE_CHAIN ] = true;
	special_behaviours_[ pH_MODE_EXCEPTIONS ] = pH_mode();
	special_behaviours_[ ONLY_RNA_AT_RNA_POSITIONS ] = true;
	special_behaviours_[ KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS ] = true;
	special_behaviours_[ KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_FOR_EXISTING_BASE_TYPE ] = true;
	//special_behaviours_[ KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_FOR_EXISTING_BASE_TYPE ] = true;
	special_behaviours_[ KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES ] = true;
	special_behaviours_[ KEEP_EXISTING_DISULFIDES ] = true;
	special_behaviours_[ NO_METAPATCHES ] = true;
	special_behaviours_[ ALLOW_ALTERNATE_BACKBONE_MATCHING ] = false;

	if ( restore_pre_talaris_behaviour() || icoor_05_2009() ) {
		special_behaviours_[ ONLY_RNA_AT_RNA_POSITIONS ] = false;
	}

	debug_assert( special_behaviours_.size() == (static_cast<core::Size>(END_OF_BEHAVIOUR_LIST) - 1) );
}


// Protected Methods ///////////////////////////////////////////////////////////

/// @brief Set whether design is only limited to protein/peptoid/dna/saccharide positions, or can happen
/// everywhere.
/// @details The DefaultPackerPalette has this set to true; everything else should set it to false.
void
PackerPalette::set_only_design_protein_peptoid_dna_saccharide(
	bool const setting
) {
	debug_assert( special_behaviours_.count( ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ) );
	special_behaviours_[ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES] = setting;
}

/// @brief Set whether we're forcing the existing base type.
/// @details Defaults to false.  The NoDesignPackerPalette sets this to true.
void
PackerPalette::set_force_existing_base_type(
	bool const setting
) {
	debug_assert( special_behaviours_.count( FORCE_EXISTING_BASE_TYPE ) );
	special_behaviours_[ FORCE_EXISTING_BASE_TYPE ] = setting;
}

/// @brief Set whether design is only limited to polymer positions, or can happen
/// everywhere.
/// @details The DefaultPackerPalette has this set to true; everything else should set it to false.
void
PackerPalette::set_only_design_polymer_residues(
	bool const setting
) {
	debug_assert( special_behaviours_.count( ONLY_DESIGN_POLYMER_RESIDUES ) );
	special_behaviours_[ONLY_DESIGN_POLYMER_RESIDUES] = setting;
}

/// @brief   Return a list of base ResidueType names to be added to the PackerPalette by properties.
/// @author  Labonte <JWLabonte@jhu.edu>
utility::vector1< std::string >
PackerPalette::get_base_names_by_properties( utility::vector1< core::chemical::ResidueProperty > const & properties )
{
	using namespace core::chemical;

	ResidueTypeCOPs base_types(
		ResidueTypeFinder( *residue_type_set_ ).properties( properties ).get_possible_base_residue_types( true, true ) );
	utility::vector1< std::string > base_names;
	for ( auto const & base_type : base_types ) {
		base_names.push_back( base_type->name() );
	}

	return base_names;
}


// Private Methods /////////////////////////////////////////////////////////////

/// @brief Given a residue type pointer, get a raw pointer to the base type.
/// @details Avoids calls to get_self_ptr().
core::chemical::ResidueType const *
PackerPalette::get_base_type_raw_ptr(
	core::chemical::ResidueTypeCOP const & restype
) const {
	if ( restype->is_base_type() ) return restype.get();
	return restype->get_base_type_cop().get();
}

/// @brief Get a residue type with a variant removed.
core::chemical::ResidueTypeCOP
PackerPalette::get_residue_type_cop_with_variant_removed(
	core::chemical::ResidueTypeCOP const & restype,
	core::chemical::VariantType const vartype_to_remove
) const {
	if ( !restype->has_variant_type( vartype_to_remove ) ) return restype;

	runtime_assert( residue_type_set_ != nullptr );
	core::chemical::ResidueTypeFinder finder( *residue_type_set_ );
	finder.base_type( restype->get_base_type_cop() ).variants( restype->variant_type_enums(), restype->custom_variant_types(), false ).disallow_variants( utility::vector1< core::chemical::VariantType >( { vartype_to_remove } ), false );
	return finder.get_representative_type();
}

/// @brief Set the special behaviours to be a map of each behaviour type to FALSE.
/// @details All behaviours default to false, and must be explicitly set to true.
void
PackerPalette::initialize_special_behaviours() {
	special_behaviours_.clear();
	for ( core::Size i=1; i < static_cast<core::Size>(END_OF_BEHAVIOUR_LIST); ++i ) {
		special_behaviours_[ static_cast<SpecialPackerPaletteBehaviour>(i) ] = false;
	}
}

/// @brief Given the existing residue, decide what to do with its base type.
/// @details Depending on what's defined in the special_behaviours_ list, we might keep the base type,
/// keep the base type with terminal variants, keep the base type with terminal and side-chain variants,
/// etc.
void
PackerPalette::decide_what_to_do_with_existing_type(
	core::conformation::Residue const & existing_residue,
	std::list< core::chemical::ResidueTypeCOP> & residue_type_list ) const
{
	debug_assert( special_behaviours_.size() == (static_cast<core::Size>(END_OF_BEHAVIOUR_LIST) - 1) );

	if ( !residue_type_set_ ) { //If the ResidueTypeSetCOP is nullptr, keep the existing type only.
		utility::vector1 < core::chemical::ResidueTypeCOP > types_to_add;
		types_to_add.push_back( existing_residue.type_ptr() );
		add_residue_types_to_list( types_to_add, residue_type_list, existing_residue.seqpos() );
		return;
	}

	core::chemical::ResidueTypeCOP existing_type_pointer(
		( special_behaviours_.at( STRIP_VIRTUAL_SIDE_CHAIN ) && existing_residue.type().has_variant_type( core::chemical::VIRTUAL_SIDE_CHAIN ) ) ?
		get_residue_type_cop_with_variant_removed( existing_residue.type_ptr(), core::chemical::VIRTUAL_SIDE_CHAIN ) :
		existing_residue.type_ptr()
	);

	debug_assert( special_behaviours_.count( FORCE_EXISTING_BASE_TYPE ) );
	if ( special_behaviours_.at( FORCE_EXISTING_BASE_TYPE ) ) { //If we're definitely keeping the existing base type, do so here.
		residue_type_list.clear();
		residue_type_list.push_back( existing_type_pointer );
		return;
	}

	//Ensure that we have settings for relevant behaviours:
	debug_assert( special_behaviours_.count( ONLY_DESIGN_POLYMER_RESIDUES ) );
	debug_assert( special_behaviours_.count( ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ) );
	debug_assert( special_behaviours_.count( ALL_DNA_TYPES_ON ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_BASE_TYPE ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_FOR_EXISTING_BASE_TYPE ) );
	debug_assert( special_behaviours_.count( EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ) );
	debug_assert( special_behaviours_.count( STRIP_VIRTUAL_SIDE_CHAIN) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_DISULFIDES ) );
	debug_assert( special_behaviours_.count( pH_MODE_EXCEPTIONS ) );
	debug_assert( special_behaviours_.count( NO_METAPATCHES ) );
	debug_assert( special_behaviours_.count( ALLOW_ALTERNATE_BACKBONE_MATCHING ) );

	if ( ! special_behaviours_.at( KEEP_EXISTING_BASE_TYPE ) ) return;

	// KEEP THIS FIRST.  If we're not designing non-polymer residues, and if this is a non-polymer residue,
	// add the current type at this position and exit.
	if ( ( special_behaviours_.at( ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ) &&
			!existing_residue.type().is_protein() && !existing_residue.type().is_peptoid() &&
			!existing_residue.type().is_DNA() )
			||
			( special_behaviours_.at( ONLY_DESIGN_POLYMER_RESIDUES ) && !existing_residue.type().is_polymer() )
			) {
		utility::vector1 < core::chemical::ResidueTypeCOP > types_to_add;
		types_to_add.push_back( existing_type_pointer );
		add_residue_types_to_list( types_to_add, residue_type_list, existing_residue.seqpos() );
		return;
	}

	// KEEP THIS SECOND.  Do nothing if we're not keeping the existing base type:
	if ( !special_behaviours_.at( KEEP_EXISTING_BASE_TYPE ) ) return;

	// If we're not designing disulfides, and this is a disulfide, add the unmodified residue type:
	if ( special_behaviours_.at( KEEP_EXISTING_DISULFIDES ) && existing_residue.type().is_disulfide_bonded() ) {
		utility::vector1 < core::chemical::ResidueTypeCOP > types_to_add;
		types_to_add.push_back( existing_type_pointer );
		add_residue_types_to_list( types_to_add, residue_type_list, existing_residue.seqpos() );
		return;
	}

	if ( special_behaviours_.at( ALL_DNA_TYPES_ON ) && existing_residue.type().is_DNA() ) {
		core::chemical::ResidueTypeFinder finder( *residue_type_set_ ); //Inefficient, but hopefully we can deprecate this behaviour soon.
		finder.variants( existing_residue.type().variant_type_enums() ).base_property( core::chemical::DNA );
		if ( special_behaviours_.at( EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ) ) {
			finder.variant_exceptions( utility::tools::make_vector1( core::chemical::ADDUCT_VARIANT ), false );
		}
		utility::vector1< core::chemical::ResidueTypeCOP > dna_types_to_add( finder.get_all_possible_residue_types( false ) );
		add_residue_types_to_list( dna_types_to_add, residue_type_list, existing_residue.seqpos() );
		return;
	}

	//If we reach here, we're ready to start a residue type search.
	utility::vector1< core::chemical::VariantType > variants;
	utility::vector1< std::string > special_variants;
	utility::vector1< core::chemical::VariantType > exception_variants;

	//Is this pH mode?
	if ( pH_mode() && special_behaviours_.at( pH_MODE_EXCEPTIONS ) ) {
		utility::vector1< core::chemical::VariantType > pH_exceptions( core::chemical::pH_mode_exceptions() );
		for ( core::Size i(1), imax(pH_exceptions.size()); i<=imax; ++i ) {
			exception_variants.push_back( pH_exceptions[i] );
		}
	}

	//Keeping the terminal variant types:
	if ( special_behaviours_.at( KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS ) ) {
		utility::vector1 < core::chemical::VariantType > present_terminal_varianttypes;
		utility::vector1 < std::string > present_on_the_fly_types;
		get_types_on_residue( existing_residue, terminal_types_, present_terminal_varianttypes, present_on_the_fly_types );
		for ( core::Size i(1), imax( present_terminal_varianttypes.size() ); i<=imax; ++i ) {
			variants.push_back( present_terminal_varianttypes[i] );
		}
		//Currently, we assume that on-the-fly types are not terminal types.
		/*for( core::Size i(1), imax( present_on_the_fly_types.size() ); i<=imax; ++i ) {
		special_variants.push_back( present_on_the_fly_types[i] );
		}*/
	}

	if ( special_behaviours_.at(KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES) ) {
		utility::vector1 < core::chemical::VariantType > present_nonterminal_varianttypes;
		utility::vector1 < std::string > present_on_the_fly_types;
		get_types_on_residue( existing_residue, non_terminal_types_, present_nonterminal_varianttypes, present_on_the_fly_types );
		for ( auto const & present_nonterminal_variant : present_nonterminal_varianttypes ) {
			variants.push_back(present_nonterminal_variant);
		}
		for ( auto const & on_the_fly_type : present_on_the_fly_types ) {
			special_variants.push_back( on_the_fly_type );
		}
	}

	if ( special_behaviours_.at( EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ) && existing_residue.type().is_DNA() ) {
		//Need this again in case ALL_DNA_TYPES_ON is false.
		exception_variants.push_back( core::chemical::ADDUCT_VARIANT );
	}

	if ( special_behaviours_.at( STRIP_VIRTUAL_SIDE_CHAIN ) && existing_residue.type().has_variant_type( core::chemical::VIRTUAL_SIDE_CHAIN ) ) {
		exception_variants.push_back( core::chemical::VIRTUAL_SIDE_CHAIN );
	}

	// TODO at some point: Could add an alternative behaviour that adds all combinations of present/absent existing
	// nonterminal variant types.  For example, if the existing type is N-methyl-phosphotyrosone, we'd add the base
	// type, the N-methyl type, the phospho type, and the N-methyl-phosho type.
	// Yes, please!  ~ Labonte

	utility::vector1< core::chemical::ResidueTypeCOP > const types_to_add(
		residue_type_set_->get_all_types_with_variants_by_basetype(
		existing_residue.type().is_base_type() ? existing_residue.type_ptr() : existing_residue.type().get_base_type_cop(),
		variants,
		special_variants,
		exception_variants,
		special_behaviours_.at( NO_METAPATCHES )
		)
	);
	add_residue_types_to_list( types_to_add, residue_type_list, existing_residue.seqpos() );

	if ( residue_type_list.empty() ) {
		residue_type_list.push_back( existing_type_pointer );
	}
}

/// @brief   Given the existing residue and a candidate base type in the base types list, decide what to
///          do with the candidate base type.
/// @details If the candidate type is the existing type, then this calls
///          decide_what_to_do_with_existing_type() and marks existing_type_processed
///          as true; otherwise, existing_type_processed is left unchanged.
/// @note    If we ever decide to allow design of one class/family of backbone to another, we will need to have this
///          function overwritten by the child class.
void
PackerPalette::decide_what_to_do_with_base_type(
	core::conformation::Residue const & existing_residue,
	std::list< core::chemical::ResidueTypeCOP> &residue_type_list,
	core::chemical::ResidueTypeCOP const & candidate_base_type,
	bool &existing_type_processed
) const {

	debug_assert( special_behaviours_.size() == (static_cast<core::Size>(END_OF_BEHAVIOUR_LIST) - 1) );

	// If the ResidueTypeSetCOP is null, then we can't design:
	if ( !residue_type_set_ ) {
		return;
	}

	// Confirm that there are settings (true or false) for relevant options.
	debug_assert( special_behaviours_.count( ONLY_DESIGN_POLYMER_RESIDUES ) );
	debug_assert( special_behaviours_.count( ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ) );
	debug_assert( special_behaviours_.count( ONLY_ALPHA_AA_AT_ALPHA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_BETA_AA_AT_BETA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_GAMMA_AA_AT_GAMMA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_DNA_AT_DNA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_OLIGOUREA_AT_OLIGOUREA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_ARAMID_AT_ARAMID_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_SACCHARIDES_AT_SACCHARIDE_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_LIGAND_AT_LIGAND_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_MATCHING_LIGAND_NAMES ) );
	debug_assert( special_behaviours_.count( EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ) );
	debug_assert( special_behaviours_.count( ONLY_RNA_AT_RNA_POSITIONS ) );
	debug_assert( special_behaviours_.count( pH_MODE_EXCEPTIONS ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES ) );
	debug_assert( special_behaviours_.count( KEEP_EXISTING_DISULFIDES ) );
	debug_assert( special_behaviours_.count( NO_METAPATCHES ) );
	debug_assert( special_behaviours_.count( ALLOW_ALTERNATE_BACKBONE_MATCHING ) );

	// If this is the existing type, then decide what to do with the existing type.  KEEP THIS FIRST.
	if ( get_base_type_raw_ptr( candidate_base_type ) == get_base_type_raw_ptr( existing_residue.type_ptr() ) ) {
		decide_what_to_do_with_existing_type( existing_residue, residue_type_list );
		existing_type_processed = true;
		return;
	}

	//BACKBONE TYPE CHECKS:
	//Non-polymer:
	if ( special_behaviours_.at( ONLY_DESIGN_POLYMER_RESIDUES ) && !existing_residue.type().is_polymer() ) {
		return; //Do nothing if this isn't a polymer residue and we're not designing non-polymer residues.
	}

	if ( special_behaviours_.at( ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES ) &&
			! existing_residue.type().is_protein() &&
			! existing_residue.type().is_peptoid() &&
			! existing_residue.type().is_aramid() &&
			! existing_residue.type().is_DNA() &&
			! existing_residue.type().is_carbohydrate() ) {
		return; //Do nothing if this isn't a protein/peptoid/DNA/sugar residue and if we're not designing non-protein/non-peptoid residues.
	}

	//Alpha amino acids and/or peptoids:
	// AMW: now with the NCAA default packer palette, only-alpha-at-alpha-positions should not
	// permit alpha-to-peptoid by default. That ought to require a separate TBD behavior.
	if (
			!( restore_pre_talaris_behaviour() || icoor_05_2009() ) &&
			special_behaviours_.at( ONLY_ALPHA_AA_AT_ALPHA_POSITIONS ) &&
			( ( existing_residue.type().is_alpha_aa() || existing_residue.type().is_peptoid() ) &&
			!( candidate_base_type->is_alpha_aa() || candidate_base_type->is_peptoid() ) )
			) {
		return;
	} else if ( //Special case -- old pre-talaris behaviour (which was also before the is_alpha_aa() property).
			(restore_pre_talaris_behaviour() || icoor_05_2009() ) &&
			special_behaviours_.at( ONLY_ALPHA_AA_AT_ALPHA_POSITIONS ) &&
			( existing_residue.type().is_protein() || existing_residue.type().is_peptoid() ) &&
			!( candidate_base_type->is_protein() || candidate_base_type->is_peptoid() )
			) {
		return;
	}

	//Beta amino acids:
	if ( special_behaviours_.at( ONLY_BETA_AA_AT_BETA_POSITIONS ) &&
			( existing_residue.type().is_beta_aa() ) &&
			!( candidate_base_type->is_beta_aa() )
			) {
		return;
	}
	//Gamma amino acids:
	if ( special_behaviours_.at( ONLY_GAMMA_AA_AT_GAMMA_POSITIONS ) &&
			( existing_residue.type().is_gamma_aa() ) &&
			!( candidate_base_type->is_gamma_aa() )
			) {
		return;
	}

	//DNA:
	if ( special_behaviours_.at( ONLY_DNA_AT_DNA_POSITIONS ) &&
			( existing_residue.type().is_DNA() ) &&
			!( candidate_base_type->is_DNA() )
			) {
		return;
	}

	//Oligourea:
	if ( special_behaviours_.at( ONLY_OLIGOUREA_AT_OLIGOUREA_POSITIONS ) &&
			existing_residue.type().is_oligourea() &&
			!( candidate_base_type->is_oligourea() ) ) {
		return;
	}

	//Aramids:
	if ( existing_residue.type().is_aramid() && special_behaviours_.at( ONLY_ARAMID_AT_ARAMID_POSITIONS ) ) {
		if ( !aramid_backbones_are_compatible( existing_residue.type(), *candidate_base_type ) ) {
			return;
		}
	}

	// Carbohydrates:
	if ( existing_residue.type().is_carbohydrate() ) {
		if ( special_behaviours_.at( ONLY_SACCHARIDES_AT_SACCHARIDE_POSITIONS ) &&
				! candidate_base_type->is_carbohydrate()  ) {
			return;
		}

		// to design a saccharide residue, we must first check that the same type of backbone exists between the
		// existing residue and any candidate residues.
		if ( ! saccharide_backbones_are_compatible( existing_residue.type(), *candidate_base_type ) ) {
			return;
		}
	}

	//Ligands:
	if ( special_behaviours_.at( ONLY_LIGAND_AT_LIGAND_POSITIONS ) ) {
		if ( existing_residue.type().is_ligand() && !( candidate_base_type->is_ligand() ) ) {
			return;
		}
		if ( !( existing_residue.type().is_ligand() ) && candidate_base_type->is_ligand() ) {
			return;
		}
		if ( special_behaviours_.at( ONLY_MATCHING_LIGAND_NAMES ) && existing_residue.type().is_ligand() && candidate_base_type->is_ligand() && existing_residue.type().base_name() != candidate_base_type->base_name() ) {
			return;
		}
	}

	//RNA:
	if ( special_behaviours_.at( ONLY_RNA_AT_RNA_POSITIONS ) &&
			( existing_residue.type().is_RNA() ) &&
			!( candidate_base_type->is_RNA() )
			) {
		return;
	}
	//Disulfides:
	if ( special_behaviours_.at( KEEP_EXISTING_DISULFIDES ) && existing_residue.type().is_disulfide_bonded() ) {
		return;
	}

	//If we reach here, we're ready to start a residue type search.
	utility::vector1< core::chemical::VariantType > variants;
	utility::vector1< std::string > special_variants;
	utility::vector1< core::chemical::VariantType > exception_variants;

	//Is this pH mode?
	if ( pH_mode() && special_behaviours_.at( pH_MODE_EXCEPTIONS ) ) {
		utility::vector1< core::chemical::VariantType > pH_exceptions( core::chemical::pH_mode_exceptions() );
		for ( core::Size i(1), imax(pH_exceptions.size()); i<=imax; ++i ) {
			exception_variants.push_back( pH_exceptions[i] );
		}
	}

	//Keeping the terminal variant types:
	if ( special_behaviours_.at( KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS ) ) {
		utility::vector1 < core::chemical::VariantType > present_terminal_varianttypes;
		utility::vector1 < std::string > present_on_the_fly_types;
		get_types_on_residue( existing_residue, terminal_types_, present_terminal_varianttypes, present_on_the_fly_types );
		for ( core::Size i(1), imax( present_terminal_varianttypes.size() ); i<=imax; ++i ) {
			variants.push_back( present_terminal_varianttypes[i] );
		}
		//Currently, we assume that on-the-fly types are not terminal types.
		/*for( core::Size i(1), imax( present_on_the_fly_types.size() ); i<=imax; ++i ) {
		special_variants.push_back( present_on_the_fly_types[i] );
		}*/
	}

	if ( special_behaviours_.at(KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES) ) {
		utility::vector1 < core::chemical::VariantType > present_nonterminal_varianttypes;
		utility::vector1 < std::string > present_on_the_fly_types;
		get_types_on_residue( existing_residue, non_terminal_types_, present_nonterminal_varianttypes, present_on_the_fly_types );
		for ( auto const & present_nonterminal_variant : present_nonterminal_varianttypes ) {
			variants.push_back(present_nonterminal_variant);
		}
		for ( auto const & on_the_fly_type : present_on_the_fly_types ) {
			special_variants.push_back( on_the_fly_type );
		}
	}

	if ( special_behaviours_.at( EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS ) && existing_residue.type().is_DNA() ) {
		exception_variants.push_back( core::chemical::ADDUCT_VARIANT );
	}

	if ( special_behaviours_.at( STRIP_VIRTUAL_SIDE_CHAIN ) && existing_residue.type().has_variant_type( core::chemical::VIRTUAL_SIDE_CHAIN ) ) {
		exception_variants.push_back( core::chemical::VIRTUAL_SIDE_CHAIN );
	}

	// TODO at some point: Could add an alternative behaviour that adds all combinations of present/absent existing
	// nonterminal variant types.  For example, if the existing type is N-methyl-phosphotyrosone, we'd add the base
	// type, the N-methyl type, the phospho type, and the N-methyl-phosho type.
	// Yes, please!  ~ Labonte

	utility::vector1< core::chemical::ResidueTypeCOP > const types_to_add(
		residue_type_set_->get_all_types_with_variants_by_basetype(
		candidate_base_type,
		variants,
		special_variants,
		exception_variants,
		special_behaviours_.at( NO_METAPATCHES )
		)
	);
	add_residue_types_to_list( types_to_add, residue_type_list, existing_residue.seqpos() );
}

/// @brief Set up the list of variant types that don't modify termini.
/// @details This is just all of the VariantTypes that exist that are NOT in the terminal_types_ list.
void
PackerPalette::initialize_non_terminal_types() {
	non_terminal_types_.clear();
	core::Size const n_terminal_types( terminal_types_.size() );
	non_terminal_types_.reserve( static_cast<core::Size>( core::chemical::N_VARIANTS ) - n_terminal_types );

	for ( core::Size i=1; i<=static_cast<core::Size>( core::chemical::N_VARIANTS ); ++i ) {

		bool add_to_list( true );
		for ( core::Size j=1; j<=n_terminal_types; ++j ) {
			if ( static_cast<core::Size>( terminal_types_[j] ) == i ) {
				add_to_list=false;
				break;
			}
		}

		if ( add_to_list ) {
			non_terminal_types_.push_back( static_cast< core::chemical::VariantType >( i ) );
		}
	}
	debug_assert( non_terminal_types_.size() == ( static_cast<core::Size>( core::chemical::N_VARIANTS ) - n_terminal_types ) );
}

/// @brief Given a list of VariantTypes and a residue, populate two new lists with the VariantTypes present on the
/// residue and the custom VariantTypes present.
/// @details present_types and present_on_the_fly_types are cleared and populated by this operation.  The
/// present_on_the_fly_types vector is a vector of strings representing types generated on the fly, for which enum don't
/// exist.  (The enzdes and metalloprotein code generate custom ResidueTypes in memory, for example.)
void
PackerPalette::get_types_on_residue(
	core::conformation::Residue const & residue,
	utility::vector1 < core::chemical::VariantType > const & types,
	utility::vector1 < core::chemical::VariantType > & present_types,
	utility::vector1 < std::string > &present_on_the_fly_types
) const {
	present_types.clear();
	present_on_the_fly_types.clear();
	for ( core::Size i=1, imax=types.size(); i<=imax; ++i ) { //Loop through everything in the types list.
		if ( residue.type().has_variant_type( types[i] ) ) {
			present_types.push_back( types[i] );
		}
	}
	present_on_the_fly_types = residue.type().custom_variant_types();
}

/// @brief Given some ResidueTypes to add to the ResidueType list, add the ones that are not already in the list.
/// @details Determines equivalency by pointer address comparison.
void
PackerPalette::add_residue_types_to_list(
	utility::vector1< core::chemical::ResidueTypeCOP > const &types_to_add,
	std::list< core::chemical::ResidueTypeCOP > &residue_type_list,
	core::Size const seqpos
) const {
	for ( core::Size i=1, imax=types_to_add.size(); i<=imax; ++i ) {
		if ( special_behaviours_.at( STRIP_VIRTUAL_SIDE_CHAIN ) && types_to_add[i]->has_variant_type( core::chemical::VIRTUAL_SIDE_CHAIN ) ) {
			continue;
		}
		bool skip_to_next(false);
		for ( std::list< core::chemical::ResidueTypeCOP >::iterator itr=residue_type_list.begin();
				itr != residue_type_list.end();
				++itr
				) { //Check that we haven't already added this type.  We'll just use pointer comparisons for speed and simplicity.
			if ( types_to_add[i].get() == (*itr).get() ) {
				skip_to_next=true;
				break;
			}
		}
		if ( skip_to_next ) continue; //If the pointers match, skip to the next

		if ( TR.Debug.visible() ) {
			TR.Debug << "Adding " << types_to_add[ i ]->name() << " at position " << seqpos << "." << std::endl;
		}
		residue_type_list.push_back( types_to_add[ i ] );
	}
}


/// @brief    Do these two ResidueTypes have compatible backbones?
/// @details  Two residues must have superimposable main chains for one to be a candidate for design for the other.
/// @author   Labonte <JWLabonte@jhu.edu>
bool
PackerPalette::saccharide_backbones_are_compatible(
	chemical::ResidueType const & existing_type,
	chemical::ResidueType const & candidate_type ) const
{
	if ( ! special_behaviours_.at( ALLOW_ALTERNATE_BACKBONE_MATCHING ) ) {
		if ( existing_type.mainchain_atoms().size() != candidate_type.mainchain_atoms().size() ) {
			return false;
		}
	}

	core::chemical::carbohydrates::CarbohydrateInfoCOP existing_info( existing_type.carbohydrate_info() ),
		candidate_info( candidate_type.carbohydrate_info() );

	if ( existing_info->anomeric_carbon() != candidate_info->anomeric_carbon() ) {
		// TODO: I think that I could actually design an aldose into a ketose provided that the ring size and length
		// of the main chain are the same, because Residue uses atom indices to place new coordinates, but I'm going
		// to start with the easy cases first.  ~ Labonte
		// Further, such a case should probably be considered a case of an alternate backbone anyhow.
		return false;
	}
	if ( existing_info->ring_size() != candidate_info->ring_size() ) {
		return false;
	}
	if ( existing_info->mainchain_glycosidic_bond_acceptor() !=
			candidate_info->mainchain_glycosidic_bond_acceptor() ) {
		if ( ! existing_type.is_upper_terminus() ) {  // We can be less stringent at the non-reducing end.
			return false;
		}
	}
	if ( ( existing_info->anomer() == candidate_info->anomer() &&
			existing_info->stereochem() == candidate_info->stereochem() ) ||
			( existing_info->anomer() != candidate_info->anomer() &&
			existing_info->stereochem() != candidate_info->stereochem() ) ) {
		// (Alpha for a D-sugar is axial, but beta for an L-sugar is also axial, so these also match.)
	} else {
		return false;
	}
	// Compare the stereochemistry of the On atoms, where n is the linkage site, by comparing internal coordinates.
	// Note: This check only works because the carbohydrate .params files have been consistently generated with
	// idealized torsion angles.
	if ( existing_type.icoor( existing_info->mainchain_glycosidic_bond_acceptor() + 1 ).phi() !=
			candidate_type.icoor( candidate_info->mainchain_glycosidic_bond_acceptor() + 1 ).phi() ) {
		if ( ! existing_type.is_upper_terminus() ) {  // We can be less stringent at the non-reducing end.
			return false;
		}
	}
	return true;
}

/// @brief Do these two ResidueTypes have compatible aramid backbones?
/// @note Returns false if either residue is not an aramid, or if they are both aramids
/// but of incompatible types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
bool
PackerPalette::aramid_backbones_are_compatible(
	core::chemical::ResidueType const &existing_type,
	core::chemical::ResidueType const &candidate_type
) const {
	using namespace core::chemical;

	if ( !existing_type.is_aramid() || !candidate_type.is_aramid() ) return false;

	static const utility::fixedsizearray1< ResidueProperty, 12 > aramid_types ({ ORTHO_ARAMID, META_ARAMID, PARA_ARAMID,
		PRE_METHYLENE_ORTHO_ARAMID, PRE_METHYLENE_META_ARAMID, PRE_METHYLENE_PARA_ARAMID,
		POST_METHYLENE_ORTHO_ARAMID, POST_METHYLENE_META_ARAMID, POST_METHYLENE_PARA_ARAMID,
		PRE_METHYLENE_POST_METHYLENE_ORTHO_ARAMID, PRE_METHYLENE_POST_METHYLENE_META_ARAMID, PRE_METHYLENE_POST_METHYLENE_PARA_ARAMID
		});

	for ( ResidueProperty const & type : aramid_types ) {
		if ( existing_type.has_property( type ) != candidate_type.has_property( type ) ) {
			return false;
		}
	}

	return true;
}

} // palette
} // pack
} // core

#ifdef    SERIALIZATION

/// @brief "Save" function for serialization.
template< class Archive >
void
core::pack::palette::PackerPalette::save( Archive & arc ) const {
	arc( CEREAL_NVP( restore_pre_talaris_behaviour_ ) ); // _Bool
	arc( CEREAL_NVP( icoor_05_2009_ ) ); // _Bool
	arc( CEREAL_NVP( pH_mode_ ) ); // _Bool
	serialize_base_residue_type_names( arc, base_residue_type_names_); // utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >
	core::chemical::serialize_residue_type_set( arc, residue_type_set_ );
	arc( CEREAL_NVP( variant_type_names_ ) ); // utility::vector1 < std::string >
	serialize_behaviours_map( arc, special_behaviours_ ); // std::map < SpecialPackerPaletteBehaviour, bool >
	serialize_VariantType_vector( arc, terminal_types_ ); // utility::vector1 < core::chemical::VariantType >
	serialize_VariantType_vector( arc, non_terminal_types_ ); // utility::vector1 < core::chemical::VariantType >
}

/// @brief "Load" function for serialization.
template< class Archive >
void
core::pack::palette::PackerPalette::load( Archive & arc ) {
	arc( restore_pre_talaris_behaviour_ ); // _Bool
	arc( icoor_05_2009_ ); // _Bool
	arc( pH_mode_ ); // _Bool
	deserialize_base_residue_type_names( arc, base_residue_type_names_); // utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >
	core::chemical::deserialize_residue_type_set( arc, residue_type_set_ );
	arc( variant_type_names_ ); // utility::vector1 < std::string >
	deserialize_behaviours_map( arc, special_behaviours_ ); // std::map < SpecialPackerPaletteBehaviour, bool >
	deserialize_VariantType_vector( arc, terminal_types_ ); // utility::vector1 < core::chemical::VariantType >
	deserialize_VariantType_vector( arc, non_terminal_types_ ); // utility::vector1 < core::chemical::VariantType >
}

/// @brief Given a utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >,
/// serialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::serialize_base_residue_type_names(
	Archive &archive,
	utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > > const &vect
) const {
	archive( vect.size() );
	for ( core::Size i=1, imax=vect.size(); i<=imax; ++i ) {
		archive( vect[i].first );
		core::chemical::serialize_residue_type(archive, vect[i].second );
	}
}

/// @brief Given a utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >,
/// deserialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::deserialize_base_residue_type_names(
	Archive &archive,
	utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > > &vect
) const {
	core::Size vectsize;
	archive( vectsize );
	vect.resize(vectsize);
	for ( core::Size i=1; i<=vectsize; ++i ) {
		std::string resname;
		core::chemical::ResidueTypeCOP restype;
		archive(resname);
		core::chemical::deserialize_residue_type(archive, restype );
		vect[i] = std::pair< std::string, core::chemical::ResidueTypeCOP >( resname, restype );
	}
}


/// @brief Given a std::map of SpecialPackerPaletteBahaviour enum to bool,
/// serialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::serialize_behaviours_map(
	Archive &arc,
	std::map < SpecialPackerPaletteBehaviour, bool > const &mymap
) const {
	arc( mymap.size() );
	for ( std::map< SpecialPackerPaletteBehaviour, bool >::const_iterator it( mymap.begin() ); it!=mymap.end(); ++it ) {
		arc( static_cast<core::Size>( it->first ) );
		arc( it->second );
	}
}

/// @brief Given a std::map of SpecialPackerPaletteBahaviour enum to bool,
/// deserialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::deserialize_behaviours_map(
	Archive &arc,
	std::map < SpecialPackerPaletteBehaviour, bool > &mymap
) const {
	core::Size mapsize;
	arc( mapsize );
	mymap.clear();
	for ( core::Size i=1; i<=mapsize; ++i ) {
		core::Size key; bool value;
		arc(key);
		arc(value);
		mymap[static_cast<SpecialPackerPaletteBehaviour>(key)] = value;
	}
}


/// @brief Given a vector of VariantTypes, serialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::serialize_VariantType_vector(
	Archive &archive,
	utility::vector1 < core::chemical::VariantType > const &vect
) const {
	archive( vect.size() );
	for ( core::Size i=1, imax=vect.size(); i<=imax; ++i ) {
		archive( static_cast<core::Size>( vect[i] ) );
	}
}

/// @brief Given a vector of VariantTypes, deserialize it.
template < class Archive >
void
core::pack::palette::PackerPalette::deserialize_VariantType_vector(
	Archive &archive,
	utility::vector1 < core::chemical::VariantType > &vect
) const {
	core::Size vectsize;
	archive( vectsize );
	vect.resize( vectsize );
	for ( core::Size i=1; i<=vectsize; ++i ) {
		core::Size curentry;
		archive( curentry );
		vect[i] = static_cast<core::chemical::VariantType>(curentry);
	}
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::palette::PackerPalette );
CEREAL_REGISTER_TYPE( core::pack::palette::PackerPalette )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_palette_PackerPalette )
#endif // SERIALIZATION
