// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/PatchOperation.cc
/// @brief  Polymorphic classes representing the contents of a residue-type patch file
/// @author Phil Bradley

// Unit headers
#include <core/chemical/PatchOperation.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>

// Project headers
#include <core/chemical/Atom.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/icoor_support.hh>

// Numeric headers
#include <numeric/conversions.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <core/chemical/MutableICoorRecord.hh> // AUTO IWYU For MutableICoorRecord


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @details Auto-generated virtual destructor
PatchOperation::~PatchOperation() = default;

static basic::Tracer tr( "core.chemical" );
static basic::Tracer TR_PatchOperations( "core.chemical.PatchOperations.hh" );

DeleteAtom::DeleteAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
DeleteAtom::apply( MutableResidueType & rsd ) const
{

	if ( !rsd.has( atom_name_ )  ) {
		TR_PatchOperations.Debug << "DeleteAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "DeleteAtom::apply: deleting: " << atom_name_ << std::endl;
		}
		rsd.delete_atom( atom_name_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteAtom::name() const {
	return "DeleteAtom";
}

SetBackboneHeavyatom::SetBackboneHeavyatom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
SetBackboneHeavyatom::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetBackboneHeavyatom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ <<
			std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetBackboneHeavyatom::apply: " << atom_name_ << std::endl;
		}
		rsd.set_backbone_heavyatom( atom_name_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetBackboneHeavyatom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetBackboneHeavyatom::name() const {
	return "SetBackboneHeavyatom";
}


SetDisulfideAtomName::SetDisulfideAtomName( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
SetDisulfideAtomName::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetDisulfideAtomName::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ <<
			std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetDisulfideAtomName::apply: " << atom_name_ << std::endl;
		}
		rsd.set_disulfide_atom_name( atom_name_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetDisulfideAtomName").
/// @author Andy Watkins (amw579@stanford.edu).
std::string
SetDisulfideAtomName::name() const {
	return "SetDisulfideAtomName";
}

/// @brief constructor the type of connection is "LOWER" or "UPPER"
SetPolymerConnectAtom::SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in ) :
	atom_name_( atom_name_in )
{
	if ( upper_lower_in == "LOWER" ) {
		upper_lower_ = -1;
	} else if ( upper_lower_in == "UPPER" ) {
		upper_lower_ = 1;
	} else {
		utility_exit_with_message( "SetPolymerConnectAtom: unrecognized switch "+upper_lower_in );
	}
}

bool
SetPolymerConnectAtom::apply( MutableResidueType & rsd ) const
{
	if ( atom_name_ != "NONE" && ! rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetPolymerConnectAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	}

	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "SetPolymerConnectAtom::apply: " <<
			atom_name_ << ' ' << upper_lower_ << std::endl;
	}

	if ( upper_lower_ == -1 ) {
		rsd.set_lower_connect_atom( atom_name_ );
		if ( atom_name_ == "NONE" ) {
			clean_up_dangling_connect( rsd, ICoordAtomIDType::POLYMER_LOWER );
		}
	} else {
		debug_assert( upper_lower_ == 1 );
		rsd.set_upper_connect_atom( atom_name_ );
		if ( atom_name_ == "NONE" ) {
			clean_up_dangling_connect( rsd, ICoordAtomIDType::POLYMER_UPPER );
		}
	}

	return false;
}

bool
SetPolymerConnectAtom::changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const
{
	if ( atom_name_ == "NONE" ) {
		if ( upper_lower_ == -1 ) {
			return rsd_type.lower_connect_id() != 0 && ( atom == "LOWER" || (
					rsd_type.has( atom_name_ ) && rsd_type.atom_index( atom ) == rsd_type.lower_connect_atom() ) );
		} else {
			return rsd_type.upper_connect_id() != 0 && ( atom == "UPPER" || (
					rsd_type.has( atom_name_ ) && rsd_type.atom_index( atom ) == rsd_type.upper_connect_atom() ) );
		}
	} else if ( atom == "LOWER" ) {
		return upper_lower_ == -1;
	} else if ( atom == "UPPER" ) {
		return upper_lower_ != -1;
	} else {
		return rsd_type.has( atom ) && rsd_type.has( atom_name_ ) && rsd_type.atom_index( atom ) == rsd_type.atom_index( atom_name_ );
	}
}


/// @brief Return the name of this PatchOperation ("SetPolymerConnectAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetPolymerConnectAtom::name() const {
	return "SetPolymerConnectAtom";
}

AddConnect::AddConnect( std::string const & connect_atom,
	Real const phi, Real const theta, Real const d,
	std::string const & parent_atom,
	std::string const & angle_atom,
	std::string const & torsion_atom
):
	connect_atom_( connect_atom ),
	phi_( phi ), theta_( theta ), d_( d ),
	parent_atom_ (  parent_atom ),
	angle_atom_  (   angle_atom ),
	torsion_atom_( torsion_atom )
{}


bool
AddConnect::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( connect_atom_ ) ||
			!rsd.has(  parent_atom_ ) ||
			!rsd.has(   angle_atom_ ) ||
			!rsd.has( torsion_atom_ ) ) return true; // failure!

	Size const connid( rsd.add_residue_connection( connect_atom_ ) );
	rsd.set_icoor( "CONN"+ObjexxFCL::string_of( connid ), phi_, theta_, d_, parent_atom_, angle_atom_, torsion_atom_ );
	return false;
}

bool
AddConnect::changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const
{
	return (atom == connect_atom_ ) && rsd_type.has( atom ) && rsd_type.has( connect_atom_ );
}

/// @brief Return the name of this PatchOperation ("AddConnect").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddConnect::name() const {
	return "AddConnect";
}

AddProperty::AddProperty( std::string const & property_in ):
	property_( property_in ),
	property_enum_( ResidueProperties::get_property_from_string( property_in ) /*Will default to NO_PROPERTY if this is an on-the-fly property*/ )
{}

bool
AddProperty::apply( MutableResidueType & rsd ) const
{
	runtime_assert_string_msg( property_enum_ != core::chemical::NO_PROPERTY, "Error in AddProperty::apply(): Could not parse \"" + property_ + "\" as a valid residue type property." );
	rsd.add_property( property_enum_ );
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "AddProperty::apply: " << property_ << std::endl;
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("AddProperty").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddProperty::name() const {
	return "AddProperty";
}

DeleteProperty::DeleteProperty( std::string const & property_in ):
	property_name_( property_in ),
	property_( ResidueProperties::get_property_from_string( property_in ) /*Will default to NO_PROPERTY if this is an on-the-fly property*/ )
{}

bool
DeleteProperty::apply( MutableResidueType & rsd ) const
{
	runtime_assert_string_msg( property_ != core::chemical::NO_PROPERTY, "Error in DeleteProperty::apply(): Could not parse \"" + property_name_ + "\" as a valid residue type property." );
	rsd.delete_property( property_ );
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "DeleteProperty::apply: " << property_name_ << std::endl;
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteProperty").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteProperty::name() const {
	return "DeleteProperty";
}

// DeleteVariantType //////////////////////////////////////////////////////////
DeleteVariantType::DeleteVariantType( std::string const & variant_in ) :
	variant_str_( variant_in ),
	variant_( ResidueProperties::get_variant_from_string( variant_in ) )
{}

// DeleteVariantType //////////////////////////////////////////////////////////
DeleteVariantType::DeleteVariantType( VariantType const variant_in ) :
	variant_str_( ResidueProperties::get_string_from_variant(variant_in ) ),
	variant_( variant_in )
{}

/// @return  false on success
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
DeleteVariantType::apply( MutableResidueType & rsd ) const
{
	if ( variant_ != NO_VARIANT ) {
		rsd.remove_variant_type( variant_ );
	} else {
		rsd.remove_variant_type( variant_str_); //For on-the-fly types.
	}
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("DeleteVariantType").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteVariantType::name() const {
	return "DeleteVariantType";
}

// AddChi ////////////////////////////////////////////////////////////////////
// Constructor for when the chi index is specified.
AddChi::AddChi(Size const & chino_in,
	std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & atom3_in,
	std::string const & atom4_in) :
	no_index_(false),
	chino_( chino_in ),
	atom1_( atom1_in ),
	atom2_( atom2_in ),
	atom3_( atom3_in ),
	atom4_( atom4_in )
{}

// Constructor for when the chi index is not specified.
/// @details In this case, the new chi will be added to the end of the list when this operation is applied.
/// @author  Labonte
AddChi::AddChi(std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & atom3_in,
	std::string const & atom4_in) :
	no_index_(true),
	chino_(),
	atom1_(atom1_in),
	atom2_(atom2_in),
	atom3_(atom3_in),
	atom4_(atom4_in)
{}

bool
AddChi::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) || !rsd.has( atom3_ ) || !rsd.has( atom4_ ) ) {
		TR_PatchOperations.Debug << "AddChi::apply failed: " << rsd.name() << " is missing atom(s) " <<
			atom1_ << ' ' << rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) <<
			atom3_ << ' ' << rsd.has( atom3_ ) << ' ' << atom4_ << ' ' << rsd.has( atom4_ ) << std::endl;
		return true; // failure
	} else {
		if ( no_index_ ) {
			rsd.add_chi(atom1_, atom2_, atom3_, atom4_);
		} else {
			rsd.add_chi(chino_, atom1_ , atom2_, atom3_, atom4_);
		}
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "AddChi::apply: " <<
				atom1_ << ' ' << atom2_ << ' ' << atom3_ << ' ' << atom4_ << std::endl;
		}
		return false;  // success
	}
}

/// @brief Return the name of this PatchOperation ("AddChi").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddChi::name() const {
	return "AddChi";
}

// AddProtonChi //////////////////////////////////////////////////////////////
AddProtonChi::AddProtonChi(
	Size const & chino_in, utility::vector1<core::Real> const & samples, utility::vector1<core::Real> const & extrasamples
):
	chino_( chino_in ), samples_(samples), extrasamples_(extrasamples)
{}

bool
AddProtonChi::apply( MutableResidueType & rsd ) const
{
	rsd.set_proton_chi( chino_, samples_, extrasamples_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("AddProtonChi").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddProtonChi::name() const {
	return "AddProtonChi";
}

RedefineChi::RedefineChi(Size const & chino_in,
	std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & atom3_in,
	std::string const & atom4_in):
	chino_( chino_in ), atom1_( atom1_in ), atom2_( atom2_in ), atom3_( atom3_in ), atom4_( atom4_in )
{}

bool
RedefineChi::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) || !rsd.has( atom3_ ) || !rsd.has( atom4_ ) ) {
		TR_PatchOperations.Debug << "RedefineChi::apply failed: " <<
			rsd.name() << " is missing atom(s) " <<
			atom1_ << ' ' << rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) <<
			atom3_ << ' ' << rsd.has( atom3_ ) << ' ' << atom4_ << ' ' << rsd.has( atom4_ ) << std::endl;
		return true; // failure
	} else {
		rsd.redefine_chi( chino_, atom1_ , atom2_, atom3_, atom4_ );
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "RedefineChi::apply: " <<
				atom1_ << ' ' << atom2_ << ' ' << atom3_ << ' ' << atom4_ << std::endl;
		}
		return false;
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("RedefineChi").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
RedefineChi::name() const {
	return "RedefineChi";
}

bool
DeleteTerminalChi::apply( MutableResidueType & rsd ) const
{
	if ( rsd.nchi() < 1 ) {
		TR_PatchOperations.Debug << "DeleteTerminalChi::apply failed: " <<
			rsd.name() << " has only " << rsd.nchi() << " chis!" << std::endl;
		return true; // failure
	} else {
		rsd.delete_terminal_chi();
		return false;
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteTerminalChi").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteTerminalChi::name() const {
	return "DeleteTerminalChi";
}


DeleteMetalbindingAtom::DeleteMetalbindingAtom(
	std::string const & atom_name
):
	atom_name_( atom_name )
{}

bool
DeleteMetalbindingAtom::apply( MutableResidueType & rsd ) const
{
	rsd.delete_metalbinding_atom( atom_name_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteMetalbindingAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteMetalbindingAtom::name() const {
	return "DeleteMetalbindingAtom";
}

DeleteActCoordAtom::DeleteActCoordAtom(
	std::string const & atom_name
):
	atom_name_( atom_name )
{}

bool
DeleteActCoordAtom::apply( MutableResidueType & rsd ) const
{
	rsd.delete_actcoord_atom( atom_name_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteActCoordAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteActCoordAtom::name() const {
	return "DeleteActCoordAtom";
}

// AddChiRotamer /////////////////////////////////////////////////////////////
// Constructor for when the chi index is specified
AddChiRotamer::AddChiRotamer(Size const & chino_in, Angle const & mean_in, Angle const & sdev_in) :
	no_index_(false),
	chino_(chino_in),
	mean_(mean_in),
	sdev_(sdev_in)
{}

/// @details In this case, the rotamer sample will be added to the last chi when this operation is applied.
/// @author  Labonte
AddChiRotamer::AddChiRotamer(Angle const & mean_in, Angle const & sdev_in) :
	no_index_(true),
	chino_(),
	mean_(mean_in),
	sdev_(sdev_in)
{}

bool
AddChiRotamer::apply( MutableResidueType & rsd ) const
{
	if ( no_index_ ) {
		rsd.add_chi_rotamer_to_last_chi(mean_ , sdev_);
	} else {
		rsd.add_chi_rotamer(chino_, mean_ , sdev_);
	}
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "AddChiRotamer::apply: " <<
			chino_ << ' ' << mean_ << ' ' << sdev_ << std::endl;
	}
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("AddChiRotamer").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddChiRotamer::name() const {
	return "AddChiRotamer";
}

// ClearChiRotamers ///////////////////////////////////////////////////////////

ClearChiRotamers::ClearChiRotamers( core::uint const chi_no_in ) :
	chi_no_( chi_no_in )
{}

/// @return  true on failure
bool
ClearChiRotamers::apply( MutableResidueType & rsd ) const
{
	if ( chi_no_ == 0 || chi_no_ > rsd.nchi() ) {
		TR_PatchOperations.Debug << "ClearChiRotamers::apply failed: " <<
			rsd.name() << " has no chi " << chi_no_ << std::endl;
		return true; // failure
	} else {
		rsd.clear_chi_rotamers( chi_no_ ) ;
		return false;  // success
	}
}

/// @brief Return the name of this PatchOperation ("ClearChiRotamers").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ClearChiRotamers::name() const {
	return "ClearChiRotamers";
}

// AddAtom ///////////////////////////////////////////////////////////////////
AddAtom::AddAtom(
	std::string const & atom_name_in,
	std::string const & atom_type_name_in,
	std::string const & mm_atom_type_name_in,
	Real const charge
):
	atom_name_( atom_name_in ),
	atom_type_name_( atom_type_name_in ),
	mm_atom_type_name_( mm_atom_type_name_in ),
	charge_( charge )
{}

/// @brief Return the name of this PatchOperation ("AddAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddAtom::name() const {
	return "AddAtom";
}

bool
AddAtom::apply( MutableResidueType & rsd ) const
{
	if ( rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "AddAtom::apply failed: " <<
			rsd.name() << " already has an atom named  '" << atom_name_ << "'." << std::endl;
		return true; // Can't add an atom if this residue already has one of the same name.
	}
	rsd.add_atom( atom_name_, atom_type_name_, mm_atom_type_name_, charge_ );
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "AddAtom::apply: " <<
			atom_name_ << ' ' << atom_type_name_ << ' ' << mm_atom_type_name_ << ' ' << charge_ << std::endl;
	}
	return false;
}


// AddAtomAlias ///////////////////////////////////////////////////////////////

AddAtomAlias::AddAtomAlias(
	std::string const & rosetta_atom_name_in,
	utility::vector1< std::string > const & aliases_in,
	std::string const & preferred_alias_in ) :
	rosetta_atom_name_( rosetta_atom_name_in ),
	aliases_( aliases_in ),
	preferred_alias_( preferred_alias_in )
{}

/// @return  true on failure
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
AddAtomAlias::apply( MutableResidueType & rsd ) const
{
	for ( auto alias : aliases_ ) {
		rsd.add_atom_alias( rosetta_atom_name_, alias );
	}
	rsd.add_canonical_atom_alias( rosetta_atom_name_, preferred_alias_ );
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("AddAtomAlias").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddAtomAlias::name() const {
	return "AddAtomAlias";
}


// AddBond ////////////////////////////////////////////////////////////////////

AddBond::AddBond(
	std::string const & atom1_in,
	std::string const & atom2_in
):
	atom1_( atom1_in ),
	atom2_( atom2_in )
{}

/// @brief Return the name of this PatchOperation ("AddBond").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddBond::name() const { return "AddBond"; }

bool
AddBond::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) ) {
		TR_PatchOperations.Debug << "AddBond::apply failed: " <<
			rsd.name() << " is missing atom(s) " <<
			atom1_ << ' ' << rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) << std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "AddBond::apply: " << atom1_ << ' ' << atom2_ << std::endl;
		}
		rsd.add_bond( atom1_, atom2_ );
	}
	return false;
}


// AddBondType ////////////////////////////////////////////////////////////////

AddBondType::AddBondType(
	std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & bond_type_in ) :
	atom1_( atom1_in ),
	atom2_( atom2_in ),
	bond_type_( bond_type_in )
{}

/// @return  true on failure
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
AddBondType::apply( MutableResidueType & rsd ) const
{
	rsd.add_bond( atom1_, atom2_, convert_to_BondName( bond_type_ ) );
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("AddBondType").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddBondType::name() const {
	return "AddBondType";
}

// ChangeBondType /////////////////////////////////////////////////////////////

ChangeBondType::ChangeBondType(
	std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & old_bond_type_in,
	std::string const & new_bond_type_in ) :
	atom1_( atom1_in ),
	atom2_( atom2_in ),
	old_bond_type_( old_bond_type_in ),
	new_bond_type_( new_bond_type_in )
{}

/// @return  true on failure
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
ChangeBondType::apply( MutableResidueType & rsd ) const
{
	rsd.change_bond_type(
		atom1_, atom2_, convert_to_BondName( new_bond_type_ ) );
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("ChangeBondType").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ChangeBondType::name() const {
	return "ChangeBondType";
}

// SetAtomicCharge ////////////////////////////////////////////////////////////

SetAtomicCharge::SetAtomicCharge(
	std::string const & atom_name_in,
	Real const charge_in
):
	atom_name_( atom_name_in ),
	charge_( charge_in )
{}

bool
SetAtomicCharge::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetAtomicCharge::apply failed: " <<
			rsd.name() << " is missing atom: " << atom_name_ << std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetAtomicCharge::apply: " << atom_name_ << ' ' << charge_ << std::endl;
		}
		rsd.atom( atom_name_ ).charge( charge_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetAtomicCharge").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetAtomicCharge::name() const {
	return "SetAtomicCharge";
}

// SetFormalCharge ////////////////////////////////////////////////////////////

SetFormalCharge::SetFormalCharge( std::string const & atom_name_in, core::SSize charge_in ) :
	atom_name_( atom_name_in ),
	charge_( charge_in )
{}

/// @return  true on failure
bool
SetFormalCharge::apply( MutableResidueType & rsd ) const
{
	if ( ! rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Error << "SetFormalCharge::apply() failed: " << rsd.name() << " is missing atom: " <<
			atom_name_ << std::endl;
		return true;  // failure
	}
	rsd.atom( atom_name_ ).formal_charge( charge_ );
	return false;  // success
}

/// @brief Return the name of this PatchOperation ("SetFormalCharge").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetFormalCharge::name() const {
	return "SetFormalCharge";
}

// SetNetFormalCharge ////////////////////////////////////////////////////////////
/// @brief Constructor
/// @details Note that this deliberately takes a signed int.
SetNetFormalCharge::SetNetFormalCharge( signed int const charge_in ) :
	charge_( charge_in )
{}

/// @return  true on failure
bool
SetNetFormalCharge::apply( MutableResidueType & rsd ) const
{
	rsd.net_formal_charge( charge_ );
	return false;  // success always
}

/// @brief Return the name of this PatchOperation ("SetNetFormalCharge").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetNetFormalCharge::name() const {
	return "SetNetFormalCharge";
}

// SetAtomType ////////////////////////////////////////////////////////////////

SetAtomType::SetAtomType(
	std::string const & atom_name_in,
	std::string const & atom_type_name_in
):
	atom_name_( atom_name_in ),
	atom_type_name_( atom_type_name_in )
{}

bool
SetAtomType::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetAtomType::apply failed: " <<
			rsd.name() << " is missing atom: " << atom_name_ << std::endl;
		return true; // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetAtomType::apply: " << atom_name_ << ' ' << atom_type_name_ << std::endl;
		}
		rsd.set_atom_type( atom_name_, atom_type_name_ );
		if ( atom_type_name_ == "VIRT" ) {
			rsd.atom( atom_name_ ).charge( 0.0 );
			rsd.atom( atom_name_ ).is_virtual( true );
		}
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetAtomType").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetAtomType::name() const {
	return "SetAtomType";
}

// Set_AA ////////////////////////////////////////////////////////////////////
bool
Set_AA::apply( MutableResidueType & rsd ) const
{
	rsd.aa( aa_ );
	return false;
}

// SetIO_String //////////////////////////////////////////////////////////////

SetIO_String::SetIO_String(
	std::string const & name3,
	char const name1
):
	name3_( name3 ),
	name1_( name1 )
{}

bool
SetIO_String::apply( MutableResidueType & rsd ) const
{
	rsd.name3( name3_ );
	rsd.name1( name1_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("SetIO_String").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetIO_String::name() const {
	return "SetIO_String";
}

SetInterchangeabilityGroup_String::SetInterchangeabilityGroup_String(
	std::string const & intgrp
):
	intgrp_( intgrp )
{}

bool
SetInterchangeabilityGroup_String::apply( MutableResidueType & rsd ) const
{
	rsd.interchangeability_group( intgrp_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("SetInterchangeabilityGroup_String").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetInterchangeabilityGroup_String::name() const {
	return "SetInterchangeabilityGroup_String";
}

// SetMMAtomType //////////////////////////////////////////////////////////////

SetMMAtomType::SetMMAtomType(
	std::string const & atom_name_in,
	std::string const & mm_atom_type_name_in
):
	atom_name_( atom_name_in ),
	mm_atom_type_name_( mm_atom_type_name_in )
{}

bool
SetMMAtomType::apply( MutableResidueType & rsd ) const
{
	if ( ! rsd.has( atom_name_ ) ) {
		if ( TR_PatchOperations.Debug.visible() ) {
			TR_PatchOperations.Debug << "SetMMAtomType::apply failed: " <<
				rsd.name() << " is missing atom: " << atom_name_ << std::endl;
		}
		return true;  // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetMMAtomType::apply: " <<
				atom_name_ << ' ' << mm_atom_type_name_ << std::endl;
		}
		rsd.atom( rsd.atom_vertex( atom_name_ ) ).mm_name( mm_atom_type_name_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetMMAtomType").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetMMAtomType::name() const {
	return "SetMMAtomType";
}

// SetICoor //////////////////////////////////////////////////////////////////

// helper function
std::string
expand_icoor_atom_name( std::string name, MutableResidueType const & rsd )
{
	std::string const nconn_tag( "%LASTCONN" );
	Size pos( name.find( nconn_tag ) );
	if ( pos < name.size() ) {
		name.replace( pos, nconn_tag.size(), ObjexxFCL::string_of( rsd.n_possible_residue_connections() ) );
	}
	return name;
}

// helper function
void
check_residue_has_atom( std::string const & atom_name, core::chemical::MutableResidueType const & rsd, std::string const & location ) {
	// The strings here might be connection designations, so avoid the inevitable failure from that.
	if ( string_to_icoord_type( atom_name ) == ICoordAtomIDType::INTERNAL &&
			! rsd.has( atom_name ) ) {
		utility_exit_with_message( location + " atom `" + atom_name + "` not found in residue " + rsd.name() );
	}
}

SetICoor::SetICoor(
	std::string const & atom_in,
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	std::string const & stub1_in,
	std::string const & stub2_in,
	std::string const & stub3_in
):
	atom_( atom_in ),
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub1_( stub1_in ),
	stub2_( stub2_in ),
	stub3_( stub3_in )
{}

/// @brief Return the name of this PatchOperation ("SetICoor").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetICoor::name() const { return "SetICoor"; }

bool
SetICoor::apply( MutableResidueType & rsd ) const
{
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "SetICoor::apply: " <<
			atom_ << ' ' << stub1_ << ' ' << stub2_ << ' ' << stub3_ << std::endl;
	}
	std::string const atom ( expand_icoor_atom_name( atom_ , rsd ) );
	std::string const stub1( expand_icoor_atom_name( stub1_, rsd ) );
	std::string const stub2( expand_icoor_atom_name( stub2_, rsd ) );
	std::string const stub3( expand_icoor_atom_name( stub3_, rsd ) );
	check_residue_has_atom( atom, rsd, "SetICoor patch:");
	check_residue_has_atom( stub1, rsd, "SetICoor patch:");
	check_residue_has_atom( stub2, rsd, "SetICoor patch:");
	check_residue_has_atom( stub3, rsd, "SetICoor patch:");
	bool const rebuild_icoor_xyz( true );
	rsd.set_icoor( atom, phi_, theta_, d_, stub1, stub2, stub3, rebuild_icoor_xyz );
	return false;
}

ChangeAncestory::ChangeAncestory(
	std::string const & target_atom,
	Ancestor which_ancestor,
	std::string const & ancestor_name
) :
	atom_( target_atom ),
	which_ancestor_( which_ancestor ),
	ancestor_name_( ancestor_name )
{}

/// @brief change the ancestory, but leave the icoors intact.
bool
ChangeAncestory::apply( MutableResidueType & rsd ) const
{
	VD const atvd( rsd.atom_vertex( atom_ ));
	MutableICoorRecordCOP const aticoor( rsd.icoor( atvd ) );

	if ( aticoor == nullptr ) {
		utility_exit_with_message("ChangeAncestory patch called on a ResidueType and atom that doesn't yet have a proper Icoor record");
	}

	std::string pa, gp, gg;

	try {
		pa = which_ancestor_ == anc_parent           ? ancestor_name_ : aticoor->stub_atom1();
		gp = which_ancestor_ == anc_grandparent      ? ancestor_name_ : aticoor->stub_atom2();
		gg = which_ancestor_ == anc_greatgrandparent ? ancestor_name_ : aticoor->stub_atom3();
	} catch ( utility::excn::Exception const & excn ) {
		std::ostringstream oss;
		oss << "Failed to apply the ChangeAncestory patch (" << atom_ << ", ";
		switch ( which_ancestor_ ) {
		case anc_parent : oss << "parent"; break;
		case anc_grandparent : oss << "grandparent"; break;
		case anc_greatgrandparent : oss << "greatgrandparent"; break;
		}
		oss << ", " << ancestor_name_ << ") to residue type " << rsd.name() << " when constructing ancetor's ICoorAtomID\n";
		oss << "Message from ICoorAtomID constructor: " << excn.msg() << "\n";

		utility_exit_with_message( oss.str() );
	}

	rsd.set_icoor( atom_, aticoor->phi(), aticoor->theta(), aticoor->d(), pa, gp, gg, true /*rebuild_xyz*/ );
	return false;
}

/// @brief Return the name of this PatchOperation ("ChangeAncestory").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ChangeAncestory::name() const {
	return "ChangeAncestory";
}

// ResetBondLength ////////////////////////////////////////////////////////////

ResetBondLength::ResetBondLength( std::string const & atm_in, core::Distance d_in ) :
	atm_( atm_in ),
	d_( d_in )
{}

/// @return  true on failure
bool
ResetBondLength::apply( MutableResidueType & rsd ) const
{
	if ( ! rsd.has( atm_ ) ) {
		TR_PatchOperations.Debug << "ResetBondLength::apply failed: " <<
			rsd.name() << " is missing atom " << atm_ << std::endl;
		return true; // failure
	} else {
		rsd.reset_bond_distance_to_atom( atm_, d_ );
		return false;  // success
	}
}

/// @brief Return the name of this PatchOperation ("ResetBondLength").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ResetBondLength::name() const {
	return "ResetBondLength";
}

// PrependMainchainAtom ///////////////////////////////////////////////////////

PrependMainchainAtom::PrependMainchainAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
PrependMainchainAtom::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "PrependMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		utility::vector1< VD > const & old_mainchain_atoms( rsd.mainchain_atoms() );
		utility::vector1< VD > new_mainchain_atoms;
		new_mainchain_atoms.push_back( rsd.atom_vertex( atom_name_ ) );
		for ( Size i = 1; i <= old_mainchain_atoms.size(); ++i ) {
			new_mainchain_atoms.push_back( old_mainchain_atoms[i] );
		}
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("PrependMainchainAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
PrependMainchainAtom::name() const {
	return "PrependMainchainAtom";
}

AppendMainchainAtom::AppendMainchainAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
AppendMainchainAtom::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "AppendMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		utility::vector1< VD > new_mainchain_atoms( rsd.mainchain_atoms() );
		new_mainchain_atoms.push_back( rsd.atom_vertex( atom_name_ ) );
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("AppendMainchainAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AppendMainchainAtom::name() const {
	return "AppendMainchainAtom";
}

ReplaceMainchainAtom::ReplaceMainchainAtom( std::string const & target, std::string const & new_atom ) :
	target_( target ),
	new_atom_( new_atom )
{}

bool
ReplaceMainchainAtom::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( target_ ) ) {
		TR_PatchOperations.Debug << "ReplaceMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << target_ << std::endl;
		return true; // failure
	} else {
		utility::vector1< VD > new_mainchain_atoms( rsd.mainchain_atoms() );
		for ( Size i = 1; i <= new_mainchain_atoms.size(); ++i ) {
			std::string mainchain_atom_name = rsd.atom_name( new_mainchain_atoms[ i ] );
			mainchain_atom_name.erase( std::remove_if( mainchain_atom_name.begin(),
				mainchain_atom_name.end(),
				::isspace ),
				mainchain_atom_name.end() );

			if ( mainchain_atom_name == target_ ) {
				new_mainchain_atoms[ i ] = rsd.atom_vertex( new_atom_ );
				break;
			}
		}
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceMainchainAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceMainchainAtom::name() const {
	return "ReplaceMainchainAtom";
}

SetNbrAtom::SetNbrAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
SetNbrAtom::apply( MutableResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "SetNbrAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		rsd.nbr_atom( atom_name_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetNbrAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetNbrAtom::name() const {
	return "SetNbrAtom";
}

SetNbrRadius::SetNbrRadius( Real const & radius ) :
	radius_( radius )
{}

bool
SetNbrRadius::apply( MutableResidueType & rsd ) const
{
	rsd.nbr_radius( radius_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("SetNbrRadius").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetNbrRadius::name() const {
	return "SetNbrRadius";
}

/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
bool
SetAllAtomsRepulsive::apply( MutableResidueType & rsd ) const {

	for ( VD atm: rsd.all_atoms() ) {
		if ( rsd.atom(atm).is_hydrogen() ) {
			rsd.set_atom_type( atm, "REPLS" );
		} else {
			rsd.set_atom_type( atm, "HREPS" );
		}
	}

	return false;
}

/// @brief Return the name of this PatchOperation ("SetAllAtomsRepulsive").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetAllAtomsRepulsive::name() const {
	return "SetAllAtomsRepulsive";
}

SetOrientAtom::SetOrientAtom(bool force_nbr_atom_orient):
	force_nbr_atom_orient_(force_nbr_atom_orient)
{}

bool
SetOrientAtom::apply( MutableResidueType & rsd ) const
{
	rsd.force_nbr_atom_orient( force_nbr_atom_orient_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("SetOrientAtom").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetOrientAtom::name() const {
	return "SetOrientAtom";
}

/// @brief Constructor.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
RemoveRotamerSpecifications::RemoveRotamerSpecifications() {}

/// @brief Strip all RotamerSpecifications from the ResidueType.
/// @return Success ("false" -- i.e. no error thrown) or failure ("true").
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
RemoveRotamerSpecifications::apply( MutableResidueType & rsd ) const {
	using namespace core::chemical::rotamers;

	if ( rsd.rotamer_library_specification() ) {
		rsd.strip_rotamer_library_specification();
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("RemoveRotamerSpecifications").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
RemoveRotamerSpecifications::name() const {
	return "RemoveRotamerSpecifications";
}

/// @brief Constructor.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
RamaPreproFilename::RamaPreproFilename(
	std::string const &non_prepro_file,
	std::string const &prepro_file
) :
	non_prepro_file_( non_prepro_file ),
	prepro_file_( prepro_file )
{}

/// @brief Set the RamaPrepro library paths in the residue type.
/// @return Success ("false" -- i.e. no error thrown) or failure ("true").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
RamaPreproFilename::apply(
	MutableResidueType & rsd
) const {
	rsd.set_rama_prepro_map_file_name( non_prepro_file_, false);
	rsd.set_rama_prepro_map_file_name( prepro_file_, true);
	return false;
}

/// @brief Return the name of this PatchOperation ("RamaPreproFilename").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
RamaPreproFilename::name() const {
	return "RamaPreproFilename";
}

/// @brief Constructor.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
RamaPreproResname::RamaPreproResname(
	std::string const &resname_in
) :
	resname_( resname_in )
{}

/// @brief Set the RamaPrepro reference string in the residue type.
/// @return Success ("false" -- i.e. no error thrown) or failure ("true").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
RamaPreproResname::apply(
	MutableResidueType & rsd
) const {
	rsd.set_rama_prepro_mainchain_torsion_potential_name( resname_, false );
	rsd.set_rama_prepro_mainchain_torsion_potential_name( resname_, true );
	return false;
}

/// @brief Return the name of this PatchOperation ("RamaPreproResname").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
RamaPreproResname::name() const {
	return "RamaPreproResname";
}

NCAARotLibPath::NCAARotLibPath( std::string const & path_in ) :
	path_( path_in )
{}

/// @brief set the NCAA rotamer library path in the residue type
bool
NCAARotLibPath::apply( MutableResidueType & rsd ) const
{
	using namespace core::chemical::rotamers;

	NCAARotamerLibrarySpecificationOP ncaa_libspec;
	if ( rsd.rotamer_library_specification() ) {
		NCAARotamerLibrarySpecificationCOP old_libspec(
			utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification const >(
			rsd.rotamer_library_specification() ) );
		if ( ! old_libspec ) {
			tr.Error << "Found existing rotamer specification " << rsd.rotamer_library_specification()->keyname() <<
				" when attempting to patch NCAA rotlib path." << std::endl;
			utility_exit_with_message( "Cannot have multiple rotamer specifications." );
		}
		ncaa_libspec = utility::pointer::make_shared< NCAARotamerLibrarySpecification >( *old_libspec );
	} else {
		ncaa_libspec = NCAARotamerLibrarySpecificationOP(
			new core::chemical::rotamers::NCAARotamerLibrarySpecification );
	}
	ncaa_libspec->ncaa_rotlib_path( path_ );
	rsd.rotamer_library_specification( ncaa_libspec );
	return false;
}

/// @brief Return the name of this PatchOperation ("NCAARotLibPath").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
NCAARotLibPath::name() const {
	return "NCAARotLibPath";
}

/// @brief Constructor.
NCAARotLibBBTorsions::NCAARotLibBBTorsions( utility::vector1< core::Size > const &indices_in ) :
	indices_( indices_in )
{}

/// @brief Set the mainchain torsion indices that a noncanonical rotamer library depends upon.
bool
NCAARotLibBBTorsions::apply( MutableResidueType & rsd ) const
{
	using namespace core::chemical::rotamers;

	NCAARotamerLibrarySpecificationOP ncaa_libspec;
	if ( rsd.rotamer_library_specification() ) {
		NCAARotamerLibrarySpecificationCOP old_libspec(
			utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification const >(
			rsd.rotamer_library_specification() ) );
		if ( ! old_libspec ) {
			tr.Error << "Found existing rotamer specification " << rsd.rotamer_library_specification()->keyname() <<
				" when attempting to patch NCAA rotlib path." << std::endl;
			utility_exit_with_message( "Cannot have multiple rotamer specifications." );
		}
		ncaa_libspec = utility::pointer::make_shared< NCAARotamerLibrarySpecification >( *old_libspec );
	} else {
		ncaa_libspec = NCAARotamerLibrarySpecificationOP(
			new core::chemical::rotamers::NCAARotamerLibrarySpecification );
	}
	ncaa_libspec->clear_rotamer_bb_torsion_indices();

	for ( core::Size i(1), imax(indices_.size()); i<=imax; ++i ) {
		ncaa_libspec->add_rotamer_bb_torsion_index(indices_[i]);
	}

	rsd.rotamer_library_specification(ncaa_libspec);

	return false;
}

/// @brief Return the name of this PatchOperation ("NCAARotLibBBTorsions").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
NCAARotLibBBTorsions::name() const {
	return "NCAARotLibBBTorsions";
}

/// @brief Constructor
/// @author Vikram K. Mulligan (vmullig@uw.edu).
NCAARotLibNumRotamerBins::NCAARotLibNumRotamerBins(
	utility::vector1< core::Size > const &binsizes_in
) :
	binsizes_( binsizes_in )
{}

/// @brief Set the number of rotamer bins per chi for an NCAA.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
NCAARotLibNumRotamerBins::apply(
	MutableResidueType & rsd
) const {
	using namespace core::chemical::rotamers;

	NCAARotamerLibrarySpecificationOP ncaa_libspec( utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification >( rsd.rotamer_library_specification_nonconst() ) );
	runtime_assert_string_msg( ncaa_libspec, "A patch may only set the number of NCAA rotamer bins (\"NCAA_ROTLIB_NUM_ROTAMER_BINS\") for a residue type that already has an NCAA rotamer library specified, either in the params file or in patches." );
	ncaa_libspec->ncaa_rotlib_n_bin_per_rot( binsizes_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("NCAARotLibNumRotamerBins").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
NCAARotLibNumRotamerBins::name() const {
	return "NCAARotLibNumRotamerBins";
}

/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
bool
ConnectSulfurAndMakeVirtualProton::apply( MutableResidueType & rsd ) const {
	std::string const & disulfide_atom_name = rsd.get_disulfide_atom_name();
	std::string const & CB_equivalent = rsd.atom_name( rsd.atom_base( rsd.atom_vertex( disulfide_atom_name ) ) );
	std::string const & CA_equivalent = rsd.atom_name( rsd.atom_base( rsd.atom_vertex( CB_equivalent ) ) );
	AddConnect ad(
		disulfide_atom_name,
		180, 68.374, 1.439,
		disulfide_atom_name,
		CB_equivalent,
		CA_equivalent
	);
	bool x = ad.apply( rsd );

	// Now I need to grab the proton sitting on the disulfide atom
	utility::vector1< VD > const & bonded_H( rsd.bonded_hydrogens( rsd.atom_vertex( disulfide_atom_name ) ) );
	if ( bonded_H.empty() ) {
		utility_exit_with_message("Cannot apply ConnectSulfurAndMakeVirtualProton to " + rsd.name() + ": atom " + disulfide_atom_name + " has no hydrogens.");
	}
	if ( bonded_H.size() != 1 ) {
		tr.Warning << "In ConnectSulfurAndMakeVirtualProton patch, " << rsd.name() << " disulfide atom " << disulfide_atom_name << " has more than one hydrogen." << std::endl;
	}
	VD HG_atm = bonded_H[1];

	rsd.set_atom_type( HG_atm, "VIRT" );
	rsd.atom( HG_atm ).mm_name( "VIRT" );
	rsd.atom( HG_atm ).charge( 0.0 );
	rsd.atom( HG_atm ).is_virtual( true );

	return x;
}

bool
ConnectSulfurAndMakeVirtualProton::changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const
{
	std::string const & disulfide_atom_name = rsd_type.get_disulfide_atom_name();
	return rsd_type.has( atom ) && rsd_type.atom_index( atom ) == rsd_type.atom_index( disulfide_atom_name );
}

/// @brief Return the name of this PatchOperation ("ConnectSulfurAndMakeVirtualProton").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ConnectSulfurAndMakeVirtualProton::name() const {
	return "ConnectSulfurAndMakeVirtualProton";
}

/// @brief Constructor.
SetBaseName::SetBaseName(
	std::string const & new_base_name
) :
	new_base_name_(new_base_name)
{}

/// @brief Apply this patch to generate a new base residue type.
bool
SetBaseName::apply(
	MutableResidueType & rsd
) const {
	rsd.name( new_base_name_ );
	rsd.base_name( new_base_name_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("SetBaseName").
std::string
SetBaseName::name() const {
	return "SetBaseName";
}

/// @brief This patch operaton DOES result in a new base residue type.
bool
SetBaseName::generates_base_residue_type() const {
	return true;
}

bool
ChiralFlipNaming::apply( MutableResidueType & rsd ) const {

	//Special-case logic for peptoids (VKM 15 Nov 2017):
	if ( rsd.is_peptoid() ) {
		runtime_assert_string_msg( rsd.is_s_peptoid(), "Error in ChiralFlipNaming patch operation: this operation can only be applied to peptoids with side-chains with S chirality." );
		int peptoid_number( std::atoi( rsd.base_name().c_str() ) );
		runtime_assert_string_msg( peptoid_number > 600 && peptoid_number < 700, "Error in ChiralFlipNaming patch operation: S-chirality peptoids must have 3-letter codes between 600 and 699" );
		runtime_assert_string_msg( peptoid_number % 2 == 1, "Error in ChiralFlipNaming patch operation: S-chirality peptoids must have odd-numbered 3-letter codes starting with 6." );
		++peptoid_number;
		std::stringstream new_num;
		new_num << peptoid_number;
		rsd.name3( new_num.str() );
		rsd.interchangeability_group( new_num.str() );
		rsd.name( new_num.str() ); //VariantTypes will be added to this in the Patch::apply() function.
		rsd.base_name( new_num.str() );
		return false;
	}

	//Flip the backbone aa:
	if ( core::chemical::is_canonical_L_aa_excluding_gly( rsd.backbone_aa() ) ) {
		rsd.backbone_aa( core::chemical::get_D_equivalent( rsd.backbone_aa() ) );
	} else if ( core::chemical::is_canonical_D_aa( rsd.backbone_aa() ) ) {
		rsd.backbone_aa( core::chemical::get_L_equivalent( rsd.backbone_aa() ) );
	}

	if ( rsd.aa() <= aa_tyr ) {
		rsd.aa( get_D_equivalent( rsd.aa() ) );

		if ( rsd.aa() == aa_dal ) {
			rsd.name3( "DAL" );
			rsd.interchangeability_group( "DAL" );
		} else if ( rsd.aa() == aa_dcs ) {
			//( igroup == "DCYS" || igroup == "CYS" ) {
			rsd.name3( "DCS" );
			rsd.interchangeability_group( "DCS" );
		} else if ( rsd.aa() == aa_das ) {
			//( igroup == "DASP" || igroup == "ASP" ) {
			rsd.name3( "DAS" );
			rsd.interchangeability_group( "DAS" );
		} else if ( rsd.aa() == aa_dgu ) {
			//( igroup == "DGLU" || igroup == "GLU" ) {
			rsd.name3( "DGU" );
			rsd.interchangeability_group( "DGU" );
		} else if ( rsd.aa() == aa_dph ) {
			//( igroup == "DPHE" || igroup == "PHE" ) {
			rsd.name3( "DPH" );
			rsd.interchangeability_group( "DPH" );
		} else if ( rsd.aa() == aa_dhi ) {
			//( igroup == "DHIS" || igroup == "HIS" ) {
			rsd.name3( "DHI" );
			rsd.interchangeability_group( "DHI" );
		} else if ( rsd.aa() == aa_dil ) {
			//( igroup == "DILE" || igroup == "ILE" ) {
			rsd.name3( "DIL" );
			rsd.interchangeability_group( "DIL" );
		} else if ( rsd.aa() == aa_dly ) {
			//( igroup == "DLYS" || igroup == "LYS" ) {
			rsd.name3( "DLY" );
			rsd.interchangeability_group( "DLY" );
		} else if ( rsd.aa() == aa_dle ) {
			//( igroup == "DLEU" || igroup == "LEU" ) {
			rsd.name3( "DLE" );
			rsd.interchangeability_group( "DLE" );
		} else if ( rsd.aa() == aa_dme ) {
			//( igroup == "DMET" || igroup == "MET" ) {
			rsd.name3( "DME" );
			rsd.interchangeability_group( "DME" );
		} else if ( rsd.aa() == aa_dan ) {
			//( igroup == "DASN" || igroup == "ASN" ) {
			rsd.name3( "DAN" );
			rsd.interchangeability_group( "DAN" );
		} else if ( rsd.aa() == aa_dpr ) {
			//( igroup == "DPRO" || igroup == "PRO" ) {
			if ( rsd.name3() == "HYP" ) {
				rsd.name3( "DHY" );
				rsd.interchangeability_group( "DHY" );
			} else if ( rsd.name3() == "0AZ" ) {
				rsd.name3( "D0A" );
				rsd.interchangeability_group( "D0A" );
			} else {
				rsd.name3( "DPR" );
				rsd.interchangeability_group( "DPR" );
			}
		} else if ( rsd.aa() == aa_dgn ) {
			//( igroup == "DGLN" || igroup == "GLN" ) {
			rsd.name3( "DGN" );
			rsd.interchangeability_group( "DGN" );
		} else if ( rsd.aa() == aa_dar ) {
			//( igroup == "DARG" || igroup == "ARG" ) {
			rsd.name3( "DAR" );
			rsd.interchangeability_group( "DAR" );
		} else if ( rsd.aa() == aa_dse ) {
			//( igroup == "DSER" || igroup == "SER" ) {
			rsd.name3( "DSE" );
			rsd.interchangeability_group( "DSE" );
		} else if ( rsd.aa() == aa_dth ) {
			//( igroup == "DTHR" || igroup == "THR" ) {
			rsd.name3( "DTH" );
			rsd.interchangeability_group( "DTH" );
		} else if ( rsd.aa() == aa_dva ) {
			//( igroup == "DVAL" || igroup == "VAL" ) {
			rsd.name3( "DVA" );
			rsd.interchangeability_group( "DVA" );
		} else if ( rsd.aa() == aa_dtr ) {
			//( igroup == "DTRP" || igroup == "TRP" ) {
			rsd.name3( "DTR" );
			rsd.interchangeability_group( "DTR" );
		} else if ( rsd.aa() == aa_dty ) {
			//( igroup == "DTYR" || igroup == "TYR" ) {
			rsd.name3( "DTY" );
			rsd.interchangeability_group( "DTY" );
		}
	} else if ( rsd.aa() == na_rad ) {
		rsd.aa( na_lra );
		rsd.name3( " 0A" );
		rsd.interchangeability_group( " 0A" );
	} else if ( rsd.aa() == na_rcy ) {
		rsd.aa( na_lrc );
		rsd.name3( " 0C" );
		rsd.interchangeability_group( " 0C" );
	} else if ( rsd.aa() == na_rgu ) {
		rsd.aa( na_lrg );
		rsd.name3( " 0G" );
		rsd.interchangeability_group( " 0G" );
	} else if ( rsd.aa() == na_ura ) {
		rsd.aa( na_lur );
		rsd.name3( " 0U" );
		rsd.interchangeability_group( " 0U" );
	} else {
		// e.g. DHLU, DC01--will be caught later and set to DAL etc. for canonicals
		rsd.interchangeability_group( "D" + rsd.name().substr( 0, rsd.name().find(":") ) );
	}
	//tr << rsd.aa() << std::endl;
	return false;
}

/// @brief Return the name of this PatchOperation ("ChiralFlipNaming").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ChiralFlipNaming::name() const {
	return "ChiralFlipNaming";
}

/// @brief This patch operaton DOES result in a new base residue type.
bool
ChiralFlipNaming::generates_base_residue_type() const {
	return true;
}

bool
ChiralFlipAtoms::apply( MutableResidueType & rsd ) const {

	// So I think I need to separately process shadow atoms.
	// The reason is because they're 0 from something
	// so I think maybe I should delete and re-add them to a residue?

	//tr << "In chiralflipatoms apply " << std::endl;

	// Flip each
	//rsd.debug_dump_icoor();
	if ( rsd.natoms() < 4 ) return true;

	//std::string base_atom_name;
	// The rest of the Patch takes care of props, and it's really out of theme.
	if ( rsd.is_l_aa() || rsd.is_s_peptoid() ) {
		//rsd.delete_property( "L_AA" );
		//rsd.add_property( "D_AA" );
		//base_atom_name = "N";
	} else if ( rsd.is_d_rna() ) {
		//rsd.delete_property( "D_RNA" );
		//rsd.add_property( "L_RNA" );
		//base_atom_name = "P";
	} else {
		utility_exit_with_message( "For some reason, calling ChiralFlipAtoms on a residue type that's not an L-amino acid, not a D-nucleic acid, and not a peptoid with a chiral sidechain." );
	}

	for ( VD atm: rsd.all_atoms() ) {
		MutableICoorRecordCOP icoor1( rsd.icoor( atm ) );

		if ( icoor1 == nullptr ) {
			tr << "In ChiralFlipAtoms, cannot flip ICoor for atom " << rsd.atom_name( atm ) << " as its ICoor record isn't set yet." << std::endl;
			continue; // Or should this be an error?
		}

		if ( icoor1->stub_atom1() == icoor1->stub_atom2() ) { continue; } // Don't flip root

		// RM: There used to be some atom name mangling here, but I don't think it's needed any longer, given we're keeping the atoms as-is

		rsd.set_icoor( rsd.atom_name( atm ),
			-1.0*icoor1->phi(),
			icoor1->theta(),
			icoor1->d(),
			icoor1->stub_atom1(),
			icoor1->stub_atom2(),
			icoor1->stub_atom3(),
			true );
	}
	rsd.fill_ideal_xyz_from_icoor();

	return false;
}

/// @brief Return the name of this PatchOperation ("ChiralFlipAtoms").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ChiralFlipAtoms::name() const {
	return "ChiralFlipAtoms";
}

/// @brief replace proton with trifluoromethyl
bool
ReplaceProtonWithTrifluoromethyl::apply( MutableResidueType & rsd ) const {
	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	//std::cout << "Proton is " << rsd.atom_name( proton_index ) << std::endl;
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	//std::cout << "Adding atom C" << atom_ << " to " << atom_ << ", bonded to " << rsd.atom_name( icoor.stub_atom1().atomno() ) << rsd.atom_name( icoor.stub_atom2().atomno() ) << rsd.atom_name( icoor.stub_atom3().atomno() ) << std::endl;
	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.511005,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	// these codes with
	rsd.add_atom( "E"+atom_, "F", "F3", -0.20 );
	rsd.add_atom( "F"+atom_, "F", "F3", -0.20 );
	rsd.add_atom( "G"+atom_, "F", "F3", -0.20 );
	rsd.add_bond( "C"+atom_, "E"+atom_ );
	rsd.add_bond( "C"+atom_, "F"+atom_ );
	rsd.add_bond( "C"+atom_, "G"+atom_ );

	rsd.set_icoor( "E"+atom_,
		88.407090/180.0*3.14159,
		70.5/180.0*3.14159,
		1.319661,
		"C"+atom_,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		true );

	rsd.set_icoor( "F"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.319661,
		"C"+atom_,
		icoor->stub_atom1(),
		"E"+atom_,
		true );

	rsd.set_icoor( "G"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.319661,
		"C"+atom_,
		icoor->stub_atom1(),
		"F"+atom_,
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithTrifluoromethyl").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithTrifluoromethyl::name() const {
	return "ReplaceProtonWithTrifluoromethyl";
}

/// @brief replace proton with methyl
bool
ReplaceProtonWithMethyl::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.511005,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	// these codes with
	rsd.add_atom( "X"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_atom( "Y"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_atom( "Z"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_bond( "C"+atom_, "X"+atom_ );
	rsd.add_bond( "C"+atom_, "Y"+atom_ );
	rsd.add_bond( "C"+atom_, "Z"+atom_ );

	rsd.set_icoor( "X"+atom_,
		88.407090/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		true );

	rsd.set_icoor( "Y"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		icoor->stub_atom1(),
		"X"+atom_,
		true );

	rsd.set_icoor( "Z"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		icoor->stub_atom1(),
		"Y"+atom_,
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithMethyl").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithMethyl::name() const {
	return "ReplaceProtonWithMethyl";
}

/// @brief replace proton with methoxy
bool
ReplaceProtonWithMethoxy::apply( MutableResidueType & rsd ) const {
	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "O"+atom_, "OH", "OE", -0.61 );
	rsd.add_bond( "O"+atom_, atom_ );
	rsd.set_icoor( "O"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.375964, // actually that's the tyr dist
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true  );


	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.22 );
	rsd.add_bond( "C"+atom_, "O"+atom_ );
	rsd.set_icoor( "C"+atom_,
		-177.688373/180.0*3.14159,
		65.681781/180.0*3.14159,
		1.393527, // aliphatic context.
		"O"+atom_,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		true );

	// these codes with
	rsd.add_atom( "X"+atom_, "Hapo", "HA", 0.14 );
	rsd.add_atom( "Y"+atom_, "Hapo", "HA", 0.14 );
	rsd.add_atom( "Z"+atom_, "Hapo", "HA", 0.14 );
	rsd.add_bond( "C"+atom_, "X"+atom_ );
	rsd.add_bond( "C"+atom_, "Y"+atom_ );
	rsd.add_bond( "C"+atom_, "Z"+atom_ );

	rsd.set_icoor( "X"+atom_,
		88.407090/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		"O"+atom_,
		icoor->stub_atom1(),
		true );

	rsd.set_icoor( "Y"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		"O"+atom_,
		"X"+atom_,
		true );

	rsd.set_icoor( "Z"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		"O"+atom_,
		"Y"+atom_,
		true );

	// fake proton chi
	rsd.add_chi( icoor->stub_atom2(), icoor->stub_atom1(), "O"+atom_, "C"+atom_ );

	rsd.add_chi_rotamer( rsd.nchi(),   0, 10 );
	rsd.add_chi_rotamer( rsd.nchi(),  60, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 120, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 180, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 240, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 300, 10 );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithMethoxy").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithMethoxy::name() const {
	return "ReplaceProtonWithMethoxy";
}

/// @brief replace proton with ethyl
bool
ReplaceProtonWithEthyl::apply( MutableResidueType & rsd ) const {
	// What is the proton to be replaced?
	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "C"+atom_, "CH2", "CT2", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.511005,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true  );

	rsd.add_atom( "A"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_atom( "B"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_bond( "C"+atom_, "A"+atom_ );
	rsd.add_bond( "C"+atom_, "B"+atom_ );

	rsd.add_atom( "D"+atom_, "CH3", "CT3", -0.27 );
	rsd.add_bond( "D"+atom_, "C"+atom_ );
	rsd.set_icoor( "D"+atom_,
		3.14159,
		70.5/180.0*3.14159,
		1.527165, // semi-arbitrary--CC bond from ethanolamine
		"C"+atom_,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		true );

	rsd.set_icoor( "A"+atom_,
		120/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		icoor->stub_atom1(),
		"D"+atom_,
		true );

	rsd.set_icoor( "B"+atom_,
		120/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		icoor->stub_atom1(),
		"A"+atom_,
		true );

	// these codes with
	rsd.add_atom( "X"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_atom( "Y"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_atom( "Z"+atom_, "Hapo", "HA", 0.09 );
	rsd.add_bond( "D"+atom_, "X"+atom_ );
	rsd.add_bond( "D"+atom_, "Y"+atom_ );
	rsd.add_bond( "D"+atom_, "Z"+atom_ );

	rsd.set_icoor( "X"+atom_,
		88.407090/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"D"+atom_,
		"C"+atom_,
		icoor->stub_atom1(),
		true );

	rsd.set_icoor( "Y"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"D"+atom_,
		"C"+atom_,
		"X"+atom_,
		true );

	rsd.set_icoor( "Z"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"D"+atom_,
		"C"+atom_,
		"Y"+atom_,
		true );

	// fake proton chi
	rsd.add_chi( icoor->stub_atom2(), icoor->stub_atom1(), "C"+atom_, "D"+atom_ );

	rsd.add_chi_rotamer( rsd.nchi(),   0, 10 );
	rsd.add_chi_rotamer( rsd.nchi(),  60, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 120, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 180, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 240, 10 );
	rsd.add_chi_rotamer( rsd.nchi(), 300, 10 );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithEthyl").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithEthyl::name() const {
	return "ReplaceProtonWithEthyl";
}

/// @brief replace proton with chlorine
bool
ReplaceProtonWithChlorine::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "CL"+atom_.substr( 1 ), "Cl", "CL", -0.10 );
	rsd.add_bond( "CL"+atom_.substr( 1 ), atom_ );
	rsd.set_icoor( "CL"+atom_.substr( 1 ),
		icoor->phi(),
		icoor->theta(),
		1.810005,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithChlorine").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithChlorine::name() const {
	return "ReplaceProtonWithChlorine";
}

/// @brief replace proton with fluorine
bool
ReplaceProtonWithFluorine::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "F"+atom_, "F", "F1", -0.22 );
	rsd.add_bond( "F"+atom_, atom_ );
	rsd.set_icoor( "F"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.332965,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithFluorine").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithFluorine::name() const {
	return "ReplaceProtonWithFluorine";
}

/// @brief replace proton with bromine
bool
ReplaceProtonWithBromine::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "BR"+atom_.substr( 1 ), "Br", "BR", -0.07 );
	rsd.add_bond( "BR"+atom_.substr( 1 ), atom_ );
	rsd.set_icoor( "BR"+atom_.substr( 1 ),
		icoor->phi(),
		icoor->theta(),
		1.903609,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithBromine").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithBromine::name() const {
	return "ReplaceProtonWithBromine";
}

/// @brief replace proton with iodine
bool
ReplaceProtonWithIodine::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "I"+atom_, "I", "I", -0.08 );
	rsd.add_bond( "I"+atom_, atom_ );
	rsd.set_icoor( "I"+atom_,
		icoor->phi(),
		icoor->theta(),
		2.142205,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithIodine").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithIodine::name() const {
	return "ReplaceProtonWithIodine";
}

/// @brief replace proton with hydroxyl
bool
ReplaceProtonWithHydroxyl::apply( MutableResidueType & rsd ) const {

	VD atm( rsd.atom_vertex( atom_ ) );

	utility::vector1< VD > const & bond_H( rsd.bonded_hydrogens( atm ) );
	if ( bond_H.empty() ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( atm ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	// //RM: I don't think this should be the responsibility of this patch.
	//for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
	// // Am I an aromatic carbon atom with a free valence?
	// // Can also be specified in params, but add if not already set!
	// rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
	// if ( rsd.atom( ii ).mm_name() == "CA" ) {
	//  AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
	//  for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
	//   if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
	//    rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
	//   }
	//  }
	// }
	//}

	VD proton_index = bond_H[1]; // Are we just assuming there's only the one?
	MutableICoorRecordCOP icoor = rsd.icoor( proton_index );
	debug_assert( icoor != nullptr );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "O"+atom_, "OH", "OH1", -0.54 );
	rsd.add_bond( "O"+atom_, atom_ );
	rsd.set_icoor( "O"+atom_,
		icoor->phi(),
		icoor->theta(),
		1.375964,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		icoor->stub_atom3(),
		true );

	rsd.add_atom( "H"+atom_, "Hpol", "H", 0.43 );
	rsd.add_bond( "H"+atom_, "O"+atom_ );

	rsd.set_icoor( "H"+atom_,
		0.000000*3.14159/180.000000,
		70.600000*3.14159/180.000000,
		0.960239,
		"O"+atom_,
		icoor->stub_atom1(),
		icoor->stub_atom2(),
		true );

	rsd.add_chi( icoor->stub_atom2(), icoor->stub_atom1(), "O"+atom_, "H"+atom_ );

	utility::vector1< Real > chi_samples;
	chi_samples.push_back( 0 );
	chi_samples.push_back( 60 );
	chi_samples.push_back( 120 );
	chi_samples.push_back( 180 );
	chi_samples.push_back( 240 );
	chi_samples.push_back( 300 );
	utility::vector1< Real > extra;
	chi_samples.push_back( 20 );
	rsd.set_proton_chi( rsd.nchi(), chi_samples, extra );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	//VD nbr_atom = rsd.nbr_vertex();
	//rsd.nbr_radius( find_nbr_dist( rsd, nbr_atom ) );
	//rsd.nbr_atom( nbr_atom );

	//rename_atoms( rsd, false );
	//rsd.assign_internal_coordinates();
	//rsd.autodetermine_chi_bonds();

	return false;
}

/// @brief Return the name of this PatchOperation ("ReplaceProtonWithHydroxyl").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
ReplaceProtonWithHydroxyl::name() const {
	return "ReplaceProtonWithHydroxyl";
}

bool
VirtualizeSidechain::apply( MutableResidueType & rsd ) const {
	// TODO: The original version of this also virtualized backbone hydrogens.
	// I'm currently keeping this behavior, but we probably shouldn't.
	for ( VD atm: rsd.all_atoms() ) {
		if ( rsd.is_backbone_heavyatom( atm ) ) {
			continue;
		}
		rsd.set_atom_type( atm, "VIRT" );
		rsd.atom( atm ).mm_name( "VIRT" );
		rsd.atom( atm ).charge( 0.0 );
		rsd.atom( atm ).is_virtual( true );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("VirtualizeSidechain").
/// @author Andrew M. Watkins (amw579@stanford.edu).
std::string
VirtualizeSidechain::name() const {
	return "VirtualizeSidechain";
}


bool
AddConnectAndTrackingVirt::apply( MutableResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!

	// Provide unique variant name
	std::string res_varname( "MP-" + atom_ + "-METAL_CONNECT" );
	Size count=0;
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for ResidueType: " + rsd.name() );
		}
		++count;
		if ( count == 1 ) {
			if ( !rsd.properties().has_custom_variant_types() || ! rsd.has_variant_type( res_varname ) ) break;
		} else {
			res_varname = "MP-" + atom_ + "-METAL_CONNECT" + std::to_string( count );
			if ( !rsd.properties().has_custom_variant_types() || ! rsd.has_variant_type( res_varname ) ) break;
		}
	}
	rsd.enable_custom_variant_types();
	rsd.add_variant_type( res_varname );


	Size const con_res = rsd.add_residue_connection( atom_ );
	tr.Trace << "Forming connection number " << con_res << " for " << rsd.name() << std::endl;
	Size virtcount = rsd.n_virtual_atoms();

	if ( con_res > rsd.n_polymeric_residue_connections() // no underflow!!
			&& virtcount < con_res-rsd.n_polymeric_residue_connections() ) {

		// Add a virt.  First, find a unique name for it:
		std::string virtname = "V" + ObjexxFCL::string_of( ++virtcount );
		while ( rsd.has( virtname )  ) {
			virtname = "V" + ObjexxFCL::string_of( ++virtcount );
		}

		tr.Trace << "Adding " << virtname << std::endl;
		//Create a new virt atom, using the unique name found above:
		rsd.add_atom( virtname, "VIRT", "VIRT", 0.0 );
		rsd.add_bond( virtname, atom_ );
		rsd.set_icoor( virtname, 1.0, 1.0, 1.37, rsd.atom_name( rsd.all_atoms()[1] ), ( "V" + ObjexxFCL::string_of( virtcount-2 ) ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), true );
	}
	//icoor = rsd.icoor( rsd.atom_index( virtname ) );

	// Okay, give CONN the icoor of the selected virt.
	if ( virtcount <= 3 ) { // like FE
		rsd.set_icoor( "CONN"+ObjexxFCL::string_of( con_res ), 0, 1.0, 1.37, rsd.atom_name( rsd.all_atoms()[1] ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), ( "V" + ObjexxFCL::string_of( virtcount ) ), true );
	} else {
		rsd.set_icoor( "CONN"+ObjexxFCL::string_of( con_res ), 10, 1.0, 1.37, rsd.atom_name( rsd.all_atoms()[1] ), ( "V" + ObjexxFCL::string_of( virtcount-2 ) ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), true );
	}
	return false;
}

bool
AddConnectAndTrackingVirt::changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const
{
	return rsd_type.has( atom ) && rsd_type.has( atom_ ) && rsd_type.atom_index( atom ) == rsd_type.atom_index( atom_ );
}

/// @brief Return the name of this PatchOperation ("AddConnectAndTrackingVirt").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddConnectAndTrackingVirt::name() const {
	return "AddConnectAndTrackingVirt";
}

bool
AddConnectDeleteChildProton::apply( MutableResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!
	if (  rsd.atom( atom_ ).is_hydrogen() ) return true; // Can't delete child proton from a proton
	rsd.add_metapatch_connect( atom_ );
	return false;
}

bool
AddConnectDeleteChildProton::changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const
{
	return rsd_type.has( atom_ ) && !rsd_type.atom_is_hydrogen( rsd_type.atom_index(atom_) ) && rsd_type.has( atom ) && rsd_type.atom_index( atom ) == rsd_type.atom_index( atom_ );
}

/// @brief Return the name of this PatchOperation ("AddConnectDeleteChildProton").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
AddConnectDeleteChildProton::name() const {
	return "AddConnectDeleteChildProton";
}

bool
DeleteChildProton::apply( MutableResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!
	rsd.delete_child_proton( atom_ );
	return false;
}

/// @brief Return the name of this PatchOperation ("DeleteChildProton").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
DeleteChildProton::name() const {
	return "DeleteChildProton";
}

bool
VirtualizeAll::apply( MutableResidueType & rsd ) const {
	for ( VD atm: rsd.all_atoms() ) {
		rsd.set_atom_type( atm, "VIRT" );
		rsd.atom( atm ).mm_name( "VIRT" );
		rsd.atom( atm ).charge( 0.0 );
		rsd.atom( atm ).is_virtual( true );
	}
	rsd.net_formal_charge( 0 );
	return false;
}

/// @brief Return the name of this PatchOperation ("VirtualizeAll").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
VirtualizeAll::name() const {
	return "VirtualizeAll";
}

// SetVirtualShadow //////////////////////////////////////////////////////////////

SetVirtualShadow::SetVirtualShadow(
	std::string const & shadower,
	std::string const & shadowee
):
	shadower_( shadower ),
	shadowee_( shadowee )
{}

bool
SetVirtualShadow::apply( MutableResidueType & rsd ) const
{
	if ( ! rsd.has( shadower_ ) ) {
		if ( TR_PatchOperations.Debug.visible() ) {
			TR_PatchOperations.Debug << "SetVirtualShadow::apply failed: " <<
				rsd.name() << " is missing atom: " << shadower_ << std::endl;
		}
		return true;  // failure
	} else  if ( ! rsd.has( shadowee_ ) ) {
		if ( TR_PatchOperations.Debug.visible() ) {
			TR_PatchOperations.Debug << "SetVirtualShadow::apply failed: " <<
				rsd.name() << " is missing atom: " << shadowee_ << std::endl;
		}
		return true;  // failure
	} else {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetVirtualShadow::apply: " <<
				shadower_ << ' ' << shadowee_ << std::endl;
		}
		rsd.set_shadowing_atom( shadower_, shadowee_ );
	}
	return false;
}

/// @brief Return the name of this PatchOperation ("SetVirtualShadow").
/// @author Vikram K. Mulligan (vmullig@uw.edu).
std::string
SetVirtualShadow::name() const {
	return "SetVirtualShadow";
}


// RenameAtom //////////////////////////////////////////////////////////////

bool
RenameAtom::apply( MutableResidueType & rsd ) const
{
	rsd.rename_atom( rsd.atom_vertex( old_name_ ), new_name_ );
	return false;
}

} // chemical
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ChangeBondType::ChangeBondType() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ChangeBondType::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( old_bond_type_ ) ); // std::string
	arc( CEREAL_NVP( new_bond_type_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ChangeBondType::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( old_bond_type_ ); // std::string
	arc( new_bond_type_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ChangeBondType );
CEREAL_REGISTER_TYPE( core::chemical::ChangeBondType )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetNbrRadius::SetNbrRadius() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetNbrRadius::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( radius_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetNbrRadius::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( radius_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetNbrRadius );
CEREAL_REGISTER_TYPE( core::chemical::SetNbrRadius )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddBond::AddBond() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddBond::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddBond::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddBond );
CEREAL_REGISTER_TYPE( core::chemical::AddBond )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetPolymerConnectAtom::SetPolymerConnectAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetPolymerConnectAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( upper_lower_ ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetPolymerConnectAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( upper_lower_ ); // int
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetPolymerConnectAtom );
CEREAL_REGISTER_TYPE( core::chemical::SetPolymerConnectAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ChangeAncestory::ChangeAncestory() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ChangeAncestory::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
	arc( CEREAL_NVP( which_ancestor_ ) ); // enum core::chemical::Ancestor
	arc( CEREAL_NVP( ancestor_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ChangeAncestory::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
	arc( which_ancestor_ ); // enum core::chemical::Ancestor
	arc( ancestor_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ChangeAncestory );
CEREAL_REGISTER_TYPE( core::chemical::ChangeAncestory )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithMethoxy::ReplaceProtonWithMethoxy() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithMethoxy::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithMethoxy::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithMethoxy );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithMethoxy )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetInterchangeabilityGroup_String::SetInterchangeabilityGroup_String() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetInterchangeabilityGroup_String::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( intgrp_ ) ); // const std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetInterchangeabilityGroup_String::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( intgrp_ ); // const std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetInterchangeabilityGroup_String );
CEREAL_REGISTER_TYPE( core::chemical::SetInterchangeabilityGroup_String )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::RedefineChi::RedefineChi() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::RedefineChi::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( chino_ ) ); // Size
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( atom3_ ) ); // std::string
	arc( CEREAL_NVP( atom4_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::RedefineChi::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( chino_ ); // Size
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( atom3_ ); // std::string
	arc( atom4_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::RedefineChi );
CEREAL_REGISTER_TYPE( core::chemical::RedefineChi )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetICoor::SetICoor() {}



/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetICoor::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
	arc( CEREAL_NVP( phi_ ) ); // Real
	arc( CEREAL_NVP( theta_ ) ); // Real
	arc( CEREAL_NVP( d_ ) ); // Real
	arc( CEREAL_NVP( stub1_ ) ); // std::string
	arc( CEREAL_NVP( stub2_ ) ); // std::string
	arc( CEREAL_NVP( stub3_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetICoor::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
	arc( phi_ ); // Real
	arc( theta_ ); // Real
	arc( d_ ); // Real
	arc( stub1_ ); // std::string
	arc( stub2_ ); // std::string
	arc( stub3_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetICoor );
CEREAL_REGISTER_TYPE( core::chemical::SetICoor )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithEthyl::ReplaceProtonWithEthyl() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithEthyl::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithEthyl::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithEthyl );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithEthyl )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetNbrAtom::SetNbrAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetNbrAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetNbrAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetNbrAtom );
CEREAL_REGISTER_TYPE( core::chemical::SetNbrAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AppendMainchainAtom::AppendMainchainAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AppendMainchainAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AppendMainchainAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AppendMainchainAtom );
CEREAL_REGISTER_TYPE( core::chemical::AppendMainchainAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithHydroxyl::ReplaceProtonWithHydroxyl() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithHydroxyl::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithHydroxyl::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithHydroxyl );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithHydroxyl )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetIO_String::SetIO_String() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetIO_String::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( name3_ ) ); // const std::string
	arc( CEREAL_NVP( name1_ ) ); // const char
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetIO_String::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( name3_ ); // const std::string
	arc( name1_ ); // const char
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetIO_String );
CEREAL_REGISTER_TYPE( core::chemical::SetIO_String )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddChi::AddChi() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddChi::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( no_index_ ) ); // _Bool
	arc( CEREAL_NVP( chino_ ) ); // Size
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( atom3_ ) ); // std::string
	arc( CEREAL_NVP( atom4_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddChi::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( no_index_ ); // _Bool
	arc( chino_ ); // Size
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( atom3_ ); // std::string
	arc( atom4_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddChi );
CEREAL_REGISTER_TYPE( core::chemical::AddChi )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetAtomType::SetAtomType() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetAtomType::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_type_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetAtomType::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( atom_type_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetAtomType );
CEREAL_REGISTER_TYPE( core::chemical::SetAtomType )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ChiralFlipAtoms::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ChiralFlipAtoms::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ChiralFlipAtoms );
CEREAL_REGISTER_TYPE( core::chemical::ChiralFlipAtoms )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithFluorine::ReplaceProtonWithFluorine() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithFluorine::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithFluorine::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithFluorine );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithFluorine )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetBackboneHeavyatom::SetBackboneHeavyatom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetBackboneHeavyatom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetBackboneHeavyatom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetBackboneHeavyatom );
CEREAL_REGISTER_TYPE( core::chemical::SetBackboneHeavyatom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetDisulfideAtomName::SetDisulfideAtomName() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetDisulfideAtomName::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetDisulfideAtomName::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetDisulfideAtomName );
CEREAL_REGISTER_TYPE( core::chemical::SetDisulfideAtomName )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetAtomicCharge::SetAtomicCharge() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetAtomicCharge::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( charge_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetAtomicCharge::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( charge_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetAtomicCharge );
CEREAL_REGISTER_TYPE( core::chemical::SetAtomicCharge )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::PatchOperation::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::PatchOperation::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::PatchOperation );
CEREAL_REGISTER_TYPE( core::chemical::PatchOperation )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetBaseName::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( new_base_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetBaseName::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( new_base_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetBaseName );
CEREAL_REGISTER_TYPE( core::chemical::SetBaseName )

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ChiralFlipNaming::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ChiralFlipNaming::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ChiralFlipNaming );
CEREAL_REGISTER_TYPE( core::chemical::ChiralFlipNaming )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddProperty::AddProperty() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddProperty::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( property_ ) ); // std::string
	arc( CEREAL_NVP( property_enum_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddProperty::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( property_ ); // std::string
	arc( property_enum_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddProperty );
CEREAL_REGISTER_TYPE( core::chemical::AddProperty )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteProperty::DeleteProperty() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteProperty::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( property_name_ ) );
	arc( CEREAL_NVP( property_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteProperty::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( property_name_ );
	arc( property_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteProperty );
CEREAL_REGISTER_TYPE( core::chemical::DeleteProperty )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithTrifluoromethyl::ReplaceProtonWithTrifluoromethyl() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithTrifluoromethyl::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithTrifluoromethyl::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithTrifluoromethyl );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithTrifluoromethyl )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddChiRotamer::AddChiRotamer() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddChiRotamer::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( no_index_ ) ); // _Bool
	arc( CEREAL_NVP( chino_ ) ); // Size
	arc( CEREAL_NVP( mean_ ) ); // Real
	arc( CEREAL_NVP( sdev_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddChiRotamer::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( no_index_ ); // _Bool
	arc( chino_ ); // Size
	arc( mean_ ); // Real
	arc( sdev_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddChiRotamer );
CEREAL_REGISTER_TYPE( core::chemical::AddChiRotamer )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddBondType::AddBondType() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddBondType::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( bond_type_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddBondType::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( bond_type_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddBondType );
CEREAL_REGISTER_TYPE( core::chemical::AddBondType )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteTerminalChi::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteTerminalChi::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteTerminalChi );
CEREAL_REGISTER_TYPE( core::chemical::DeleteTerminalChi )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddConnect::AddConnect() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddConnect::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( phi_ ) ); // core::Real
	arc( CEREAL_NVP( theta_ ) ); // core::Real
	arc( CEREAL_NVP( d_ ) ); // core::Real
	arc( CEREAL_NVP( connect_atom_ ) ); // const std::string
	arc( CEREAL_NVP( parent_atom_ ) ); // const std::string
	arc( CEREAL_NVP( angle_atom_ ) ); // const std::string
	arc( CEREAL_NVP( torsion_atom_ ) ); // const std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddConnect::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( phi_ ); // core::Real
	arc( theta_ ); // core::Real
	arc( d_ ); // core::Real
	arc( connect_atom_ ); // const std::string
	arc( parent_atom_ ); // const std::string
	arc( angle_atom_ ); // const std::string
	arc( torsion_atom_ ); // const std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddConnect );
CEREAL_REGISTER_TYPE( core::chemical::AddConnect )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddConnectAndTrackingVirt::AddConnectAndTrackingVirt() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddConnectAndTrackingVirt::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddConnectAndTrackingVirt::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddConnectAndTrackingVirt );
CEREAL_REGISTER_TYPE( core::chemical::AddConnectAndTrackingVirt )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetFormalCharge::SetFormalCharge() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetFormalCharge::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( charge_ ) ); // core::SSize
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetFormalCharge::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( charge_ ); // core::SSize
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetFormalCharge );
CEREAL_REGISTER_TYPE( core::chemical::SetFormalCharge )

/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetNetFormalCharge::SetNetFormalCharge() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetNetFormalCharge::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( charge_ ) ); // signed int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetNetFormalCharge::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( charge_ ); // signed int
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetNetFormalCharge );
CEREAL_REGISTER_TYPE( core::chemical::SetNetFormalCharge )

/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddConnectDeleteChildProton::AddConnectDeleteChildProton() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddConnectDeleteChildProton::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddConnectDeleteChildProton::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddConnectDeleteChildProton );
CEREAL_REGISTER_TYPE( core::chemical::AddConnectDeleteChildProton )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddProtonChi::AddProtonChi() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddProtonChi::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( chino_ ) ); // Size
	arc( CEREAL_NVP( samples_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( extrasamples_ ) ); // utility::vector1<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddProtonChi::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( chino_ ); // Size
	arc( samples_ ); // utility::vector1<core::Real>
	arc( extrasamples_ ); // utility::vector1<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddProtonChi );
CEREAL_REGISTER_TYPE( core::chemical::AddProtonChi )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ClearChiRotamers::ClearChiRotamers() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ClearChiRotamers::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( chi_no_ ) ); // core::uint
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ClearChiRotamers::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( chi_no_ ); // core::uint
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ClearChiRotamers );
CEREAL_REGISTER_TYPE( core::chemical::ClearChiRotamers )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteVariantType::DeleteVariantType() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteVariantType::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( variant_str_ ) );
	arc( CEREAL_NVP( variant_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteVariantType::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( variant_str_ );
	arc( variant_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteVariantType );
CEREAL_REGISTER_TYPE( core::chemical::DeleteVariantType )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetMMAtomType::SetMMAtomType() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetMMAtomType::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( mm_atom_type_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetMMAtomType::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( mm_atom_type_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetMMAtomType );
CEREAL_REGISTER_TYPE( core::chemical::SetMMAtomType )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithMethyl::ReplaceProtonWithMethyl() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithMethyl::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithMethyl::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithMethyl );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithMethyl )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithChlorine::ReplaceProtonWithChlorine() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithChlorine::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithChlorine::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithChlorine );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithChlorine )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteMetalbindingAtom::DeleteMetalbindingAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteMetalbindingAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteMetalbindingAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteMetalbindingAtom );
CEREAL_REGISTER_TYPE( core::chemical::DeleteMetalbindingAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteActCoordAtom::DeleteActCoordAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteActCoordAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteActCoordAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteActCoordAtom );
CEREAL_REGISTER_TYPE( core::chemical::DeleteActCoordAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::NCAARotLibPath::NCAARotLibPath() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::NCAARotLibPath::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( path_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::NCAARotLibPath::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( path_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::NCAARotLibPath );
CEREAL_REGISTER_TYPE( core::chemical::NCAARotLibPath )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddAtomAlias::AddAtomAlias() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddAtomAlias::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( rosetta_atom_name_ ) ); // std::string
	arc( CEREAL_NVP( aliases_ ) ); // utility::vector1< std::string >
	arc( CEREAL_NVP( preferred_alias_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddAtomAlias::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( rosetta_atom_name_ ); // std::string
	arc( aliases_ ); // utility::vector1< std::string >
	arc( preferred_alias_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddAtomAlias );
CEREAL_REGISTER_TYPE( core::chemical::AddAtomAlias )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ResetBondLength::ResetBondLength() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ResetBondLength::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atm_ ) ); // std::string
	arc( CEREAL_NVP( d_ ) ); // core::Distance
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResetBondLength::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atm_ ); // std::string
	arc( d_ ); // core::Distance
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResetBondLength );
CEREAL_REGISTER_TYPE( core::chemical::ResetBondLength )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetOrientAtom::SetOrientAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetOrientAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( force_nbr_atom_orient_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetOrientAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( force_nbr_atom_orient_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetOrientAtom );
CEREAL_REGISTER_TYPE( core::chemical::SetOrientAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteChildProton::DeleteChildProton() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteChildProton::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteChildProton::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteChildProton );
CEREAL_REGISTER_TYPE( core::chemical::DeleteChildProton )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetAllAtomsRepulsive::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetAllAtomsRepulsive::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetAllAtomsRepulsive );
CEREAL_REGISTER_TYPE( core::chemical::SetAllAtomsRepulsive )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ConnectSulfurAndMakeVirtualProton::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ConnectSulfurAndMakeVirtualProton::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ConnectSulfurAndMakeVirtualProton );
CEREAL_REGISTER_TYPE( core::chemical::ConnectSulfurAndMakeVirtualProton )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceMainchainAtom::ReplaceMainchainAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceMainchainAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( target_ ) ); // std::string
	arc( CEREAL_NVP( new_atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceMainchainAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( target_ ); // std::string
	arc( new_atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceMainchainAtom );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceMainchainAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::DeleteAtom::DeleteAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::DeleteAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::DeleteAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::DeleteAtom );
CEREAL_REGISTER_TYPE( core::chemical::DeleteAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithIodine::ReplaceProtonWithIodine() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithIodine::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithIodine::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithIodine );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithIodine )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::ReplaceProtonWithBromine::ReplaceProtonWithBromine() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithBromine::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ReplaceProtonWithBromine::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ReplaceProtonWithBromine );
CEREAL_REGISTER_TYPE( core::chemical::ReplaceProtonWithBromine )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::PrependMainchainAtom::PrependMainchainAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::PrependMainchainAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::PrependMainchainAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::PrependMainchainAtom );
CEREAL_REGISTER_TYPE( core::chemical::PrependMainchainAtom )


/// @brief Default constructor required by cereal to deserialize this class
core::chemical::AddAtom::AddAtom() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AddAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_type_name_ ) ); // std::string
	arc( CEREAL_NVP( mm_atom_type_name_ ) ); // std::string
	arc( CEREAL_NVP( charge_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AddAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( atom_name_ ); // std::string
	arc( atom_type_name_ ); // std::string
	arc( mm_atom_type_name_ ); // std::string
	arc( charge_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AddAtom );
CEREAL_REGISTER_TYPE( core::chemical::AddAtom )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_PatchOperation )

/// @brief Default constructor required by cereal to deserialize this class
core::chemical::SetVirtualShadow::SetVirtualShadow() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::SetVirtualShadow::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( CEREAL_NVP( shadower_ ) ); // std::string
	arc( CEREAL_NVP( shadowee_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::SetVirtualShadow::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::PatchOperation >( this ) );
	arc( shadower_ ); // std::string
	arc( shadowee_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::SetVirtualShadow );
CEREAL_REGISTER_TYPE( core::chemical::SetVirtualShadow )


#endif // SERIALIZATION
