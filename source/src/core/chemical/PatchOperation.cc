// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/PatchOperation.cc
/// @brief  Polymorphic classes representing the contents of a residue-type patch file
/// @author Phil Bradley

// Unit headers
#include <core/chemical/PatchOperation.hh>

// Package headers
#include <core/chemical/ResidueType.hh>

// Project headers
#include <core/chemical/Atom.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/Bond.hh>

// Numeric headers
#include <numeric/conversions.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>


namespace core {
namespace chemical {

/// @details Auto-generated virtual destructor
PatchOperation::~PatchOperation() {}

static THREAD_LOCAL basic::Tracer tr( "core.chemical" );
static THREAD_LOCAL basic::Tracer TR_PatchOperations( "core.chemical.PatchOperations.hh" );

DeleteAtom::DeleteAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
DeleteAtom::apply( ResidueType & rsd ) const
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

SetBackboneHeavyatom::SetBackboneHeavyatom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
SetBackboneHeavyatom::apply( ResidueType & rsd ) const
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
SetPolymerConnectAtom::apply( ResidueType & rsd ) const
{
	if ( atom_name_ == "NONE" || rsd.has( atom_name_ ) ) {
		if ( TR_PatchOperations.Trace.visible() ) {
			TR_PatchOperations.Trace << "SetPolymerConnectAtom::apply: " <<
				atom_name_ << ' ' << upper_lower_ << std::endl;
		}
		if ( upper_lower_ == -1 ) {
			rsd.set_lower_connect_atom( atom_name_ );
		} else {
			debug_assert( upper_lower_ == 1 );
			rsd.set_upper_connect_atom( atom_name_ );
		}
	} else {
		TR_PatchOperations.Debug << "SetPolymerConnectAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	}
	return false;
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
AddConnect::apply( ResidueType & rsd ) const
{
	if ( !rsd.has( connect_atom_ ) ||
			!rsd.has(  parent_atom_ ) ||
			!rsd.has(   angle_atom_ ) ||
			!rsd.has( torsion_atom_ ) ) return true; // failure!

	Size const connid( rsd.add_residue_connection( connect_atom_ ) );
	rsd.set_icoor( "CONN"+ObjexxFCL::string_of( connid ), phi_, theta_, d_, parent_atom_, angle_atom_, torsion_atom_ );
	return false;
}

AddProperty::AddProperty( std::string const & property_in ):
	property_( property_in )
{}

bool
AddProperty::apply( ResidueType & rsd ) const
{
	rsd.add_property( property_ );
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "AddProperty::apply: " << property_ << std::endl;
	}
	return false;
}

DeleteProperty::DeleteProperty( std::string const & property_in ):
	property_( property_in )
{}

bool
DeleteProperty::apply( ResidueType & rsd ) const
{
	rsd.delete_property( property_ );
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "DeleteProperty::apply: " << property_ << std::endl;
	}
	return false;
}


// DeleteVariantType //////////////////////////////////////////////////////////
DeleteVariantType::DeleteVariantType( std::string const & variant_in ) :
	variant_( variant_in )
{}

/// @return  false on success
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
DeleteVariantType::apply( ResidueType & rsd ) const
{
	rsd.remove_variant_type( variant_ );
	return false;  // success
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
AddChi::apply( ResidueType & rsd ) const
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


// AddProtonChi //////////////////////////////////////////////////////////////
AddProtonChi::AddProtonChi(
	Size const & chino_in, utility::vector1<core::Real> const & samples, utility::vector1<core::Real> const & extrasamples
):
	chino_( chino_in ), samples_(samples), extrasamples_(extrasamples)
{}

bool
AddProtonChi::apply( ResidueType & rsd ) const
{
	rsd.set_proton_chi( chino_, samples_, extrasamples_ );
	return false;
}


RedefineChi::RedefineChi(Size const & chino_in,
	std::string const & atom1_in,
	std::string const & atom2_in,
	std::string const & atom3_in,
	std::string const & atom4_in):
	chino_( chino_in ), atom1_( atom1_in ), atom2_( atom2_in ), atom3_( atom3_in ), atom4_( atom4_in )
{}

bool
RedefineChi::apply( ResidueType & rsd ) const
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

bool
DeleteTerminalChi::apply( ResidueType & rsd ) const
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

DeleteMetalbindingAtom::DeleteMetalbindingAtom(
	std::string const & atom_name
):
	atom_name_( atom_name )
{}

bool
DeleteMetalbindingAtom::apply( ResidueType & rsd ) const
{
	rsd.delete_metalbinding_atom( atom_name_ );
	return false;
}

DeleteActCoordAtom::DeleteActCoordAtom(
	std::string const & atom_name
):
	atom_name_( atom_name )
{}

bool
DeleteActCoordAtom::apply( ResidueType & rsd ) const
{
	rsd.delete_act_coord_atom( atom_name_ );
	return false;
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
AddChiRotamer::apply( ResidueType & rsd ) const
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


// ClearChiRotamers ///////////////////////////////////////////////////////////

ClearChiRotamers::ClearChiRotamers( core::uint const chi_no_in ) :
	chi_no_( chi_no_in )
{}

/// @return  true on failure
bool
ClearChiRotamers::apply( ResidueType & rsd ) const
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


// AddAtom ///////////////////////////////////////////////////////////////////
#if defined(WIN32) && !defined(WIN_PYROSETTA)
AddAtomWIN32::AddAtomWIN32(
#else
AddAtom::AddAtom(
#endif
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

bool
#if defined(WIN32) && !defined(WIN_PYROSETTA)
AddAtomWIN32::apply( ResidueType & rsd ) const
#else
AddAtom::apply( ResidueType & rsd ) const
#endif
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

AddAtomAlias::AddAtomAlias( std::string const & rosetta_atom_name_in, std::string const & alias_in ) :
	rosetta_atom_name_( rosetta_atom_name_in ),
	alias_( alias_in )
{}

/// @return  true on failure
/// @remarks Because ResidueType is not throwing exceptions, this will never return true.  Failure will lead to exits
/// from ResidueType. ~Labonte
bool
AddAtomAlias::apply( ResidueType & rsd ) const
{
	rsd.add_atom_alias( rosetta_atom_name_, alias_ );
	return false;  // success
}


// AddBond ////////////////////////////////////////////////////////////////////

AddBond::AddBond(
	std::string const & atom1_in,
	std::string const & atom2_in
):
	atom1_( atom1_in ),
	atom2_( atom2_in )
{}

bool
AddBond::apply( ResidueType & rsd ) const
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
AddBondType::apply( ResidueType & rsd ) const
{
	rsd.add_bond( atom1_, atom2_, convert_to_BondName( bond_type_ ) );
	return false;  // success
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
ChangeBondType::apply( ResidueType & rsd ) const
{
	rsd.change_bond_type(
			atom1_, atom2_, convert_to_BondName( old_bond_type_ ), convert_to_BondName( new_bond_type_ ) );
	return false;  // success
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
SetAtomicCharge::apply( ResidueType & rsd ) const
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


// SetFormalCharge ////////////////////////////////////////////////////////////

SetFormalCharge::SetFormalCharge( std::string const & atom_name_in, core::SSize charge_in ) :
	atom_name_( atom_name_in ),
	charge_( charge_in )
{}

/// @return  true on failure
bool
SetFormalCharge::apply( ResidueType & rsd ) const
{
	if ( ! rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Error << "SetFormalCharge::apply() failed: " << rsd.name() << " is missing atom: " <<
			atom_name_ << std::endl;
		return true;  // failure
	}
	rsd.atom( atom_name_ ).formal_charge( charge_ );
	return false;  // success
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
SetAtomType::apply( ResidueType & rsd ) const
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


// SetIO_String //////////////////////////////////////////////////////////////

SetIO_String::SetIO_String(
	std::string const & name3,
	char const name1
):
	name3_( name3 ),
	name1_( name1 )
{}

bool
SetIO_String::apply( ResidueType & rsd ) const
{
	rsd.name3( name3_ );
	rsd.name1( name1_ );
	return false;
}

SetInterchangeabilityGroup_String::SetInterchangeabilityGroup_String(
	std::string const & intgrp
):
	intgrp_( intgrp )
{}

bool
SetInterchangeabilityGroup_String::apply( ResidueType & rsd ) const
{
	rsd.interchangeability_group( intgrp_ );
	return false;
}

AppendInterchangeabilityGroup_String::AppendInterchangeabilityGroup_String(
	std::string const & intgrp_addendum
):
	intgrp_addendum_( intgrp_addendum )
{}

bool
AppendInterchangeabilityGroup_String::apply( ResidueType & rsd ) const
{
	std::string new_intgrp = rsd.interchangeability_group() + intgrp_addendum_;
	rsd.interchangeability_group( new_intgrp );
	return false;
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
SetMMAtomType::apply( ResidueType & rsd ) const
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
		rsd.set_mm_atom_type( atom_name_, mm_atom_type_name_ );
	}
	return false;
}


// SetICoor //////////////////////////////////////////////////////////////////

// helper function
std::string
expand_icoor_atom_name( std::string name, ResidueType const & rsd )
{
	std::string const nconn_tag( "%NCONN" );
	Size pos( name.find( nconn_tag ) );
	if ( pos < name.size() ) {
		//std::cout << "name before replace: " << name << std::endl;
		name.replace( pos, nconn_tag.size(), ObjexxFCL::string_of( rsd.n_residue_connections() ) );
		//std::cout << "name after replace: " << name << std::endl;
	}
	return name;
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

bool
SetICoor::apply( ResidueType & rsd ) const
{
	if ( TR_PatchOperations.Trace.visible() ) {
		TR_PatchOperations.Trace << "SetICoor::apply: " <<
			atom_ << ' ' << stub1_ << ' ' << stub2_ << ' ' << stub3_ << std::endl;
	}
	std::string const atom ( expand_icoor_atom_name( atom_ , rsd ) );
	std::string const stub1( expand_icoor_atom_name( stub1_, rsd ) );
	std::string const stub2( expand_icoor_atom_name( stub2_, rsd ) );
	std::string const stub3( expand_icoor_atom_name( stub3_, rsd ) );
	//  bool const rebuild_icoor_xyz( ICoorAtomID( stub1, rsd ).is_internal() &&
	//                 ICoorAtomID( stub2, rsd ).is_internal() &&
	//                 ICoorAtomID( stub3, rsd ).is_internal() );
	bool const rebuild_icoor_xyz( true );
	rsd.set_icoor( atom, phi_, theta_, d_, stub1, stub2, stub3, rebuild_icoor_xyz );
	return false;
}


// ResetBondLength ////////////////////////////////////////////////////////////

ResetBondLength::ResetBondLength( std::string const & atm_in, core::Distance d_in ) :
	atm_( atm_in ),
	d_( d_in )
{}

/// @return  true on failure
bool
ResetBondLength::apply( ResidueType & rsd ) const
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


// PrependMainchainAtom ///////////////////////////////////////////////////////

PrependMainchainAtom::PrependMainchainAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
PrependMainchainAtom::apply( ResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "PrependMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		AtomIndices const & old_mainchain_atoms( rsd.mainchain_atoms() );
		AtomIndices new_mainchain_atoms;
		new_mainchain_atoms.push_back( rsd.atom_index( atom_name_ ) );
		for ( Size i = 1; i <= old_mainchain_atoms.size(); ++i ) {
			new_mainchain_atoms.push_back( old_mainchain_atoms[i] );
		}
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}

AppendMainchainAtom::AppendMainchainAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
AppendMainchainAtom::apply( ResidueType & rsd ) const
{
	if ( !rsd.has( atom_name_ ) ) {
		TR_PatchOperations.Debug << "AppendMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << atom_name_ << std::endl;
		return true; // failure
	} else {
		AtomIndices new_mainchain_atoms( rsd.mainchain_atoms() );
		new_mainchain_atoms.push_back( rsd.atom_index( atom_name_ ) );
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}

ReplaceMainchainAtom::ReplaceMainchainAtom( std::string const & target, std::string const & new_atom ) :
	target_( target ),
	new_atom_( new_atom )
{}

bool
ReplaceMainchainAtom::apply( ResidueType & rsd ) const
{
	if ( !rsd.has( target_ ) ) {
		TR_PatchOperations.Debug << "ReplaceMainchainAtom::apply failed: " <<
			rsd.name() << " is missing atom " << target_ << std::endl;
		return true; // failure
	} else {
		AtomIndices new_mainchain_atoms( rsd.mainchain_atoms() );
		for ( Size i = 1; i <= new_mainchain_atoms.size(); ++i ) {
			std::string mainchain_atom_name = rsd.atom_name( new_mainchain_atoms[ i ] );
			mainchain_atom_name.erase( std::remove_if( mainchain_atom_name.begin(),
				mainchain_atom_name.end(),
				::isspace ),
				mainchain_atom_name.end() );

			if ( mainchain_atom_name == target_ ) {
				new_mainchain_atoms[ i ] = rsd.atom_index( new_atom_ );
				break;
			}
		}
		rsd.set_mainchain_atoms( new_mainchain_atoms );
	}
	return false;
}


SetNbrAtom::SetNbrAtom( std::string const & atom_name_in ) :
	atom_name_( atom_name_in )
{}

bool
SetNbrAtom::apply( ResidueType & rsd ) const
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

SetNbrRadius::SetNbrRadius( Real const & radius ) :
	radius_( radius )
{}

bool
SetNbrRadius::apply( ResidueType & rsd ) const
{
	rsd.nbr_radius( radius_ );
	return false;
}


SetOrientAtom::SetOrientAtom(bool force_nbr_atom_orient):
	force_nbr_atom_orient_(force_nbr_atom_orient)
{}

bool
SetOrientAtom::apply( ResidueType & rsd ) const
{
	rsd.force_nbr_atom_orient( force_nbr_atom_orient_ );
	return false;
}

NCAARotLibPath::NCAARotLibPath( std::string const & path_in ) :
	path_( path_in )
{}

/// @brief set the NCAA rotamer library path in the residue type
bool
NCAARotLibPath::apply( ResidueType & rsd ) const
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
		ncaa_libspec = NCAARotamerLibrarySpecificationOP( new NCAARotamerLibrarySpecification( *old_libspec ) );
	} else {
		ncaa_libspec = NCAARotamerLibrarySpecificationOP(
			new core::chemical::rotamers::NCAARotamerLibrarySpecification );
	}
	ncaa_libspec->ncaa_rotlib_path( path_ );
	rsd.rotamer_library_specification( ncaa_libspec );
	return false;
}

/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
bool
ConnectSulfurAndMakeVirtualProton::apply( ResidueType & rsd ) const {
	std::string disulfide_atom_name = rsd.get_disulfide_atom_name();
	std::string CB_equivalent = rsd.atom_name( rsd.atom_base( rsd.atom_index( disulfide_atom_name ) ) );
	std::string CA_equivalent = rsd.atom_name( rsd.atom_base( rsd.atom_index( CB_equivalent ) ) );
	AddConnect ad(
		disulfide_atom_name,
		180, 68.374, 1.439,
		disulfide_atom_name,
		CB_equivalent,
		CA_equivalent
	);
	bool x = ad.apply( rsd );

	// Now I need to grab the proton sitting on the disulfide atom
	std::string HG_name = rsd.atom_name( rsd.attached_H_begin( rsd.atom_index( disulfide_atom_name ) ) );
	rsd.set_atom_type( HG_name, "VIRT" );
	rsd.set_mm_atom_type( HG_name, "VIRT" );
	rsd.atom( HG_name ).charge( 0.0 );
	rsd.atom( HG_name ).is_virtual( true );

	return x;
}
/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
bool
SetAllAtomsRepulsive::apply( ResidueType & rsd ) const {

	for ( Size i = 1; i <= rsd.natoms(); ++i ) {
		if ( rsd.atom_is_hydrogen( i ) ) {
			rsd.set_atom_type( rsd.atom_name( i ), "REPLS" );
		} else {
			rsd.set_atom_type( rsd.atom_name( i ), "HREPS" );
		}
	}

	return false;
}

bool
ChiralFlipNaming::apply( ResidueType & rsd ) const {

	if ( rsd.aa() <= aa_tyr ) {
		rsd.aa( get_D_equivalent( rsd.aa() ) );

		if ( rsd.aa() == aa_dal ) {
			rsd.name3( "DAL" );
			rsd.interchangeability_group( "DAL" );
			//std::cout << "igroup is " << rsd.interchangeability_group() << std::endl;
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
			rsd.name3( "DPR" );
			rsd.interchangeability_group( "DPR" );
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

	} else {
		// e.g. DHLU, DC01--will be caught later and set to DAL etc. for canonicals
		rsd.interchangeability_group( "D" + rsd.name().substr( 0, rsd.name().find(":") ) );
	}
	//tr << rsd.aa() << std::endl;
	return false;
}


bool
ChiralFlipAtoms::apply( ResidueType & rsd ) const {

	// So I think I need to separately process shadow atoms.
	// The reason is because they're 0 from something
	// so I think maybe I should delete and re-add them to a residue?

	//tr << "In chiralflipatoms apply " << std::endl;

	// Flip each
	//rsd.debug_dump_icoor();
	rsd.delete_property( "L_AA" );
	rsd.add_property( "D_AA" );

	if ( rsd.natoms() < 4 ) return true;
	for ( core::Size ii = 2; ii <= rsd.natoms(); ++ii ) {
		//pretty_print_atomicoor(std::cout, rsd.icoor( ii ), rsd );

		AtomICoor icoor1 = rsd.icoor( ii );

		std::string n1 = icoor1.stub_atom1().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom1().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom1().atomno() == 0 ? " N  " : rsd.atom( icoor1.stub_atom1().atomno() ).name() ) );
		std::string n2 = icoor1.stub_atom2().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom2().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom2().atomno() == 0 ? " N  " : rsd.atom( icoor1.stub_atom2().atomno() ).name() ) );
		std::string n3 = icoor1.stub_atom3().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom3().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom3().atomno() == 0 ? " N  " : rsd.atom( icoor1.stub_atom3().atomno() ).name() ) );

		rsd.set_icoor( rsd.atom( ii ).name(),
			-1.0*icoor1.phi(),
			icoor1.theta(),
			icoor1.d(),
			n1, n2, n3 );

		//pretty_print_atomicoor(std::cout, rsd.icoor( ii ), rsd );
	}

	rsd.finalize();
	rsd.fill_ideal_xyz_from_icoor();

	rsd.finalize();
	//rsd.debug_dump_icoor();
	// AMW: The unpatched type does not typically set these because they are unneeded
	// Since this is now a D-aa and it needs to hijack other rotamers, it needs them
	// (Maybe it would use them without; who knows! Don't think so, though.)
	//is backbone aa ness handled elsewhere?
	//rsd.backbone_aa( rsd.name() );

	return false;
}

PatchOperationOP
patch_operation_from_patch_file_line( std::string const & line ) {
	using numeric::conversions::radians;
	std::istringstream l( line );
	std::string tag, atom1, atom2, atom3, atom4, atom_name, atom_alias, atom_type_name, mm_atom_type_name, bond_type,
		property, variant, dummy;
	Real charge, mean, sdev, radius;
	Size chino;
	SSize formal_charge;
	l >> tag;
	if ( l.fail() || tag[0] == '#' ) return 0;
	if ( tag == "ADD_ATOM" ) {
		if ( line.size() < 25 ) return 0;
		atom_name = line.substr( 9,4); l >> tag;
		l >> atom_type_name; // = line.substr( 14,4);
		l >> mm_atom_type_name; // = line.substr( 19,4);
		l >> charge;
		if ( l.fail() ) return 0;
#if defined(WIN32) && !defined(WIN_PYROSETTA)
		return PatchOperationOP( new AddAtomWIN32( atom_name, atom_type_name, mm_atom_type_name, charge ) );
#else
		return PatchOperationOP( new AddAtom( atom_name, atom_type_name, mm_atom_type_name, charge ) );
#endif
	} else if ( tag == "DELETE_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new DeleteAtom( atom_name ) );

	} else if ( tag == "ADD_ATOM_ALIAS" ) {
		l >> atom_name >> atom_alias;
		if ( l.fail() ) {
			return NULL;
		}
		return PatchOperationOP( new AddAtomAlias( atom_name, atom_alias ) );

	} else if ( tag == "SET_BACKBONE_HEAVYATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new SetBackboneHeavyatom( atom_name ) );

	} else if ( tag == "SET_IO_STRING" ) { // 13 character tag
		// NOTE - USE FIXED WIDTH IO SINCE NAME3 CAN CONTAIN INTERNAL WHITESPACE (EG DNA,RNA)
		if ( line.size() < 19 ) return 0;
		std::string const three_letter_code( line.substr(14,3) ), one_letter_code( line.substr(18,1) );
		return PatchOperationOP( new SetIO_String( three_letter_code, one_letter_code[0] ) );

	} else if ( tag == "SET_INTERCHANGEABILITY_GROUP" ) {
		std::string intgrp;
		l >> intgrp;
		return PatchOperationOP( new SetInterchangeabilityGroup_String( intgrp ) );

	} else if ( tag == "APPEND_INTERCHANGEABILITY_GROUP" ) {
		std::string intgrp_addendum;
		l >> intgrp_addendum;
		return PatchOperationOP( new AppendInterchangeabilityGroup_String( intgrp_addendum ) );

	} else if ( tag == "NBR_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new SetNbrAtom( atom_name ) );

	} else if ( tag == "NBR_RADIUS" ) {
		Real radius;
		l >> radius;
		return PatchOperationOP( new SetNbrRadius( radius ) );

	} else if ( tag == "ADD_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new AddProperty( property ) );

	} else if ( tag == "DELETE_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new DeleteProperty( property ) );

	} else if ( tag == "DELETE_VARIANT_TYPE" ) {
		l >> variant;
		if ( l.fail() ) { return 0; }
		return PatchOperationOP( new DeleteVariantType( variant ) );

		// Added by Andy M. Chen in June 2009
		// This is needed for adding new side-chain torsions, which occurs in certain PTMs.
	} else if ( tag == "ADD_CHI" ) {
		if ( line.substr(8, 3) == "N+1" ) {
			l >> dummy >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) return 0;
			return PatchOperationOP( new AddChi(atom1, atom2, atom3, atom4) );
		} else {
			l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) return 0;
			return PatchOperationOP( new AddChi(chino, atom1, atom2, atom3, atom4) );
		}

		// Added by Andy M. Chen in June 2009
		// This is needed for PTMs, which often result in one or more extra chi angles.
	} else if ( tag == "ADD_PROTON_CHI" ) {
		Size chino, nsamples, nextra_samples;
		std::string dummy;
		l >> chino;
		l >> dummy; // should be "SAMPLES"
		l >> nsamples;
		utility::vector1< Real > samples( nsamples );
		for ( Size ii = 1; ii <= nsamples; ++ii ) {
			l >> samples[ ii ];
		}
		l >> dummy; // should be "EXTRA"
		l >> nextra_samples;
		utility::vector1< Real > extra_samples( nextra_samples );
		for ( Size ii = 1; ii <= nextra_samples; ++ii ) {
			l >> extra_samples[ ii ];
		}
		if ( l.fail() ) return 0;
		return PatchOperationOP( new AddProtonChi( chino, samples, extra_samples ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "REDEFINE_CHI" ) {
		l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new RedefineChi( chino, atom1, atom2, atom3, atom4 ) );

	} else if ( tag == "DELETE_TERMINAL_CHI" ) {
		return PatchOperationOP( new DeleteTerminalChi() );
	} else if ( tag == "DELETE_METALBINDING_ATOM" ) {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) return 0;

		return PatchOperationOP( new DeleteMetalbindingAtom( atom_name ) );
	} else if ( tag == "DELETE_ACT_COORD_ATOM" ) {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) return 0;

		return PatchOperationOP( new DeleteActCoordAtom( atom_name ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_CHI_ROTAMER" ) {
		if ( line.substr(16, 1) == "N" ) {
			l >> dummy >> mean >> sdev;
			if ( l.fail() ) return 0;
			return PatchOperationOP( new AddChiRotamer(mean, sdev) );
		} else {
			l >> chino >> mean >> sdev;
			if ( l.fail() ) return 0;
			return PatchOperationOP( new AddChiRotamer(chino, mean, sdev) );
		}

	} else if ( tag == "CLEAR_CHI_ROTAMERS" ) {
		l >> chino;
		if ( l.fail() ) { return 0; }
		return PatchOperationOP( new ClearChiRotamers( chino ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_BOND" ) {
		l >> atom1 >> atom2;
		if ( l.fail() ) return 0;
		return PatchOperationOP( new AddBond( atom1, atom2 ) );

	} else if ( tag == "ADD_BOND_TYPE" ) {
		l >> atom1 >> atom2 >> bond_type;
		if ( l.fail() ) {
			return NULL;
		}
		return PatchOperationOP( new AddBondType( atom1, atom2, bond_type ) );

	} else if ( tag == "CHANGE_BOND_TYPE" ) {
		std::string old_bond_type;
		l >> atom1 >> atom2 >> old_bond_type >> dummy >> bond_type;
		if ( l.fail() ) {
			return 0;
		}
		return PatchOperationOP( new ChangeBondType( atom1, atom2, old_bond_type, bond_type ) );
		
	} else if ( tag == "ADD_CONNECT" ) {
		std::string connect_atom;
		l >> connect_atom;
		if ( l.fail() ) return 0;
		l >> tag;
		if ( l.fail() ) {
			return PatchOperationOP( new AddConnect(
				connect_atom, 0.0, 0.0, 0.0, connect_atom, connect_atom, connect_atom ) );
		} else {
			Real phi, theta, d;
			std::string parent_atom, angle_atom, torsion_atom;
			l >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;
			if ( l.fail() || tag != "ICOOR" ) {
				utility_exit_with_message( "bad line in patchfile: "+line );
			}
			return PatchOperationOP( new AddConnect(
				connect_atom, radians(phi), radians(theta), d, parent_atom, angle_atom, torsion_atom ) );
		}

	} else if ( tag == "SET_ATOM_TYPE" ) {
		l >> atom_name >> atom_type_name;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetAtomType( atom_name, atom_type_name ) );

	} else if ( tag == "SET_MM_ATOM_TYPE" ) {
		runtime_assert( l.good() );
		l >> atom_name >> mm_atom_type_name;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetMMAtomType( atom_name, mm_atom_type_name ) );

	} else if ( tag == "SET_FORMAL_CHARGE" ) {
		l >> atom_name >> formal_charge;
		if ( l.fail() ) {
			return NULL;
		}
		return PatchOperationOP( new SetFormalCharge( atom_name, formal_charge ) );

	} else if ( tag == "SET_ATOMIC_CHARGE" ) {
		l >> atom_name >> charge;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetAtomicCharge( atom_name, charge ) );

	} else if ( tag == "SET_POLYMER_CONNECT" ) {
		l >> tag >> atom_name; // tag should be "UPPER" or "LOWER"
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetPolymerConnectAtom( atom_name, tag ) );

	} else if ( tag == "SET_ICOOR" ) {
		Real phi,theta,d;
		std::string stub1, stub2, stub3;
		l >> atom_name >> phi >> theta >> d >> stub1 >> stub2 >> stub3;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetICoor( atom_name, radians(phi), radians(theta), d, stub1, stub2, stub3 ) );

	} else if ( tag == "RESET_BOND_LENGTH" ) {
		core::Distance d;
		l >> atom_name >> d;
		if ( l.fail() ) {
			utility_exit_with_message( "bad line in patchfile: " + line );
		}
		return PatchOperationOP( new ResetBondLength( atom_name, d ) );

	} else if ( tag == "PREPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new PrependMainchainAtom( atom_name ) );

	} else if ( tag == "APPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new AppendMainchainAtom( atom_name ) );

	}  else if ( tag == "REPLACE_MAINCHAIN_ATOM" ) {
		std::string target, new_atom;
		l >> target >> new_atom;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new ReplaceMainchainAtom( target, new_atom ) );

	} else if ( tag == "NCAA_ROTLIB_PATH" ) {
		std::string path;
		l >> path;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new NCAARotLibPath( path ) );

	} else if ( tag == "SET_NBR_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetNbrAtom( atom_name ) );

	} else if ( tag == "SET_NBR_RADIUS" ) {
		l >> radius;
		if ( l.fail() ) utility_exit_with_message( line );
		return PatchOperationOP( new SetNbrRadius( radius ) );

	} else if ( tag == "SET_ORIENT_ATOM" ) {
		l >> tag;
		if ( l.fail() ) utility_exit_with_message( line );
		if ( tag == "NBR" ) {
			return PatchOperationOP( new SetOrientAtom(true) );
		} else if ( tag == "DEFAULT" ) {
			return PatchOperationOP( new SetOrientAtom(false) );
		} else {
			tr.Warning << "Unknown SET_ORIENT ATOM tag: " << tag << std::endl;
			return 0;
		}
	} else if ( tag == "SET_ALL_ATOMS_REPULSIVE" ) {
		return PatchOperationOP( new SetAllAtomsRepulsive() );
	} else if ( tag == "CONNECT_SULFUR_AND_MAKE_VIRTUAL_PROTON" ) {
		return PatchOperationOP( new ConnectSulfurAndMakeVirtualProton() );
	} else if ( tag == "CHIRAL_FLIP_NAMING" ) {
		return PatchOperationOP( new ChiralFlipNaming );//( atom1, atom2 ) );
	} else if ( tag == "CHIRAL_FLIP_ATOMS" ) {
		//std::string atom1, atom2;
		//l >> atom1 >> atom2;
		return PatchOperationOP( new ChiralFlipAtoms );//( atom1, atom2 ) );
	}
	tr.Warning << "patch_operation_from_patch_file_line: bad line: " << line << std::endl;

	return 0;
}

} // chemical
} // core
