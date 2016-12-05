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

// Project headers
#include <core/chemical/Atom.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/bond_support.hh>
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
PatchOperation::~PatchOperation() = default;

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
		name.replace( pos, nconn_tag.size(), ObjexxFCL::string_of( rsd.n_possible_residue_connections() ) );
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
ChangeAncestory::apply( ResidueType & rsd ) const
{
	Size const atind( rsd.atom_index( atom_ ));
	AtomICoor const aticoor( rsd.icoor( atind ));

	ICoorAtomID pa, gp, gg;

	try {
		pa = which_ancestor_ == anc_parent           ? ICoorAtomID( ancestor_name_, rsd ) : aticoor.stub_atom1();
		gp = which_ancestor_ == anc_grandparent      ? ICoorAtomID( ancestor_name_, rsd ) : aticoor.stub_atom2();
		gg = which_ancestor_ == anc_greatgrandparent ? ICoorAtomID( ancestor_name_, rsd ) : aticoor.stub_atom3();
	} catch ( utility::excn::EXCN_Base const & excn ) {
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

	rsd.set_icoor( atom_, aticoor.phi(), aticoor.theta(), aticoor.d(), pa, gp, gg, true /*rebuild_xyz*/ );
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
	rsd.finalize();
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
	rsd.finalize();
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


/// @brief Constructor.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
RemoveRotamerSpecifications::RemoveRotamerSpecifications() {}

/// @brief Strip all RotamerSpecifications from the ResidueType.
/// @return Success ("false" -- i.e. no error thrown) or failure ("true").
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
RemoveRotamerSpecifications::apply( ResidueType & rsd ) const {
	using namespace core::chemical::rotamers;

	if ( rsd.rotamer_library_specification() ) {
		rsd.strip_rotamer_library_specification();
	}
	return false;
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
	ResidueType & rsd
) const {
	rsd.set_rama_prepro_map_file_name( non_prepro_file_, false);
	rsd.set_rama_prepro_map_file_name( prepro_file_, true);
	return false;
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
	ResidueType & rsd
) const {
	rsd.set_rama_prepro_mainchain_torsion_potential_name( resname_, false );
	rsd.set_rama_prepro_mainchain_torsion_potential_name( resname_, true );
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
	ResidueType & rsd
) const {
	using namespace core::chemical::rotamers;

	NCAARotamerLibrarySpecificationOP ncaa_libspec( utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification >( rsd.rotamer_library_specification_nonconst() ) );
	runtime_assert_string_msg( ncaa_libspec, "A patch may only set the number of NCAA rotamer bins (\"NCAA_ROTLIB_NUM_ROTAMER_BINS\") for a residue type that already has an NCAA rotamer library specified, either in the params file or in patches." );
	ncaa_libspec->ncaa_rotlib_n_bin_per_rot( binsizes_ );
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

	} else if ( rsd.aa() == na_rad ) {
		rsd.name3( " 0A" );
		rsd.interchangeability_group( " 0A" );
	} else if ( rsd.aa() == na_rcy ) {
		rsd.name3( " 0C" );
		rsd.interchangeability_group( " 0C" );
	} else if ( rsd.aa() == na_rgu ) {
		rsd.name3( " 0G" );
		rsd.interchangeability_group( " 0G" );
	} else if ( rsd.aa() == na_ura ) {
		rsd.name3( " 0U" );
		rsd.interchangeability_group( " 0U" );
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
	if ( rsd.natoms() < 4 ) return true;
	
	std::string base_atom_name;
	// The rest of the Patch takes care of props, and it's really out of theme.
	if ( rsd.is_l_aa() ) {
		//rsd.delete_property( "L_AA" );
		//rsd.add_property( "D_AA" );
		base_atom_name = " N  ";
	} else if ( rsd.is_d_rna() ) {
		//rsd.delete_property( "D_RNA" );
		//rsd.add_property( "L_RNA" );
		base_atom_name = " P  ";
	} else {
		utility_exit_with_message( "For some reason, calling ChiralFlipAtoms on a non-AA, non-NA" );
	}

	for ( core::Size ii = 2; ii <= rsd.natoms(); ++ii ) {

		AtomICoor icoor1 = rsd.icoor( ii );

		std::string n1 = icoor1.stub_atom1().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom1().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom1().atomno() == 0 ? base_atom_name : rsd.atom( icoor1.stub_atom1().atomno() ).name() ) );
		std::string n2 = icoor1.stub_atom2().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom2().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom2().atomno() == 0 ? base_atom_name : rsd.atom( icoor1.stub_atom2().atomno() ).name() ) );
		std::string n3 = icoor1.stub_atom3().is_polymer_upper() ? "UPPER" :
			( icoor1.stub_atom3().is_polymer_lower() ? "LOWER" :
			( icoor1.stub_atom3().atomno() == 0 ? base_atom_name : rsd.atom( icoor1.stub_atom3().atomno() ).name() ) );

		rsd.set_icoor( rsd.atom( ii ).name(),
			-1.0*icoor1.phi(),
			icoor1.theta(),
			icoor1.d(),
			n1, n2, n3,
			true );
	}
	rsd.fill_ideal_xyz_from_icoor();

	return false;
}

/// @brief replace proton with trifluoromethyl
bool
ReplaceProtonWithTrifluoromethyl::apply( ResidueType & rsd ) const {

	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	//std::cout << "Proton is " << rsd.atom_name( proton_index ) << std::endl;
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	//std::cout << "Adding atom C" << atom_ << " to " << atom_ << ", bonded to " << rsd.atom_name( icoor.stub_atom1().atomno() ) << rsd.atom_name( icoor.stub_atom2().atomno() ) << rsd.atom_name( icoor.stub_atom3().atomno() ) << std::endl;
	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.511005,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		true );

	rsd.set_icoor( "F"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.319661,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		"E"+atom_,
		true );

	rsd.set_icoor( "G"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.319661,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		"F"+atom_,
		true );

	rosetta_recharge_fullatom( rsd );
	find_bonds_in_rings( rsd );
	return false;
}

/// @brief replace proton with methyl
bool
ReplaceProtonWithMethyl::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.511005,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		true );

	rsd.set_icoor( "Y"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		"X"+atom_,
		true );

	rsd.set_icoor( "Z"+atom_,
		-120.0/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
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

/// @brief replace proton with methoxy
bool
ReplaceProtonWithMethoxy::apply( ResidueType & rsd ) const {
	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "O"+atom_, "OH", "OE", -0.61 );
	rsd.add_bond( "O"+atom_, atom_ );
	rsd.set_icoor( "O"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.375964, // actually that's the tyr dist
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
		true  );


	rsd.add_atom( "C"+atom_, "CH3", "CT3", -0.22 );
	rsd.add_bond( "C"+atom_, "O"+atom_ );
	rsd.set_icoor( "C"+atom_,
		-177.688373/180.0*3.14159,
		65.681781/180.0*3.14159,
		1.393527, // aliphatic context.
		"O"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
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
		rsd.atom_name( icoor.stub_atom1().atomno() ),
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
	rsd.add_chi( rsd.atom_name( icoor.stub_atom2().atomno() ), rsd.atom_name( icoor.stub_atom1().atomno() ), "O"+atom_, "C"+atom_ );

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

/// @brief replace proton with ethyl
bool
ReplaceProtonWithEthyl::apply( ResidueType & rsd ) const {
	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "C"+atom_, "CH2", "CT2", -0.27 );
	rsd.add_bond( "C"+atom_, atom_ );
	rsd.set_icoor( "C"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.511005,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		true );

	rsd.set_icoor( "A"+atom_,
		120/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		"D"+atom_,
		true );

	rsd.set_icoor( "B"+atom_,
		120/180.0*3.14159,
		70.5/180.0*3.14159,
		1.083573,
		"C"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
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
		rsd.atom_name( icoor.stub_atom1().atomno() ),
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
	rsd.add_chi( rsd.atom_name( icoor.stub_atom2().atomno() ), rsd.atom_name( icoor.stub_atom1().atomno() ), "C"+atom_, "D"+atom_ );

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

/// @brief replace proton with chlorine
bool
ReplaceProtonWithChlorine::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "CL"+atom_.substr( 1 ), "Cl", "CL", -0.10 );
	rsd.add_bond( "CL"+atom_.substr( 1 ), atom_ );
	rsd.set_icoor( "CL"+atom_.substr( 1 ),
		icoor.phi(),
		icoor.theta(),
		1.810005,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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

/// @brief replace proton with fluorine
bool
ReplaceProtonWithFluorine::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "F"+atom_, "F", "F1", -0.22 );
	rsd.add_bond( "F"+atom_, atom_ );
	rsd.set_icoor( "F"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.332965,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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

/// @brief replace proton with bromine
bool
ReplaceProtonWithBromine::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "BR"+atom_.substr( 1 ), "Br", "BR", -0.07 );
	rsd.add_bond( "BR"+atom_.substr( 1 ), atom_ );
	rsd.set_icoor( "BR"+atom_.substr( 1 ),
		icoor.phi(),
		icoor.theta(),
		1.903609,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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

/// @brief replace proton with iodine
bool
ReplaceProtonWithIodine::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "I"+atom_, "I", "I", -0.08 );
	rsd.add_bond( "I"+atom_, atom_ );
	rsd.set_icoor( "I"+atom_,
		icoor.phi(),
		icoor.theta(),
		2.142205,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
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

/// @brief replace proton with hydroxyl
bool
ReplaceProtonWithHydroxyl::apply( ResidueType & rsd ) const {

	// What is the proton to be replaced?
	if ( rsd.number_bonded_hydrogens( rsd.atom_index( atom_ ) ) == 0 ) {
		return true;
	}

	// Now the valence isn't free.
	rsd.atom( rsd.atom_index( atom_ ) ).set_property( chemical::AROMATIC_CARBON_WITH_FREE_VALENCE, false );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( rsd.atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( rsd.bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( rsd.atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					rsd.atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	Size proton_index = rsd.attached_H_begin( rsd.atom_index( atom_ ) );
	AtomICoor icoor = rsd.icoor( proton_index );
	rsd.delete_atom( proton_index );

	rsd.add_atom( "O"+atom_, "OH", "OH1", -0.54 );
	rsd.add_bond( "O"+atom_, atom_ );
	rsd.set_icoor( "O"+atom_,
		icoor.phi(),
		icoor.theta(),
		1.375964,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		rsd.atom_name( icoor.stub_atom3().atomno() ),
		true );

	rsd.add_atom( "H"+atom_, "Hpol", "H", 0.43 );
	rsd.add_bond( "H"+atom_, "O"+atom_ );

	rsd.set_icoor( "H"+atom_,
		0.000000*3.14159/180.000000,
		70.600000*3.14159/180.000000,
		0.960239,
		"O"+atom_,
		rsd.atom_name( icoor.stub_atom1().atomno() ),
		rsd.atom_name( icoor.stub_atom2().atomno() ),
		true );

	rsd.add_chi( rsd.atom_name( icoor.stub_atom2().atomno() ), rsd.atom_name( icoor.stub_atom1().atomno() ), "O"+atom_, "H"+atom_ );

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

bool
AddConnectDeleteChildProton::apply( ResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!
	rsd.add_metapatch_connect( atom_ );
	return false;
}

bool
DeleteChildProton::apply( ResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!
	rsd.delete_child_proton( atom_ );
	return false;
}

bool
VirtualizeAll::apply( ResidueType & rsd ) const {
	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		rsd.set_atom_type( rsd.atom_name( ii ), "VIRT" );
		rsd.set_mm_atom_type( rsd.atom_name( ii ), "VIRT" );
		rsd.atom( ii ).charge( 0.0 );
		rsd.atom( ii ).is_virtual( true );
	}
	return false;
}

bool
AddConnectAndTrackingVirt::apply( ResidueType & rsd ) const {
	if ( !rsd.has( atom_ ) ) return true; // failure!

	// Provide unique variant name
	std::string res_varname( atom_ + "-METAL_CONNECT" );
	Size count=0;
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for ResidueType: " + rsd.name() );
		}
		++count;
		if ( count == 1 ) {
			if ( ! rsd.has_variant_type( res_varname ) ) break;
		} else {
			res_varname = atom_ + "-METAL_CONNECT" + utility::to_string( count );
			if ( ! rsd.has_variant_type( res_varname ) ) break;
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
		rsd.set_atom_base( virtname, atom_ );
		rsd.set_icoor( virtname, 1.0, 1.0, 1.37, rsd.atom_name( 1 ), ( "V" + ObjexxFCL::string_of( virtcount-2 ) ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), true );
	}
	//icoor = rsd.icoor( rsd.atom_index( virtname ) );

	// Okay, give CONN the icoor of the selected virt.
	if ( virtcount <= 3 ) { // like FE
		rsd.set_icoor( "CONN"+ObjexxFCL::string_of( con_res ), 0, 1.0, 1.37, rsd.atom_name( 1 ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), ( "V" + ObjexxFCL::string_of( virtcount ) ), true );
	} else {
		rsd.set_icoor( "CONN"+ObjexxFCL::string_of( con_res ), 10, 1.0, 1.37, rsd.atom_name( 1 ), ( "V" + ObjexxFCL::string_of( virtcount-2 ) ), ( "V" + ObjexxFCL::string_of( virtcount-1 ) ), true );
	}
	return false;
}

PatchOperationOP
patch_operation_from_patch_file_line(
	std::string const & line,
	std::map< std::string, Real > const & atomic_charge_reassignments
) {
	using numeric::conversions::radians;
	std::istringstream l( line );
	std::string tag, atom1, atom2, atom3, atom4, atom_name, atom_alias, atom_type_name, mm_atom_type_name, bond_type,
		property, variant, dummy;
	Real charge, mean, sdev, radius;
	Size chino;
	SSize formal_charge;
	l >> tag;
	if ( l.fail() || tag[0] == '#' ) return nullptr;
	if ( tag == "ADD_ATOM" ) {
		if ( line.size() < 25 ) return nullptr;
		atom_name = line.substr( 9,4); l >> tag;
		l >> atom_type_name; // = line.substr( 14,4);
		l >> mm_atom_type_name; // = line.substr( 19,4);
		l >> charge;
		if ( l.fail() ) return nullptr;

		//fd let command line override charge
		if ( atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) ) != atomic_charge_reassignments.end() ) {
			tr.Trace << "reassigning patch atomic charge " << atom_name << " atomtype: " << atom_type_name << " --> " <<
				atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second << std::endl;
			// note that we set charge and also parse_charge, so this will over-ride the parse_charge if those are the charges we are using
			charge = atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second;
		}

#if defined(WIN32) && !defined(WIN_PYROSETTA)
		return PatchOperationOP( new AddAtomWIN32( atom_name, atom_type_name, mm_atom_type_name, charge ) );
#else
		return PatchOperationOP( new AddAtom( atom_name, atom_type_name, mm_atom_type_name, charge ) );
#endif
	} else if ( tag == "DELETE_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new DeleteAtom( atom_name ) );

	} else if ( tag == "ADD_ATOM_ALIAS" ) {
		l >> atom_name >> atom_alias;
		if ( l.fail() ) {
			return nullptr;
		}
		return PatchOperationOP( new AddAtomAlias( atom_name, atom_alias ) );

	} else if ( tag == "SET_BACKBONE_HEAVYATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new SetBackboneHeavyatom( atom_name ) );

	} else if ( tag == "SET_IO_STRING" ) { // 13 character tag
		// NOTE - USE FIXED WIDTH IO SINCE NAME3 CAN CONTAIN INTERNAL WHITESPACE (EG DNA,RNA)
		if ( line.size() < 19 ) return nullptr;
		std::string const three_letter_code( line.substr(14,3) ), one_letter_code( line.substr(18,1) );
		return PatchOperationOP( new SetIO_String( three_letter_code, one_letter_code[0] ) );

	} else if ( tag == "SET_INTERCHANGEABILITY_GROUP" ) {
		std::string intgrp;
		l >> intgrp;
		return PatchOperationOP( new SetInterchangeabilityGroup_String( intgrp ) );

	} else if ( tag == "NBR_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new SetNbrAtom( atom_name ) );

	} else if ( tag == "NBR_RADIUS" ) {
		Real radius;
		l >> radius;
		return PatchOperationOP( new SetNbrRadius( radius ) );

	} else if ( tag == "ADD_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new AddProperty( property ) );

	} else if ( tag == "DELETE_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new DeleteProperty( property ) );

	} else if ( tag == "DELETE_VARIANT_TYPE" ) {
		l >> variant;
		if ( l.fail() ) { return nullptr; }
		return PatchOperationOP( new DeleteVariantType( variant ) );

		// Added by Andy M. Chen in June 2009
		// This is needed for adding new side-chain torsions, which occurs in certain PTMs.
	} else if ( tag == "ADD_CHI" ) {
		if ( line.substr(8, 3) == "N+1" ) {
			l >> dummy >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) return nullptr;
			return PatchOperationOP( new AddChi(atom1, atom2, atom3, atom4) );
		} else {
			l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) return nullptr;
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
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new AddProtonChi( chino, samples, extra_samples ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "REDEFINE_CHI" ) {
		l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new RedefineChi( chino, atom1, atom2, atom3, atom4 ) );

	} else if ( tag == "DELETE_TERMINAL_CHI" ) {
		return PatchOperationOP( new DeleteTerminalChi() );
	} else if ( tag == "DELETE_METALBINDING_ATOM" ) {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) return nullptr;

		return PatchOperationOP( new DeleteMetalbindingAtom( atom_name ) );
	} else if ( tag == "DELETE_ACT_COORD_ATOM" ) {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) return nullptr;

		return PatchOperationOP( new DeleteActCoordAtom( atom_name ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_CHI_ROTAMER" ) {
		if ( line.substr(16, 1) == "N" ) {
			l >> dummy >> mean >> sdev;
			if ( l.fail() ) return nullptr;
			return PatchOperationOP( new AddChiRotamer(mean, sdev) );
		} else {
			l >> chino >> mean >> sdev;
			if ( l.fail() ) return nullptr;
			return PatchOperationOP( new AddChiRotamer(chino, mean, sdev) );
		}

	} else if ( tag == "CLEAR_CHI_ROTAMERS" ) {
		l >> chino;
		if ( l.fail() ) { return nullptr; }
		return PatchOperationOP( new ClearChiRotamers( chino ) );

		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_BOND" ) {
		l >> atom1 >> atom2;
		if ( l.fail() ) return nullptr;
		return PatchOperationOP( new AddBond( atom1, atom2 ) );

	} else if ( tag == "ADD_BOND_TYPE" ) {
		l >> atom1 >> atom2 >> bond_type;
		if ( l.fail() ) {
			return nullptr;
		}
		return PatchOperationOP( new AddBondType( atom1, atom2, bond_type ) );

	} else if ( tag == "CHANGE_BOND_TYPE" ) {
		std::string old_bond_type;
		l >> atom1 >> atom2 >> old_bond_type >> dummy >> bond_type;
		if ( l.fail() ) {
			return nullptr;
		}
		return PatchOperationOP( new ChangeBondType( atom1, atom2, old_bond_type, bond_type ) );

	} else if ( tag == "ADD_CONNECT" ) {
		std::string connect_atom;
		l >> connect_atom;
		if ( l.fail() ) return nullptr;
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
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetAtomType( atom_name, atom_type_name ) );

	} else if ( tag == "SET_MM_ATOM_TYPE" ) {
		runtime_assert( l.good() );
		l >> atom_name >> mm_atom_type_name;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetMMAtomType( atom_name, mm_atom_type_name ) );

	} else if ( tag == "SET_FORMAL_CHARGE" ) {
		l >> atom_name >> formal_charge;
		if ( l.fail() ) {
			return nullptr;
		}
		return PatchOperationOP( new SetFormalCharge( atom_name, formal_charge ) );

	} else if ( tag == "SET_ATOMIC_CHARGE" ) {
		l >> atom_name >> charge;

		//fd let command line override charge
		if ( atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) ) != atomic_charge_reassignments.end() ) {
			tr.Trace << "reassigning patch atomic charge " << atom_name << " atomtype: " << atom_type_name << " --> " <<
				atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second << std::endl;
			// note that we set charge and also parse_charge, so this will over-ride the parse_charge if those are the charges we are using
			charge = atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second;
		}


		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetAtomicCharge( atom_name, charge ) );

	} else if ( tag == "SET_POLYMER_CONNECT" ) {
		l >> tag >> atom_name; // tag should be "UPPER" or "LOWER"
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetPolymerConnectAtom( atom_name, tag ) );

	} else if ( tag == "SET_ICOOR" ) {
		Real phi,theta,d;
		std::string stub1, stub2, stub3;
		l >> atom_name >> phi >> theta >> d >> stub1 >> stub2 >> stub3;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetICoor( atom_name, radians(phi), radians(theta), d, stub1, stub2, stub3 ) );
	} else if ( tag == "SET_ANCESTOR" ) {
		std::string atom_name, which_anc, anc_atom_name;
		Ancestor anc( anc_parent );
		l >> atom_name >> which_anc >> anc_atom_name;
		if ( which_anc == "PARENT" ) {
			anc = anc_parent;
		} else if ( which_anc == "GRANDPARENT" ) {
			anc = anc_grandparent;
		} else if ( which_anc == "GREATGRANDPARENT" ) {
			anc = anc_greatgrandparent;
		} else {
			std::ostringstream oss;
			oss << "While reading the SET_ANCESTOR patch operation line, did not find PARENT, GRANDPARENT, or GREATGRANDPARENT as the third string on that line\n";
			throw utility::excn::EXCN_Msg_Exception( oss.str() );
		}
		return PatchOperationOP( new ChangeAncestory( atom_name, anc, anc_atom_name ));
	} else if ( tag == "RESET_BOND_LENGTH" ) {
		core::Distance d;
		l >> atom_name >> d;
		if ( l.fail() ) {
			utility_exit_with_message( "bad line in patchfile: " + line );
		}
		return PatchOperationOP( new ResetBondLength( atom_name, d ) );

	} else if ( tag == "PREPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new PrependMainchainAtom( atom_name ) );

	} else if ( tag == "APPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new AppendMainchainAtom( atom_name ) );

	}  else if ( tag == "REPLACE_MAINCHAIN_ATOM" ) {
		std::string target, new_atom;
		l >> target >> new_atom;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceMainchainAtom( target, new_atom ) );
	} else if ( tag == "RAMA_PREPRO_FILENAME" ) {
		std::string file1, file2;
		l >> file1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		l >> file2;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new RamaPreproFilename( file1, file2 ) );
	} else if ( tag == "RAMA_PREPRO_RESNAME" ) {
		std::string name;
		l >> name;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new RamaPreproResname( name ) );
	} else if ( tag == "REMOVE_ROTAMER_SPECIFICATIONS" ) {
		return PatchOperationOP( new RemoveRotamerSpecifications() );
	} else if ( tag == "NCAA_ROTLIB_PATH" ) {
		std::string path;
		l >> path;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new NCAARotLibPath( path ) );
	} else if ( tag == "NCAA_ROTLIB_NUM_ROTAMER_BINS" ) {
		core::Size n_rots(0);
		utility::vector1<core::Size> n_bins_per_rot;
		l >> n_rots;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		n_bins_per_rot.resize( n_rots );
		for ( Size i = 1; i <= n_rots; ++i ) {
			Size bin_size(0);
			l >> bin_size;
			if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
			n_bins_per_rot[i] = bin_size;
		}
		return PatchOperationOP( new NCAARotLibNumRotamerBins( n_bins_per_rot ) );
	} else if ( tag == "SET_NBR_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new SetNbrAtom( atom_name ) );

	} else if ( tag == "SET_NBR_RADIUS" ) {
		l >> radius;
		if ( l.fail() ) utility_exit_with_message("bad line in patchfile: " + line );
		return PatchOperationOP( new SetNbrRadius( radius ) );

	} else if ( tag == "SET_ORIENT_ATOM" ) {
		l >> tag;
		if ( l.fail() ) utility_exit_with_message("bad line in patchfile: " +  line );
		if ( tag == "NBR" ) {
			return PatchOperationOP( new SetOrientAtom(true) );
		} else if ( tag == "DEFAULT" ) {
			return PatchOperationOP( new SetOrientAtom(false) );
		} else {
			tr.Warning << "Unknown SET_ORIENT ATOM tag: " << tag << std::endl;
			return nullptr;
		}
	} else if ( tag == "SET_ALL_ATOMS_REPULSIVE" ) {
		return PatchOperationOP( new SetAllAtomsRepulsive() );
	} else if ( tag == "CONNECT_SULFUR_AND_MAKE_VIRTUAL_PROTON" ) {
		return PatchOperationOP( new ConnectSulfurAndMakeVirtualProton() );
	} else if ( tag == "REPLACE_PROTON_WITH_CHLORINE" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithChlorine( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_FLUORINE" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithFluorine( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_BROMINE" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithBromine( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_IODINE" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithIodine( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_METHYL" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithMethyl( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_TRIFLUOROMETHYL" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithTrifluoromethyl( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_HYDROXYL" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithHydroxyl( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_METHOXY" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithMethoxy( atom1 ) );
	} else if ( tag == "REPLACE_PROTON_WITH_ETHYL" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new ReplaceProtonWithEthyl( atom1 ) );
	} else if ( tag == "ADD_CONNECT_AND_DELETE_CHILD_PROTON" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new AddConnectDeleteChildProton( atom1 ) );
	} else if ( tag == "DELETE_CHILD_PROTON" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new DeleteChildProton( atom1 ) );
	} else if ( tag == "ADD_CONNECT_AND_TRACKING_VIRT" ) {
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad line in patchfile: " + line );
		return PatchOperationOP( new AddConnectAndTrackingVirt( atom1 ) );
	} else if ( tag == "CHIRAL_FLIP_NAMING" ) {
		return PatchOperationOP( new ChiralFlipNaming );//( atom1, atom2 ) );
	} else if ( tag == "CHIRAL_FLIP_ATOMS" ) {
		//std::string atom1, atom2;
		//l >> atom1 >> atom2;
		return PatchOperationOP( new ChiralFlipAtoms );//( atom1, atom2 ) );
	} else if ( tag == "VIRTUALIZE_ALL" ) {
		return PatchOperationOP( new VirtualizeAll );//( atom1, atom2 ) );
	}
	tr.Warning << "patch_operation_from_patch_file_line: bad line: " << line << std::endl;

	return nullptr;
}

} // chemical
} // core
