// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/PatchOperationFactory.cc
/// @brief  Generate PatchOperations from patchfile lines.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/PatchOperationFactory.hh>

// Package headers
#include <core/chemical/PatchOperation.hh>

#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption

#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.PatchOperationFactory" );

/// @brief add a PatchOperation prototype, using its default type name as the map key
void
PatchOperationFactory::factory_register( PatchOperationCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const & type( creator->keyname() );
	if ( creator_map_.count( type ) != 0 ) {
		utility_exit_with_message( "PatchOperationFactory::factory_register already has a mover creator with name \"" + type + "\".  Conflicting PatchOperation names" );
	}
	creator_map_[ type ] = creator;
}

PatchOperationOP
PatchOperationFactory::newPatchOperation( std::string const & tag, std::string const & line, std::istream & input, std::map< std::string, Real > const & atomic_charge_reassignments ) const
{
	CreatorMap::const_iterator iter( creator_map_.find( tag ) );
	if ( iter != creator_map_.end() ) {
		debug_assert( iter->second != nullptr );
		return iter->second->create_operation(line, input, atomic_charge_reassignments);
	} else {
		TR.Warning << "Unable to understand Patch operation type " << tag << std::endl;
		TR.Warning << "  Availible operation types: " << '\n';
		for ( auto const & entry: creator_map_ ) {
			TR.Warning << entry.first << ", ";
		}
		TR.Warning << std::endl;
		return nullptr; // Matches original behavior of patch_operation_from_patch_file_line()
	}
}


///////////////////////////////////////////////////////////////////////////////
// The actual Operations.

class AddAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & atomic_charge_reassignments) const override {
		std::string atom_name, atom_type_name, mm_atom_type_name;
		core::Real charge;

		if ( line.size() < 25 ) {
			utility_exit_with_message( "Could not parse an ADD_ATOM line.  Note that the first 5 characters following the \"ADD_ATOM\" tag are whitespace-sensitive, since atom names have different whitespace paddings in PDB files.\nLine: " + line );
			return nullptr;
		}
		l >> atom_name; //Needed to advance cursor, although overwritten by next line.
		atom_name = line.substr( 9,4);
		l >> atom_type_name; // = line.substr( 14,4);
		l >> mm_atom_type_name; // = line.substr( 19,4);
		l >> charge;
		if ( l.fail() ) {
			utility_exit_with_message( "Could not parse an ADD_ATOM line.  Note that the first 5 characters following the \"ADD_ATOM\" tag are whitespace-sensitive, since atom names have different whitespace paddings in PDB files.\nLine: " + line );
			return nullptr;
		}

		//fd let command line override charge
		if ( atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) ) != atomic_charge_reassignments.end() ) {
			TR.Trace << "reassigning patch atomic charge " << atom_name << " atomtype: " << atom_type_name << " --> " <<
				atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second << std::endl;
			// note that we set charge and also parse_charge, so this will over-ride the parse_charge if those are the charges we are using
			charge = atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second;
		}

		return utility::pointer::make_shared< AddAtom >( atom_name, atom_type_name, mm_atom_type_name, charge );
	}
};

static PatchOperationRegistrator< AddAtomCreator > reg_AddAtomCreator;

///////////////

class DeleteAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_ATOM"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse DELETE_ATOM patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< DeleteAtom >( atom_name );
	}
};

static PatchOperationRegistrator< DeleteAtomCreator > reg_DeleteAtomCreator;

///////////////

class AddAtomAliasCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_ATOM_ALIAS"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name, atom_alias;
		utility::vector1< std::string > atom_aliases;
		std::string const atom_name_full( line.substr( 15, 4 ) );  // Rosetta atom, including whitespace
		std::string const primary_alias( line.substr( 20, 4) );  // The preferred alias used in most cases by the PDB.
		l >> atom_name >> atom_alias;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_ATOM_ALIAS patch operation." );
			return nullptr;
		}
		while ( ! l.fail() ) {
			atom_aliases.push_back( atom_alias );
			l >> atom_alias;
		}
		return utility::pointer::make_shared< AddAtomAlias >( atom_name_full, atom_aliases, primary_alias );
	}
};

static PatchOperationRegistrator< AddAtomAliasCreator > reg_AddAtomAliasCreator;

///////////////

class SetBackboneHeavyatomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_BACKBONE_HEAVYATOM"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse SET_BACKBONE_HEAVYATOM patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetBackboneHeavyatom >( atom_name );
	}
};

static PatchOperationRegistrator< SetBackboneHeavyatomCreator > reg_SetBackboneHeavyatomCreator;

///////////////

class SetDisulfideAtomNameCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_DISULFIDE_ATOM_NAME"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse SET_DISULFIDE_ATOM_NAME patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetDisulfideAtomName >( atom_name );

	}
};

static PatchOperationRegistrator< SetDisulfideAtomNameCreator > reg_SetDisulfideAtomNameCreator;

///////////////

class SetAACreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_AA"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string aa;
		l >> aa;
		return utility::pointer::make_shared< Set_AA >( aa );
	}
};

static PatchOperationRegistrator< SetAACreator > reg_SetAACreator;

///////////////

class SetBaseNameCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_BASE_NAME"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string newname;

		l >> newname;
		runtime_assert_string_msg( !newname.empty(), "Could not parse patch operation line \"" + line + "\"." );
		return utility::pointer::make_shared< SetBaseName >( newname );
	}
};

static PatchOperationRegistrator< SetBaseNameCreator > reg_SetBaseNameCreator;

///////////////

class SetIOStringCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_IO_STRING"; }
	PatchOperationOP create_operation(std::string const & line, std::istream &, std::map< std::string, Real > const & ) const override {
		// NOTE - USE FIXED WIDTH IO SINCE NAME3 CAN CONTAIN INTERNAL WHITESPACE (EG DNA,RNA)
		if ( line.size() < 19 ) {
			utility_exit_with_message( "Failed to parse SET_IO_STRING patch operation.  Note that this is a whitespace-sensitive line.\nLine: " + line );
			return nullptr;
		}
		std::string const three_letter_code( line.substr(14,3) ), one_letter_code( line.substr(18,1) );
		return utility::pointer::make_shared< SetIO_String >( three_letter_code, one_letter_code[0] );
	}
};

static PatchOperationRegistrator< SetIOStringCreator > reg_SetIOStringCreator;

///////////////

class SetInterchangeabilityGroupCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_INTERCHANGEABILITY_GROUP"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string intgrp;
		l >> intgrp;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse SET_INTERCHANGEABILITY_GROUP patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetInterchangeabilityGroup_String >( intgrp );
	}
};

static PatchOperationRegistrator< SetInterchangeabilityGroupCreator > reg_SetInterchangeabilityGroupCreator;

///////////////

class NbrAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "NBR_ATOM"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse NBR_ATOM patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetNbrAtom >( atom_name );
	}
};

static PatchOperationRegistrator< NbrAtomCreator > reg_NbrAtomCreator;

///////////////

class NbrRadiusCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "NBR_RADIUS"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		Real radius;
		l >> radius;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse NBR_RADIUS patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetNbrRadius >( radius );
	}
};

static PatchOperationRegistrator< NbrRadiusCreator > reg_NbrRadiusCreator;

///////////////

class AddPropertyCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_PROPERTY"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string property;
		l >> property;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_PROPERTY patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< AddProperty >( property );
	}
};

static PatchOperationRegistrator< AddPropertyCreator > reg_AddPropertyCreator;

///////////////

class DeletePropertyCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_PROPERTY"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string property;
		l >> property;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse DELETE_PROPERTY patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< DeleteProperty >( property );
	}
};

static PatchOperationRegistrator< DeletePropertyCreator > reg_DeletePropertyCreator;

///////////////

class DeleteVariantTypeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_VARIANT_TYPE"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string variant;
		l >> variant;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse DELETE_VARIANT_TYPE patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< DeleteVariantType >( variant );
	}
};

static PatchOperationRegistrator< DeleteVariantTypeCreator > reg_DeleteVariantTypeCreator;

///////////////

class AddChiCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_CHI"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1, atom2, atom3, atom4, dummy;
		core::Size chino;
		if ( line.substr(8, 3) == "N+1" ) {
			l >> dummy >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) {
				utility_exit_with_message( "Failed to parse ADD_CHI patch operation." );
				return nullptr;
			}
			return utility::pointer::make_shared< AddChi >(atom1, atom2, atom3, atom4);
		} else {
			l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
			if ( l.fail() ) {
				utility_exit_with_message( "Failed to parse ADD_CHI patch operation." );
				return nullptr;
			}
			return utility::pointer::make_shared< AddChi >(chino, atom1, atom2, atom3, atom4);
		}
	}
};

static PatchOperationRegistrator< AddChiCreator > reg_AddChiCreator;

///////////////

class AddProtonChiCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_PROTON_CHI"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		Size chino, nsamples, nextra_samples;
		std::string dummy2;
		l >> chino;
		l >> dummy2; // should be "SAMPLES"
		l >> nsamples;
		utility::vector1< Real > samples( nsamples );
		for ( Size ii = 1; ii <= nsamples; ++ii ) {
			l >> samples[ ii ];
		}
		l >> dummy2; // should be "EXTRA"
		l >> nextra_samples;
		utility::vector1< Real > extra_samples( nextra_samples );
		for ( Size ii = 1; ii <= nextra_samples; ++ii ) {
			l >> extra_samples[ ii ];
		}
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_PROTON_CHI patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< AddProtonChi >( chino, samples, extra_samples );
	}
};

static PatchOperationRegistrator< AddProtonChiCreator > reg_AddProtonChiCreator;

///////////////

class RedefineChiCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REDEFINE_CHI"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1, atom2, atom3, atom4;
		core::Size chino;
		l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse REDEFINE_CHI patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< RedefineChi >( chino, atom1, atom2, atom3, atom4 );
	}
};

static PatchOperationRegistrator< RedefineChiCreator > reg_RedefineChiCreator;

///////////////

class DeleteTerminalChiCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_TERMINAL_CHI"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< DeleteTerminalChi >();
	}
};

static PatchOperationRegistrator< DeleteTerminalChiCreator > reg_DeleteTerminalChiCreator;

///////////////

class DeleteMetalBindingAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_METALBINDING_ATOM"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse DELETE_METALBINDING_ATOM patch operation." );
			return nullptr;
		}

		return utility::pointer::make_shared< DeleteMetalbindingAtom >( atom_name );
	}
};

static PatchOperationRegistrator< DeleteMetalBindingAtomCreator > reg_DeleteMetalBindingAtomCreator;

///////////////

class DeleteActCoordAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_ACT_COORD_ATOM"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse DELETE_ACT_COORD_ATOM patch operation." );
			return nullptr;
		}

		return utility::pointer::make_shared< DeleteActCoordAtom >( atom_name );
	}
};

static PatchOperationRegistrator< DeleteActCoordAtomCreator > reg_DeleteActCoordAtomCreator;

///////////////

class AddChiRotamerCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_CHI_ROTAMER"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string dummy;
		Size chino;
		Real mean, sdev;
		if ( line.substr(16, 1) == "N" ) {
			l >> dummy >> mean >> sdev;
			if ( l.fail() ) {
				utility_exit_with_message( "Failed to parse ADD_CHI_ROTAMER patch operation." );
				return nullptr;
			}
			return utility::pointer::make_shared< AddChiRotamer >(mean, sdev);
		} else {
			l >> chino >> mean >> sdev;
			if ( l.fail() ) {
				utility_exit_with_message( "Failed to parse ADD_CHI_ROTAMER patch operation." );
				return nullptr;
			}
			return utility::pointer::make_shared< AddChiRotamer >(chino, mean, sdev);
		}
	}
};

static PatchOperationRegistrator< AddChiRotamerCreator > reg_AddChiRotamerCreator;

///////////////

class ClearChiRotamersCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "CLEAR_CHI_ROTAMERS"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		Size chino;
		l >> chino;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse CLEAR_CHI_ROTAMERS patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< ClearChiRotamers >( chino );
	}
};

static PatchOperationRegistrator< ClearChiRotamersCreator > reg_ClearChiRotamersCreator;

///////////////

class AddBondCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_BOND"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1, atom2;
		l >> atom1 >> atom2;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_BOND patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< AddBond >( atom1, atom2 );
	}
};

static PatchOperationRegistrator< AddBondCreator > reg_AddBondCreator;

///////////////

class AddBondTypeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_BOND_TYPE"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1, atom2, bond_type;
		l >> atom1 >> atom2 >> bond_type;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_BOND_TYPE patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< AddBondType >( atom1, atom2, bond_type );
	}
};

static PatchOperationRegistrator< AddBondTypeCreator > reg_AddBondTypeCreator;

///////////////

class ChangeBondTypeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "CHANGE_BOND_TYPE"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1, atom2, bond_type;
		std::string old_bond_type, dummy;
		l >> atom1 >> atom2 >> old_bond_type >> dummy >> bond_type;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse CHANGE_BOND_TYPE patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< ChangeBondType >( atom1, atom2, old_bond_type, bond_type );
	}
};

static PatchOperationRegistrator< ChangeBondTypeCreator > reg_ChangeBondTypeCreator;

///////////////

class AddConnectCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_CONNECT"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		using numeric::conversions::radians;
		std::string tag;
		std::string connect_atom;
		l >> connect_atom;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse ADD_CONNECT patch operation." );
			return nullptr;
		}
		l >> tag;
		if ( l.fail() ) {
			return utility::pointer::make_shared< AddConnect >(
				connect_atom, 0.0, 0.0, 0.0, connect_atom, connect_atom, connect_atom );
		} else {
			Real phi, theta, d;
			std::string parent_atom, angle_atom, torsion_atom;
			l >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;
			if ( l.fail() || tag != "ICOOR" ) {
				utility_exit_with_message( "bad line in patchfile: "+line );
			}
			return utility::pointer::make_shared< AddConnect >(
				connect_atom, radians(phi), radians(theta), d, parent_atom, angle_atom, torsion_atom );
		}
	}
};

static PatchOperationRegistrator< AddConnectCreator > reg_AddConnectCreator;

///////////////

class SetAtomTypeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ATOM_TYPE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name, atom_type_name;
		l >> atom_name >> atom_type_name;
		if ( l.fail() ) utility_exit_with_message( "bad SET_ATOM_TYPE line in patchfile: " + line );
		return utility::pointer::make_shared< SetAtomType >( atom_name, atom_type_name );
	}
};

static PatchOperationRegistrator< SetAtomTypeCreator > reg_SetAtomTypeCreator;

///////////////

class SetMMAtomTypeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_MM_ATOM_TYPE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name, mm_atom_type_name;
		runtime_assert( l.good() );
		l >> atom_name >> mm_atom_type_name;
		if ( l.fail() ) utility_exit_with_message( "bad SET_MM_ATOM_TYPE line in patchfile: " + line );
		return utility::pointer::make_shared< SetMMAtomType >( atom_name, mm_atom_type_name );
	}
};

static PatchOperationRegistrator< SetMMAtomTypeCreator > reg_SetMMAtomTypeCreator;

///////////////

class SetFormalChargeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_FORMAL_CHARGE"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		SSize formal_charge;
		l >> atom_name >> formal_charge;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse SET_FORMAL_CHARGE patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetFormalCharge >( atom_name, formal_charge );
	}
};

static PatchOperationRegistrator< SetFormalChargeCreator > reg_SetFormalChargeCreator;

///////////////

class SetNetFormalChargeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_NET_FORMAL_CHARGE"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
		signed int net_formal_charge;
		l >> net_formal_charge;
		if ( l.fail() ) {
			utility_exit_with_message( "Failed to parse SET_NET_FORMAL_CHARGE patch operation." );
			return nullptr;
		}
		return utility::pointer::make_shared< SetNetFormalCharge >( net_formal_charge );
	}
};

static PatchOperationRegistrator< SetNetFormalChargeCreator > reg_SetNetFormalChargeCreator;

///////////////

class SetAtomicChargeCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ATOMIC_CHARGE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & atomic_charge_reassignments) const override {
		std::string atom_name;
		Real charge;
		l >> atom_name >> charge;

		//fd let command line override charge
		if ( atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) ) != atomic_charge_reassignments.end() ) {
			TR.Trace << "reassigning patch atomic charge " << atom_name << " --> " <<
				atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second << std::endl;
			// note that we set charge and also parse_charge, so this will over-ride the parse_charge if those are the charges we are using
			charge = atomic_charge_reassignments.find( ObjexxFCL::stripped( atom_name ) )->second;
		}

		if ( l.fail() ) utility_exit_with_message( "bad SET_ATOMIC_CHARGE line in patchfile: " + line );
		return utility::pointer::make_shared< SetAtomicCharge >( atom_name, charge );
	}
};

static PatchOperationRegistrator< SetAtomicChargeCreator > reg_SetAtomicChargeCreator;

///////////////

class SetPolymerConnectCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_POLYMER_CONNECT"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string tag, atom_name;
		l >> tag >> atom_name; // tag should be "UPPER" or "LOWER"
		if ( l.fail() ) utility_exit_with_message( "bad SET_POLYMER_CONNECT line in patchfile: " + line );
		return utility::pointer::make_shared< SetPolymerConnectAtom >( atom_name, tag );
	}
};

static PatchOperationRegistrator< SetPolymerConnectCreator > reg_SetPolymerConnectCreator;

///////////////

class SetICOORCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ICOOR"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		using numeric::conversions::radians;
		Real phi,theta,d;
		std::string atom_name, stub1, stub2, stub3;
		l >> atom_name >> phi >> theta >> d >> stub1 >> stub2 >> stub3;
		if ( l.fail() ) utility_exit_with_message( "bad SET_ICOOR line in patchfile: " + line );
		return utility::pointer::make_shared< SetICoor >( atom_name, radians(phi), radians(theta), d, stub1, stub2, stub3 );
	}
};

static PatchOperationRegistrator< SetICOORCreator > reg_SetICOORCreator;

///////////////

class SetAncestorCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ANCESTOR"; }
	PatchOperationOP create_operation(std::string const &, std::istream & l, std::map< std::string, Real > const & ) const override {
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
			utility_exit_with_message( "While reading the SET_ANCESTOR patch operation line, did not find PARENT, GRANDPARENT, or GREATGRANDPARENT as the third string on that line." );
		}
		return utility::pointer::make_shared< ChangeAncestory >( atom_name, anc, anc_atom_name );
	}
};

static PatchOperationRegistrator< SetAncestorCreator > reg_SetAncestorCreator;

///////////////

class ResetBondLengthCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "RESET_BOND_LENGTH"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		core::Distance d;
		l >> atom_name >> d;
		if ( l.fail() ) {
			utility_exit_with_message( "bad RESET_BOND_LENGTH line in patchfile: " + line );
		}
		return utility::pointer::make_shared< ResetBondLength >( atom_name, d );
	}
};

static PatchOperationRegistrator< ResetBondLengthCreator > reg_ResetBondLengthCreator;

///////////////

class PrependMainchainAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "PREPEND_MAINCHAIN_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad PREPEND_MAINCHAIN_ATOM line in patchfile: " + line );
		return utility::pointer::make_shared< PrependMainchainAtom >( atom_name );
	}
};

static PatchOperationRegistrator< PrependMainchainAtomCreator > reg_PrependMainchainAtomCreator;

///////////////

class AppendMainchainAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "APPEND_MAINCHAIN_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad APPEND_MAINCHAIN_ATOM line in patchfile: " + line );
		return utility::pointer::make_shared< AppendMainchainAtom >( atom_name );
	}
};

static PatchOperationRegistrator< AppendMainchainAtomCreator > reg_AppendMainchainAtomCreator;

///////////////

class ReplaceMainchainAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_MAINCHAIN_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string target, new_atom;
		l >> target >> new_atom;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_MAINCHAIN_ATOM line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceMainchainAtom >( target, new_atom );
	}
};

static PatchOperationRegistrator< ReplaceMainchainAtomCreator > reg_ReplaceMainchainAtomCreator;

///////////////

class RamaPreproFilenameCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "RAMA_PREPRO_FILENAME"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string file1, file2;
		l >> file1;
		if ( l.fail() ) utility_exit_with_message( "bad RAMA_PREPRO_FILENAME line in patchfile: " + line );
		l >> file2;
		if ( l.fail() ) utility_exit_with_message( "bad RAMA_PREPRO_FILENAME line in patchfile: " + line );
		return utility::pointer::make_shared< RamaPreproFilename >( file1, file2 );
	}
};

static PatchOperationRegistrator< RamaPreproFilenameCreator > reg_RamaPreproFilenameCreator;

///////////////

class RamaPreproResnameCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "RAMA_PREPRO_RESNAME"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string name;
		l >> name;
		if ( l.fail() ) utility_exit_with_message( "bad RAMA_PREPRO_RESNAME line in patchfile: " + line );
		return utility::pointer::make_shared< RamaPreproResname >( name );
	}
};

static PatchOperationRegistrator< RamaPreproResnameCreator > reg_RamaPreproResnameCreator;

///////////////

class RemoveRotamerSpecificationsCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REMOVE_ROTAMER_SPECIFICATIONS"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< RemoveRotamerSpecifications >();
	}
};

static PatchOperationRegistrator< RemoveRotamerSpecificationsCreator > reg_RemoveRotamerSpecificationsCreator;

///////////////

class NCAARotlibPathCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "NCAA_ROTLIB_PATH"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string path;
		l >> path;
		if ( l.fail() ) utility_exit_with_message( "bad NCAA_ROTLIB_PATH line in patchfile: " + line );
		return utility::pointer::make_shared< NCAARotLibPath >( path );
	}
};

static PatchOperationRegistrator< NCAARotlibPathCreator > reg_NCAARotlibPathCreator;

///////////////

class NCAARotlibBBTorsionsCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "NCAA_ROTLIB_BB_TORSIONS"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		utility::vector1< core::Size > ncaa_rotlib_bb_torsions;
		while ( !l.eof() ) {
			core::Size index;
			l >> index;
			if ( l.fail() ) utility_exit_with_message( "Bad NCAA_ROTLIB_BB_TORSIONS line in patchfile: " + line );
			ncaa_rotlib_bb_torsions.push_back(index);
		}
		return utility::pointer::make_shared< NCAARotLibBBTorsions >( ncaa_rotlib_bb_torsions );
	}
};

static PatchOperationRegistrator< NCAARotlibBBTorsionsCreator > reg_NCAARotlibBBTorsionsCreator;

///////////////

class NCAARotlibNumRotamerBinsCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "NCAA_ROTLIB_NUM_ROTAMER_BINS"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		core::Size n_rots(0);
		utility::vector1<core::Size> n_bins_per_rot;
		l >> n_rots;
		if ( l.fail() ) utility_exit_with_message( "bad NCAA_ROTLIB_NUM_ROTAMER_BINS line in patchfile: " + line );
		n_bins_per_rot.resize( n_rots );
		for ( Size i = 1; i <= n_rots; ++i ) {
			Size bin_size(0);
			l >> bin_size;
			if ( l.fail() ) utility_exit_with_message( "bad NCAA_ROTLIB_NUM_ROTAMER_BINS line in patchfile: " + line );
			n_bins_per_rot[i] = bin_size;
		}
		return utility::pointer::make_shared< NCAARotLibNumRotamerBins >( n_bins_per_rot );
	}
};

static PatchOperationRegistrator< NCAARotlibNumRotamerBinsCreator > reg_NCAARotlibNumRotamerBinsCreator;

///////////////

class SetNbrAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_NBR_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom_name;
		l >> atom_name;
		if ( l.fail() ) utility_exit_with_message( "bad SET_NBR_ATOM line in patchfile: " + line );
		return utility::pointer::make_shared< SetNbrAtom >( atom_name );
	}
};

static PatchOperationRegistrator< SetNbrAtomCreator > reg_SetNbrAtomCreator;

///////////////

class SetNbrRadiusCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_NBR_RADIUS"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		Real radius;
		l >> radius;
		if ( l.fail() ) utility_exit_with_message("bad SET_NBR_RADIUS line in patchfile: " + line );
		return utility::pointer::make_shared< SetNbrRadius >( radius );
	}
};

static PatchOperationRegistrator< SetNbrRadiusCreator > reg_SetNbrRadiusCreator;

///////////////

class SetOrientAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ORIENT_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string tag;
		l >> tag;
		if ( l.fail() ) utility_exit_with_message("bad SET_ORIENT_ATOM line in patchfile: " +  line );
		if ( tag == "NBR" ) {
			return utility::pointer::make_shared< SetOrientAtom >(true);
		} else if ( tag == "DEFAULT" ) {
			return utility::pointer::make_shared< SetOrientAtom >(false);
		} else {
			TR.Warning << "Unknown SET_ORIENT ATOM tag: " << tag << std::endl;
			return nullptr;
		}
	}
};

static PatchOperationRegistrator< SetOrientAtomCreator > reg_SetOrientAtomCreator;

///////////////

class SetAllAtomsRepulsiveCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_ALL_ATOMS_REPULSIVE"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< SetAllAtomsRepulsive >();
	}
};

static PatchOperationRegistrator< SetAllAtomsRepulsiveCreator > reg_SetAllAtomsRepulsiveCreator;

///////////////

class ConnectSulfurAndMakeVirtualProtonCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "CONNECT_SULFUR_AND_MAKE_VIRTUAL_PROTON"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< ConnectSulfurAndMakeVirtualProton >();
	}
};

static PatchOperationRegistrator< ConnectSulfurAndMakeVirtualProtonCreator > reg_ConnectSulfurAndMakeVirtualProtonCreator;

///////////////

class ReplaceProtonWithChlorineCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_CHLORINE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_CHLORINE line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithChlorine >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithChlorineCreator > reg_ReplaceProtonWithChlorineCreator;

///////////////

class ReplaceProtonWithFluorineCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_FLUORINE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_FLUORINE line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithFluorine >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithFluorineCreator > reg_ReplaceProtonWithFluorineCreator;

///////////////

class ReplaceProtonWithBromineCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_BROMINE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_BROMINE line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithBromine >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithBromineCreator > reg_ReplaceProtonWithBromineCreator;

///////////////

class ReplaceProtonWithIodineCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_IODINE"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_IODINE line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithIodine >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithIodineCreator > reg_ReplaceProtonWithIodineCreator;

///////////////

class ReplaceProtonWithMethylCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_METHYL"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_METHYL line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithMethyl >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithMethylCreator > reg_ReplaceProtonWithMethylCreator;

///////////////

class ReplaceProtonWithTrifluoromethylCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_TRIFLUOROMETHYL"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_TRIFLUOROMETHYL line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithTrifluoromethyl >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithTrifluoromethylCreator > reg_ReplaceProtonWithTrifluoromethylCreator;

///////////////

class ReplaceProtonWithHydroxylCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_HYDROXYL"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_HYDROXYL line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithHydroxyl >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithHydroxylCreator > reg_ReplaceProtonWithHydroxylCreator;

///////////////

class ReplaceProtonWithMethoxyCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_METHOXY"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_METHOXY line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithMethoxy >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithMethoxyCreator > reg_ReplaceProtonWithMethoxyCreator;

///////////////

class ReplaceProtonWithEthylCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "REPLACE_PROTON_WITH_ETHYL"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad REPLACE_PROTON_WITH_ETHYL line in patchfile: " + line );
		return utility::pointer::make_shared< ReplaceProtonWithEthyl >( atom1 );
	}
};

static PatchOperationRegistrator< ReplaceProtonWithEthylCreator > reg_ReplaceProtonWithEthylCreator;

///////////////

class AddConnectAndDeleteChildProtonCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_CONNECT_AND_DELETE_CHILD_PROTON"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad ADD_CONNECT_AND_DELETE_CHILD_PROTON line in patchfile: " + line );
		return utility::pointer::make_shared< AddConnectDeleteChildProton >( atom1 );
	}
};

static PatchOperationRegistrator< AddConnectAndDeleteChildProtonCreator > reg_AddConnectAndDeleteChildProtonCreator;

///////////////

class DeleteChildProtonCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "DELETE_CHILD_PROTON"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad DELETE_CHILD_PROTON line in patchfile: " + line );
		return utility::pointer::make_shared< DeleteChildProton >( atom1 );
	}
};

static PatchOperationRegistrator< DeleteChildProtonCreator > reg_DeleteChildProtonCreator;

///////////////

class AddConnectAndTrackingVirtCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "ADD_CONNECT_AND_TRACKING_VIRT"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		std::string atom1;
		l >> atom1;
		if ( l.fail() ) utility_exit_with_message( "bad ADD_CONNECT_AND_TRACKING_VIRT line in patchfile: " + line );
		return utility::pointer::make_shared< AddConnectAndTrackingVirt >( atom1 );
	}
};

static PatchOperationRegistrator< AddConnectAndTrackingVirtCreator > reg_AddConnectAndTrackingVirtCreator;

///////////////

class ChiralFlipNamingCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "CHIRAL_FLIP_NAMING"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< ChiralFlipNaming >();//( atom1, atom2 ) );
	}
};

static PatchOperationRegistrator< ChiralFlipNamingCreator > reg_ChiralFlipNamingCreator;

///////////////

class ChiralFlipAtomsCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "CHIRAL_FLIP_ATOMS"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		//std::string atom1, atom2;
		//l >> atom1 >> atom2;
		return utility::pointer::make_shared< ChiralFlipAtoms >();//( atom1, atom2 ) );
	}
};

static PatchOperationRegistrator< ChiralFlipAtomsCreator > reg_ChiralFlipAtomsCreator;

///////////////

class VirtualizeSidechainCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "VIRTUALIZE_SIDECHAIN"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< VirtualizeSidechain >();//( atom1, atom2 ) );
	}
};

static PatchOperationRegistrator< VirtualizeSidechainCreator > reg_VirtualizeSidechainCreator;

///////////////

class VirtualizeAllCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "VIRTUALIZE_ALL"; }
	PatchOperationOP create_operation(std::string const &, std::istream &, std::map< std::string, Real > const & ) const override {
		return utility::pointer::make_shared< VirtualizeAll >();//( atom1, atom2 ) );
	}
};

static PatchOperationRegistrator< VirtualizeAllCreator > reg_VirtualizeAllCreator;

///////////////

class SetVirtualShadowCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "SET_VIRTUAL_SHADOW"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		runtime_assert( l.good() );
		std::string shadower, shadowee;
		l >> shadower >> shadowee;
		if ( l.fail() ) utility_exit_with_message( "bad SET_VIRTUAL_SHADOW line in patchfile: " + line );
		return utility::pointer::make_shared< SetVirtualShadow >( shadower, shadowee );
	}
};

static PatchOperationRegistrator< SetVirtualShadowCreator > reg_SetVirtualShadowCreator;

///////////////

class RenameAtomCreator: public PatchOperationCreator
{
public:
	std::string keyname() const override { return "RENAME_ATOM"; }
	PatchOperationOP create_operation(std::string const & line, std::istream & l, std::map< std::string, Real > const & ) const override {
		runtime_assert( l.good() );
		//std::string old_name, new_name;
		//l >> old_name >> new_name;
		std::string old_name = line.substr( 12, 4 ); // Rosetta atom
		std::string new_name = line.substr( 17, 4 );

		if ( l.fail() ) utility_exit_with_message( "bad RENAME_ATOM line in patchfile: " + line );
		return utility::pointer::make_shared< RenameAtom >( old_name, new_name );
	}
};

static PatchOperationRegistrator< RenameAtomCreator > reg_RenameAtomCreator;

///////////////

} //namespace chemical
} //namespace core
