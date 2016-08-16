// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/PatchOperation.hh
/// @brief  Polymorphic classes representing the contents of a residue-type patch file
/// @author Phil Bradley

#ifndef INCLUDED_core_chemical_PatchOperation_hh
#define INCLUDED_core_chemical_PatchOperation_hh

// Unit headers
#include <core/chemical/PatchOperation.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <map>

namespace core {
namespace chemical {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  A single operation that needs to be applied in a residue patch
class PatchOperation : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PatchOperation();

	/// @brief Returns true to signal failure, false to indicate success.
	virtual
	bool
	apply( ResidueType & rsd ) const = 0;

	/// @brief Which atom, if any, is added. Used for fast matching of ResidueType/Patches to PDB residues.
	virtual
	std::string
	adds_atom(){ return ""; }

	/// @brief Which atom, if any, is deleted. Used for fast matching of ResidueType/Patches to PDB residues.
	virtual
	std::string
	deletes_atom(){ return ""; }

	/// @brief Which property, if any, is added.
	virtual
	std::string
	adds_property(){ return ""; }

	/// @brief Which property, if any, is deleted.
	virtual
	std::string
	deletes_property(){ return ""; }

	/// @brief Which variant, if any, is deleted.
	virtual
	std::string
	deletes_variant(){ return ""; }

	/// @brief Generates a new aa
	virtual
	bool
	may_change_aa(){ return false; }

	/// @brief Generates name3.
	virtual
	std::string
	generates_name3(){ return ""; }

	/// @brief Generates interchangeability_group.
	virtual
	std::string
	generates_interchangeability_group(){ return ""; }

	/// @brief Generates base residue -- legacy for D_AA -- do not use otherwise.
	virtual
	bool
	generates_base_residue(){ return false; }

	/// @brief Special -- does this apply to 'minimal', placeholder types? Generally true, unless updating aa or name3.
	virtual
	bool
	applies_to_placeholder() const { return false; }

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete an atom
class DeleteAtom : public PatchOperation {
public:

	/// @brief constructor
	DeleteAtom( std::string const & atom_name_in );

	/// @brief delete an atom from ResidueType rsd
	bool
	apply( ResidueType & rsd ) const;

	virtual
	std::string
	deletes_atom(){ return atom_name_; }

private:
	/// name of the atom to be deleted
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as backbone heavy atom
class SetBackboneHeavyatom : public PatchOperation {
public:
	/// @brief constructor
	SetBackboneHeavyatom( std::string const & atom_name_in );

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( ResidueType & rsd ) const;

private:
	// name of the atom to be set
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as polymer connection
class SetPolymerConnectAtom : public PatchOperation {
public:

	/// @brief constructor the type of connection is "LOWER" or "UPPER"
	SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in );

	/// @brief set an atom in ResidueType rsd as a polymer connection atom
	bool
	apply( ResidueType & rsd ) const;

private:
	// "NONE" to delete the connection by setting its atomno to ZERO
	std::string atom_name_;
	// -1 for lower connection, 1 for upper connection
	int upper_lower_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AddConnect : public PatchOperation {
public:

	// c-tor
	AddConnect( std::string const & connect_atom,
		Real const phi, Real const theta, Real const d,
		std::string const & parent_atom,
		std::string const & angle_atom,
		std::string const & torsion_atom
	);

	/// add a property
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string const connect_atom_;
	Real const phi_;
	Real const theta_;
	Real const d_;
	std::string const  parent_atom_;
	std::string const   angle_atom_;
	std::string const torsion_atom_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a property to ResidueType
class AddProperty : public PatchOperation {
public:

	/// @brief constructor
	AddProperty( std::string const & property_in );

	/// @brief add a property
	bool
	apply( ResidueType & rsd ) const;

	/// @brief Which property, if any, is added.
	virtual
	std::string
	adds_property(){ return property_; }

private:
	/// property to be added
	std::string property_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete a property from ResidueType
///    Added by Andy M. Chen in June 2009
///    This is needed for deleting properties, which occurs in certain PTM's (e.g. methylation)
class DeleteProperty : public PatchOperation {
public:

	/// @brief constructor
	DeleteProperty( std::string const & property_in );

	/// @brief delete a property
	bool
	apply( ResidueType & rsd ) const;

	/// @brief Which property, if any, is deleted.
	virtual
	std::string
	deletes_property(){ return property_; }

private:
	// property to be added
	std::string property_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for deleting a VariantType from a ResidueType.
/// @author  Labonte <JWLabonte@jhu.edu>
class DeleteVariantType : public PatchOperation {
public:
	// Constructor
	DeleteVariantType( std::string const & variant_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

	/// @brief Which variant, if any, is deleted.
	virtual
	std::string
	deletes_variant(){ return variant_; }

private:
	std::string variant_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief   Add a chi angle to ResidueType.
/// @author  Added by Andy M. Chen in June 2009
/// @details This is needed for PTMs, which often result in one or more extra chi angles.
class AddChi : public PatchOperation {
public:
	/// @brief Constructor for when the chi index is specified.
	AddChi(Size const & chino_in,
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & atom3_in,
		std::string const & atom4_in);

	/// @brief Constructor for when the chi index is not specified.
	AddChi(std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & atom3_in,
		std::string const & atom4_in);

	/// @brief Add a chi angle.
	bool apply(ResidueType & rsd) const;

private:
	bool no_index_;  // indicates that no chi index is provided, and the new chi will be added to the end of the list
	Size chino_;

	// atoms defining the added chi angle
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};


class AddProtonChi : public PatchOperation {
public:

	/// constructor
	AddProtonChi(
		Size const & chino_in,
		utility::vector1<core::Real> const & samples,
		utility::vector1<core::Real> const & extrasamples
	);

	/// add a proton chi angle
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	utility::vector1<core::Real> samples_;
	utility::vector1<core::Real> extrasamples_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Redefine a chi angle
///    Added by Andy M. Chen in June 2009
///    This is needed for certain PTMs
class RedefineChi : public PatchOperation {
public:
	/// constructor
	RedefineChi(
		Size const & chino_in,
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & atom3_in,
		std::string const & atom4_in
	);

	/// redefine a chi angle
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete terminal chi angle
///    Added by Andrew M. Watkins in April 2015
class DeleteTerminalChi : public PatchOperation {
public:
	/// constructor
	DeleteTerminalChi() {};

	/// redefine a chi angle
	bool
	apply( ResidueType & rsd ) const;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete a metal binding atom
///    Added by Andrew M. Watkins in April 2015
class DeleteMetalbindingAtom : public PatchOperation {
public:
	/// constructor
	DeleteMetalbindingAtom( std::string const & atom_name );

	bool
	apply( ResidueType & rsd ) const;
private:
	std::string atom_name_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete an act coord atom
///    Added by Andrew M. Watkins in April 2015
class DeleteActCoordAtom : public PatchOperation {
public:
	/// constructor
	DeleteActCoordAtom( std::string const & atom_name );

	bool
	apply( ResidueType & rsd ) const;
private:
	std::string atom_name_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief   Add a rotamer sample to a chi angle of the ResidueType.
/// @author  Added by Andy M. Chen in June 2009
/// @details This is needed for PTMs.
class AddChiRotamer : public PatchOperation {
public:
	/// @brief Constructor for when the chi index is specified
	AddChiRotamer(Size const & chino_in, Angle const & mean_in, Angle const & sdev_in);

	/// @brief Constructor for when the chi index is not specified
	AddChiRotamer(Angle const & mean_in, Angle const & sdev_in);

	/// @brief Add a rotamer sample.
	bool apply(ResidueType & rsd) const;

private:
	bool no_index_;  // indicates that no chi index is provided, and the rotamer sample will be added to the last chi
	Size chino_;
	Real mean_;
	Real sdev_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for clearing all rotamer bins from the chi of a ResidueType.
/// @note    This is useful if one has redefined a chi.
/// @author  Labonte <JWLabonte@jhu.edu>
class ClearChiRotamers : public PatchOperation {
public:
	// Constructor
	ClearChiRotamers( core::uint const chi_no_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

private:
	core::uint chi_no_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add an atom to ResidueType
#if defined(WIN32) && !defined(WIN_PYROSETTA)
class AddAtomWIN32 : public PatchOperation {
#else
class AddAtom : public PatchOperation {
#endif
public:

	/// constructor
#if defined(WIN32) && !defined(WIN_PYROSETTA)
	AddAtomWIN32(
#else
	AddAtom(
#endif
		std::string const & atom_name_in,
		std::string const & atom_type_name_in,
		std::string const & mm_atom_type_name_in,
		Real const charge
	);

	/// add an atom
	bool
	apply( ResidueType & rsd ) const;

	virtual
	std::string
	adds_atom(){ return atom_name_; }

private:
	std::string atom_name_;
	std::string atom_type_name_;
	std::string mm_atom_type_name_;
	Real const charge_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for adding an atom alias to a ResidueType.
/// @note    See residue_io.cc for a description of atom aliases.
/// @remarks Atom aliases were graciously added to Rosetta by Rocco Moretti.
/// @author  Labonte <JWLabonte@jhu.edu>
class AddAtomAlias : public PatchOperation {
public:
	// Constructor
	AddAtomAlias( std::string const & rosetta_atom_name_in, std::string const & alias_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

	virtual
	std::string
	adds_atom(){ return alias_; }

private:
	std::string rosetta_atom_name_;
	std::string alias_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a bond to ResidueType
class AddBond : public PatchOperation {
public:

	/// constructor
	AddBond(
		std::string const & atom1_in,
		std::string const & atom2_in
	);

	/// add a bond
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a bond is added
	std::string atom1_;
	std::string atom2_;

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for adding a specific type of bond to a ResidueType.
/// @note    See residue_io.cc for a description of bond types.
/// @author  Labonte <JWLabonte@jhu.edu>
class AddBondType : public PatchOperation {
public:
	// Constructor
	AddBondType(
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & bond_type_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

private:
	std::string atom1_;
	std::string atom2_;
	std::string bond_type_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for changing the bond type of a given bond.
/// @note    See residue_io.cc for a description of bond types.
/// @author  Labonte <JWLabonte@jhu.edu>
class ChangeBondType : public PatchOperation {
public:
	// Constructor
	ChangeBondType(
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & old_bond_type_in,
		std::string const & new_bond_type_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( ResidueType & rsd ) const;

private:
	std::string atom1_;
	std::string atom2_;
	std::string old_bond_type_;
	std::string new_bond_type_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom's charge
class SetAtomicCharge : public PatchOperation {
public:

	/// constructor
	SetAtomicCharge(
		std::string const & atom_name_in,
		Real const charge_in
	);

	/// set an atom's charge
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's charge
	Real charge_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for setting the formal charge of a ResidueType's atom.
/// @author  Labonte <JWLabonte@jhu.edu>
class SetFormalCharge : public PatchOperation {
public:
	// Constructor
	SetFormalCharge( std::string const & atom_name_in, core::SSize charge_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
	core::SSize charge_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's chemical type
class SetAtomType : public PatchOperation {
public:
	/// constructor
	SetAtomType(
		std::string const & atom_name_in,
		std::string const & atom_type_name_in
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string atom_type_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set residue's name1 and name3
class SetIO_String : public PatchOperation {
public:
	/// constructor
	SetIO_String(
		std::string const & name3,
		char const name1
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

	/// @brief Generates name3.
	virtual
	std::string
	generates_name3(){ return name3_; }

	virtual
	bool
	applies_to_placeholder() const { return true; }

private:
	std::string const name3_;
	char const name1_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the interchangeability_group string for a ResidueType
class SetInterchangeabilityGroup_String : public PatchOperation {
public:
	/// constructor
	SetInterchangeabilityGroup_String(
		std::string const & intgrp
	);

	bool
	apply( ResidueType & rsd ) const;

	std::string
	generates_interchangeability_group(){ return intgrp_; }

private:
	std::string const intgrp_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's MM chemical type
class SetMMAtomType : public PatchOperation {
public:
	/// constructor
	SetMMAtomType(
		std::string const & atom_name_in,
		std::string const & mm_atom_type_name_in
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string mm_atom_type_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom's AtomICoord
class SetICoor : public PatchOperation {
public:
	/// constructor
	SetICoor(
		std::string const & atom_in,
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub1_in,
		std::string const & stub2_in,
		std::string const & stub3_in
	);

	/// set an atom's AtomICoord
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_;
	/// AtomICoord data
	Real phi_;
	Real theta_;
	Real d_;
	std::string stub1_;
	std::string stub2_;
	std::string stub3_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Change the parent, grandparent, or great-grandparent of an atom
enum Ancestor { anc_parent, anc_grandparent, anc_greatgrandparent };
class ChangeAncestory : public PatchOperation {
public:
	/// constructor
	ChangeAncestory(
		std::string const & target_atom,
		Ancestor which_ancestor,
		std::string const & ancestor_name
	);

	/// @brief change the ancestory, but leave the icoors intact.
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_;
	Ancestor which_ancestor_;
	std::string ancestor_name_;
};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for resetting the length of a bond within a ResidueType.
/// @note    This is useful for when an atom is rehybridized within a patch file.
/// @author  Labonte <JWLabonte@jhu.edu>
class ResetBondLength : public PatchOperation {
public:
	// Constructor
	ResetBondLength( std::string const & atm_in, core::Distance d_in );

	/// @brief  Apply this patch to the given ResidueType.
	virtual bool apply( ResidueType & rsd ) const;

private:
	std::string atm_;
	core::Distance d_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom before the first mainchain atom
class PrependMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	PrependMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the first mainchain atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom after the last mainchain atom
class AppendMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	AppendMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace a mainchain atom
class ReplaceMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	ReplaceMainchainAtom( std::string const & target, std::string const & new_atom );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string target_;
	std::string new_atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor atom
class SetNbrAtom : public PatchOperation {
public:
	/// @brief constructor
	SetNbrAtom( std::string const & atom_name_in );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetNbrRadius : public PatchOperation {
public:
	/// @brief constructor
	SetNbrRadius( Real const & radius );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const;

private:
	Real radius_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetAllAtomsRepulsive : public PatchOperation {
public:
	/// @brief constructor
	SetAllAtomsRepulsive() {};

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set orient atom selection mode.
class SetOrientAtom : public PatchOperation {
public:
	SetOrientAtom(bool force_nbr_atom_orient);

	bool
	apply( ResidueType & rsd ) const;

private:
	bool force_nbr_atom_orient_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the path to a rotamer library for an NCAA that is not in dunbrack
class NCAARotLibPath : public PatchOperation {
public:
	/// @brief constructor
	NCAARotLibPath( std::string const & path_in );

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string path_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
class ConnectSulfurAndMakeVirtualProton : public PatchOperation {
public:
	/// @brief constructor
	ConnectSulfurAndMakeVirtualProton() {};

	bool
	apply( ResidueType & rsd ) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipNaming : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipNaming() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};

	virtual
	bool
	applies_to_placeholder() const { return true; }

	/// @brief Generates a new aa
	virtual
	bool
	may_change_aa(){ return true; }

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const;

private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipAtoms : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipAtoms() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};
	/// @brief set the NCAA rotamer library path in the residue type

	bool
	apply( ResidueType & rsd ) const;

private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithMethoxy: public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethoxy( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithEthyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithEthyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithMethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with trifluoromethyl
class ReplaceProtonWithTrifluoromethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithTrifluoromethyl( std::string atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with chlorine
class ReplaceProtonWithChlorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithChlorine( std::string atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with fluorine
class ReplaceProtonWithFluorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithFluorine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with bromine
class ReplaceProtonWithBromine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithBromine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with iodine
class ReplaceProtonWithIodine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithIodine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with hydroxyl
class ReplaceProtonWithHydroxyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithHydroxyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect and tracking virt to the atom
class AddConnectAndTrackingVirt : public PatchOperation {
public:
	/// @brief constructor
	AddConnectAndTrackingVirt( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect to the atom, delete child proton
class AddConnectDeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	AddConnectDeleteChildProton( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete child proton
class DeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	DeleteChildProton( std::string const & atom ): atom_( atom ) {};

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_;
};

/// @brief  Virtual constructor, returns 0 if no match
PatchOperationOP
patch_operation_from_patch_file_line(
	std::string const & line,
	std::map< std::string, Real > const & atomic_charge_reassignments
);

} // chemical
} // core

#endif
