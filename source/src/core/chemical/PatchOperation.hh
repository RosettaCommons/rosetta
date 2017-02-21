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
#include <utility>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  A single operation that needs to be applied in a residue patch
class PatchOperation : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~PatchOperation() override;

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete an atom
class DeleteAtom : public PatchOperation {
public:

	/// @brief constructor
	DeleteAtom( std::string const & atom_name_in );

	/// @brief delete an atom from ResidueType rsd
	bool
	apply( ResidueType & rsd ) const override;


	std::string
	deletes_atom() override{ return atom_name_; }

private:
	/// name of the atom to be deleted
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as backbone heavy atom
class SetBackboneHeavyatom : public PatchOperation {
public:
	/// @brief constructor
	SetBackboneHeavyatom( std::string const & atom_name_in );

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	// name of the atom to be set
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetBackboneHeavyatom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as polymer connection
class SetPolymerConnectAtom : public PatchOperation {
public:

	/// @brief constructor the type of connection is "LOWER" or "UPPER"
	SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in );

	/// @brief set an atom in ResidueType rsd as a polymer connection atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	// "NONE" to delete the connection by setting its atomno to ZERO
	std::string atom_name_;
	// -1 for lower connection, 1 for upper connection
	int upper_lower_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetPolymerConnectAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	std::string connect_atom_;
	Real phi_;
	Real theta_;
	Real d_;
	std::string  parent_atom_;
	std::string   angle_atom_;
	std::string torsion_atom_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnect();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a property to ResidueType
class AddProperty : public PatchOperation {
public:

	/// @brief constructor
	AddProperty( std::string const & property_in );

	/// @brief add a property
	bool
	apply( ResidueType & rsd ) const override;

	/// @brief Which property, if any, is added.

	std::string
	adds_property() override{ return property_; }

private:
	/// property to be added
	std::string property_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddProperty();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

	/// @brief Which property, if any, is deleted.

	std::string
	deletes_property() override{ return property_; }

private:
	// property to be added
	std::string property_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteProperty();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for deleting a VariantType from a ResidueType.
/// @author  Labonte <JWLabonte@jhu.edu>
class DeleteVariantType : public PatchOperation {
public:
	// Constructor
	DeleteVariantType( std::string const & variant_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( ResidueType & rsd ) const override;

	/// @brief Which variant, if any, is deleted.

	std::string
	deletes_variant() override{ return variant_; }

private:
	std::string variant_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteVariantType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply(ResidueType & rsd) const override;

private:
	bool no_index_;  // indicates that no chi index is provided, and the new chi will be added to the end of the list
	Size chino_;

	// atoms defining the added chi angle
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	utility::vector1<core::Real> samples_;
	utility::vector1<core::Real> extrasamples_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddProtonChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	RedefineChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete a metal binding atom
///    Added by Andrew M. Watkins in April 2015
class DeleteMetalbindingAtom : public PatchOperation {
public:
	/// constructor
	DeleteMetalbindingAtom( std::string const & atom_name );

	bool
	apply( ResidueType & rsd ) const override;
private:
	std::string atom_name_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteMetalbindingAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete an act coord atom
///    Added by Andrew M. Watkins in April 2015
class DeleteActCoordAtom : public PatchOperation {
public:
	/// constructor
	DeleteActCoordAtom( std::string const & atom_name );

	bool
	apply( ResidueType & rsd ) const override;
private:
	std::string atom_name_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteActCoordAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply(ResidueType & rsd) const override;

private:
	bool no_index_;  // indicates that no chi index is provided, and the rotamer sample will be added to the last chi
	Size chino_;
	Real mean_;
	Real sdev_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddChiRotamer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply( ResidueType & rsd ) const override;

private:
	core::uint chi_no_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ClearChiRotamers();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add an atom to ResidueType
class AddAtom : public PatchOperation {
public:

	/// constructor
	AddAtom(
		std::string const & atom_name_in,
		std::string const & atom_type_name_in,
		std::string const & mm_atom_type_name_in,
		Real const charge
	);

	/// add an atom
	bool
	apply( ResidueType & rsd ) const override;


	std::string
	adds_atom() override{ return atom_name_; }

private:
	std::string atom_name_;
	std::string atom_type_name_;
	std::string mm_atom_type_name_;
	Real charge_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply( ResidueType & rsd ) const override;


	std::string
	adds_atom() override{ return alias_; }

private:
	std::string rosetta_atom_name_;
	std::string alias_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddAtomAlias();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atoms between which a bond is added
	std::string atom1_;
	std::string atom2_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddBond();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply( ResidueType & rsd ) const override;

private:
	std::string atom1_;
	std::string atom2_;
	std::string bond_type_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddBondType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply( ResidueType & rsd ) const override;

private:
	std::string atom1_;
	std::string atom2_;
	std::string old_bond_type_;
	std::string new_bond_type_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ChangeBondType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's charge
	Real charge_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetAtomicCharge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for setting the formal charge of a ResidueType's atom.
/// @author  Labonte <JWLabonte@jhu.edu>
class SetFormalCharge : public PatchOperation {
public:
	// Constructor
	SetFormalCharge( std::string const & atom_name_in, core::SSize charge_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( ResidueType & rsd ) const override;

private:
	std::string atom_name_;
	core::SSize charge_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetFormalCharge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string atom_type_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetAtomType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

	/// @brief Generates name3.

	std::string
	generates_name3() override{ return name3_; }


	bool
	applies_to_placeholder() const override { return true; }

private:
	std::string name3_;
	char name1_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetIO_String();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

	std::string
	generates_interchangeability_group() override{ return intgrp_; }

private:
	std::string intgrp_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetInterchangeabilityGroup_String();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string mm_atom_type_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetMMAtomType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

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
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetICoor();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( ResidueType & rsd ) const override;

private:
	/// atom's name
	std::string atom_;
	Ancestor which_ancestor_;
	std::string ancestor_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ChangeAncestory();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply( ResidueType & rsd ) const override;

private:
	std::string atm_;
	core::Distance d_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ResetBondLength();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom before the first mainchain atom
class PrependMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	PrependMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the first mainchain atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	PrependMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom after the last mainchain atom
class AppendMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	AppendMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AppendMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace a mainchain atom
class ReplaceMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	ReplaceMainchainAtom( std::string const & target, std::string const & new_atom );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string target_;
	std::string new_atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor atom
class SetNbrAtom : public PatchOperation {
public:
	/// @brief constructor
	SetNbrAtom( std::string const & atom_name_in );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetNbrAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetNbrRadius : public PatchOperation {
public:
	/// @brief constructor
	SetNbrRadius( Real const & radius );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const override;

private:
	Real radius_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetNbrRadius();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetAllAtomsRepulsive : public PatchOperation {
public:
	/// @brief constructor
	SetAllAtomsRepulsive() {};

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const override;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set orient atom selection mode.
class SetOrientAtom : public PatchOperation {
public:
	SetOrientAtom(bool force_nbr_atom_orient);

	bool
	apply( ResidueType & rsd ) const override;

private:
	bool force_nbr_atom_orient_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetOrientAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Remove existing rotamer specifications (of any type).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RemoveRotamerSpecifications : public PatchOperation {
public:
	/// @brief Constructor.
	///
	RemoveRotamerSpecifications();

	/// @brief Strip all RotamerSpecifications from the ResidueType.
	///
	bool
	apply( ResidueType & rsd ) const;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the filenames for RamaPrePro scoring tables.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RamaPreproFilename : public PatchOperation {
public:

	/// @brief Constructor.
	///
	RamaPreproFilename( std::string const &non_prepro_file, std::string const &prepro_file );

	/// @brief Set the RamaPrepro library paths in the residue type.
	///
	bool
	apply( ResidueType & rsd ) const override;

private:

	/// @brief The rama table file to use for scoring residues that do NOT precede proline.
	///
	std::string non_prepro_file_;

	/// @brief The rama table file to use for scoring residues that DO precede proline.
	///
	std::string prepro_file_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the residue name for RamaPrePro scoring tables.
/// @details This is the name in the scoring table AND the reference string used to look up the table.  Should be
/// unique.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RamaPreproResname : public PatchOperation {
public:

	/// @brief Constructor.
	///
	RamaPreproResname( std::string const &resname_in );

	/// @brief Set the RamaPrepro reference string in the residue type.
	///
	bool
	apply( ResidueType & rsd ) const override;

private:

	/// @brief The name listed for the amino acid type in the scoring table.
	/// @details Must be unique.
	std::string resname_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the path to a rotamer library for an NCAA that is not in dunbrack
class NCAARotLibPath : public PatchOperation {
public:
	/// @brief constructor
	NCAARotLibPath( std::string const & path_in );

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string path_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NCAARotLibPath();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the number of rotamer bins per chi for an NCAA that is not in dunbrack.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
class NCAARotLibNumRotamerBins : public PatchOperation {
public:
	/// @brief Constructor
	///
	NCAARotLibNumRotamerBins( utility::vector1< core::Size > const &binsizes_in );

	/// @brief Set the number of rotamer bins per chi for an NCAA.
	///
	bool
	apply( ResidueType & rsd ) const;

private:
	utility::vector1 < core::Size > binsizes_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
class ConnectSulfurAndMakeVirtualProton : public PatchOperation {
public:
	/// @brief constructor
	ConnectSulfurAndMakeVirtualProton() {};

	bool
	apply( ResidueType & rsd ) const override;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipNaming : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipNaming() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};


	bool
	applies_to_placeholder() const override { return true; }

	/// @brief Generates a new aa

	bool
	may_change_aa() override{ return true; }

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const override;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipAtoms : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipAtoms() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};
	/// @brief set the NCAA rotamer library path in the residue type

	bool
	apply( ResidueType & rsd ) const override;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithMethoxy: public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethoxy( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithMethoxy();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithEthyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithEthyl( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithEthyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithMethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethyl( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithMethyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with trifluoromethyl
class ReplaceProtonWithTrifluoromethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithTrifluoromethyl( std::string atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithTrifluoromethyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with chlorine
class ReplaceProtonWithChlorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithChlorine( std::string atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithChlorine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with fluorine
class ReplaceProtonWithFluorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithFluorine( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithFluorine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with bromine
class ReplaceProtonWithBromine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithBromine( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithBromine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with iodine
class ReplaceProtonWithIodine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithIodine( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithIodine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with hydroxyl
class ReplaceProtonWithHydroxyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithHydroxyl( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithHydroxyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect and tracking virt to the atom
class AddConnectAndTrackingVirt : public PatchOperation {
public:
	/// @brief constructor
	AddConnectAndTrackingVirt( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnectAndTrackingVirt();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect to the atom, delete child proton
class AddConnectDeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	AddConnectDeleteChildProton( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnectDeleteChildProton();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete child proton
class DeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	DeleteChildProton( std::string  atom ): atom_(std::move( atom )) {};

	bool
	apply( ResidueType & rsd ) const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteChildProton();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief virtualize all
class VirtualizeAll: public PatchOperation {
public:
	/// @brief constructor
	VirtualizeAll() {};

	bool
	apply( ResidueType & rsd ) const;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set virtual shadow atoms
class SetVirtualShadow : public PatchOperation {
public:
	/// constructor
	SetVirtualShadow(
		std::string const & shadower_,
		std::string const & shadowee_
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const override;

private:
	/// atom's name
	std::string shadower_;
	/// atom's type name
	std::string shadowee_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetVirtualShadow();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Virtual constructor, returns 0 if no match
PatchOperationOP
patch_operation_from_patch_file_line(
	std::string const & line,
	std::map< std::string, Real > const & atomic_charge_reassignments
);

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_PatchOperation )
#endif // SERIALIZATION


#endif
