// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Phil Bradley
/// see patch.cc to understand whats going on

#ifndef INCLUDED_core_chemical_PatchOperation_hh
#define INCLUDED_core_chemical_PatchOperation_hh


// // Unit headers
#include <core/chemical/PatchOperation.fwd.hh>

// // Package headers
#include <core/chemical/ResidueType.hh>

//Tracer header
#include <basic/Tracer.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>




namespace core {
namespace chemical {

static basic::Tracer TR_PatchOperations("core.chemical.PatchOperations.hh");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief  A single operation that needs to be applied in a residue patch
class PatchOperation : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PatchOperation();

	/// @brief  Returns TR_PatchOperationsUE to signal failure
	virtual
	bool
	apply( ResidueType & rsd ) const = 0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete an atom
class DeleteAtom : public PatchOperation {
public:

	/// @brief constructor
	DeleteAtom( std::string const & atom_name_in ) :
		atom_name_( atom_name_in )
	{}

	/// @brief delete an atom from ResidueType rsd
	bool
	apply( ResidueType & rsd ) const
	{

		if ( !rsd.has( atom_name_ )  ) {
			TR_PatchOperations.Debug << "DeleteAtom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ << std::endl;
			return true; // failure
		} else {
			//std::cout << "DeleteAtom::apply: deleting: " << atom_name_ << std::endl;
			rsd.delete_atom( atom_name_ );
		}
		return false;
	}

private:
	/// name of the atom to be deleted
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as backbone heavy atom
class SetBackboneHeavyatom : public PatchOperation {
public:
	/// @brief constructor
	SetBackboneHeavyatom( std::string const & atom_name_in ) :
		atom_name_( atom_name_in )
	{}

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "SetBackboneHeavyatom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ <<
				std::endl;
			return true; // failure
		} else {
			//std::cout << "SetBackboneHeavyatom::apply: " << atom_name_ << std::endl;
			rsd.set_backbone_heavyatom( atom_name_ );
		}
		return false;
	}

private:
	// name of the atom to be set
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as polymer connection
class SetPolymerConnectAtom : public PatchOperation {
public:

	/// @brief constructor the type of connection is "LOWER" or "UPPER"
	SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in ) :
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

	/// @brief set an atom in ResidueType rsd as a polymer connection atom
	bool
	apply( ResidueType & rsd ) const
	{
		if ( atom_name_ == "NONE" || rsd.has( atom_name_ ) ) {
			//std::cout << "SetPolymerConnectAtom::apply: " << atom_name_ << ' ' << upper_lower_ << std::endl;
			if ( upper_lower_ == -1 ) {
				rsd.set_lower_connect_atom( atom_name_ );
			} else {
				assert( upper_lower_ == 1 );
				rsd.set_upper_connect_atom( atom_name_ );
			}
		} else {
			TR_PatchOperations.Debug << "SetPolymerConnectAtom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ <<
				std::endl;
			return true; // failure
		}
		return false;
	}

private:
	/// "NONE" to delete the connection by setting its atomno to ZERO
	std::string atom_name_;
	/// -1 for lower connection, 1 for upper connection
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
							):
		connect_atom_( connect_atom ),
		phi_( phi ), theta_( theta ), d_( d ),
		parent_atom_ (  parent_atom ),
		angle_atom_  (   angle_atom ),
		torsion_atom_( torsion_atom )
	{}


	/// add a property
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( connect_atom_ ) ||
				 !rsd.has(  parent_atom_ ) ||
				 !rsd.has(   angle_atom_ ) ||
				 !rsd.has( torsion_atom_ ) ) return true; // failure!

		Size const connid( rsd.add_residue_connection( connect_atom_ ) );
		rsd.set_icoor( "CONN"+ObjexxFCL::string_of( connid ), phi_, theta_, d_, parent_atom_, angle_atom_, torsion_atom_ );
		return false;
	}

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

	/// constructor
	AddProperty( std::string const & property_in ):
		property_( property_in )
	{}

	/// add a property
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.add_property( property_ );
		//std::cout << "AddProperty::apply: " << property_ << std::endl;
		return false;
	}

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

	/// constructor
	DeleteProperty( std::string const & property_in ):
		property_( property_in )
	{}

	/// add a property
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.delete_property( property_ );
		//std::cout << "DeleteProperty::apply: " << property_ << std::endl;
		return false;
	}

private:
	/// property to be added
	std::string property_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Add a chi angle to ResidueType
///    Added by Andy M. Chen in June 2009
///    This is needed for PTM's, which often result in one or more extra chi angles

class AddChi : public PatchOperation {
public:

	/// constructor
	AddChi(
		Size const & chino_in, std::string const & atom1_in, std::string const & atom2_in,
		std::string const & atom3_in, std::string const & atom4_in
	):
		chino_( chino_in ), atom1_( atom1_in ), atom2_( atom2_in ), atom3_( atom3_in ), atom4_( atom4_in )
	{}

	/// add a chi angle
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) || !rsd.has( atom3_ ) || !rsd.has( atom4_ ) )
		{
			TR_PatchOperations.Debug << "AddChi::apply failed: " << rsd.name() << " is missing atom(s) " << atom1_ << ' '
				<< rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) << atom3_ << ' '
				<< rsd.has( atom3_ ) << atom4_ << ' ' << rsd.has( atom4_ ) << std::endl;
			return true; // failure
		}
		else
		{
			rsd.add_chi( chino_, atom1_ , atom2_, atom3_, atom4_ );
			//std::cout << "AddChi::apply: " << atom1_ << ' ' << atom2_
			//	<< ' ' << atom3_ << ' ' << atom4_ << std::endl;
			return false;
		}
		return false;
	}


private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Redefine a chi angle
///    Added by Andy M. Chen in June 2009
///    This is needed for certain PTM's

class RedefineChi : public PatchOperation {
public:

	/// constructor
	RedefineChi(
		Size const & chino_in, std::string const & atom1_in, std::string const & atom2_in,
		std::string const & atom3_in, std::string const & atom4_in
	):
		chino_( chino_in ), atom1_( atom1_in ), atom2_( atom2_in ), atom3_( atom3_in ), atom4_( atom4_in )
	{}

	/// redefine a chi angle
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) || !rsd.has( atom3_ ) || !rsd.has( atom4_ ) )
		{
			TR_PatchOperations.Debug << "RedefineChi::apply failed: " << rsd.name() << " is missing atom(s) " << atom1_ << ' '
				<< rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) << atom3_ << ' '
				<< rsd.has( atom3_ ) << atom4_ << ' ' << rsd.has( atom4_ ) << std::endl;
			return true; // failure
		}
		else
		{
			rsd.redefine_chi( chino_, atom1_ , atom2_, atom3_, atom4_ );
			//std::cout << "RedefineChi::apply: " << atom1_ << ' ' << atom2_
			//	<< ' ' << atom3_ << ' ' << atom4_ << std::endl;
			return false;
		}
		return false;
	}


private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Add a rotamer sample to a chi angle of the ResidueType
///    Added by Andy M. Chen in June 2009
///    This is needed for PTM's

class AddChiRotamer : public PatchOperation {
public:

	/// constructor
	AddChiRotamer(
		Size const & chino_in, Real const & mean_in, Real const & sdev_in
	):
		chino_( chino_in ), mean_( mean_in ), sdev_( sdev_in )
	{}

	/// add a rotamer sample
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.add_chi_rotamer( chino_, mean_ , sdev_ );
		//std::cout << "AddChiRotamer::apply: " << chino_ << ' ' << mean_
		//	<< ' ' << sdev_ << std::endl;
		return false;
	}


private:
	/// atoms between which a chi angle is added
	Size chino_;
	Real mean_;
	Real sdev_;
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
	):
		atom_name_( atom_name_in ),
		atom_type_name_( atom_type_name_in ),
		mm_atom_type_name_( mm_atom_type_name_in ),
		charge_( charge )
	{}

	/// add an atom
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.add_atom( atom_name_, atom_type_name_, mm_atom_type_name_, charge_ );
		//std::cout << "AddAtom::apply: " << atom_name_ << ' ' << atom_type_name_ << ' ' << mm_atom_type_name_ << ' ' <<
		//	charge_ << std::endl;
		return false;
	}

private:
	std::string atom_name_;
	std::string atom_type_name_;
	std::string mm_atom_type_name_;
	Real const charge_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a bond to ResidueType
class AddBond : public PatchOperation {
public:

	/// constructor
	AddBond(
		std::string const & atom1_in,
		std::string const & atom2_in
	):
		atom1_( atom1_in ),
		atom2_( atom2_in )
	{}

	/// add a bond
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom1_ ) || !rsd.has( atom2_ ) ) {
			TR_PatchOperations.Debug << "AddBond::apply failed: " << rsd.name() << " is missing atom(s) " << atom1_ << ' ' <<
				rsd.has( atom1_ ) << ' ' << atom2_ << ' ' << rsd.has( atom2_ ) << std::endl;
			return true; // failure
		} else {
			//std::cout << "AddBond::apply: " << atom1_ << ' ' << atom2_ << std::endl;
			rsd.add_bond( atom1_, atom2_ );
		}
		return false;
	}

private:
	/// atoms between which a bond is added
	std::string atom1_;
	std::string atom2_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom's charge
class SetAtomicCharge : public PatchOperation {
public:

	/// constructor
	SetAtomicCharge(
		std::string const & atom_name_in,
		Real const charge_in
	):
		atom_name_( atom_name_in ),
		charge_( charge_in )
	{}

	/// set an atom's charge
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's charge
	Real charge_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's chemical type
class SetAtomType : public PatchOperation {
public:
	/// constructor
	SetAtomType(
		std::string const & atom_name_in,
		std::string const & atom_type_name_in
	):
		atom_name_( atom_name_in ),
		atom_type_name_( atom_type_name_in )
	{}

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "SetAtomType::apply failed: " << rsd.name() << " is missing atom: " << atom_name_ << std::endl;
			return true; // failure
		} else {
			//std::cout << "SetAtomType::apply: " << atom_name_ << ' ' << atom_type_name_ << std::endl;
			rsd.set_atom_type( atom_name_, atom_type_name_ );
		}
		return false;
	}

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string atom_type_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's chemical type
class SetIO_String : public PatchOperation {
public:
	/// constructor
	SetIO_String(
		std::string const & name3,
		char const name1
	):
		name3_( name3 ),
		name1_( name1 )
	{}

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.name3( name3_ );
		rsd.name1( name1_ );
		return false;
	}

private:
	std::string const name3_;
	char const name1_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's MM chemical type
class SetMMAtomType : public PatchOperation {
public:
	/// constructor
	SetMMAtomType(
		std::string const & atom_name_in,
		std::string const & mm_atom_type_name_in
	):
		atom_name_( atom_name_in ),
		mm_atom_type_name_( mm_atom_type_name_in )
	{}

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "SetAtomType::apply failed: " << rsd.name() << " is missing atom: " << atom_name_ << std::endl;
			return true; // failure
		} else {
			//std::cout << "SetAtomType::apply: " << atom_name_ << ' ' << mm_atom_type_name_ << std::endl;
			rsd.set_mm_atom_type( atom_name_, mm_atom_type_name_ );
		}
		return false;
	}

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
	):
		atom_( atom_in ),
		phi_( phi_in ),
		theta_( theta_in ),
		d_( d_in ),
		stub1_( stub1_in ),
		stub2_( stub2_in ),
		stub3_( stub3_in )
	{}

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
/// @brief add a mainchain atom before the first mainchain atom
class PrependMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	PrependMainchainAtom( std::string const & atom_name_in ) :
		atom_name_( atom_name_in )
	{}

	/// @brief set an atom to be the first mainchain atom
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "PrependMainchainAtom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ <<	std::endl;
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

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom after the last mainchain atom
class AppendMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	AppendMainchainAtom( std::string const & atom_name_in ) :
		atom_name_( atom_name_in )
	{}

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "AppendMainchainAtom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ <<	std::endl;
			return true; // failure
		} else {
			AtomIndices new_mainchain_atoms( rsd.mainchain_atoms() );
			new_mainchain_atoms.push_back( rsd.atom_index( atom_name_ ) );
			rsd.set_mainchain_atoms( new_mainchain_atoms );
		}
		return false;
	}

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor atom
class SetNbrAtom : public PatchOperation {
public:
	/// @brief constructor
	SetNbrAtom( std::string const & atom_name_in ) :
		atom_name_( atom_name_in )
	{}

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const
	{
		if ( !rsd.has( atom_name_ ) ) {
			TR_PatchOperations.Debug << "SetNbrAtom::apply failed: " << rsd.name() << " is missing atom " << atom_name_ <<	std::endl;
			return true; // failure
		} else {
			rsd.nbr_atom( atom_name_ );
		}
		return false;
	}

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetNbrRadius : public PatchOperation {
public:
	/// @brief constructor
	SetNbrRadius( Real const & radius ) :
		radius_( radius )
	{}

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.nbr_radius( radius_ );
		return false;
	}

private:
	Real radius_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set orient atom selection mode.
class SetOrientAtom : public PatchOperation {
	public:
		SetOrientAtom(bool force_nbr_atom_orient):
			force_nbr_atom_orient_(force_nbr_atom_orient)
		{}

		bool
		apply( ResidueType & rsd ) const
		{
			rsd.force_nbr_atom_orient( force_nbr_atom_orient_ );
			return false;
		}

	private:
		bool force_nbr_atom_orient_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the path to a rotamer library for an NCAA that is not in dunbrack
class NCAARotLibPath : public PatchOperation {
public:
	/// @brief constructor
	NCAARotLibPath( std::string const & path_in ) :
		path_( path_in )
	{}

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const
	{
		rsd.set_ncaa_rotlib_path( path_ );
		rsd.set_use_ncaa_rotlib( true );
		return false;
	}

private:
	std::string path_;
};

/// @brief  Virtual constructor, returns 0 if no match
PatchOperationOP
patch_operation_from_patch_file_line( std::string const & line );



} // chemical
} // core



#endif
