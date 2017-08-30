// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/ResidueInformation.cc
/// @brief  Method definitions for ResidueInformation
/// @author Sergey Lyskov


// Unit headers
#include <core/io/ResidueInformation.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/NomenclatureManager.hh>

// Numeric headers
#include <numeric/xyzVector.hh>


namespace core {
namespace io {

ResidueInformation::ResidueInformation() :
	resName_( "" ),
	chainID_( ' ' ),
	resSeq_( 0 ),
	iCode_( ' ' ),
	terCount_( 0 ),
	atoms_(),
	xyz_(),
	temps_(),
	segmentID_( "    " ),
	rosetta_resName_( "" )
{}

ResidueInformation::ResidueInformation( AtomInformation const & ai ) :
	chainID_( ai.chainID ),
	resSeq_( ai.resSeq ),
	iCode_( ai.iCode ),
	terCount_( ai.terCount ),
	atoms_(),
	xyz_(),
	temps_(),
	segmentID_( ai.segmentID )
{
	resName( ai.resName ); // Needed to properly set name-dependent variables
}

bool
ResidueInformation::operator==( ResidueInformation const & that) const
{
	return
		resName_  == that.resName_ &&
		chainID_  == that.chainID_ &&
		resSeq_   == that.resSeq_  &&
		iCode_    == that.iCode_   &&
		terCount_ == that.terCount_;
}

bool
ResidueInformation::operator!=( ResidueInformation const & that) const
{
	return ! ( *this == that );
}

std::string const & ResidueInformation::resName() const { return resName_; }
char ResidueInformation::chainID()  const { return chainID_; }
int  ResidueInformation::resSeq()   const { return resSeq_; }
char ResidueInformation::iCode()    const { return iCode_; }
int  ResidueInformation::terCount() const { return terCount_; }
std::string const & ResidueInformation::segmentID() const { return segmentID_; }

std::string const & ResidueInformation::rosetta_resName() const { return rosetta_resName_; }

void ResidueInformation::resName(  std::string const & setting )
{
	resName_ = setting;
	rosetta_resName_ = core::io::NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( setting ).first;

	for ( Size ii = 1 ; ii <= atoms_.size(); ++ii ) atoms_[ ii ].resName = setting;
}
void ResidueInformation::chainID(  char setting ) { chainID_ = setting; }
void ResidueInformation::resSeq(   int setting ) { resSeq_ = setting; }
void ResidueInformation::iCode(    char setting ) { iCode_ = setting; }
void ResidueInformation::terCount( int setting ) { terCount_ = setting; }
void ResidueInformation::set_xyz( std::string const &atomname, Vector const &vect ) { xyz_[atomname] = vect; }
void ResidueInformation::set_temp( std::string const &atomname, core::Real const &val ) { temps_[atomname] = static_cast<double>(val); }
void ResidueInformation::segmentID( std::string const & setting ) { segmentID_ = setting; }

utility::vector1< AtomInformation > const & ResidueInformation::atoms() const { return atoms_; }
void ResidueInformation::append_atom( AtomInformation const & new_atom ) {
	atoms_.push_back( new_atom );
	if ( ! xyz_.count( new_atom.name ) ) {
		xyz_[ new_atom.name ] = Vector( new_atom.x, new_atom.y, new_atom.z );
		temps_[ new_atom.name ] = new_atom.temperature;
	}
}

void ResidueInformation::rename_atom( std::string const & orig_name, std::string const & new_name )
{
	if ( xyz_.find( orig_name ) == xyz_.end() ) {
		utility_exit_with_message( "Failed to find atom \"" + orig_name + "\" in ResidueInformation when trying to rename it to \"" + new_name );
	}
	for ( Size ii = 1; ii <= atoms_.size(); ++ii ) {
		if ( atoms_[ ii ].name == orig_name ) {
			atoms_[ ii ].name = new_name;
			// Do not insert this atom into the xyz_ and temps_ maps if the
			// the atom would displace an older atom
			if ( ! xyz_.count( new_name ) ) {
				xyz_[ new_name ] = xyz_[ orig_name ];
				temps_[ new_name ] = atoms_[ ii ].temperature;
			}
			xyz_.erase( orig_name );
			temps_.erase( orig_name );
			return;
		}
	}
}

/// @details Append all the atoms from the source residue, updating their
/// resName so that they all have the same resName, and so that the xyz_
/// and temps_ maps are updated.
void ResidueInformation::append_atoms( ResidueInformation const & source )
{
	for ( core::Size ii = 1; ii <= source.atoms_.size(); ++ii ) {
		AtomInformation const & iiat = source.atoms_[ ii ];
		if ( ! xyz_.count( iiat.name ) ) {
			AtomInformation new_at( iiat );
			new_at.resName = resName_;
			atoms_.push_back( new_at );
			xyz_[ iiat.name ] = Vector( iiat.x, iiat.y, iiat.z );
			temps_[ iiat.name ] = iiat.temperature;
		}
	}
}

std::map< std::string, Vector > const & ResidueInformation::xyz() const {
	return xyz_;
}

//< map of names to B-factors;  redundant but used a lot in reader
std::map< std::string, double > const & ResidueInformation::temps() const {
	return temps_;
}

std::string ResidueInformation::resid() const {
	std::string buf;
	buf.resize(7);
	// This is horribly hacky. Is this necessary?
	sprintf(&buf[0], "%4d%c%c", resSeq_, iCode_, chainID_ );
	buf.resize(6);
	return buf;
}

} // namespace io
} // namespace core
