// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pose/ResidueIndexDescription.cc
/// @brief  Classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/pose/ResidueIndexDescription.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/selection.hh>
#include <core/pose/chains_util.hh>

// Utility headers
#include <utility>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>
#include <ostream>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

std::string RID_Source::source_string() const {
	return "";
}

std::string
RID_FileSource::source_string() const {
	return  "line " + utility::to_string( linenum_ ) + " of the file named " + fname_;
}

ResidueIndexDescription::~ResidueIndexDescription() = default;

std::string
ResidueIndexDescription::source_string() const {
	if ( source_ == nullptr ) {
		return "";
	} else {
		return "from " + source_->source_string() + " ";
	}
}

core::Size
ResidueIndexDescription::do_error( bool no_error, std::string const & msg ) const {
	if ( !no_error ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, msg);
	}
	return 0;
}

std::ostream & operator<<( std::ostream & out, ResidueIndexDescription const & rid ) {
	rid.show(out);
	return out;
}

ResidueIndexDescriptionPoseNum::~ResidueIndexDescriptionPoseNum() = default;

core::Size
ResidueIndexDescriptionPoseNum::resolve_index(
	core::pose::Pose const & pose,
	bool no_error /*= false */
) const
{
	if ( pose_index_ == 0 ) {
		return do_error( no_error, "Residue index description " + source_string() + " is asking for non-existent residue zero.");
	}
	// If no_error is set, too-big numbers should be returned, rather than the zero from do_error().
	if ( pose_index_ > pose.size() && !no_error ) {
		return do_error( no_error, "Residue index description " + source_string() + "exceeds the number of residues in the Pose: Asked for residue " +
			utility::to_string( pose_index_ ) + " pose has only " + utility::to_string( pose.size() ) + " residues." );
	}
	return pose_index_;
}

void
ResidueIndexDescriptionPoseNum::show( std::ostream & out ) const {
	out << "Residue " << pose_index_;
}

ResidueIndexDescriptionCOP
make_rid_posenum( core::Size resnum ) {
	return ResidueIndexDescriptionCOP( new ResidueIndexDescriptionPoseNum( resnum ) );
}

ResidueIndexDescriptionPDB::~ResidueIndexDescriptionPDB() = default;

core::Size
ResidueIndexDescriptionPDB::resolve_index(
	core::pose::Pose const & pose,
	bool no_error /*= false */
) const
{
	pose::PDBInfoCOP info = pose.pdb_info();
	if ( info == nullptr ) {
		return do_error( no_error, "Cannot parse PDB residue " + source_string() + utility::to_string( resindex_ ) + insertion_code_ + ' ' + "on chain " + chain_ + " because the pose doesn't have PDB information");
	}
	core::Size resid = pose.pdb_info()->pdb2pose().find( chain_, resindex_, insertion_code_ );
	if ( resid == 0 ) {
		return do_error( no_error, "PDB residue " + source_string() + utility::to_string( resindex_ ) + insertion_code_ + ' ' + "on chain " + chain_ + " is not present in the input Pose." );
	}
	return resid;
}

void
ResidueIndexDescriptionPDB::show( std::ostream & out ) const {
	out << "Residue " << resindex_ << insertion_code_ << " chain " << chain_;
}

ResidueIndexDescriptionRefPose::~ResidueIndexDescriptionRefPose() = default;

core::Size
ResidueIndexDescriptionRefPose::resolve_index(
	core::pose::Pose const & pose,
	bool /*no_error = false */
) const
{
	return get_resnumber_from_reference_pose( refpose_name_, refpose_number_, refpose_offset_, pose );
}

void
ResidueIndexDescriptionRefPose::show( std::ostream & out ) const {
	out << "Residue offset " << refpose_offset_ << " from refpose " << refpose_name_ << " " << refpose_number_;
}

ResidueIndexDescriptionLastResidue::~ResidueIndexDescriptionLastResidue() = default;

core::Size
ResidueIndexDescriptionLastResidue::resolve_index(
	core::pose::Pose const & pose,
	bool /*no_error = false */
) const
{
	// TODO: consider if we should ignore virtual residues (e.g. virtual root) or accomodate symmetry.
	return pose.size();
}

void
ResidueIndexDescriptionLastResidue::show( std::ostream & out ) const {
	out << "Last residue in pose";
}

ResidueIndexDescriptionChainEnd::~ResidueIndexDescriptionChainEnd() = default;

core::Size
ResidueIndexDescriptionChainEnd::resolve_index(
	core::pose::Pose const & pose,
	bool no_error /*=false */
) const
{
	if ( chain_no_ != core::Size(-1) ) { // Chain of -1 means to use chain letter instead.
		if ( chain_no_ == 0 ) {
			return do_error( no_error, "Attempted to access non-existent chain 0 in residue index description " + source_string() );
		}
		if ( chain_no_ > pose.num_chains() ) {
			return do_error( no_error, "Attempted to access chain " + std::to_string( chain_no_ ) + " in residue index description " + source_string() + ": Pose only has " + std::to_string( pose.num_chains() ) + " chains." );
		}
		if ( chain_start_ ) {
			return pose.conformation().chain_begin( chain_no_ );
		} else {
			return pose.conformation().chain_end( chain_no_ );
		}
	} else if ( chain_letter_ ) { // Chain letter of 0x00 (not '0') means invalid chain letter.
		utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain_letter_, pose); // sorted order
		if ( chain_ids.size() == 0 ) {
			return do_error( no_error, "Attempted to access non-existent chain `" + std::to_string(chain_letter_) + "` in residue index description " + source_string() );
		}
		if ( chain_start_ ) {
			return pose.conformation().chain_begin( chain_ids[1] ); // Get the first residue which is part of this chain.
		} else {
			return pose.conformation().chain_end( chain_ids[ chain_ids.size() ] ); // Get last residue which is part of chain
		}
	} else {
		return do_error( no_error, "Attempted to resolve an invalid ResidueIndexDescriptionChainEnd " + source_string() );
	}

}

void
ResidueIndexDescriptionChainEnd::show( std::ostream & out ) const {
	out << ( chain_start_ ? "Starting":"Ending" ) << " residue of chain ";
	if ( chain_no_ != core::Size(-1) ) {
		out << chain_no_;
	} else if ( chain_letter_ ) {
		out << chain_letter_;
	} else {
		out << "???";
	}
}

// stupid visual studio has a different definition of std::isalpha than the C++ standard
// so I have to write my own version of isalpha...
bool
character_is_USA_letter( char c )
{
	return ( 'A' <= c && c <= 'Z' ) || ( 'a' <= c && c <= 'z' );
}

void
parse_PDBnum_icode(
	std::string const & token,
	std::string const & fname,
	Size const lineno,
	int & PDBnum,
	char & icode
)
{
	std::istringstream PDBnum_s;
	if ( character_is_USA_letter( *token.rbegin() ) ) {
		PDBnum_s.str(token.substr(0, token.length() - 1));
		icode = *token.rbegin();
	} else {
		PDBnum_s.str(token);
		icode = ' ';
	}

	char remaining;
	if ( !(PDBnum_s >> PDBnum) || PDBnum_s.get(remaining) ) {
		std::stringstream err_msg;
		err_msg
			<< "Problem in parse_PDBNum_icode. On line " << lineno
			<< " of the file '" << fname << "', "
			<< "the token '" << token << "' "
			<< "is not a valid <PDBNUM>[<ICODE>] identifier.";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg.str() );
	}

}

}
}

#ifdef    SERIALIZATION

template< class Archive >
void
core::pose::RID_Source::save( Archive & ) const {}

template< class Archive >
void
core::pose::RID_Source::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::RID_Source );
CEREAL_REGISTER_TYPE( core::pose::RID_Source )


template< class Archive >
void
core::pose::RID_FileSource::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::RID_Source >( this ) );
	arc( CEREAL_NVP( fname_ ) ); // std::string
	arc( CEREAL_NVP( linenum_ ) ); // core::Size
}

template< class Archive >
void
core::pose::RID_FileSource::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::RID_Source >( this ) );
	arc( fname_ ); // std::string
	arc( linenum_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::RID_FileSource );
CEREAL_REGISTER_TYPE( core::pose::RID_FileSource )


template< class Archive >
void
core::pose::ResidueIndexDescription::save( Archive & arc ) const {
	arc( CEREAL_NVP( source_ ) ); // ResidueIndexDescriptionCOP
}

template< class Archive >
void
core::pose::ResidueIndexDescription::load( Archive & arc ) {
	core::pose::RID_SourceOP source;
	arc( source );
	source_ = source;
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescription );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescription )


template< class Archive >
void
core::pose::ResidueIndexDescriptionPoseNum::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( CEREAL_NVP( pose_index_ ) ); // core::Size
}

template< class Archive >
void
core::pose::ResidueIndexDescriptionPoseNum::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( pose_index_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescriptionPoseNum );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescriptionPoseNum )


template< class Archive >
void
core::pose::ResidueIndexDescriptionPDB::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( CEREAL_NVP( chain_ ) ); // char
	arc( CEREAL_NVP( resindex_ ) ); // int
	arc( CEREAL_NVP( insertion_code_ ) ); // char
}

template< class Archive >
void
core::pose::ResidueIndexDescriptionPDB::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( chain_ ); // char
	arc( resindex_ ); // int
	arc( insertion_code_ ); // char
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescriptionPDB );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescriptionPDB )


template< class Archive >
void
core::pose::ResidueIndexDescriptionRefPose::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( CEREAL_NVP( refpose_name_ ) ); // std::string
	arc( CEREAL_NVP( refpose_number_ ) ); // core::Size
	arc( CEREAL_NVP( refpose_offset_ ) ); // signed long
}

template< class Archive >
void
core::pose::ResidueIndexDescriptionRefPose::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( refpose_name_ ); // std::string
	arc( refpose_number_ ); // core::Size
	arc( refpose_offset_ ); // signed long
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescriptionRefPose );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescriptionRefPose )


template< class Archive >
void
core::pose::ResidueIndexDescriptionLastResidue::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
}

template< class Archive >
void
core::pose::ResidueIndexDescriptionLastResidue::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescriptionLastResidue );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescriptionLastResidue )


template< class Archive >
void
core::pose::ResidueIndexDescriptionChainEnd::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( CEREAL_NVP( chain_start_ ) ); // bool
	arc( CEREAL_NVP( chain_no_ ) ); // core::Size
	arc( CEREAL_NVP( chain_letter_ ) ); // char
}

template< class Archive >
void
core::pose::ResidueIndexDescriptionChainEnd::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::ResidueIndexDescription >( this ) );
	arc( chain_start_ ); // bool
	arc( chain_no_ ); // core::Size
	arc( chain_letter_ ); // char
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::ResidueIndexDescriptionChainEnd );
CEREAL_REGISTER_TYPE( core::pose::ResidueIndexDescriptionChainEnd )


CEREAL_REGISTER_DYNAMIC_INIT( core_pose_ResidueIndexDescription )
#endif // SERIALIZATION
