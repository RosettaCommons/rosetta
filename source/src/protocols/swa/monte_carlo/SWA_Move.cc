// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/SWA_Move.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <utility/tools/make_vector1.hh>
#include <map>
#include <iostream>
#include <utility/string_util.hh>

using utility::make_tag_with_dashes;

namespace protocols {
namespace swa {
namespace monte_carlo {


	//Constructor
	SWA_Move::SWA_Move( MoveElement const & move_element,
											utility::vector1< Attachment > const & attachments,
											MoveType const & move_type ):
		utility::pointer::ReferenceCount(),
		move_element_( move_element ),
		attachments_( attachments ),
		move_type_( move_type )
	{
	}

	//Constructor
	SWA_Move::SWA_Move( MoveElement const & move_element,
											Attachment const & attachment,
											MoveType const & move_type ):
		utility::pointer::ReferenceCount(),
		move_element_( move_element ),
		attachments_( utility::tools::make_vector1( attachment ) ),
		move_type_( move_type )
	{
	}

	//Constructor
	SWA_Move::SWA_Move( Size const res,
											utility::vector1< Attachment > const & attachments,
											MoveType const & move_type ):
		utility::pointer::ReferenceCount(),
		move_element_( utility::tools::make_vector1( res ) ),
		attachments_( attachments ),
		move_type_( move_type )
	{
	}

	//Constructor
	SWA_Move::SWA_Move( Size const res,
											Attachment const & attachment,
											MoveType const & move_type ):
		utility::pointer::ReferenceCount(),
		move_element_( utility::tools::make_vector1( res ) ),
		attachments_( utility::tools::make_vector1( attachment ) ),
		move_type_( move_type )
	{
	}

	SWA_Move::SWA_Move():
		utility::pointer::ReferenceCount(),
		move_type_( NO_ADD_OR_DELETE )
	{
	}


	//Destructor
	SWA_Move::~SWA_Move()
	{}

	SWA_Move::SWA_Move( SWA_Move const & src ):
		ReferenceCount( src )
	{
		*this = src;
	}

	SWA_Move &
	SWA_Move::operator=( SWA_Move const & src )
	{
		move_element_ = src.move_element_;
		attachments_ = src.attachments_;
		move_type_ = src.move_type_;
		return *this;
	}

	Size
	SWA_Move::moving_res() const{
			runtime_assert( move_element_.size() == 1 );
			return move_element_[ 1 ];
		}

	Size
	SWA_Move::attached_res() const{
		runtime_assert( attachments_.size() == 1 );
		return attachments_[ 1 ].attached_res();
	}


	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////

	//Constructor
	Attachment::Attachment():
		utility::pointer::ReferenceCount(),
		attached_res_( 0 ),
		attachment_type_( NO_ATTACHMENT )
	{
	}

	//Constructor
	Attachment::Attachment( Size const & attached_res,
													 AttachmentType const attachment_type ):
		utility::pointer::ReferenceCount(),
		attached_res_( attached_res ),
		attachment_type_( attachment_type )
	{
	}

	//Destructor
	Attachment::~Attachment()
	{}

	Attachment::Attachment( Attachment const & src ):
		ReferenceCount( src )
	{
		*this = src;
	}

	Attachment &
	Attachment::operator=( Attachment const & src )
	{
		attached_res_ = src.attached_res_;
		attachment_type_ = src.attachment_type_;
		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	// following will be deprecated after we implement attachments
	// std::string
	// to_string( MovingResidueCase const & moving_residue_case ){

	// 	static bool init( false );
	// 	static std::map< MovingResidueCase, std::string> moving_residue_case_name;

	// 	if ( !init ){
	// 		moving_residue_case_name[ NO_CASE ] = "NO_CASE";
	// 		moving_residue_case_name[ CHAIN_TERMINUS_5PRIME ] = "CHAIN_TERMINUS_5PRIME";
	// 		moving_residue_case_name[ CHAIN_TERMINUS_3PRIME ] = "CHAIN_TERMINUS_3PRIME";
	// 		moving_residue_case_name[ INTERNAL ] = "INTERNAL";
	// 		moving_residue_case_name[ FLOATING_BASE ] = "FLOATING_BASE";
	// 		init = true;
	// 	}

	// 	return moving_residue_case_name[ moving_residue_case ];
	// }

	///////////////////////////////////////////////////////////////////////////////////////////
	std::string
	to_string( MoveType const & move_type ){

		static bool init( false );
		static std::map< MoveType, std::string> move_type_name;

		if ( !init ){
			move_type_name[ NO_ADD_OR_DELETE ] = "NO_ADD_OR_DELETE";
			move_type_name[ ADD ]    = "ADD";
			move_type_name[ DELETE ] = "DELETE";
			move_type_name[ RESAMPLE ] = "RESAMPLE";
			move_type_name[ RESAMPLE_INTERNAL_LOCAL ] = "RESAMPLE_INTERNAL_LOCAL";
			init = true;
		}

		return move_type_name[ move_type ];
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	std::string
	to_string( AttachmentType const & attachment_type ){

		static bool init( false );
		static std::map< AttachmentType, std::string> attachment_type_name;

		if ( !init ){
			attachment_type_name[ NO_ATTACHMENT ] = "NO_ATTACHMENT";
			attachment_type_name[ ATTACHED_TO_PREVIOUS ] = "ATTACHED_TO_PREVIOUS";
			attachment_type_name[ ATTACHED_TO_NEXT ] = "ATTACHED_TO_NEXT";
			attachment_type_name[ SKIP_BULGE_TO_PREVIOUS ] = "SKIP_BULGE_TO_PREVIOUS";
			attachment_type_name[ SKIP_BULGE_TO_NEXT ] = "SKIP_BULGE_TO_NEXT";
		}

		return attachment_type_name[ attachment_type ];
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	std::ostream &
	operator <<( std::ostream & os, SWA_Move const & swa_move )
	{
		os << "Choice " << to_string( swa_move.move_type() )	<<
			" res " << make_tag_with_dashes( swa_move.move_element() );
		if ( swa_move.attachments().size() == 0 ) {
			os << " with no attachments";
		} else {
			os << " with attachments ";
			for ( Size n = 1; n <= swa_move.attachments().size(); n++ ) os << " " << swa_move.attachments()[n];
		}
		return os;
	}


	/////////////////////////////////////////////////////////////////////////////////////////
	std::ostream &
	operator <<( std::ostream & os, Attachment const & attachment )
	{
		os << to_string( attachment.attachment_type() ) << ":" << attachment.attached_res();
		return os;
	}

} //monte_carlo
} //swa
} //protocols
