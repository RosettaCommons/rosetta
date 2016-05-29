// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMove.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <utility/tools/make_vector1.hh>
#include <map>
#include <iostream>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.StepWiseMove" );

using utility::make_tag_with_dashes;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

std::map< MoveType, std::string> move_type_name;
std::map< AttachmentType, std::string> attachment_type_name;

//Constructor
StepWiseMove::StepWiseMove( MoveElement const & move_element,
	utility::vector1< Attachment > const & attachments,
	MoveType const & move_type ):
	utility::pointer::ReferenceCount(),
	move_element_( move_element ),
	attachments_( attachments ),
	move_type_( move_type )
{
}

//Constructor
StepWiseMove::StepWiseMove( MoveElement const & move_element,
	utility::vector1< Attachment > const & attachments,
	MoveType const & move_type,
	std::string const submotif_tag ):
	utility::pointer::ReferenceCount(),
	move_element_( move_element ),
	attachments_( attachments ),
	move_type_( move_type ),
	submotif_tag_( submotif_tag )
{
}

//Constructor
StepWiseMove::StepWiseMove( MoveElement const & move_element,
	Attachment const & attachment,
	MoveType const & move_type ):
	utility::pointer::ReferenceCount(),
	move_element_( move_element ),
	attachments_( utility::tools::make_vector1( attachment ) ),
	move_type_( move_type )
{
}

//Constructor
StepWiseMove::StepWiseMove( Size const res,
	utility::vector1< Attachment > const & attachments,
	MoveType const & move_type ):
	utility::pointer::ReferenceCount(),
	move_element_( utility::tools::make_vector1( res ) ),
	attachments_( attachments ),
	move_type_( move_type )
{
}

//Constructor
StepWiseMove::StepWiseMove( Size const res,
	Attachment const & attachment,
	MoveType const & move_type ):
	utility::pointer::ReferenceCount(),
	move_element_( utility::tools::make_vector1( res ) ),
	attachments_( utility::tools::make_vector1( attachment ) ),
	move_type_( move_type )
{
}

StepWiseMove::StepWiseMove():
	utility::pointer::ReferenceCount(),
	move_type_( NO_MOVE )
{
}

StepWiseMove::StepWiseMove( utility::vector1< std::string > swa_move_string_vector,
	core::pose::full_model_info::FullModelParametersCOP full_model_parameters /* = 0 */ ):
	utility::pointer::ReferenceCount(),
	move_type_( NO_MOVE )
{
	using namespace ObjexxFCL;
	Size const num_strings = swa_move_string_vector.size();
	if ( num_strings == 0 ) return; // blank.
	//  runtime_assert( num_strings >= 4 );

	Size n( 1 );
	move_type_ = move_type_from_string( swa_move_string_vector[ n ] ); n++;

	bool string_is_ok( false );
	while ( n <= num_strings ) {
		//   std::vector <int> ints = ints_of( swa_move_string_vector[ n ], string_is_ok );
		std::vector< int > resnum;
		std::vector< char > chains;
		string_is_ok = utility::get_resnum_and_chain_from_one_tag( swa_move_string_vector[ n ], resnum, chains );
		if ( !string_is_ok ) break;
		for ( Size i = 0; i < resnum.size(); i++ ) {
			if ( full_model_parameters ) {
				move_element_.push_back( full_model_parameters->conventional_to_full( resnum[i], chains[i] ) );
			} else {
				move_element_.push_back( resnum[i] );
			}
		}
		n++;
	}

	//  runtime_assert( n != num_strings );
	while ( n <= num_strings ) {
		AttachmentType attachment_type = attachment_type_from_string( swa_move_string_vector[ n ] );
		n++;

		if ( attachment_type == SUBMOTIF ) {
			submotif_tag_ = swa_move_string_vector[ n ];
			n++;
		} else {
			//   Size attachment_res = int_of( swa_move_string_vector[ n ] );
			// attachment res should be one integer.
			std::vector< int > resnum;
			std::vector< char > chains;
			string_is_ok = utility::get_resnum_and_chain_from_one_tag( swa_move_string_vector[ n ], resnum, chains );
			runtime_assert( string_is_ok );
			runtime_assert( resnum.size() == 1 );
			Size attachment_res = resnum[0];
			if ( full_model_parameters ) attachment_res = full_model_parameters->conventional_to_full( resnum[0], chains[0] );
			n++;

			attachments_.push_back( Attachment( attachment_res, attachment_type ) );
		}
	}

	if ( move_type_ == RESAMPLE_INTERNAL_LOCAL ) {
		runtime_assert( attachments_.size() == 2 );
		runtime_assert( attachments_[ 1 ].attachment_type() == BOND_TO_PREVIOUS );
		runtime_assert( attachments_[ 2 ].attachment_type() == BOND_TO_NEXT );
	}
}


//Destructor
StepWiseMove::~StepWiseMove()
{}

StepWiseMove::StepWiseMove( StepWiseMove const & src ):
	ReferenceCount( src )
{
	*this = src;
}

Size
StepWiseMove::moving_res() const {
	runtime_assert( move_element_.size() == 1 );
	return move_element_[ 1 ];
}

Size
StepWiseMove::attached_res() const {
	runtime_assert( attachments_.size() == 1 );
	return attachments_[ 1 ].attached_res();
}

AttachmentType
StepWiseMove::attachment_type() const {
	runtime_assert( attachments_.size() == 1 );
	return attachments_[ 1 ].attachment_type();
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

///////////////////////////////////////////////////////////////////////////////////////////
void
initialize_move_type_name(){
	static bool init( false );
	if ( !init ) {
		move_type_name[ NO_MOVE ] = "NO_MOVE";
		move_type_name[ ADD ]    = "ADD";
		move_type_name[ DELETE ] = "DELETE";
		move_type_name[ FROM_SCRATCH ] = "FROM_SCRATCH";
		move_type_name[ RESAMPLE ] = "RESAMPLE";
		move_type_name[ RESAMPLE_INTERNAL_LOCAL ] = "RESAMPLE_INTERNAL_LOCAL";
		move_type_name[ ADD_SUBMOTIF ]    = "ADD_SUBMOTIF";
		move_type_name[ ADD_LOOP_RES ]    = "ADD_LOOP_RES";
		move_type_name[ DELETE_LOOP_RES ] = "DELETE_LOOP_RES";
		init = true;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
std::string
to_string( MoveType const & move_type ){
	initialize_move_type_name();
	return move_type_name[ move_type ];
}
///////////////////////////////////////////////////////////////////////////////////////////
MoveType
move_type_from_string( std::string const & name ){
	initialize_move_type_name();
	MoveType move_type( NO_MOVE );
	for ( std::map< MoveType, std::string>::const_iterator it = move_type_name.begin(),
			end = move_type_name.end(); it != end; ++it ) {
		if ( it->second == name ) {
			move_type = it->first;
		}
	}
	runtime_assert( move_type );
	return move_type;
}

///////////////////////////////////////////////////////////////////////////////////////////
void
initialize_attachment_type_name(){
	static bool init( false );
	if ( !init ) {
		attachment_type_name[ NO_ATTACHMENT ] = "NO_ATTACHMENT";
		attachment_type_name[ BOND_TO_PREVIOUS ] = "BOND_TO_PREVIOUS";
		attachment_type_name[ BOND_TO_NEXT ] = "BOND_TO_NEXT";
		attachment_type_name[ JUMP_TO_PREV_IN_CHAIN ] = "JUMP_TO_PREV_IN_CHAIN";
		attachment_type_name[ JUMP_TO_NEXT_IN_CHAIN ] = "JUMP_TO_NEXT_IN_CHAIN";
		attachment_type_name[ JUMP_DOCK ] = "JUMP_DOCK";
		attachment_type_name[ SUBMOTIF ] = "SUBMOTIF";
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
std::string
to_string( AttachmentType const & attachment_type ){
	initialize_attachment_type_name();
	return attachment_type_name[ attachment_type ];
}

///////////////////////////////////////////////////////////////////////////////////////////
AttachmentType
attachment_type_from_string( std::string const & name ){
	initialize_attachment_type_name();
	AttachmentType attachment_type( NO_ATTACHMENT );
	for ( std::map< AttachmentType, std::string>::const_iterator it = attachment_type_name.begin(),
			end = attachment_type_name.end(); it != end; ++it ) {
		if ( it->second == name ) {
			attachment_type = it->first;
		}
	}
	runtime_assert( attachment_type );
	return attachment_type;
}

/////////////////////////////////////////////////////////////////////////////////////////
std::ostream &
operator <<( std::ostream & os, StepWiseMove const & swa_move )
{
	os << "Choice " << to_string( swa_move.move_type() ) <<
		" res " << make_tag_with_dashes( swa_move.move_element() );
	if ( swa_move.attachments().size() == 0 ) {
		os << " with no attachments";
	} else {
		os << " with attachments ";
		for ( Size n = 1; n <= swa_move.attachments().size(); n++ ) os << " " << swa_move.attachments()[n];
	}
	if ( swa_move.submotif_tag().size() > 0 ) {
		os << " submotif_tag " << swa_move.submotif_tag();
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

/////////////////////////////////////////////////////////////////////////////////////////
bool
operator==( Attachment const & a, Attachment const & b ) {
	return ( ( a.attached_res_     == b.attached_res_ ) &&
		( a.attachment_type_  == b.attachment_type_ ) );
}

/////////////////////////////////////////////////////////////////////////////////////////
bool
operator==( StepWiseMove const & a, StepWiseMove const & b ) {

	return ( ( a.move_element_ == b.move_element_ ) &&
		( a.attachments_  == b.attachments_ ) &&
		( a.move_type_    == b.move_type_ ) &&
		( a.submotif_tag_ == b.submotif_tag_ ) );

}


} //mover
} //monte_carlo
} //stepwise
} //protocols
