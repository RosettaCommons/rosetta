// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMove.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_StepWiseMove_HH
#define INCLUDED_protocols_stepwise_monte_carlo_StepWiseMove_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.fwd.hh>
#include <core/types.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <ostream>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

typedef utility::vector1< core::Size >  MoveElement;
typedef utility::vector1< Attachment> Attachments;

// If you add something here, update to_string( AttachmenType ) in StepWiseMove.cc
enum AttachmentType { NO_ATTACHMENT = 0,
	BOND_TO_PREVIOUS,
	BOND_TO_NEXT,
	JUMP_TO_PREV_IN_CHAIN,
	JUMP_TO_NEXT_IN_CHAIN,
	JUMP_DOCK,
	SUBMOTIF, // special, used to fill submotif_tag_
	LAST_ATTACHMENT_TYPE };

// If you add something here, update to_string( MoveType ) in StepWiseMove.cc
enum MoveType { NO_MOVE = 0,
	ADD, // create new suite or jump
	DELETE, // split existing suite or jump
	FROM_SCRATCH, // create new residues (several)
	RESAMPLE, // resample existing suite or jump between two partitions
	RESAMPLE_INTERNAL_LOCAL, // ERRASER-style, KIC move.
	ADD_SUBMOTIF, // may be deprecated in factor of FROM_SCRATCH
	ADD_LOOP_RES, // add residue to design loop
	DELETE_LOOP_RES, // delete residue from design loop
	LAST_ADD_OR_DELETE_CHOICE };

std::string to_string( AttachmentType const & attachment_type );

std::string to_string( MoveType const & move_type_name );

MoveType
move_type_from_string( std::string const & name );

AttachmentType
attachment_type_from_string( std::string const & name );

/////////////////////////////////////////////////////////////////////
class Attachment: public utility::pointer::ReferenceCount {

public:

	//constructor
	Attachment();

	Attachment( Size const & attachment_res,
		AttachmentType const attachment_type );

	Attachment( Attachment const & src );

	friend
	bool
	operator==( Attachment const & a, Attachment const & b );

	//destructor
	~Attachment();

public:

	void set_attached_res( Size const & setting ){ attached_res_ = setting; }
	Size attached_res() const { return attached_res_; }

	void set_attachment_type( AttachmentType const & setting ){ attachment_type_ = setting; }
	AttachmentType attachment_type() const{ return attachment_type_; }

private:

	Size attached_res_;
	AttachmentType attachment_type_;

};

/////////////////////////////////////////////////////////////////////
class StepWiseMove: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWiseMove();

	StepWiseMove( MoveElement const & move_element,
		Attachments const & attachments,
		MoveType const & move_type );

	StepWiseMove( MoveElement const & move_element,
		Attachments const & attachments,
		MoveType const & move_type,
		std::string const & submotif_tag );

	StepWiseMove( MoveElement const & move_element,
		Attachment const & attachment,
		MoveType const & move_type );

	StepWiseMove( Size const moving_res,
		Attachments const & attachments,
		MoveType const & move_type );

	StepWiseMove( Size const moving_res,
		Attachment const & attachment,
		MoveType const & move_type );

	StepWiseMove( StepWiseMove const & src );

	StepWiseMove( utility::vector1< std::string > swa_move_string_vector,
		core::pose::full_model_info::FullModelParametersCOP full_model_parameters = 0 /* to convert resnum, chain to Rosetta res */ );

	friend
	bool
	operator==( StepWiseMove const & a, StepWiseMove const & b );

	friend
	bool
	operator!=( StepWiseMove const & a, StepWiseMove const & b ) {
		return !( a == b );
	}

	//destructor
	~StepWiseMove();

public:

	void set_move_element( MoveElement const & setting ){ move_element_ = setting; }
	MoveElement move_element() const{ return move_element_; }

	Size moving_res() const;

	Size attached_res() const;

	AttachmentType attachment_type() const;

	void set_attachments( Attachments const & setting ){ attachments_ = setting; }
	Attachments attachments() const{ return attachments_; }

	void set_move_type( MoveType const & setting ){ move_type_ = setting; }
	MoveType move_type() const{ return move_type_; }

	bool is_jump() { return ( attachment_type() == JUMP_TO_NEXT_IN_CHAIN ||
		attachment_type() == JUMP_TO_PREV_IN_CHAIN ||
		attachment_type() == JUMP_DOCK ); }

	std::string const & submotif_tag() const { return submotif_tag_; }

private:

	MoveElement move_element_;
	Attachments attachments_;
	MoveType move_type_;

	//kind of ad hoc. only filled in submotif creation.
	std::string submotif_tag_;

};

std::ostream &
operator <<( std::ostream & os, StepWiseMove const & swa_move );

std::ostream &
operator <<( std::ostream & os, Attachment const & attachment );

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
