// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/SWA_Move.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_monte_carlo_SWA_Move_HH
#define INCLUDED_protocols_swa_monte_carlo_SWA_Move_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/monte_carlo/SWA_Move.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>
#include <ostream>

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	typedef utility::vector1< core::Size >  MoveElement;
	typedef utility::vector1< Attachment> Attachments;

	enum AttachmentType{ NO_ATTACHMENT = 0, ATTACHED_TO_PREVIOUS, ATTACHED_TO_NEXT,
											 SKIP_BULGE_TO_PREVIOUS, SKIP_BULGE_TO_NEXT, LAST_ATTACHMENT_TYPE };

	enum MoveType{ NO_ADD_OR_DELETE = 0, ADD, DELETE, RESAMPLE, RESAMPLE_INTERNAL_LOCAL, LAST_ADD_OR_DELETE_CHOICE };

	std::string to_string( AttachmentType const & attachment_type );

	std::string to_string( MoveType const & move_type_name );

	/////////////////////////////////////////////////////////////////////
	class Attachment: public utility::pointer::ReferenceCount {

	public:

		//constructor
		Attachment();

		Attachment( Size const & attachment_res,
								AttachmentType const attachment_type );

		Attachment( Attachment const & src );

		Attachment &
		operator=( Attachment const & src );

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
	class SWA_Move: public utility::pointer::ReferenceCount {

	public:

		//constructor
		SWA_Move();

		SWA_Move( MoveElement const & move_element,
							Attachments const & attachments,
							MoveType const & move_type );

		SWA_Move( Size const moving_res,
							Attachments const & attachments,
							MoveType const & move_type );

		SWA_Move( Size const moving_res,
							Attachment const & attachments,
							MoveType const & move_type );

		SWA_Move( SWA_Move const & src );

		SWA_Move &
		operator=( SWA_Move const & src );

		//destructor
		~SWA_Move();

	public:

		void set_move_element( MoveElement const & setting ){ move_element_ = setting; }
		MoveElement move_element() const{ return move_element_; }

		Size moving_res() const;

		Size attached_res() const;

		void set_attachments( Attachments const & setting ){ attachments_ = setting; }
		Attachments attachments() const{ return attachments_; }

		void set_move_type( MoveType const & setting ){ move_type_ = setting; }
		MoveType move_type() const{ return move_type_; }

	private:

		MoveElement move_element_;
		Attachments attachments_;
		MoveType move_type_;

	};

	std::ostream &
	operator <<( std::ostream & os, SWA_Move const & swa_move );

	std::ostream &
	operator <<( std::ostream & os, Attachment const & attachment );

} //monte_carlo
} //swa
} //protocols

#endif
