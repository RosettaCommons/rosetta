// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/conformation/signals/LengthEvent.hh
/// @brief  signal for a change in length of residues in a Conformation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_signals_LengthEvent_hh
#define INCLUDED_core_conformation_signals_LengthEvent_hh

// type headers
#include <core/types.hh>

// unit headers
#include <core/conformation/signals/LengthEvent.fwd.hh>

// package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/Conformation.fwd.hh>

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/exit.hh>

#include <iosfwd>

namespace core {
namespace conformation {
namespace signals {


/// @brief signals a change in length of residues in a Conformation
/// @remarks When accessing residue information, take care as to which
///  data member you choose.  For almost all situations the ResidueCAP
///  'residue' should be used instead of the Conformation.  See remarks
///  below.
struct LengthEvent : public GeneralEvent {


	// typedefs
	typedef core::Size Size;
	typedef GeneralEvent Super;


	/// @brief the type of length change
	enum Tag {
		EMPTY, // empty event, e.g. for default constructor
		INVALIDATE, // for safety, e.g. if during copy assignment something radical occurs
		RESIDUE_APPEND,
		RESIDUE_PREPEND,
		RESIDUE_DELETE
	};


	/// @brief default constructor
	inline
	LengthEvent() :
		Super(),
		tag( EMPTY ),
		position( 0 ),
		conformation_size( 0 ),
		length_change( 0 ),
		residue()
	{}


	/// @brief constructor
	/// @param t type of length change
	/// @param pos residue position
	LengthEvent(
		Conformation const * conf,
		Tag const t,
		Size const & pos,
		int const & len_chg,
		Residue const * res
	);

	/// @brief copy constructor
	inline
	LengthEvent( LengthEvent const & rval ) :
		Super( rval ),
		tag( rval.tag ),
		position( rval.position ),
		conformation_size( rval.conformation_size ),
		length_change( rval.length_change ),
		residue( rval.residue )
	{
#ifndef NDEBUG
		check_consistency();
#endif
	}


	/// @brief default destructor
	inline
	virtual
	~LengthEvent() {}


	/// @brief copy assignment
	inline
	LengthEvent &
	operator =( LengthEvent const & rval ) {
#ifndef NDEBUG
		rval.check_consistency();
#endif
		if ( this != &rval ) {
			Super::operator =( rval );
			tag = rval.tag;
			position = rval.position;
			conformation_size = rval.conformation_size;
			length_change = rval.length_change;
			residue = rval.residue;
		}
		return *this;
	}

	void
	check_consistency() const;

	/// @brief tag indicating type of length change
	Tag tag;

	/// @brief residue position where the event happened
	Size position;

	/// @brief The length of the conformation after the length event
	/// @details This is needed for things like regenerating mappings
	/// after multiple applications of events which may invalidate the confromation.
	Size conformation_size;

	/// @brief overall length change of the conformation
	int length_change;

	/// @brief direct access to residue
	/// @remarks Almost always want to use this to access the residue instead of
	///  the conformation.  Calling Conformation::residue() can cause an internal
	///  update/re-sync inside Pose, which may have consequences if you're depending
	///  upon multiple residue operations to be setup (such as bond angle/length
	///  changes) prior to an internal update.
	Residue const * residue;


};

std::ostream &
operator<<( std::ostream & out, LengthEvent const & event );

} // namespace signals
} // namespace conformation
} // namespace core


#endif /* INCLUDED_core_conformation_signals_LengthEvent_HH */
