// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_moves_symmetry_SetupForSymmetryMover_hh

#define INCLUDED_protocols_moves_symmetry_SetupForSymmetryMover_hh

// Unit headers
#include <protocols/moves/symmetry/SetupForSymmetryMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers

namespace protocols {
namespace moves {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////

class SetupForSymmetryMover : public Mover
{
public:

	// default constructor
	SetupForSymmetryMover();

	SetupForSymmetryMover( std::string const & );

	~SetupForSymmetryMover();

	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new SetupForSymmetryMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose );
	virtual void parse_my_tag( 
			utility::tag::TagPtr const tag,
			moves::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

	virtual std::string get_name() const;

	// setter
	void slide_into_contact(bool val) { slide_ = val; }

private:
	bool slide_;
	std::string symmdef_file_;
};

///////////////

class ExtractAsymmetricUnitMover : public Mover
{
public:

	// default constructor
	ExtractAsymmetricUnitMover();

	~ExtractAsymmetricUnitMover();

	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new ExtractAsymmetricUnitMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose );
	virtual void parse_my_tag( 
			utility::tag::TagPtr const tag,
			moves::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

	virtual std::string get_name() const;
};

} // symmetry
} // moves
} // rosetta
#endif
