// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/CstInfoMover.hh
/// @brief A Mover to output information about constraints
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_CstInfoMover_hh
#define INCLUDED_protocols_simple_moves_CstInfoMover_hh

// Unit headers
#include <protocols/simple_moves/CstInfoMover.fwd.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace simple_moves {

///@brief A Mover to output information about constraints
class CstInfoMover : public protocols::moves::Mover {

public:

	CstInfoMover();

	// copy constructor
	CstInfoMover( CstInfoMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~CstInfoMover();

	void cst_file( std::string const & setting ) { cst_file_ = setting; }
	void prefix( std::string const & setting ) { prefix_ = setting; }
	void recursive( bool setting ) { recursive_ = setting; }

	std::string const & cst_file() const { return cst_file_; }
	std::string const & prefix() const { return prefix_; }
	bool recursive() const { return recursive_; }

	virtual void
	apply( core::pose::Pose & pose );

	void
	add_info_for_csts( core::scoring::constraints::ConstraintCOPs const & cstlist, core::pose::Pose & pose, std::string const & tag ) const;

	core::scoring::constraints::ConstraintCOPs
	get_constraints_from_file( std::string const & filename, core::pose::Pose const & pose ) const;

	core::scoring::constraints::ConstraintCOPs
	get_constraints_from_pose( core::pose::Pose const & pose ) const;

public:
	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//CstInfoMover & operator=( CstInfoMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

private:

	std::string cst_file_;
	std::string prefix_;
	bool recursive_;

};

std::ostream &operator<< (std::ostream &os, CstInfoMover const &mover);


} //protocols
} //simple_moves


#endif //protocols/simple_moves_CstInfoMover_hh







