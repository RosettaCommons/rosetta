// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_movers/CstInfoMover.hh
/// @brief A Mover to output information about constraints
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_constraint_movers_CstInfoMover_hh
#define INCLUDED_protocols_constraint_movers_CstInfoMover_hh

// Unit headers
#include <protocols/constraint_movers/CstInfoMover.fwd.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace constraint_movers {

///@brief A Mover to output information about constraints
class CstInfoMover : public protocols::moves::Mover {

public:

	CstInfoMover();

	// copy constructor
	CstInfoMover( CstInfoMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~CstInfoMover();

	void cst_file( std::string const & setting ) { cst_file_ = setting; }
	void dump_cst_file( std::string const & setting ) { dump_cst_file_ = setting; }
	void prefix( std::string const & setting ) { prefix_ = setting; }
	void recursive( bool setting ) { recursive_ = setting; }

	std::string const & cst_file() const { return cst_file_; }
	std::string const & dump_cst_file() const { return dump_cst_file_; }
	std::string const & prefix() const { return prefix_; }
	bool recursive() const { return recursive_; }

	void
	apply( core::pose::Pose & pose ) override;

	void
	add_info_for_csts( core::scoring::constraints::ConstraintCOPs const & cstlist, core::pose::Pose & pose, std::string const & tag ) const;

	core::scoring::constraints::ConstraintCOPs
	get_constraints_from_file( std::string const & filename, core::pose::Pose const & pose ) const;

	core::scoring::constraints::ConstraintCOPs
	get_constraints_from_pose( core::pose::Pose const & pose ) const;

public:
	void
	show( std::ostream & output=std::cout ) const override;


	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//CstInfoMover & operator=( CstInfoMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	std::string cst_file_;
	std::string dump_cst_file_;
	std::string prefix_;
	bool recursive_;

};

std::ostream &operator<< (std::ostream &os, CstInfoMover const &mover);


} //protocols
} //constraint_movers


#endif //protocols/constraint_movers_CstInfoMover_hh
