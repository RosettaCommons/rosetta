// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/OligomericAverageDegreeFilter.hh
/// @brief  header file for OligomericAverageDegreeFilter class
/// @author Neil King (neilking@u.washington.edu)


#ifndef INCLUDED_protocols_matdes_OligomericAverageDegreeFilter_hh
#define INCLUDED_protocols_matdes_OligomericAverageDegreeFilter_hh

// Unit Headers
#include <protocols/matdes/OligomericAverageDegreeFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace matdes {

class OligomericAverageDegreeFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	OligomericAverageDegreeFilter();

	// @brief constructor with arguments
	OligomericAverageDegreeFilter( core::pack::task::TaskFactoryOP task_factory, core::Real const t, core::Real const d, bool jump_set, core::Size jump, std::string dof_names, bool mcomp );

	// @brief copy constructor
	OligomericAverageDegreeFilter( OligomericAverageDegreeFilter const & rval );

	~OligomericAverageDegreeFilter() override;


public:// virtual constructor


	// @brief make clone
	protocols::filters::FilterOP clone() const override;

	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override;


public:// accessor

	// @brief get name of this filter

public:// setters

	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void threshold( core::Real const t );
	void distance_threshold( core::Real const d );
	void jump_set( bool const jump_set );
	void jump_id( core::Size const jump );
	void sym_dof_names( std::string const & dof_names );
	void write2pdb( bool const write );
	void verbose( bool const verb );
	void multicomp( bool const multicomp );

public:// getters
	core::pack::task::TaskFactoryOP task_factory() const;
	core::Real threshold() const;
	core::Real distance_threshold() const;
	bool jump_set() const;
	core::Size jump_id() const;
	std::string sym_dof_names() const;
	bool write2pdb() const;
	bool verbose() const;
	bool multicomp() const;

public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & ) override;

public:// virtual main operation

	// @brief returns true if the given pose passes the filter, false otherwise.
	bool apply( core::pose::Pose const & pose ) const override;

	/// @brief
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	/// @brief calc oligomeric AverageDegree
	core::Real compute( core::pose::Pose const & pose, bool const & verbose, bool const & write ) const;
	void write_to_pdb( core::pose::Pose const & pose, core::Size const residue, std::string const & residue_name, core::Size const neighbors ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	core::pack::task::TaskFactoryOP task_factory_;
	core::Real threshold_;
	core::Real distance_threshold_;
	bool jump_set_;
	core::Size jump_id_;
	std::string sym_dof_names_;
	bool write2pdb_, verbose_, multicomp_;

};

} // matdes
} // protocols

#endif
