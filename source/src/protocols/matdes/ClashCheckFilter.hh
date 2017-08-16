// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/ClashCheckFilter.hh
/// @brief  header file for ClashCheckFilter class
/// @author Neil King (neilking@u.washington.edu)


#ifndef INCLUDED_protocols_matdes_ClashCheckFilter_hh
#define INCLUDED_protocols_matdes_ClashCheckFilter_hh

// Unit Headers
#include <protocols/matdes/ClashCheckFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers

namespace protocols {
namespace matdes {

class ClashCheckFilter : public protocols::filters::Filter {
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
	ClashCheckFilter();

	// @brief constructor with arguments
	ClashCheckFilter( core::pack::task::TaskFactoryOP task_factory, core::Real const c, std::string s, core::Size const n, core::Size const t, bool const v, bool const w );

	// @brief copy constructor
	ClashCheckFilter( ClashCheckFilter const & rval );

	~ClashCheckFilter() override;


public:// virtual constructor


	// @brief make clone
	protocols::filters::FilterOP clone() const override;

	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override;


public:// accessor

	// @brief get name of this filter

public:// setters

	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void clash_dist( core::Real const c );
	void sym_dof_names( std::string const & s );
	void nsub_bblock( core::Size const n );
	void threshold( core::Size const t );
	void verbose( bool const v );
	void write( bool const w );

public:// getters
	core::pack::task::TaskFactoryOP task_factory() const;
	core::Real clash_dist() const;
	std::string sym_dof_names() const;
	core::Size nsub_bblock() const;
	core::Size threshold() const;
	bool verbose() const;
	bool write() const;

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
	core::Size compute( core::pose::Pose const & pose, bool const & v, bool const & w ) const;

	void write_to_pdb( core::pose::Pose const & pose, std::string const & residue_name, core::Size const residue, std::string const & atom_name ) const;
	void write_pymol_string_to_pdb( std::string const & pymol_selection ) const;

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
	core::Real clash_dist_;
	std::string sym_dof_names_;
	core::Size nsub_bblock_;
	core::Size threshold_;
	bool verbose_;
	bool write_;

};

} // matdes
} // protocols

#endif
