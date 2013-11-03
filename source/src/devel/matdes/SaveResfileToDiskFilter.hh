// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/SaveResfileToDiskFilter.hh
/// @brief  header file for SaveResfileToDiskFilter class
/// @author Neil King (neilking@u.washington.edu)


#ifndef INCLUDED_devel_matdes_SaveResfileToDiskFilter_hh
#define INCLUDED_devel_matdes_SaveResfileToDiskFilter_hh

// Unit Headers
#include <devel/matdes/SaveResfileToDiskFilter.fwd.hh>

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

namespace devel {
namespace matdes {

class SaveResfileToDiskFilter : public protocols::filters::Filter {
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
	SaveResfileToDiskFilter();

	// @brief constructor with arguments
	SaveResfileToDiskFilter( core::pack::task::TaskFactoryOP task_factory,
			utility::vector1<core::Size> selected_resis,
			bool designable_only,
			std::string resfile_name,
			std::string resfile_suffix,
			std::string resfile_prefix,
			std::string resfile_general_property,
			std::string selected_resis_property
	); 

	// @brief copy constructor
	SaveResfileToDiskFilter( SaveResfileToDiskFilter const & rval );

	virtual ~SaveResfileToDiskFilter();


public:// virtual constructor

	// @brief make clone
	virtual protocols::filters::FilterOP clone() const;

	// @brief make fresh instance
	virtual protocols::filters::FilterOP fresh_instance() const;


public:// accessor

	// @brief get name of this filter
	virtual std::string name() const { return "SaveResfileToDisk"; }

public:// setters

	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void selected_resis( utility::vector1<core::Size> const r );
	void designable_only( bool const d );
	void resfile_name( std::string const n );
	void resfile_suffix( std::string const s );
	void resfile_prefix( std::string const p );
	void resfile_general_property( std::string const g );
	void selected_resis_property( std::string const g );

public:// getters
	core::pack::task::TaskFactoryOP task_factory() const;
	utility::vector1< core::Size > selected_resis() const;
	bool designable_only() const;
	std::string resfile_name() const;
	std::string resfile_suffix() const;
	std::string resfile_prefix() const;
	std::string resfile_general_property() const;
	std::string selected_resis_property() const;

public:// parser

	virtual void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & );


public:// virtual main operation

	// @brief returns true if the given pose passes the filter, false otherwise.
  virtual bool apply( core::pose::Pose const & pose ) const;

	// @brief public functions
	utility::vector1< core::Size > select_residues( Pose const & pose ) const;
	void write_resfile( Pose const & pose, utility::vector1< core::Size > const & selected_residues ) const;

	/// @brief
//  virtual core::Real report_sm( core::pose::Pose const & pose ) const;
//  virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;


private:

  core::pack::task::TaskFactoryOP task_factory_;
	utility::vector1< core::Size > selected_resis_;
	bool designable_only_;
	std::string resfile_name_;
	std::string resfile_suffix_;
	std::string resfile_prefix_;
	std::string resfile_general_property_;
	std::string selected_resis_property_;

};

} // matdes
} // devel

#endif
