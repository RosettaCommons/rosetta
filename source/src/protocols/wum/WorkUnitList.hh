// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/WorkUnitList.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_WorkUnitList_hh
#define INCLUDED_protocols_wum_WorkUnitList_hh


#include <core/types.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace wum {

class WorkUnitList;

/// @ brief just a prettier name for a pointer to a work unitlist.
typedef  const WorkUnitList* WorkUnitListCAP;

/// @brief WOrkUnitList is a store for WorkUnitCOPs. THe purpose of this class is to store all the possible WorkUnits
///  that a protocol might need. When needed these are cloned and then used. THis class acts like a library of WorkUnit blueprints.
class WorkUnitList{
public:
	WorkUnitList(){}

	/// @brief Add a WorkUnit to the list, each workunit must be named with a string that is used later to retrieve it !
	void register_work_unit( const std::string &name, WorkUnitBaseOP the_work_unit );

	/// @brief Return a COP to a workunit with a given name. If multiple WUs were registered with the same name, the first is returned.
	const WorkUnitBaseCOP get_work_unit( const std::string &name ) const;

	/// @brief Return an OP to a workunit with a given name but clone it.
	const WorkUnitBaseOP  get_work_unit_clone( const std::string &name ) const { return get_work_unit( name )->clone(); }

	/// @brief Return an OP to a workunit with the same name as the one given as a parameter.
	const WorkUnitBaseCOP get_work_unit( const WorkUnitBase &wu ) const { return get_work_unit( wu.get_wu_type() ); }

	/// @brief Return an OP to a workunit with the same name as the one given as a parameter, but as a clone
	const WorkUnitBaseOP get_work_unit_clone( const WorkUnitBase &wu ) const { return get_work_unit( wu.get_wu_type() )->clone(); }

	void merge( const WorkUnitList & source );
protected:

	/// An STL map is used to associate the WUs with strings.
	std::map< std::string, WorkUnitBaseCOP > work_unit_list_;
};


}
}

#endif

