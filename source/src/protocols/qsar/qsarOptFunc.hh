// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/qsarOptFunc.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_qsarOptFunc_HH
#define INCLUDED_protocols_qsar_qsarOptFunc_HH

#include <core/optimization/Multifunc.hh>
#include <protocols/qsar/qsarOptFunc.fwd.hh>
//external headers
#include <cppdb/frontend.h>

#include <map>
#include <list>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

namespace protocols {
namespace qsar {

struct qsarOptData
{
	bool activity;
	std::map<std::string,core::Real> score_map;
	std::string tag;
};


class qsarOptFunc : public core::optimization::Multifunc
{

public:

	qsarOptFunc(
		utility::sql_database::sessionOP db_session,
		core::optimization::Multivec const & initial_values,
		std::map<std::string,core::Size> const & grid_indices);

	virtual ~qsarOptFunc() {}


	void setup_data_map();

	void set_initial_values(core::optimization::Multivec const & initial_values);

	// func
	virtual
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	// dfunc
	virtual
	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const;

	/// @brief Error state reached -- derivative does not match gradient
	virtual
	void
	dump( core::optimization::Multivec const & vars, core::optimization::Multivec const & vars2 ) const;

private:

	qsarOptData get_struct_data(core::Size const & struct_id);

	std::list<qsarOptData> data_map_;
	core::optimization::Multivec initial_values_;
	std::map<std::string,core::Size> grid_indices_;
	cppdb::statement score_selection_;
	cppdb::statement struct_id_selection_;
	cppdb::statement tag_activity_selection_;
	core::Real cutoff_;

};


}
}


#endif /* QSAROPTFUNC_HH_ */
