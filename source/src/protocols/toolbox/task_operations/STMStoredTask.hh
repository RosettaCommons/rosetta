// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/STMStoredTask.hh
/// @brief
/// @author Neil King (neilking@uw.edu)


#ifndef INCLUDED_protocols_toolbox_task_operations_STMStoredTask_hh
#define INCLUDED_protocols_toolbox_task_operations_STMStoredTask_hh

#include <protocols/toolbox/task_operations/STMStoredTask.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/types.hh>
#include <map>

namespace protocols {
namespace toolbox {
namespace task_operations {

class STMStoredTask: public basic::datacache::CacheableData {
public:

	// constructors
	STMStoredTask();
	STMStoredTask( const STMStoredTask& rval );
	virtual basic::datacache::CacheableDataOP clone() const;
	virtual basic::datacache::CacheableDataOP fresh_instance() const;

	// setter
	void set_task( core::pack::task::PackerTaskOP const task, std::string const task_name );

	// getters
	core::pack::task::PackerTaskOP get_task( std::string task_name ) const;
	bool has_task( std::string const task_name ) const;

private:
#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const {
			// manually do this because we don't support OP of abstract class
			int size = tasks_.size();
			ar & size;
			for( std::map< std::string, core::pack::task::PackerTaskOP >::const_iterator itr = tasks_.begin(); itr != tasks_.end(); itr++ ) {
				ar & itr->first;
				core::pack::task::PackerTask_ * tmp = static_cast<core::pack::task::PackerTask_ *> (itr->second.get());
				ar & tmp;
			}
	}
	template<class Archive>
	void load(Archive & ar, const unsigned int version) {
			int size;
			ar & size; //tasks_.size();
			for( int i = 0; i < size; i++ ) {
				std::string tmp_string;
				core::pack::task::PackerTask_ * tmpptr;
				ar & tmp_string;
				ar & tmpptr;
				core::pack::task::PackerTaskOP tmpop(tmpptr);
				tasks_[tmp_string] = tmpop;
			}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

	std::map< std::string, core::pack::task::PackerTaskOP > tasks_;

};

} // task_operations
} // toolbox
} // protocols

#endif
