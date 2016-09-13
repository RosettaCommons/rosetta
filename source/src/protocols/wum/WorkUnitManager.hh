// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/WorkUnitManager.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_WorkUnitManager_hh
#define INCLUDED_protocols_wum_WorkUnitManager_hh

#include <protocols/wum/WorkUnitManager.fwd.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/WorkUnitList.hh>

#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <string>
#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace wum {

class WorkUnitQueue;
class WorkUnitQueue_Swapped;
class WorkUnitManager;

enum WUM_MPI_TAGS {
	WUM_MPI_REQUEST_WU=101,     // A slave wants a new job
	WUM_MPI_SEND_WU,        // A slave wants to send the master a new job
	WUM_MPI_DATA_BLOCK,     // Tag for the data blocks
	WUM_MPI_SPINDOWN,           // A Master wants to shut down a slave ready for shut down
	WUM_MPI_end
};

class WorkUnitQueue {
public:
	typedef std::list < WorkUnitBaseOP >::iterator  iterator;
	typedef std::list < WorkUnitBaseOP >::const_iterator  const_iterator;

public:
	WorkUnitQueue():memory_limit_(0) {};

	virtual ~WorkUnitQueue() {};

	virtual core::Size  size() const { return wus_.size(); }

	virtual void add( WorkUnitBaseOP new_wu )       { if ( is_under_memory_limit() ) wus_.push_back( new_wu ); }
	virtual void push_back( WorkUnitBaseOP new_wu ) { if ( is_under_memory_limit() ) wus_.push_back( new_wu ); }
	virtual void push_front( WorkUnitBaseOP new_wu ){ if ( is_under_memory_limit() ) wus_.push_front( new_wu ); }

	virtual WorkUnitBaseOP &next();
	virtual WorkUnitBaseOP pop_next();
	virtual iterator erase( iterator i );

	iterator       begin()       { return wus_.begin(); }
	const_iterator begin() const { return wus_.begin(); }
	iterator       end()         { return wus_.end(); }
	const_iterator end() const   { return wus_.end(); }

	virtual void clear() { wus_.clear(); }

	/// @brief return total memory foot print in bytes
	core::Size mem_foot_print() const;

	/// @brief report number of total structures, and memory foot prints
	void mem_stats( core::Size &n_structs, core::Size &structs_memory, core::Size &WU_memory ) const;

	void set_memory_limit( core::Size memory_limit ) {
		memory_limit_ = memory_limit;
	}

	bool is_under_memory_limit() const {
		if ( memory_limit_ == 0 ) return true;
		if ( mem_foot_print() < memory_limit_ ) return true;
		return false;
	}

protected:
	std::list < WorkUnitBaseOP > wus_;

private:
	core::Size memory_limit_;
};


// as above but uses a disk-swap to prevent overflows
class WorkUnitQueue_Swapped: public WorkUnitQueue {
public:
	WorkUnitQueue_Swapped( WorkUnitManager *wum, std::string  swap_file, core::Size memory_limit ):
		WorkUnitQueue(),
		swap_file_(std::move( swap_file )),
		memory_limit_( memory_limit ),
		wum_( wum )
	{
	}

	~WorkUnitQueue_Swapped() {};

	//virtual core::Size  size();

	void add( WorkUnitBaseOP new_wu ) override;

	const std::string &swap_file() const { return swap_file_; }

protected:

	virtual void add_to_swap( WorkUnitBaseOP new_wu );

private:
	void set_swap_file( const std::string &swap_file ){ swap_file_ = swap_file; }
	std::string swap_file_;

	WorkUnitQueue swap_buffer_;

	// parameters
	core::Size  memory_limit_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size  read_swap_limit_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size  clean_swap_limit_;
	core::Size  max_swap_buffer_size_;
	// variables
	core::Size  n_swap_total_;
	core::Size  n_swap_dead_;

private:
	WorkUnitManager *wum_;
};


class WorkUnitManager: public utility::pointer::ReferenceCount {
public:
	friend class WorkUnitQueue_Swapped;

	WorkUnitManager(){
	}

	~WorkUnitManager() override= default;

	// @brief Main loop
	void virtual go()=0;

	void register_work_units( const protocols::wum::WorkUnitList &work_unit_list );

public:
	typedef WorkUnitQueue::iterator       iterator;
	typedef WorkUnitQueue::const_iterator const_iterator;

	WorkUnitQueue& outbound(){ return outbound_wus_; }
	WorkUnitQueue& inbound() { return inbound_wus_; }
	const WorkUnitQueue& outbound() const { return outbound_wus_; }
	const WorkUnitQueue& inbound()  const { return inbound_wus_; }

protected:
	const protocols::wum::WorkUnitList &work_unit_list() const { return work_unit_list_; }
	protocols::wum::WorkUnitList &work_unit_list() { return work_unit_list_; }

protected:
	void write_queues_to_file( const std::string &prefix = "default" ) const;
	void write_work_unit( const WorkUnitBaseOP &wu, std::ostream &out ) const;
	void write_queue( const WorkUnitQueue &the_queue, std::ostream &out ) const;

	void read_queues_from_file( const std::string &prefix = "default" );
	bool read_work_unit( WorkUnitBaseOP &qualified_wu,  std::istream &in );
	void read_queue( WorkUnitQueue &the_queue, std::istream &in  );

	/// @brief return total memory foot print in bytes
	core::Size mem_foot_print() const {
		return inbound().mem_foot_print() + outbound().mem_foot_print();
	};

private:
	protocols::wum::WorkUnitList work_unit_list_;
private:
	WorkUnitQueue inbound_wus_;
	WorkUnitQueue outbound_wus_;


};


}
}

#endif

