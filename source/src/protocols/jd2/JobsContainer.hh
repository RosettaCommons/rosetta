// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobsContainer.hh
/// @brief  header file for Job classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_JobsContainer_hh
#define INCLUDED_protocols_jd2_JobsContainer_hh

//unit headers
#include <protocols/jd2/JobsContainer.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//C++ headers
#include <string>
#include <list>
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
	namespace jd2 {

		/// @details The JobsContainer class contains a list of JobsOPs.
		/// It permits the list to be versatile -- for example, if one wishes
		/// to load only a subset of jobs into memory at any given time.
		class JobsContainer : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< JobsContainer >
		{
		public:

			/// @brief JobsContainer constructor.
			///
			JobsContainer();

			/// @brief JobsContainer destructor.
			///
			virtual ~JobsContainer();
			
			/// @brief Get self const owning pointers.
			///
			inline JobsContainerCOP get_self_ptr() const { return shared_from_this(); }

			/// @brief Get self owning pointers.
			///
			inline JobsContainerOP get_self_ptr() { return shared_from_this(); }
			
			/// @brief Get a specific job, by number.
			/// @details Should work even if jobs have been deleted, since this uses a
			/// map instead of an array.
			JobOP operator [](core::Size const index);
		
			/// @brief Get a specific job, by number (const-access).
			/// @details Should work even if jobs have been deleted, since this uses a
			/// map instead of an array.
			JobCOP operator [](core::Size const index) const;
			
			/// @brief Assignment operator.
			///
			JobsContainer &
			operator=( JobsContainer const & src );
			
			/// @brief Get the total number of jobs.
			/// @details Might not be the number in the joblist_ map, if we're only holding
			/// a subset in memory at any given time.
			core::Size size() const { return total_jobs_; }
			
			/// @brief Add a job to the list of jobs.
			///
			void push_back( JobOP new_job );
			
			/// @brief Set the total number of jobs.
			/// @details This overrides whatever the length of the joblist_ is.
			void set_total_jobs( core::Size const total_jobs );
			
			/// @brief Clear the jobs list
			///
			void clear();
			
			/// @brief Access the last element.
			///
			JobOP back();
			
			/// @brief Randomize the order of elements (the map keys)
			///
			void shuffle();
			
			/// @brief Erase an element in the jobs list
			///
			void erase( core::Size const index );
			
			/// @brief Does the job with the given index exist in the currently-loaded list of jobs?
			///
			bool has_job( core::Size const index ) const;
			
			/// @brief Can the job with the given index be deleted?
			///
			bool can_be_deleted ( core::Size const index ) const;
			
			/// @brief Return the index of the highest job currently loaded in memory.
			/// 
			core::Size highest_job_index() const { return highest_job_index_; }
			
			/// @brief Set the job inputter.
			/// @details Needed if the job list is to be updated by the job inputter.
			void set_job_inputter( JobInputterOP jobinputter ) { job_inputter_ = jobinputter; return; }

			/// @brief Get a list of job indices currently in memory.
			/// @details The output vector is cleared and populated with the current job indices in memory.
			void get_loaded_job_indices( utility::vector1 < core::Size > &output ) const;
			
			/// @brief Mark all jobs currently in memory as deletable.
			/// @details This will result in all jobs being purged the next time a higher-index job is requeted.
			void set_all_current_jobs_as_deletable();
			
			/// @brief Set whether or not this JobsContainer forces job purging when higher-numbered jobs are requested.
			/// @details If true, and if a job inputter is used that limits the number of jobs in memory at any given time,
			/// then when a higher job index is requested than exists in memory, the JobsContainer will force job purging
			/// until that job exists in memory.  If false (the default), then only deletable jobs will be purged, and the
			/// job inputter can attempt to add jobs up to the requested job (which might fail, if the requested job's number
			/// is too high).
			void set_force_job_purging( bool const val=true ) { force_job_purging_ = val; return; }
			
			/// @brief Get whether or not this JobsContainer forces job purging when higher-numbered jobs are requested.
			/// @details If true, and if a job inputter is used that limits the number of jobs in memory at any given time,
			/// then when a higher job index is requested than exists in memory, the JobsContainer will force job purging
			/// until that job exists in memory.  If false (the default), then only deletable jobs will be purged, and the
			/// job inputter can attempt to add jobs up to the requested job (which might fail, if the requested job's number
			/// is too high).
			bool force_job_purging() const { return force_job_purging_; }
			
		private:
		
			/// @brief The list of owning pointers to jobs, with the job number as the key.
			///
			std::map < core::Size, JobOP > joblist_;
			
			/// @brief The highest-index job in the joblist_ map.
			/// @details Starts at 0 and counts up as elements are added.
			core::Size highest_job_index_;
			
			/// @brief The total number of jobs.
			/// @details  By default, this will be the size of the joblist_, unless overridden.
			core::Size total_jobs_;
			
			/// @brief Has the total number of jobs been overridden?
			/// @details Default false.
			bool total_jobs_set_;
			
			/// @brief The JobInputter that populated this object (or which is otherwise associated with it).
			///
			JobInputterOP job_inputter_;
			
			/// @brief If true, and if a job inputter is used that limits the number of jobs in memory at any given time,
			/// then when a higher job index is requested than exists in memory, the JobsContainer will force job purging
			/// until that job exists in memory.  If false (the default), then only deletable jobs will be purged, and the
			/// job inputter can attempt to add jobs up to the requested job (which might fail, if the requested job's number
			/// is too high).
			bool force_job_purging_;

		}; // JobsContainer class


	} // namespace jd2
} // namespace protocols


#endif //INCLUDED_protocols_jd2_JobsContainer_hh
