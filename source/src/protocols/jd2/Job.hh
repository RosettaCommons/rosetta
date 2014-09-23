// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/Job.hh
/// @brief  header file for Job classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_Job_hh
#define INCLUDED_protocols_jd2_Job_hh

//unit headers
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobOutputterObserver.hh>
// AUTO-REMOVED #include <protocols/jd2/InnerJob.fwd.hh>

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

/// @details The Job class is directly used by the job distributor.  It contains
/// an index into nstruct and a lightweight pointer to an InnerJob (to share
/// those heavier classes).  It also directly contains job-associated output
/// data. The idea here is that your scores, etc are associated with the job,
/// not with a particular pose within that running job.  Thus, the Job object is
/// a better place to store them than the pose.
class Job : public utility::pointer::ReferenceCount
{
public:

	/// TODO these should be maps rather than lists.
	typedef std::list< std::string > Strings;
	typedef std::map< std::string, std::string > StringStringPairs;
	typedef std::map< std::string, core::Real > StringRealPairs;

	Job( InnerJobOP inner_job, core::Size nstruct_index );

	/// @brief returns a copy of this object whose "output fields" are zeroed
	/// out.  Used by the JobDistributor in cases where the job fails and must be
	/// retried to prevent accumulation of Job state after a failure.  This
	/// implementation was chosen over a clear_all_output function to prevent
	/// mover A from deleting mover B's hard work!  You probably should not be
	/// trying to call this function.  The exception:  If you want an
	/// intermediate-output pose (not the final pose) to not have the aggregated
	/// accessory data in the "real" Job object.
	JobOP copy_without_output() const;

	virtual ~Job();

	///@brief Note: only compare if the pointers to the poses are to the
	///same location
	friend
	bool
	operator == ( Job const & a, Job const & b );

	friend
	bool
	operator != ( Job const & a, Job const & b );


	virtual
	void
	show( std::ostream & out ) const;

	friend
	std::ostream &
	operator << ( std::ostream & out, const Job & job );


	///@brief access to inner-job ... use is discouraged - use sparingly!
	/// --- DO NOT use my_job->inner_job()->get_pose()
	/// INSTEAD use my_job->get_pose()
	InnerJobCOP inner_job() const;

	///@brief return the input tag (a short string, generally)
	std::string const & input_tag() const;

	///@brief nonconst access is intended only for the JobInputter to load poses into the InnerJob, and the Parser to add constraints, and the JobDistributor to delete completed inputs (recycle memory)
	InnerJobOP inner_job_nonconst();

	///get_pose : will return
	// 1)  a pose saved in Job-Object ... if available
	// 2)  re-route the call to the JobInputter::pose_from_job

	///@brief return a COP to the input pose
	core::pose::PoseCOP get_pose() const;

	///@brief in-place copy of input pose
	void get_pose( core::pose::Pose& ) const;

	core::Size nstruct_index() const;
	///@brief
	core::Size nstruct_max() const;

	///////////////////////////THIS SECTION OF FUNCTIONS IS MEANT TO BE USED BY MOVERS/////////////////////////////
	//It is safe to call these functions even in the absence of a job distributor - it will store the data in a dummy object.  You are not making your protocol dependent on JD2 by using these functions (although this extra output might get "lost" if you do not emit it by another method like the Tracers.) -- SML 10/20/11
	//functions for loading output info into the job
	///@brief add an output string
	void add_string( std::string const & string_in );

	///@brief add output strings
	void add_strings( Strings const & );

	///@brief add a string/string pair
	void add_string_string_pair( std::string const & string1, std::string const & string2 );

	///@brief add a string/real pair
	void add_string_real_pair( std::string const & string_in, core::Real const real_in );


	////////////////////THIS SECTION OF FUNCTIONS IS FORBIDDEN FOR USE BY MOVERS//////////////////////////////////
	//If Movers try to use these functions, those Movers become tied to JD2 and may fail if JD2 is not present (because there will be no data in these strings).  To prevent this, DO NOT CALL these functions from within Movers. SML 10/20/11

	//functions for returning output info from the job.  You get iterators so that this interface can stay constant as the underlying implementation changes
	Strings::const_iterator output_strings_begin() const;
	Strings::const_iterator output_strings_end() const;

	StringStringPairs::const_iterator output_string_string_pairs_begin() const;
	StringStringPairs::const_iterator output_string_string_pairs_end() const;

	StringRealPairs::const_iterator output_string_real_pairs_begin() const;
	StringRealPairs::const_iterator output_string_real_pairs_end() const;

	////////////////////////END SECTION//////////////////////////////////////////////////////////////////

	//there are no functions for deleting output info.  This is on purpose - use copy_without_output instead

	Strings get_strings() const { return long_strings_; }

	StringStringPairs get_string_string_pairs() const { return string_string_pairs_; };

	StringRealPairs get_string_real_pairs() const { return string_real_pairs_; };

	void set_status_prefix( std::string prefix ) {
		status_prefix_ = prefix;
	}

	std::string const& status_prefix() const {
		return status_prefix_;
	}

	bool completed() const {
		return completed_;
	}

	bool to_do() const {
		return !completed_ && !bad();
	}

	bool bad() const;

	void set_completed(bool value = true)  {
		completed_ = value;
	}

	void set_bad(bool value = true);

	void add_output_observer( JobOutputterObserverAP an_observer );
	void remove_output_observer( JobOutputterObserverAP old_observer );
	void call_output_observers( core::pose::Pose const & pose );

private:
	//InnerJobCOP inner_job() const;
	//bookkeeping data
	///@brief a pointer to the "heavy" InnerJob which maintains the starting pose for the job (shared across nstruct)
	InnerJobOP inner_job_;
	///@brief which nstruct is this?
	core::Size const nstruct_index_;

	std::string status_prefix_;

	/// the following block of data makes the Job class pretty heavy. If we create thousands of Jobs this can be a problem...
	// put these into an extra class and store a pointer ... NULL if nothing has been stored yet?
	// cleanup after job is written to output... OL 6/2/09
	//Storage units for output data
	///@brief used for arbitrary string data (stuff you've preformatted).  Intended to be appended to the end of a PDB or dumped to a tracer if not in a PDB output mode.
	Strings long_strings_;
	///@brief string-string pairs.  Inserted into SCORE: lines in scorefiles/silentfiles.
	StringStringPairs string_string_pairs_;
	///@brief string-real pairs (scoretype/score pairs).  Inserted into SCORE: lines in scorefiles/silentfiles
	StringRealPairs string_real_pairs_;

	///@brief container of evaluators
	//Oliver??

	bool completed_;

#ifndef PTR_MODERN
	typedef std::set< JobOutputterObserverAP > JobOutputterObservers;
#else
	typedef std::set< JobOutputterObserverAP, utility::pointer::owner_less< JobOutputterObserverAP > > JobOutputterObservers;
#endif
	JobOutputterObservers output_observers_;

}; // Job

extern JobCOP const JD2_BOGUS_JOB;

} // namespace jd2
} // namespace protocols



#endif //INCLUDED_protocols_jd2_Job_HH
