// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/SilentFileJobOutputter.hh
/// @brief  header file for SilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Oliver Lange olange@u.washington.edu


#ifndef INCLUDED_protocols_jd2_SilentFileJobOutputter_hh
#define INCLUDED_protocols_jd2_SilentFileJobOutputter_hh

//unit headers
#include <protocols/jd2/SilentFileJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//utility headers
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

//C++ headers
#include <string>

namespace protocols {
namespace jd2 {

/// @details this is a implementation of JobOutputter for silent-file-based output.
// todo:
// read designated silent-file in beginning and implement method
// job_completed()
class SilentFileJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

	typedef protocols::jd2::FileJobOutputter parent;

	SilentFileJobOutputter();
	virtual ~SilentFileJobOutputter();

	/// @brief this function flushes any internal buffers - see parent class for explanation
	virtual void flush();

	//////////////////////////////creating output functions/////////////////////////////////////////
	/// @brief this function outputs the final result of a job.
	virtual
	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag );

	/// @brief this function is intended for saving
	/// mid-protocol poses; for example the final centroid
	/// structure in a combined centroid/fullatom protocol.
	/// --->these go to file silent_filename+tag
	virtual
	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it
	/// belongs here since it needs to check the same output
	/// locations as the class normally writes to.  This class
	/// checks wherever output goes to see if the job's
	/// expected output already exists (on disk or whatever).
	/// This is the most basic form of checkpointing.
	virtual
	bool job_has_completed( JobCOP job );

public: // accessors

	/// @brief this is the master function for determining the
	/// unique output identifier for a job
	virtual
	std::string output_name( JobCOP job );

	virtual
	std::string filename( JobCOP ) const {
		return silent_file_;
	}

	void
	set_silent_file_name( utility::file::FileName name );

	void
	set_forced_silent_struct_type( std::string const& );

	void
	set_write_separate_scorefile( bool write_separate_scorefile );

	void set_write_no_structures( bool value = true ) {
		bWriteNoStructures_ = value;
		bWriteIntermediateStructures_ = !value && bWriteIntermediateStructures_;
	}

	//////////////////////////////////////scorefile functions/////////////////////////////////////
protected:
	//called by final_- and other_pose methods
	core::io::silent::SilentStructOP dump_pose(
		utility::file::FileName const & filename,
		JobCOP job,
		core::pose::Pose const & pose,
		bool bWriteScoreOnly,
		int copy_count = -1, /* if 0 or positive attach as postfix to job-tag ONLY used in other_pose output*/
		std::string const & suffix = "" /* appended to tag */
	);

private: // methods

	void add_silent_struct(
		core::io::silent::SilentStructOP ss,
		utility::file::FileName const & fn
	);

	/// @brief called by the constructor to set filename and options
	void set_defaults();

	void read_done_jobs();

	/// @brief this function is called at the end of job distribution to flush the buffers
	void write_all_structs();

private: // members

	// write intermediate files ( from calls to other_pose )
	bool bWriteIntermediateFiles_;

	// write also structural information to intermediate files
	// ( from calls to other_pose )
	bool bWriteIntermediateStructures_;

	/// @brief toggle to switch off writing of structures
	bool bWriteNoStructures_;

	/// @brief whether to write a separate scorefile that contains
	///the scorelines from the silent file
	bool write_separate_scorefile_;

	//  file name for silent-output
	utility::file::FileName silent_file_;

	// list of tags already written
	utility::vector1< std::string > silent_file_tags_;

	// used for buffering output.
	core::Size n_to_buffer_;

	// when n_to_buffer_ is reached write all structures with a probability of ... (default 1)
	core::Real random_flush_frequency_;

	//utility::vector1< core::io::silent::SilentStructOP > saved_structs_;
	utility::vector1< std::pair< core::io::silent::SilentStructOP, utility::file::FileName > > saved_structs_;

	/// @brief override choice of silent_struct_type
	std::string forced_silent_struct_type_;
}; // SilentFileJobOutputter

// This is necessary because some protocols append prefixes
// to the tags and thus a simple string comparison will not
// recognise that S_1234_4 and C_S_1234_4 are the same.
class CompareTags: public std::unary_function<std::string, bool > {
public:
	CompareTags( const std::string & querytag ): querytag_(querytag) {}

	bool operator () ( const std::string & compare_tag ) const {
		// Strings match if all the characters of the shorter string match all of the last characters of the other.
		core::Size offset1 = compare_tag.find("S_");
		if ( offset1 == std::string::npos ) offset1 = 0;
		core::Size offset2 = querytag_.find("S_");
		if ( offset2 == std::string::npos ) offset2 = 0;

		return ( compare_tag.substr(offset1) == querytag_.substr(offset2) );
	}
private:
	const std::string querytag_;
};

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_SilentFileJobOutputter_HH
