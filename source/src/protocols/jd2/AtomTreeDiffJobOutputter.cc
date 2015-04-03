// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/AtomTreeDiffJobOutputter.cc
/// @brief  AtomTreeDiffJobOutputter
/// @author Gordon Lemmon (gordon.h.lemmon@vanderbilt.edu); Rocco Moretti (rmoretti@u.washington.edu)

///Unit headers
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/Job.fwd.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <core/svn_version.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>


///C++ headers
//#include <string>
//#include <sstream>

static thread_local basic::Tracer TR( "protocols.jd2.AtomTreeDiffJobOutputter" );

namespace protocols {
namespace jd2 {

using namespace basic::options;
using namespace basic::options::OptionKeys;

AtomTreeDiffJobOutputter::AtomTreeDiffJobOutputter():
	bb_precision_(6),
	sc_precision_(4),
	bondlen_precision_(2),
	use_input_(false)
{

	set_precision( option[ out::file::atom_tree_diff_bb ], option[ out::file::atom_tree_diff_sc ], option[ out::file::atom_tree_diff_bl ] );

	// Add directory path and prefix/suffix (if specified) to plain file name:
	outfile_name_= basic::options::option[ OptionKeys::out::file::atom_tree_diff ]();
	{
		/// TODO abstract this section to a generic file for all job outputters
		utility::file::FileName outfile(outfile_name_);
		std::ostringstream oss;
		oss << basic::options::option[ basic::options::OptionKeys::out::prefix ]() << outfile.base()
		    << basic::options::option[ basic::options::OptionKeys::out::suffix ]();
		outfile.base( oss.str() );
		if ( basic::options::option[ basic::options::OptionKeys::out::path::pdb ].user() ){
			outfile.path(basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path());
			outfile.vol(basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().vol());
		} else if( basic::options::option[ basic::options::OptionKeys::out::path::all ].user() ){
			outfile.path(basic::options::option[ basic::options::OptionKeys::out::path::all ]().path());
			outfile.vol(basic::options::option[ basic::options::OptionKeys::out::path::all ]().vol());
		}else{
			outfile.path("");
			outfile.vol("");
		}

		if( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ] && outfile.ext() != "gz" ) {
			outfile.ext( ".gz" ); // else use input extension
		}
		outfile_name_ = outfile.name();
	}

	//Looks like we want to only have the file open when we're actually writing, moved this to dump_pose

	/*
	if( utility::file::file_exists(outfile_name) ) {
		atom_tree_diff_.read_file(outfile_name);
		out_.open_append( outfile_name.c_str() );
	} else {
		out_.open( outfile_name.c_str() );
		if( basic::options::option[ basic::options::OptionKeys::run::version ] ) {
			out_ << "# Mini-Rosetta version " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url() << "\n";
		}
	}
	if ( !out_.good() ) {
		utility_exit_with_message( "Unable to open file: " + outfile_name + "\n" );
	}
	*/
}

AtomTreeDiffJobOutputter::~AtomTreeDiffJobOutputter(){}

void
AtomTreeDiffJobOutputter::final_pose( JobOP job, core::pose::Pose const & pose, std::string const & /*tag*/ ){

	call_output_observers( pose, job );

	std::map< std::string, core::Real > scores;
	for( Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end());
			 it != end; ++it) {
		// Need to use operator[], as we want later (more recent) items in the string_real_pairs list to overwrite the earlier ones
		// (e.g. if we've recomputed a value that was read in from the inital silent file.)
		// Standard constructors/copy functions ignores items which are already in map.
		scores[it->first] = it->second;
	}

	std::string output= output_name(job);
	dump_pose( output, pose, scores, job);

}

/// @brief Appends pose to atom tree diff file
void
AtomTreeDiffJobOutputter::dump_pose(
	std::string const & tag,
	core::pose::Pose const & pose,
	std::map< std::string, core::Real > scores,
	JobCOP job
){

	if( utility::file::file_exists(outfile_name_) ) {
		//atom_tree_diff_.read_file(outfile_name_);
		out_.open_append( outfile_name_.c_str() );
	} else {
		out_.open( outfile_name_.c_str() );
		if( basic::options::option[ basic::options::OptionKeys::run::version ] ) {
			out_ << "# Mini-Rosetta version " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url() << "\n";
		}
	}
	if ( !out_.good() ) {
		utility_exit_with_message( "Unable to open file: " + outfile_name_ + "\n" );
	}

	if ( use_input_ && job ) {
		// Use the input pose as the reference pose, matching the original ligand docking application.
		std::string const ref_tag( job->input_tag() );
		TR<< "ref_tag: "<< ref_tag<< std::endl;
		if( last_ref_tag_ != ref_tag ){
			TR<< "writing reference pose"<< std::endl;
			std::map< std::string, core::Real > temp_scores;
			temp_scores["is_reference_pose"]= 1;
			job->get_pose(last_ref_pose_); // deep copy for reference pose
			last_ref_tag_= ref_tag; // deep copy for reference pose
			// Yes, we're outputting the input pose, but using the output tag -- we need the extra numbers on the end so we can slice them off when reading
			core::import_pose::atom_tree_diffs::dump_reference_pose(out_, "%REF%_"+tag, temp_scores, last_ref_pose_);
		}
	} else {
		// Use the current structure as the reference pose, if we don't already have an applicable one.
		core::Size end= tag.find_last_of('_');
		std::string const ref_tag= tag.substr(0, end);
		TR<< "ref_tag: "<< ref_tag<< std::endl;
		if( last_ref_tag_ != ref_tag ){
			TR<< "writing reference pose"<< std::endl;
			std::map< std::string, core::Real > temp_scores;
			temp_scores["is_reference_pose"]= 1;
			core::import_pose::atom_tree_diffs::dump_reference_pose(out_, "%REF%_"+tag, temp_scores, pose);
			last_ref_tag_= ref_tag;// deep copy for reference pose
			last_ref_pose_= pose;// deep copy for reference pose
		}
	}
	if( used_tags_.find(tag) != used_tags_.end() ){
		basic::Error() << "Tag " << tag << " already exists in silent file; writing structure anyway..." << std::endl;
	}
	core::import_pose::atom_tree_diffs::dump_atom_tree_diff(out_, tag, scores, last_ref_pose_, pose, bb_precision_, sc_precision_, bondlen_precision_);
	used_tags_.insert(tag);
	// Can't flush compressed streams -- results in file truncation
	if ( out_.uncompressed() ) out_.flush();
	out_.close();
}


void
AtomTreeDiffJobOutputter::other_pose(
  JobOP job,
	core::pose::Pose const & pose,
	std::string const & /*tag*/,
	int /*copy_count*/, /*default -1 */
	bool /*score_only*/ /*default false*/
 ){
	call_output_observers( pose, job );
	// do something with this function later
	return;
}

bool
AtomTreeDiffJobOutputter::job_has_completed( JobCOP job ){
	// was the job completed beforehand ( in the silent file before this app even started) ?
	if( used_tags_.find( output_name( job )) != used_tags_.end() ) return true;
	// did WE complete the job later ?
	if( job->completed() ) return true;

	return false;
}

std::string
AtomTreeDiffJobOutputter::output_name( JobCOP job ){
	return affixed_numbered_name( job );
}

void
AtomTreeDiffJobOutputter::set_precision(int bb_precision, int sc_precision, int bondlen_precision ) {
  bb_precision_ = bb_precision;
  sc_precision_ = sc_precision;
  bondlen_precision_ = bondlen_precision;
}

/// @brief use input as reference pose?
void
AtomTreeDiffJobOutputter::use_input_for_ref(bool use_input) {
	// Reset reference tag, in case we're changing course midstream
	last_ref_tag_ = "";
	use_input_ = use_input;
}

//CREATOR SECTION
std::string
AtomTreeDiffJobOutputterCreator::keyname() const
{
        return "AtomTreeDiffJobOutputter";
}

protocols::jd2::JobOutputterOP
AtomTreeDiffJobOutputterCreator::create_JobOutputter() const {
        return protocols::jd2::JobOutputterOP( new AtomTreeDiffJobOutputter );
}

}//jd2
}//protocols
