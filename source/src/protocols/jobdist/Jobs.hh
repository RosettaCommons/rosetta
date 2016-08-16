// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jobdist/Jobs.hh
///
/// @brief  Objects representing inputs for the JobDistributor
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_jobdist_Jobs_hh
#define INCLUDED_protocols_jobdist_Jobs_hh

#include <protocols/jobdist/Jobs.fwd.hh>

#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <sstream>


namespace protocols {
namespace jobdist {

static THREAD_LOCAL basic::Tracer jobs_tracer( "protocol.jobdist.jobs.hh" );

/// @brief Each Job object describes a particular input to Rosetta.
///
/// @details
/// Ordinarily, an input is just a single PDB file.
/// In other cases, e.g. docking, input might be a pair of PDB files.
/// In still other cases, input might be a list of tags to extract from a silent file.
/// And so on.
/// One BasicJob class is defined below that should meet most needs,
/// but if you need to carry additional information, you can subclass it.
///
/// Contains a map of strings to strings that can be used to store arbitrary
/// extra data if you're too lazy to make a subclass.
/// Templated get/set functions accept anything that can be
/// read from and written to a C++ stream with >> and << operators.
/// All data is stored as strings internally.
/// Some numeric data may experience slight round-off error,
/// though I think we're storing enough digits to avoid that.
///
/// Example:
///   BasicJob jobdata("1abc.pdb", 10);
///   jobdata.set("cmp_to", "1xyz.pdb");
///   jobdata.set("score", 42.42);
///
///   std::string cmp_to;
///   core::Real score;
///   if( jobdata.get("cmp_to", cmp_to) ) std::cout << "cmp_to = " << cmp_to << std::endl;
///   else std::cout << "Unable to recover cmp_to!" << std::endl;
///   if( jobdata.get("score", score) ) std::cout << "score = " << score << std::endl;
///   else std::cout << "Unable to recover score!" << std::endl;
///
/// Although Jobs are small objects, it is common for Rosetta to loop many times over the same input.
/// In order to avoid creating (size of -s and -l) * (nstruct) JobData objects,
/// each Job also has an nstruct counter.
///
/// So, each Job represents a unique input,
/// and each tuple of (Job, struct_n) represents a unique output.
///
class BasicJob : public utility::pointer::ReferenceCount
{
public:

	/// @brief You MUST ensure that input_tag is a UNIQUE identifier for this Job!
	BasicJob(std::string input_tag, std::string native_tag, int nstruct=1):
		input_id_(input_tag),
		native_id_(native_tag),
		nstruct_(nstruct),
		preserve_whole_input_tag_(false)
	{}

	virtual ~BasicJob() {}

	/// @brief The number of times this job should be repeated.
	virtual int nstruct() const
	{ return nstruct_; }

	/// @brief The tag supplied at create time.
	virtual std::string input_tag() const
	{ return input_id_; }

	/// @brief The tag supplied at create time.
	virtual std::string native_tag() const
	{ return native_id_; }

	/// @brief The tag supplied at create time plus the nstruct number (from 1 to nstruct inclusive).
	virtual std::string output_tag(int struct_n) const;

	/// @brief The tag supplied at create time plus the nstruct number (from 1 to nstruct inclusive).
	virtual std::string output_file_name() const {
		return output_file_name_;
	};

	/// @brief The tag supplied at create time plus the nstruct number (from 1 to nstruct inclusive).
	void set_output_file_name( std::string str ) {
		output_file_name_ = str;
	};

	/// @brief Extracts named value.  Works for anything that deserializes from string.  Returns false on error.
	template <typename T>
	bool get(std::string const & key, T & value)
	{
		if ( extra_data_.find(key) == extra_data_.end() ) return false;
		std::istringstream is( extra_data_[key] );
		is >> value;
		return !is.fail();
	}

	/// @brief Specialization for strings.
	bool get(std::string const & key, std::string & value)
	{
		if ( extra_data_.find(key) == extra_data_.end() ) return false;
		value = extra_data_[key];
		return true;
	}

	/// @brief Set named value.  Works for anything that serializes to a string.
	template <typename T>
	void set(std::string const & key, T const & value)
	{
		std::ostringstream os;
		os.precision(20); // More than enough even for doubles (?)
		os << value;
		extra_data_[key] = os.str();
	}

	/// @brief Specialization for strings.
	void set(std::string const & key, std::string const & value)
	{ extra_data_[key] = value; }

	void set_preserve_whole_input_tag( bool setting ){ preserve_whole_input_tag_ = setting; }
protected:

	std::string input_id_;
	std::string native_id_; /// name of the native
	int nstruct_;
	std::map< std::string, std::string > extra_data_;
	std::string output_file_name_;
	bool preserve_whole_input_tag_;
}; // BasicJob


} // namespace jobdist
} // namespace protocols

#endif // INCLUDED_protocols_jobdist_Jobs_HH
