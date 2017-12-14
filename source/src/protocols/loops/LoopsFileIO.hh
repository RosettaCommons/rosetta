// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopsFileIO.hh
/// @brief This class exists to handle the reading and writing of loops files.
/// @author Brian D. Weitzner

#ifndef INCLUDED_protocols_loops_LoopsFileIO_HH
#define INCLUDED_protocols_loops_LoopsFileIO_HH

// Unit header
#include <protocols/loops/LoopsFileIO.fwd.hh>

// Package headers
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/ResidueIndexDescription.hh>

// Utility headers
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <fstream>
#include <string>

namespace protocols {
namespace loops {

/// @brief Checks if there is a problem with the beginning and ending residues defined
/// in a loops file.
void
validate_loop_start_stop(
	bool prohibit_single_residue_loops,
	core::Size start,
	core::Size stop,
	std::string const & filename,
	core::Size linecount
);


class LoopFromFileData
{
public:
	LoopFromFileData();

	LoopFromFileData(
		core::pose::ResidueIndexDescriptionFromFile const & start_res,
		core::pose::ResidueIndexDescriptionFromFile const & cutpoint_res,
		core::pose::ResidueIndexDescriptionFromFile const & end_res,
		core::Real skip_rate,
		bool extended,
		bool prohibit_single_residue_loops = true
	);

	/// constructed the other way around (for the the PoseNumberedLoopReader)
	LoopFromFileData(
		SerializedLoop const & loop,
		std::string const & fname, // file from which this loop was created
		bool prohibit_single_residue_loops = true
	);

	/// @brief loop-index resolution function: construct a SerializedLoopData object
	/// by possibly retrieving data from a Pose.  This function also performs the
	/// loop-index checks performed by the PoseNumberedLoopFileReader.
	SerializedLoop
	resolve_as_serialized_loop_from_pose( core::pose::Pose const & pose ) const;

	core::pose::ResidueIndexDescriptionFromFile const & start_res() const { return start_res_; }
	void start_res( core::pose::ResidueIndexDescriptionFromFile const & setting ) { start_res_ = setting; }

	core::pose::ResidueIndexDescriptionFromFile const & cutpoint_res() const { return cutpoint_res_; }
	void cutpoint_res( core::pose::ResidueIndexDescriptionFromFile const & setting ) { cutpoint_res_ = setting; }

	core::pose::ResidueIndexDescriptionFromFile const & end_res() const { return end_res_; }
	void end_res( core::pose::ResidueIndexDescriptionFromFile const & setting ) { end_res_ = setting; }

	core::Real skip_rate() const { return skip_rate_; }
	void skip_rate( core::Real setting ) { skip_rate_ = setting; }

	bool extended() const { return extended_; }
	void extended( bool setting ) { extended_ = setting; }

	bool prohibit_single_residue_loops() const { return prohibit_single_residue_loops_; }
	void prohibit_single_residue_loops( bool setting ) { prohibit_single_residue_loops_ = setting; }

private:
	core::pose::ResidueIndexDescriptionFromFile start_res_;
	core::pose::ResidueIndexDescriptionFromFile cutpoint_res_;
	core::pose::ResidueIndexDescriptionFromFile end_res_;
	core::Real skip_rate_;
	bool extended_;
	bool prohibit_single_residue_loops_;
};

class LoopsFileData : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~LoopsFileData() override;

	LoopsOP resolve_loops( core::pose::Pose const & pose ) const;
	SerializedLoopList resolve_as_serialized_loops( core::pose::Pose const & pose ) const;
	core::Size size() const;
	void resize( core::Size new_size );
	void push_back( LoopFromFileData const & loop );
	void insert_loop_at_index( LoopFromFileData const & loop, core::Size i );

	LoopFromFileData const & operator[] ( core::Size const i ) const;

private:
	utility::vector1< LoopFromFileData > loops_file_data_;
};

/// @brief This class ensures that the Loops object that is needed to run any of the various
/// forms of loop modeling is correctly initialized from a Pose.  If the residues specified
/// from a loops file have not been resolved into the residue indices for a Pose, then
/// this class will die with an assertion failure.
class GuardedLoopsFromFile : public utility::pointer::ReferenceCount
{
public:
	/// @brief default ctor; sets the object in an "in charge" state.
	GuardedLoopsFromFile();

	/// @brief constructor from a loops-file-data object: sets this object in an "in charge" state.
	GuardedLoopsFromFile( LoopsFileData const & lfd );

	/// @brief constructor from loops pointer: sets this object in a "not in charge" state.
	GuardedLoopsFromFile( LoopsOP loops );

	/// @brief constructor from a loops object: set this object in an "in charge" state.
	GuardedLoopsFromFile( Loops const & loops );

	/// @brief copy constructor; takes it's "in charge" state from src.
	GuardedLoopsFromFile( GuardedLoopsFromFile const & src );

	/// @brief copy constructor; takes it's "in charge" state from src. Set the copy as not-in-charge.
	GuardedLoopsFromFile( GuardedLoopsFromFile const & src, bool  );

	/// @brief virtual dstor
	~GuardedLoopsFromFile() override;

	/// @brief assignment operator; takes it's "in charge" state from rhs
	GuardedLoopsFromFile & operator = ( GuardedLoopsFromFile const & rhs );

	/// @brief set to "in charge" state.
	void in_charge( bool setting );

	/// @brief get "in charge" state.
	bool in_charge() const;

	/// @brief Resolve the loop indices, and mark the state as resolved, so that calls to loops() will succeed.
	/// This function will re-resolve loop indices with a new pose, which may be important if the same loop_file_data_
	/// is being applied to a pose which has different PDB indices.  If I am not in charge, this is a no-op.
	void resolve_loop_indices( core::pose::Pose const & );

	//// @brief Resolve the loop indices if they have not yet been resolved.  If I am not in charge, this is a no-op.
	void resolve_loop_indices_once( core::pose::Pose const & );

	/// @brief request the LoopsCOP pointer; asserts that the loop indices
	/// have been resolved or that "I am not in charge".
	LoopsCOP loops() const;

	/// @brief request the LoopsOP pointer; asserts that the loop indices
	/// have been resolved or that "I am not in charge".
	LoopsOP loops();

	/// @brief set the loops owning pointer object directly
	void set_loops_pointer( LoopsOP setting );

	/// @brief set the loops to copy the contents of settings into the existing Loops object
	void loops( Loops const & setting );

	/// @brief set the LoopsFileData object directly
	void loops( LoopsFileData const & setting );

	/// @brief read access to the LoopsFileData
	LoopsFileData const & loops_file_data() const;

private:
	bool in_charge_;
	bool pose_has_resolved_loop_indices_;
	bool rely_on_loopfile_indices_;
	LoopsFileData loops_file_data_;
	LoopsOP loops_;
};

///////////////////////////////////////////////////////////////////////////
// a class for reading a loops file where the format of that file might
// be one of many formats
class LoopsFileIO : public utility::pointer::ReferenceCount {

public:

	//constructor
	LoopsFileIO();

	//copy constructor
	LoopsFileIO( const LoopsFileIO & src );

	// destructor
	~LoopsFileIO() override;

	friend std::ostream & operator<<( std::ostream & os, const LoopsFileIO & loops_file_io );

	/// @brief Return an "unresolved" list of loops specified in a file
	/// which can be turned into a "resolved" list of residue indices in
	/// a particular Pose by giving each LoopFromFileData object access
	/// to that Pose.
	/// Note: prohibit_single_residue_loops used to be called "strict_looprelax_checks_"
	/// which was decidedly opaque
	LoopsFileDataOP read_loop_file(
		std::string const & filename,
		bool prohibit_single_residue_loops = true
	);

	LoopsFileDataOP read_loop_file_stream(
		std::istream & loopfstream,
		std::string const & filename,
		bool prohibit_single_residue_loops = true
	);


}; // LoopsFileIO


/// This is the main legacy loop-reading function, which will read the pose-numbered
/// file.  This functionality is used by a great many number of places often having
/// nothing to do with representing a loop, so it will persist.
class PoseNumberedLoopFileReader {
public:
	PoseNumberedLoopFileReader();

	SerializedLoopList
	read_pose_numbered_loops_file(
		std::istream & is,
		std::string const & filename,
		bool strict_looprelax_checks = true // i.e. prohibit_single_residue_loops
	);

	/// @brief if the input stream has had some number of lines already removed from it,
	/// indicate how many.
	void set_linecount_offset( core::Size );

	/// @brief For code that relys on reading loop-file-formatted ranges if residues
	/// but which really ought to use
	void hijack_loop_reading_code_set_loop_line_begin_token( std::string const & token );

private:
	std::string loop_line_begin_token_; // usually "LOOP"
	core::Size linecount_offset_;

};

/// The following enumerators are used for parsing JSON formatted loop files
enum ResidueIdentifier {
	start=1,
	stop,
	cut_point,
	number_of_residue_identifiers=cut_point
};

enum LoopConfiguration {
	extras=number_of_residue_identifiers + 1,
	resSeq,
	iCode,
	chainID,
	skip_rate,
	extend,
	use_pose_numbering,
	number_of_configuration_keywords=use_pose_numbering
};

class JSONFormattedLoopsFileReader {
public:

	/// @brief if the input stream has had some number of lines already removed from it,
	/// indicate how many.
	void set_linecount_offset( core::Size );

	LoopsFileDataOP
	read_loop_file(
		std::istream & is,
		std::string const & filename,
		bool prohibit_single_residue_loops
	);

private: // methods
	LoopsFileDataOP
	parse_json_formatted_data(
		utility::json_spirit::mValue & json_data,
		bool prohibit_single_residue_loops,
		std::string const & filename
	);

	void
	ensure_all_fields_are_valid( utility::json_spirit::mValue & json_data, std::string const & filename );

	core::pose::ResidueIndexDescriptionFromFile
	parse_json_residue_info(
		utility::json_spirit::mValue & json_loop_data,
		ResidueIdentifier residue_identifier,
		std::string const & filename,
		core::Size & approximate_linenumber
	);

	void
	parse_configuration_options(
		utility::json_spirit::mValue & json_loop_data,
		LoopFromFileData & loop
	);

	void setup_residue_type_map();
	std::string name_from_residue_identifier( ResidueIdentifier residue_identifier );
	std::string name_from_loop_configuration( LoopConfiguration loop_configuration );

private:
	core::Size linecount_offset_;

	static bool initialized_;
	static utility::vector1< std::string > valid_loop_file_keys_;

};

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopsFileIO_HH
