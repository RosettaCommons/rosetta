// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CheckPointer
/// @brief Application-level code for Abrelax, fold_cst and JumpingFoldCst protocols
/// @details
///    use -help to see options
///    usage of class:
///
/// @author Mike Tyka
/// @author Oliver Lange

#ifndef INCLUDED_protocols_checkpoint_CheckPointer_hh
#define INCLUDED_protocols_checkpoint_CheckPointer_hh

// Unit Headers

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#ifdef WIN32
#include <protocols/moves/MonteCarlo.hh>
#endif

#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>
#include <ostream>


namespace protocols {
namespace checkpoint {

class FileBuffer {
public:
	FileBuffer( std::string const & filename, bool gzipped = false ):
		filename_( filename ),
		gzipped_( gzipped )
	{
	}


	void set_contents(const std::string &contents ){ contents_ = contents; }

	void dump();

	core::Size size(){ return contents_.length(); }

private:
	std::string filename_;
	bool gzipped_;
	std::string contents_;
};


class CheckPointer : public utility::pointer::ReferenceCount {
public:
	CheckPointer( std::string const& type );

	virtual ~CheckPointer() {
		clear_checkpoints();
	}

	void clear_checkpoints();

	void flush_checkpoints();

	bool recover_checkpoint(
		core::pose::Pose &pose,
		moves::MonteCarlo *mc,
		std::string const& current_tag,
		std::string const& id,
		bool fullatom = false,
		bool foldtree = false
	);

	bool recover_checkpoint(
		core::pose::Pose &pose,
		moves::MonteCarloOP mc,
		std::string const& current_tag,
		std::string const& id,
		bool fullatom = false,
		bool foldtree = false
	)
	{
		return recover_checkpoint( pose, &(*mc), current_tag, id, fullatom, foldtree );
	}

	bool recover_checkpoint(
		core::pose::Pose &pose,
		std::string const& current_tag,
		std::string const& id,
		bool fullatom = false,
		bool foldtree = false
	)
	{
		return recover_checkpoint( pose, NULL, current_tag, id, fullatom, foldtree );
	}

	void checkpoint(
		core::pose::Pose &pose,
		moves::MonteCarlo *mc,
		std::string const& current_tag,
		std::string const& id,
		bool foldtree = false
	);

	void checkpoint(
		core::pose::Pose &pose,
		moves::MonteCarloOP mc,
		std::string const& current_tag,
		std::string const& id,
		bool foldtree = false
	){
		checkpoint( pose, &(*mc), current_tag, id, foldtree );
	}

	void checkpoint(
		core::pose::Pose &pose,
		std::string const& current_tag,
		std::string const& id,
		bool foldtree = false
	){
		checkpoint( pose, NULL, current_tag, id, foldtree );
	}


	std::string const& type() const {
		return type_;
	}

	void set_type( const std::string &new_type) {
		type_ = new_type;
	}


	/// print checksum data
	void debug( const std::string &tag, const std::string &label, core::Real data1, core::Real data2=0.0, core::Real data3=0.0 ) const;

	void set_disabled( bool value = true ){ disabled_ = value; }
	bool get_disabled() const { return disabled_; }

	core::Size get_checkpoint_recoveries() const { return count_checkpoint_recoveries_; }
private:
	std::string type_;
	std::vector< std::string > checkpoint_ids_;

	bool disabled_;
	bool delete_checkpoints_;
	core::Size count_checkpoint_recoveries_;


	core::Size file_buffer_size();

	std::vector < FileBuffer > file_buffer;


};

}
}

#endif
