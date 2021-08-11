// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.hh
/// @brief a region of a Pose with secondary structures of known DSSP at either teminus
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_pose_sewing_data_storage_PoseWithTerminalSegmentsOfKnownDSSP_hh
#define INCLUDED_protocols_pose_sewing_data_storage_PoseWithTerminalSegmentsOfKnownDSSP_hh

#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <string>

namespace protocols {
namespace pose_sewing {
namespace data_storage {

/// @brief a region of a Pose with secondary structures of known DSSP at either teminus
class PoseWithTerminalSegmentsOfKnownDSSP : public utility::VirtualBase {

public:

	PoseWithTerminalSegmentsOfKnownDSSP();
	PoseWithTerminalSegmentsOfKnownDSSP(PoseWithTerminalSegmentsOfKnownDSSP const &);

	virtual ~PoseWithTerminalSegmentsOfKnownDSSP();

	PoseWithTerminalSegmentsOfKnownDSSPOP
	clone() const;

	void
	set_N_term_DSSP (char);

	char
	get_N_term_DSSP() const;

	void
	set_C_term_DSSP (char);

	char
	get_C_term_DSSP() const;

	void
	set_N_term_length ( core::Size );

	core::Size
	get_N_term_length () const;

	void
	set_C_term_length ( core::Size );

	core::Size
	get_C_term_length () const;

	void
	set_filename( std::string const & filename);

	std::string
	get_filename() const;

	void
	set_segfile_path( std::string const & segfile_path);

	std::string
	get_segfile_path() const;

	void
	set_secstruct( std::string const & secstruct);

	std::string
	get_secstruct() const;

	void
	set_source_pose ( core::pose::PoseCOP );

	void
	store_source_pose_for_segment( std::string const & src_pose, bool store_mmTF = true, bool output_elements = true );

	core::pose::PoseOP
	create_source_pose_for_segment(bool add_elements);

	void
	store_reference_pdb( std::string );

	///@brief Get the source pose by returning an already-loaded pose
	/// or loading a new pose.  By default we clone this new pose and store it.
	///  Set clone_if_new to false if this is not desired (for speed) - such as is done
	///  during SewAnythingAddMover.
	core::pose::PoseCOP
	get_source_pose(bool clone_if_new=true);

	///@brief Get the source pose by returning an already-loaded pose
	/// or loading a new pose.  By default we clone this new pose and store it.
	///  Set clone_if_new to false if this is not desired (for speed) - such as is done
	///  during SewAnythingAddMover.
	///
	/// This gives you an OP - the one just loaded or a clone of the stored one.
	core::pose::PoseOP
	get_source_pose_op(bool clone_if_new=true);

private:

	char N_term_DSSP_;
	char C_term_DSSP_;
	core::Size N_term_length_;
	core::Size C_term_length_;
	std::string filename_;
	std::string secstruct_;
	std::string segfile_path_ = "unknown";
	core::pose::PoseCOP source_pose_;

};


} //protocols
} //pose_sewing
} //data_storage



#endif //INCLUDED_protocols_pose_sewing_data_storage_PoseWithTerminalSegmentsOfKnownDSSP_hh





