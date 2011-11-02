// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/read_patchdock.hh
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_read_patchdock_hh
#define INCLUDED_protocols_protein_interface_design_read_patchdock_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#ifdef WIN32
#include <string>
#endif

// Project Headers
// C++ headers
namespace protocols {
namespace protein_interface_design {

struct Transformation
{
  typedef core::Real Real;

	Real alpha, beta, gamma; // Euler angles
	numeric::xyzVector< Real > translation; // translation
};

class PatchdockReader : public utility::pointer::ReferenceCount
{
public:
	PatchdockReader();
	~PatchdockReader();
	/// @brief if no native is read
	void read_poses( core::pose::Pose & input_pose, std::string & input_tag );

/// @brief reads input and native poses from file. If patchdock flags are used will read the patchdock transformation
/// and transform the input pose accordingly
	void read_poses( core::pose::Pose & input_pose, core::pose::Pose & native_pose, std::string & input_tag, std::string & native_tag );

	void read_patchdock( std::string & input_tag, std::string & native_tag );
	Transformation read_patchdock_entry();

	void clear_internals();
	core::Size number_of_patchdock_entries();
private:
	std::string patchdock_fname_;
	core::Size patchdock_entry_num_;
	core::pose::PoseOP saved_input_pose_, saved_native_pose_;
	std::string saved_input_tag_, saved_native_tag_;
	utility::vector1< Transformation > saved_transformations_;
};
}
}

#endif /*INCLUDED_READ_PATCHDOCK_H_*/

