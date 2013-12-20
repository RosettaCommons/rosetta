// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/DomainAssemblyJobInputter.hh
/// @brief  Declaration of DomainAssemblyJobInputter which creates a fusion pose from two PDB files 
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

#ifndef INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_hh
#define INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <devel/domain_assembly/DomainAssemblyJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

//utility headers

namespace devel {
namespace domain_assembly {

///@details This is the simplest implementation of JobInputter, which reads from -s/-l and PDB files.
class DomainAssemblyJobInputter : public protocols::jd2::JobInputter
{
public:

	DomainAssemblyJobInputter();

	virtual ~DomainAssemblyJobInputter();

	///@brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference
 	virtual void pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job );

	///@brief this function determines what jobs exist using the da options file name and nstruct
	virtual void fill_jobs( protocols::jd2::Jobs & jobs );

	/// @brief Return the type of input source that the DomainAssemblyJobInputter is currently
	///  using.
	/// @return Always <em>PDB_FILE</em>.
	virtual protocols::jd2::JobInputterInputSource::Enum input_source() const;

}; // DomainAssemblyJobInputter

} // namespace domain_assembly
} // namespace devel

#endif //INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_HH
