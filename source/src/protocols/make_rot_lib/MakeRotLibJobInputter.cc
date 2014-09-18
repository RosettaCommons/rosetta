// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJobInputter.cc
/// @brief  Implementation file for MakeRotLibJobInputter class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/make_rot_lib/MakeRotLibJobInputter.hh>
#include <protocols/make_rot_lib/MakeRotLibJobInputterCreator.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>
#include <protocols/make_rot_lib/MakeRotLibJob.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/JobInputter.hh>

// core headers
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

///Utility headers
#include <utility/vector1.hh>

// basic headers
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/make_rot_lib.OptionKeys.gen.hh>
#include <basic/options/util.hh>

// c++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace make_rot_lib {

static thread_local basic::Tracer TR( "protocols.make_rot_lib.MakeRotLibJobInputter" );

protocols::make_rot_lib::MakeRotLibJobInputter::MakeRotLibJobInputter() :
	jd2::JobInputter()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

  TR << "Instantiate MakeRotLibJobInputter" << std::endl;

	mrlod_ = new MakeRotLibOptionsData( option[ OptionKeys::make_rot_lib::options_file ].value() );
}

protocols::make_rot_lib::MakeRotLibJobInputter::~MakeRotLibJobInputter(){}

///@details This function will first see if the pose already exists in the Job.  If not, it will read it into the pose reference, and hand a COP cloned from that pose to the Job. If the pose pre-exists it just copies the COP's pose into it.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.
void
protocols::make_rot_lib::MakeRotLibJobInputter::pose_from_job( core::pose::Pose & pose, jd2::JobOP job)
{
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::conformation;

  TR << "MakeRotLibJobInputter::pose_from_job" << std::endl;

  if( !job->inner_job()->get_pose() ){
    TR << "filling pose from Job " << job->input_tag() << std::endl;

		// get correct patch based on polymer type
		std::string patch_name;
		if ( mrlod_->get_polymer_type() == PEPTIDE ) {
			TR << "Making a PROTEIN rotamer library..." << std::endl;
			patch_name = "_p:MethylatedCtermProteinFull_p:AcetylatedNtermProteinFull";
		} else if ( mrlod_->get_polymer_type() == PEPTOID ) {
			TR << "Making a PEPTOID rotamer library..." << std::endl;
			patch_name = "_p:AcetylatedNtermDimethylatedCtermPeptoidFull";
		} else {
			utility_exit_with_message("The MakeRotLib protocol only works for peptides and peptoids.");
		}

		// create fully patched name
		std::stringstream fullname;
		fullname << mrlod_->get_name() << patch_name;

		// make single residue pose
		ResidueTypeSetCAP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueType const & RT( RTS->name_map( fullname.str() ) );
		Residue R( RT, true );
		pose.append_residue_by_jump( R, 1 );

		// load pose in to inner job
    load_pose_into_job(pose, job);

  } else {
    TR << "filling pose from saved copy " << job->input_tag() << std::endl;
    pose = *(job->inner_job()->get_pose());
  }

}

///@details this function determines what jobs exist from the make rot lib options data. Each job calculates the rotamers for one omg, phi, psi, eps bin.
void protocols::make_rot_lib::MakeRotLibJobInputter::fill_jobs( jd2::Jobs & jobs )
{
  TR << "MakeRotLibJobInputter::fill_jobs" << std::endl;

  jobs.clear();

	// nstruct will always be 1
	core::Size nstruct( 1 );

	// make single InnerJob that will be used in each job
	jd2::InnerJobOP ij( new jd2::InnerJob( mrlod_->get_name(), nstruct ) );

	// sanity check
	runtime_assert( mrlod_->get_omg_range().step != 0 ); // results in infinite loop
	runtime_assert( mrlod_->get_phi_range().step != 0 ); // results in infinite loop
	runtime_assert( mrlod_->get_psi_range().step != 0 ); // results in infinite loop
	runtime_assert( mrlod_->get_eps_range().step != 0 ); // results in infinite loop

	runtime_assert( mrlod_->get_omg_range().low <= mrlod_->get_omg_range().high ); // results in no loop
	runtime_assert( mrlod_->get_phi_range().low <= mrlod_->get_phi_range().high ); // results in no loop
	runtime_assert( mrlod_->get_psi_range().low <= mrlod_->get_psi_range().high ); // results in no loop
	runtime_assert( mrlod_->get_eps_range().low <= mrlod_->get_eps_range().high ); // results in no loop

	// create a job for all combinitorial combinations of omg, phi, psi and eps
	for ( core::Real o(mrlod_->get_omg_range().low ); o <= mrlod_->get_omg_range().high; o += mrlod_->get_omg_range().step ) {
		for ( core::Real h(mrlod_->get_phi_range().low ); h <= mrlod_->get_phi_range().high; h += mrlod_->get_phi_range().step ) {
			for ( core::Real s(mrlod_->get_psi_range().low ); s <= mrlod_->get_psi_range().high; s += mrlod_->get_psi_range().step ) {
				for ( core::Real e(mrlod_->get_eps_range().low ); e <= mrlod_->get_eps_range().high; e += mrlod_->get_eps_range().step ) {
					TR << "pushing the omg: " << o << " phi: "<< h << " psi: " << s << " eps: " << e << " bin" << std::endl;
					jobs.push_back( jd2::JobOP( new MakeRotLibJob( ij, nstruct, o, h, s, e, mrlod_ ) ) );
				}
			}
		}
	}

}

/// @brief Return the type of input source that the MakeRotLibJobInputter is currently using.
jd2::JobInputterInputSource::Enum MakeRotLibJobInputter::input_source() const
{
  return jd2::JobInputterInputSource::MAKE_ROT_LIB;
}

//CREATOR SECTION
std::string
MakeRotLibJobInputterCreator::keyname() const
{
  return "MakeRotLibJobInputter";
}

protocols::jd2::JobInputterOP
MakeRotLibJobInputterCreator::create_JobInputter() const
{
  return new MakeRotLibJobInputter;
}

}//make_rot_lib
}//protocols
