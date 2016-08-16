// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/make_rot_lib/MakeRotLibJobInputter.cc
/// @brief Implementation file for MakeRotLibJobInputter class
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
#include <utility/file/FileName.hh>

// basic headers
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/make_rot_lib.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/util.hh>

// c++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace make_rot_lib {

static THREAD_LOCAL basic::Tracer TR( "protocols.make_rot_lib.MakeRotLibJobInputter" );

protocols::make_rot_lib::MakeRotLibJobInputter::MakeRotLibJobInputter() :
	jd2::JobInputter()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Instantiate MakeRotLibJobInputter" << std::endl;

	// add ACE and NME residues even if they're not uncommented!
	// Don't do this! They are uncommented!!
	//option[ in::file::extra_res_fa ].push_back( utility::file::FileName( std::string( "terminal/ACE.params" ) ) );
	//option[ in::file::extra_res_fa ].push_back( utility::file::FileName( std::string( "terminal/NME.params" ) ) );

	mrlod_ = MakeRotLibOptionsDataOP( new MakeRotLibOptionsData( option[ OptionKeys::make_rot_lib::options_file ].value() ) );
}

protocols::make_rot_lib::MakeRotLibJobInputter::~MakeRotLibJobInputter(){}

/// @details This function will first see if the pose already exists in the Job. If not, it will read it into the pose reference, and hand a COP cloned from that pose to the Job. If the pose pre-exists it just copies the COP's pose into it. The Job object (within its InnerJob) contains a PoseCOP. This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.
void
protocols::make_rot_lib::MakeRotLibJobInputter::pose_from_job( core::pose::Pose & pose, jd2::JobOP job)
{
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys::make_rot_lib;

	TR << "MakeRotLibJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
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
		if ( option[ use_terminal_residues ].value() ) {
			std::string fullname = mrlod_->get_name();
			//TR << "fullname is " << fullname << std::endl;
			// make single residue pose
			if ( option[ use_terminal_residues ].value( ) ) {
				ResidueTypeSetCOP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
				ResidueType const & ACE( RTS->name_map( "ACE" ) );
				ResidueType const & RT( RTS->name_map( fullname ) );
				ResidueType const & NME( RTS->name_map( "NME" ) );
				Residue N( ACE, true );
				Residue R( RT, true );
				Residue C( NME, true );
				pose.append_residue_by_jump( N, 1 );
				pose.append_residue_by_bond( R, true );
				pose.append_residue_by_bond( C, true );
			}
		} else {
			std::stringstream fullname;
			fullname << mrlod_->get_name() << patch_name;

			// make single residue pose
			ResidueTypeSetCOP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
			ResidueType const & RT( RTS->name_map( fullname.str() ) );
			Residue R( RT, true );
			pose.append_residue_by_jump( R, 1 );
		}

		// load pose in to inner job
		load_pose_into_job(pose, job);

	} else {
		TR << "filling pose from saved copy " << job->input_tag() << std::endl;
		pose = *(job->inner_job()->get_pose());
	}

}

/// @details this function determines what jobs exist from the make rot lib options data. Each job calculates the rotamers for one omg, phi, psi, eps bin.
void protocols::make_rot_lib::MakeRotLibJobInputter::fill_jobs( jd2::JobsContainer & jobs )
{
	TR << "MakeRotLibJobInputter::fill_jobs" << std::endl;

	jobs.clear();

	// nstruct will always be 1
	core::Size nstruct( 1 );
	bool semirotameric = mrlod_->get_semirotameric();
	core::Size n_bb = mrlod_->get_n_bb();
	utility::vector1 < core::Size > bb_ids = mrlod_->get_bb_ids();
	// make single InnerJob that will be used in each job
	jd2::InnerJobOP ij( new jd2::InnerJob( mrlod_->get_name(), nstruct ) );

	// sanity check
	runtime_assert( mrlod_->get_omg_range().step != 0 ); // results in infinite loop
	for ( core::Size i = 1; i <= n_bb; ++i ) {
		runtime_assert( mrlod_->get_bb_range( i ).step != 0 ); // results in infinite loop
	}
	runtime_assert( mrlod_->get_eps_range().step != 0 ); // results in infinite loop
	//TR << "Sanity checks on step size nonzero passed." << std::endl;

	runtime_assert( mrlod_->get_omg_range().low <= mrlod_->get_omg_range().high ); // results in no loop
	for ( core::Size i = 1; i <= n_bb; ++i ) {
		runtime_assert( mrlod_->get_bb_range( i ).low <= mrlod_->get_bb_range( i ).high ); // results in no loop
	}
	runtime_assert( mrlod_->get_eps_range().low <= mrlod_->get_eps_range().high ); // results in no loop
	//TR << "Sanity checks on low less than or equal to high passed." << std::endl;

	// initialize bb loop indices
	utility::vector1 < core::Real > i;
	utility::vector1 < core::Real > maxes;
	utility::vector1 < core::Real > steps;
	i.resize(n_bb+1);
	maxes.resize(n_bb+1);
	steps.resize(n_bb+1);
	// Each array contains one more than n_bb because we need an extra counter
	// for the arbitrary depth for loop.
	for ( core::Size tmp = 1; tmp <= n_bb + 1; ++tmp ) {
		if ( tmp == n_bb + 1 ) {
			i[ tmp ] = 0;
			maxes[ tmp ] = 1;
			steps[ tmp ] = 1;
		} else {
			i[ tmp ] = mrlod_->get_bb_range( tmp ).low;
			maxes[ tmp ] = mrlod_->get_bb_range( tmp ).high;
			steps[ tmp ] = mrlod_->get_bb_range( tmp ).step;
		}
	}
	//TR << "initialized loop indices " << std::endl;

	// create a job for all combinitorial combinations of omg, phi, psi and eps
	for ( core::Real o(mrlod_->get_omg_range().low ); o <= mrlod_->get_omg_range().high; o += mrlod_->get_omg_range().step ) {
		//for ( core::Real h(mrlod_->get_phi_range().low ); h <= mrlod_->get_phi_range().high; h += mrlod_->get_phi_range().step ) {
		// for ( core::Real s(mrlod_->get_psi_range().low ); s <= mrlod_->get_psi_range().high; s += mrlod_->get_psi_range().step ) {
		core::Size p = 1;
		while ( i [ n_bb + 1 ] == 0 ) {
			utility::vector1< core::Real > bbs;
			bbs.resize( n_bb );
			for ( core::Size init_i = 1; init_i <= n_bb; ++init_i ) {
				bbs[ init_i ] = i[ init_i ];
			}
			for ( core::Real e( mrlod_->get_eps_range().low ); e <= mrlod_->get_eps_range().high; e += mrlod_->get_eps_range().step ) {
				TR << "pushing the omg: " << o;
				for ( core::Size push_i = 1; push_i <= n_bb; ++push_i ) {
					TR << " bb" << bb_ids[ push_i ] << ": " << bbs[ push_i ];
				}
				TR << " eps: " << e << " bin" << std::endl;
				jobs.push_back( jd2::JobOP( new MakeRotLibJob( ij, nstruct, o, bbs, bb_ids, e, mrlod_, semirotameric ) ) );
			}

			i[ 1 ] += steps[ 1 ];
			while ( i[ p ] > maxes[ p ] ) {
				if ( p <= n_bb ) i[ p ] = mrlod_->get_bb_range( p ).low;
				else             i[ p ] = 0;
				p = p + 1;
				i[ p ] += steps[ p ];
				if ( i[ p ] <= maxes[ p ] ) p = 1;
			}
		}
	}
	TR << "pushed all jobs" << std::endl;

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
	return protocols::jd2::JobInputterOP( new MakeRotLibJobInputter );
}

}//make_rot_lib
}//protocols
