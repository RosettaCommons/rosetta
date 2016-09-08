// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
//#include <protocols/FlexPepDocking/FlexPepDockingProtocol.hh>

#include <numeric/constants.hh>
#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/util.hh>//option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/scoring/Interface.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;


class RMSByResStatistics : public protocols::moves::Mover
{
public:
	RMSByResStatistics() {}

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const
	{ return std::string("RMSByResStatistics"); }
};


void
RMSByResStatistics::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;
	using core::Real;
	protocols::jd2::JobOP cur_job
		( protocols::jd2::JobDistributor::get_instance()->current_job() );
	core::pose::Pose native = *this->get_native_pose();
	Size nres = pose.size();
	ObjexxFCL::FArray1D_bool res_subset( nres, false );
	ObjexxFCL::FArray1D_bool pep_subset( nres, false );
	ObjexxFCL::FArray1D_bool pep_intrf_subset( nres, false );

	// calculate native interface
	core::Size rb_jump_ = 1;
	protocols::scoring::InterfaceOP native_interface
		= new protocols::scoring::Interface();
	native.update_residue_neighbors();
	native_interface->calculate( native );

	Size i = 1;
	while(pose.chain(i)==1) // skip to chain B
		i++;
  for( ; i <= nres ; ++i)
    {
			res_subset(i) = true;
			res_subset(i-1) = false;
			pep_subset(i) = true;
			if(native_interface->is_interface(i))
				pep_intrf_subset(i)=true;

			Real res_rms =
				rmsd_no_super_subset( pose, native, res_subset, is_polymer_heavyatom );
			std::ostringstream ss;
			ss << i;
			ss << pose.residue(i).name1();
			cur_job->add_string_real_pair( ss.str(), res_rms );

			Real res_rms_bb =
				rmsd_no_super_subset( pose, native, res_subset, is_protein_backbone );
			std::ostringstream ss_bb;
			ss_bb << "bb_";
			ss_bb << i;
			ss_bb << pose.residue(i).name1();
			cur_job->add_string_real_pair( ss_bb.str(), res_rms_bb );

			Real res_rms_ca_cb =
				rmsd_no_super_subset( pose, native, res_subset, is_protein_CA_or_CB );
			std::ostringstream ss_ca_cb;
			ss_ca_cb << "CA_CB_";
			ss_ca_cb << i;
			ss_ca_cb << pose.residue(i).name1();
			cur_job->add_string_real_pair( ss_ca_cb.str(), res_rms_ca_cb );
		}

	Real pep_rms_bb =
		rmsd_no_super_subset( pose, native, pep_subset, is_protein_backbone );
	cur_job->add_string_real_pair( "rmsBB", pep_rms_bb );
	Real pep_intrf_rms_bb =
		rmsd_no_super_subset( pose, native, pep_intrf_subset, is_protein_backbone );
	cur_job->add_string_real_pair( "irmsBB", pep_intrf_rms_bb );


}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
  using namespace core;
  using namespace std;
	using basic::options::option;
	using namespace basic::options::OptionKeys;

  init(argc, argv);
	// initialize calculator
	protocols::moves::MoverOP rmsCalculator = new RMSByResStatistics();
	// read native pose: (TODO: look how this should be handled in Job Distributor 2)
	core::pose::PoseOP native_pose = new core::pose::Pose();
	if (option[ in::file::native ].user()) {
	  if ( option[ in::file::centroid_input ].user() ) {
	    core::import_pose::centroid_pose_from_pdb( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
	  } else {
	    core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
	  }
	}
	rmsCalculator->set_native_pose( native_pose );
	// run
	protocols::jd2::JobDistributor::get_instance()->go( rmsCalculator );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
