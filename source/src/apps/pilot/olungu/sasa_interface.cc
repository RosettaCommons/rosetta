// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /apps/pilot/olungu/sasa_interface.cc
/// @brief calculate metrics over antibody peritopes
/// @author Oana Lungu


#include <basic/MetricValue.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <devel/init.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/metrics.hh>

#include <string>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("");

//namespaces
using namespace core;
using namespace core::pose;
using namespace utility;
using namespace protocols::antibody;

int
main( int argc, char* argv[] ){

	try{

		devel::init(argc, argv);

		//-s read in PDB
		core::pose::Pose pose;
		std::string pdbname(basic::options::option[ basic::options::OptionKeys::in::file::s ].value()[1]);
		core::import_pose::pose_from_pdb( pose, pdbname );
		protocols::antibody::AntibodyInfo abinfo( pose );

		std::pair<ParatopeMetric< core::Real >, ParatopeMetric<core::Real> > sasa_result = paratope_sasa(pose, abinfo);
		//std::cout << loop << " CDR_sasa " << loop_sasa << std::endl;
		//std::cout << loop <<" CDR hydrophobic sasa " << hydrop_loop_sasa << std::endl;
		//std::cout << "peritope_sasa " << sasa_result.first.paratope << std::endl;
		//std::cout << "peritope hydrophobic sasa " << hydrophobic_sasa << std::endl;
		//std::cout << "Total SASA is: " << Real mr.value() << std::endl;
		//std::cout << "Total hydrophobic SASA is: " << hSASA << std::endl;

		core::SSize total_charge = pose_charge( pose );
		std::cout << "Total charge " << total_charge << std::endl;
		ParatopeMetric< core::SSize>p_charge = paratope_charge( pose, abinfo );
		std::cout << "peritope charge " << p_charge.paratope << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

