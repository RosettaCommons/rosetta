// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /apps/pilot/olungu/sasa_interface.cc
/// @brief
/// @author

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh> //necessary for fold tree tricks
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
//#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/antibody/AntibodyUtil.hh>
#include <basic/MetricValue.hh>
#include <devel/init.hh>

// C++ Headers
#include <sstream>
#include <set>
#include <string>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

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

	devel::init(argc, argv);

    //-s read in PDB
    core::pose::Pose pose;
    std::string pdbname(basic::options::option[ basic::options::OptionKeys::in::file::s ].value()[1]);
    core::import_pose::pose_from_pdb( pose, pdbname );

		 // Build atom subsets for computing SASA
    AntibodyInfo ab_info = AntibodyInfo(pose);
    std::cout << ab_info;
    
    
    
		id::AtomID_Map< bool > atom_subset;
		atom_subset.clear();
		core::pose::initialize_atomid_map( atom_subset, pose, false );

// Register a "SasaResiCalculator", giving it the name "sasa_resi"
//		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculator;
//		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa_resi", sasa_calculator );

		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator2 = new core::pose::metrics::simple_calculators::SasaCalculator;
		 core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator2);
//
//  for (core::Size i=1; i<=CDRNameEnum_total; ++i){
//		CDRNameEnum loop = static_cast<CDRNameEnum>(i);
//		Size loop_start = ab_info.get_CDR_loop(loop).start();
//		Size loop_end = ab_info.get_CDR_loop(loop).stop();
////    Size CDR_start = ab_info.get_CDR_start(loop);
////    Size CDR_end = ab_info.get_CDR_end(loop);
//    std::cout << "CDR start " << CDR_start << std::endl;
//    std::cout << "CDR end " << CDR_end << std::endl;
//
//		for (Size ii=loop_start; ii<=loop_end; ++ii){
//				 core::conformation::Residue const & irsd( pose.residue( ii ) );
//				 std::cout << "residue " << irsd.name3() << ii << std::endl;
//
//				 for ( Size ia = 1; ia <= irsd.natoms(); ++ia ) {
//					 id::AtomID const iid( ia, ii );
//					 atom_subset[ iid ] = true;
//				 }
//		}
//	}
// 	// Use this Calculator to report the total sasa for a given Pose
//	// Get the bound SASA, and the SASA for chain1/chain2 alone
//	id::AtomID_Map< Real > atom_sasa ;
//	utility::vector1< Real > residue_sasa ;
//	core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];
//	//Real resi = core::scoring::calc_per_atom_sasa( pose, atom_sasa, residue_sasa, probe_radius, false, atom_subset );
//
//	//std::cout << "SASA of peritope is: " << resi << std::endl;     // report the value//
//
//
//basic::MetricValue<Real> mr;
//pose.metric("sasa","total_sasa",mr);
//
//std::cout << "Total SASA is: " << mr.value() << std::endl;     // report the value
//
//    basic::MetricValue<utility::vector1< Real > > sasaresi;
//    pose.metric("sasa","residue_sasa",sasaresi);
//    
//    std::cout << "sasaresi : " << sasaresi.value() << std::endl;
//    pose.print_metric("sasa","residue_sasa");
//    
}
