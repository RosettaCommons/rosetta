// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file stuff to compare designs against the native pdb
/// @brief
/// @author Florian Richter, floric@u.washington.edu

// Unit headers
#include <protocols/enzdes/DesignVsNativeComparison.hh>
#include <protocols/enzdes/enzdes_util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
//#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/MetricValue.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <core/scoring/constraints/ConstraintSet.hh> //temporary include

// Utility Headers

#include <basic/Tracer.hh>


// option key includes

#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.enzdes.DesignVsNativeComparison" );

namespace protocols {
namespace enzdes {


DesignVsNativeComparison::DesignVsNativeComparison(){
	native_poses_.clear();
}

DesignVsNativeComparison::~DesignVsNativeComparison(){}

void
DesignVsNativeComparison::compare_to_native(
	core::pose::Pose const & pose,
	utility::vector1< std::pair< std::string, std::string > > const & calculators,
	core::scoring::ScoreFunctionCOP scorefxn,
	utility::vector1< core::io::silent::SilentEnergy > & silent_Es)
{

	//first, extract the pdb code from the pose input tag
	std::string pdb_code = enzutil::get_pdb_code_from_pose_tag( pose );

	std::map<  std::string, core::pose::PoseOP >::iterator map_it = native_poses_.find( pdb_code );

	if( map_it == native_poses_.end() ){

		std::string native_path = basic::options::option[basic::options::OptionKeys::enzdes::compare_native].value();

		std::string native_filename = native_path + pdb_code + ".pdb";
		core::pose::PoseOP new_native( new core::pose::Pose() );

		core::import_pose::pose_from_pdb( *new_native, native_filename );

		//need to score the native, in case the calculators rely on cached energies
		//also need to remove non protein residues, just to make sure
		toolbox::pose_manipulation::remove_non_protein_residues( *new_native );

		(*scorefxn)(*new_native);

		native_poses_.insert( std::pair< std::string, core::pose::PoseOP >( pdb_code, new_native ) );

		map_it = native_poses_.find( pdb_code );

	}

	core::pose::PoseOP native = map_it->second;

	core::Real native_totalE = native->energies().total_energies().dot( native->energies().weights() );

	//remove the ligand from the pose, so the comparison makes sense
	core::pose::PoseOP pureprotpose( new core::pose::Pose( pose ) );

	toolbox::pose_manipulation::remove_non_protein_residues( *pureprotpose );

	//arrghh, bug, when scoring a pose where one of the deleted residues was constrained, program crashes
	//because the deleted residue survives in the long range energy containers:(
	//temporary hack workaround: put an empty constraint set into the pureprotpose
	pureprotpose->constraint_set( core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() ) );
	core::Real pose_totalE = (*scorefxn)(*pureprotpose );

	//ok, we have our native, now use each of the calculators to compare
	//and of course we're also looking at the total energy
	silent_Es.push_back( core::io::silent::SilentEnergy( "NaCo_dTotE", pose_totalE - native_totalE, 1, 12) );

	basic::MetricValue< core::Real > mval_real;
	basic::MetricValue< core::Size > mval_size;

	core::Real poseval, nativeval;

	for( utility::vector1< std::pair< std::string, std::string > >::const_iterator calc_it = calculators.begin();
			 calc_it != calculators.end(); ++calc_it ){


		if( calc_it->first == "hbond_pm" || calc_it->first == "burunsat_pm" || calc_it->first == "NLconts_pm" || calc_it->second == "total_pos_charges" || calc_it->second == "total_neg_charges" ){
			pureprotpose->metric( calc_it->first, calc_it->second, mval_size );
			poseval = mval_size.value();

			native->metric( calc_it->first, calc_it->second, mval_size );
			nativeval = mval_size.value();
		}
		else{
			pureprotpose->metric( calc_it->first, calc_it->second, mval_real );
			poseval = mval_real.value();

			native->metric( calc_it->first, calc_it->second, mval_real );
			nativeval = mval_real.value();
		}

		//std::cerr << "For " << calc_it->first << ", design val is " << poseval << " and nativeval is " << nativeval << std::endl;

		//now create a silent energy out of the calculated value
		core::Real native_diff = poseval - nativeval;
		std::string se_name = "NaCo_"+calc_it->first;
		if( calc_it->first == "charges_pm" ) se_name = "NaCo_" + calc_it->second;
		int width = std::max( 10, (int) se_name.length() + 3 );

		//core::io::SilentEnergy se( se_name , native_diff, 1, width );
		silent_Es.push_back( core::io::silent::SilentEnergy ( se_name , native_diff, 1, width ) );

	} //loop over calculators


} //compare_to_native


}//enzdes
}//protocols
