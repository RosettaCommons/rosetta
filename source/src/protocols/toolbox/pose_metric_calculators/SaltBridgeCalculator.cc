// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// How many salt bridge interactions are there?
///
/// @details
/// Not much detailed here. Iterate through the oxygens of acidic residues and compare that to
/// the distance of the polar hydrogens in basic residues. Default distance is 3.2A.
/// Wait, you want to know how to use this? Well, within your protocol, you need to do the following:
/// First, create the calculator. To do this, see below:
/// core::pose::metrics::PoseMetricCalculatorOP sb_calculator = new protocols::toolbox::pose_metric_calculators::SaltBridgeCalculator();
/// Then you must register this so that the pose understands it. See below:
/// core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sb_metric", sb_calculator );
/// To actually get the metric, you have to print it. For example:
/// core::pose::Pose pose;
/// pose.print_metric("sb_metric", "salt_bridge")
/// Where sb_metric is the name that it is registered under and "salt_bridge" is the key, seen below.
///
///
///
/// @author
/// Steven Combs
///
/////////////////////////////////////////////////////////////////////////
#include <protocols/toolbox/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.SaltBridgeCalculator" );

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {


/// @brief default constructor sets distance_cutoff to 3.2. This is what is usually defined as a Hbond between heavy atom and Hydrogen
SaltBridgeCalculator::SaltBridgeCalculator() :
	distance_cutoff_(3.2),
	salt_bridge_total_(0)
{

}


/// @brief constructur where you define what the distance cutoff is for the salt bridge
SaltBridgeCalculator::SaltBridgeCalculator(core::Real dist_cutoff) :
	distance_cutoff_(dist_cutoff),
	salt_bridge_total_(0)
{

}


void SaltBridgeCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const{
	if ( key == "salt_bridge" ) {
		basic::check_cast( valptr, &salt_bridge_total_, "salt_bridge expects to return a Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( salt_bridge_total_ );

	} else {
		basic::Error() << "SaltBridgeCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}
}


std::string SaltBridgeCalculator::print( std::string const & key ) const{
	if ( key == "salt_bridge" ) {
		return utility::to_string(salt_bridge_total_);
	}
	basic::Error() << "SaltBridgeCalculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


/// @brief not sure why they name this function recompute as you are actually computing the metric. Whateva
void SaltBridgeCalculator::recompute(core::pose::Pose const & pose){
	salt_bridge_total_ = 0;
	//start iterating through the residues
	for ( core::Size res_num1=1; res_num1 <= pose.n_residue(); ++res_num1 ) {
		//assign the number to a residue based on the seqpos
		core::conformation::Residue acceptor(pose.residue(res_num1));
		//continue only if the residue is either aspartic or glutamic acid
		if ( acceptor.name3() == "ASP" || acceptor.name3() =="GLU" ) {
			for ( core::Size res_num2=1; res_num2 <= pose.n_residue(); ++res_num2 ) {
				core::conformation::Residue donate(pose.residue(res_num2));
				//only continue if this residue is a his, lys or arg
				if ( donate.name3() == "HIS" || donate.name3() == "LYS" || donate.name3() =="ARG" ) {
					//set up a flag that will stop us from double counting salt bridges
					bool get_out_of_loop=false;
					//start iteration through acceptor heavy atoms
					for (
							core::chemical::AtomIndices::const_iterator
							anum  = acceptor.accpt_pos().begin(),
							anume = acceptor.accpt_pos().end(); anum != anume; ++anum ) {
						core::Size const acc_atm( *anum );


						for
							(
									core::chemical::AtomIndices::const_iterator
									don_num = donate.Hpos_polar_sc().begin(),
									don_nume = donate.Hpos_polar_sc().end();
									don_num != don_nume; ++don_num
									) {
							core::Size const don_atm(*don_num);

							//get the distance between the donor residue and polar hydrogen sidechain
							core::Real distance(acceptor.xyz(acc_atm).distance(donate.xyz(don_atm)) );

							if ( get_out_of_loop ==false && distance < distance_cutoff_ ) {
								++salt_bridge_total_;
								get_out_of_loop=true;
							}

						}


					}


				}
			}
		}

	}


}


}
}
}
