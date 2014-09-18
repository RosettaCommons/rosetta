// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file protocols/pmut_scan/AlterSpecDisruptionDriver.cc
/// @brief A protocol that tries to find interface-disrupting mutations as phase 1 of an altered-specificity protocol
/// @author Steven Lewis smlewi@gmail.com

// Unit headers
#include <protocols/pmut_scan/AlterSpecDisruptionDriver.hh>
#include <protocols/pmut_scan/Mutant.hh>

//project Headers
#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace pmut_scan {

static thread_local basic::Tracer TR( "protocols.pmut_scan.AlterSpecDisruptionDriver" );

///
/// @begin AlterSpecDisruptionDriver::AlterSpecDisruptionDriver
///
/// @brief
/// Main constructor for the class.
///
AlterSpecDisruptionDriver::AlterSpecDisruptionDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string list_file, bool output_mutant_structures ) :
	PointMutScanDriver(pdb_file_names, double_mutant_scan, list_file, output_mutant_structures), IAM_(NULL)
{
	IAM_ = new protocols::analysis::InterfaceAnalyzerMover(
		1, //interface_jump
		true, //output to Tracer - we are ignoring this output and will collect from getters
		get_scorefxn(), //pass in the scorefunction we are using
		false, //do not compute packstat
		true, //do repack input - we already did, but this will keep the PackerTasks similar
		true,  //do repack separated pose
		false); //do not bother with JD2 tracer name hookups

	//IAM_->set_use_resfile(false);
	IAM_->set_use_centroid_dG(false);
	IAM_->set_compute_packstat(false); //should be redundant
	IAM_->set_compute_interface_sc(false);
	IAM_->set_compute_separated_sasa(false);
	IAM_->set_compute_interface_energy(true); //this is the one value we'll need
	IAM_->set_calc_hbond_sasaE(false);
	IAM_->set_compute_interface_delta_hbond_unsat(false);
	IAM_->set_skip_reporting(true);

	interface_set_.clear();

}


///
/// @begin AlterSpecDisruptionDriver::~AlterSpecDisruptionDriver
///
/// @brief
/// Destructor.
///
AlterSpecDisruptionDriver::~AlterSpecDisruptionDriver() {}

///@details calculate dG of binding by scoring, separating, repacking, and rescoring

core::Energy AlterSpecDisruptionDriver::score(core::pose::Pose & pose) {

	//It would be good if there were a way to test an invariant here - primarily that we have an interface of two chains!
	//runtime_assert(pose.num_chains() >= 2);
	if(!(pose.conformation().num_chains() >= 2)){
		throw utility::excn::EXCN_Msg_Exception("AlterSpecDisruptionDriver requires two chains to calculate an interface to disrupt.");
	}

	//no clear way to get IAM to use the same TaskFactory
	IAM_->set_compute_interface_energy(true); //speed
	IAM_->apply(pose);
	//NOTICE THIS IS NEGATED!  The parent code is expecting to stabilize things; we want to destabilize them.
	return -IAM_->get_separated_interface_energy();

}

///@brief offers a chance for child classes to inject mutant selection logic
bool AlterSpecDisruptionDriver::reject_mutant( Mutant const & mutant, core::pose::Pose const & pose ){

	if(reject_on_chains(mutant)) return true;

	if(reject_on_interface(mutant, pose)) return true;

	return false;

}

///@brief reject Mutant based on chain IDs of constituent mutations
bool AlterSpecDisruptionDriver::reject_on_chains( Mutant const & mutant ){

	if(mutant.n_mutations() >= 2){ //if there is more than one mutation, check that they are on the same chain.  Disruptions on different chains are likely unrecoverable clashes.
		typedef utility::vector1< MutationData >::const_iterator mut_iter;

		char const chain(mutant.mutations_begin()->pdb_chain());
		for ( mut_iter it(mutant.mutations_begin()), end(mutant.mutations_end()); it != end; ++it ) {
			if(chain != it->pdb_chain()) {
				//TR << "skipping Mutant " << mutant << " because it has a mutation on chain " << it->pdb_chain() << " which does not match the first mutation's chain " << chain << std::endl;
				return true;
			} //if skipping
		} //for all mutations
	} //if we have more than one mutation

	return false;
}

///@brief reject Mutant based on interface-ness of constituent mutations.  This function accumulates state - on first call it determines what the interface is, and assumes that interface is the same for all subsequent calls!  If you are trying to use different Poses this may cause problems.
bool AlterSpecDisruptionDriver::reject_on_interface( Mutant const & mutant, core::pose::Pose const & pose )
{

	if( interface_set_.empty() ){
		IAM_->set_compute_interface_energy(false); //speed

		core::pose::Pose copy(pose);
		IAM_->apply(copy);
		interface_set_ = IAM_->get_interface_set();

		// TR << "interface set:";
		// for(std::set< core::Size >::const_iterator it(interface_set_.begin()), end(interface_set_.end()); it != end; ++it)
		// 	TR << " " << *it;
		// TR << std::endl;

		IAM_->set_compute_interface_energy(true); //speed
	}

	typedef utility::vector1< MutationData >::const_iterator mut_iter;

	std::set< core::Size >::iterator const interface_end(interface_set_.end());
	for( mut_iter it(mutant.mutations_begin()), end(mutant.mutations_end()); it != end; ++it ) {
		core::Size const resid(it->pose_resnum());
		if(interface_set_.find(resid) == interface_end) { //if a mutation is not in the interface set, skip this Mutant
			//TR << "skipping Mutant " << mutant << " because it contains position " << resid << " which is not in the interface" << std::endl;
			return true;
		}
	}

	return false;
}



} // namespace pmut_scan
} // namespace protocols

