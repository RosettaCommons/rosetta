// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Nikolas Sgourakis


//#include <basic/options/option.hh>
//#include <core/scoring/ResidualDipolarCoupling.hh>
#include <protocols/scoring/ResidualDipolarCouplingRigidSegments.hh>
#include <iostream>
#include <devel/init.hh>
#include <protocols/jd2/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/moves/Mover.hh>


using namespace protocols::moves;
using namespace protocols::scoring;

class RDCScoreMover : public Mover {
public:
	RDCScoreMover( ResidualDipolarCouplingRigidSegmentsOP );

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "RDCScoreMover"; }


private:

	ResidualDipolarCouplingRigidSegmentsOP rdcrs_;

};


RDCScoreMover::RDCScoreMover(ResidualDipolarCouplingRigidSegmentsOP rdcrs ):
	rdcrs_(rdcrs)
{}

void RDCScoreMover::apply( core::pose::Pose& pose ) {


	//	core::scoring::constraints::add_constraints_from_cmdline( pose, *sfxn_ );
  // print jobname and score to screen
	std::cout  << protocols::jd2::current_output_name() << ' ' << rdcrs_->compute_total_score(pose) << ' ' <<rdcrs_->compute_pairwise_score()
   <<std::endl;
}


int main( int argc, char * argv []   ){

	try {

	using namespace protocols::scoring;
	using namespace protocols::jd2;

	using namespace core;

  devel::init( argc, argv );


	  ResidualDipolarCouplingRigidSegmentsOP my_rdcs = new ResidualDipolarCouplingRigidSegments ;
		std::cout << *my_rdcs;
	//  cout <<"hello";


	try{
		RDCScoreMover* scoremover = new RDCScoreMover( my_rdcs );

			/*		if ( option[ in::file::keep_input_scores ]() ){
			scoremover->set_keep_input_scores();
		}
		if ( option[ rescore::skip]() ){
			scoremover->set_skip_scoring();
			} */


		protocols::jd2::JobDistributor::get_instance()->go( scoremover );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
	}

 	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
