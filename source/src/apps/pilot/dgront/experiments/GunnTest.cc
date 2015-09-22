// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 *
 *  Created on: Apr 21, 2009
 *      Author: dgront
 */

#include <core/types.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/abinitio/GunnCost.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


static THREAD_LOCAL basic::Tracer tr( "GunnTest" );

void register_options() {

  OPT( in::file::native );
  OPT( in::file::s );
}


class GunnTest : public protocols::moves::Mover {
public:

    GunnTest() {

        using namespace core;
        using namespace basic::options;
	using namespace basic::options::OptionKeys;

	frag_size_ = 9;
	native = new core::pose::Pose;
	core::import_pose::pose_from_pdb(*native, option[ in::file::native ]());

	for(Size i=1;i<=native->total_residue()-frag_size_+1;i++) {
	    data_p.push_back( new protocols::abinitio::GunnTuple );
	    data_n.push_back( new protocols::abinitio::GunnTuple );
	}
	computeGunnTuples(*native,frag_size_,data_n);
    }

    virtual ~GunnTest() {}

    virtual void apply( core::pose::Pose & pose ) {

	computeGunnTuples(pose,frag_size_,data_p);
	for(Size i=1;i<=data_p.size();i++) {
	    protocols::abinitio::GunnTuple & t = *data_p[i];
	    std::cout<<t.q1<<" ";
	    std::cout<<t.q2<<" ";
	    std::cout<<t.q3<<" ";
	    std::cout<<t.q4<<" ";
	    std::cout<<t.q5<<" ";
	    std::cout<<t.q6<<" ";
	    std::cout<<gun.score_tuple(*data_n[i],t)<<std::endl;
	}
    }

    void computeGunnTuples(core::pose::Pose & pose,Size frag_size,utility::vector1<protocols::abinitio::GunnTupleOP> result) {

	for(Size i=1;i<=pose.total_residue()-frag_size + 1;i++) {
	    protocols::abinitio::GunnTuple & t = *result[i];
	    gun.compute_gunn( pose, i, i+frag_size_-1, t);
	}
    }

    Size frag_size_;
    utility::vector1<protocols::abinitio::GunnTupleOP> data_n;
    utility::vector1<protocols::abinitio::GunnTupleOP> data_p;
    core::pose::PoseOP native;
    protocols::abinitio::GunnCost gun;
};


int main(int argc, char * argv[]) {
try {
    using namespace protocols;
    using namespace protocols::jobdist;
    using namespace protocols::moves;

  register_options();
  devel::init(argc, argv);

    GunnTest test;
    not_universal_main( test );
} catch ( utility::excn::EXCN_Base const & e ) {
                          std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                              }

  return 0;
}
