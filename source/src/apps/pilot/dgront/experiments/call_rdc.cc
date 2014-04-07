// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief

// libRosetta headers
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <protocols/moves/Mover.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

// RDC and scoring
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/excn/Exceptions.hh>

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(in::file::silent);
  OPT(in::file::rdc);
  OPT(in::file::s);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
}

int
main( int argc, char * argv [] ) {
    try {
        using namespace basic::options;
        using namespace basic::options::OptionKeys;

	devel::init(argc, argv);
	register_options();

//------------- Read the native pose  ----------
    core::pose::Pose native_pose;
    if ( option[ in::file::native ].user() )
        core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]().name() );


    core::scoring::ResidualDipolarCoupling rdc;
    rdc.read_RDC_file();
    std::cout << "native:\t" << rdc.compute_dipscore(native_pose) << 	std::endl;

//------------- Read the pose for scoring ----------
    core::pose::Pose fa_pose;
    utility::vector1<utility::file::FileName> s = option[in::file::s]();
    for(core::Size i=1;i<=s.size();i++) {
	core::import_pose::pose_from_pdb( fa_pose, s[i].name());
	std::cout << s[i].name() << ":\t" << rdc.compute_dipscore(fa_pose) << std::endl;
    }

//------------- Now, let's find a best-scoring fragment based on RDC ----------
    core::Size Nmer_size = 35;
    core::Size frag_from = 10;
    core::scoring::ResidualDipolarCoupling::RDC_lines rdc_data_on_fragment;
    core::scoring::ResidualDipolarCoupling::RDC_lines data = rdc.get_RDC_data();
    for ( core::scoring::ResidualDipolarCoupling::RDC_lines::const_iterator it = data.begin(); it != data.end(); ++it ) {
	if ( it->res1() >= frag_from && it->res1() <= frag_from + Nmer_size - 1 ) {

	    core::scoring::RDC clone(it->res1() - frag_from + 1, it->atom1(), it->res2() - frag_from + 1,
			it->atom2(), it->Jdipolar());
	    rdc_data_on_fragment.push_back( clone );

//	    std::cerr << "Copied: "<< *it << " as " << clone << std::endl;
	}
    }
    core::scoring::ResidualDipolarCoupling rdc_on_Nmer( rdc_data_on_fragment );

    core::pose::Pose frag_pose;
    core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
    std::string seq = native_pose.sequence().substr(frag_from-1,Nmer_size);
    core::pose::make_pose_from_sequence(frag_pose, seq,*rsd_set);
    core::scoring::store_RDC_in_pose( &rdc_on_Nmer, frag_pose );

    for(core::Size ires=1;ires< native_pose.total_residue() - Nmer_size - 1;ires++) {
	for(core::Size j=0;j<Nmer_size;j++) {
    	    frag_pose.set_phi( j+1, native_pose.phi(ires+j) );
    	    frag_pose.set_psi( j+1, native_pose.psi(ires+j) );
    	    frag_pose.set_omega( j+1, native_pose.omega(ires+j) );
	}
	std::cout << ires << " : " << rdc_on_Nmer.compute_dipscore(frag_pose) << std::endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}

