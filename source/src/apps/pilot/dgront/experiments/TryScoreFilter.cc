// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Dominik Gront

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_filters/PDDFScoreFilter.hh>

#include <devel/init.hh>

#include <core/types.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/fragment/FragmentIO.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <core/fragment/FragSet.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/assembly.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <string>
#include <sstream>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>



basic::Tracer TR("TryScoreFilter");

OPT_KEY( String, score_name )
OPT_KEY( Real, cutoff )

using namespace core;

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(in::file::native);
  OPT(in::file::s);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
  OPT(score::saxs::ref_pddf);
  OPT(score::saxs::custom_ff);

  NEW_OPT( score_name, "name of the filtered score type", "" );
  NEW_OPT( cutoff, "cutoff value of the score",50 );
}


class TryScoreFilter {

public:

    /// @brief Constructor parses command line and sets up the job
    TryScoreFilter() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace std;

	the_filter_ = new protocols::simple_filters::PDDFScoreFilter();

//------------- set up sampling  -----------------
	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_bb(true);

	std::string frag_large_file, frag_small_file;
	if (option[ in::file::fragA ].user())
	    frag_large_file  = option[ in::file::fragA ]();
	else
	    frag_large_file  = option[ in::file::frag9 ]();
	if (option[ in::file::fragB ].user())
	    frag_small_file  = option[ in::file::fragB ]();
	else
	    frag_small_file  = option[ in::file::frag3 ]();
 	fragset_large_ = core::fragment::FragmentIO(option[ OptionKeys::abinitio::number_9mer_frags ] ).read_data( frag_large_file );
	fragset_small_ = core::fragment::FragmentIO(option[ OptionKeys::abinitio::number_3mer_frags ] ).read_data( frag_small_file );

//	core::scoring::ScoreFunctionOP scorefxn_ = core::scoring::getScoreFunction();
//	scorefxn_->set_weight( core::scoring::pddf_score, 1.0 );

	setup_pose();

	abinitio_ = new protocols::abinitio::ClassicAbinitio( fragset_large_, fragset_small_ ,movemap_ );
	abinitio_->init( cen_pose_ );
	abinitio_->return_centroid( true );
    }

    /// @brief
    void run() {

	for(Size i=1;i<=10;i++) {
	    bool done=false;
	    do {
		abinitio_->apply( cen_pose_ );
	        done = the_filter_->apply( cen_pose_ );
	        TR << the_filter_->recent_score()<<" : "<<done<<std::endl;
    	    } while( !done);

	    std::stringstream ss;
	    ss<<"pddf_"<<i<<".pdb";
	    cen_pose_.dump_pdb(ss.str());
	}
    }

private:
    protocols::abinitio::ClassicAbinitioOP abinitio_;
    protocols::simple_filters::PDDFScoreFilterOP the_filter_;
    std::string sequence_;
    core::pose::Pose cen_pose_;
    core::kinematics::MoveMapOP  movemap_;
    core::fragment::FragSetOP fragset_large_;
    core::fragment::FragSetOP fragset_small_;

    void setup_pose() {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

        core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("centroid");

	std::string sequence;
	if (option[in::file::fasta].user()) {
	    sequence = core::sequence::read_fasta_file(option[in::file::fasta]()[1])[1]->sequence();
		core::pose::make_pose_from_sequence(cen_pose_, sequence,*rsd_set);
	}

	if (option[in::file::s].user()) {

		core::pose::PoseOP tmp_pose(new core::pose::Pose);
		std::string fn = option[in::file::s](1);
		core::import_pose::pose_from_pdb(*tmp_pose, fn);

		sequence = tmp_pose->sequence();
		core::pose::make_pose_from_sequence(cen_pose_, sequence,*rsd_set);
		for (Size i = 1; i <= sequence.size(); ++i) {
			cen_pose_.set_phi(i, tmp_pose->phi(i));
			cen_pose_.set_psi(i, tmp_pose->psi(i));
			cen_pose_.set_omega(i, tmp_pose->omega(i));
		}
	}
    }
};




////////////////////////////////////////////////////////
int main( int argc, char * argv [] ) {
try {
    protocols::abinitio::ClassicAbinitio::register_options();
    protocols::abinitio::AbrelaxApplication::register_options();
    register_options();
    devel::init( argc, argv );

    TryScoreFilter job;
    job.run();
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }

    return 0;
}

