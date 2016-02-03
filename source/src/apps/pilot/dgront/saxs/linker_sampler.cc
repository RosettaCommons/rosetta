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
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>

#include <protocols/moves/Mover.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>


#include <core/fragment/FragmentIO.hh>

#include <core/fragment/FragSet.hh>

#include <protocols/evaluation/TimeEvaluator.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "linker_sampler" );

OPT_1GRP_KEY( IntegerVector, assembly, link_centers )
OPT_1GRP_KEY( IntegerVector, assembly, link_lengths )

using namespace core;

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(in::file::fasta);
  OPT(in::file::native);
  OPT(out::nstruct);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
  OPT(score::saxs::q_min);
  OPT(score::saxs::q_max);
  OPT(score::saxs::q_step);
  OPT(score::saxs::custom_ff);
  OPT(score::saxs::ref_cen_spectrum);
  OPT(score::saxs::ref_fa_spectrum);

  OPT(relax::fast);

  NEW_OPT( assembly::link_centers, "the center between the domains",0);
  NEW_OPT( assembly::link_lengths, "the size of each linker",0);
}


class DomainAssemblerNDocker {

public:

    /// @brief Constructor parses command line and sets up the job
    DomainAssemblerNDocker() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace std;

        core::chemical::ResidueTypeSetCAP rsd_set_cen = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
        core::chemical::ResidueTypeSetCAP rsd_set_fa = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::import_pose::pose_from_file( cen_pose_, option[in::file::s]()[1].name(), core::import_pose::PDB_file);
	core::util::switch_to_residue_type_set(cen_pose_, core::chemical::CENTROID);


//------------- set up sampling  -----------------
	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_bb(false);
	loop_pos_ = option[assembly::link_centers]();
	loop_len_ = option[assembly::link_lengths]();

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

	has_native_ = false;
        if ( option[ in::file::native ].user() ) {
    	    core::import_pose::pose_from_file( native_pose_, option[ in::file::native ]().name() , core::import_pose::PDB_file);
    	    has_native_ = true;
	}
    }

// ------------- Run the protocol ------------------
    void run() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	for(Size i=0;i<(Size)option[ out::nstruct ];i++) {
// ------------- Set loop conformation to -150,150,180 -----------------
	    for(Size j=1;j<=loop_pos_.size();j++) {
		Size pos = loop_pos_[j];
		Size len = loop_len_[j];
		for(Size i=pos-len;i<=pos+len;i++) {
		    movemap_->set_bb(i,true);
		    cen_pose_.set_phi(i, -150.0);
		    cen_pose_.set_psi(i, 150.0);
		    cen_pose_.set_omega(i, 180.0);
		}
	    }
//	    cen_pose_.dump_pdb("start_"+ObjexxFCL::lead_zero_string_of(i+1,5)+".pdb");

// ------------- Run Abinitio protocol --------------------------
	    protocols::abinitio::AbrelaxApplication abrelax_app;
	    protocols::evaluation::TimeEvaluatorOP run_time_( NULL );
	    abrelax_app.add_evaluation( run_time_ = new protocols::evaluation::TimeEvaluator );
	    if ( run_time_ ) run_time_->reset();
	    do {
		protocols::abinitio::ClassicAbinitio abinitio( fragset_large_, fragset_small_ ,movemap_ );
		abinitio.init( cen_pose_ );
		abinitio.apply( cen_pose_ );
	    } while( !abrelax_app.check_filters( cen_pose_ ) );

// ------------- Switch to all-atom --------------------------
	    core::pose::Pose fa_pose(cen_pose_);
	    core::util::switch_to_residue_type_set( fa_pose, core::chemical::FA_STANDARD);
//	    cen_pose_.conformation().detect_bonds();//apl fix this !

// ------------- Relax  --------------------------
	    if ( option[ OptionKeys::relax::fast].user() ) {
	    // core::scoring::ScoreFunctionOP scorefxn_fa( core::scoring::get_score_function() );
    		core::scoring::ScoreFunctionOP scorefxn_fa = core::scoring::get_score_function();
		protocols::moves::MoverOP relaxer = new protocols::relax::FastRelax( scorefxn_fa );
		relaxer->apply( fa_pose );
	    }
// ------------- output  --------------------------
    	    if ( !option[ out::nooutput ]() ) {
        	core::io::silent::SilentStructOP ss
            	    = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
        	ss->fill_struct( fa_pose );
        	ss->set_decoy_tag("S_"+ObjexxFCL::lead_zero_string_of(i+1,5));

        	if ( has_native_ ) {
            	    core::Real CA_rmsd = core::scoring::CA_rmsd( native_pose_, fa_pose );
            	    core::Real CA_maxsub = core::scoring::CA_maxsub( native_pose_, fa_pose );
            	    core::Real CA_gdtmm = core::scoring::CA_gdtmm( native_pose_, fa_pose );
            	    ss->add_energy( "rmsd", CA_rmsd );
            	    ss->add_energy( "maxsub", CA_maxsub );
            	    ss->add_energy( "gdtmm", CA_gdtmm );
        	}
        	ss->add_energy( "time", run_time_->apply(fa_pose) );
        	sfd_out.write_silent_struct( *ss, option[ out::file::silent ]() );
    	    }
        }
    }

private:

    utility::vector1<Size> loop_pos_;
    utility::vector1<Size> loop_len_;
    core::io::silent::SilentFileData sfd_out;
    std::string sequence_;
    core::pose::Pose cen_pose_;
    core::pose::Pose fa_pose_;
    core::kinematics::MoveMapOP  movemap_;
    core::fragment::FragSetOP fragset_large_;
    core::fragment::FragSetOP fragset_small_;
    core::pose::Pose native_pose_;
    bool has_native_;
};


////////////////////////////////////////////////////////
int main( int argc, char * argv [] ) {
    try {
    protocols::abinitio::ClassicAbinitio::register_options();
    protocols::abinitio::AbrelaxApplication::register_options();
    register_options();
    devel::init( argc, argv );

    DomainAssemblerNDocker job;
    job.run();
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
    return 0;
}

