// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Mike Tyka
/// @author Dominik Gront

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/relax/ClassicRelax.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <devel/init.hh>

#include <core/types.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/rms_util.hh>

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

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>


static thread_local basic::Tracer TR( "protocols::DommainAssemblerNDocker" );

OPT_1GRP_KEY( String, assembly, cfg )

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
  OPT(score::saxs::ref_spectrum);

  NEW_OPT( assembly::cfg, "input config file", "" );
}


class DomainAssemblerNDocker {

public:

    /// @brief Constructor parses command line and sets up the job
    DomainAssemblerNDocker() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace std;

//-------------set up sequence -----------------
	if (option[in::file::fasta].user())
	    sequence_ = core::sequence::read_fasta_file(option[in::file::fasta]()[1])[1]->sequence();
	else
	    utility_exit_with_message("You must give a sequence using in::file::fasta flag !");

//------------- set up sampling  -----------------
	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_bb(true);

//------------- set up domains -----------------
	if(!option[assembly::cfg].user())
	    utility_exit_with_message("You must give a config file that defines domains !");
	setup_cen_domains(utility::io::izstream( option[ assembly::cfg ]() ) );

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
    	    core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
    	    core::import_pose::pose_from_pdb( native_pose_, *rsd_set, option[ in::file::native ]() );
    	    has_native_ = true;
	}
    }

    /// @brief
    void run() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	for(Size i=0;i<(Size)option[ out::nstruct ];i++) {

	protocols::abinitio::AbrelaxApplication abrelax_app;
	do {
	    protocols::abinitio::ClassicAbinitio abinitio( fragset_large_, fragset_small_ ,movemap_ );
	    abinitio.init( cen_pose_ );
	    abinitio.apply( cen_pose_ );
	} while( !abrelax_app.check_filters( cen_pose_ ) );

        // output
        if ( !option[ out::nooutput ]() ) {
            core::io::silent::SilentStructOP ss
                = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
            ss->fill_struct( cen_pose_ );
            ss->set_decoy_tag("S_0000001");

            if ( has_native_ ) {
                core::Real CA_rmsd
                        = core::scoring::CA_rmsd( native_pose_, cen_pose_ );
                core::Real CA_maxsub
                        = core::scoring::CA_maxsub( native_pose_, cen_pose_ );
                core::Real CA_gdtmm
                        = core::scoring::CA_gdtmm( native_pose_, cen_pose_ );
                ss->add_energy( "rmsd", CA_rmsd );
                ss->add_energy( "maxsub", CA_maxsub );
                ss->add_energy( "gdtmm", CA_gdtmm );
            }
            sfd_out.write_silent_struct( *ss, option[ out::file::silent ]() );
        }
        }
    }

private:

    core::io::silent::SilentFileData sfd_out;
    std::string sequence_;
    core::pose::Pose cen_pose_;
    core::kinematics::MoveMapOP  movemap_;
    core::fragment::FragSetOP fragset_large_;
    core::fragment::FragSetOP fragset_small_;
    core::pose::Pose native_pose_;
    bool has_native_;

    /// @brief reads PDB files for domains and copies dihedrals onto a blank pose
    void setup_cen_domains(std::istream & input) {

	std::string line;
        core::chemical::ResidueTypeSetCAP rsd_set_cen = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
        core::chemical::ResidueTypeSetCAP rsd_set_fa = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	core::pose::make_pose_from_sequence(cen_pose_, sequence_,
			*rsd_set_cen);
	TR.Debug << "Creating an empty pose with "<<cen_pose_.total_residue()<<" residues"<<std::endl;
	for (Size pos = 1; pos <= cen_pose_.total_residue(); pos++) {
	    if (!cen_pose_.residue(pos).is_protein())
		continue;
	    cen_pose_.set_phi(pos, -150);
	    cen_pose_.set_psi(pos, 150);
	    cen_pose_.set_omega(pos, 180);
	}

	while( getline( input, line ) ) {
	    if ( line.substr(0,1) == "#" ) continue;
    	    std::istringstream line_stream( line );
    	    if(line.length() < 3)
    	      continue;
	    TR.Debug << "Processing line " << line << std::endl;

	    std::string fname;
	    Size resi;
	    Size resj;
	    line_stream >> fname >> resi >> resj;
	    core::pose::Pose tmp_pose;
	    core::import_pose::pose_from_pdb( tmp_pose, *rsd_set_fa,  fname);

	    for(Size ir = 1; ir <= tmp_pose.total_residue(); ir++ ) {
		cen_pose_.set_phi( ir + resi, tmp_pose.phi( ir ) );
		cen_pose_.set_psi( ir + resi, tmp_pose.psi( ir ) );
		cen_pose_.set_omega( ir + resi, tmp_pose.omega( ir ) );
		movemap_->set_bb(ir + resi,false);
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

    DomainAssemblerNDocker job;
    job.run();
} catch ( utility::excn::EXCN_Base const & e ) {
                         std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                            }
    return 0;
}

