// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/jeliazkov/calc_irms_despite_mismatch
/// @brief pilot app for comparing interface rmsds where partner residues might not be equivalent
/// sample usage: calc_irms_despite_mismatch.macosclangrelease
//                    -l 1ahw.list -native 1ahw.pdb -aln_chain G -inteface_jump 1 -out:file:score_only score.sc
/// @author Jeliazko Jeliazkov (jeliazkov@jhu.edu)


#include <protocols/jd2/internal_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <devel/init.hh>

#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

#include <core/pose/PDBInfo.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>

/// my headers
#include <string>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


static basic::Tracer TR( "apps.pilot.jeliazkov.cal_irms_despite_mismatch" );

OPT_KEY( String, aln_chain )
OPT_KEY( Integer, interface_jump )

class CalcIrms : public protocols::moves::Mover {
public:
	// default constructor
	CalcIrms(){

		using namespace basic::options;

		if ( not option[ OptionKeys::interface_jump].user() ) {
			utility_exit_with_message("No interface jump given by -inteface_jump flag! Exiting!");
		}

		if ( not  option[ OptionKeys::aln_chain].user() ) {
			utility_exit_with_message("No chain given by -aln_chain flag! Exiting!");
			return;
		}

		aln_chain_ = option[ OptionKeys::aln_chain ];
		interface_jump_ = option[ OptionKeys::interface_jump ];

		if ( aln_chain_.size() > 1 ) {
			// exit gracefully?
			utility_exit_with_message("Multiple chains given for alignment! Not possible at the moment!");
			return;
		}


	};
	// destructor
	~CalcIrms() override= default;

	void apply( core::pose::Pose & model_pose ) override
	{

		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace protocols::analysis;

		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

		// sequence-based mapping
		id::SequenceMapping seq_map = sequence::map_seq1_seq2(
			utility::pointer::make_shared< sequence::Sequence >(native_pose_),
			utility::pointer::make_shared< sequence::Sequence >(model_pose)
		);

		// for debugging, what's in the chain
		for ( core::Size i=1; i <= native_pose_.size(); ++i ) {
			if ( is_aln_chain_[i] ) {
				TR << native_pose_.pdb_info()->pose2pdb(i) << std::endl;
				TR << "Native residue: " << native_pose_.sequence()[i] << std::endl;
				TR << "Model  residue: " << model_pose.sequence()[seq_map[i]] << std::endl;
			}
		}

		// backbone align on overlapping set, but score model first (just in case)
		scorefxn_->score(model_pose);
		Real chain_rms = rmsd_with_super_subset( native_pose_, model_pose, is_aln_chain_, seq_map, is_protein_backbone_including_O );
		TR << "Aligned chain " << aln_chain_ << " with RMS: " << chain_rms << std::endl;

		// for debugging, what's at the interface
		for ( core::Size i=1; i <= native_pose_.size(); ++i ) {
			if ( is_native_interface_(i) ) {
				if ( native_pose_.sequence()[i] != model_pose.sequence()[seq_map[i]] ) {
					TR << "Warning: sequence mismatch at interface: " << std::endl;
					TR << native_pose_.pdb_info()->pose2pdb(i) << std::endl;
					TR << "Native residue: " << native_pose_.sequence()[i] << std::endl;
					TR << "Model  residue: " << model_pose.sequence()[seq_map[i]] << std::endl;
				}
			}
		}

		Real irms_ns( 0.0 );
		Real irms_s( 0.0 );

		irms_ns = rmsd_no_super_subset( native_pose_, model_pose, is_native_interface_, seq_map, is_heavyatom );
		irms_s = rmsd_with_super_subset( native_pose_, model_pose, is_native_interface_, seq_map, is_heavyatom );


		setPoseExtraScore(model_pose, "irms_ns", irms_ns);
		setPoseExtraScore(model_pose, "irms_s", irms_s);

		// also get interaction energy from docking metrics -- vector1_int of length 1, with value 2 (I think)
		//TR << "What am I?? " << utility::vector1_int(1,2) << std::endl;
		Real Isc = protocols::docking::calc_interaction_energy(model_pose, scorefxn_, utility::vector1_int(1,interface_jump_));
		setPoseExtraScore(model_pose, "I_sc", Isc);

		// also also get interaction energy by repacking away?
		InterfaceAnalyzerMover iam = InterfaceAnalyzerMover(  interface_jump_, false, nullptr, false, false, true );
		iam.apply_const(model_pose);

		setPoseExtraScore(model_pose, "dG", iam.get_separated_interface_energy());

		return;
	} // apply

	std::string get_name() const override { return "CalcIrms"; }


	protocols::moves::MoverOP
	fresh_instance() const override {
		return utility::pointer::make_shared< CalcIrms >();
	}


	bool
	reinitialize_for_each_job() const override { return false; }


	bool
	reinitialize_for_new_input() const override { return false; }

	void
	setup_native() {

		using namespace basic::options;
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		import_pose::pose_from_file( native_pose_, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file );

		scorefxn_ = get_score_function();
		scorefxn_->score(native_pose_);

		// calculate using interface from crystal -- assumes those residues are present in model!
		protocols::scoring::Interface interface( interface_jump_ );
		interface.distance( 10.0 ); // from protocols/docking/util or something
		interface.calculate( native_pose_ );

		is_native_interface_ = ObjexxFCL::FArray1D_bool( native_pose_.size(), false );
		is_aln_chain_ = ObjexxFCL::FArray1D_bool( native_pose_.size(), false );

		for ( Size i = 1; i <= native_pose_.size(); ++i ) {
			// note () access to this array implies 1-index, whereas [] implies 0
			if ( interface.is_interface( i ) ) is_native_interface_( i ) = true;
			if ( native_pose_.pdb_info()->chain(i) == aln_chain_[0] ) is_aln_chain_( i ) = true;
		}

		TR << "The alignment residues: " << std::endl;
		for ( Size i = 1; i <= native_pose_.size(); ++i ) {
			TR << is_aln_chain_(i) << " ,";
		}
		TR << std::endl;
	}

private:
	// Things we should only set once, while we load *many* models
	core::pose::Pose native_pose_;
	ObjexxFCL::FArray1D_bool is_native_interface_;
	ObjexxFCL::FArray1D_bool is_aln_chain_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string aln_chain_;
	core::Size interface_jump_;

};

using CalcIrmsOP = utility::pointer::shared_ptr<CalcIrms>;


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		protocols::jd2::register_options();

		NEW_OPT( aln_chain, "Chain for alignment.", "" );
		NEW_OPT( interface_jump, "Jump number for interface metric calculations.", 1 );

		// initialize core
		devel::init(argc, argv);

		CalcIrmsOP calc_irms = utility::pointer::make_shared< CalcIrms >();
		calc_irms->setup_native();
		protocols::jd2::JobDistributor::get_instance()->go( calc_irms );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
