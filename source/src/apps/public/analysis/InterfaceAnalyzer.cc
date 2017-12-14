// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/analysis/InterfaceAnalyzer.cc
/// @brief Q&D protocol to run InterfaceAnalyzerMover as protocol
/// @author Steven Lewis, Bryan Der, Ben Stranges

// Unit Headers
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/Mover.hh>

#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.public.analysis.InterfaceAnalyzer" );

//define local options
basic::options::IntegerOptionKey const jumpnum("jumpnum");
basic::options::BooleanOptionKey const compute_packstat("compute_packstat");
basic::options::BooleanOptionKey const tracer_data_print("tracer_data_print");
basic::options::StringVectorOptionKey const fixedchains( "fixedchains" );
basic::options::BooleanOptionKey const compute_interface_sc("compute_interface_sc"); // sc == shape complementarity
basic::options::BooleanOptionKey const pack_input("pack_input");
basic::options::BooleanOptionKey const pack_separated("pack_separated");
basic::options::BooleanOptionKey const use_jobname("use_jobname");
basic::options::BooleanOptionKey const add_regular_scores_to_scorefile("add_regular_scores_to_scorefile");
basic::options::BooleanOptionKey const use_resfile("use_resfile");
basic::options::StringOptionKey const interface("interface");

// mover deffinition
class IAMover : public protocols::moves::Mover {
public:

	IAMover();

	~IAMover() override = default;

	void apply( core::pose::Pose& pose ) override;


	std::string
	get_name() const override {
		return "IAMover";
	}

	bool reinitialize_for_each_job() const override { return true; }

	bool reinitialize_for_new_input() const override { return true; }
	protocols::moves::MoverOP fresh_instance() const override {return protocols::moves::MoverOP( new IAMover );}


	void assign_IA_mover(core::pose::Pose & pose);

private:
	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;
	core::scoring::ScoreFunctionOP scorefxn_;
};

IAMover::IAMover() : scorefxn_(core::scoring::get_score_function()) {}

//assign the correct constructor for the mover, and figure out the multichain assignment for that ctor
void IAMover::assign_IA_mover(core::pose::Pose & pose){

	core::Size const num_chains(pose.conformation().num_chains());

	//shared booleans - making local copies for readability
	bool const tracer = basic::options::option[ tracer_data_print ].value();
	bool const comp_packstat = basic::options::option[ compute_packstat ].value();
	bool const pack_in = basic::options::option[ pack_input ].value();
	bool const pack_sep = basic::options::option[ pack_separated ].value();
	bool const jobname = basic::options::option[ use_jobname ].value();

	//if 2 chains, or no multichain ctor, use interface_jump constructor
	if ( (num_chains <= 2) ) {
		TR << "Computing interface between two chains in pose" << std::endl;
		core::Size const interface_jump = basic::options::option[ jumpnum ].value();
		IAM_ = protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(
			interface_jump,
			tracer,
			scorefxn_,
			comp_packstat,
			pack_in,
			pack_sep,
			jobname
			) );
	} else if ( basic::options::option[fixedchains].active() ) {
		utility::vector1<std::string> fixed_chains_string (basic::options::option[fixedchains].value());
		//parse the fixed chains to figure out pose chain nums
		std::set< int > fixed_chains; //This is a set of the CHAIN IDs, not residue ids
		TR << "Fixed chains are: " ;
		for ( core::Size j = 1; j <= fixed_chains_string.size(); ++j ) {
			char this_chain (fixed_chains_string[ j ][0]);
			for ( core::Size i = 1; i<=pose.size(); ++i ) {
				if ( pose.pdb_info()->chain( i ) == this_chain ) {
					fixed_chains.insert( pose.chain(i) );
					break; //once we know something about the chain we can skip - we just need the chain id
				}
			}
			TR << this_chain << ", ";
		}
		TR << "these will be moved together." << std::endl;

		IAM_ = protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(
			fixed_chains,
			tracer,
			scorefxn_,
			comp_packstat,
			pack_in,
			pack_sep,
			jobname
			) );
	} else if ( basic::options::option[interface].active() ) {
		std::string dock_chains = basic::options::option[interface].value();
		TR << "Using interface definition: "<<dock_chains <<std::endl;
		IAM_ = protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(
			dock_chains,
			tracer,
			scorefxn_,
			comp_packstat,
			pack_in,
			pack_sep,
			jobname
			) );
	} else {
		utility_exit_with_message("More than two chains present but no -fixedchains or -interface declared.   Aborting.");
	}
	//IAM_->set_use_resfile(basic::options::option[use_resfile].value());

	if ( ! basic::options::option[ compute_interface_sc ] ) {
		IAM_->set_compute_interface_sc( false );
	}

	return;
} //end assign_IA_mover

///begin apply
void IAMover::apply( core::pose::Pose & pose ) {

	//check to make sure there are enough chains
	if ( pose.conformation().num_chains() < 2 ) {
		TR.Error << "pose has only one chain, skipping" << std::endl;
		set_last_move_status( protocols::moves::FAIL_BAD_INPUT);
		return;
	}

	//fill the interface analyzer mover
	assign_IA_mover( pose );

	//now apply and get cool data and stuff
	IAM_->apply(pose);
	//flesh out scores for scorefile, if desired
	if ( basic::options::option[add_regular_scores_to_scorefile].value() ) ( *scorefxn_)(pose);

	return;
}//end apply

int
main( int argc, char* argv[] )
{
	try {

		using basic::options::option;
		option.add( jumpnum, "jump between chains of interface" ).def(1);
		option.add( compute_packstat, "compute packstat (of interface residues only)" ).def(false);
		option.add( tracer_data_print, "print to tracer, not scorefile" ).def(false);
		option.add( fixedchains, "Which chain(s) is/are moved away from the others, for 3 or more chains" );
		option.add( compute_interface_sc, "Compute the shape complementarity score").def(true);
		option.add( pack_input, "pack the input pose").def(false);
		option.add( add_regular_scores_to_scorefile, "adds default scores to scorefile").def(false);
		option.add( pack_separated, "pack the separated chains at the separated dG phase").def(false);
		option.add( use_jobname, "appended _0001 job name for output instead of input pose name").def(false);
		option.add( use_resfile, "use a resfile during the packing stages").def(false);
		option.add( interface, "dock_chains interface definition, optional, ex LH_A.  Can handle any number of chains. ");
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new IAMover() ));

		TR << "************************d**o**n**e**************************************" << std::endl;
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
