// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/sequence/SequenceProfile.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Utility headers
#include <utility/file/FileName.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <utility/excn/Exceptions.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "main" );


class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, profile )
OPT_KEY( File, vall )
//OPT_1GRP_KEY( File, cluster, out )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( profile, "blast .checkpoint sequence profile","");
	NEW_OPT( vall, "vall","");
	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::xyz);
	OPT(in::file::fasta);
	OPT(in::file::pssm);
	OPT(in::file::checkpoint);
	OPT(in::file::talos_phi_psi);
	OPT(in::file::torsion_bin_probs);
	OPT(in::path::database);

	OPT(frags::scoring::config);
	OPT(frags::scoring::profile_score);
	OPT(frags::ss_pred);
	OPT(frags::n_frags);
	OPT(frags::n_candidates);
	OPT(frags::frag_sizes);
	OPT(frags::write_ca_coordinates);
	OPT(frags::allowed_pdb);
	OPT(frags::denied_pdb);
	OPT(frags::describe_fragments);
	OPT(frags::keep_all_protocol);
	OPT(frags::bounded_protocol);
	OPT(frags::quota_protocol);
	OPT(frags::picking::selecting_rule);
	OPT(frags::picking::selecting_scorefxn);
	OPT(frags::picking::quota_config_file);
	OPT(frags::picking::query_pos);

	OPT(constraints::cst_file);

	OPT(out::file::frag_prefix);

	OPT(frags::nonlocal_pairs);
	OPT(frags::nonlocal::min_contacts_per_res);
	OPT(frags::contacts::min_seq_sep);
	OPT(frags::contacts::dist_cutoffs);
	OPT(frags::contacts::centroid_distance_scale_factor);
	OPT(frags::contacts::type);

}

class FragmentBasedDesign : public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef protocols::frag_picker::FragmentPicker FragmentPicker;
	typedef protocols::frag_picker::FragmentPickerOP FragmentPickerOP;
	typedef core::sequence::SequenceProfile SequenceProfile;
	typedef core::sequence::SequenceProfileOP SequenceProfileOP;
public:
	FragmentBasedDesign(utility::file::FileName seqprof, utility::file::FileName vall) {
		// SequenceProfileOP profile = new SequenceProfile();
		//  picker_ = new FragmentPicker();
		//  picker_->parse_command_line();
		//  profile->read_from_checkpoint(seqprof);
		//  picker_->set_query_seq(profile);
		//picker_->read_vall(vall);
	}

	virtual void apply(Pose & pose) {

		picker_ = new FragmentPicker();
		picker_->parse_command_line();
		picker_->add_query_ss(pose.secstruct(), "bla");
		std::string sequence = pose.sequence();
		picker_->set_query_seq(sequence);
		picker_->pick_candidates();
	}

	virtual std::string get_name() const {return "FragmentBasedDesign";}
	//virtual protocols::moves::MoverOP fresh_instance() const { return new FragmentBasedDesign();}

	FragmentPickerOP picker_;
};

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main( int argc, char** argv ) {
	try {

		ThisApplication::register_options();
		devel::init( argc, argv );
		// mover
		protocols::moves::MoverOP protocol;
		protocol = new FragmentBasedDesign( option[profile], option[vall]);

		// run
		protocols::jd2::JobDistributor::get_instance()->go( protocol );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}
