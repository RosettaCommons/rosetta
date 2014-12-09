// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include <devel/init.hh>

// utility headers
#include <utility/vector1.hh>

#include <core/sequence/Sequence.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>

#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/io/mpistream.hh>


static thread_local basic::Tracer trace( "picker" );

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void register_options() {

	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::xyz);
	OPT(in::file::fasta);
	OPT(in::file::pssm);
	OPT(in::file::checkpoint);
	OPT(in::file::talos_phi_psi);
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
	OPT(frags::picking::quota_config_file);
	OPT(frags::picking::query_pos);
//	OPT(frags::p_value_selection);
	OPT(in::file::torsion_bin_probs);
	OPT(out::file::frag_prefix);
	OPT(constraints::cst_file);
	OPT(in::path::database);

//	OPT(rdc::correct_NH_length);
//	OPT(rdc::reduced_couplings);
}

int main(int argc, char * argv[]) {
        try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	register_options();
	devel::init(argc, argv);

	if (option[in::file::native].user()) {
		trace.Debug << option[in::file::native]() << std::endl;
	}

	//---------- Set up a picker.
	FragmentPickerOP pickIt;
	if(option[frags::p_value_selection]() == true)
	    pickIt = new FragmentPicker("PValuedFragmentScoreManager");
	else
	    pickIt = new FragmentPicker();
	pickIt->parse_command_line();
	trace << "After setup; size of a query is: " << pickIt->size_of_query()
			<< std::endl;

	//-------- Trata ta ta, tra ta... (picking fragment candidates)
	trace << "Picking candidates" << std::endl;

	if (option[frags::picking::quota_config_file].user() || option[frags::quota_protocol].user() ) {
	    trace << "Running quota protocol" << std::endl;
	    pickIt->quota_protocol();
	} else {
	    if (option[frags::keep_all_protocol].user()) {
		trace << "Running keep-all protocol" << std::endl;
		pickIt->keep_all_protocol();
	    } else {
		trace << "Running bounded protocol" << std::endl;
		pickIt->bounded_protocol();
	    }
	}

	basic::prof_show();
        } catch ( utility::excn::EXCN_Base const & e ) {
                                  std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                      }
            return 0;

}

