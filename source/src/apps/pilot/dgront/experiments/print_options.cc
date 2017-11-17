// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief A very simple case of SAXS score that reads a PDB and prints the score value
/// @author Dominik Gront

#include <devel/init.hh>
#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

}


////////////////////////////////////////////////////////
int main( int argc, char * argv [] ) {
    try {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    register_options();
    devel::init( argc, argv );

// grep OPT fragment_picker.cc | sed 's/OPT//' | tr -d '; \t' | awk '{printf "&%s,\n",$i}'

    static const OptionKey* opts[] = {
&(in::file::native),
&(in::file::s),
&(in::file::xyz),
&(in::file::fasta),
&(in::file::pssm),
&(in::file::checkpoint),
&(in::file::talos_phi_psi),
&(in::file::torsion_bin_probs),
&(in::path::database),
&(frags::scoring::config),
&(frags::scoring::profile_score),
&(frags::ss_pred),
&(frags::n_frags),
&(frags::n_candidates),
&(frags::frag_sizes),
&(frags::write_ca_coordinates),
&(frags::allowed_pdb),
&(frags::denied_pdb),
&(frags::describe_fragments),
&(frags::keep_all_protocol),
&(frags::bounded_protocol),
&(frags::quota_protocol),
&(frags::picking::selecting_rule),
&(frags::picking::quota_config_file),
&(frags::picking::query_pos),
&(constraints::cst_file),
&(out::file::frag_prefix)
    };

    for(int i=0;i<28;i++) {
          basic::options::Option *o = &option[*opts[i]];
//	  std::cout << o->name() << " " << o->description() << std::endl;

	std::cout << "\t<tr>\n";
	std::cout << "\t\t<td width=30%>" << o->name() << "</td>\n";
	std::cout << "\t\t<td width=60%>" << o->description() << "</td>\n";
	std::cout << "\t\t<td width=30%>" << o->name() << "</td>\n";
	std::cout << "\t</tr>\n";
    }
    } catch (utility::excn::Exception const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
    return 0;
}

