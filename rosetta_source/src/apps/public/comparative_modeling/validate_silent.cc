// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/comparative_modeling/validate_silent.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <fstream>
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>

using namespace std;

void read_target_sequence(const string& filename, string* target_sequence) {
  assert(target_sequence);

  ifstream in(filename.c_str());
  if (!in.is_open()) {
    utility_exit_with_message("Failed to open input silent file");
  }

	// SEQUENCE: MNDDVDIQQSYPFSIETMPVPKKLKVGETAEIRCQLH
	string line;
  getline(in, line);
  in.close();

	// MNDDVDIQQSYPFSIETMPVPKKLKVGETAEIRCQLH
	*target_sequence = line.substr(10);
}

int main(int argc, char* argv[]) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
  devel::init(argc, argv);

  if (!option[in::file::silent].user()) {
    utility_exit_with_message("Failed to provide required argument -in:file:silent");
  }

  if (!option[out::file::silent].user()) {
    utility_exit_with_message("Failed to provide required argument -out:file:silent");
  }

  const string input_file = option[in::file::silent]()[1];
  const string output_file = option[out::file::silent]();

  string target_sequence;
  read_target_sequence(input_file, &target_sequence);

	SilentFileData sfd_in, sfd_out;
	sfd_in.read_file(input_file);

	size_t num_good = 0, num_bad = 0;

	utility::vector1<string> tags = sfd_in.tags();
	for (utility::vector1<string>::const_iterator i = tags.begin(); i != tags.end(); ++i) {
		SilentStructOP decoy = sfd_in[*i];
		string sequence = decoy->sequence().one_letter_sequence();

		if (sequence == target_sequence) {
			sfd_out.write_silent_struct(*decoy, output_file, false);
			++num_good;
		} else {
			cerr << "Removed tag: " << *i << " seq: " << sequence << endl;
			++num_bad;
		}
	}

	cout << "good: " << num_good << " bad: " << num_bad << endl;
}
