// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/util.cc
/// @brief Utility functions for the tensorflow manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <basic/tensorflow_manager/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "basic.tensorflow_manager.util" );

#include <utility/string_util.hh>

namespace basic {
namespace tensorflow_manager {

/// @brief Writes the following message:
/// The <module_name> requires compilation with Tensorflow support.  To compile with Tensorflow
/// support...
/// (Followed by instructions terminating in a carriage return character.)
std::string
get_tensorflow_compilation_instructions(
	std::string const & module_name,
	bool const use_single_quotations
) {
	std::string outstring(
		"The " + module_name + " requires compilation with Tensorflow support.  To compile with Tensorflow "
		"support:\n\n"
		"1.  Download the Tensorflow 1.15 precompiled libraries for your operating system from one of the "
		"following.  (Note that GPU versions require CUDA drivers; see https://www.tensorflow.org/install/lang_c "
		"for more information.)\n"
		"\tLinux/CPU: https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-1.15.0.tar.gz\n"
		"\tLinux/GPU: https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-gpu-linux-x86_64-1.15.0.tar.gz\n"
		"\tWindows/CPU: https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-windows-x86_64-1.15.0.zip\n"
		"\tWindows/GPU: https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-gpu-windows-x86_64-1.15.0.zip\n"
		"\tMacOS/CPU: https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-darwin-x86_64-1.15.0.tar.gz\n"
		"\tMacOS/GPU: None available.\n\n"
		"2.  Unzip/untar the archive into a suitable directory (~/mydir/ is used here as an example), and add the following environment variables:\n"
		"\tLinux, Windows:\n\t\tLIBRARY_PATH=$LIBRARY_PATH:~/mydir/lib\n\t\tLD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/mydir/lib\n"
		"\tMacOS:\n\t\tLIBRARY_PATH=$LIBRARY_PATH:~/mydir/lib\n\t\tDYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/mydir/lib\n\n"
		"3.  Edit your user.settings file (Rosetta/main/source/tools/build/user.settings), and uncomment (i.e. remove the "
		"octothorp from the start of) the following lines:\n"
		"\timport os\n"
		"\t\t'program_path'  : os.environ['PATH'].split(':'),\n"
		"\t\t'ENV' : os.environ,\n\n"
		"4.  Compile Rosetta, appending extras=tensorflow (for CPU-only) or extras=tensorflow_gpu (for GPU) to your scons "
		"command.  For example:\n"
		"\t./scons.py -j 8 mode=release extras=tensorflow bin\n"
	);
	if ( use_single_quotations == false ) {
		utility::replace_in( outstring, "'", "\"" );
	}
	return outstring;
}

} //tensorflow_manager
} //basic


