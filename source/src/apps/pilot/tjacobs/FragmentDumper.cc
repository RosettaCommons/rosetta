// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentDumper.cc
///
/// @brief

/// @author Tim Jacobs

// Core headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Basic headers
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

//utility
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>

//C++ headers
#include <iostream>
#include <string>
#include <stdlib.h>

namespace FragmentDumper {
	basic::options::IntegerOptionKey const num_frags_to_print( "num_frags_to_print" );
}

int
main( int argc, char * argv [] )
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::fragment;
	
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	using namespace utility::file;
	
	// initialize devel
	option.add( FragmentDumper::num_frags_to_print, "The number of fragments to dump to pdbs");
	devel::init(argc, argv);
	
	//get fragments
	utility::vector1<core::fragment::FragSetOP> frag_sets;
	if(option[ OptionKeys::loops::frag_files ].user())
	{
		FileVectorOption frag_files( option[ OptionKeys::loops::frag_files ] );
		for ( Size i = 1; i <= frag_files.size(); ++i ) {
			frag_sets.push_back(FragmentIO().read_data( frag_files[i] ));
		}
	}
	else
	{
		utility_exit_with_message("Must provide fragments files with loops::frag_files");
	}
	
	core::Size num_frags_to_print = option[FragmentDumper::num_frags_to_print].def(10);
	
	cout << "Number of fragment sets: " << frag_sets.size() << endl;
		
	core::Size num_printed_frags=0;
	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	for(Size frag_sets_index=1; frag_sets_index<=frag_sets.size() ; ++frag_sets_index)
	{
		fragment::FragSetOP cur_frag_set=frag_sets[frag_sets_index];
		cout << "Number of fragment in set " << frag_sets_index << ": " << cur_frag_set->size() << endl;
		for(fragment::FrameIterator it = cur_frag_set->begin(); it!=cur_frag_set->end(); ++it)
		{
			for(Size i = 1; i<=it->nr_frags(); ++i)
			{
				core::pose::Pose fragment_pose;
				it->fragment_as_pose(i, fragment_pose, restype_set);
		
				std::stringstream filename;
				filename << "fragment_" << num_printed_frags << ".pdb";
				fragment_pose.dump_pdb(filename.str());
				
				++num_printed_frags;
				if(num_printed_frags >= num_frags_to_print)
				{
					//we are done
					return 0;
				}
			}
		}
	}
	utility_exit_with_message("Printed all fragments");
	return 0;
}
