// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @author Steven Combs
/// @brief a simple application to generate statistics for orbital interactions.
/// Command line: generate_orbital_statistics -database <database> -l <pdb_list> -orbitals:Hpol -orbitals:sc_stats
/// note that you can do -orbitals:Hpol for polar hydrogen - orbital interactions or -orbitals:Haro for aromatic hydrogen - orbital interactions

//TODO add stat calc for bb orbitals and orbital-orbital interactions

#include <basic/options/util.hh>
#include <devel/init.hh>


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/orbitals/OrbitalsStatistics.hh>

#include <numeric/histograms/TwoDHistogram.hh>
#include <core/chemical/orbitals/OrbitalType.hh>

#include <utility/file/FileName.hh>


utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> > add_histograms_together
(
		utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> >  histogram_vector,
		utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> >  orbital_histogram

){
	utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> > new_vector_histogram;
	numeric::histograms::TwoDHistogram<core::Size, core::SSize> histogram;

	//put new data
	for(core::Size x=1; x<= histogram_vector.size(); ++x){
		new_vector_histogram.push_back(histogram);
	}

	//iterate through the vector containing all the histograms
	for(core::Size x=1; x<= histogram_vector.size(); ++x){
		//iterate through the elments of the vector to add
		for(core::Size i=0; i <= 100; ++i){
			for(core::SSize j=-10; j<= 10; ++j){
				std::pair<core::Size, core::SSize> pair(i,j);
				core::Size new_counts = histogram_vector[x].lookup_counts(pair) + orbital_histogram[x].lookup_counts(pair);
				new_vector_histogram[x].insert_data(pair, new_counts);
			}
		}
	}
	return new_vector_histogram;
}


int main( int argc, char * argv [] )
{
    try {
	devel::init(argc, argv);
	//using an old job manager because I cant figure out how to use the new job
	//manager with this code. Taken from Ian's old job manager
	std::vector< utility::file::FileName > pdb_file_names;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].active() )
		pdb_file_names = basic::options::option[ basic::options::OptionKeys::in::file::s ]().vector(); // make a copy (-s)
	std::vector< utility::file::FileName > list_file_names;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::l ].active() )
		list_file_names = basic::options::option[ basic::options::OptionKeys::in::file::l ]().vector(); // make a copy (-l)

	for(std::vector< utility::file::FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			pdb_file_names.push_back( utility::file::FileName(line) );
		}
		data.close();
	}
	//end parsing proteins from files


	////////////////////start set up of OrbitalsStatistics  class/////////////////////////////////////////

	//create instance of orbitalsstatics outside of loop function. Only want 1 instance
	core::scoring::orbitals::OrbitalsStatistics orbital_stats;
	//make the histogram the same as the starting histogram when OrbitalsStatistics is initialized
	numeric::histograms::TwoDHistogram<core::Size, core::SSize> histogram(orbital_stats.get_2D_histogram());

	//a vector that will contain all the histograms.
	utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> > histogram_vector;

	//get number of histograms from construction
	core::Size number_of_histograms(orbital_stats.get_number_of_histograms());

	//add histograms to the previously created vector. Notice that the number created is based
	//upon number_of_histograms_ which is constructed in the
	for(core::Size x=1; x<= number_of_histograms; ++x){
		histogram_vector.push_back(histogram);
	}

	////////////////////end set up of OrbitalsStatistics class/////////////////////////////////////////

	//begin iterating through list of files and composing statistics
	for(std::vector< utility::file::FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i) {
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose, i->name(), core::import_pose::PDB_file);

		if(basic::options::option[basic::options::OptionKeys::orbitals::sc_stats]){
			orbital_stats.sc_H_orbital(pose); //where the magic happens

		}
		if(basic::options::option[basic::options::OptionKeys::orbitals::bb_stats]){
			orbital_stats.bb_stats(pose); //where the magic happens

		}


		if(i+1 != i_end){
			histogram_vector=add_histograms_together(histogram_vector, orbital_stats.get_histogram_vector());

		}else{
			histogram_vector=add_histograms_together(histogram_vector, orbital_stats.get_histogram_vector());


			//iterate through the vector containing all the histograms
				for(core::Size x=1; x<= histogram_vector.size(); ++x){
					//open up stream to output files
					std::ofstream combined_total_histogram;
					std::stringstream convert_ss;
					convert_ss << x << ".histogram";
					combined_total_histogram.open(convert_ss.str().c_str());

					//iterate through the elments of the vector to add
					for(core::Size i=0; i <= 100; ++i){
						for(core::SSize j=-10; j<= 10; ++j){
							std::pair<core::Size, core::SSize> pair(i,j);
							combined_total_histogram << i << " " << j << " " << histogram_vector[x].lookup_counts(pair) << std::endl;
						}
					}
				}


		}//end else statement


	}

    }

catch ( utility::excn::EXCN_Base const & e ) {
                          std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                              }
    return 0;


}
