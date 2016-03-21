// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/jd2/EnsembleJobInputter.cc
/// @brief A Job Inputter for distributing a job based on a set of input structures that make up an ensemble
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/EnsembleJobInputter.hh>
#include <protocols/jd2/EnsembleJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <utility/io/util.hh>
#include <utility/string_util.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.EnsembleJobInputter" );


namespace protocols {
namespace jd2 {

protocols::jd2::EnsembleJobInputter::EnsembleJobInputter(){
	TR << "Instantiate EnsembleJobInputter" << std::endl;
}

protocols::jd2::EnsembleJobInputter::~EnsembleJobInputter(){}


/// @details this function determines what jobs exist from -s/-l
void protocols::jd2::EnsembleJobInputter::fill_jobs( JobsContainer & jobs ){
	
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "EnsembleJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	utility::vector1< std::string > const inputs( basic::options::start_files() );
	core::Size const nstruct( get_nstruct() );
	
	bool seed_weights = option[ basic::options::OptionKeys::jd2::seed_ensemble_weights_file].user() || option[ basic::options::OptionKeys::jd2::seed_ensemble_weights ].user();
	bool seeded_ensemble = option[ basic::options::OptionKeys::jd2::seed_ensemble ]() || seed_weights;
	bool grid_ensemble = option[ basic::options::OptionKeys::jd2::grid_ensemble ]();
	
	if (seeded_ensemble && grid_ensemble){
		utility_exit_with_message("Cannot pass both seed_ensemble and grid_ensemble.  Please choose either or!");
	}
	
	if (nstruct < inputs.size()){
		TR << "Passed nstruct is lower than the number of inputs.  Using ensemble as seeds!" << std::endl;
		seeded_ensemble = true;
		grid_ensemble = false;
	}
	
	
	utility::vector1< core::Size > nstruct_for_inputs;
	
	if (seeded_ensemble){
		nstruct_for_inputs.resize(inputs.size(), 0);
		
		if (seed_weights){
			utility::vector1< core::Real > weights;
			if (option[ basic::options::OptionKeys::jd2::seed_ensemble_weights_file].user() ){
				weights = read_get_weights( option[ basic::options::OptionKeys::jd2::seed_ensemble_weights_file](), inputs );
			}
			else {
				weights = option[ basic::options::OptionKeys::jd2::seed_ensemble_weights]();
			}
			
			if (weights.size() != inputs.size()){
				utility_exit_with_message(" Number of weights must equal the number of inputs! ");
			}
			TR.Debug << utility::to_string(weights) << std::endl;
			numeric::random::WeightedSampler sampler;
			sampler.weights( weights );
			for (core::Size i = 1; i <= nstruct; ++i){
				nstruct_for_inputs[ sampler.random_sample(numeric::random::rg()) ] += 1;
			}
			
		}
		else {
			for ( core::Size i = 1; i <= nstruct; ++i ){
				nstruct_for_inputs[ numeric::random::rg().random_range( 1, inputs.size() )] += 1;
			}
		}
		
	} else if (grid_ensemble){
		core::Size n = nstruct/inputs.size();
		core::Size remainder = nstruct % inputs.size();
		nstruct_for_inputs.resize(inputs.size(), n);
		
		if (remainder != 0){
			nstruct_for_inputs[inputs.size()]+=remainder;
			TR << " Grid search not even! Adding " << remainder << " to nstruct of last ensemble!" << std::endl;
		}
		
	} else {
		utility_exit_with_message(" We should never be here.  Grid ensemble and seeded ensemble are both false!");
	}
	
	core::Size job_index = 1;
	for (core::Size i = 1; i <= inputs.size(); ++i ){
		if (nstruct_for_inputs[i] == 0) continue;
		InnerJobOP ijob( new InnerJob( inputs[i], nstruct_for_inputs[i] ) );
		
		for ( core::Size x = 1; x <= nstruct_for_inputs[i]; ++x ) {
			jobs.push_back( JobOP( new Job( ijob, job_index ) ) );
			TR.Debug << "pushing " << inputs[i] << " nstruct index " << x << std::endl;
			job_index+=1;
			
		}
		TR << "pushed " << inputs[i] << " as nstruct " << nstruct_for_inputs[i] << std::endl;
	}
	
} //fill_jobs

utility::vector1< core::Real >
EnsembleJobInputter::read_get_weights(std::string const & filename, utility::vector1< std::string > const & inputs) const {
	
	std::map< std::string, core::Real > inputs_to_weights;
	inputs_to_weights["ALL"] = 0; //Allows us to specify an ALL, and then only upweight certain models.
	
	TR << "Reading weights: " << filename << std::endl;
	utility::vector1< std::string > const lines( utility::io::get_lines_from_file_data( filename ) );
	if ( lines.size() == 0 ) {
		utility_exit_with_message( "Weight file must contain data!" );
	}
	
	for (core::Size i = 1; i <= lines.size(); ++i){
		std::string fname;
		core::Real weight;
		std::istringstream ss( lines[ i ] );
		
		ss >> fname >> weight;
		inputs_to_weights[ fname ] = weight;
		
	}
	
	utility::vector1< core::Real > weights(inputs.size(), inputs_to_weights["ALL"]); //If weight not given, we defualt to zero.
	
	
	//Match and make the weights.
	
	for (std::map<std::string, core::Real>::const_iterator it = inputs_to_weights.begin(); it != inputs_to_weights.end(); ++it){
	
		for ( core::Size i = 1; i <= inputs.size(); ++i){
			if (inputs[ i ].find(it->first) != std::string::npos) {
    			weights[ i ] = it->second;
				break;
			}
		}
	}
	return weights;
	
}


//CREATOR SECTION
std::string
EnsembleJobInputterCreator::keyname() const
{
	return "EnsembleJobInputter";
}

protocols::jd2::JobInputterOP
EnsembleJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new EnsembleJobInputter );
}



} //protocols
} //jd2






