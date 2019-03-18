// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCData.cc
/// @brief   Implementation of class RDCData
/// @details last Modified: 07/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/rdc/RDCData.hh>

// Package headers
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/io/nmr/AtomSelection.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/rdc/parameters.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// C++ headers
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <numeric>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

static basic::Tracer TR( "core.scoring.nmr.rdc.RDCData" );

/// @brief constructor with filename
RDCData::RDCData(
	std::string const & filename,
	pose::Pose const & pose
) :
	basic::datacache::CacheableData(),
	number_alignment_media_(0)
{
	register_options();
	init_rdc_data_from_file(filename, pose);
}

/// @brief copy constructor
RDCData::RDCData(RDCData const & other) :
	basic::datacache::CacheableData(other),
	number_alignment_media_(other.number_alignment_media_)
{
	// Make deep copy of PCSMultiSet vector
	rdc_multiset_vec_.clear();
	rdc_multiset_vec_.resize(other.rdc_multiset_vec_.size());
	for ( Size i = 1; i <= rdc_multiset_vec_.size(); ++i ) {
		rdc_multiset_vec_[i] = RDCMultiSetOP( new RDCMultiSet( *(other.rdc_multiset_vec_[i]) ) );
	}
}

/// @brief assignment operator
RDCData &
RDCData::operator=(RDCData const & rhs) {
	if ( this != & rhs ) {
		basic::datacache::CacheableData::operator=(rhs);
		// Make deep copy of RDCMultiSet vector
		rdc_multiset_vec_.clear();
		rdc_multiset_vec_.resize(rhs.rdc_multiset_vec_.size());
		for ( Size i = 1; i <= rdc_multiset_vec_.size(); ++i ) {
			rdc_multiset_vec_[i] = RDCMultiSetOP( new RDCMultiSet( *(rhs.rdc_multiset_vec_[i]) ) );
		}
		number_alignment_media_ = rhs.number_alignment_media_;
	}
	return *this;
}

/// @brief destructor
RDCData::~RDCData() { }

basic::datacache::CacheableDataOP
RDCData::clone() const {
	return basic::datacache::CacheableDataOP( new RDCData( *this ) );
}

/// @brief register options
void
RDCData::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::rdc::multiset_weights);
}

/// @brief compute the overall RDC score and scores for the individual alignment media
Real
RDCData::compute_score_all_media(
	pose::Pose const & pose,
	utility::vector1<Real> & scores_all_media,
	utility::vector1< RDCTensorCOP > & tensors_all_media
)
{
	Real total_score(0);
	for ( Size i = 1; i <= number_alignment_media_; ++i ) {
		// register updated coordinates for atoms with an assigned RDC value from the pose
		rdc_multiset_vec_[i]->update_spin_coordinates( pose );

		// calculate total score and scores for each alignment medium
		Real individual_score(0);
		// Split behavior depending on if we are going to solve the tensor or
		// calculate the score from fixed tensor values.
		if ( rdc_multiset_vec_[i]->tensor_fixed() ) {
			individual_score = rdc_multiset_vec_[i]->compute_rdc_values_and_score_from_tensor();
		} else {
			if ( rdc_multiset_vec_[i]->get_computation_type() == RDCMultiSet::SVD ) {
				rdc_multiset_vec_[i]->update_matrix_A();
				individual_score = rdc_multiset_vec_[i]->solve_tensor_and_compute_score_by_svd();
			} else {
				individual_score = rdc_multiset_vec_[i]->solve_tensor_and_compute_score_by_nls();
			}
		}
		rdc_multiset_vec_[i]->set_atom_derivatives( pose );
		tensors_all_media.push_back( rdc_multiset_vec_[i]->get_tensor_const() );
		scores_all_media.push_back( individual_score * rdc_multiset_vec_[i]->get_weight());
		total_score += individual_score * rdc_multiset_vec_[i]->get_weight();
	}
	return total_score;
}

void
RDCData::show(std::ostream & tracer) const {
	auto sum_calc = [](Size const & a, RDCMultiSetOP const & b) { return a + b->get_total_number_rdc(); };
	Size total_rdc = std::accumulate(rdc_multiset_vec_.begin(), rdc_multiset_vec_.end(), 0, sum_calc);
	tracer << "   * * * RDCData Summary Report * * *   " << std::endl;
	tracer << "No Alignment Media: " << number_alignment_media_ << std::endl;
	tracer << "Total No RDCs:      " << total_rdc << std::endl;
	for ( Size i = 1; i <= number_alignment_media_; ++i ) {
		rdc_multiset_vec_[i]->show(tracer);
	}
}

/// @brief utility function used during construction of RDCData object
void
RDCData::init_rdc_data_from_file(
	std::string const & filename,
	pose::Pose const & pose
)
{
	using namespace core::io::nmr;
	using ObjexxFCL::uppercased;
	using utility::string2Real;
	using utility::to_string;

	std::ifstream infile;
	std::string line;
	utility::vector1< RDCMultiSetOP > multiset_vec;
	std::map< std::string, std::string > multiset_params_map;
	utility::vector1< std::string > experiment_filenames, singleset_weighting_schemes;
	utility::vector1< Real > singleset_weights;
	Size line_number(0);
	bool multisetblock(false);

	// initialize multiset weights from cml options
	utility::vector1<Real> multiset_weights;
	if ( basic::options::option[ basic::options::OptionKeys::nmr::rdc::multiset_weights ].user() ) {
		multiset_weights = basic::options::option[ basic::options::OptionKeys::nmr::rdc::multiset_weights ]();
	}

	TR.Info << "Opening file '" << utility::file_basename(filename) <<"' " << std::endl;
	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open RDC input file." );
	}

	while ( std::getline(infile, line) ) {
		++line_number;

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;

		} else if ( uppercased(line.substr(0,8)) == "MULTISET" ) {
			// Start of RDCMultiSet group
			multisetblock = true;
			continue;

		} else if ( uppercased(line.substr(0,3)) == "END" ) {
			// End of RDCMultiSet group
			// Create RDCMultiSetOP and set members
			// clear vector of experiment names, weights, weighting schemes and multiset_params map afterwards
			if ( !experiment_filenames.empty() && !singleset_weighting_schemes.empty()
					&& !singleset_weights.empty() && !multiset_params_map.empty() ) {

				std::string alignment_medium, computation_type;
				utility::vector1<Real> tensor_vals;

				// Check that all parameter for this alignment medium have been provided and convert them to the appropriate type
				// Mandatory are 'alignment_medium', 'computation_type' and 'alignment_tensor'
				// Optional are 'averaging' and 'fixed_tensor'
				if ( multiset_params_map.find("alignment_medium") != multiset_params_map.end() ) {
					alignment_medium = multiset_params_map["alignment_medium"];
				} else {
					utility_exit_with_message("ERROR in creating RDCData. The alignment medium name has not been provided.");
				}

				if ( multiset_params_map.find("computation_type") != multiset_params_map.end() ) {
					computation_type = uppercased(multiset_params_map["computation_type"]);
					if ( computation_type != "SVD"   && computation_type != "NLS"  && computation_type != "NLSDA" &&
							computation_type != "NLSR"  && computation_type != "NLSDAR" ) {
						utility_exit_with_message("ERROR in creating RDCData. Invalid computation type encountered. Possible computation types are \"SVD\", \"NLS\", \"NLSDA\", \"NLSR\" and \"NLSDAR\".");
					}
				} else {
					utility_exit_with_message("ERROR in creating RDCData. The computation type has not been provided.");
				}

				if ( multiset_params_map.find("alignment_tensor") != multiset_params_map.end() ) {
					tensor_vals = read_rdc_tensor_values_from_string(multiset_params_map["alignment_tensor"]);
					if ( convert_string_to_normalization_type(basic::options::option[ basic::options::OptionKeys::nmr::rdc::normalization_type ]()) == NORM_TYPE_CH ) {
						tensor_vals.push_back(rdc_D_max(RDC_TYPE_CAHA, basic::options::option[ basic::options::OptionKeys::nmr::rdc::correct_sign ]()));
					} else {
						tensor_vals.push_back(rdc_D_max(RDC_TYPE_NH, basic::options::option[ basic::options::OptionKeys::nmr::rdc::correct_sign ]()));
					}
					// tensor_vals[1] == alpha
					// tensor_vals[2] == beta
					// tensor_vals[3] == gamma
					// tensor_vals[4] == Da
					// tensor_vals[5] == R
					// tensor_vals[6] == Dmax
				} else {
					utility_exit_with_message("ERROR in creating RDCData. Alignment tensor values have not been provided.");
				}

				// Create RDCMultiSet, set tensor and additional parameter
				RDCMultiSetOP multiset_ptr( new RDCMultiSet(experiment_filenames, alignment_medium, pose) );
				multiset_ptr->set_computation_type(computation_type);
				RDCTensorOP tensor( new RDCTensor );
				tensor->set_tensor_in_pas(tensor_vals);
				multiset_ptr->set_tensor(tensor);

				if ( multiset_params_map.find("averaging") != multiset_params_map.end() ) {
					std::string averaging_type = uppercased(multiset_params_map["averaging"]);
					if ( averaging_type != "SUM" && averaging_type != "MEAN" ) {
						utility_exit_with_message("ERROR in creating RDCData. Invalid averaging type encountered. Possible averaging types are \"SUM\" and \"MEAN\".");
					}
					multiset_ptr->set_averaging_type(averaging_type);
				}

				if ( multiset_params_map.find("fixed_tensor") != multiset_params_map.end() ) {
					if ( uppercased(multiset_params_map["fixed_tensor"]) == "TRUE" ) {
						multiset_ptr->fix_tensor();
					}
				}

				runtime_assert_msg(multiset_ptr->get_number_experiments() == singleset_weights.size(),
					"ERROR in creating RDCData from input file. The number of experiments is inconsistent with the size of the SingleSet weights vector.");
				runtime_assert_msg(multiset_ptr->get_number_experiments() == singleset_weighting_schemes.size(),
					"ERROR in creating RDCData from input file. The number of experiments is inconsistent with the size of the SingleSet weighting scheme vector.");

				for ( Size i = 1; i <= multiset_ptr->get_number_experiments(); ++i ) {
					multiset_ptr->get_rdc_singleset_vec()[i]->set_weight(singleset_weights[i]);
					multiset_ptr->get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme(singleset_weighting_schemes[i]);
				}
				multiset_ptr->update_single_rdc_weighting();
				multiset_vec.push_back(multiset_ptr);

			} else {
				TR.Warning << "ERROR in creating RDCData. No RDC experiments or alignment parameter are available to create RDCMultiSet at line "
					<< line_number << " in RDC input file" << std::endl;
				TR.Warning << "Skip line " << line_number << std::endl;
			}
			multisetblock = false;
			experiment_filenames.clear();
			singleset_weights.clear();
			singleset_weighting_schemes.clear();
			multiset_params_map.clear();
			continue;

		} else if ( multisetblock ) {
			// Within RDCMultiSet group
			// Read alignment tensor parameter as key-value pairs
			// The following parameter are mandatory
			if ( line.find("alignment_medium") != std::string::npos ) {
				read_key_value_pair_from_line(line, "alignment_medium", multiset_params_map, line_number);

			} else if ( line.find("computation_type") != std::string::npos ) {
				read_key_value_pair_from_line(line, "computation_type", multiset_params_map, line_number);

			} else if ( line.find("alignment_tensor") != std::string::npos ) {
				read_key_value_pair_from_line(line, "alignment_tensor", multiset_params_map, line_number);

				// Parameter that the user may provide optionally
			} else if ( line.find("averaging") != std::string::npos ) {
				read_key_value_pair_from_line(line, "averaging", multiset_params_map, line_number);

			} else if ( line.find("fixed_tensor") != std::string::npos ) {
				read_key_value_pair_from_line(line, "fixed_tensor", multiset_params_map, line_number);

				// Read parameter related to one dataset and store them in a vector
			} else if ( line.find("dataset") != std::string::npos ) {
				Size idx_key = line.find("dataset");
				Size idx_equal_sign = line.find("=", idx_key + 7);
				if ( idx_equal_sign == std::string::npos ) {
					utility_exit_with_message("ERROR: No equal sign found after parameter \"dataset\" in NMR input file. Please provide input parameter \" key = value \".");
				} else {
					utility::vector1<std::string>  dataset_params = read_rdc_dataset_params_list(utility::trim(line.substr(idx_equal_sign+1)));
					std::string rdc_datafile     = dataset_params[1];
					Real weight                  = string2Real(dataset_params[2]);
					std::string weighting_scheme = uppercased(dataset_params[3]);

					// Check validity of input
					if ( weighting_scheme != "CONST" && weighting_scheme != "SIGMA" && weighting_scheme != "OBSIG" ) {
						utility_exit_with_message("ERROR in creating RDCData. No valid single value weighting type in line " + to_string(line_number)
							+ ". Possible weighting types are: \"CONST\", \"SIGMA\" and \"OBSIG\".");
					}
					experiment_filenames.push_back(rdc_datafile);
					singleset_weights.push_back(weight);
					singleset_weighting_schemes.push_back(weighting_scheme);
				}
			}
		}
	}
	if ( multiset_weights.size() != multiset_vec.size() ) {
		TR.Warning << "Number of RDCMultiSet weights is inconsistent with the number of RDCMultiSets in RDC input file. Set all weights to 1 instead." << std::endl;
		TR.Warning << "Make sure to provide weights on command line." << std::endl;
		for ( Size i = 1; i <= multiset_vec.size(); ++i ) {
			multiset_vec[i]->set_weight(1.0);
		}
	} else {
		for ( Size i = 1; i <= multiset_vec.size(); ++i ) {
			multiset_vec[i]->set_weight(multiset_weights[i]);
		}
	}
	infile.close();

	// Finally, set the private data member of RDCData
	rdc_multiset_vec_ = multiset_vec;
	number_alignment_media_ = multiset_vec.size();

	TR.Info << "Finished reading RDC data for " << number_alignment_media_ << " alignment media from file " << utility::file_basename(filename) << std::endl;
}

Size
RDCData::get_total_number_rdc() const {
	auto sum = [](Size const & a, RDCMultiSetOP const & b) { return a + b->get_total_number_rdc(); };
	return std::accumulate(rdc_multiset_vec_.begin(), rdc_multiset_vec_.end(), 0, sum);
}

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core
