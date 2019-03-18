// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSData.cc
/// @brief   Implementation of class PCSData
/// @details last Modified: 06/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pcs/PCSData.hh>

// Package headers
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableData.hh>
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
#include <algorithm>
#include <numeric>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "core.scoring.nmr.pcs.PCSData" );

/// @brief construct with filename
PCSData::PCSData(
	std::string const & filename,
	pose::Pose const & pose
) :
	basic::datacache::CacheableData(),
	number_tags_(0)
{
	register_options();
	init_from_cml();
	init_pcs_data_from_file(filename, pose);
}

/// @brief copy constructor
PCSData::PCSData(PCSData const & other) :
	basic::datacache::CacheableData(other),
	number_tags_(other.number_tags_),
	optimize_tensors_(other.optimize_tensors_)
{
	// Make deep copy of PCSMultiSet vector
	pcs_multiset_vec_.clear();
	pcs_multiset_vec_.resize(other.pcs_multiset_vec_.size());
	for ( Size i = 1; i <= pcs_multiset_vec_.size(); ++i ) {
		pcs_multiset_vec_[i] = PCSMultiSetOP( new PCSMultiSet( *(other.pcs_multiset_vec_[i]) ) );
	}
}

/// @brief assignment operator
PCSData &
PCSData::operator=(PCSData const & rhs) {
	if ( this != & rhs ) {
		basic::datacache::CacheableData::operator=(rhs);
		// Make deep copy of PCSMultiSet vector
		pcs_multiset_vec_.clear();
		pcs_multiset_vec_.resize(rhs.pcs_multiset_vec_.size());
		for ( Size i = 1; i <= pcs_multiset_vec_.size(); ++i ) {
			pcs_multiset_vec_[i] = PCSMultiSetOP( new PCSMultiSet( *(rhs.pcs_multiset_vec_[i]) ) );
		}
		number_tags_ = rhs.number_tags_;
		optimize_tensors_ = rhs.optimize_tensors_;
	}
	return *this;
}

/// @brief destructor
PCSData::~PCSData() { }

basic::datacache::CacheableDataOP
PCSData::clone() const {
	return basic::datacache::CacheableDataOP( new PCSData( *this ) );
}

/// @brief register options
void
PCSData::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pcs::multiset_weights);
	option.add_relevant(OptionKeys::nmr::pcs::optimize_tensor);
}

void
PCSData::init_from_cml() {
	using namespace basic::options;
	optimize_tensors_ = option[ basic::options::OptionKeys::nmr::pcs::optimize_tensor ]();
}

/// @brief compute the overall PCS score and individual scores for each tagging site
Real
PCSData::compute_score_all_tags(
	pose::Pose & pose,
	utility::vector1<Real> & scores_all_tags,
	utility::vector1< utility::vector1< PCSTensorCOP > > & tensors_all_lanthanides
)
{
	tensors_all_lanthanides.resize(number_tags_);
	scores_all_tags.resize(number_tags_);

	// calculate total score and scores for each tagging site
	Real total_score(0);
	for ( Size i = 1; i <= number_tags_; ++i ) {
		Real individual_score = pcs_multiset_vec_[i]->compute_score(pose, tensors_all_lanthanides[i]);
		total_score += individual_score * pcs_multiset_vec_[i]->get_weight();
		scores_all_tags[i] = individual_score * pcs_multiset_vec_[i]->get_weight();
	}
	return total_score;
}

void
PCSData::show(std::ostream & tracer) const {
	auto sum_calc = [](Size const & a, PCSSingleSetOP const & b) { return a + b->get_number_pcs(); };
	Size total_pcs(0);
	for ( Size i = 1; i <= number_tags_; ++i ) {
		total_pcs += std::accumulate(std::begin(pcs_multiset_vec_[i]->get_pcs_singleset_vec()),
			std::end(pcs_multiset_vec_[i]->get_pcs_singleset_vec()), 0, sum_calc);
	}
	tracer << "   * * * PCSData Summary Report * * *   " << std::endl;
	tracer << "No Spinlabel Sites: " << number_tags_ << std::endl;
	tracer << "Total No PCSs:      " << total_pcs << std::endl;
	for ( Size i = 1; i <= number_tags_; ++i ) {
		pcs_multiset_vec_[i]->show(tracer);
	}
}

/// @brief utility function used during construction of PCSData object
void
PCSData::init_pcs_data_from_file(
	std::string const & filename,
	pose::Pose const & pose
)
{
	using namespace io::nmr;
	using id::AtomID;
	using id::NamedAtomID;
	using utility::to_string;
	using utility::string2Real;
	using ObjexxFCL::uppercased;

	std::ifstream infile;
	std::string line;
	utility::vector1< PCSSingleSetOP > singleset_vec;
	utility::vector1< PCSMultiSetOP > multiset_vec;
	std::map< std::string, std::string > multiset_params_map;
	pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	Size line_number(0);
	bool multisetblock(false);

	// initialize multiset weights from cml options
	utility::vector1<Real> multiset_weights;
	if ( basic::options::option[ basic::options::OptionKeys::nmr::pcs::multiset_weights ].user() ) {
		multiset_weights = basic::options::option[ basic::options::OptionKeys::nmr::pcs::multiset_weights ]();
	}

	TR.Info << "Opening file '" << utility::file_basename(filename) <<"' " << std::endl;
	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open PCS input file " + filename );
	}

	while ( std::getline(infile, line) ) {
		++line_number;

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;
		} else if ( uppercased(line.substr(0,8)) == "MULTISET" ) {
			// Start of PCSMultiSet group
			multisetblock = true;
			continue;

		} else if ( uppercased(line.substr(0,3)) == "END" ) {
			// End of PCSMultiSet group
			// Create PCSMultiSetOP and set members
			// Clear singleset vector and multiset_params map afterwards
			if ( !singleset_vec.empty() && !multiset_params_map.empty() ) {

				char chain_id;
				Size sl_pos;

				// Check that all parameter for the grid search have been provided and convert them to the appropriate types
				if ( multiset_params_map.find("spinlabel_position") != multiset_params_map.end() ) {
					sl_pos = utility::string2Size(multiset_params_map["spinlabel_position"]);
				} else {
					utility_exit_with_message("ERROR in creating PCSData. The spinlabel site residue number for the PCSMultiSet has not been provided.");
				}

				if ( multiset_params_map.find("chain_id") != multiset_params_map.end() ) {
					runtime_assert_msg(multiset_params_map["chain_id"].size() == 1, "ERROR in creating PCSData. The tagging site chain ID for the PCSMultiSet must be a single character.");
					chain_id = std::toupper(multiset_params_map["chain_id"].c_str()[0]);
				} else {
					utility_exit_with_message("ERROR in creating PCSData. The tagging site chain ID for the PCSMultiSet has not been provided.");
				}

				// Try to convert tagging site number from pdb to pose numbering
				if ( pdbinfo ) {
					TR.Trace << "Converting PCS spinlabel residue " + to_string(sl_pos) + " " + to_string(chain_id) + " in PCS input file from pdb to pose numbering." << std::endl;
					sl_pos = pdbinfo->pdb2pose(chain_id, sl_pos);
					if ( sl_pos == 0 ) { // Residue not found
						utility_exit_with_message( "ERROR: Cannot convert PCS spinlabel residue " + to_string(sl_pos) + " " + to_string(chain_id)
							+ " in PCS input file from pdb to pose numbering. Residue number was not found." );
					}
				} else { // Assume pose numbering instead
					TR.Warning << "Cannot convert PCS spinlabel residue in input file from pdb to pose numbering. No PDBInfo object. Using pose numbering instead." << std::endl;
				}

				// Create PCSMultiSet object and put it into the multiset vector
				PCSMultiSetOP multiset_ptr;

				if ( multiset_params_map.find("spinlabel_type") != multiset_params_map.end() ) {
					if ( multiset_params_map.find("gridsearch") != multiset_params_map.end() ) {
						TR.Warning << "Spinlabel and gridsearch option have been both set in PCSMultiSet, but only one method can be used for PCS calculation. ";
						TR.Warning << "Using spinlabel representation and skipping gridsearch." << std::endl;
					}
					std::string sl_type = uppercased(multiset_params_map["spinlabel_type"]);

					NMRSpinlabelOP spinlabel_ptr( new NMRSpinlabel("fa_standard", sl_type) );
					if ( !spinlabel_ptr ) {
						utility_exit_with_message("ERROR in creating PCSData. Could not create NMRSpinlabel object. Check if spinlabel ResidueType exists in Rosetta database.");
					}
					multiset_ptr = PCSMultiSetOP( new PCSMultiSet(singleset_vec, sl_pos, spinlabel_ptr) );
				} else {
					if ( multiset_params_map.find("gridsearch") != multiset_params_map.end() ) {
						utility::vector1<std::string> gridsearch_params = read_gridsearch_values_from_string(multiset_params_map["gridsearch"]);
						std::string grid_atom1 = uppercased(gridsearch_params[1]);
						std::string grid_atom2 = uppercased(gridsearch_params[2]);
						Real distance_to_atom1 = string2Real(gridsearch_params[3]);
						Real grid_stepsize     = string2Real(gridsearch_params[4]);
						Real grid_min_radius   = string2Real(gridsearch_params[5]);
						Real grid_max_radius   = string2Real(gridsearch_params[6]);

						AtomID atom1(named_atom_id_to_atom_id(NamedAtomID(grid_atom1, sl_pos), pose));
						AtomID atom2(named_atom_id_to_atom_id(NamedAtomID(grid_atom2, sl_pos), pose));
						NMRGridSearchOP gridsearch_ptr = NMRGridSearchOP( new NMRGridSearch(atom1, atom2, pose, distance_to_atom1, grid_stepsize, grid_min_radius, grid_max_radius) );
						if ( !gridsearch_ptr ) {
							utility_exit_with_message("ERROR in creating PCSData. Could not create NMRGridSearch object. Check gridsearch parameter.");
						}
						multiset_ptr = PCSMultiSetOP( new PCSMultiSet(singleset_vec, sl_pos, gridsearch_ptr) );
					} else {
						utility_exit_with_message("ERROR in creating PCSData. No gridsearch parameter or spinlabel type has been provided for the PCSMultiSet.");
					}
				}

				if ( multiset_params_map.find("fixed_tensor") != multiset_params_map.end() ) {
					if ( uppercased(multiset_params_map["fixed_tensor"]) == "TRUE" ) {
						multiset_ptr->fix_tensors();
					}
				}
				multiset_vec.push_back(multiset_ptr);

			} else {
				TR.Warning << "ERROR in creating PCSData. No PCS datasets, gridsearch or spinlabel parameter are available to create PCSMultiSet at line "
					<< line_number << " in PCS input file" << std::endl;
				TR.Warning << "Skip line " << line_number << std::endl;
			}
			multisetblock = false;
			singleset_vec.clear();
			multiset_params_map.clear();
			continue;

		} else if ( multisetblock ) {
			// Within the PCSMultiSet group
			// Read mandatory parameter
			if ( line.find("spinlabel_position") != std::string::npos ) {
				read_key_value_pair_from_line(line, "spinlabel_position", multiset_params_map, line_number);

			} else if ( line.find("chain_id") != std::string::npos ) {
				read_key_value_pair_from_line(line, "chain_id", multiset_params_map, line_number);

			} else if ( line.find("gridsearch") != std::string::npos ) {
				read_key_value_pair_from_line(line, "gridsearch", multiset_params_map, line_number);

			} else if ( line.find("spinlabel_type") != std::string::npos ) {
				read_key_value_pair_from_line(line, "spinlabel_type", multiset_params_map, line_number);

				// fixed_tensor is obligatory and defaults to false
			} else if ( line.find("fixed_tensor") != std::string::npos ) {
				read_key_value_pair_from_line(line, "fixed_tensor", multiset_params_map, line_number);

				// Read parameter related to one dataset and create PCSSingleSet
			} else if ( line.find("dataset") != std::string::npos ) {
				Size idx_key = line.find("dataset");
				Size idx_equal_sign = line.find("=", idx_key + 7);
				if ( idx_equal_sign == std::string::npos ) {
					utility_exit_with_message("ERROR: No equal sign found after parameter \"dataset\" in NMR input file. Please provide input parameter \" key = value \".");
				} else {
					utility::vector1<std::string>  dataset_params = read_pcs_dataset_params_list(utility::trim(line.substr(idx_equal_sign+1)));
					std::string pcs_datafile     = dataset_params[1];
					std::string metal_ion_label  = dataset_params[2];
					Real weight                  = string2Real(dataset_params[3]);
					std::string weighting_scheme = uppercased(dataset_params[4]);
					std::string ave_type         = uppercased(dataset_params[5]);
					std::string computation_type = uppercased(dataset_params[6]);
					utility::vector1<Real> tensor_vals(8);
					tensor_vals[1] = string2Real(dataset_params[12]);  // alpha
					tensor_vals[2] = string2Real(dataset_params[13]);  // beta
					tensor_vals[3] = string2Real(dataset_params[14]);  // gamma
					tensor_vals[4] = string2Real(dataset_params[7]);  // xM
					tensor_vals[5] = string2Real(dataset_params[8]);  // yM
					tensor_vals[6] = string2Real(dataset_params[9]);  // zM
					tensor_vals[7] = string2Real(dataset_params[10]);  // Xax
					tensor_vals[8] = string2Real(dataset_params[11]);  // Xrh

					// Check validity of input
					if ( computation_type != "SVD"   && computation_type != "NLS"   &&
							computation_type != "NLSAX" && computation_type != "NLSRH" && computation_type != "NLSAXRH" ) {
						utility_exit_with_message("ERROR in creating PCSData. No valid computation type in line " + to_string(line_number)
							+ ". Possible computation types are: \"SVD\", \"NLS\", \"NLSAX\", \"NLSRH\" and \"NLSAXRH\".");
					}
					if ( weighting_scheme != "CONST" && weighting_scheme != "SIGMA" && weighting_scheme != "OBSIG" ) {
						utility_exit_with_message("ERROR in creating PCSData. No valid single value weighting type in line " + to_string(line_number)
							+ ". Possible weighting types are: \"CONST\", \"SIGMA\" and \"OBSIG\".");
					}
					if ( ave_type != "SUM" && ave_type != "MEAN" ) {
						utility_exit_with_message("ERROR in creating PCSData. No valid averaging type in in line " + to_string(line_number)
							+ ". Possible averaging types are: \"MEAN\" and \"SUM\".");
					}

					// Create PCSSingleSet from the experiment parameters and put it into the singleset vector
					PCSSingleSetOP singleset_ptr( new PCSSingleSet( pcs_datafile, metal_ion_label, pose, weight, weighting_scheme, computation_type) );
					singleset_ptr->set_averaging_type(ave_type);
					PCSTensorOP tensor( new PCSTensor );
					tensor->set_tensor_in_pas(tensor_vals);
					singleset_ptr->set_tensor(tensor);
					singleset_vec.push_back(singleset_ptr);
				}
			}
		}
	}
	if ( multiset_weights.size() != multiset_vec.size() ) {
		TR.Warning << "Number of PCSMultiSet weights is inconsistent with the number of PCSMultiSets in PCS input file. Set all weights to 1 instead." << std::endl;
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

	// Finally, set the private data member of PCSData
	pcs_multiset_vec_ = multiset_vec;
	number_tags_ = multiset_vec.size();

	TR.Info << "Finished reading PCS data for " << number_tags_ << " tagging site(s) from file " << utility::file_basename(filename) << std::endl;
}

Size
PCSData::get_total_number_pcs() const {
	auto sum = [](Size const & a, PCSMultiSetOP const & b) { return a + b->get_total_number_pcs(); };
	return std::accumulate(pcs_multiset_vec_.begin(), pcs_multiset_vec_.end(), 0, sum);
}

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core
