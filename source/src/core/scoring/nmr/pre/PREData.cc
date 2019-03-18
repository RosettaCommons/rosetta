// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREData.cc
/// @brief   Implementation of class PREData
/// @details last Modified: 10/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pre/PREData.hh>

// Package headers
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
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
namespace pre {

static basic::Tracer TR( "core.scoring.nmr.pre.PREData" );

/// @brief constructor with filename
PREData::PREData(
	std::string const & filename,
	pose::Pose const & pose
) :
	basic::datacache::CacheableData(),
	number_spinlabel_sites_(0)
{
	register_options();
	init_pre_data_from_file(filename, pose);
}

/// @brief copy constructor
PREData::PREData(PREData const & other) :
	basic::datacache::CacheableData(other),
	number_spinlabel_sites_(other.number_spinlabel_sites_)
{
	// Make deep copy of PREMultiSet vector
	pre_multiset_vec_.clear();
	pre_multiset_vec_.resize(other.pre_multiset_vec_.size());
	for ( Size i = 1; i <= pre_multiset_vec_.size(); ++i ) {
		pre_multiset_vec_[i] = PREMultiSetOP( new PREMultiSet( *(other.pre_multiset_vec_[i]) ) );
	}
}

/// @brief assignment operator
PREData &
PREData::operator=(PREData const & rhs) {
	if ( this != & rhs ) {
		basic::datacache::CacheableData::operator=(rhs);
		// Make deep copy of PREMultiSet vector
		pre_multiset_vec_.clear();
		pre_multiset_vec_.resize(rhs.pre_multiset_vec_.size());
		for ( Size i = 1; i <= pre_multiset_vec_.size(); ++i ) {
			pre_multiset_vec_[i] = PREMultiSetOP( new PREMultiSet( *(rhs.pre_multiset_vec_[i]) ) );
		}
		number_spinlabel_sites_ = rhs.number_spinlabel_sites_;
	}
	return *this;
}

/// @brief destructor
PREData::~PREData() {}

basic::datacache::CacheableDataOP
PREData::clone() const {
	return basic::datacache::CacheableDataOP( new PREData( *this ));
}

/// @brief compute the overall PRE score and individual scores for each spinlabel site
Real
PREData::compute_score_all_spinlabel(
	pose::Pose & pose,
	utility::vector1<Real> & individual_scores
)
{
	using WeightCoordVector = PREMultiSet::WeightCoordVector;

	individual_scores.resize(number_spinlabel_sites_);

	// calculate total score and scores for each spinlabel site
	Real total_score(0);
	for ( Size i = 1; i <= number_spinlabel_sites_; ++i ) {
		WeightCoordVector para_ion_position;
		Real individual_score = pre_multiset_vec_[i]->find_para_ion_position_and_compute_pre_score(pose, para_ion_position);
		pre_multiset_vec_[i]->set_atom_derivatives(pose, para_ion_position);
		total_score += individual_score * pre_multiset_vec_[i]->get_weight();
		individual_scores[i] = individual_score * pre_multiset_vec_[i]->get_weight();
	}
	return total_score;
}

void
PREData::show(std::ostream & tracer) const {
	auto sum_calc = [](Size const & a, PREMultiSetOP const & b) { return a + b->get_total_number_pre(); };
	Size total_pre = std::accumulate(pre_multiset_vec_.begin(), pre_multiset_vec_.end(), 0, sum_calc);
	tracer << "   * * * PREData Summary Report * * *   " << std::endl;
	tracer << "No Spinlabel sites: " << number_spinlabel_sites_ << std::endl;
	tracer << "Total No PREs:      " << total_pre << std::endl;
	for ( Size i = 1; i <= number_spinlabel_sites_; ++i ) {
		pre_multiset_vec_[i]->show(tracer);
	}
}

/// @brief register options
void
PREData::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pre::multiset_weights);
}

/// @brief utility function used during construction of PREData object
void
PREData::init_pre_data_from_file(
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
	utility::vector1< PRESingleSetOP > singleset_vec;
	utility::vector1< PREMultiSetOP > multiset_vec;
	std::map< std::string, std::string > multiset_params_map;
	pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	Size line_number(0);
	bool multisetblock(false);

	// initialize multiset weights from cml options
	utility::vector1<Real> multiset_weights;
	if ( basic::options::option[ basic::options::OptionKeys::nmr::pre::multiset_weights ].user() ) {
		multiset_weights = basic::options::option[ basic::options::OptionKeys::nmr::pre::multiset_weights ]();
	}

	TR.Info << "Opening file '" << utility::file_basename(filename) <<"' " << std::endl;
	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open PRE input file." );
	}

	while ( std::getline(infile, line) ) {
		++line_number;

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;

		} else if ( uppercased(line.substr(0,8)) == "MULTISET" ) {
			// Start of PREMultiSet group
			multisetblock = true;
			continue;

		} else if ( uppercased(line.substr(0,3)) == "END" ) {
			// End of PREMultiSet group
			// Create PREMultiSetOP and set members
			// Clear singleset vector and multiset_params map afterwards
			if ( !singleset_vec.empty() && !multiset_params_map.empty() ) {

				Size spinlabel_position;
				char chain_id;
				std::string ion_type;
				Real temperature;

				// Check that all parameter for the spinlabel and experimental conditions have been provided and convert them to the appropriate type
				// Mandatory are 'spinlabel_position', 'chain_id', 'spinlabel_type' or 'gridsearch', 'ion_type', 'temperature' and 'dataset',
				// and one of either option 'protein_mass' or 'taur'. Optional are 'tauc_min', 'tauc_max' and 'averaging'.
				if ( multiset_params_map.find("spinlabel_position") != multiset_params_map.end() ) {
					spinlabel_position = utility::string2Size(multiset_params_map["spinlabel_position"]);
				} else {
					utility_exit_with_message("ERROR in creating PREData. The spinlabel position for the PREMultiSet has not been provided.");
				}

				if ( multiset_params_map.find("chain_id") != multiset_params_map.end() ) {
					runtime_assert_msg(multiset_params_map["chain_id"].size() == 1, "ERROR in creating PREData. The spinlabel site chain ID for the PREMultiSet must be a single character.");
					chain_id = std::toupper(multiset_params_map["chain_id"].c_str()[0]);
				} else {
					utility_exit_with_message("ERROR in creating PREData. The spinlabel site chain ID for the PREMultiSet has not been provided.");
				}

				if ( multiset_params_map.find("ion_type") != multiset_params_map.end() ) {
					ion_type = multiset_params_map["ion_type"];
				} else {
					utility_exit_with_message("ERROR in creating PREData. The spinlabel paramagnetic ion type for the PREMultiSet has not been provided.");
				}

				if ( multiset_params_map.find("temperature") != multiset_params_map.end() ) {
					temperature = string2Real(multiset_params_map["temperature"]);
				} else {
					utility_exit_with_message("ERROR in creating PREData. The temperature for the PREMultiSet has not been provided.");
				}

				// Try to convert tagging site number from pdb to pose numbering
				if ( pdbinfo ) {
					TR.Debug << "Converting PRE spinlabel site residue " + to_string(spinlabel_position) + " " + to_string(chain_id) + " in PRE input file from pdb to pose numbering." << std::endl;
					spinlabel_position = pdbinfo->pdb2pose(chain_id, spinlabel_position);
					if ( spinlabel_position == 0 ) { // Residue not found
						utility_exit_with_message( "ERROR: Cannot convert PRE spinlabel site residue residue " + to_string(spinlabel_position) + " "
							+ to_string(chain_id) + " in PRE input file from pdb to pose numbering. Residue number was not found." );
					}
				} else { // Assume pose numbering instead
					TR.Warning << "Cannot convert PRE spinlabel site residue in input file from pdb to pose numbering. No PDBInfo object. Using pose numbering instead." << std::endl;
				}

				// Finally create PREMultiSet with spinlabel or gridsearch
				PREMultiSetOP multiset_ptr;

				if ( multiset_params_map.find("spinlabel_type") != multiset_params_map.end() ) {
					if ( multiset_params_map.find("gridsearch") != multiset_params_map.end() ) {
						TR.Warning << "Gridsearch and spinlabel option have been both set in PREMultiSet, but only one method can be used for PRE calculation. ";
						TR.Warning << "Using spinlabel representation and skipping gridsearch." << std::endl;
					}
					std::string spinlabel_type = uppercased(multiset_params_map["spinlabel_type"]);
					NMRSpinlabelOP spinlabel_ptr( new NMRSpinlabel("fa_standard", spinlabel_type) );
					if ( !spinlabel_ptr ) {
						utility_exit_with_message("ERROR in creating PREData. Could not create NMRSpinlabel object. Check if spinlabel ResidueType exists in Rosetta database.");
					}
					multiset_ptr = PREMultiSetOP( new PREMultiSet(singleset_vec, pose, spinlabel_position, ion_type, spinlabel_ptr) );
				} else {
					if ( multiset_params_map.find("gridsearch") != multiset_params_map.end() ) {
						utility::vector1<std::string> gridsearch_params = read_gridsearch_values_from_string(multiset_params_map["gridsearch"]);
						std::string grid_atom1 = uppercased(gridsearch_params[1]);
						std::string grid_atom2 = uppercased(gridsearch_params[2]);
						Real distance_to_atom1 = string2Real(gridsearch_params[3]);
						Real grid_stepsize     = string2Real(gridsearch_params[4]);
						Real grid_min_radius   = string2Real(gridsearch_params[5]);
						Real grid_max_radius   = string2Real(gridsearch_params[6]);

						AtomID atom1(named_atom_id_to_atom_id(NamedAtomID(grid_atom1, spinlabel_position), pose));
						AtomID atom2(named_atom_id_to_atom_id(NamedAtomID(grid_atom2, spinlabel_position), pose));
						NMRGridSearchOP gridsearch_ptr = NMRGridSearchOP( new NMRGridSearch(atom1, atom2, pose, distance_to_atom1, grid_stepsize, grid_min_radius, grid_max_radius) );
						if ( !gridsearch_ptr ) {
							utility_exit_with_message("ERROR in creating PREData. Could not create NMRGridSearch object. Check gridsearch parameter.");
						}
						multiset_ptr = PREMultiSetOP( new PREMultiSet(singleset_vec, pose, spinlabel_position, ion_type, gridsearch_ptr) );
					} else {
						utility_exit_with_message("ERROR in creating PREData. No gridsearch parameter or spinlabel type has been provided for the PREMultiSet.");
					}
				}
				multiset_ptr->set_temperature(temperature);
				if ( multiset_params_map.find("protein_mass") != multiset_params_map.end() ) {
					Real protein_mass = string2Real(multiset_params_map["protein_mass"]);
					multiset_ptr->set_protein_mass(protein_mass);
				} else if ( multiset_params_map.find("taur") != multiset_params_map.end() ) {
					Real taur = string2Real(multiset_params_map["taur"]);
					multiset_ptr->set_tau_r(taur);
				} else {
					utility_exit_with_message("ERROR in creating PREData. Neither the rotational correlation time nor the protein mass for the PREMultiSet have not been provided.");
				}

				if ( multiset_params_map.find("tauc_min") != multiset_params_map.end() &&
						multiset_params_map.find("tauc_max") != multiset_params_map.end() ) {
					// correlation time is provided in nsec in input file
					Real tauc_min = string2Real(multiset_params_map["tauc_min"]) * 1.0e-9;
					Real tauc_max = string2Real(multiset_params_map["tauc_max"]) * 1.0e-9;
					if ( tauc_min > tauc_max ) {
						utility_exit_with_message("ERROR in creating PREData. Value of \"tauc_min\" is greater than \"tauc_max\".");
					}
					multiset_ptr->set_tau_c_limits(tauc_min, tauc_max);
				} else {
					// Raise exit if only one of the two limits for tauc has been provided
					if ( multiset_params_map.find("tauc_min") != multiset_params_map.end() &&
							multiset_params_map.find("tauc_max") == multiset_params_map.end() ) {
						utility_exit_with_message("ERROR in creating PREData. Value of \"tauc_min\" but not of \"tauc_max\" provided for PREMultiSet.");
					} else if ( multiset_params_map.find("tauc_min") == multiset_params_map.end() &&
							multiset_params_map.find("tauc_max") != multiset_params_map.end() ) {
						utility_exit_with_message("ERROR in creating PREData. Value of \"tauc_max\" but not of \"tauc_min\" provided for PREMultiSet.");
					}
				}
				if ( multiset_params_map.find("averaging") != multiset_params_map.end() ) {
					std::string averaging_type = uppercased(multiset_params_map["averaging"]);
					runtime_assert( averaging_type == "SUM"  || averaging_type == "MEAN" );
					multiset_ptr->set_averaging_type(averaging_type);
				}
				multiset_vec.push_back(multiset_ptr);

			} else {
				TR.Warning << "Error in creating PREData. No PRE experiments or spinlabel parameter are available to create PREMultiSet at line "
					<< line_number << " in PRE input file" << std::endl;
				TR.Warning << "Skip line " << line_number << std::endl;
			}
			multisetblock = false;
			singleset_vec.clear();
			multiset_params_map.clear();
			continue;

		} else if ( multisetblock ) {
			// Within the PREMultiSet block
			// Read spinlabel and experimental parameter as key-value pairs
			if ( line.find("spinlabel_position") != std::string::npos ) {
				read_key_value_pair_from_line(line, "spinlabel_position", multiset_params_map, line_number);

			} else if ( line.find("chain_id") != std::string::npos ) {
				read_key_value_pair_from_line(line, "chain_id", multiset_params_map, line_number);

			} else if ( line.find("gridsearch") != std::string::npos ) {
				read_key_value_pair_from_line(line, "gridsearch", multiset_params_map, line_number);

			} else if ( line.find("spinlabel_type") != std::string::npos ) {
				read_key_value_pair_from_line(line, "spinlabel_type", multiset_params_map, line_number);

			} else if ( line.find("ion_type") != std::string::npos ) {
				read_key_value_pair_from_line(line, "ion_type", multiset_params_map, line_number);

			} else if ( line.find("protein_mass") != std::string::npos ) {
				read_key_value_pair_from_line(line, "protein_mass", multiset_params_map, line_number);

			} else if ( line.find("taur") != std::string::npos ) {
				read_key_value_pair_from_line(line, "taur", multiset_params_map, line_number);

			} else if ( line.find("temperature") != std::string::npos ) {
				read_key_value_pair_from_line(line, "temperature", multiset_params_map, line_number);

			} else if ( line.find("tauc_min") != std::string::npos ) {
				read_key_value_pair_from_line(line, "tauc_min", multiset_params_map, line_number);

			} else if ( line.find("tauc_max") != std::string::npos ) {
				read_key_value_pair_from_line(line, "tauc_max", multiset_params_map, line_number);

			} else if ( line.find("averaging") != std::string::npos ) {
				read_key_value_pair_from_line(line, "averaging", multiset_params_map, line_number);

			} else if ( line.find("dataset") != std::string::npos ) {
				Size idx_key = line.find("dataset");
				Size idx_equal_sign = line.find("=", idx_key + 7);
				if ( idx_equal_sign == std::string::npos ) {
					utility_exit_with_message("ERROR: No equal sign found after parameter \"dataset\" in NMR input file. Please provide input parameter \" key = value \".");
				} else {
					utility::vector1<std::string>  dataset_params = read_pre_dataset_params_list(utility::trim(line.substr(idx_equal_sign+1)));
					std::string pre_datafile     = dataset_params[1];
					Real weight                  = string2Real(dataset_params[2]);
					std::string weighting_scheme = uppercased(dataset_params[3]);
					std::string pre_rate_type    = uppercased(dataset_params[4]);
					Real field_strength          = string2Real(dataset_params[5]);

					// Check validity of input
					if ( weighting_scheme != "CONST" && weighting_scheme != "SIGMA" && weighting_scheme != "OBSIG" ) {
						utility_exit_with_message("ERROR in creating PREData. No valid single value weighting type in line " + to_string(line_number)
							+ ". Possible weighting types are: \"CONST\", \"SIGMA\" and \"OBSIG\".");
					}
					if ( pre_rate_type != "R1" && pre_rate_type != "R2" ) {
						utility_exit_with_message("ERROR in creating PREData. No valid PRE rate type in line " + to_string(line_number)
							+ ". Possible weighting types are: \"R1\" and \"R2\".");
					}
					PRESingleSetOP singleset_ptr( new PRESingleSet( pre_datafile, pose, weight, pre_rate_type, weighting_scheme) );
					singleset_ptr->set_field_strength(field_strength);
					singleset_vec.push_back(singleset_ptr);
				}
			}
		}
	}
	if ( multiset_weights.size() != multiset_vec.size() ) {
		TR.Warning << "Number of PREMultiSet weights is inconsistent with the number of PREMultiSets in PRE input file. Set all weights to 1 instead." << std::endl;
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

	// Finally, set the private data member of PREData
	pre_multiset_vec_ = multiset_vec;
	number_spinlabel_sites_ = multiset_vec.size();

	TR.Info << "Finished reading PRE data for " << number_spinlabel_sites_ << " spinlabel site(s) from file " << utility::file_basename(filename) << std::endl;
}

Size
PREData::get_total_number_pre() const {
	auto sum = [](Size const & a, PREMultiSetOP const & b) { return a + b->get_total_number_pre(); };
	return std::accumulate(pre_multiset_vec_.begin(), pre_multiset_vec_.end(), 0, sum);
}

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core
