// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/ParaNMRScoreMover.cc
/// @brief   Implementation file of ParaNMRScoreMover.
/// @details last Modified: 11/15/18
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Project headers
#include <protocols/nmr/ParaNMRScoreMover.hh>
#include <protocols/nmr/ParaNMRScoreMoverCreator.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSData.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/rdc/RDCData.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/scoring/nmr/pre/PREData.hh>
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/rdc/parameters.hh>

// Core headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <ObjexxFCL/format.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// C++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

namespace protocols {
namespace nmr {

static basic::Tracer TR( "protocols.nmr.ParaNMRScoreMover" );

ParaNMRScoreMover::ParaNMRScoreMover() :
	protocols::moves::Mover(),
	verbose_(false),
	write_calc_values_(false),
	write_tensor_info__(false)
{
	sfxn_ = core::scoring::get_score_function();
}

ParaNMRScoreMover::ParaNMRScoreMover( ParaNMRScoreMover const & other ) :
	protocols::moves::Mover( other ),
	sfxn_(other.sfxn_),
	verbose_(other.verbose_),
	write_calc_values_(other.write_calc_values_),
	write_tensor_info__(other.write_tensor_info__)
{}

ParaNMRScoreMover &
ParaNMRScoreMover::operator=( ParaNMRScoreMover const & rhs ) {
	if ( this == &rhs ) {
		return *this;
	}
	return *( new ParaNMRScoreMover( *this ) );
}

ParaNMRScoreMover::~ParaNMRScoreMover() {}

std::string
ParaNMRScoreMoverCreator::keyname() const {
	return ParaNMRScoreMover::mover_name();
}

std::string
ParaNMRScoreMover::mover_name() {
	return "ParaNMRScoreMover";
}

std::string
ParaNMRScoreMover::get_name() const {
	return mover_name();
}

protocols::moves::MoverOP
ParaNMRScoreMoverCreator::create_mover() const {
	return protocols::moves::MoverOP(new ParaNMRScoreMover);
}

protocols::moves::MoverOP
ParaNMRScoreMover::clone() const {
	return protocols::moves::MoverOP(new ParaNMRScoreMover(*this));
}

protocols::moves::MoverOP
ParaNMRScoreMover::fresh_instance() const {
	return protocols::moves::MoverOP(new ParaNMRScoreMover);
}

void
ParaNMRScoreMover::apply( core::pose::Pose & pose ) {
	using namespace core::scoring;
	using namespace core::pose::datacache;
	using namespace core::scoring::nmr;

	// At least one NMR score term must have nonzero weight in the scorefunction.
	if ( sfxn_->has_zero_weight( nmr_pcs ) && sfxn_->has_zero_weight( nmr_rdc ) && sfxn_->has_zero_weight( nmr_pre ) ) {
		utility_exit_with_message("PCS, RDC and PRE score terms have zero weight. At least one NMR score term must have nonzero weight to score the input pose with NMR data.");
	}

	// Score input pose. This also initializes the NMR data and puts them in the pose's datacache.
	(*sfxn_)(pose);

	// Look for the individual NMR data types in the datacache.
	// Calculate scores and Q-factors separated by individual experiments/alignment media/spinlabel sites and add them to the scorefile.
	if ( pose.data().has( CacheableDataType::NMR_PCS_DATA ) ) {
		if ( sfxn_->has_zero_weight( nmr_pcs ) ) {
			TR.Warning << "Weight for 'nmr_pcs' not set in scorefunction. Input pose will not be scored with PCS data." << std::endl;
		} else {
			pcs::PCSDataOP pcs_data = utility::pointer::static_pointer_cast< pcs::PCSData > ( pose.data().get_ptr( CacheableDataType::NMR_PCS_DATA ) );
			add_pcs_scores_to_scorefile(pose, *pcs_data, (*sfxn_)[ nmr_pcs ]);
		}
	}
	if ( pose.data().has( CacheableDataType::NMR_RDC_DATA ) ) {
		if ( sfxn_->has_zero_weight( nmr_rdc ) ) {
			TR.Warning << "Weight for 'nmr_rdc' not set in scorefunction. Input pose will not be scored with RDC data." << std::endl;
		} else {
			rdc::RDCDataOP rdc_data = utility::pointer::static_pointer_cast< rdc::RDCData > ( pose.data().get_ptr( CacheableDataType::NMR_RDC_DATA ) );
			add_rdc_scores_to_scorefile(pose, *rdc_data, (*sfxn_)[ nmr_rdc ]);
		}
	}
	if ( pose.data().has( CacheableDataType::NMR_PRE_DATA ) ) {
		if ( sfxn_->has_zero_weight( nmr_pre ) ) {
			TR.Warning << "Weight for 'nmr_pre' not set in scorefunction. Input pose will not be scored with PRE data." << std::endl;
		} else {
			pre::PREDataOP pre_data = utility::pointer::static_pointer_cast< pre::PREData > ( pose.data().get_ptr( CacheableDataType::NMR_PRE_DATA ) );
			add_pre_scores_to_scorefile(pose, *pre_data, (*sfxn_)[ nmr_pre ]);
		}
	}
}

void
ParaNMRScoreMover::add_pcs_scores_to_scorefile(
	core::pose::Pose & pose,
	core::scoring::nmr::pcs::PCSData & pcs_data,
	core::Real const pcs_wght
) const {
	using namespace core::scoring::nmr::pcs;
	using namespace core::scoring::nmr;
	using namespace ObjexxFCL::format;
	using namespace utility::file;
	using core::scoring::nmr::convert_rdc_type_to_string;
	using utility::trim;
	using utility::strip;

	// Get the current job and its input tag
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string input_tag = job->input_tag();

	TR.Info << "Adding PCS score(s) and Q-factor(s) to scorefile." << std::endl;

	core::Size n_tags = pcs_data.get_number_tags();
	utility::vector1< PCSMultiSetOP > const & multiset_vec = pcs_data.get_pcs_multiset_vec();
	core::Real total_score(0);
	core::Real total_qfac(0);
	core::Real sqdiff_total(0);
	core::Real sqmeas_total(0);

	// Output experimental vs calculated values
	std::ofstream outfile;
	if ( write_calc_values_ ) {
		FileName tag_filename(input_tag);
		std::string filename = tag_filename.base() + std::string("_pcs_pred.txt");
		outfile.open(filename.c_str(), std::ios::out );
		if ( !outfile.is_open() ) {
			TR.Warning << "Unable to open output file " << filename << ". Skip writing of experimental vs. calculated PCS values to file." << std::endl;
		}
	}

	// Output tensor info
	std::fstream tensor_info;
	if ( write_tensor_info__ ) {
		tensor_info.open("pcs_tensor.info", std::ios::in | std::ios::out | std::ios::app );
		if ( !tensor_info.is_open() ) {
			TR.Warning << "Unable to open output file pcs_tensor.info. Skip writing of PCS tensor info to file." << std::endl;
		} else {
			std::string line;
			std::string header = LJ(8, "Position") + RJ(6, "Metal") + RJ(12, "Experiments") + RJ(5, "PCSs") + RJ(10, "Xax") + RJ(10, "Xrh")
				+ RJ(10, "alpha") + RJ(10, "beta") + RJ(10, "gamma") + RJ(10, "xM") + RJ(10, "yM") + RJ(10, "zM") + LJ(12," Description");
			while ( std::getline(tensor_info, line) && line.find(header) == std::string::npos )
					tensor_info.seekg(0, std::ios::end);
			if ( tensor_info.eof() || tensor_info.tellg() == 0 ) {
				tensor_info.clear();
				tensor_info << header << std::endl;
			}
		}
	}

	// Loop over all tagging sites
	for ( core::Size i = 1; i <= n_tags; ++i ) {
		core::Size n_metals = multiset_vec[i]->get_number_metal_ions();
		utility::vector1< PCSSingleSetOP > const & singleset_vec = multiset_vec[i]->get_pcs_singleset_vec();
		core::Real score_per_tag(0);
		core::Real sqdiff_per_tag(0);
		core::Real sqmeas_per_tag(0);

		// Loop over Lanthanides per tagging site
		for ( core::Size j = 1; j <= n_metals; ++j ) {
			core::Size n_pcs = singleset_vec[j]->get_number_pcs();
			utility::vector1< PCSSingle > const & single_pcs_vec = singleset_vec[j]->get_single_pcs_vec();
			core::Real score_per_metal(0);
			core::Real sqdiff(0);
			core::Real sqmeas(0);

			if ( write_calc_values_ && outfile.is_open() ) {
				outfile << "#################################################################" << std::endl;
				outfile << "#Dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
				outfile << "#" << std::right << std::setw(7) << "Residue" << std::setw(5) << "Atom"
					<< std::setw(14) << "Exp_PCS" << std::setw(12) << "Calc_PCS" << std::setw(12)
					<< "Dev_PCS" << "     " << "Dev > Err" << std::endl;
				outfile << "#################################################################" << std::endl;
			}

			// Loop over PCSs for this lanthanide
			for ( core::Size k = 1; k <= n_pcs; ++k ) {
				core::Real diff = single_pcs_vec[k].get_pcs_exp() - single_pcs_vec[k].get_pcs_calc();
				score_per_metal += diff * diff * single_pcs_vec[k].get_weight();
				sqdiff += diff*diff;
				sqmeas += single_pcs_vec[k].get_pcs_exp() * single_pcs_vec[k].get_pcs_exp();

				// Print experimental vs calculated values
				if ( write_calc_values_ && outfile.is_open() ) {
					outfile << " ( " << std::right;
					for ( auto const & id : single_pcs_vec[k].get_protein_spins() ) {
						outfile << std::setw(5) << id.rsd() << " " << std::setw(4) << trim(pose.residue(id.rsd()).atom_name(id.atomno())) << " ";
					}
					outfile << ") ";
					outfile << std::setw(11) << std::fixed << std::setprecision(3) << single_pcs_vec[k].get_pcs_exp() * singleset_vec[j]->get_scaling_factor();
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << single_pcs_vec[k].get_pcs_calc() * singleset_vec[j]->get_scaling_factor();
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << diff * singleset_vec[j]->get_scaling_factor();
					outfile << "      " << std::left << ((diff > single_pcs_vec[k].get_pcs_err()) ? "True  * " : "False   ");
					outfile << std::endl;
				}
			} // PCSs per lanthanide

			// Calculate scores and q-factors and add them to scorefile
			//score_per_tag += std::sqrt(score_per_metal) * singleset_vec[j]->get_weight();
			score_per_tag += score_per_metal * singleset_vec[j]->get_weight();
			core::Real qfac_per_metal = std::sqrt(sqdiff/sqmeas);
			sqdiff_per_tag += sqdiff;
			sqmeas_per_tag += sqmeas;

			if ( verbose_ ) {
				std::stringstream ss;
				ss << "pcs_score_" << multiset_vec[i]->get_tag_residue_number() << "_" << singleset_vec[j]->get_metal_ion_label();
				std::string tag = ss.str();
				// Score reported in scorefile = raw score * singleset weight * multiset weight * nmr_pcs weight of sfxn
				job->add_string_real_pair(tag,score_per_metal * singleset_vec[j]->get_weight() * multiset_vec[i]->get_weight() * pcs_wght);
				ss.str("");
				ss.clear();
				ss << "pcs_qfac_" << multiset_vec[i]->get_tag_residue_number() << "_" << singleset_vec[j]->get_metal_ion_label();
				tag = ss.str();
				job->add_string_real_pair(tag,qfac_per_metal);
			}

			// Print PCS tensor
			TR.Info << "PCS tensor for dataset: " << singleset_vec[j]->get_dataset_name() << ", Tagging site: " << multiset_vec[i]->get_tag_residue_number()
				<< ", Lanthanide: " << singleset_vec[j]->get_metal_ion_label() << std::endl;
			PCSTensorOP tensor(singleset_vec[j]->get_tensor());
			tensor->diagonalize_tensor();
			tensor->reorder_tensor();
			tensor->show_tensor_stats(TR.Info);

			// Output PCS tensor params
			if ( write_tensor_info__ && tensor_info.is_open() ) {
				tensor_info << LJ(8, multiset_vec[i]->get_tag_residue_number());
				tensor_info << std::right;
				tensor_info << A(6, singleset_vec[j]->get_metal_ion_label())
					<< I(12, 1)
					<< I(5, singleset_vec[j]->get_number_pcs())
					<< F(10, 3, tensor->get_ax())
					<< F(10, 3, tensor->get_rh())
					<< F(10, 3, tensor->get_alpha())
					<< F(10, 3, tensor->get_beta())
					<< F(10, 3, tensor->get_gamma())
					<< F(10, 3, tensor->get_metal_center().x())
					<< F(10, 3, tensor->get_metal_center().y())
					<< F(10, 3, tensor->get_metal_center().z());
				tensor_info << std::left << " " << input_tag << std::endl;
			}

		} // Lanthanides per tagging site
		total_score += score_per_tag * multiset_vec[i]->get_weight();
		sqdiff_total += sqdiff_per_tag;
		sqmeas_total += sqmeas_per_tag;

	} // all tagging sites
	total_qfac = std::sqrt(sqdiff_total/sqmeas_total);

	// Calculate total score and total q-factor
	job->add_string_real_pair("pcs_score_total", total_score * pcs_wght);
	job->add_string_real_pair("pcs_qfac_total", total_qfac);

	// Close output files
	if ( outfile.is_open() ) {
		outfile.close();
	}
	if ( tensor_info.is_open() ) {
		tensor_info.close();
	}
}

void
ParaNMRScoreMover::add_rdc_scores_to_scorefile(
	core::pose::Pose & pose,
	core::scoring::nmr::rdc::RDCData & rdc_data,
	core::Real const rdc_wght
) const {
	using namespace core::scoring::nmr::rdc;
	using namespace core::scoring::nmr;
	using namespace ObjexxFCL::format;
	using namespace utility::file;
	using core::scoring::nmr::convert_rdc_type_to_string;
	using utility::trim;
	using utility::strip;

	// Get the current job and its input tag
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string input_tag = job->input_tag();

	TR.Info << "Adding RDC score(s) and Q-factor(s) to scorefile." << std::endl;

	core::Size n_media = rdc_data.get_number_alignment_media();
	utility::vector1< RDCMultiSetOP > const & multiset_vec = rdc_data.get_rdc_multiset_vec();
	core::Real total_score(0);
	core::Real total_qfac(0);
	core::Real sqdiff_total(0);
	core::Real sqmeas_total(0);

	// Output experimental vs calculated values
	std::ofstream outfile;
	if ( write_calc_values_ ) {
		FileName tag_filename(input_tag);
		std::string filename = tag_filename.base() + std::string("_pcs_pred.txt");
		outfile.open(filename.c_str(), std::ios::out );
		if ( !outfile.is_open() ) {
			TR.Warning << "Unable to open output file " << filename << ". Skip writing of experimental vs. calculated RDC values to file." << std::endl;
		}
	}

	// Output tensor info
	std::fstream tensor_info;
	if ( write_tensor_info__ ) {
		tensor_info.open("rdc_tensor.info", std::ios::in | std::ios::out | std::ios::app );
		if ( !tensor_info.is_open() ) {
			TR.Warning << "Unable to open output file rdc_tensor.info. Skip writing of RDC tensor info to file." << std::endl;
		} else {
			std::string line;
			std::string header = LJ(10, "Medium") + RJ(12, "Experiments") + RJ(5, "RDCs") + RJ(10, "Da") + RJ(10, "R")
				+ RJ(12, "Aa") + RJ(12, "Ar") + RJ(10, "alpha") + RJ(10, "beta") + RJ(10, "gamma") + LJ(12," Description");
			while ( std::getline(tensor_info, line) && line.find(header) == std::string::npos )
					tensor_info.seekg(0, std::ios::end);
			if ( tensor_info.eof() ) {
				tensor_info.clear();
				tensor_info << header << std::endl;
			}
		}
	}

	// Loop over all alignment media
	for ( core::Size i = 1; i <= n_media; ++i ) {
		core::Size n_experiments = multiset_vec[i]->get_number_experiments();
		utility::vector1< RDCSingleSetOP > const & singleset_vec = multiset_vec[i]->get_rdc_singleset_vec();
		core::Real score_per_medium(0);
		core::Real sqdiff_per_medium(0);
		core::Real sqmeas_per_medium(0);
		RDC_NORM_TYPE norm_type = multiset_vec[i]->get_normalization_type();

		// Loop of experiments per alignment medium
		for ( core::Size j = 1; j <= n_experiments; ++j ) {
			core::Size n_rdc = singleset_vec[j]->get_number_rdc();
			utility::vector1< RDCSingle > const & single_rdc_vec = singleset_vec[j]->get_single_rdc_vec();
			core::Real score_per_experiment(0);
			core::Real sqdiff(0);
			core::Real sqmeas(0);
			core::Real scaling_factor(1.0);
			if ( norm_type == NORM_TYPE_NH ) {
				scaling_factor = rdc_scaling_factor_toNH(singleset_vec[j]->get_rdc_type());
			} else if ( norm_type == NORM_TYPE_CH ) {
				scaling_factor = rdc_scaling_factor_toCH(singleset_vec[j]->get_rdc_type());
			} else if ( norm_type == NORM_TYPE_NONE ) {
				scaling_factor=1.0;
			}

			if ( write_calc_values_ && outfile.is_open() ) {
				outfile << "#####################################################################################" << std::endl;
				outfile << "#Dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
				outfile << "#" << std::right << std::setw(8) << "ResidueA" << std::setw(6) << "AtomA" << "    " << std::setw(8) << "ResidueB" << std::setw(6) << "AtomB"
					<< std::setw(14) << "Exp_RDC" << std::setw(12) << "Calc_RDC" << std::setw(12) << "Dev_RDC"
					<< "     " << "Dev > Err" << std::endl;
				outfile << "#####################################################################################" << std::endl;
			}

			// Loop over RDCs for this experiment
			for ( core::Size k = 1; k <= n_rdc; ++k ) {
				core::Real diff = single_rdc_vec[k].get_rdc_exp() - single_rdc_vec[k].get_rdc_calc();
				score_per_experiment += diff * diff * single_rdc_vec[k].get_weight();
				sqdiff += diff * diff;
				sqmeas += single_rdc_vec[k].get_rdc_exp() * single_rdc_vec[k].get_rdc_exp();

				// Print experimental vs calculated values
				if ( write_calc_values_ && outfile.is_open() ) {
					outfile << " ( " << std::right;
					for ( auto const & id : single_rdc_vec[k].get_spinsAB() ) {
						outfile << std::setw(6) << id.first.rsd() << " " << std::setw(5) << trim(pose.residue(id.first.rsd()).atom_name(id.first.atomno())) << " ";
					}
					outfile << ")  ( ";
					for ( auto const & id : single_rdc_vec[k].get_spinsAB() ) {
						outfile << std::setw(6) << id.second.rsd() << " " << std::setw(5) << trim(pose.residue(id.second.rsd()).atom_name(id.second.atomno())) << " ";
					}
					outfile << ") ";
					outfile << std::setw(11) << std::fixed << std::setprecision(3) << single_rdc_vec[k].get_rdc_exp() * scaling_factor;
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << single_rdc_vec[k].get_rdc_calc() * scaling_factor;
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << diff * scaling_factor;
					outfile << "      " << std::left << ((diff > single_rdc_vec[k].get_rdc_err()) ? "True  * " : "False   ");
					outfile << std::endl;
				}
			} // RDCs for this experiment

			// Calculate scores and q-factors and add them to scorefile
			//score_per_medium = std::sqrt(score_per_experiment) * singleset_vec[j]->get_weight();
			score_per_medium = score_per_experiment * singleset_vec[j]->get_weight();
			core::Real qfac_per_experiment = std::sqrt(sqdiff/sqmeas);
			sqdiff_per_medium += sqdiff;
			sqmeas_per_medium += sqmeas;

			if ( verbose_ ) {
				std::stringstream ss;
				ss << "rdc_score_" << multiset_vec[i]->get_alignment_medium_label() << "_"
					<< convert_rdc_type_to_string(singleset_vec[j]->get_rdc_type());
				std::string tag = ss.str();
				// Score reported in scorefile = raw score * singleset weight * multiset weight * nmr_rdc weight in sfxn
				job->add_string_real_pair(tag,score_per_experiment * singleset_vec[j]->get_weight() * multiset_vec[i]->get_weight() * rdc_wght);
				ss.str("");
				ss.clear();
				ss << "rdc_qfac_" << multiset_vec[i]->get_alignment_medium_label() << "_"
					<< convert_rdc_type_to_string(singleset_vec[j]->get_rdc_type());
				tag = ss.str();
				job->add_string_real_pair(tag,qfac_per_experiment);
			}
		} // experiments per alignment medium

		total_score += score_per_medium * multiset_vec[i]->get_weight();
		sqdiff_total += sqdiff_per_medium;
		sqmeas_total += sqmeas_per_medium;

		// Print RDC tensor
		TR.Info << "RDC tensor for alignment medium: " << multiset_vec[i]->get_alignment_medium_label() << std::endl;
		RDCTensorOP tensor(utility::pointer::const_pointer_cast< RDCTensor >(multiset_vec[i]->get_tensor()));
		tensor->set_Dmax(rdc_D_max(RDC_TYPE_NH));
		tensor->diagonalize_tensor();
		tensor->reorder_tensor();
		tensor->show_tensor_stats(TR.Info);

		// Output RDC tensor params
		if ( write_tensor_info__ && tensor_info.is_open() ) {
			tensor_info << LJ(10, multiset_vec[i]->get_alignment_medium_label());
			tensor_info << std::right;
			tensor_info << I(12, multiset_vec[i]->get_number_experiments())
				<< I(5, multiset_vec[i]->get_total_number_rdc())
				<< F(10, 3, tensor->get_Da())
				<< F(10, 3, tensor->get_R())
				<< E(12, 3, tensor->get_ax())
				<< E(12, 3, tensor->get_rh())
				<< F(10, 3, tensor->get_alpha())
				<< F(10, 3, tensor->get_beta())
				<< F(10, 3, tensor->get_gamma());
			tensor_info << std::left << " " << input_tag << std::endl;
		}
	} // alignment media
	total_qfac = std::sqrt(sqdiff_total/sqmeas_total);

	// Calculate total score and total q-factor
	job->add_string_real_pair("rdc_score_total", total_score * rdc_wght);
	job->add_string_real_pair("rdc_qfac_total", total_qfac);

	// Close output files
	if ( outfile.is_open() ) {
		outfile.close();
	}
	if ( tensor_info.is_open() ) {
		tensor_info.close();
	}
}

void
ParaNMRScoreMover::add_pre_scores_to_scorefile(
	core::pose::Pose & pose,
	core::scoring::nmr::pre::PREData & pre_data,
	core::Real const pre_wght
) const {
	using namespace core::scoring::nmr::pre;
	using namespace core::scoring::nmr;
	using namespace ObjexxFCL::format;
	using namespace utility::file;
	using utility::trim;
	using utility::strip;

	// Get the current job and its input tag
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string input_tag = job->input_tag();

	TR.Info << "Adding PRE score(s) and Q-factor(s) to scorefile." << std::endl;

	core::Size n_spinlabel = pre_data.get_number_spinlabel_sites();
	utility::vector1< PREMultiSetOP > const & multiset_vec = pre_data.get_pre_multiset_vec();
	core::Real total_score(0);
	core::Real total_qfac(0);
	core::Real sqdiff_total(0);
	core::Real sqmeas_total(0);

	// Output experimental vs calculated values
	std::ofstream outfile;
	if ( write_calc_values_ ) {
		FileName tag_filename(input_tag);
		std::string filename = tag_filename.base() + std::string("_pcs_pred.txt");
		outfile.open(filename.c_str(), std::ios::out );
		if ( !outfile.is_open() ) {
			TR.Warning << "Unable to open output file " << filename << ". Skip writing of experimental vs. calculated PRE values to file." << std::endl;
		}
	}

	// Output spinlabel info
	std::fstream pre_info;
	if ( write_tensor_info__ ) {
		pre_info.open("pre_spinlabel.info", std::ios::in | std::ios::out | std::ios::app );
		if ( !pre_info.is_open() ) {
			TR.Warning << "Unable to open output file pre_spinlabel.info. Skip writing of PRE spinlabel info to file." << std::endl;
		} else {
			std::string line;
			std::string header = LJ(8, "Position") + RJ(12, "Experiments") + RJ(5, "PREs") + RJ(10, "Spinlabel")
				+ RJ(8, "Radical") + RJ(5, "Size") + RJ(12, "tau_r") + RJ(12, "tau_c") + RJ(12, "tau_t") + LJ(12," Description");
			while ( std::getline(pre_info, line) && line.find(header) == std::string::npos )
					pre_info.seekg(0, std::ios::end);
			if ( pre_info.eof() ) {
				pre_info.clear();
				pre_info << header << std::endl;
			}
		}
	}

	// Loop over number of spinlabel sites
	for ( core::Size i = 1; i <= n_spinlabel; ++i ) {
		core::Size n_experiments = multiset_vec[i]->get_number_experiments();
		utility::vector1< PRESingleSetOP > const & singleset_vec = multiset_vec[i]->get_pre_singleset_vec();
		core::Real score_per_sl_site(0);
		core::Real sqdiff_per_sl_site(0);
		core::Real sqmeas_per_sl_site(0);

		// Loop over number of experiments for this spinlabel site
		for ( core::Size j = 1; j <= n_experiments; ++j ) {
			core::Size n_pre = singleset_vec[j]->get_number_pre();
			utility::vector1< PRESingle > const & single_pre_vec = singleset_vec[j]->get_pre_single_vec();
			core::Real score_per_experiment(0);
			core::Real sqdiff(0);
			core::Real sqmeas(0);

			if ( write_calc_values_ ) {
				outfile << "#################################################################" << std::endl;
				outfile << "#Dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
				outfile << "#" << std::right << std::setw(7) << "Residue" << std::setw(5) << "Atom"
					<< std::setw(14) << "Exp_PRE" << std::setw(12) << "Calc_PRE" << std::setw(12)
					<< "Dev_PRE" << "     " << "Dev > Err" << std::endl;
				outfile << "#################################################################" << std::endl;
			}

			// Loop over PREs for this experiment
			for ( core::Size k = 1; k <= n_pre; ++k ) {
				core::Real diff = single_pre_vec[k].get_pre_exp() - single_pre_vec[k].get_pre_calc();
				score_per_experiment += diff * diff * single_pre_vec[k].get_weight();
				sqdiff += diff * diff;
				sqmeas += single_pre_vec[k].get_pre_exp() * single_pre_vec[k].get_pre_exp();

				// Print experimental vs calculated values
				if ( write_calc_values_ ) {
					outfile << " ( " << std::right;
					for ( auto const & id : single_pre_vec[k].get_protein_spins() ) {
						outfile << std::setw(5) << id.rsd() << " " << std::setw(4) << trim(pose.residue(id.rsd()).atom_name(id.atomno())) << " ";
					}
					outfile << ") ";
					outfile << std::setw(11) << std::fixed << std::setprecision(3) << single_pre_vec[k].get_pre_exp() * singleset_vec[j]->get_scaling_factor();
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << single_pre_vec[k].get_pre_calc() * singleset_vec[j]->get_scaling_factor();
					outfile << std::setw(12) << std::fixed << std::setprecision(3) << diff * singleset_vec[j]->get_scaling_factor();
					outfile << "      " << std::left << ((diff > single_pre_vec[k].get_pre_err()) ? "True  * " : "False   ");
					outfile << std::endl;
				}
			} // PREs for this experiment

			// Calculate scores and q-factors and add them to scorefile
			//score_per_sl_site += std::sqrt(score_per_experiment) * singleset_vec[j]->get_weight();
			score_per_sl_site += score_per_experiment * singleset_vec[j]->get_weight();
			core::Real qfac_per_experiment = std::sqrt(sqdiff/sqmeas);
			sqdiff_per_sl_site += sqdiff;
			sqmeas_per_sl_site += sqmeas;

			if ( verbose_ ) {
				std::stringstream ss;
				ss << "pre_score_" << multiset_vec[i]->get_spinlabel_site_rsd() << "_" << singleset_vec[j]->pre_rate_type_to_string()
					<< "_" << convert_pre_spin_type_to_string(single_pre_vec[1].get_pre_spin_type());
				std::string tag = ss.str();
				// Score reported in scorefile = raw score * multiset weight * singleset weight * nmr_pre weight in sfxn
				job->add_string_real_pair(tag,score_per_experiment * singleset_vec[j]->get_weight() * multiset_vec[i]->get_weight() * pre_wght);
				ss.str("");
				ss.clear();
				ss << "pre_qfac_" << multiset_vec[i]->get_spinlabel_site_rsd() << "_" << singleset_vec[j]->pre_rate_type_to_string()
					<< "_" << convert_pre_spin_type_to_string(single_pre_vec[1].get_pre_spin_type());
				tag = ss.str();
				job->add_string_real_pair(tag,qfac_per_experiment);
			}
		} // experiments for this spinlabel site
		total_score += score_per_sl_site * multiset_vec[i]->get_weight();
		sqdiff_total += sqdiff_per_sl_site;
		sqmeas_total += sqmeas_per_sl_site;

		// Print PRESpinlabel and experiment info
		TR.Info << "PRE spinlabel and experiment info for position: " << multiset_vec[i]->get_spinlabel_site_rsd() << std::endl;
		multiset_vec[i]->show(TR.Info);

		// Output PRESpinlabel and experiment info
		if ( write_tensor_info__ && pre_info.is_open() ) {
			pre_info << LJ(8, multiset_vec[i]->get_spinlabel_site_rsd());
			pre_info << std::right;
			pre_info << I(12, multiset_vec[i]->get_number_experiments())
				<< I(5,  multiset_vec[i]->get_total_number_pre());
			if ( multiset_vec[i]->get_spinlabel() ) {
				pre_info << A(10, multiset_vec[i]->get_spinlabel()->get_code())
					<< A(8,  multiset_vec[i]->get_spinlabel()->get_radical_atom())
					<< I(5,  multiset_vec[i]->get_spinlabel()->get_current_ensemble_size());
			} else {
				pre_info << A(10, "None")
					<< A(8,  "None")
					<< I(5,  "None");
			}
			pre_info << E(12, 3, multiset_vec[i]->get_tau_r())
				<< E(12, 3, multiset_vec[i]->get_tau_c())
				<< E(12, 3, multiset_vec[i]->get_tau_t());
			pre_info << std::left << " " << input_tag << std::endl;
		}
	} // all spinlabel sites
	total_qfac = std::sqrt(sqdiff_total/sqmeas_total);

	// Calculate total score and total q-factor
	job->add_string_real_pair("pre_score_total", total_score * pre_wght);
	job->add_string_real_pair("pre_qfac_total", total_qfac);

	// Close output files
	if ( pre_info.is_open() ) {
		pre_info.close();
	}
	if ( outfile.is_open() ) {
		outfile.close();
	}
}

/// @brief Parse tags of XML script
void
ParaNMRScoreMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & //pose//
)
{
	try {
		if ( tag->hasOption("scorefxn") ) {
			sfxn_ = protocols::rosetta_scripts::parse_score_function(tag, "scorefxn", datamap);
		}
		if ( tag->hasOption("verbose") ) {
			verbose_ = tag->getOption< bool >("verbose", false);
		}
		if ( tag->hasOption("output_exp_calc") ) {
			write_calc_values_ = tag->getOption< bool >("output_exp_calc", false);
		}
		if ( tag->hasOption("write_tensor_info_") ) {
			write_tensor_info__ = tag->getOption< bool >("write_tensor_info_", false);
		}
	} catch ( utility::excn::RosettaScriptsOptionError & excn ) {
		TR << "caught exception " << excn.msg() << std::endl;
	}
}

/// @brief Create XML schema definition for ParaNMRScoreMover
void
ParaNMRScoreMover::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd) {
	using namespace utility::tag;

	// Basic attributes
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Scorefunction used by this Mover for scoring with paramagnetic NMR data." )
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "Write separate score and Q-factor values for each spin-label site, alignment medium or metal ion to scorefile.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "output_exp_calc", xsct_rosetta_bool, "Write a table of experimental vs. calculated NMR values for each NMR dataset to prediction file.", "false")
		+ XMLSchemaAttribute::attribute_w_default( "write_tensor_info_", xsct_rosetta_bool, "Write tensor and/or spinlabel info for each NMR dataset to info file.", "false");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"This Mover scores the input pose with PCS, RDC or PRE data, writes the scores and Q-factors to the scorefile, and optionally creates a file of experimental vs. calculated NMR values and a file with tensor infos.", attlist );
}

void
ParaNMRScoreMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ParaNMRScoreMover::provide_xml_schema( xsd );
}

} // nmr
} // protocols
