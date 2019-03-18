// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/rdc/RDCEnergy.cc
/// @brief   Implementation of class RDCEnergy
/// @details last Modified: 08/03/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/rdc/RDCEnergy.hh>
#include <protocols/nmr/rdc/RDCEnergyCreator.hh>

// Project headers
#include <core/io/nmr/AtomSelection.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>
#include <core/scoring/nmr/rdc/RDCData.hh>
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/rdc/parameters.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/constants.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace protocols {
namespace nmr {
namespace rdc {

static basic::Tracer TR( "protocols.nmr.rdc.RDCEnergy" );

/// @brief Instantiate a new RDCEnergy
core::scoring::methods::EnergyMethodOP
RDCEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const {
	return core::scoring::methods::EnergyMethodOP( new RDCEnergy );
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
core::scoring::ScoreTypes
RDCEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::nmr_rdc );
	return sts;
}

/// @brief default constructor
RDCEnergy::RDCEnergy() :
	core::scoring::methods::WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new RDCEnergyCreator ) )
{ }

/// @brief copy constructor
RDCEnergy::RDCEnergy(RDCEnergy const & other) :
	core::scoring::methods::WholeStructureEnergy( other ),
	atom_id_to_rdc_xyz_deriv_map_(other.atom_id_to_rdc_xyz_deriv_map_)
{ }

/// @brief destructor
RDCEnergy::~RDCEnergy() { }

core::scoring::methods::EnergyMethodOP
RDCEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new RDCEnergy );
}

/// @brief Return RDCData from pose. Create RDCData if not present and attach them to the pose.
core::scoring::nmr::rdc::RDCData &
RDCEnergy::get_rdc_data_from_pose( core::pose::Pose & pose ) const {
	using namespace core::pose::datacache;
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::rdc;

	if ( pose.data().has( CacheableDataType::NMR_RDC_DATA ) ) {
		return *( utility::pointer::static_pointer_cast< RDCData > ( pose.data().get_ptr( CacheableDataType::NMR_RDC_DATA ) ) );
	} else {
		std::string inputfile;
		if ( basic::options::option[ basic::options::OptionKeys::nmr::rdc::input_file ].active() ) {
			inputfile = basic::options::option[ basic::options::OptionKeys::nmr::rdc::input_file ]();
		} else {
			utility_exit_with_message( "No RDC input file given. You must provide the input file on the command line \"-nmr::rdc::input_file filename\"." );
		}
		RDCDataOP rdc_data_all_ptr(utility::pointer::static_pointer_cast< RDCData >( NMRDataFactory::get_instance()->get_nmr_data("RDC", inputfile, pose) ) );
		pose.data().set( CacheableDataType::NMR_RDC_DATA, rdc_data_all_ptr );
		return *rdc_data_all_ptr;
	}
}

/// @brief Calculate the RDC score and write it into the pose's total energy EnergyMap
void
RDCEnergy::finalize_total_energy(
	Pose & pose,
	ScoreFunction const & /*sxfn*/,
	EnergyMap & totals
) const
{
	totals[ core::scoring::nmr_rdc ] = calcualate_total_score_and_tensors( pose );
}

/// @brief Calculate the total RDC score from RDCData retrieved from the pose
core::Real
RDCEnergy::calcualate_total_score_and_tensors(Pose & pose) const {
	using namespace core::scoring::nmr::rdc;
	RDCData & rdc_data_all = get_rdc_data_from_pose(pose);

	// Some basic setup
	Real total_score(0);
	Size number_media(rdc_data_all.get_number_alignment_media());
	utility::vector1< RDCMultiSetOP > & rdc_multiset_vec = rdc_data_all.get_rdc_multiset_vec();

	for ( Size i = 1; i <= number_media; ++i ) {
		rdc_multiset_vec[i]->update_spin_coordinates( pose );
		Real multiset_score(0);

		// Split behavior depending on if we are going to solve the tensor or
		// calculate the score from fixed tensor values.
		if ( rdc_multiset_vec[i]->tensor_fixed() ) {
			multiset_score = rdc_multiset_vec[i]->compute_rdc_values_and_score_from_tensor();
		} else {
			if ( rdc_multiset_vec[i]->get_computation_type() == RDCMultiSet::SVD ) {
				rdc_multiset_vec[i]->update_matrix_A();
				multiset_score = rdc_multiset_vec[i]->solve_tensor_and_compute_score_by_svd();
			} else {
				multiset_score = rdc_multiset_vec[i]->solve_tensor_and_compute_score_by_nls();
			}
		}
		rdc_multiset_vec[i]->set_atom_derivatives( pose );
		total_score += multiset_score * rdc_multiset_vec[i]->get_weight();

		TR.Trace << "RDC Score for alignment medium " << rdc_multiset_vec[i]->get_alignment_medium_label() << " " << multiset_score
			<< " (unweighted) or " << multiset_score * rdc_multiset_vec[i]->get_weight() << " (weighted)." << std::endl;
		RDCTensorCOP tensor = rdc_multiset_vec[i]->get_tensor_const();
		TR.Trace << "RDC Tensor for alignment medium " << rdc_multiset_vec[i]->get_alignment_medium_label() << ": " << std::endl;
		rdc_multiset_vec[i]->get_computation_type() == RDCMultiSet::SVD ? tensor->show_tensor_stats(TR.Trace, false) : tensor->show_tensor_stats(TR.Trace, true);
	}

	TR.Trace << "Total RDC score " << total_score << std::endl;
	return total_score;
}

/// @brief Called at the beginning of atom tree minimization, this method
///        allows the derived class the opportunity to initialize pertinent data
///        that will be used during minimization.
///        Here, the function creates and updates the atom_id_to_rdc_xyz_deriv_map_ which
///        is needed by the eval_atom_derivative() function.
void
RDCEnergy::setup_for_minimizing(
	Pose & pose,
	ScoreFunction const & /*sxfn*/,
	core::kinematics::MinimizerMapBase const & /*minmap*/
) const
{
	using namespace core::scoring::nmr::rdc;

	RDCData & rdc_data_all = get_rdc_data_from_pose(pose);

	// Set to zero because we are using operator+=
	Vector fij_spinA(0 ,0 ,0);
	Vector fij_spinB(0 ,0 ,0);
	atom_id_to_rdc_xyz_deriv_map_.default_value(fij_spinA);
	atom_id_to_rdc_xyz_deriv_map_.fill_with(fij_spinA);

	utility::vector1< RDCMultiSetOP > & multiset_vec = rdc_data_all.get_rdc_multiset_vec();
	Size number_media(rdc_data_all.get_number_alignment_media());
	utility::vector1<RDCSingle>::const_iterator iter;

	for ( Size i = 1; i <= number_media; ++i ) {
		utility::vector1< RDCSingleSetOP > const & singleset_vec = multiset_vec[i]->get_rdc_singleset_vec();
		Size number_experiments(multiset_vec[i]->get_number_experiments());
		core::conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
		if ( multiset_vec[i]->symmetric_rdc_calc() && core::pose::symmetry::is_symmetric(pose) ) {
			syminfo_ptr = core::pose::symmetry::symmetry_info( pose );
		}

		for ( Size j = 1; j <= number_experiments; ++j ) {
			utility::vector1<RDCSingle> const & single_rdc_vec = singleset_vec[j]->get_single_rdc_vec();

			for ( iter = single_rdc_vec.begin(); iter != single_rdc_vec.end(); ++iter ) {
				utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > const & spinsAB = iter->get_spinsAB();
				utility::vector1< std::pair< Vector, Vector > > const & atom_derivatives = iter->get_atom_derivatives();
				runtime_assert_msg(spinsAB.size() == atom_derivatives.size(), "ERROR in setup of atom tree minimization from RDC derivatives. Vector of AtomIDs and of derivatives have unequal length.");
				for ( Size k = 1; k <= spinsAB.size(); ++k ) {
					fij_spinA = atom_id_to_rdc_xyz_deriv_map_.get(spinsAB[k].first);  // returns (0,0,0) if atom is not present in map, otherwise the corresponding derivative
					fij_spinB = atom_id_to_rdc_xyz_deriv_map_.get(spinsAB[k].second);
					fij_spinA += atom_derivatives[k].first;                    // accumulates the derivative for that particular atom (for all rdc data)
					fij_spinB += atom_derivatives[k].second;
					atom_id_to_rdc_xyz_deriv_map_.set(spinsAB[k].first, fij_spinA);   // sets the accumulated derivative and resizes the map if necessary
					atom_id_to_rdc_xyz_deriv_map_.set(spinsAB[k].second, fij_spinB);
				}
				// Do another pass to set RDC derivatives for symmetric residues
				if ( syminfo_ptr ) {
					for ( Size k = 1; k <= spinsAB.size(); ++k ) {
						utility::vector1< Size >spinA_symm_rsd = syminfo_ptr->bb_clones(spinsAB[k].first.rsd());
						utility::vector1< Size >spinB_symm_rsd = syminfo_ptr->bb_clones(spinsAB[k].second.rsd());
						runtime_assert(spinA_symm_rsd.size() == spinB_symm_rsd.size());

						for ( Size l = 1, l_end = spinA_symm_rsd.size(); l <= l_end; ++l ) {
							core::id::AtomID spinA_symm(spinsAB[k].first.atomno(), spinA_symm_rsd[l]);
							core::id::AtomID spinB_symm(spinsAB[k].second.atomno(), spinB_symm_rsd[l]);

							// return the derivative of that atom, if not present in map, return the default value which is (0,0,0)
							Vector fij_spinA_symm = atom_id_to_rdc_xyz_deriv_map_.get(spinA_symm);
							Vector fij_spinB_symm = atom_id_to_rdc_xyz_deriv_map_.get(spinB_symm);

							// rotate derivative vector into the frame of the symmetric subunit
							// and accumulate the contribution from all rdc experiments
							fij_spinA_symm += core::scoring::nmr::apply_vector_rotation(pose, atom_derivatives[k].first, spinsAB[k].first.rsd(), spinA_symm_rsd[l]);
							fij_spinB_symm += core::scoring::nmr::apply_vector_rotation(pose, atom_derivatives[k].second, spinsAB[k].second.rsd(), spinB_symm_rsd[l]);

							// set the accumulated derivative and resize the map if necessary
							atom_id_to_rdc_xyz_deriv_map_.set(spinA_symm, fij_spinA_symm);
							atom_id_to_rdc_xyz_deriv_map_.set(spinB_symm, fij_spinB_symm);
						}
					}
				}
			}
		}
	}
}

/// @brief Evaluate the xyz derivative of the RDC for an atom in the pose.
///        Called during the atomtree derivative calculation, atom_tree_minimize.cc,
///        through the ScoreFunction::eval_atom_derivative intermediary.
///        F1 and F2 should not be zeroed, rather, setup_for_minimizing() accumulates the RDC
///        contribution from the xyz derivatives of atom id
void
RDCEnergy::eval_atom_derivative(
	core::id::AtomID const & id,
	Pose const & pose,
	core::kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( !atom_id_to_rdc_xyz_deriv_map_.has(id) ) {
		// only checks if id residue number is in range of AtomID_Map
		return;
	}
	Vector fij = atom_id_to_rdc_xyz_deriv_map_.get(id); //  returns (0,0,0) if id is not present in map, otherwise the accumulated derivative of the RDC for this atom

	Vector atom_x = pose.xyz(id);
	Vector const f2( -fij );
	Vector const atom_y = atom_x - f2;   // a "fake" atom in the direction of the gradient
	Vector const f1( atom_x.cross( atom_y ) );

	F1 += weights[ core::scoring::nmr_rdc ] * f1;
	F2 += weights[ core::scoring::nmr_rdc ] * f2;
}


/// @brief Indicate in the context-graphs-required list which
///        context-graphs this energy method requires that the Pose
///        maintain when doing neighbor evaluation. Context graphs are allowed.
void
RDCEnergy::indicate_required_context_graphs(utility::vector1< bool > &) const { }

/// @brief Return the version of the energy method
core::Size
RDCEnergy::version() const { return 1; }

/// @brief register options
void
RDCEnergy::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::rdc::show_info);
}

/// @brief show additional information of the energy method
void
RDCEnergy::show_additional_info(
	std::ostream & tracer,
	Pose & pose,
	bool verbose
) const
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::rdc;

	if ( basic::options::option[ basic::options::OptionKeys::nmr::rdc::show_info ].user() ) {
		verbose = basic::options::option[ basic::options::OptionKeys::nmr::rdc::show_info ]();
	}

	// Some basic setup
	RDCData & rdc_data_all = get_rdc_data_from_pose(pose);
	Real rdc_score = calcualate_total_score_and_tensors(pose);
	Size number_media(rdc_data_all.get_number_alignment_media());
	utility::vector1< Real > scores_all_media(number_media, 0.0); // set all scores to zero because we are going to use +=
	utility::vector1< RDCMultiSetOP > & multiset_vec = rdc_data_all.get_rdc_multiset_vec();

	tracer << " * * * * * * *       RDCEnergy Info       * * * * * * * " << std::endl;
	tracer << "Total Score: " << rdc_score << std::endl;
	tracer << "Is norm of " << number_media << " alignment media." << std::endl;

	// Back-calculate RDCs from determined tensor params
	// Show calc vs exp RDC and tensor params
	// Loop over RDCMultiSets
	Real sum_all_dev_square(0);
	Real sum_all_exp_square(0);
	Size total_number_rdcs(0);
	Size total_number_experiments(0);
	for ( Size i = 1; i <= number_media; ++i ) {
		Real sum_dev_square(0);
		Real sum_exp_square(0);
		Size number_experiments(multiset_vec[i]->get_number_experiments());
		utility::vector1< RDCSingleSetOP > const & singleset_vec = multiset_vec[i]->get_rdc_singleset_vec();

		// Tensor is for all experiments the same
		RDCTensor tensor( *(multiset_vec[i]->get_tensor() ) );
		Real Dmax_ = multiset_vec[i]->get_normalization_type() == NORM_TYPE_CH ? rdc_D_max(RDC_TYPE_CAHA, multiset_vec[i]->correct_sign()) :
			rdc_D_max(RDC_TYPE_NH, multiset_vec[i]->correct_sign()) ;
		tensor.set_Dmax(Dmax_);
		tensor.diagonalize_tensor();
		tensor.reorder_tensor();

		Real Da(tensor.get_Da());
		Real R(tensor.get_R());
		Real alpha(tensor.get_alpha());
		Real beta(tensor.get_beta());
		Real gamma(tensor.get_gamma());

		multiset_vec[i]->update_spin_coordinates(pose);
		utility::vector1< utility::vector1< RDCMultiSet::SpinPairCoordinates > > const & spin_coordinates = multiset_vec[i]->get_spin_coordinates();
		ObjexxFCL::FArray1D<Real> const & rdc_values = multiset_vec[i]->get_rdc_values();
		ObjexxFCL::FArray1D<Real> const & rdc_single_weights = multiset_vec[i]->get_rdc_single_weights();
		Size number_rdcs_per_medium(multiset_vec[i]->get_total_number_rdc());
		Size num_subunits(1);
		if ( multiset_vec[i]->symmetric_rdc_calc() && core::pose::symmetry::is_symmetric(pose) ) {
			core::conformation::symmetry::SymmetryInfoCOP syminfo_ptr = core::pose::symmetry::symmetry_info( pose );
			num_subunits = syminfo_ptr->subunits();
		}

		// In case that we calculate RDCs and automatically take into account symmetry
		// the total length of these vectors is number of RDCs per medium times the number of subunits
		runtime_assert_msg(number_rdcs_per_medium*num_subunits == spin_coordinates.size(), "ERROR in RDCEnergy's show_additional_info() function. Number of RDC and length of spin coordinate vector are not the same");
		runtime_assert_msg(number_rdcs_per_medium*num_subunits == rdc_values.size(), "ERROR in RDCEnergy's show_additional_info() function. Number of RDC and length of RDC value vector are not the same");
		runtime_assert_msg(number_rdcs_per_medium*num_subunits == rdc_single_weights.size(), "ERROR in RDCEnergy's show_additional_info() function. Number of RDC and length of RDC weights vector are not the same");

		Vector euler_angles(alpha, beta, gamma);
		numeric::xyzMatrix< Real > rotM(rotation_matrix_from_euler_angles(euler_angles, tensor.get_euler_convention()));

		RDC_NORM_TYPE norm_type = multiset_vec[i]->get_normalization_type();
		bool correct_sign = multiset_vec[i]->correct_sign();

		// convert Da and R to tensor eigenvalues
		// Da is scaled to NH
		Real A_zz = (4.0 * Da) / (3.0 * Dmax_);
		Real A_xx = Da/Dmax_ * (-2.0/3.0 + R);
		Real A_yy = Da/Dmax_ * (-2.0/3.0 - R);
		Size index_offset(0);
		Real Dmax;

		// Loop over single RDC experiments per alignment medium
		for ( Size j = 1; j <= number_experiments; ++j ) {
			Real score_single_experiment(0);
			Size number_rdcs_per_experiment(singleset_vec[j]->get_number_rdc());
			// RDC prefactor
			Dmax = rdc_D_max(singleset_vec[j]->get_rdc_type(), correct_sign);

			// if RDCs are normalized to N-H use Dmax(NH) otherwise Dmax(CH)
			if ( norm_type == NORM_TYPE_NH || norm_type == NORM_TYPE_NONE ) {
				Dmax /= rdc_scaling_factor_toNH(singleset_vec[j]->get_rdc_type());
			} else if ( norm_type == NORM_TYPE_CH ) {
				Dmax /= rdc_scaling_factor_toCH(singleset_vec[j]->get_rdc_type());
			}
			tracer << "RDC dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
			if ( verbose ) {
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * * * * * * *                exp vs. calc RDC               * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " SpinA_Resid_Atom     SpinB_Resid_Atom     Obs_RDC     Scal_Obs_RDC     Calc_RDC     Scal_Calc_RDC     Deviation     Abs_Deviation " << std::endl;
				tracer << " --------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
			}
			// loop over no subunits
			for ( Size su = 1; su <= num_subunits; ++su ) {

				// Calculate all RDCs for this experiment
				for ( Size k = 1, no_rdcs = number_rdcs_per_experiment; k <= no_rdcs; ++k ) {
					Real calc_rdc(0);
					Size n_eq_spins(spin_coordinates[index_offset + no_rdcs*(su-1) + k].size());

					// loop over equivalent spins
					for ( Size m = 1; m <= n_eq_spins; ++m ) {

						// vector between spins A and B
						Real x(spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].second.x()
							- spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].first.x());
						Real y(spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].second.y()
							- spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].first.y());
						Real z(spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].second.z()
							- spin_coordinates[index_offset + no_rdcs*(su-1) + k][m].first.z());

						if ( norm_type == NORM_TYPE_NH || norm_type == NORM_TYPE_NONE ) {
							// Scale bond vector to length of NH bond (1.041 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.041 / d);
							y *= (1.041 / d);
							z *= (1.041 / d);
						} else if ( norm_type == NORM_TYPE_CH ) {
							// Scale bond vector to length of CAHA bond (1.107 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.107 / d);
							y *= (1.107 / d);
							z *= (1.107 / d);
						}

						// transformed vector after rotation
						Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
						Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
						Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);

						calc_rdc += 1.5 * Dmax * (A_xx * (x_t * x_t) + A_yy * (y_t * y_t) + A_zz * (z_t * z_t));
					}
					if ( multiset_vec[i]->get_averaging_type() == core::scoring::nmr::MEAN ) { calc_rdc /= n_eq_spins; }
					Real rdc_dev(calc_rdc - rdc_values(index_offset + no_rdcs*(su-1) + k));
					sum_dev_square += rdc_dev*rdc_dev;
					sum_exp_square += rdc_values(index_offset + no_rdcs*(su-1) + k)*rdc_values(index_offset + no_rdcs*(su-1) + k);
					Real rdc_abs_dev(std::abs(rdc_dev));
					score_single_experiment += rdc_dev * rdc_dev * rdc_single_weights(index_offset + no_rdcs*(su-1) + k);

					if ( verbose ) {
						std::ostringstream converter;

						for ( auto const & spin_pair : singleset_vec[j]->get_single_rdc_vec()[k].get_spinsAB() ) {
							converter << std::setw(4) << spin_pair.first.rsd() << " " << pose.residue(spin_pair.first.rsd()).atom_name(spin_pair.first.atomno()) << " ";
						}
						tracer << " ( " << std::right << converter.str() << ")";
						converter.clear();
						converter.str("");

						for ( auto const & spin_pair : singleset_vec[j]->get_single_rdc_vec()[k].get_spinsAB() ) {
							converter << std::setw(4) << spin_pair.second.rsd() << " " << pose.residue(spin_pair.second.rsd()).atom_name(spin_pair.second.atomno()) << " ";
						}
						tracer << " ( " << std::right << converter.str() << ")";

						tracer << std::setw(12) << std::fixed << std::setprecision(3);
						if ( norm_type == NORM_TYPE_NH ) {
							tracer << rdc_values(index_offset + no_rdcs*(su-1) + k) * rdc_scaling_factor_toNH(singleset_vec[j]->get_rdc_type());
						} else if ( norm_type == NORM_TYPE_CH ) {
							tracer << rdc_values(index_offset + no_rdcs*(su-1) + k) * rdc_scaling_factor_toCH(singleset_vec[j]->get_rdc_type());
						} else {
							tracer << rdc_values(index_offset + no_rdcs*(su-1) + k);
						}
						tracer << std::setw(17) << std::fixed << std::setprecision(3) << rdc_values(index_offset + no_rdcs*(su-1) + k);
						tracer << std::setw(13) << std::fixed << std::setprecision(3);
						if ( norm_type == NORM_TYPE_NH ) {
							tracer << calc_rdc * rdc_scaling_factor_toNH(singleset_vec[j]->get_rdc_type());
						} else if ( norm_type == NORM_TYPE_CH ) {
							tracer << calc_rdc * rdc_scaling_factor_toCH(singleset_vec[j]->get_rdc_type());
						} else {
							tracer << calc_rdc;
						}
						tracer << std::setw(18) << std::fixed << std::setprecision(3) << calc_rdc;
						tracer << std::setw(14) << std::fixed << std::setprecision(3) << rdc_dev;
						tracer << std::setw(18) << std::fixed << std::setprecision(3) << rdc_abs_dev << std::endl;
					} // number RDCs per experiment
				} // number subunits
			}
			if ( verbose ) {
				tracer << " --------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
			}
			score_single_experiment /= num_subunits; // Since we have summed up over number of subunits we need to divide here again.
			//scores_all_media[i] += std::sqrt(score_single_experiment) * singleset_vec[j]->get_weight();
			scores_all_media[i] += score_single_experiment * singleset_vec[j]->get_weight();
			tracer << "Weighted score for experiment " << singleset_vec[j]->get_dataset_name() << " " << score_single_experiment << std::endl;

			// Increment the index offset
			index_offset += number_rdcs_per_experiment * num_subunits;
		}
		tracer << "Overall weighted score for alignment medium " << multiset_vec[i]->get_alignment_medium_label() << ": " << scores_all_media[i] << std::endl;
		tensor.show_tensor_stats(tracer, true);
		Real Q_fac_per_medium(std::sqrt(sum_dev_square/sum_exp_square));
		Real Rmsd_per_medium(std::sqrt(sum_dev_square/ (number_rdcs_per_medium * num_subunits)));
		tracer << "RDC Q-factor and Rmsd for alignment medium " << multiset_vec[i]->get_alignment_medium_label() << ": " << Q_fac_per_medium << ", " << Rmsd_per_medium << std::endl;
		sum_all_dev_square += sum_dev_square;
		sum_all_exp_square += sum_exp_square;
		total_number_rdcs += number_rdcs_per_medium * num_subunits;
		total_number_experiments += number_experiments;
	}
	Real total_Q_fac(std::sqrt(sum_all_dev_square/sum_all_exp_square));
	Real total_Rmsd(std::sqrt(sum_all_dev_square/total_number_rdcs));
	tracer << "Total RDC Q-factor and Rmsd for " << number_media << " alignment media and " << total_number_experiments
		<< " experiments " << total_Q_fac << ", " << total_Rmsd << std::endl;
}

} // namespace rdc
} // namespace nmr
} // namespace protocols
