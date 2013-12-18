// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/PseudocontactShiftEnergy.cc
 ///
 /// @brief
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified June 2009
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergyCreator.hh>

// Package headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftData.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftInput.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftTensor.hh>
#include <protocols/scoring/methods/pcs/GridSearchIterator.hh>
#include <protocols/scoring/methods/pcs/TensorsOptimizer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/PCS.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/scoring/EnergyMap.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Objexx headers
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs {

/// @details This must return a fresh instance of the CarbonHBondEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
PseudocontactShiftEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return new PCS_Energy;
}

core::scoring::ScoreTypes
PseudocontactShiftEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::pcs );
	return sts;
}

basic::Tracer TR_PCS_Energy("protocols.scoring.methods.pcs.PCS_Energy");

void
PCS_Energy::indicate_required_context_graphs( utility::vector1< bool > & ) const{
}

PCS_Energy &
PCS_Energy::operator=(PCS_Energy const & other){
	std::cerr << "Error, == operator not correctly implemented in the PCS_Energy" << std::endl;
	utility_exit_with_message("Exiting");
	if ( this != &other ) {
	}
	return *this;
}

PCS_Energy::PCS_Energy(PCS_Energy const & src ):
	parent( src )
{
}

PCS_Energy::~PCS_Energy(){
}

/// c-tor
PCS_Energy::PCS_Energy() :
	parent( new PseudocontactShiftEnergyCreator )
{}

/// clone
core::scoring::methods::EnergyMethodOP
PCS_Energy::clone() const{
	return new PCS_Energy;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
PCS_Energy::finalize_total_energy(core::pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {
	//using namespace conformation;
	totals[ core::scoring::pcs ] = calculate_pcs_score( pose, false );

} // finalize_total_energy



PCS_data &
PCS_Energy::PCS_data_from_pose(core::pose::Pose & pose) const{

	bool have_exclusions_changed = PCS_Energy_parameters_manager::get_instance()->has_exclude_residues_vector_changed();

	if ( (!have_exclusions_changed) &&
			 ( pose.data().has( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_DATA ) ) ){
		return *( static_cast< PCS_data * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_DATA )() ) );
	}

	bool has_exclude_residues = PCS_Energy_parameters_manager::get_instance()->has_exclude_residues_vector();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1<std::string> vec_filename;
	utility::vector1<core::Real> vec_weight;
	vec_filename = PCS_Energy_parameters_manager::get_instance()->get_vector_filename();
	vec_weight = PCS_Energy_parameters_manager::get_instance()->get_vector_weight();

	TR_PCS_Energy << "Initialization of PCS_data" << std::endl;

	if(vec_filename.size() == 0){
		utility_exit_with_message("Missing input file for PCS. Review your setup file");
	}

	PCS_data_input pcs_d_i = PCS_data_input_manager::get_instance()->get_input_data(vec_filename, vec_weight);

	PCS_dataOP pcs_d;
	if (has_exclude_residues) {
		utility::vector1< bool > exclude_residues;
		exclude_residues = PCS_Energy_parameters_manager::get_instance()->get_vector_exclude_residues();
		pcs_d = new PCS_data(pcs_d_i, exclude_residues);
	} else {
		pcs_d = new PCS_data(pcs_d_i);
	}

	if ( have_exclusions_changed ) {
		PCS_Energy_parameters_manager::get_instance()->exclude_residues_vector_is_current();
	}

	pose.data().set( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_DATA, pcs_d );

	return *pcs_d;
}


void
PCS_Energy::dump_PCS_info(
	utility::vector1<PCS_tensor> const & vec_tensor,
	numeric::xyzVector< core::Real > const & best_coo,
	PCS_data const & pcs_d
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	core::Size i, j;
	static core::Size n_rescore(1);
	utility::vector1<core::Real> A(5, 0);

	if( option[ basic::options::OptionKeys::PCS::write_extra ].user() ){

		std::string file_dump (option[ basic::options::OptionKeys::PCS::write_extra ]());

		std::ofstream myfile;
		if(n_rescore == 1){
			myfile.open (file_dump.c_str(), std::ios::out);
			myfile << "# Tensor: Xxx Xxy Xxz Xyy Xyz x y z" << std::endl;
			myfile << "# Spins: res_num atom_name PCS_exp PCS_calc PCS_dev PCS_abs_dev" << std::endl;
		}
		else{
			myfile.open (file_dump.c_str(), std::ios::app);
		}
		if (!myfile.is_open ()){
			std::cerr << "Unable to open the file '" << file_dump  <<"'" << std::endl;
			utility_exit();
		}


		const utility::vector1<core::Real> & X_all(pcs_d.get_X_all());
		const utility::vector1<core::Real> & Y_all(pcs_d.get_Y_all());
		const utility::vector1<core::Real> & Z_all(pcs_d.get_Z_all());

		for (i = 1 ; i <= pcs_d.get_n_lanthanides(); ++i){
			PCS_data_per_lanthanides const & PCS_d_p_l (pcs_d.get_pcs_data_per_lanthanides_all()[i]);
			utility::vector1<PCS_line_data> const & PCS_d_l_a_s (pcs_d.get_PCS_data_line_all_spin());

			myfile << "#" << PCS_d_p_l.get_filename() << " RESCORE NUMBER " << n_rescore << std::endl;

			core::Real Xxx(vec_tensor[i].chi_xx());
			core::Real Xxy(vec_tensor[i].chi_xy());
			core::Real Xxz(vec_tensor[i].chi_xz());
			core::Real Xyy(vec_tensor[i].chi_yy());
			core::Real Xyz(vec_tensor[i].chi_yz());
			const utility::vector1<core::Size> & A_index( PCS_d_p_l.get_A_index());
			const ObjexxFCL::FArray1D< core::Real > & fstyle_b(PCS_d_p_l.get_fstyle_b());

			myfile << "# Tensor: " << std::setw(10) << Xxx << " " << Xxy << " "  << Xxz << " " << Xyy << " " << Xyz << " " << best_coo.x() << " " << best_coo.y() << " " << best_coo.z() << std::endl;

			for ( j = 1; j <= PCS_d_p_l.get_n_pcs(); ++j){
				core::Real PCS_exp (fstyle_b(j));
				core::Size idx (A_index[j]);
				core::Real x (X_all[idx]);
				core::Real y (Y_all[idx]);
				core::Real z (Z_all[idx]);
				core::Size res_num(PCS_d_l_a_s[idx].residue_num());
				std::string atom_name(PCS_d_l_a_s[idx].atom_name());
				fill_A_line(A, best_coo.x(), best_coo.y(), best_coo.z(), x, y, z);
				core::Real PCS_calc(A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz);
				core::Real PCS_dev (PCS_exp - PCS_calc);
				core::Real PCS_abs_dev (fabs(PCS_exp - PCS_calc));
				myfile << res_num <<" " << atom_name << std::setw(10) << PCS_exp << "  " << PCS_calc<< "  " << PCS_dev<< "  " << PCS_abs_dev<< "  " << std::endl;
			}
		}
		n_rescore ++;
		myfile.close();
	}
}


core::Real
PCS_Energy::calculate_pcs_score(core::pose::Pose & pdb, bool print_to_tracer) const{

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


  utility::vector1<PCS_tensor> vec_tensor;
	utility::vector1<core::Real> vec_score;

	core::Real pcs_score_total;
	core::Size i;


	core::Real pcs_weight (PCS_Energy_parameters_manager::get_instance()->get_pcs_weight());

	PCS_data &pcs_d = PCS_data_from_pose(pdb);

	if (pcs_weight == 0){
		return 0;
	}


	//alloc best score and vector and coordinate
	for ( i = 1; i <= pcs_d.get_n_lanthanides(); ++i){
		PCS_tensor PCS_t = PCS_tensor(0, 0, 0, 0, 0, ((pcs_d.get_pcs_data_per_lanthanides_all())[i]).get_filename());
		vec_tensor.push_back(PCS_t);
	}
	vec_score.resize(pcs_d.get_n_lanthanides());

	numeric::xyzVector< core::Real > best_coo;

	//call to calculate the tensors and the score
	pcs_score_total = calculate_scores_and_tensors_from_pose_and_PCS_data(vec_score, vec_tensor, best_coo, pdb, pcs_d);

	bool minimize_best_tensor;
	minimize_best_tensor = PCS_Energy_parameters_manager::get_instance()->get_minimize_best_tensor();

	print_to_tracer = false;

	if(minimize_best_tensor){
		if(print_to_tracer){//Only called in PCS_main at the moment, quick flag.

			TR_PCS_Energy << "*** Before minimization of the tensor ***" << std::endl;
			TR_PCS_Energy << "Score: " << pcs_score_total << std::endl;
			TR_PCS_Energy << "Sum of: ";
			for(i = 1; i <= vec_score.size(); ++i){
				TR_PCS_Energy	<< vec_score[i] << " ";
			}
			TR_PCS_Energy	<< std::endl;
			TR_PCS_Energy << "Score weighted: " << pcs_weight * pcs_score_total << std::endl;
			TR_PCS_Energy << "Tensors found:" << std::endl;
			for(i = 1; i <= vec_tensor.size(); ++i){
				TR_PCS_Energy	<< vec_tensor[i] << std::endl;
			}
			TR_PCS_Energy << "Lanthanide position: " << best_coo.x() << " " << best_coo.y() << " " << best_coo.z() << std::endl;
			}

		core::Real optimized_score (minimize_tensors_from_PCS_data(vec_tensor, best_coo, pcs_d));
		//		std::cerr << pcs_score_total << " -> " << optimized_score << std::endl;
		core::Real tolerance( 0.001);
		if((pcs_score_total + tolerance) < optimized_score){
			TR_PCS_Energy << "Warning, optimized score has a higher value than starting position. Problem with minimizer?" << std::endl;
			TR_PCS_Energy << pcs_score_total << " -> " << optimized_score << std::endl;
			//utility_exit_with_message("PROBLEM WITH MINIMIZER");
		}

		pcs_score_total = optimized_score;

		if(print_to_tracer){//Only called in PCS_main at the moment, quick flag.
			TR_PCS_Energy << "*** After minimization of the tensor ***" << std::endl;
			TR_PCS_Energy << "Score: " << optimized_score << std::endl;
			TR_PCS_Energy << "Sum of: NOT AVAILABLE";
			TR_PCS_Energy	<< std::endl;
			TR_PCS_Energy << "Score weighted: " << pcs_weight * optimized_score << std::endl;
			TR_PCS_Energy << "Tensors found:" << std::endl;
			for(i = 1; i <= vec_tensor.size(); ++i){
				TR_PCS_Energy	<< vec_tensor[i] << std::endl;
			}
			TR_PCS_Energy << "Lanthanide position: " << best_coo.x() << " " << best_coo.y() << " " << best_coo.z() << std::endl;
			}
	}


	dump_PCS_info(vec_tensor, best_coo, pcs_d);

	return (pcs_weight * pcs_score_total);
}

core::Real
PCS_Energy::minimize_tensors_from_PCS_data(	utility::vector1<PCS_tensor> & vec_best_tensor,
																						numeric::xyzVector< core::Real > & best_coo,
																						PCS_data const & pcs_d
																						) const{

	core::Size i;

	utility::vector1<core::Real> vect_to_opt;
	vect_to_opt.push_back(best_coo.x());
	vect_to_opt.push_back(best_coo.y());
	vect_to_opt.push_back(best_coo.z());

	if(vec_best_tensor.size() !=  pcs_d.get_n_lanthanides()){
		utility_exit_with_message("n_lanthanides and vec_best_tensor size differs in minimize_tensors_from_PCS_data");
	}

	for(i = 1; i <= vec_best_tensor.size(); ++i){
		vect_to_opt.push_back(vec_best_tensor[i].chi_xx());
		vect_to_opt.push_back(vec_best_tensor[i].chi_xy());
		vect_to_opt.push_back(vec_best_tensor[i].chi_xz());
		vect_to_opt.push_back(vec_best_tensor[i].chi_yy());
		vect_to_opt.push_back(vec_best_tensor[i].chi_yz());
	}

	TensorsOptimizer tensors_opt(pcs_d);
	//	optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone_atol", 0.0000001, true, false, false );
	core::optimization::MinimizerOptions options( "dfpmin", 0.00001, true, false, false );
	core::optimization::Minimizer minimizer(tensors_opt, options );

	core::Real optimized_cost(minimizer.run( vect_to_opt ));

	best_coo.assign(vect_to_opt[1], vect_to_opt[2], vect_to_opt[3]);

	for(i = 1; i <= vec_best_tensor.size(); ++i){
		vec_best_tensor[i].reset_tensor((core::Real)vect_to_opt[3 + 5*(i-1) + 1],
																		(core::Real)vect_to_opt[3 + 5*(i-1) + 2],
																		(core::Real)vect_to_opt[3 + 5*(i-1) + 3],
																		(core::Real)vect_to_opt[3 + 5*(i-1) + 4],
																		(core::Real)vect_to_opt[3 + 5*(i-1) + 5]);
	}
	return (optimized_cost);
}


//This will be called for each new pose
core::Real
PCS_Energy::calculate_scores_and_tensors_from_pose_and_PCS_data(	utility::vector1<core::Real> & vec_best_score,
																																	utility::vector1<PCS_tensor> & vec_best_tensor,
																																	numeric::xyzVector< core::Real > & best_coo,
																																	core::pose::Pose const & pdb,
																																	PCS_data & pcs_d) const{

	//	using namespace basic::options;
	//	using namespace basic::options::OptionKeys;
	core::Real x, y, z;
	core::Real best_score, score; //, score2;
	core::Size i;
	core::Size size_of;
	utility::vector1<core::Real> vec_score_temp;
	utility::vector1<PCS_tensor> vec_tensor_temp;

	//ref to do atomic best switch in order to avoid copy stuff
	utility::vector1<core::Real> * vec_score_ref_current;
	utility::vector1<PCS_tensor> * vec_tensor_ref_current;

	//PCS_Energy_parameters_manager::get_instance()->print_grid_param();

	//some basic checking...
	size_of = vec_best_score.size();
	if((size_of != vec_best_tensor.size())||
		 (size_of != pcs_d.get_n_lanthanides())){
		std::cerr << "Problem in calculate_scores_and_tensors_from_pose_and_PCS_data function" << std::endl;
		std::cerr << "n_lanthanides =  " << pcs_d.get_n_lanthanides();
		std::cerr << "vec_best_tensor.size() = " << vec_best_tensor.size();
		std::cerr << "vec_best_score.size() = " << vec_best_score.size() << std::endl;
		utility_exit_with_message("Exiting");
	}

	for( i = 1; i <= size_of; i++){
		vec_tensor_temp.push_back(vec_best_tensor[i]);
	}

	vec_score_temp.resize(size_of);

	//Must be after the resize statment!
	vec_score_ref_current = & vec_score_temp;
	vec_tensor_ref_current = & vec_tensor_temp;

	pcs_d.update_X_Y_Z_all(pdb);

	//	TR_PCS_Energy << "Reading grid search paramaters" << std::endl;

	core::Real grid_edge (PCS_Energy_parameters_manager::get_instance()->get_grid_edge());
	core::Real grid_step (PCS_Energy_parameters_manager::get_instance()->get_grid_step());
	core::Real grid_small_cutoff (PCS_Energy_parameters_manager::get_instance()->get_grid_small_cutoff());
	core::Real grid_large_cutoff (PCS_Energy_parameters_manager::get_instance()->get_grid_large_cutoff());
	core::Real grid_cone_angle_cutoff (PCS_Energy_parameters_manager::get_instance()->get_grid_cone_angle_cutoff());
	std::string grid_atom_name_1 (PCS_Energy_parameters_manager::get_instance()->get_grid_atom_name_1());
	std::string grid_atom_name_2 (PCS_Energy_parameters_manager::get_instance()->get_grid_atom_name_2());
	core::Size grid_residue_num_1 (PCS_Energy_parameters_manager::get_instance()->get_grid_residue_num_1());
	core::Size grid_residue_num_2 (PCS_Energy_parameters_manager::get_instance()->get_grid_residue_num_2());
	core::Real grid_k_vector (PCS_Energy_parameters_manager::get_instance()->get_grid_k_vector());

	if(grid_residue_num_1 > pdb.total_residue()){
		std::cerr << "Error: Couldn't find residue " << grid_residue_num_1 << std::endl;
		std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
		utility_exit_with_message("Can't define gridsearchiterator");
	}
	if( grid_residue_num_2> pdb.total_residue()){
		std::cerr << "Error: Couldn't find residue " << grid_residue_num_1 << std::endl;
		std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
		utility_exit_with_message("Can't define gridsearchiterator");
	}
	if( ! pdb.residue(grid_residue_num_1).has(grid_atom_name_1)){
		std::cerr << "Error: Couldn't find the atom " << grid_atom_name_1 << " in residue " << grid_residue_num_1 << std::endl;
		std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
		utility_exit_with_message("Can't define gridsearchiterator");
	}
	if( ! pdb.residue(grid_residue_num_2).has(grid_atom_name_2)){
		std::cerr << "Error: Couldn't find the atom " << grid_atom_name_2 << " in residue " <<  grid_residue_num_2<< std::endl;
		std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
		utility_exit_with_message("Can't define gridsearchiterator");
	}

	numeric::xyzVector< core::Real > coo1 = pdb.residue(grid_residue_num_1).atom(grid_atom_name_1).xyz();
	numeric::xyzVector< core::Real > coo2 = pdb.residue(grid_residue_num_2).atom(grid_atom_name_2).xyz();

	GridSearchIterator grid_it(coo1, coo2, grid_k_vector, grid_edge, grid_step, grid_small_cutoff, grid_large_cutoff, grid_cone_angle_cutoff);

	best_coo.assign(coo1.x(), coo1.y(), coo1.z());
	best_score = 999999999999999999999999999.9; //std::numeric_limits::infinity();x
	bool test_at_least_one_iteration = false;

	while(grid_it.next_center(x, y, z) == true){
		test_at_least_one_iteration = true;
		//TR_PCS_Energy << "trying x= " << x << "y= " << y << "z= " << z << std::endl;
		//		std::cout  << x << " " << y << " " << z << " SCANNER " << std::endl;
		pcs_d.update_matrix_A_all(x, y, z);

		score = 0;
		for(i = 1; i <= pcs_d.get_n_lanthanides(); ++i){
			(*vec_score_ref_current)[i] = pcs_d.get_pcs_data_per_lanthanides_all()[i].calculate_tensor_and_cost_with_svd((*vec_tensor_ref_current)[i]);
			//score += (*vec_score_ref_current)[i];
			score += (*vec_score_ref_current)[i] * (*vec_score_ref_current)[i];

			if (score > best_score){ // if a single lanthanide already give a worse score, no need to look for other lanthanides
				continue;
			}
		}

		if ( score < best_score){
			best_score = score;

			best_coo.assign(x, y, z);

			//atomic switch
			if((vec_score_ref_current != &vec_score_temp) && (vec_score_ref_current != &vec_best_score)){ //test to make sure...
				std::cerr << "Problem in calculate_scores_and_tensors_from_pose_and_PCS_data function" << std::endl;
				std::cerr << "The atomic switch is not working (1)" << std::endl;
				utility_exit_with_message("Exiting");
			}
			if(vec_tensor_ref_current == &vec_best_tensor){
				vec_tensor_ref_current = &vec_tensor_temp;
				vec_score_ref_current = &vec_score_temp;
			}
			else{
				if(vec_tensor_ref_current != &vec_tensor_temp){ //test to make sure...
					std::cerr << "Problem in calculate_scores_and_tensors_from_pose_and_PCS_data function" << std::endl;
					std::cerr << "The atomic switch is not working (2)" << std::endl;
					utility_exit_with_message("Exiting");
				}
				vec_tensor_ref_current = &vec_best_tensor;
				vec_score_ref_current = &vec_best_score;
			}
			//end atomic switch
		}
	} //while


	if(test_at_least_one_iteration == false){
		std::cerr << "The description of the grid search given is too restrictive" << std::endl;
		utility_exit_with_message("Exiting");
	}

	if( vec_tensor_ref_current == &vec_best_tensor ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_tensor[i].copy_from_ref(vec_tensor_temp[i]);
		}
	}

	if( vec_score_ref_current == &vec_best_score ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_score[i] = vec_score_temp[i];
		}
	}
	return (sqrt(best_score));
}

core::Size
PCS_Energy::version() const
{
	return 1; // Initial versioning
}

PCS_Energy_parameters_manager::PCS_Energy_parameters_manager(){
	/*
//Do I need to initialize to some values?? In principle no.
	grid_edge_ = 20;
	grid_step_ = 2;
	grid_small_cutoff_ = 4;
	grid_large_cutoff_ = 	10;
	grid_cone_angle_cutoff_ = 180;
	grid_atom_name_1_ = "CA";
	grid_atom_name_2_ = "CB";
	grid_residue_num_1_ = 68;
	grid_residue_num_2_ = 68;
	grid_k_vector_ = 1;
	minimize_best_tensor_ = true;
	pcs_weight_ = 100;
	*/

	utility::vector1< bool > tmp(false, 0);
	vec_exclude_residues_ = tmp;
	vec_exclude_residues_exists_ = false;
	vec_exclude_residues_changed_ = false;
}


//rvernon -> partial PCS score machinery in development

// The input exclude residues vector is a list of residue numbers
// This function converts that into a vector of bools, where the index numbers are residues,
// so excluded residues return "true" and non excluded residues return "false"
//
// This is sized up to the largest excluded residue, so when doing the exclusion check you also have
// to say "excluded = false" for residues outside of the array.
void
PCS_Energy_parameters_manager::set_vector_exclude_residues(utility::vector1< core::Size > const vec_exclude) {

	if (vec_exclude.size() > 0) {

		core::Size largest_n = 0;

		for (core::Size i = 1; i <= vec_exclude.size(); ++i) {
			if ( vec_exclude[i] > largest_n ) {
				largest_n = vec_exclude[i];
			}
		}

		utility::vector1< bool > temp (largest_n, false);

		for (core::Size i = 1; i <= vec_exclude.size(); ++i) {
			temp[vec_exclude[i]] = true;
		}

		if ( vec_exclude_residues_ == temp ) {
			vec_exclude_residues_changed_ = false;
		} else {
			vec_exclude_residues_changed_ = true;
		}

		vec_exclude_residues_ = temp;
		vec_exclude_residues_exists_ = true;
	}
}

void
PCS_Energy_parameters_manager::remove_vector_exclude_residues() {
	utility::vector1< bool > temp (false, 0);
	vec_exclude_residues_ = temp;
	vec_exclude_residues_exists_ = false;
	vec_exclude_residues_changed_ = true;
}

bool
PCS_Energy_parameters_manager::has_exclude_residues_vector() {
	return vec_exclude_residues_exists_;
}

bool
PCS_Energy_parameters_manager::has_exclude_residues_vector_changed() {
	return vec_exclude_residues_changed_;
}


utility::vector1< bool >
PCS_Energy_parameters_manager::get_vector_exclude_residues() {
	return vec_exclude_residues_;
}

void
PCS_Energy_parameters_manager::exclude_residues_vector_is_current() {
	vec_exclude_residues_changed_ = false;
}


//rvernon

void
PCS_Energy_parameters_manager::set_vector_name_and_weight(utility::vector1<std::string> const vec_filename,
																													utility::vector1<core::Real> const vec_individual_weight){

	vec_filename_ = vec_filename;
	vec_individual_weight_ = vec_individual_weight;
}

void
PCS_Energy_parameters_manager::set_grid_param(core::Real const grid_edge,
																							core::Real const grid_step,
																							core::Real const grid_small_cutoff,
																							core::Real const grid_large_cutoff,
																							core::Real const grid_cone_angle_cutoff,
																							std::string const grid_atom_name_1,
																							std::string const grid_atom_name_2,
																							core::SSize const grid_residue_num_1,
																							core::SSize const grid_residue_num_2,
																							core::Real const grid_k_vector,
																							bool const minimize_best_tensor,
																							core::Real const pcs_weight
																							){

	if((grid_residue_num_1 < 0)||(grid_residue_num_2 < 0)){
		utility_exit_with_message("Residue num negative. Please review your setup files");
	}
	core::Size grid_residue_num_1_positif(grid_residue_num_1);
	core::Size grid_residue_num_2_positif(grid_residue_num_2);

	grid_edge_ = grid_edge;
	grid_step_ = grid_step;
	grid_small_cutoff_ = grid_small_cutoff;
	grid_large_cutoff_ = 	grid_large_cutoff;
	grid_cone_angle_cutoff_ = grid_cone_angle_cutoff;
	grid_atom_name_1_ = grid_atom_name_1;
	grid_atom_name_2_ = grid_atom_name_2;
	grid_residue_num_1_ = grid_residue_num_1_positif;
	grid_residue_num_2_ = grid_residue_num_2_positif;
	grid_k_vector_ = grid_k_vector;
	minimize_best_tensor_ = minimize_best_tensor;
	pcs_weight_ = pcs_weight;
}

void
PCS_Energy_parameters_manager::print_grid_param() const{
	std::cout <<
		grid_edge_<< " " <<
		grid_step_<<" " <<
		grid_small_cutoff_<<" " <<
		grid_large_cutoff_<<" " <<
		grid_cone_angle_cutoff_<<" " <<
		grid_atom_name_1_<<" " <<
		grid_atom_name_2_<<" " <<
		grid_residue_num_1_<<" " <<
		grid_residue_num_2_<<" " <<
		grid_k_vector_<<" " <<
		minimize_best_tensor_<<" " <<
		pcs_weight_<<" " <<
		std::endl;
}


//	void
//	PCS_Energy_parameters_manager::print_grid_param() const;

core::Real
PCS_Energy_parameters_manager::get_grid_edge() const{
	return grid_edge_;
}

core::Real
PCS_Energy_parameters_manager::get_grid_step() const{
	return grid_step_;
}

core::Real
PCS_Energy_parameters_manager::get_grid_small_cutoff() const{
	return grid_small_cutoff_;
}

core::Real
PCS_Energy_parameters_manager::get_grid_large_cutoff() const{
	return grid_large_cutoff_;
}

core::Real
PCS_Energy_parameters_manager::get_grid_cone_angle_cutoff() const{
	return grid_cone_angle_cutoff_;
}

std::string
PCS_Energy_parameters_manager::get_grid_atom_name_1() const{
	return grid_atom_name_1_;
}

std::string
PCS_Energy_parameters_manager::get_grid_atom_name_2() const{
	return grid_atom_name_2_;
}

core::Size
PCS_Energy_parameters_manager::get_grid_residue_num_1() const{
	return grid_residue_num_1_;
}

core::Size
PCS_Energy_parameters_manager::get_grid_residue_num_2() const{
	return grid_residue_num_2_;
}

core::Real
PCS_Energy_parameters_manager::get_grid_k_vector() const{
	return grid_k_vector_;
}

bool
PCS_Energy_parameters_manager::get_minimize_best_tensor() const{
	return minimize_best_tensor_;
}

core::Real
PCS_Energy_parameters_manager::get_pcs_weight() const{
	return pcs_weight_;
}

utility::vector1<std::string> const &
PCS_Energy_parameters_manager::get_vector_filename() const{
	return vec_filename_;
}

utility::vector1<core::Real> const &
PCS_Energy_parameters_manager::get_vector_weight() const{
	return vec_individual_weight_;
}

PCS_Energy_parameters_manager * PCS_Energy_parameters_manager::instance_( 0 );

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex PCS_Energy_parameters_manager::singleton_mutex_;

std::mutex & PCS_Energy_parameters_manager::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
PCS_Energy_parameters_manager * PCS_Energy_parameters_manager::get_instance()
{
	boost::function< PCS_Energy_parameters_manager * () > creator = boost::bind( &PCS_Energy_parameters_manager::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

PCS_Energy_parameters_manager *
PCS_Energy_parameters_manager::create_singleton_instance()
{
	return new PCS_Energy_parameters_manager;
}


} // pcs
} // methods
} // scoring
} // protocols
