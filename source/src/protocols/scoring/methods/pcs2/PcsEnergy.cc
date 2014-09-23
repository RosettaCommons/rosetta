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
 /// @file protocols/scoring/methods/pcs2/PcsEnergy.cc
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
 /// @last_modified February 2010
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsEnergy.hh>
#include <protocols/scoring/methods/pcs2/PcsEnergyCreator.hh>

// Package headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenterManagerSingleton.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>
#include <protocols/scoring/methods/pcs2/PcsInputCenterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsTensor.hh>
// AUTO-REMOVED #include <protocols/scoring/methods/pcs2/GridSearchIterator.hh>
#include <protocols/scoring/methods/pcs2/GridSearchIteratorCA.hh>
#include <protocols/scoring/methods/pcs2/TensorsOptimizer.hh>
#include <protocols/scoring/methods/pcs2/TensorsOptimizerSvd.hh>
#include <protocols/scoring/methods/pcs2/TensorsOptimizerFix.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/PCS.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/optimization/Minimizer.hh>

#include <core/scoring/EnergyMap.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <fstream>

#include <protocols/scoring/methods/pcs2/PcsInputLine.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

/// @details This must return a fresh instance of the CarbonHBondEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
PcsEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new PcsEnergy );
}

core::scoring::ScoreTypes
PcsEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::pcs2 );
	return sts;
}

static thread_local basic::Tracer TR_PcsEnergy( "protocols.scoring.methods.pcs.PcsEnergy" );

void
PcsEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const{
}

PcsEnergy &
PcsEnergy::operator=(PcsEnergy const & other){
	std::cerr << "Error, = operator not correctly implemented in the PcsEnergy" << std::endl;
	utility_exit_with_message("Exiting");
	if ( this != &other ) {
	}
	return *this;
}

PcsEnergy::PcsEnergy(PcsEnergy const & src ):
	parent( src )
{
	//	TR_PcsEnergy << " () called" << std::endl;
}

PcsEnergy::~PcsEnergy(){
}

/// c-tor
PcsEnergy::PcsEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new PcsEnergyCreator ) )
{
	//	TR_PcsEnergy << "constructor called" << std::endl;
}

/// clone
core::scoring::methods::EnergyMethodOP
PcsEnergy::clone() const{
	//	TR_PcsEnergy << " clone called" << std::endl;
	return core::scoring::methods::EnergyMethodOP( new PcsEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
PcsEnergy::finalize_total_energy(core::pose::Pose & pose,
																 core::scoring::ScoreFunction const &,
																 core::scoring::EnergyMap & totals
																 ) const {

	core::Size i_multi_data;
	/*
	core::Size static kkk(1);

	kkk++;
	if(kkk > 50){
		utility_exit_with_message("EXITING AFTER 50 PCS evals");
	}
	*/

	GridSearchIteratorCA grid_it(pose);
	PcsDataCenterManager & pcs_d_c_m = PCS_multi_data_from_pose(pose); //OK
	//PcsDataCenterManagerSingleton & pcs_d_c_m_s = PCS_multi_data_from_noone(); //DEOSNT SPEED UP

	totals[ core::scoring::pcs2 ] = 0;



	//for ( i_multi_data = 1; i_multi_data <= pcs_d_c_m_s.get_n_multi_data(); ++i_multi_data){
	for ( i_multi_data = 1; i_multi_data <= pcs_d_c_m.get_n_multi_data(); ++i_multi_data){
		//PcsDataCenter & pcs_d_c = (pcs_d_c_m_s.get_PCS_data_all())[i_multi_data];
		PcsDataCenter & pcs_d_c = (pcs_d_c_m.get_PCS_data_all())[i_multi_data];
		pcs_d_c.update_X_Y_Z_all(pose);
		totals[ core::scoring::pcs2 ] += calculate_pcs_score_on_PCS_data_center_CA( pose, false, pcs_d_c, i_multi_data, grid_it);
	}



	TR_PcsEnergy << "Score found: " << totals[core::scoring::pcs2] << std::endl;
}




PcsDataCenterManager &
PcsEnergy::PCS_multi_data_from_pose(core::pose::Pose & pose) const{
	core::Size i_multi_data;
	core::Size n_multi_data;


	n_multi_data = PcsEnergyParameterManager::get_instance()->get_n_multi_data();

	if (
			 ( pose.data().has( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_MULTI_DATA ) ) ){

		//std::cerr << "PcsDataCenterManager was cached" << std::endl<< std::endl;
		return *( utility::pointer::static_pointer_cast< protocols::scoring::methods::pcs2::PcsDataCenterManager > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_MULTI_DATA ) ) );
	}

	//	TR_PcsEnergy << "PcsDataCenterManager was NOT cached" << std::endl;



	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	PcsDataCenterManagerOP pcs_d_c_m_OP;

	pcs_d_c_m_OP = PcsDataCenterManagerOP( new PcsDataCenterManager() );

	for(i_multi_data = 1; i_multi_data <= n_multi_data; ++i_multi_data ){

		utility::vector1<std::string> vec_filename;
		utility::vector1<core::Real> vec_weight;
		vec_filename = PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_vector_filename();
		vec_weight = PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_vector_weight();

		//		TR_PcsEnergy << "Initialization of PcsDataCenter of paramagnetic center " << i_multi_data << " / " << n_multi_data << std::endl;

		core::Size start(PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_include_only_start());

		core::Size end(PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_include_only_end());

		core::Real individual_scale(PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_individual_scale());

		core::Size i2;
		//		TR_PcsEnergy << "File(s) to open: ";
		for( i2 = 1; i2 <= vec_filename.size() ; i2++ ){
			//TR_PcsEnergy << vec_filename[i2] << " " ;
		}
		//TR_PcsEnergy << std::endl;

		if(vec_filename.size() == 0){
			utility_exit_with_message("Missing input file for PCS. Review your setup file");
		}

		PcsInputCenter pcs_i_c = PcsInputCenterManager::get_instance()->get_PcsInputCenter_for(vec_filename, vec_weight);

		//TR_PcsEnergy << pcs_i_c << std::endl;

		PcsDataCenterOP pcs_d_c_OP;
		//			PcsDataCenter pcs_d (pcs_i_c, exclude_residues);
		//pcs_d = PcsDataCenter(pcs_i_c, exclude_residues);
		//pcs_d = PcsDataCenter(pcs_i_c);
		pcs_d_c_OP = PcsDataCenterOP( new PcsDataCenter(pcs_i_c, start, end, individual_scale) );

		//(*pcs_d_c_m_OP).get_PCS_data_all().push_back(pcs_d);
		(*pcs_d_c_m_OP).get_PCS_data_all().push_back(*pcs_d_c_OP);
		//pose.data().set( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_DATA, pcs_d );
	}

	//	(*pcs_d_c_m_OP).set_n_multi_data(n_multi_data);

	pose.data().set( core::pose::datacache::CacheableDataType::PSEUDOCONTACT_SHIFT_MULTI_DATA, pcs_d_c_m_OP );

	//	TR_PcsEnergy << *(PcsInputCenterManager::get_instance());

	return *pcs_d_c_m_OP;
}

PcsDataCenterManagerSingleton &
PcsEnergy::PCS_multi_data_from_noone() const{

	PcsEnergyParameterManager & pcs_e_m =	*(PcsEnergyParameterManager::get_instance());
	PcsDataCenterManagerSingleton & pcs_d_c_m_s = *(PcsDataCenterManagerSingleton::get_instance(pcs_e_m));
	return(pcs_d_c_m_s);
}


void
PcsEnergy::dump_PCS_info(
	utility::vector1<PcsTensor> const & vec_tensor,
	numeric::xyzVector< core::Real > const & best_coo,
	PcsDataCenter const & pcs_d_c
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	core::Size i, j;
	static core::Size n_rescore(1);
	//


	if( option[ basic::options::OptionKeys::PCS::write_extra ].user() ){
		utility::vector1<core::Real> A(5, 0);

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


		const utility::vector1<core::Real> & X_all(pcs_d_c.get_X_all());
		const utility::vector1<core::Real> & Y_all(pcs_d_c.get_Y_all());
		const utility::vector1<core::Real> & Z_all(pcs_d_c.get_Z_all());

		for (i = 1 ; i <= pcs_d_c.get_n_lanthanides(); ++i){
			PcsDataLanthanide const & PCS_d_p_l (pcs_d_c.get_pcs_data_per_lanthanides_all()[i]);
			utility::vector1<PcsInputLine> const & PCS_d_l_a_s (pcs_d_c.get_PCS_data_line_all_spin());

			myfile << "#" << PCS_d_p_l.get_filename() << " RESCORE NUMBER " << n_rescore << std::endl;

			core::Real Xxx(vec_tensor[i].get_chi_xx());
			core::Real Xxy(vec_tensor[i].get_chi_xy());
			core::Real Xxz(vec_tensor[i].get_chi_xz());
			core::Real Xyy(vec_tensor[i].get_chi_yy());
			core::Real Xyz(vec_tensor[i].get_chi_yz());
			const utility::vector1<core::Size> & A_index( PCS_d_p_l.get_A_index());
			//			const ObjexxFCL::FArray1D< core::Real > & fstyle_b(PCS_d_p_l.get_fstyle_b());
			const utility::vector1< core::Real > & cstyle_b(PCS_d_p_l.get_cstyle_b());

			myfile << "# Tensor: " << std::setw(10) << Xxx << " " << Xxy << " "  << Xxz << " " << Xyy << " " << Xyz << " " << best_coo.x() << " " << best_coo.y() << " " << best_coo.z() << std::endl;

			for ( j = 1; j <= PCS_d_p_l.get_n_pcs(); ++j){
				core::Real PCS_exp (cstyle_b[j]);
				core::Size idx (A_index[j]);
				core::Real x (X_all[idx]);
				core::Real y (Y_all[idx]);
				core::Real z (Z_all[idx]);
				core::Size res_num(PCS_d_l_a_s[idx].get_residue_num());
				std::string atom_name(PCS_d_l_a_s[idx].get_atom_name());
				core::Real r5;
				r5 = fill_A_line_fast(A, best_coo.x(), best_coo.y(), best_coo.z(), x, y, z);
				core::Real PCS_calc(A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz);
				PCS_calc = PCS_calc / r5 * FACT_USI_PRECALC_FOR_A_3;
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
PcsEnergy::calculate_pcs_score_on_PCS_data_center_CA(core::pose::Pose & pdb,
																										 bool,
																										 PcsDataCenter & pcs_d_c,
																										 core::Size i_multi_data,
																										 GridSearchIteratorCA & grid_it) const{

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


  utility::vector1<PcsTensor> vec_tensor;
	utility::vector1<core::Real> vec_score;

	core::Real pcs_score_total;
	core::Size i;


	core::Real pcs_weight_center (PcsEnergyParameterManager::get_instance()->get_PcsEnergyParameter_for(i_multi_data).get_pcs_weight());

	if (pcs_weight_center == 0){
		return 0;
	}

	//alloc best score and vector and coordinate
	for ( i = 1; i <= pcs_d_c.get_n_lanthanides(); ++i){
		PcsTensor PCS_t = PcsTensor(0, 0, 0, 0, 0, ((pcs_d_c.get_pcs_data_per_lanthanides_all())[i]).get_filename());
		vec_tensor.push_back(PCS_t);
	}
	vec_score.resize(pcs_d_c.get_n_lanthanides());

	numeric::xyzVector< core::Real > best_coo;

	//call to calculate the tensors and the score
	pcs_score_total = CA_search_scores_and_tensors(vec_score, vec_tensor, best_coo, pdb, pcs_d_c, i_multi_data, grid_it);

	//	pcs_score_total = CA_search_scores_and_tensors_with_svd(vec_score, vec_tensor, best_coo, pdb, pcs_d_c, i_multi_data, grid_it);

	//	std::cerr << "x y z " << best_coo.x() << " " << best_coo.y() << " " << best_coo.z() << std::endl;

	dump_PCS_info(vec_tensor, best_coo, pcs_d_c);

	//This weight is the weight for the center, defined for each stage.
	return (pcs_weight_center * pcs_score_total);
}


core::Real
PcsEnergy::minimize_tensors_from_PCS_data(	utility::vector1<PcsTensor> & vec_best_tensor,
																						numeric::xyzVector< core::Real > & best_coo,
																						PcsDataCenter const & /*pcs_d_c*/,
																						core::optimization::Minimizer & minimizer,
																						utility::vector1<core::Real> & vect_to_opti
																						) const{

	core::Size i;

	core::Real optimized_cost(minimizer.run( vect_to_opti ));

	best_coo.assign(vect_to_opti[1], vect_to_opti[2], vect_to_opti[3]);

	for(i = 1; i <= vec_best_tensor.size(); ++i){
		vec_best_tensor[i].reset_tensor((core::Real)vect_to_opti[3 + 5*(i-1) + 1],
																		(core::Real)vect_to_opti[3 + 5*(i-1) + 2],
																		(core::Real)vect_to_opti[3 + 5*(i-1) + 3],
																		(core::Real)vect_to_opti[3 + 5*(i-1) + 4],
																		(core::Real)vect_to_opti[3 + 5*(i-1) + 5]);
	}
	return (optimized_cost);
}

core::Real
PcsEnergy::minimize_tensors_from_PCS_data_with_svd(	utility::vector1<PcsTensor> & /*vec_best_tensor*/,
																										numeric::xyzVector< core::Real > & best_coo,
																										PcsDataCenter const & /*pcs_d_c*/,
																										core::optimization::Minimizer & minimizer,
																										utility::vector1<core::Real> & vect_to_opti
																										) const{

	//	core::Size i;

	core::Real optimized_cost(minimizer.run( vect_to_opti ));

	best_coo.assign(vect_to_opti[1], vect_to_opti[2], vect_to_opti[3]);

	return (optimized_cost);
}


core::Real
PcsEnergy::minimize_tensors_fix_from_PCS_data(	utility::vector1<PcsTensor> & vec_best_tensor,
																								PcsDataCenter const & pcs_d_c/*,
																								core::Real xM,
																								core::Real yM,
																								core::Real zM*/
																								) const{

	core::Size i;

	utility::vector1<core::Real> vect_to_opt;

	if(vec_best_tensor.size() !=  pcs_d_c.get_n_lanthanides()){
		utility_exit_with_message("n_lanthanides and vec_best_tensor size differs in minimize_tensors_from_PCS_data");
	}

	for(i = 1; i <= vec_best_tensor.size(); ++i){
		vect_to_opt.push_back(vec_best_tensor[i].get_chi_xx());
		vect_to_opt.push_back(vec_best_tensor[i].get_chi_xy());
		vect_to_opt.push_back(vec_best_tensor[i].get_chi_xz());
		vect_to_opt.push_back(vec_best_tensor[i].get_chi_yy());
		vect_to_opt.push_back(vec_best_tensor[i].get_chi_yz());
	}
	//FIX
	TensorsOptimizerFix tensors_opt_fix(pcs_d_c);
	//	optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone_atol", 0.0000001, true, false, false );
	//	optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	core::optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	core::optimization::Minimizer minimizer(tensors_opt_fix, options );

	core::Real optimized_cost(minimizer.run( vect_to_opt ));


	for(i = 1; i <= vec_best_tensor.size(); ++i){
		vec_best_tensor[i].reset_tensor((core::Real)vect_to_opt[5*(i-1) + 1],
																		(core::Real)vect_to_opt[5*(i-1) + 2],
																		(core::Real)vect_to_opt[5*(i-1) + 3],
																		(core::Real)vect_to_opt[5*(i-1) + 4],
																		(core::Real)vect_to_opt[5*(i-1) + 5]);
	}
	return (optimized_cost);
}




//This will be called for each new pose
core::Real
PcsEnergy::CA_search_scores_and_tensors(utility::vector1<core::Real> & vec_best_score,
																				utility::vector1<PcsTensor> & vec_best_tensor,
																				numeric::xyzVector< core::Real > & best_coo,
																				core::pose::Pose const & /* pdb*/,
																				PcsDataCenter & pcs_d_c,
																				core::Size,
																				GridSearchIteratorCA & grid_it) const{

	core::Real x, y, z;
	core::Real best_score/*, score*/;  // redefined later ~Labonte
	core::Size i;
	core::Size size_of;
	utility::vector1<core::Real> vec_score_temp;
	utility::vector1<PcsTensor> vec_tensor_temp;

	//ref to do atomic best switch in order to avoid copy stuff
	utility::vector1<core::Real> * vec_score_ref_current;
	utility::vector1<PcsTensor> * vec_tensor_ref_current;

	//CA_SEARCH
	TensorsOptimizer tensors_opt(pcs_d_c);
	//optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone_atol", 0.0000001, true, false, false );
	//	optimization::MinimizerOptions options( "dfpmin", 0.00001, true, false, false );
	//GOOD optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	//	optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	core::optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	core::optimization::Minimizer minimizer(tensors_opt, options );

	//some basic checking...
	size_of = vec_best_score.size();
	if((size_of != vec_best_tensor.size())||
		 (size_of != pcs_d_c.get_n_lanthanides())){
		std::cerr << "Problem in CA_search_scores_and_tensors function" << std::endl;
		std::cerr << "n_lanthanides =  " << pcs_d_c.get_n_lanthanides();
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

	grid_it.reset();
	best_coo.assign(0, 0, 0);

	core::Size i_junk;
	utility::vector1<PcsTensor> vec_tensor_junk;
	utility::vector1<core::Real> vec_score_junk;
	for ( i_junk = 1; i_junk <= pcs_d_c.get_n_lanthanides(); ++i_junk){
		PcsTensor PCS_t_junk = PcsTensor(0, 0, 0, 0, 0, ((pcs_d_c.get_pcs_data_per_lanthanides_all())[i_junk]).get_filename());
		vec_tensor_junk.push_back(PCS_t_junk);
	}

	numeric::xyzVector< core::Real >  temp_coo;
	best_score = 999999999999999999999999999.9;
	bool test_at_least_one_iteration = false;
	core::Size n_best_found(0);
	core::Real epsilon(0.001);

	utility::vector1<core::Real> vect_to_opti;

	vect_to_opti.push_back(0);
	vect_to_opti.push_back(0);
	vect_to_opti.push_back(0);
	core::Size ii;
	for(ii = 1; ii <= pcs_d_c.get_n_lanthanides(); ++ii){
		vect_to_opti.push_back(0);
		vect_to_opti.push_back(0);
		vect_to_opti.push_back(0);
		vect_to_opti.push_back(0);
		vect_to_opti.push_back(0);
	}

	core::Size n_identical(1);
	PcsTensor pcs_t(0, 0, 0, 0, 0, "");
	utility::vector1<PcsDataLanthanide> & pcs_d_l_vec(pcs_d_c.get_pcs_data_per_lanthanides_all());

	while(grid_it.next_center(x, y, z) == true){
		//std::cerr << "Trying " << x << " " << y << " " << z << std::endl;
		if(n_best_found > n_identical){
			TR_PcsEnergy << "Found " << n_best_found << " the same value -> best score" << std::endl;
			break;
		}

		test_at_least_one_iteration = true;
		//		std::cerr << "Trying " << x << " " << y << " " << z << std::endl;
		temp_coo.assign(x, y, z);
		//score = 0;


		vect_to_opti[1] = x;
		vect_to_opti[2] = y;
		vect_to_opti[3] = z;

		//X_Y_Z is updated
		pcs_d_c.update_matrix_A_all_for_svd(x, y, z);

		for ( ii = 1; ii <= pcs_d_c.get_n_lanthanides(); ++ii){

			PcsDataLanthanide & pcs_d_l(pcs_d_l_vec[ii]);
			pcs_d_l.calculate_tensor_only_with_svd(pcs_t);

			vect_to_opti[3 + (ii-1)*5 + 1] = pcs_t.get_chi_xx(); //0;
			vect_to_opti[3 + (ii-1)*5 + 2] = pcs_t.get_chi_xy(); //0;
			vect_to_opti[3 + (ii-1)*5 + 3] = pcs_t.get_chi_xz(); //0;
			vect_to_opti[3 + (ii-1)*5 + 4] = pcs_t.get_chi_yy(); //0;
			vect_to_opti[3 + (ii-1)*5 + 5] = pcs_t.get_chi_yz(); //0;
		}


		core::Real score(minimize_tensors_from_PCS_data(vec_tensor_junk, temp_coo, pcs_d_c, minimizer, vect_to_opti));
		//		std::cerr << "Temp found: " << score << std::endl;
		//		std::cerr << "Score: " << score << std::endl;
		//		std::cerr << "Found " << temp_coo.x() << " " << temp_coo.y() << " " << temp_coo.z() << std::endl;

		if(fabs(score - best_score) < epsilon){
			n_best_found++;
		}

		if ( score < best_score){
			if(best_score - score > epsilon){
				n_best_found = 0;
			}
			for(i = 1; i <= pcs_d_c.get_n_lanthanides(); ++i){
				(*vec_tensor_ref_current)[i] = vec_tensor_junk[i];
				(*vec_score_ref_current)[i] = score /  pcs_d_c.get_n_lanthanides();
			}
		}

		if ( score < best_score){
			best_score = score;

			best_coo.assign(temp_coo.x(), temp_coo.y(), temp_coo.z());

			//atomic switch
			if((vec_score_ref_current != &vec_score_temp) && (vec_score_ref_current != &vec_best_score)){ //test to make sure...
				std::cerr << "Problem in grid_search_scores_and_tensors function" << std::endl;
				std::cerr << "The atomic switch is not working (1)" << std::endl;
				utility_exit_with_message("Exiting");
			}
			if(vec_tensor_ref_current == &vec_best_tensor){
				vec_tensor_ref_current = &vec_tensor_temp;
				vec_score_ref_current = &vec_score_temp;
			}
			else{
				if(vec_tensor_ref_current != &vec_tensor_temp){ //test to make sure...
					std::cerr << "Problem in grid_search_scores_and_tensors function" << std::endl;
					std::cerr << "The atomic switch is not working (2)" << std::endl;
					utility_exit_with_message("Exiting");
				}
				vec_tensor_ref_current = &vec_best_tensor;
				vec_score_ref_current = &vec_best_score;
			}
			//end atomic switch
		}
	} //while


	if(!(n_best_found > n_identical)){
		//		TR_PcsEnergy << " All point visited " << std::endl;
	}

	if(test_at_least_one_iteration == false){
		utility_exit_with_message("The description of the grid search given is too restrictive");
	}

	if( vec_tensor_ref_current == &vec_best_tensor ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_tensor[i].reset_from_ref(vec_tensor_temp[i]);
		}
	}

	if( vec_score_ref_current == &vec_best_score ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_score[i] = vec_score_temp[i];
		}
	}
	return (best_score);
}



core::Real
PcsEnergy::CA_search_scores_and_tensors_with_svd(utility::vector1<core::Real> & vec_best_score,
																								utility::vector1<PcsTensor> & vec_best_tensor,
																								numeric::xyzVector< core::Real > & best_coo,
																								core::pose::Pose const & /* pdb*/,
																								PcsDataCenter & pcs_d_c,
																								core::Size,
																								GridSearchIteratorCA & grid_it) const{

	core::Real x, y, z;
	core::Real best_score/*, score*/;  // redefined later ~Labonte
	core::Size i;
	core::Size size_of;
	utility::vector1<core::Real> vec_score_temp;
	utility::vector1<PcsTensor> vec_tensor_temp;

	//ref to do atomic best switch in order to avoid copy stuff
	utility::vector1<core::Real> * vec_score_ref_current;
	utility::vector1<PcsTensor> * vec_tensor_ref_current;

	TensorsOptimizerSvd tensors_opt_svd(pcs_d_c);
	//CA_SVD_SEARCH
	//optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone_atol", 0.0000001, true, false, false );
	//	optimization::MinimizerOptions options( "dfpmin", 0.00001, true, false, false );
	core::optimization::MinimizerOptions options( "dfpmin", 0.000001, true, false, false );
	core::optimization::Minimizer minimizer(tensors_opt_svd, options );

	//some basic checking...
	size_of = vec_best_score.size();
	if((size_of != vec_best_tensor.size())||
		 (size_of != pcs_d_c.get_n_lanthanides())){
		std::cerr << "Problem in CA_search_scores_and_tensors function" << std::endl;
		std::cerr << "n_lanthanides =  " << pcs_d_c.get_n_lanthanides();
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

	grid_it.reset();
	best_coo.assign(0, 0, 0);

	core::Size i_junk;
	utility::vector1<PcsTensor> vec_tensor_junk;
	utility::vector1<core::Real> vec_score_junk;
	for ( i_junk = 1; i_junk <= pcs_d_c.get_n_lanthanides(); ++i_junk){
		PcsTensor PCS_t_junk = PcsTensor(0, 0, 0, 0, 0, ((pcs_d_c.get_pcs_data_per_lanthanides_all())[i_junk]).get_filename());
		vec_tensor_junk.push_back(PCS_t_junk);
	}


	numeric::xyzVector< core::Real >  temp_coo;
	best_score = 999999999999999999999999999.9;
	bool test_at_least_one_iteration = false;
	core::Size n_best_found(0);
	core::Real epsilon(0.0001);


	utility::vector1<core::Real> vect_to_opti;

	vect_to_opti.push_back(0);
	vect_to_opti.push_back(0);
	vect_to_opti.push_back(0);

	core::Size n_identical(1);

	while(grid_it.next_center(x, y, z) == true){
		//		std::cerr << "Trying " << x << " " << y << " " << z << std::endl;
		if(n_best_found > n_identical){
			TR_PcsEnergy << "Found " << n_best_found << " the same value -> best score" << std::endl;
			break;
		}

		test_at_least_one_iteration = true;
		//		std::cerr << "Trying " << x << " " << y << " " << z << std::endl;
		temp_coo.assign(x, y, z);
		//score = 0;


		vect_to_opti[1] = x;
		vect_to_opti[2] = y;
		vect_to_opti[3] = z;

		core::Real score(minimize_tensors_from_PCS_data_with_svd(vec_tensor_junk, temp_coo, pcs_d_c, minimizer, vect_to_opti));

		//		std::cerr << "Score: " << score << std::endl;
		//		std::cerr << "Found " << temp_coo.x() << " " << temp_coo.y() << " " << temp_coo.z() << std::endl;

		if(fabs(score - best_score) < epsilon){
			n_best_found++;
		}

		if ( score < best_score){
			if(best_score - score > epsilon){
				n_best_found = 0;
			}
			for(i = 1; i <= pcs_d_c.get_n_lanthanides(); ++i){
				(*vec_tensor_ref_current)[i] = vec_tensor_junk[i];
				(*vec_score_ref_current)[i] = score /  pcs_d_c.get_n_lanthanides();
			}
		}

		if ( score < best_score){
			best_score = score;

			best_coo.assign(temp_coo.x(), temp_coo.y(), temp_coo.z());

			//atomic switch
			if((vec_score_ref_current != &vec_score_temp) && (vec_score_ref_current != &vec_best_score)){ //test to make sure...
				std::cerr << "Problem in grid_search_scores_and_tensors function" << std::endl;
				std::cerr << "The atomic switch is not working (1)" << std::endl;
				utility_exit_with_message("Exiting");
			}
			if(vec_tensor_ref_current == &vec_best_tensor){
				vec_tensor_ref_current = &vec_tensor_temp;
				vec_score_ref_current = &vec_score_temp;
			}
			else{
				if(vec_tensor_ref_current != &vec_tensor_temp){ //test to make sure...
					std::cerr << "Problem in grid_search_scores_and_tensors function" << std::endl;
					std::cerr << "The atomic switch is not working (2)" << std::endl;
					utility_exit_with_message("Exiting");
				}
				vec_tensor_ref_current = &vec_best_tensor;
				vec_score_ref_current = &vec_best_score;
			}
			//end atomic switch
		}
	} //while


	if(!(n_best_found > n_identical)){
		TR_PcsEnergy << " All point visited " << std::endl;
	}

	if(test_at_least_one_iteration == false){
		utility_exit_with_message("The description of the grid search given is too restrictive");
	}

	if( vec_tensor_ref_current == &vec_best_tensor ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_tensor[i].reset_from_ref(vec_tensor_temp[i]);
		}
	}

	if( vec_score_ref_current == &vec_best_score ){
 		for(i = 1; i <= vec_best_score.size(); ++i){
			vec_best_score[i] = vec_score_temp[i];
		}
	}
	return (best_score);
}
core::Size
PcsEnergy::version() const
{
	return 1; // Initial versioning
}


} // PCS
} // methods
} // scoring
} // protocols
