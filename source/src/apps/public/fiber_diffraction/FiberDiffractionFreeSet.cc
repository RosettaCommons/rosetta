// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/wojtek/FiberDiffractionFreeSet.cc
/// @brief  FiberDiffraction free set generator
/// @author Wojciech Potrzebowski and Ingemar Andre


#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

//Fiber diffraction modules
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/chemical/AA.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>

#include <numeric/random/random.hh>

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>
#include <numeric/util.hh>
//#include <math.h>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/basic_sys_util.hh>
// option key includes


using namespace core;
using namespace utility;
static basic::Tracer TR("apps.pilot.wojtek.FiberDiffractionFreeSet");

int
main( int argc, char * argv [] ) {

	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// initialize core
		devel::init(argc, argv);

		utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
		utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
		utility::vector0< utility::vector0 < int > >::iterator nvals;
		core::Size lmax, Rmax;
		core::Real a_(0), b_(0), c_(0), p_(0), radius_(0);

		core::Real res_cutoff_low_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_low ]();
		core::Real res_cutoff_high_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_high ]();

		if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user() ) {
			a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
		}

		if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user() ) {
			b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
		}

		if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user() ) {
			p_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ]();
		}

		if ( !(a_ > 0 ) ) {
			utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
		}

		if ( !(b_ > 0 ) ) {
			utility_exit_with_message("The number of turns, score::fiber_diffraction::b, must be set!");
		}
		if ( !(p_ > 0 ) ) {
			utility_exit_with_message("Helical pitch, score::fiber_diffraction::p, must be set!");
		}
		c_ = p_*a_;


		core::scoring::fiber_diffraction::getFiberDiffractionData(c_, res_cutoff_high_, res_cutoff_low_).getAllFiberData(layer_lines_I, layer_lines_R, nvals, lmax, Rmax);
		TR << "Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;


		TR<<"Calculating sampling points"<<std::endl;

		if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::radius ].user() ) {
			radius_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::radius ]();
		}

		if ( !(radius_ > 0 ) ) {
			utility_exit_with_message( "Helical radius, score::fiber_diffraction::radius, must be set!");
		}
		core::Real max_r_value(radius_);
		core::Real structure_cutoff(4*M_PI*max_r_value);
		TR<<"Maximum radius and cutoff: "<< max_r_value <<" "<<structure_cutoff<< std::endl;

		utility::vector0< utility::vector1< core::Real > > sampling_points_lE;
		utility::vector0< core::Size > sampling_points_number_l;
		utility::vector0< core::Size > lowest_bessel_orders_l;
		utility::vector0< core::Size > lowest_resolution_l;
		utility::vector0< core::Size > highest_resolution_l;
		utility::vector0< utility::vector1< core::Size > > selected_R_l;
		utility::vector0< utility::vector1< core::Real > >selected_Rinv_l;
		utility::vector0< utility::vector1< core::Real > > work_Rinv_l;
		utility::vector0< int >  work_R_l;
		utility::vector0< core::Size > free_selected_l;

		work_Rinv_l.resize(lmax+1);
		work_R_l.resize(lmax+1);
		for ( Size l=0; l<= lmax; ++l ) {
			work_R_l[l]=0;
			work_Rinv_l[l].resize(layer_lines_R[l].size());
			for ( Size i=1; i<= layer_lines_R[l].size(); ++i ) {
				work_Rinv_l[l][i]=-100.0;
			}
		}

		core::Size chi_iterations_(20);
		if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::chi_iterations ].user() ) {
			chi_iterations_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::chi_iterations ]();
		}

		core::Real free_prcnt_thr(0.2);
		//Start from lowest order and repeat util free set is below 20%
		int n_minus;
		core::Size selection_iters(0);
		core::Size independent_Rs(0);
		core::Real free_set_percent(1), free_set_percent_l(1), max_percent(0);
		utility::vector0<int>::iterator position;
		while ( free_set_percent > free_prcnt_thr && selection_iters< chi_iterations_ ) {
			core::scoring::fiber_diffraction::bessel_roots(lmax, c_, res_cutoff_high_, res_cutoff_low_, structure_cutoff, sampling_points_lE,\
				sampling_points_number_l, lowest_bessel_orders_l, highest_resolution_l, lowest_resolution_l,\
				layer_lines_R, nvals);

			core::scoring::fiber_diffraction::interpolate_sampled_to_grid(lmax, sampling_points_lE,sampling_points_number_l,\
				highest_resolution_l, lowest_resolution_l, layer_lines_R, selected_R_l, selected_Rinv_l);
			//Remove lowest_bessel_order
			//check percentage for each layer line
			max_percent = 0;
			for ( Size l=0; l <= lmax; ++l ) {
				if ( sampling_points_number_l[l]==0 ) continue;
				independent_Rs = 0;
				//free_set_percent_l = 1.0-float(work_R_l[l])/float(highest_resolution_l[l]-lowest_resolution_l[l]);
				for ( Size i=1; i<= sampling_points_number_l[l]; ++i ) {
					if ( work_Rinv_l[l][selected_R_l[l][i]]==-100.0 ) independent_Rs++;
					work_Rinv_l[l][selected_R_l[l][i]]=selected_Rinv_l[l][i];
					//Supplementing working set with neighboring points
					if ( fabs(sampling_points_lE[l][i]-layer_lines_R[l][selected_R_l[l][i]+1]) < fabs(sampling_points_lE[l][i]-layer_lines_R[l][selected_R_l[l][i]-1]) ) {
						if ( (selected_R_l[l][i]+1)<highest_resolution_l[l] ) {
							if ( work_Rinv_l[l][selected_R_l[l][i]+1]==-100.0 ) independent_Rs++;
							work_Rinv_l[l][selected_R_l[l][i]+1] = layer_lines_R[l][selected_R_l[l][i]+1];
						}
					} else {
						if ( (selected_R_l[l][i]-1)>lowest_resolution_l[l] ) {
							if ( work_Rinv_l[l][selected_R_l[l][i]-1]==-100.0 ) independent_Rs++;
							work_Rinv_l[l][selected_R_l[l][i]-1] = layer_lines_R[l][selected_R_l[l][i]-1];
						}
					}
					///////////////////////////////////////////Up and down///////////////////////////////////
					/*if ((selection_iters == chi_iterations_-1) && (free_set_percent_l>free_prcnt_thr)) {
					if (selected_R_l[l][i]+1<highest_resolution_l[l]) {
					if (work_Rinv_l[l][selected_R_l[l][i]+1]==-100.0) { independent_Rs++; TR<<"Condition sattisfied"<<std::endl;}
					work_Rinv_l[l][selected_R_l[l][i]+1] = layer_lines_R[l][selected_R_l[l][i]+1];
					}
					if (selected_R_l[l][i]-1>lowest_resolution_l[l]) {
					if (work_Rinv_l[l][selected_R_l[l][i]-1]==-100.0) { independent_Rs++; TR<<"Condition sattisfied"<<std::endl;}
					work_Rinv_l[l][selected_R_l[l][i]-1] = layer_lines_R[l][selected_R_l[l][i]-1];
					}
					}*/
					////////////////////////////////////////////////////////////////////////////////////////
				}
				work_R_l[l]+=independent_Rs;
				free_set_percent_l = 1.0-float(work_R_l[l])/float(highest_resolution_l[l]-lowest_resolution_l[l]);
				TR<<"Indep R "<<work_R_l[l]<<" : "<<free_set_percent_l<<std::endl;
				if ( free_set_percent_l>max_percent ) max_percent = free_set_percent_l;
				//If free precentage for each layer line is higher than 20% than next lowest bessel orders are taken
				if ( free_set_percent_l>free_prcnt_thr ) {
					position = std::find(nvals[l].begin(), nvals[l].end(), lowest_bessel_orders_l[l]);
					if ( position != nvals[l].end() && nvals[l].size()!=1 ) {
						nvals[l].erase(position);
					}
					n_minus = -lowest_bessel_orders_l[l];
					position = std::find(nvals[l].begin(), nvals[l].end(), n_minus);
					if ( position != nvals[l].end() && nvals[l].size()!=1 ) {
						nvals[l].erase(position);
					}
				}
			}
			free_set_percent = max_percent;
			selection_iters++;
		}

		TR<<"Generating random points"<<std::endl;
		//Final check for the lines which are above 20%
		//Add some extra points to a free set at random
		Size  random_R(1), selection_rounds(0);
		Size number_of_frees(0);
		for ( Size l=0; l <= lmax; ++l ) {
			free_set_percent_l = 1.0-float(work_R_l[l])/float(highest_resolution_l[l]-lowest_resolution_l[l]);
			number_of_frees = (highest_resolution_l[l]-lowest_resolution_l[l])-work_R_l[l]+1;
			//Randomly select
			independent_Rs = 0;
			if ( free_set_percent_l>free_prcnt_thr ) {
				selection_rounds = number_of_frees-ceil((free_prcnt_thr/free_set_percent_l)*number_of_frees);
				free_selected_l.resize(number_of_frees);

				int number_of_frees_R =0;
				for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
					if ( layer_lines_R[l][R] == 0.0 ) continue;
					if ( layer_lines_R[l][R]*layer_lines_R[l][R]+(l/c_)*(l/c_) < res_cutoff_low_*res_cutoff_low_ || \
							layer_lines_R[l][R]*layer_lines_R[l][R]+(l/c_)*(l/c_) > res_cutoff_high_*res_cutoff_high_ ) continue;
					if ( work_Rinv_l[l][R]==-100.0 ) {
						free_selected_l[number_of_frees_R]=R;
						number_of_frees_R++;
					}
				}

				for ( Size bin=0; bin<selection_rounds; ++bin ) {
					//maximum_number of trials
					Size max_trials = 0;
					while ( work_Rinv_l[l][free_selected_l[random_R]]!=-100.0 && max_trials<10 ) {
						random_R = numeric::random::rg().random_range( 1, free_selected_l.size() );
						max_trials++;
					}
					if ( work_Rinv_l[l][free_selected_l[random_R]]==-100.0 ) independent_Rs++;
					else TR<<"Random R "<<random_R<<" : "<<free_selected_l[random_R]<<" : "<<work_Rinv_l[l][free_selected_l[random_R]] <<std::endl;
					work_Rinv_l[l][free_selected_l[random_R]] = layer_lines_R[l][free_selected_l[random_R]];
				}
			}
			work_R_l[l]+=independent_Rs;
			free_set_percent_l = 1.0-float(work_R_l[l])/float(highest_resolution_l[l]-lowest_resolution_l[l]);
			TR<<"Free percentage for layer line "<<l<<" : "<<free_set_percent_l<<std::endl;
		}

		//Writting to a file
		std::string free_outfile = "free_layer_lines";
		std::string work_outfile = "work_layer_lines";
		std::ofstream out_free, out_work;
		out_free.open(free_outfile.c_str(), std::ios::out);
		out_work.open(work_outfile.c_str(), std::ios::out);
		TR<<"Writting to file... "<<std::endl;
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if ( layer_lines_R[l][R] == 0.0 ) continue;
				if ( layer_lines_R[l][R]*layer_lines_R[l][R]+(l/c_)*(l/c_) < res_cutoff_low_*res_cutoff_low_ || \
						layer_lines_R[l][R]*layer_lines_R[l][R]+(l/c_)*(l/c_) > res_cutoff_high_*res_cutoff_high_ ) continue;
				if ( work_Rinv_l[l][R]==-100.0 ) {
					out_free << layer_lines_I[l][R]<<" "<< layer_lines_R[l][R] << " " << l << std::endl;
				} else {
					out_work << layer_lines_I[l][R]<<" "<< layer_lines_R[l][R] << " " << l << std::endl;
				}
			}
		}
		out_free.close();
		out_work.close();
		TR<<"Done... "<<std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
