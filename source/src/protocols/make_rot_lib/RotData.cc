// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University

// unit headers
#include <protocols/make_rot_lib/RotData.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace MakeRotLib {

RotData::RotData( Size NumChi, Size NumCluster ) :
  phi_( 0 ),
  psi_( 0 ),
	omega_( 0 ),
	min_omega_( 0 ),
	epsilon_( 0 ),
	min_epsilon_( 0 ),
  energy_( 0 ),
  probability_( 0 ),
  num_chi_( NumChi ),
	num_clusters_( NumCluster ),
  cluster_num_( 0 ),
	// debug
	twist_( 0 ),
	inter_rep_( 0 ),
	inter_atr_( 0 ),
	intra_rep_( 0 ),
	intra_atr_( 0 ),
	solvation_( 0 )
{
  inp_chi_.assign( NumChi, 0 );
  min_chi_.assign( NumChi, 0 );
	lib_chi_val_.assign( NumChi, 0 );
  std_dev_.assign( NumChi, 0 );
  cen_dst_.assign( NumCluster, 0 );
}

} // namespace MakeRotLib
} // namespace protocols
