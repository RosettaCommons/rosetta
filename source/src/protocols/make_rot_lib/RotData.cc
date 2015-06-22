// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// unit headers
#include <protocols/make_rot_lib/RotData.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace make_rot_lib {

static thread_local basic::Tracer TR( "protocols.make_rot_lib.RotData" );

	
bool operator==( RotData & r1, RotData & r2 ) {

	if ( r1.get_num_bbs() != r2.get_num_bbs()) {
		return false;
	}
	if (r1.get_bbs() == r2.get_bbs() ) {
		return false;
	}
	for ( core::Size i = 1; i <= r1.get_num_bbs(); ++i ) {
		if (r1.get_bb_id( i ) != r2.get_bb_id( i ) ) {
			return false;
		}
	}
	if ( r1.get_omg() != r2.get_omg() ) {
		return false;
	}
	if ( r1.get_min_omg() != r2.get_min_omg() ) {
		return false;
	}
	
	if ( r1.get_eps() != r2.get_eps() ) {
		return false;
	}
	if ( r1.get_min_eps() != r2.get_min_eps() ) {
		return false;
	}
	
	
	if ( r1.get_energy() != r2.get_energy() ) {
		return false;
	}
	
	if ( r1.get_probability() != r2.get_probability() ) {
		return false;
	}
	
	if ( r1.get_num_chi() != r2.get_num_chi() ) {
		return false;
	}

	if ( r1.get_num_clusters() != r2.get_num_clusters() ) {
		return false;
	}
	
	if ( r1.get_num_chi() != r2.get_num_chi() ) {
		return false;
	}
	
	core::Size num_chi = r1.get_num_chi();
	for ( core::Size i = 1; i <= num_chi; ++i ) {
		
		if ( r1.get_inp_chi(i) != r2.get_inp_chi(i) ) {
			return false;
		}
		
		if ( r1.get_min_chi(i) != r2.get_min_chi(i) ) {
			return false;
		}
		
		if ( r1.get_lib_chi_val(i) != r2.get_lib_chi_val(i) ) {
			return false;
		}
		
		if ( r1.get_std_dev(i) != r2.get_std_dev(i) ) {
			return false;
		}
	}
	
	return true;
}

RotData::RotData( core::Size NumChi, core::Size NumCluster ) :
	num_bbs_( 2 ),
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
	solvation_( 0 ),
	semirotameric_( false )
{
	bbs_.resize( num_bbs_ );
	inp_chi_.assign( NumChi, 0 );
	min_chi_.assign( NumChi, 0 );
	lib_chi_val_.assign( NumChi, 0 );
	std_dev_.assign( NumChi, 0 );
	cen_dst_.assign( NumCluster, 0 );
	
	semi_energy_dist_.assign( 36, 0 );
	semi_prob_dist_.assign( 36, 0 );
}

RotData::RotData( core::Size NumChi, core::Size NumCluster, bool semirotameric ) :
	num_bbs_( 2 ),
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
	solvation_( 0 ),
	semirotameric_( semirotameric )
{
	bbs_.resize( num_bbs_ );
	inp_chi_.assign( NumChi, 0 );
	min_chi_.assign( NumChi, 0 );
	lib_chi_val_.assign( NumChi, 0 );
	std_dev_.assign( NumChi, 0 );
	cen_dst_.assign( NumCluster, 0 );
	
	semi_energy_dist_.assign( 37, 0 );
	semi_prob_dist_.assign( 37, 0 );
}
	
RotData::RotData( core::Size NumChi, core::Size NumBBs, core::Size NumCluster ) :
	num_bbs_( NumBBs ),
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
	solvation_( 0 ),
	semirotameric_( false )
{
    bbs_.resize( num_bbs_ );
    bb_ids_.resize( num_bbs_ );
    inp_chi_.assign( NumChi, 0 );
    min_chi_.assign( NumChi, 0 );
    lib_chi_val_.assign( NumChi, 0 );
    std_dev_.assign( NumChi, 0 );
	cen_dst_.assign( NumCluster, 0 );
	
	semi_energy_dist_.assign( 36, 0 );
	semi_prob_dist_.assign( 36, 0 );
}

RotData::RotData( core::Size NumChi, core::Size NumBBs, core::Size NumCluster, bool semirotameric ) :
	num_bbs_( NumBBs ),
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
	solvation_( 0 ),
	semirotameric_( semirotameric )
{
	bbs_.resize( num_bbs_ );
	bb_ids_.resize( num_bbs_ );
	inp_chi_.assign( NumChi, 0 );
	min_chi_.assign( NumChi, 0 );
	lib_chi_val_.assign( NumChi, 0 );
	std_dev_.assign( NumChi, 0 );
	cen_dst_.assign( NumCluster, 0 );
		
	semi_energy_dist_.assign( 36, 0 );
	semi_prob_dist_.assign( 36, 0 );
}
	

/// @brief Output function, primarily for debugging purposes.
void
RotData::show( std::ostream & out ) const {
	out << "RotData: ";
	out << num_bbs_ << " ";
	for ( core::Size bb_i = 1; bb_i <= bbs_.size(); ++bb_i ) {
		out << bbs_[ bb_i ] << " ";
	} // out << psi_ << " ";
	out << omega_ << " ";
	out << epsilon_ << " ";
	out << energy_ << " ";
	out << probability_ << " ";
	out << num_chi_ << " ";
	out << num_clusters_ << " ";
	out << cluster_num_ << " ";
	out << min_omega_ << " ";
	out << min_epsilon_;
 
	if ( semirotameric_ ) {
		out << " ";
		for ( core::Size i = 1; i <= 36; ++i ) {
			out << semi_energy_dist_[i] << " ";
		}
		for ( core::Size i = 1; i <= 36; ++i ) {
			out << semi_prob_dist_[i] << " ";
		}
		
	}
	out << std::endl;

	assert( inp_chi_.size() == num_chi_ );
	assert( min_chi_.size() == num_chi_ );
	assert( lib_chi_val_.size() == num_chi_ );
	assert( std_dev_.size() == num_chi_ );
	for( core::Size ii(1); ii <= num_chi_; ++ii ) {
		out << "chi " << ii << " " << inp_chi_[ii] << " " << min_chi_[ii] << " " << lib_chi_val_[ii] << " " << std_dev_[ii] << std::endl;
	}
	assert( cen_dst_.size() == num_clusters_ );
	out << "cen_dist";
	for( core::Size jj(1); jj <= cen_dst_.size(); ++jj ) {
		out << " " << cen_dst_[jj];
	}
	out << std::endl;
}

/// @brief input function, primarily for debugging purposes.
/// @details Return true on failure.
bool
RotData::load( std::istream & in ) {
	std::string tag;
	in >> tag;
	if( tag != "RotData:" ) {
		TR.Warning << "Expected 'RotData:', found '"<< tag << "'" << std::endl;
		return true;
	}
	in >> num_bbs_;
	for (core::Size bb_i = 1; bb_i <= num_bbs_; ++bb_i ) {
		in >> bbs_[ bb_i ];
	}
	in >> omega_;
	in >> epsilon_;
	in >> energy_;
	in >> probability_;
	in >> num_chi_;
	in >> num_clusters_;
	in >> cluster_num_;
	in >> min_omega_;
	in >> min_epsilon_;

	bbs_.resize( num_bbs_ );
	inp_chi_.resize( num_chi_ );
	min_chi_.resize( num_chi_ );
	lib_chi_val_.resize( num_chi_ );
	std_dev_.resize( num_chi_ );
	for ( core::Size ii( 1 ); ii <= num_chi_; ++ii ) {
		core::Size index;
		in >> tag;
		if ( tag != "chi" ) {
			TR.Warning << "Expected 'chi' found '" << tag << "'" << std::endl;
			return true;
		}
		in >> index;
		if ( index != ii ) {
			TR.Warning << "Expected " << ii << " found " << index << std::endl;
			return true;
		}
		in >> inp_chi_[ ii ] >> min_chi_[ ii ] >> lib_chi_val_[ ii ] >> std_dev_[ ii ];
	}
	in >> tag;
	if ( tag != "cen_dist" ) {
		TR.Warning << "Expected 'cen_dist' found '" << tag << "'" << std::endl;
		return true;
	}
	cen_dst_.resize( num_clusters_ );
	for ( core::Size jj( 1 ); jj <= num_clusters_; ++jj ) {
		in >> cen_dst_[ jj ];
	}
	return ! in.good();
}

} // namespace MakeRotLib
} // namespace protocols
