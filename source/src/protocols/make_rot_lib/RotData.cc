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
#include <basic/Tracer.hh>

namespace protocols {
namespace make_rot_lib {

static basic::Tracer TR("protocols.make_rot_lib.RotData");

RotData::RotData( core::Size NumChi, core::Size NumCluster ) :
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

/// @brief Output function, primarily for debugging purposes.
void
RotData::show( std::ostream & out ) const {
	out << "RotData: ";
	out << phi_ << " ";
	out << psi_ << " ";
	out << omega_ << " ";
	out << epsilon_ << " ";
	out << energy_ << " ";
  out << probability_ << " ";
	out << num_chi_ << " ";
	out << num_clusters_ << " ";
	out << cluster_num_ << " ";
  out << min_omega_ << " ";
  out << min_epsilon_ << std::endl;

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
	in >> phi_;
	in >> psi_;
	in >> omega_;
	in >> epsilon_;
	in >> energy_;
  in >> probability_;
	in >> num_chi_;
	in >> num_clusters_;
	in >> cluster_num_;
  in >> min_omega_;
  in >> min_epsilon_;

	inp_chi_.resize(num_chi_);
	min_chi_.resize(num_chi_);
	lib_chi_val_.resize(num_chi_);
	std_dev_.resize(num_chi_);
	for( core::Size ii(1); ii <= num_chi_; ++ii ) {
		core::Size index;
		in >> tag;
		if( tag != "chi" ) {
			TR.Warning << "Expected 'chi' found '" << tag << "'" << std::endl;
			return true;
		}
		in >> index;
		if( index != ii ) {
			TR.Warning << "Expected " << ii << " found " << index << std::endl;
			return true;
		}
		in >> inp_chi_[ii] >> min_chi_[ii] >> lib_chi_val_[ii] >> std_dev_[ii];
	}
	in >> tag;
	if( tag != "cen_dist") {
		TR.Warning << "Expected 'cen_dist' found '" << tag << "'" << std::endl;
		return true;
	}
	cen_dst_.resize(num_clusters_);
	for( core::Size jj(1); jj <= num_clusters_; ++jj ) {
		in >> cen_dst_[jj];
	}
	return ! in.good();
}

} // namespace MakeRotLib
} // namespace protocols
