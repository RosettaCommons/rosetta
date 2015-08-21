// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.cc
/// @brief  SasaCalculatorLegacy class
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/sasa.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>

// Utility headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>

#include <utility/assert.hh>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {


void InterfaceSasaDefinitionCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "delta_sasa" ) {
		basic::check_cast( valptr, &delta_sasa_, "delta_sasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( delta_sasa_ );

	} else if ( key == "frac_ch1_dsasa" ) {
		basic::check_cast( valptr, &delta_sasa_, "frac_ch1_dsasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( fraction_chain1_delta_sasa_ );

	} else if ( key == "frac_ch2_dsasa" ) {
		basic::check_cast( valptr, &delta_sasa_, "frac_ch2_dsasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( fraction_chain2_delta_sasa_ );

	} else if ( key == "delta_atom_sasa" ) {
		basic::check_cast( valptr, &atom_delta_sasa_, "delta_atom_sasa expects to return a id::AtomID_Map< Real >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< Real > > *>(valptr))->set( atom_delta_sasa_ );

	} else if ( key == "delta_residue_sasa" ) {
		basic::check_cast( valptr, &residue_delta_sasa_, "delta_residue_sasa expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue< utility::vector1< Real > > *>(valptr))->set( residue_delta_sasa_ );

	} else if ( key == "interface_atoms" ) {
		basic::check_cast( valptr, &interface_atoms_, "interface_atoms expects to return a id::AtomID_Map< bool >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< bool > > *>(valptr))->set( interface_atoms_ );

	} else if ( key == "interface_residues" ) {
		basic::check_cast( valptr, &interface_residues_, "interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( interface_residues_ );

	} else if ( key == "first_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain1_interface_residues_, "first_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain1_interface_residues_ );

	} else if ( key == "second_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain2_interface_residues_, "second_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain2_interface_residues_ );

	} else if ( key == "num_interface_residues" ) {
		basic::check_cast( valptr, &num_interface_residues_, "num_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_interface_residues_ );

	} else if ( key == "num_first_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain1_interface_residues_, "num_first_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain1_interface_residues_ );

	} else if ( key == "num_second_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain2_interface_residues_, "num_second_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain2_interface_residues_ );

	} else if ( key == "first_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch1_begin_num_, "first_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_begin_num_ );

	} else if ( key == "first_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch1_end_num_, "first_chain_last_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_end_num_ );

	} else if ( key == "second_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch2_begin_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_begin_num_ );

	} else if ( key == "second_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch2_end_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_end_num_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string InterfaceSasaDefinitionCalculator::print( std::string const & key ) const {

	if ( key == "delta_sasa" ) {
		return utility::to_string( delta_sasa_ );
	} else if ( key == "frac_ch1_dsasa" ) {
		return utility::to_string( fraction_chain1_delta_sasa_ );
	} else if ( key == "frac_ch2_dsasa" ) {
		return utility::to_string( fraction_chain2_delta_sasa_ );
	} else if ( key == "delta_atom_sasa" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "delta_residue_sasa" ) {
		return utility::to_string( residue_delta_sasa_ );
	} else if ( key == "interface_atoms" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "first_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "second_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "num_interface_residues" ) {
		return utility::to_string( num_interface_residues_ );
	} else if ( key == "num_first_chain_interface_residues" ) {
		return utility::to_string( num_chain1_interface_residues_ );
	} else if ( key == "num_second_chain_interface_residues" ) {
		return utility::to_string( num_chain2_interface_residues_ );
	} else if ( key == "first_chain_first_resnum" ) {
		return utility::to_string( ch1_begin_num_ );
	} else if ( key == "first_chain_last_resnum" ) {
		return utility::to_string( ch1_end_num_ );
	} else if ( key == "second_chain_first_resnum" ) {
		return utility::to_string( ch2_begin_num_ );
	} else if ( key == "second_chain_last_resnum" ) {
		return utility::to_string( ch2_end_num_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


void InterfaceSasaDefinitionCalculator::recompute( Pose const & this_pose ) {

	verify_chain_setup( this_pose );

	Real const probe_radius( basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] );
	Real const atom_delta_sasa_thres(0.2);   // atom must lose this much sasa to be buried by the interface
	Real const residue_delta_sasa_thres(2);   // residue must lose this much sasa to be buried by the interface

	// Build atom subsets for computing SASA
	// note: we'll compute sasa ONLY for the chains of interest, ignoring any other chains which may be present!
	id::AtomID_Map< bool > bound_atom_subset, chain1_atom_subset, chain2_atom_subset;
	bound_atom_subset.clear();
	chain1_atom_subset.clear();
	chain2_atom_subset.clear();
	core::pose::initialize_atomid_map( bound_atom_subset, this_pose, false );
	core::pose::initialize_atomid_map( chain1_atom_subset, this_pose, false );
	core::pose::initialize_atomid_map( chain2_atom_subset, this_pose, false );
	for ( Size ir = ch1_begin_num_; ir <= ch1_end_num_; ++ir ) {
		core::conformation::Residue const & irsd( this_pose.residue( ir ) );
		for ( Size ia = 1; ia <= irsd.natoms(); ++ia ) {
			id::AtomID const iid( ia, ir );
			bound_atom_subset[ iid ] = true;
			chain1_atom_subset[ iid ] = true;
		}
	}
	for ( Size ir = ch2_begin_num_; ir <= ch2_end_num_; ++ir ) {
		core::conformation::Residue const & irsd( this_pose.residue( ir ) );
		for ( Size ia = 1; ia <= irsd.natoms(); ++ia ) {
			id::AtomID const iid( ia, ir );
			bound_atom_subset[ iid ] = true;
			chain2_atom_subset[ iid ] = true;
		}
	}

	// Get the bound SASA, and the SASA for chain1/chain2 alone
	id::AtomID_Map< Real > bound_atom_sasa, chain1_atom_sasa, chain2_atom_sasa;
	utility::vector1< Real > bound_residue_sasa, chain1_residue_sasa, chain2_residue_sasa;
	Real bound_sasa = core::scoring::calc_per_atom_sasa( this_pose, bound_atom_sasa, bound_residue_sasa, probe_radius, false, bound_atom_subset );
	Real chain1_sasa = core::scoring::calc_per_atom_sasa( this_pose, chain1_atom_sasa, chain1_residue_sasa, probe_radius, false, chain1_atom_subset );
	Real chain2_sasa = core::scoring::calc_per_atom_sasa( this_pose, chain2_atom_sasa, chain2_residue_sasa, probe_radius, false, chain2_atom_subset );

	// Fill in "delta_sasa" member data, and put together a lists of interface atoms and interface residues
	delta_sasa_ = chain1_sasa + chain2_sasa - bound_sasa;

	//std::cerr << "total ch1 sasa is " << chain1_sasa << ", total ch2 sasa is " << chain2_sasa << ", bound sasa is " << bound_sasa << ", and delta sasa is " << delta_sasa_ << std::endl;
	atom_delta_sasa_.clear();
	core::pose::initialize_atomid_map( atom_delta_sasa_, this_pose, 0. );
	interface_atoms_.clear();
	core::pose::initialize_atomid_map( interface_atoms_, this_pose, false );
	interface_residues_.clear();
	chain1_interface_residues_.clear();
	chain2_interface_residues_.clear();
	residue_delta_sasa_.clear();
	residue_delta_sasa_.assign( this_pose.total_residue(), 0. );

	core::Real delta_sasa_atomch1(0.0), delta_sasa_atomch2(0.0);

	for ( Size ir = ch1_begin_num_; ir <= ch1_end_num_; ++ir ) {
		core::conformation::Residue const & irsd( this_pose.residue( ir ) );
		for ( Size ia = 1; ia <= irsd.natoms(); ++ia ) {
			id::AtomID const iid( ia, ir );
			atom_delta_sasa_[iid] = chain1_atom_sasa[ iid ] - bound_atom_sasa[ iid ];
			delta_sasa_atomch1 += atom_delta_sasa_[iid];
			if ( atom_delta_sasa_[iid] > atom_delta_sasa_thres ) {
				interface_atoms_[ iid ] = true;
			}
		}
		residue_delta_sasa_[ir] = chain1_residue_sasa[ir] - bound_residue_sasa[ir];
		if ( residue_delta_sasa_[ir] > residue_delta_sasa_thres ) {
			interface_residues_.insert( ir );
			chain1_interface_residues_.insert( ir );
		}
	}

	for ( Size ir = ch2_begin_num_; ir <= ch2_end_num_; ++ir ) {
		core::conformation::Residue const & irsd( this_pose.residue( ir ) );
		for ( Size ia = 1; ia <= irsd.natoms(); ++ia ) {
			id::AtomID const iid( ia, ir );
			atom_delta_sasa_[iid] = chain2_atom_sasa[ iid ] - bound_atom_sasa[ iid ];
			delta_sasa_atomch2 += atom_delta_sasa_[iid];
			if ( atom_delta_sasa_[iid] > atom_delta_sasa_thres ) {
				interface_atoms_[ iid ] = true;
			}
		}
		residue_delta_sasa_[ir] = chain2_residue_sasa[ir] - bound_residue_sasa[ir];
		if ( residue_delta_sasa_[ir] > residue_delta_sasa_thres ) {
			interface_residues_.insert( ir );
			chain2_interface_residues_.insert( ir );
		}
	}

	fraction_chain1_delta_sasa_ = delta_sasa_atomch1 / chain1_sasa;
	fraction_chain2_delta_sasa_ = delta_sasa_atomch2 / chain2_sasa;

	//std::cerr << "summing up atomd delta sasas, ch1 delta sasa is " << delta_sasa_atomch1 << " and ch2 delta sasa is " << delta_sasa_atomch2 << std::endl;
	num_interface_residues_ = interface_residues_.size();
	num_chain1_interface_residues_ = chain1_interface_residues_.size();
	num_chain2_interface_residues_ = chain2_interface_residues_.size();

	return;
}

} // simple_calculators
} // metrics
} // pose
} // core
