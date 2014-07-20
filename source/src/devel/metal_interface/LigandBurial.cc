// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/LigandBurial.cc
/// @brief Determines the SASA of the ligand, and number of neighbors of the ligand as a measure of burial.
/// @details
/// @author Bryan Der


#include <devel/metal_interface/LigandBurial.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/sasa.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <sstream>




//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.bder.LigandBurial");

typedef numeric::xyzVector<core::Real> point;
typedef point axis;
typedef core::pose::Pose Pose;
typedef std::set< core::Size > SetSize;

using namespace core;

namespace devel {
namespace metal_interface {

///@brief


LigandBurial::LigandBurial( core::pose::Pose const & pose, std::string ligand_3_letter_code )
	: pose_(pose), ligand_3_letter_code_(ligand_3_letter_code)
{
	if(ligand_3_letter_code == "ZN") {
	 	ligand_3_letter_code_ = " ZN";
	}
}

LigandBurial::~LigandBurial()
{
}

basic::MetricValue< std::set<Size> >
LigandBurial::get_ligand_neighbors() {
	return ligand_neighbors_;
}
core::Real
LigandBurial::get_ligand_sasa() {
	return ligand_sasa_;
}


core::Size
LigandBurial::find_ligand() {

	//make this safe for ZN

	for ( Size i(1); i <= pose_.total_residue(); ++i ) {
		//if ( pose_.residue(i).name3() == ligand_3_letter_code_ || pose_.residue(i).name3() == " ZN" ) {
		if ( pose_.residue(i).name3() == ligand_3_letter_code_ ) {
			ligand_resnum_ = i;
			return ligand_resnum_;
		}
	}


	//else
	TR << "Did not find the ligand" << std::endl;

	return 0;
}


void
LigandBurial::register_calculators() {

	using namespace pose::metrics;
	using namespace pose::metrics::simple_calculators;


	std::string const calc_stem("ligand_nbr_calc_");
	calcname_.str("");
	calcname_ << calc_stem << ligand_resnum_;

	TR << "Registering NeighborsByDistance Calculator " << calcname_.str() << std::endl;
	if( !CalculatorFactory::Instance().check_calculator_exists( calcname_.str() ) ){
		CalculatorFactory::Instance().register_calculator( calcname_.str(), new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( ligand_resnum_ ) );
	}

	TR << "Registering SASA Calculator" << std::endl;
	if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc" ) ){
		CalculatorFactory::Instance().register_calculator( "sasa_calc", new SasaCalculatorLegacy() );
	}

	return;
}


void
LigandBurial::calculate_ligand_neighbors() {

	//number of neighbors of ligand
	basic::MetricValue< Size > num_n;
	pose_.metric( calcname_.str(), "num_neighbors", num_n);
	TR << pose_.pdb_info()->name() << " " << ligand_3_letter_code_ << " ligand has " << num_n.value() << " neighbors." << std::endl;

	//list each neighbor individually
	basic::MetricValue< SetSize > ligand_neighbors;
	pose_.metric( calcname_.str(), "neighbors", ligand_neighbors );


	TR << "These are the neighbor residues: " << std::endl;
	for( std::set<Size>::const_iterator it(ligand_neighbors.value().begin()), end(ligand_neighbors.value().end()); it != end; ++it){
		TR << *it << "+";
	}
	TR << std::endl;

	ligand_neighbors_ = ligand_neighbors;

	return;
}


void
LigandBurial::calculate_ligand_sasa() {

	//assume the ligand is the second of 2 chains.
	//possible to append 1 residue to new pose?
	PoseOP chain_ligand_only = pose_.split_by_chain(2);

	using core::id::AtomID;
	utility::vector1<Real> rsd_sasa(pose_.n_residue(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose_,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose_,false);
	for(Size i = 1; i <= pose_.n_residue(); i++) {
		for(Size j = 1; j <= pose_.residue(i).nheavyatoms(); j++) {
			atom_mask[AtomID(j,i)] = true;
		}
	}
	Real const probe_radius = 1.4;
	core::scoring::calc_per_atom_sasa( pose_, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );

	TR << pose_.pdb_info()->name() << "  SASA value for " << ligand_3_letter_code_ << " ligand " << ligand_resnum_ << ": " << rsd_sasa[ligand_resnum_] << std::endl;

	ligand_sasa_ = rsd_sasa[ligand_resnum_];

	return;

	// SASA of HIZ residue alone is about 320 square-Angstroms
	// chain_ligand_only.dump_pdb("ligand_only.pdb");

}




} //metal_interface
} //devel
