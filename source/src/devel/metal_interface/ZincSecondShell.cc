// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/ZincSecondShell.cc
/// @brief Computes SASA of a zinc atom of a match grafted onto its scaffold.
/// @details The matches were generated with only a zinc atom its virtual atoms as a transition state, so the clash detection cannot discern if there is room for other zinc-coordinating residues to enter the metal site.
/// @author Bryan Der


#include <devel/metal_interface/ZincSecondShell.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>

#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <sstream>




//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.bder.ZincSecondShell");

typedef numeric::xyzVector<core::Real> point;
typedef point axis;
typedef core::pose::Pose Pose;

using namespace core;

namespace devel {
namespace metal_interface {

///@brief


ZincSecondShell::ZincSecondShell( core::pose::Pose const & pose, utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr)
	: pose_(pose), msr_(msr)
{
	first_shell_atom_ids_.clear();
	first_shell_atom_names_.clear();
	second_shell_atom_ids_.clear();
	second_shell_atom_names_.clear();
	second_shell_atom_hbond_energy_.clear();
	third_shell_atom_ids_.clear();
}

ZincSecondShell::~ZincSecondShell()
{
}

// basic::MetricValue< std::set<Size> >
// ZincSecondShell::get_ligand_neighbors() {
// 	return ligand_neighbors_;
// }
// core::Real
// ZincSecondShell::get_ligand_sasa() {
// 	return ligand_sasa_;
// }


utility::vector1< id::AtomID >
ZincSecondShell::get_first_shell_atom_ids() {
	return first_shell_atom_ids_;
}

utility::vector1< id::AtomID >
ZincSecondShell::get_second_shell_atom_ids() {
	return second_shell_atom_ids_;
}

basic::MetricValue< id::AtomID_Map< Real > >
ZincSecondShell::get_atom_sasa() {
	return atom_sasa_;
}

basic::MetricValue< id::AtomID_Map< Size > >
ZincSecondShell::get_atom_hbonds() {
	return atom_hbonds_;
}


void
ZincSecondShell::register_calculators() {
	using namespace pose::metrics;
	//using namespace protocols::toolbox::pose_metric_calculators;

	TR << "Registering sasa Calculator" << std::endl;
	if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_calc" ) ){
		CalculatorFactory::Instance().register_calculator( "sasa_calc", new pose::metrics::simple_calculators::SasaCalculatorLegacy() );
	}

	TR << "Registering hbond Calculator" << std::endl;
	if( !CalculatorFactory::Instance().check_calculator_exists( "hbond_calc" ) ){
		CalculatorFactory::Instance().register_calculator( "hbond_calc", new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
	}

	return;

}


void
ZincSecondShell::fill_first_shell_atom_ids() {

	for(Size i(2) /*zinc is 1*/; i <= msr_.size(); ++i) {
		Size resnum = msr_[i]->get_seqpos();
		std::string ligand_atom_name = msr_[i]->get_ligand_atom_name();
		first_shell_atom_names_.push_back( ligand_atom_name );

		if(ligand_atom_name == " NE2") { //histidine
			id::AtomID atomid( pose_.residue(resnum).atom_index(" NE2"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " ND1") { //histidine
			id::AtomID atomid( pose_.residue(resnum).atom_index(" ND1"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " OD2") { //aspartate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OD2"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " OD1") { //aspartate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OD1"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " OE2") { //glutamate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OE2"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " OE1") { //glutamate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OE1"), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}
		else if(ligand_atom_name == " SG ") { //cysteine
			id::AtomID atomid( pose_.residue(resnum).atom_index(" SG "), resnum );
			first_shell_atom_ids_.push_back( atomid );
		}

	}

		TR << "first shell size:  " << first_shell_atom_ids_.size() << std::endl;

		for( Size i = 1; i <= first_shell_atom_ids_.size(); ++i ) {
			TR << "  FIRST SHELL : " << first_shell_atom_names_[i] << " " << first_shell_atom_ids_[i].atomno() << " " << first_shell_atom_ids_[i].rsd() << std::endl;
		}

	return;
}


void
ZincSecondShell::fill_second_shell_atom_ids() {

	for(Size i(2) /*zinc is 1*/; i <= msr_.size(); ++i) {
		Size resnum = msr_[i]->get_seqpos();
		std::string ligand_atom_name = msr_[i]->get_ligand_atom_name();
		second_shell_atom_hbond_energy_.push_back(0.0);//initialize hbond energy

		if(ligand_atom_name == " NE2") { //histidine
			id::AtomID atomid( pose_.residue(resnum).atom_index(" ND1"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " ND1" );
		}
		else if(ligand_atom_name == " ND1") { //histidine
			id::AtomID atomid( pose_.residue(resnum).atom_index(" NE2"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " NE2" );
		}
		else if(ligand_atom_name == " OD2") { //aspartate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OD1"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " OD1" );
		}
		else if(ligand_atom_name == " OD1") { //aspartate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OD2"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " OD2" );
		}
		else if(ligand_atom_name == " OE2") { //glutamate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OE1"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " OE1" );
		}
		else if(ligand_atom_name == " OE1") { //glutamate
			id::AtomID atomid( pose_.residue(resnum).atom_index(" OE2"), resnum );
			second_shell_atom_ids_.push_back( atomid );
			second_shell_atom_names_.push_back( " OE2" );
		}
		//no second shell atom id for a Cys residue

	}

	TR << "second shell size: " << second_shell_atom_ids_.size() << std::endl;

	for( Size i = 1; i <= second_shell_atom_ids_.size(); ++i ) {
		TR << "  SECOND SHELL: " << second_shell_atom_names_[i] << " " << second_shell_atom_ids_[i].atomno() << " " << second_shell_atom_ids_[i].rsd() << std::endl;
	}

	return;
}



void
ZincSecondShell::calculate_hbonds_and_sasa(	Pose const & pose ) {

	pdbname_ = pose.pdb_info()->name();

	pose.metric( "hbond_calc", "atom_Hbonds", atom_hbonds_ );
	pose.metric( "sasa_calc", "atom_sasa", atom_sasa_);

	hbond_sasa_have_been_calculated_ = true;

	return;
}


void
ZincSecondShell::report_buried_unsat(	utility::vector1< id::AtomID > atom_ids) {

	if(!hbond_sasa_have_been_calculated_) {
		TR << "YOU CANNOT REPORT ON BURIED UNSAT STATUS UNTIL HBONDS AND SASA HAVE BEEN CALCULATED" << std::endl;
		return;
	}

	std::stringstream ss;
	ss << pdbname_ << " second_shell_summary: ";

	for( Size j = 1; j <= atom_ids.size(); ++j){
		core::id::AtomID atid( atom_ids[j] );
		Size this_atomno( atid.atomno() );
		Size this_resnum( atid.rsd() );

		//TR << "  j: " << j << " atmono: " << this_atomno << " resnum: " << this_resnum << std::endl;

		if( pose_.residue(this_resnum).atom_type( this_atomno ).is_acceptor() || pose_.residue(this_resnum).atom_type( this_atomno ).is_donor() ){
			//we have to add up the sasas for the H attached to this atom
			Real cursasa =  atom_sasa_.value()[ atid ];
			for( Size hcount = pose_.residue(this_resnum).type().attached_H_begin( this_atomno ); hcount<= pose_.residue(this_resnum).type().attached_H_end( this_atomno ); hcount++){
				cursasa = cursasa + atom_sasa_.value()[ core::id::AtomID ( hcount, this_resnum ) ];
			}

			if( cursasa == 0 /*< burial_sasa_cutoff_*/ ){
				Size satisfac_cut = satisfaction_cutoff( pose_.residue(this_resnum).type().atom_type( this_atomno ).name() );
				Size bonded_heavyatoms = pose_.residue(this_resnum).n_bonded_neighbor_all_res( this_atomno ) - pose_.residue(this_resnum).type().number_bonded_hydrogens( this_atomno );
				if( ( bonded_heavyatoms + atom_hbonds_.value()[ atid ] ) < satisfac_cut ){
					TR << pose_.pdb_info()->name() << " buried unsat: " << pose_.residue(this_resnum).type().atom_name( this_atomno ) << " of res " << this_resnum << " has " << atom_sasa_.value()[atid] << " sasa and " << cursasa << " combined sasa, and " << bonded_heavyatoms << " bonded heavyatoms, and " <<  atom_hbonds_.value()[ atid ] << " hbonds, counts as buried unsatisfied." << std::endl;
					ss << atid << " buried_unsat. ";
				}
				else {
					TR << pose_.pdb_info()->name() << " buried and hbonded:" << pose_.residue(this_resnum).type().atom_name( this_atomno ) << " of res " << this_resnum << " has a hydrogen bond." << std::endl;
					ss << atid << " hbonded. ";
				}
			}
			else {
				TR << pose_.pdb_info()->name() << " not buried: " << pose_.residue(this_resnum).type().atom_name( this_atomno ) << " of res " << this_resnum << " has " << atom_sasa_.value()[atid] << " sasa and " << cursasa << " combined sasa." << std::endl;
				ss << atid << " exposed. ";
			}
		}
	}//for

	TR << ss.str() << std::endl;

	return;
}



core::Size
ZincSecondShell::satisfaction_cutoff( std::string atom_type )
{
	if( atom_type == "OH" ) return 2;
	else if (atom_type == "OCbb") return 2;
	else if( atom_type ==  "S")	return 2;
	//everything else we expect to have 3 bonded/h-bonded neighbours to count as satisfied
	else return 3;
}



} //metal_interface
} //devel
