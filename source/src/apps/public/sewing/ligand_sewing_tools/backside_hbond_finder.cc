// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/zinc_statistic_generator.cc
/// @brief Builds zinc ions off coordinating residues in different positions to build a library of allowed zinc positions relative to the parent residue
/// @details Identify metal binding histidines and attempt to build backside hydrogen bonds (SER, THR, or ASN) to non-binding nitrogens that aren't solvent exposed
/// @author Sharon Guffy



//Devel
#include <devel/init.hh>
//Protocols
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
//Core
#include <core/util/metalloproteins_util.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
//#include <protocols/residue_selectors/LigandMetalContactSelector.hh>
#include <core/select/util/SelectResiduesByLayer.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/AtomGraph.hh>
#include <core/conformation/AtomGraphData.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

//Basic
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/ReferenceCount.hh>
//Numeric
#include <numeric/xyzVector.hh>


using basic::Error;
using basic::Warning;
//using basic::T;
static basic::Tracer TR("apps.pilot.guffysl.backside_hbond_finder");



//Define local options
namespace local {
//We'll determine solvent exposure based on the number of neighbors (rough approximation)
basic::options::RealOptionKey const atom_sasa_cutoff( "atom_sasa_cutoff" );
basic::options::RealOptionKey const neighbor_distance_cutoff( "neighbor_distance_cutoff" );
struct HbondInfo{
	core::Size resnum;
	std::string residue_type;
	utility::fixedsizearray1< core::Real, core::pack::dunbrack::DUNBRACK_MAX_SCTOR > chi_angles;
	core::Real hbond_score;
	HbondInfo( core::Size a1, std::string a2, utility::fixedsizearray1< core::Real, core::pack::dunbrack::DUNBRACK_MAX_SCTOR > a3, core::Real a4 ){
		resnum = a1;
		residue_type = a2;
		chi_angles = a3;
		hbond_score = a4;
	}
	friend bool operator<( HbondInfo const & a, HbondInfo const & b ){
		//Compare first on Hbond score
		if ( a.hbond_score < b.hbond_score ) return true;
		if ( a.hbond_score > b.hbond_score ) return false;
		//Then on resnum
		if ( a.resnum < b.resnum ) return true;
		if ( a.resnum > b.resnum ) return false;
		//Then on residue type
		if ( a.residue_type < b.residue_type ) return true;
		if ( a.residue_type > b.residue_type ) return false;
		//I really can't see it getting past hbond score, though
		return false;
	}
};
}

//Header file information
class BacksideHbondFinderMover: public protocols::moves::Mover{
public:
	//Default constructor
	BacksideHbondFinderMover();
	//Copy constructor
	BacksideHbondFinderMover (BacksideHbondFinderMover const & other);
	//Destructor
	virtual ~BacksideHbondFinderMover();
	//Set private data
	//Get private data
	//Virtual methods
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	void mutate_residue( core::pose::Pose & pose, core::Size res_position, std::string target_res_name, bool );
	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;
	core::Real get_atom_sasa_cutoff() const{ return atom_sasa_cutoff_; }
	void set_atom_sasa_cutoff( core::Size cutoff ){ atom_sasa_cutoff_ = cutoff; }
	core::Real get_neighbor_distance_cutoff() const{ return neighbor_distance_cutoff_; }
	void set_neighbor_distance_cutoff( core::Real cutoff ){ neighbor_distance_cutoff_ = cutoff; }
	void
	check_for_hbonds( core::pose::Pose & pose, core::Size his_resnum, core::Size test_resnum, utility::fixedsizearray1< core::Real, 5 > const & bb_angles, std::string res_to_mutate, std::set< local::HbondInfo > & backside_hbond_res );
	//For RosettaScripts
	//Any other methods
private:
	core::Size atom_sasa_cutoff_;
	core::Real neighbor_distance_cutoff_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::set< core::Size > off_limits_residues_;
	//We need both energy and evaluator for each
	core::scoring::methods::EnergyMethodOP vdw_energy_;
	core::scoring::hbonds::HBondEnergyOP hbond_;
	//core::scoring::hbonds::HBondEnergyOP hbond_eval_;
	//private data
};
//Owning pointers
typedef utility::pointer::shared_ptr< BacksideHbondFinderMover > BacksideHbondFinderMoverOP;
typedef utility::pointer::shared_ptr< BacksideHbondFinderMover const > BacksideHbondFinderMoverCOP;
//Define methods
BacksideHbondFinderMover::BacksideHbondFinderMover(){
	//Initialize options
	scorefxn_ = core::scoring::get_score_function();
	atom_sasa_cutoff_ = basic::options::option[ local::atom_sasa_cutoff ].value();
	neighbor_distance_cutoff_ = basic::options::option[ local::neighbor_distance_cutoff ].value();
	//Define energy methods that we'll want to use for scoring
	core::scoring::methods::EnergyMethodOptions options( scorefxn_->energy_method_options() );
	//core::scoring::etable::Etable etable( * (core::scoring::ScoringManager::get_instance()->etable( options ).lock() ) );
	//van der Waals terms (fa_rep, fa_atr, fa_sol)
	if ( options.analytic_etable_evaluation() ) {
		vdw_energy_= core::scoring::methods::EnergyMethodOP( new core::scoring::etable::AnalyticEtableEnergy(  * (core::scoring::ScoringManager::get_instance()->etable( options ).lock() ), options, false ) );
	} else {
		vdw_energy_ = core::scoring::methods::EnergyMethodOP( new core::scoring::etable::TableLookupEtableEnergy(  * (core::scoring::ScoringManager::get_instance()->etable( options ).lock() ), options, false ) );
	}
	//lk_ball energy
	//core::scoring::lkball::LK_BallEnergyOP lk_ball( new core::scoring::lkball::LK_BallEnergy(options ) );
	//hbond energy
	hbond_ = core::scoring::hbonds::HBondEnergyOP( new core::scoring::hbonds::HBondEnergy( options.hbond_options() ) );
	off_limits_residues_.clear();
}

BacksideHbondFinderMover::BacksideHbondFinderMover(BacksideHbondFinderMover const & other):
	protocols::moves::Mover(other)
{
}
BacksideHbondFinderMover::~BacksideHbondFinderMover(){}

//Virtual methods
protocols::moves::MoverOP BacksideHbondFinderMover::clone() const{
	return protocols::moves::MoverOP( new BacksideHbondFinderMover(*this) );
}
protocols::moves::MoverOP BacksideHbondFinderMover::fresh_instance() const{
	return protocols::moves::MoverOP( new BacksideHbondFinderMover);
}
std::string BacksideHbondFinderMover::get_name() const{
	return "BacksideHbondFinderMover";
}


void
BacksideHbondFinderMover::mutate_residue( core::pose::Pose & pose, core::Size res_position, std::string target_res_name, bool only_dihedrals ){
	TR.Debug << "Begin mutate residue" << std::endl;
	core::chemical::ResidueTypeSetCOP restype_set( pose.conformation().residue_type_set_for_conf() );
	core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue(
		restype_set->name_map( target_res_name ), pose.residue( res_position ), pose.conformation() );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( res_position ), *new_res, pose.conformation(), only_dihedrals );
	pose.replace_residue( res_position, *new_res, false );
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only( res_position );
	TR.Debug << "End mutate residue" << std::endl;
}

void BacksideHbondFinderMover::apply(core::pose::Pose & pose){
	off_limits_residues_.clear();
	TR.Debug << "Begin set up energy methods" << std::endl;
	//Set up the energy methods
	//WILL THIS WORK?
	scorefxn_->score( pose );
	if ( scorefxn_->energy_method_options().analytic_etable_evaluation() ) {
		std::dynamic_pointer_cast< core::scoring::etable::AnalyticEtableEnergy >( vdw_energy_ )->setup_for_scoring( pose, *scorefxn_ );
	} else {
		std::dynamic_pointer_cast< core::scoring::etable::TableLookupEtableEnergy >( vdw_energy_ )->setup_for_scoring( pose, *scorefxn_ );
	}
	hbond_->setup_for_scoring( pose, *scorefxn_ );
	//lk_ball->setup_for_scoring( pose, *scorefxn );

	//Fill the residue PointGraph
	core::conformation::PointGraphOP pose_point_graph( new core::conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pose_point_graph );
	core::conformation::find_neighbors< core::conformation::PointGraphVertexData, core::conformation::PointGraphEdgeData >
		( pose_point_graph, neighbor_distance_cutoff_ );
	//Fill the atom graph
	core::conformation::AtomGraphOP pose_atom_graph( new core::conformation::AtomGraph );
	core::conformation::atom_graph_from_conformation( pose.conformation(), *pose_atom_graph );
	if ( !pose_atom_graph->get_vertex( 1 ).data().residue_id() ) {
		TR << "AtomGraph did NOT properly initialize." << std::endl;
	} else {
		TR << "AtomGraph loaded data" << std::endl;
	}
	//For neighbor calculations, we're actually going to use the code in SelectResiduesByLayer (??? Should we really even use this anymore?)
	core::select::util::SelectResiduesByLayer layer_select( false, true, true ); //Select boundary and surface residues
	utility::vector1< core::Size > not_core = layer_select.compute( pose );
	//utility::vector1< core::Size > const & surface = layer_select.selected_surface_residues();
	//Not sure if we should consider boundary residues or not




	//Go ahead and get atom SASA values, which we will check
	core::scoring::sasa::SasaCalcOP atom_sasa_calc( new core::scoring::sasa::SasaCalc );
	atom_sasa_calc->calculate( pose );
	core::id::AtomID_Map< core::Real > atom_sasa = atom_sasa_calc->get_atom_sasa();

	TR.Debug << "Found neighbors" << std::endl;
	//I think the best thing to do here is to use residue selectors rather than trying to do everything manually with the AtomGraph
	//First find zinc ion
	std::set< core::Size > zinc_indices;
	core::select::residue_selector::ResidueNameSelectorOP zinc_selector( new core::select::residue_selector::ResidueNameSelector );
	zinc_selector->set_residue_name3( " ZN" );
	core::select::residue_selector::ResidueSubset zinc_subset = zinc_selector->apply( pose );
	for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
		if ( zinc_subset.at( resnum ) ) {
			zinc_indices.insert( resnum );
			off_limits_residues_.insert( resnum );
		}
	}
	//Then find residues bound to zinc ion
	std::set< core::Size > his_indices;

	/*
	core::select::residue_selector::LigandMetalContactSelector bonded;
	bonded.set_input_set_selector( zinc_selector );
	core::select::residue_selector::ResidueSubset his_subset = bonded.apply( pose );
	for( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ){
	if( his_subset.at( resnum ) && zinc_indices.count( resnum ) == 0 ){
	his_indices.insert( resnum );
	off_limits_residues_.insert( resnum );
	}
	}
	*/

	//This is the part where we fix the bad his tautomers
	for ( core::Size zn_resnum: zinc_indices ) {
		utility::vector1< core::id::AtomID > coordinating_atoms = core::util::find_metalbinding_atoms( pose, zn_resnum, 1.0 );
		for ( core::id::AtomID coord_id: coordinating_atoms ) {
			his_indices.insert( coord_id.rsd() );
			off_limits_residues_.insert( coord_id.rsd() );
			TR << "Atom name: " << pose.residue( coord_id.rsd() ).atom_name( coord_id.atomno() ) << std::endl;
			//If the coordinating residue is not a histidine (HIS or HIS_D), ignore it
			/*
			if( pose.residue_type( coord_id.rsd() ).name3() != "HIS" ){
			TR << "Skipping non-histidine coordinating residue" << std::endl;
			continue;
			}
			*/
			//Now to get down to the mutation
			//If the coordinating atom is ND1, the residue should be HIS
			//If the coordinating atom is NE2, the residue should be HIS_D
			if ( pose.residue( coord_id.rsd() ).atom_name( coord_id.atomno() ) == " ND1" ) {
				if ( pose.residue_type( coord_id.rsd() ).base_name() != "HIS" ) {
					//mutate TODO
					TR << "Mutating HIS_D to HIS" << std::endl;
					mutate_residue( pose, coord_id.rsd(), "HIS", false );
				}
			}
			if ( pose.residue( coord_id.rsd() ).atom_name( coord_id.atomno() ) == " NE2" ) {
				if ( pose.residue_type( coord_id.rsd() ).base_name() != "HIS_D" ) {
					//mutate TODO
					TR << "Mutating HIS to HIS_D" << std::endl;
					mutate_residue( pose, coord_id.rsd(), "HIS_D", false );
				}
			}
		}
	}




	TR.Debug << "Selected zinc and histidine residues. Preparing to loop." << std::endl;
	core::Size num_his_cleared = 0;
	core::Size hbonds_found = 0;
	core::Real total_hbond_score = 0;


	utility::vector1< core::id::AtomID > all_metalbinding_atoms;
	for ( core::Size zn_resnum: zinc_indices ) {
		utility::vector1< core::id::AtomID > coordinating_atoms = core::util::find_metalbinding_atoms( pose, zn_resnum, 1.0 );
		all_metalbinding_atoms.insert( all_metalbinding_atoms.end(), coordinating_atoms.begin(), coordinating_atoms.end() );
	}


	for ( core::Size his_resnum: his_indices ) {
		TR.Debug << "Outer loop: his_resnum " << his_resnum << std::endl;
		/*
		if( std::find( surface.begin(), surface.end(), his_resnum ) == surface.end() ){
		TR << "Histidine " << his_resnum << " is solvent-exposed!" << std::endl;
		++num_his_cleared;
		continue;
		}
		*/
		//Otherwise this is either a boundary residue or a core residue
		//I suspect that some of the boundary residues may have solvent-exposed second nitrogens. I should check for that.
		//It's also possible that some non-HIS residues will be included
		utility::vector1< core::Real > his_atom_sasa = atom_sasa[ his_resnum ];
		core::Real backside_atom_sasa = 0;
		core::id::AtomID backside_atom_id;


		//Find all of this residue's metal binding atoms
		//Each atom can only coordinate 1 metal, so if a residue is coordinating more than once, skip it & count it as satisfied (assumes no backbone coordination)
		//This will also take care of the bidentate Asp/Glu case
		core::id::AtomID metalbinding_atom;
		bool found = false;
		bool satisfied = false;
		for ( core::id::AtomID coord_id: all_metalbinding_atoms ) {
			if ( coord_id.rsd() == his_resnum ) {
				if ( found ) {
					TR << "Residue " << his_resnum << " is satisfied by coordinating with multiple atoms!" << std::endl;
					++num_his_cleared;
					satisfied = true;
					break;
				} else {
					found = true;
					metalbinding_atom = coord_id;
				}
			}
		}
		if ( satisfied ) {
			continue;
		}




		if ( pose.residue_type( his_resnum ).name3() == "HIS" ) {
			if ( metalbinding_atom.atomno() == 10 ) { //NE2 is coordinating (HIS_D), so ND1 (atom 7 ) is backside
				//get sasa for nd1
				backside_atom_sasa = his_atom_sasa.at( 7 );
				//set atom id to (his_resnum, nd1_atomnum)
				backside_atom_id = core::id::AtomID( 7, his_resnum );
			} else {
				//get sasa for ne2
				backside_atom_sasa = his_atom_sasa.at( 10 );
				//set atom id to (ne2_atomnum, his_resnum)
				backside_atom_id = core::id::AtomID( 10, his_resnum );
			}
		} else if ( pose.residue_type( his_resnum ).name3() == "ASP" ) {
			if ( metalbinding_atom.atomno() == 8 ) { //od2 is coordinating, so oD1 (atom 7 ) is backside
				//get sasa for nd1
				backside_atom_sasa = his_atom_sasa.at( 7 );
				//set atom id to (his_resnum, nd1_atomnum)
				backside_atom_id = core::id::AtomID( 7, his_resnum );
			} else {
				//get sasa for ne2
				backside_atom_sasa = his_atom_sasa.at( 8 );
				//set atom id to (ne2_atomnum, his_resnum)
				backside_atom_id = core::id::AtomID( 8, his_resnum );
			}
		} else if ( pose.residue_type( his_resnum ).name3() == "GLU" ) {

			if ( metalbinding_atom.atomno() == 9 ) { //oe2 is coordinating, so oe1 (atom 8 ) is backside
				//get sasa for nd1
				backside_atom_sasa = his_atom_sasa.at( 8 );
				//set atom id to (his_resnum, nd1_atomnum)
				backside_atom_id = core::id::AtomID( 8, his_resnum );
			} else {
				//get sasa for ne2
				backside_atom_sasa = his_atom_sasa.at( 9 );
				//set atom id to (ne2_atomnum, his_resnum)
				backside_atom_id = core::id::AtomID( 9, his_resnum );
			}
		} else { //must be a cysteine. We don't want to deal with these right now
			TR << "WARNING: The specified binding site contains cysteine, which is not supported by this protocol." << std::endl;
			++num_his_cleared;
			continue;
		}
		if ( backside_atom_sasa > atom_sasa_cutoff_ ) { //Potentially enough surface area for water to coordinate
			TR << "Backside atom of histidine " << his_resnum << " is solvent-exposed!" << std::endl;
			++num_his_cleared;
			continue;
		}

		//NOW: look for backside hydrogen bonds
		//We'll want to store the position, residue type, and chi angles for residues that pass the test
		//Easiest way is using map of LigandCoordInfo --stores local coordinates, res name, coordinating atom index, ligand atom index, and chi angles
		std::set< local::HbondInfo > backside_hbond_res;
		//The graph doesn't know who its lower neighbors are, but it knows its upper neighbors
		std::set< core::Size > neighbor_res;
		//If the resnum is his_resnum, then add all the neighbors
		for ( core::conformation::PointGraph::UpperEdgeListConstIter iter = pose_point_graph->get_vertex( his_resnum ).const_upper_edge_list_begin();
				iter != pose_point_graph->get_vertex( his_resnum ).const_upper_edge_list_end();
				++iter ) {
			neighbor_res.insert( iter->upper_vertex() );
		}
		//If the resnum is less than his_resnum (his_resnum will be the upper vertex), then add the resnum
		utility::graph::UEVertex< core::conformation::PointGraphVertexData, core::conformation::PointGraphEdgeData > current_lower_res;
		for ( core::Size i = 1; i < his_resnum; ++i ) {
			current_lower_res = pose_point_graph->get_vertex( i );
			if ( current_lower_res.edge_exists( his_resnum ) ) {
				neighbor_res.insert( i );
			}
		}
		TR.Debug << "Found neighbors." << std::endl;
		//Iterate over these residues
		for ( core::Size possible_binding_partner: neighbor_res ) {
			TR.Debug << "Begin inner loop. Residue " << possible_binding_partner << std::endl;
			if ( off_limits_residues_.count( possible_binding_partner ) != 0 ) {
				continue;
			}

			//We can check for backbone hydrogen bonds first. Not counting on them existing, but it's possible
			//No need to check for clashes since we're only looking at the backbone
			core::scoring::EnergyMap bb_emap;
			hbond_->residue_pair_energy( pose.residue( his_resnum ), pose.residue( possible_binding_partner ), pose, *scorefxn_, bb_emap );
			if ( bb_emap[ core::scoring::hbond_bb_sc ] != 0 ) {
				TR << "Backbone-sidechain hydrogen bond detected!" << std::endl;
				backside_hbond_res.insert( local::HbondInfo( possible_binding_partner, "", utility::fixedsizearray1< core::Real, core::pack::dunbrack::DUNBRACK_MAX_SCTOR >(), bb_emap[ core::scoring::hbond_bb_sc ] ) );
			}
			if ( pose.residue_type( possible_binding_partner ).name3() == "PRO" || pose.residue_type( possible_binding_partner ).name3() == "GLY" ) {
				continue; //Don't mutate prolines or glycines
			}
			utility::fixedsizearray1< core::Real, 5 > bb_angles;
			bb_angles[ 1 ] = pose.phi( possible_binding_partner );
			bb_angles[ 2 ] = pose.psi( possible_binding_partner );

			check_for_hbonds( pose, his_resnum, possible_binding_partner, bb_angles, "SER", backside_hbond_res );
			check_for_hbonds( pose, his_resnum, possible_binding_partner, bb_angles, "THR", backside_hbond_res );

			if ( std::find( not_core.begin(), not_core.end(), possible_binding_partner ) != not_core.end() ) {
				check_for_hbonds( pose, his_resnum, possible_binding_partner, bb_angles, "ASN", backside_hbond_res );
			}
		}//end iterate over 2nd shell
		//Now we've found all the possible partners
		if ( backside_hbond_res.size() == 0 ) { //If we didn't find any
			TR << "Did not find hydrogen bond partner for residue " << his_resnum << std::endl;
			continue;
		}
		//Otherwise we'll need to pick the best one
		local::HbondInfo best_hbond( 0, "", utility::fixedsizearray1< core::Real, core::pack::dunbrack::DUNBRACK_MAX_SCTOR >(), 0 );
		core::Real lowest_score = 10;
		for (  local::HbondInfo current_hbond: backside_hbond_res ) {
			if ( current_hbond.hbond_score < lowest_score ) {
				best_hbond = current_hbond;
				lowest_score = best_hbond.hbond_score;
			}
		}
		//Now mutate the residue
		if ( best_hbond.residue_type != "" ) {
			mutate_residue( pose, best_hbond.resnum, best_hbond.residue_type, true );
			//Now set the correct chi angles
			for ( core::Size chi_i = 1; chi_i <= pose.residue( best_hbond.resnum ).nchi(); ++chi_i ) {
				pose.set_chi( chi_i, best_hbond.resnum, best_hbond.chi_angles[ chi_i ] );
			}
			//Now add the residue to off_limits_residues
			//The residue only needs to be off limits if we mutated it
			off_limits_residues_.insert( best_hbond.resnum );
		}
		total_hbond_score += lowest_score;
		//Set this his as cleared
		TR << "Found hydrogen bond partner for histidine " << his_resnum << std::endl;
		++num_his_cleared;
		++hbonds_found;
	}//end iterate over his
	if ( num_his_cleared < 3 ) {
		//Set job status to fail do not retry
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		return;
	}
	//otherwise we'll want to output the pose
	TR << "Successfully completed binding site hydrogen bonds!" << std::endl;
	TR << "RESNUMS: ";
	for ( core::Size resn: off_limits_residues_ ) {
		TR << resn << ",";
	}
	TR << std::endl;
	//Now add the average Hbond score to the pose's score file
	core::Real avg_hbond_score = total_hbond_score / hbonds_found;
	core::pose::setPoseExtraScore( pose, "avg_hbond_score", avg_hbond_score );
	core::pose::setPoseExtraScore( pose, "num_hbonds", hbonds_found );
	//Things to consider:
	//We need to make sure the new Hbond residues don't clash with each other, the binding site residues, or the zinc
	//Don't allow asparagine at buried positions to reduce likelihood of BUNSAT
}

void
BacksideHbondFinderMover::check_for_hbonds( core::pose::Pose & pose, core::Size his_resnum, core::Size test_resnum, utility::fixedsizearray1< core::Real, 5 > const & bb_angles, std::string res_to_mutate, std::set< local::HbondInfo > & backside_hbond_res ){
	//First mutate the residue
	mutate_residue( pose, test_resnum, res_to_mutate, true );
	//Get the rotamer library for this residue
	core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( test_resnum ).type() );
	utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamers = std::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotlib )->get_all_rotamer_samples( bb_angles );
	for ( core::pack::dunbrack::DunbrackRotamerSampleData rotamer: rotamers ) {
		if ( rotamer.probability() < 0.05 ) {
			continue;
		}
		//Set the chi angles
		for ( core::Size chi_i = 1; chi_i <= rotamer.nchi(); ++chi_i ) {
			pose.set_chi( chi_i, test_resnum, rotamer.chi_mean()[ chi_i ] );
		}

		//Check for clashes with anything in off_limits_residues
		//What would be an appropriate energy cutoff for determining if there is a clash? We don't want to be too strict since minimization may fix minor issues
		//Setup for scoring
		core::scoring::EnergyMap emap;
		bool clash = false;
		core::Real current_rep = 0;
		for ( core::Size resn: off_limits_residues_ ) {
			current_rep = emap[ core::scoring::fa_rep ];
			//Score
			if ( scorefxn_->energy_method_options().analytic_etable_evaluation() ) {
				std::dynamic_pointer_cast< core::scoring::etable::AnalyticEtableEnergy >( vdw_energy_ )->residue_pair_energy( pose.residue( resn ), pose.residue( test_resnum ), pose, *scorefxn_, emap );
			} else {
				std::dynamic_pointer_cast< core::scoring::etable::TableLookupEtableEnergy >( vdw_energy_ )->residue_pair_energy( pose.residue( resn ), pose.residue( test_resnum ), pose, *scorefxn_, emap );
			}
			if ( emap[ core::scoring::fa_rep ] - current_rep > 5 ) {
				clash = true;
			}
		}
		if ( clash ) {
			TR << "Clash detected!" << std::endl;
			continue;
		}
		hbond_->residue_pair_energy( pose.residue( his_resnum ), pose.residue( test_resnum ), pose, *scorefxn_, emap );
		//If there is a sc-sc hbond, then that term will be nonzero. If so, add it to backside_hbond_res
		if ( emap[ core::scoring::hbond_sc ] != 0 ) {
			backside_hbond_res.insert( local::HbondInfo( test_resnum, res_to_mutate, rotamer.chi_mean(), emap[ core::scoring::hbond_sc ] ) );
		}
	}//end iterate over rotamers
} //end function



int main( int argc, char* argv[] ){
	try{
		basic::options::option.add( local::atom_sasa_cutoff, "Minimum solvent-accessible surface area for an atom to be considered exposed (in square Angstroms)" ).def( 1.4 );
		basic::options::option.add( local::neighbor_distance_cutoff, "Distance cutoff for neighbor detection (in Angstroms)" ).def( 10.0 );
		devel::init( argc, argv );
		BacksideHbondFinderMoverOP bsh_finder( new BacksideHbondFinderMover );
		protocols::jd2::JobDistributor::get_instance()->go( bsh_finder );
		return 0;
	}
catch ( utility::excn::Exception const &e ){
	std::cout << "Caught exception " << e.msg() << std::endl;
	return -1;
}
}
