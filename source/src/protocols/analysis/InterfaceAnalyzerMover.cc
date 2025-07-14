// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/analysis/InterfaceAnalyzerMover.cc
/// @brief InterfaceAnalyzerMover implementation - class for in-depth interface quality analysis
/// @author Steven Lewis, Bryan Der, Ben Stranges, Jared Adolf-Bryfogle

// Unit Headers
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMoverCreator.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/dssp/Dssp.hh>

//#include <core/scoring/Energies.hh>
//#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/util.hh>

#include <protocols/jd2/util.hh>

#include <protocols/toolbox/rigid_body/util.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator2.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/CalcInterNeighborGroup.hh>

#include <core/scoring/packstat/compute_sasa.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <core/types.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
//#include <utility/exit.hh>
#include <utility>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <sstream>
#include <set>
#include <string>

#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask


static basic::Tracer TR( "protocols.analysis.InterfaceAnalyzerMover" );
static basic::Tracer TRinterface( "protocols.analysis.InterfaceAnalyzerMover.interface_selection" );
static basic::Tracer TRhbonds( "protocols.analysis.InterfaceAnalyzerMover.missing_hbonds" );


///stupid helper function needed because ternary operator does not allow variable return types
std::ostream & which_ostream( std::ostream & ost, std::ostream & oss, bool const tracer){
	if ( tracer ) return ost;
	return oss;
}

namespace protocols {
namespace analysis {

using namespace protocols::moves;
using namespace core;
using namespace utility;




///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////InterfaceAnalyzerMover///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

InterfaceAnalyzerMover::InterfaceAnalyzerMover(
	core::Size interface_jump,
	bool const tracer,
	core::scoring::ScoreFunctionCOP sf,
	bool compute_packstat,
	bool pack_input,
	bool pack_separated,
	bool use_jobname,
	bool detect_disulfide_in_separated_pose

) :
	Mover(),
	interface_jump_(interface_jump),
	tracer_(tracer),
	compute_packstat_(compute_packstat),
	explicit_constructor_(false),
	pack_input_(pack_input),
	pack_separated_(pack_separated),
	use_jobname_(use_jobname),
	detect_disulfide_in_separated_pose_(detect_disulfide_in_separated_pose),
	included_nres_(0)
	//hbond_exposure_ratio_(0),
	//total_hb_sasa_(0),
{
	protocols::moves::Mover::type( "InterfaceAnalyzer" );
	if ( sf ) { sf_ = sf->clone();}
	set_defaults();
}

//Alternate constructor for multichain poses
//Takes a set of ints for the chain nums
InterfaceAnalyzerMover::InterfaceAnalyzerMover(
	std::set<int> const & fixed_chains,
	bool const tracer,
	core::scoring::ScoreFunctionCOP sf,
	bool compute_packstat,
	bool pack_input,
	bool pack_separated,
	bool use_jobname,
	bool detect_disulfide_in_separated_pose

) :
	Mover(),
	interface_jump_(1),
	fixed_chains_(fixed_chains),
	tracer_(tracer),
	compute_packstat_(compute_packstat),
	explicit_constructor_(true),
	pack_input_(pack_input),
	pack_separated_(pack_separated),
	use_jobname_(use_jobname),
	detect_disulfide_in_separated_pose_(detect_disulfide_in_separated_pose),
	included_nres_(0)
	//hbond_exposure_ratio_(0),
	//total_hb_sasa_(0),


{
	protocols::moves::Mover::type( "InterfaceAnalyzer" );
	if ( sf ) { sf_ = sf->clone();}
	set_defaults();
}

InterfaceAnalyzerMover::InterfaceAnalyzerMover(
	std::string dock_chains,
	const bool tracer,
	core::scoring::ScoreFunctionCOP sf,
	bool compute_packstat,
	bool pack_input,
	bool pack_separated,
	bool use_jobname,
	bool detect_disulfide_in_separated_pose

):
	Mover(),
	interface_jump_(1),
	tracer_(tracer),
	compute_packstat_(compute_packstat),
	explicit_constructor_(true),
	pack_input_(pack_input),
	pack_separated_(pack_separated),
	use_jobname_(use_jobname),
	detect_disulfide_in_separated_pose_(detect_disulfide_in_separated_pose),
	included_nres_(0)
	//hbond_exposure_ratio_(0),
	//total_hb_sasa_(0),

{
	protocols::moves::Mover::type( "InterfaceAnalyzer" );
	if ( sf ) { sf_ = sf->clone();}
	set_defaults();
	dock_chains_ = dock_chains;
}

InterfaceAnalyzerMover::~InterfaceAnalyzerMover() = default;

void
InterfaceAnalyzerMover::set_defaults() {
	calcs_ready_ = false;
	compute_interface_sc_ = true;
	compute_separated_sasa_ = true;
	compute_interface_energy_ = true;
	calc_hbond_sasaE_ = true;
	compute_interface_delta_hbond_unsat_ = true;
	skip_reporting_ = false;
	use_resfile_ = false;
	use_centroid_ = false;

}


void
InterfaceAnalyzerMover::init_per_residue_data(core::pose::Pose const & pose) {

	per_residue_data_ = PerResidueInterfaceData(); //Reinit to clear any data
	per_residue_data_.regional_avg_per_residue_SASA_int.resize(3,0);
	per_residue_data_.regional_avg_per_residue_SASA_sep.resize(3,0);
	per_residue_data_.regional_avg_per_residue_energy_int.resize(3, 0);

	per_residue_data_.regional_avg_per_residue_dG.resize(3, 0);
	per_residue_data_.regional_avg_per_residue_dSASA.resize(3, 0);
	per_residue_data_.regional_avg_per_residue_energy_sep.resize(3, 0);

	per_residue_data_.interface_residues.resize(pose.size(), false);
	per_residue_data_.complexed_energy.resize(pose.size(), 0.0);
	per_residue_data_.dG.resize(pose.size(), 0.0);
	per_residue_data_.dSASA.resize(pose.size(), 0.0);
	per_residue_data_.dSASA_sc.resize(pose.size(), 0.0);
	per_residue_data_.dhSASA.resize(pose.size(), 0.0);
	per_residue_data_.dhSASA_sc.resize(pose.size(), 0.0);
	per_residue_data_.dhSASA_rel_by_charge.resize(pose.size(), 0.0);

	per_residue_data_.SASA.resize(pose.size(), 0.0);
	per_residue_data_.dSASA_fraction.resize(pose.size(), 0.0);
	per_residue_data_.separated_energy.resize(pose.size(), 0.0);
	per_residue_data_.separated_sasa.resize(pose.size(), 0.0);
}

void
InterfaceAnalyzerMover::init_data(core::pose::Pose const & pose) {

	data_ = InterfaceData(); //Reinit to clear any data
	data_.sc_value = 0;

	data_.dSASA.resize(3, 0);
	data_.dSASA_sc.resize(3, 0);

	data_.dhSASA.resize(3, 0);
	data_.dhSASA_sc.resize(3, 0);

	data_.dhSASA_rel_by_charge.resize(3, 0);

	data_.complexed_SASA = 0;
	data_.separated_SASA = 0;

	data_.dG.resize(3, 0);
	data_.gly_dG = 0;
	data_.centroid_dG = 0;
	data_.dG_dSASA_ratio = 0;
	data_.crossterm_interface_energy = 0;
	data_.crossterm_interface_energy_dSASA_ratio = 0;
	data_.complex_total_energy.resize(3, 0);
	data_.separated_total_energy.resize(3, 0);
	data_.total_hb_E = 0;


	data_.packstat = 0;

	data_.delta_unsat_hbonds = 0;
	data_.interface_hbonds = 0;
	data_.interface_nres.resize(3, 0);
	data_.complexed_interface_score.resize(3, 0);
	data_.separated_interface_score.resize(3, 0);
	data_.aromatic_nres.resize(3, 0);
	data_.aromatic_dSASA_fraction.resize(3, 0);
	data_.aromatic_dG_fraction.resize(3, 0);

	data_.ss_sheet_nres.resize(3, 0);
	data_.ss_loop_nres.resize(3, 0);
	data_.ss_helix_nres.resize(3, 0);

	data_.interface_to_surface_fraction.resize(3, 0);
	data_.interface_residues.resize(3, vector1< bool >(pose.size(), false));

	data_.delta_metrics.resize( 0, 0 );

}

void
InterfaceAnalyzerMover::init_on_new_input(const core::pose::Pose & pose){
	//Reinit structs to make sure they are clear.

	if ( !fixed_chain_strings_.empty() ) {
		fixed_chains_.clear();

		////parse the fixed chains strings to figure out pose chain nums
		for ( std::string const & fixed_chain: fixed_chain_strings_ ) {
			for ( core::Size i = 1; i<=pose.size(); ++i ) {
				if ( pose.pdb_info()->chain( i ) == fixed_chain ) {
					fixed_chains_.insert( pose.chain( i ) );
					break; //once we know something about the chain we can skip - we just need the chain id
				}
			}
		}
	}

	init_data(pose);
	init_per_residue_data(pose);

	interface_set_.clear();
	upstream_chains_.insert(0);
	downstream_chains_.insert(0);
	score_data_.clear();
	included_nres_ = 0;
	include_residue_.clear();
	include_residue_.resize(pose.size(), true);  //By default we don't ignore any residues for whole - pose calculations.

}

void
InterfaceAnalyzerMover::setup_scorefxn() {
	//Taken from per_residue_energies.cc
	// set up ScoreFunction, make certain that hydrogen bonding energies
	// are kept in the EnergyGraph.

	if ( !sf_ ) {
		TR << "NULL scorefunction. Initialize from cmd line." << std::endl;
		sf_ = core::scoring::get_score_function();
	}
	core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( sf_->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	sf_->set_energy_method_options( *emopts );

}

void InterfaceAnalyzerMover::apply_const( core::pose::Pose const & pose){

	using namespace core;
	init_on_new_input(pose);
	core::pose::Pose complexed_pose( pose );

	//If we've specified a ligand chain, then every other chain is the fixed chain
	if ( ligand_chain_.size() != 0 ) {
		// Branched carbohydrate ligands generally have the same chain (ex. A),
		// but will have different chain ids in Rosetta (ex. 1,2)
		// Using all chain ids associated with a single chain
		// will avoid problems with branching molecules
		utility::vector1< core::Size > ligand_chain_ids =
			core::pose::get_chain_ids_from_chain( ligand_chain_, pose );
		for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
			if ( !ligand_chain_ids.has_value(i) ) {
				//char chain = core::pose::get_chain_from_chain_id(i, pose);
				fixed_chains_.insert( static_cast<int>(i) );
			}
		}
	}

	core::pose::Pose separated_pose;

	//check for multichain poses
	if ( explicit_constructor_ ) {
		//fix the foldtree to reflect the fixed chains we want
		TR << "Using explicit constructor" << std::endl;

		if ( dock_chains_.size() != 0 ) {
			setup_for_dock_chains( complexed_pose, dock_chains_ );
		}

	} else if ( pose.conformation().num_chains() > 2 ) {
		//if multi chains with wrong constructor, work but print a warning
		TR.Warning << "more than 2 chains present w/o using the explicit constructor!  Values might be over the wrong jump." << std::endl;
	} else if ( pose.conformation().num_chains() < 2 ) {
		utility_exit_with_message_status( "InterfaceAnalyzerMover: pose has only one chain, aborting analysis \n", 1 );
		//remaining code works for actual interfaces
	} else {
		TR << "Using normal constructor" << std::endl;
	}

	set_pose_info( complexed_pose );
	//register calculators here if need be
	if ( !calcs_ready_ ) {
		register_calculators();
		calcs_ready_ = true;
	}
	make_interface_set(complexed_pose);
	//If there are no residues detected at the interface, don't bother with anything else. Report everything as zero and return.
	if ( interface_set_.empty() ) {
		if ( ! skip_reporting_ ) {
			setup_score_data();
			report_data();
		}
		TR << "No Interface detected.  Skipping calculations" << std::endl;
		return;
	}

	setup_scorefxn();

	//init compexed and separated pose objects
	if ( explicit_constructor_ ) {
		//Don't use the foldtree to determine which residues to move.
		separated_pose = make_separated_pose(complexed_pose, downstream_chains_);
	} else {
		separated_pose = make_separated_pose(complexed_pose, interface_jump_);
	}

	//Redetect disulfides after separation to cleave any that are in the complex pose between the chains.
	//Causes problems with clashing disulfides during RAbD
	if ( detect_disulfide_in_separated_pose_ == false ) {
		TR << "NOT detecting disulfides in the separated pose." << std::endl;
	} else {
		TR << "Detecting disulfides in the separated pose." << std::endl;
		separated_pose.conformation().detect_disulfides();
	}

	//actual computation here
	if ( compute_separated_sasa_ ) compute_separated_sasa( complexed_pose, separated_pose );
	if ( compute_interface_energy_ ) compute_interface_energy( complexed_pose, separated_pose );
	//mut_to_gly(complexed_pose, separated_pose );  this didn't help
	if ( use_centroid_ ) calc_centroid_dG( complexed_pose, separated_pose );
	if ( calc_hbond_sasaE_ ) calc_hbond_sasaE( complexed_pose ); //get hbond E at interface
	if ( compute_packstat_ ) compute_interface_packstat( complexed_pose );
	//find the change in interface hbonds
	if ( compute_interface_delta_hbond_unsat_ ) compute_interface_delta_hbond_unsat( complexed_pose, separated_pose );
	//find the shape complementarity stats for the interface
	if ( compute_interface_sc_ ) compute_interface_sc(interface_jump_, complexed_pose);
	// calculate difference for a RealMetric in complexed vs separated pose
	if ( !delta_metrics_.empty() ) compute_delta_metric( complexed_pose, separated_pose );

	setup_score_data();

	if ( !skip_reporting_ ) {
		//always will fill a selection option to get the pymol selection
		print_pymol_selection_of_interface_residues( complexed_pose, interface_set_);

		//report to tracer or job all this cool stuff we calculated
		report_data();
	}

	return;
}

/// @details InterfaceAnalyzerMover computes various interface statistics and makes them available through getters
void InterfaceAnalyzerMover::apply( core::pose::Pose & pose )
{
	apply_const( pose );
	add_score_info_to_pose( pose );
}


void InterfaceAnalyzerMover::set_pose_info( core::pose::Pose const & pose ) {
	if ( use_jobname_ ) {
		posename_base_ = protocols::jd2::current_output_name();
	} else {
		posename_ = pose.pdb_info()->name();
		posename_base_ = posename_.base();
	}
	// Intended to be used for ONLY two chain constructor.
	if ( ! explicit_constructor_ ) {
		chain1_ = pose.residue(pose.fold_tree().upstream_jump_residue(interface_jump_)).chain();
		chain2_ = pose.residue(pose.fold_tree().downstream_jump_residue(interface_jump_)).chain();
	}
}

void
InterfaceAnalyzerMover::setup_for_dock_chains( core::pose::Pose & pose, std::string dock_chains){
	TR << "Using interface constructor" <<std::endl;
	if ( ! dock_chains.find('_') ) {
		utility_exit_with_message("Unrecognized interface: "+dock_chains+" must have side1 and side2, ex: LH_A or L_H to calculate interface data");
	}

	fixed_chains_.clear();
	fixed_chain_strings_.clear();
	vector1< std::string > chainsSP = utility::string_split( dock_chains_, '_' );
	if ( pose.conformation().num_chains() == ( chainsSP[ 1 ].length() + chainsSP[ 2 ].length() ) ) {
		for ( core::Size i = 1; i <= chainsSP[ 1 ].length(); ++i ) {
			// Setup fixed chains - and let Bens multichain code do its thing
			// Branched carbohydrate ligands or protein-conjugated glycans
			// generally have the same chain (ex. A), but
			// will have different chain ids in Rosetta (ex. 1,2)
			// Using all chain ids associated with a single chain will
			// avoid problems with branching molecules
			for ( core::Size chain_id : core::pose::get_chain_ids_from_chain( std::string{chainsSP[ 1 ].at( i - 1 )}, pose ) ) {
				fixed_chains_.insert( chain_id );
			}
		}
		return;
	} else {

		//Here, we find the chains that we want to ignore in the Mover.
		std::set< int > temp_fixed_chains;
		for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
			std::string chain = core::pose::get_chain_from_chain_id( i, pose );
			if ( chain.size() == 1 && dock_chains.find(chain) !=  std::string::npos ) {
				temp_fixed_chains.insert( i );
			} else {
				fixed_chains_.insert( i );
				ignored_chains_.insert( i );
				TR << "Ignoring chain: " <<
					core::pose::get_chain_from_chain_id( i, pose ) <<
					" (" << chain << ")" << std::endl;
			}
		}

		//Move the ignored chains far away so they arn't part of the calculations.

		//Large number to make sure These castaways will not be near mobile chains.  Axis would be find as well here, but this makes them at least 2000 A away.
		core::pose::Pose sep_pose = make_separated_pose( pose, ignored_chains_, 3000 );

		//Debugging
		//sep_pose.dump_pdb("ignored_chain_sep.pdb");

		//Give the pose a normal foldtree.  Copy data into current pose.
		core::pose::set_reasonable_fold_tree( sep_pose );

		pose = sep_pose;

		//Setup fixedchains, keeping the ignored chains as fixed.
		for ( core::Size i = 1; i <= chainsSP[ 1 ].length(); ++i ) {
			// Branched carbohydrate ligands or protein-conjugated glycans generally
			// have the same chain (ex. A), but will have different chain ids in Rosetta (ex. 1,2)
			// Using all chain ids associated with a single chain will
			// avoid problems with branching molecules
			for ( core::Size chain_id : core::pose::get_chain_ids_from_chain( std::string{chainsSP[ 1 ].at( i - 1 )}, pose ) ) {
				fixed_chains_.insert( chain_id );
			}
			return;
		}
		return;
	}
}


/// @details makes the interface sets for either constructor
//always will make a new interface set, make sure this is what you want
void InterfaceAnalyzerMover::make_interface_set( core::pose::Pose & pose ){

	upstream_chains_.clear();
	downstream_chains_.clear();
	if ( explicit_constructor_ ) {
		TR << "Making an interface set with fixed chains" << std::endl;
		make_multichain_interface_set( pose, fixed_chains_ );
	} else {
		//std::set< core::Size > interface_set ( get_interface_set() );
		TR << "Making an interface set with default calculator." << std::endl;
		basic::MetricValue< std::set< core::Size > > mv_interface_set;
		pose.metric( InterfaceNeighborDefinition_, "interface_residues", mv_interface_set );
		interface_set_ = mv_interface_set.value();
		upstream_chains_.insert( chain1_ );
		downstream_chains_.insert( chain2_ );
	}
	if ( interface_set_.empty() ) {
		TR << "NO INTERFACE FOR: " << posename_base_ << std::endl;
		return;
	}

	//Transform interface_set_ to vector1<bool> to match all the per residue data and regular data.
	data_.interface_nres[ total ] = interface_set_.size();
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		auto it =  interface_set_.find( i );
		if ( it != interface_set_.end() ) {

			per_residue_data_.interface_residues[ i ] = true;
			data_.interface_residues[ total ][ i ] = true;

			core::Size chain = pose.residue( i ).chain();
			auto it_chain =  upstream_chains_.find( chain );
			InterfaceRegion region;

			if ( it_chain != upstream_chains_.end() ) {
				region = side1;
			} else {
				region = side2;
			}
			data_.interface_nres[ region ] += 1;
			data_.interface_residues[ region ][ i ] = true;
		} else {
			per_residue_data_.interface_residues[ i ] = false;
		}
	}
}

// sets the interface set for a multichain pose
// uses everything in the fixed chains as one side, everything else is the other side
void InterfaceAnalyzerMover::make_multichain_interface_set( core::pose::Pose & pose,
	std::set<int> & fixed_chains){
	using namespace core;
	using namespace utility;
	using namespace protocols::toolbox;

	std::set<core::Size> fixed_side_residues, other_side_residues;
	std::set<core::Size> fixed_chain_nums;
	std::set<core::Size> other_chain_nums;


	//itterate over all residues determine what part of the interface they are
	//also select what chain(s) are upstream and downstream of the interface
	for ( core::Size ii = 1; ii<= pose.size(); ++ii ) {
		auto it_chain = ignored_chains_.find( pose.chain( ii ) );
		if ( it_chain != ignored_chains_.end() ) {
			include_residue_[ii] = false; // Ignore this residue in total energy/Sasa calculations
		}
		if ( fixed_chains.count( pose.chain( ii ) ) ) {
			fixed_side_residues.insert( ii );
			fixed_chain_nums.insert( pose.chain( ii ) );
		} else {
			other_side_residues.insert( ii );
			other_chain_nums.insert( pose.chain( ii ) );
		}
	}
	//now assign the correct chains
	upstream_chains_ = fixed_chain_nums;
	downstream_chains_ = other_chain_nums;
	//debugging
	//TR << "Fixed residues: " << fixed_side_residues.size() << "   Other side residues: " << other_side_residues.size() << std::endl;

	//prep a vector of a pair of these residue sets for Steven's calculator
	std::pair< std::set< core::Size >, std::set< core::Size > > side_pairs;
	side_pairs.first = fixed_side_residues;
	side_pairs.second = other_side_residues;
	InterfaceAnalyzerMover::group_set interface_definition_vector;
	interface_definition_vector.push_back( side_pairs );
	chain_groups_ = interface_definition_vector ;

	CalcInterNeighborGroup calc = CalcInterNeighborGroup(chain_groups_, basic::options::option[ basic::options::OptionKeys::pose_metrics::inter_group_neighbors_cutoff ] );
	calc.compute(pose);
	interface_set_ = calc.neighbors();

	//debugging
	TR << "Interface set residues total: " << interface_set_.size() << std::endl;

}//end make_multichain_interface_set


/// @details reorder the fold tree to allow multichain interfaces to be evaluated returns the new chain for the jump
core::Size InterfaceAnalyzerMover::reorder_foldtree_find_jump( core::pose::Pose & pose ,
	std::set< int > & fixed_chains){
	using namespace core;
	using core::kinematics::Edge;

	TR << "Reordering FoldTree to find Jump for multichain interface" << std::endl;
	core::Size mobile_jump( 1 );  //setup, will be changed later
	if ( fixed_chains.empty() ) {
		utility_exit_with_message_status( "Can't find fixed chains.  Exiting...\n", 1 );
	}
	//now get info about chains
	core::Size numchains ( pose.conformation().num_chains() );

	//Ligand chain last chain in pose, no need to reorder.
	if ( ligand_chain_.size() != 0 ) {
		// Branched carbohydrate ligands or protein-conjugated glycans generally
		// have the same chain (ex. A), but will have different chain ids in Rosetta (ex. 1,2)
		// Using all chain ids associated with a single chain will
		// avoid problems with branching molecules
		for ( core::Size ligand_chain_id : core::pose::get_chain_ids_from_chain( ligand_chain_, pose ) ) {
			if ( ligand_chain_id == numchains ) {
				core::Size newjump = core::pose::get_jump_id_from_chain_id( ligand_chain_id, pose );
				return newjump;
			}
		}
	}

	utility::vector1<core::Size> chain_starts;
	utility::vector1<core::Size> chain_ends;
	for ( core::Size ii = 1; ii <= numchains; ++ii ) {
		chain_starts.push_back( pose.conformation().chain_begin( ii ) );
		chain_ends.push_back( pose.conformation().chain_end( ii ) );
	}

	//find a non mobile chain
	core::Size anchor_chain (1); //switch later if needed
	for ( core::Size ii = 1; ii<= numchains; ++ii ) {
		if ( !fixed_chains.count( ii ) ) {
			anchor_chain = ii;
			break;
		}
	}

	//now build our fold tree
	//this will move the mobile chains together
	kinematics::FoldTree foldtree ( pose.fold_tree() );
	core::Size this_jump( 1 );
	foldtree.clear();
	bool previous_mobile (false);
	for ( core::Size ii = 1; ii <= numchains; ++ii ) {
		foldtree.add_edge( Edge( chain_starts[ ii ], chain_ends[ ii ], kinematics::Edge::PEPTIDE ) );
		if ( ii != anchor_chain ) {
			if ( fixed_chains.count( ii ) && previous_mobile ) { //not in the first mobile chain
				foldtree.add_edge( Edge( chain_starts[ ii ], chain_ends[ *( fixed_chains.begin() ) ], this_jump ) );
			} else if ( fixed_chains.count( ii ) ) { //ie in the first mobile chain
				foldtree.add_edge( Edge( chain_starts[ ii ], chain_ends[ anchor_chain ], this_jump ) );
				mobile_jump = this_jump; //sets this jump to be the new mobile one
				previous_mobile=true;
			} else foldtree.add_edge( Edge( chain_starts[ ii ], chain_ends[ anchor_chain ], this_jump ) );
			++this_jump;
		}
	}
	//now set all the new values!
	foldtree.reorder( chain_starts[ anchor_chain ] );
	pose.fold_tree( foldtree );
	//TR << "Remade foldtree:\n"<< foldtree << std::endl;
	return mobile_jump;
}//end reorder foldtree


/// @brief Creates a single pose that is separated structurally at the interface
/// @details Uses the jump number to find downstream residues and move them in the Z axis in cartesian space according to step_size (A)
core::pose::Pose InterfaceAnalyzerMover::make_separated_pose( core::pose::Pose & pose, core::Size interface_jump, core::Size step_size /* 1000*/){
	using namespace core;

	( *sf_ )( pose ); //shits, giggles, and segfault prevention
	pose::Pose separated_pose( pose );

	///Get residues after the jump.
	core::Vector translation_vector(0, 0, step_size);//Translate on Z-axis
	utility::vector1< core::Size > residues_after_jump = core::kinematics::residues_downstream_of_jump( pose.fold_tree(), interface_jump);
	protocols::toolbox::rigid_body::translate(separated_pose, translation_vector, pose, residues_after_jump);

	( *sf_ )( separated_pose );
	//debugging step to make sure the jump is right
	//separated_pose.dump_pdb("IAM_test.pdb");
	return separated_pose;
}

/// @brief Creates a single pose that is separated structurally at the interface
/// @details Moves the moving_chains in cartesian space a distance according to step_size (A) in the Z axis
core::pose::Pose InterfaceAnalyzerMover::make_separated_pose( core::pose::Pose & pose, std::set<core::Size> const & chain_nums, core::Size step_size /* 1000*/){
	using namespace core;

	( *sf_ )( pose ); //shits, giggles, and segfault prevention
	pose::Pose separated_pose( pose );

	utility::vector1< core::Size > resnums;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::Size chain_num = pose.residue( i ).chain();
		if ( chain_nums.count (chain_num ) ) {
			resnums.push_back(i);
		}
	}
	core::Vector translation_vector(0, 0, step_size);//Translate on Z-axis
	protocols::toolbox::rigid_body::translate(separated_pose, translation_vector, pose, resnums);

	( *sf_ )( separated_pose );
	//debugging step to make sure the jump is right
	//separated_pose.dump_pdb("IAM_test.pdb");
	return separated_pose;
}

/// @detail Only want to register the calculators once, thus the 'if' statement in apply
void InterfaceAnalyzerMover::register_calculators()
{
	using namespace core::pose::metrics;
	using namespace protocols::simple_pose_metric_calculators;
	//determine name
	std::ostringstream interface_jump_cast;
	interface_jump_cast << interface_jump_;
	std::string const ijump(interface_jump_cast.str());


	//TR << "Running register_calculators" << std::endl;

	//this sucks but I can't figure out a way to iterate...
	Sasa_ = "Sasa_" + ijump;
	if ( CalculatorFactory::Instance().check_calculator_exists( Sasa_ ) ) {
		TR.Warning << "In InterfaceAnalyzerMover, calculator " << Sasa_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( Sasa_, utility::pointer::make_shared< core::pose::metrics::simple_calculators::SasaCalculator2 >() );
	}

	InterfaceNeighborDefinition_ = "InterfaceNeighborDefinition_" + ijump;
	if ( CalculatorFactory::Instance().check_calculator_exists( InterfaceNeighborDefinition_ ) ) {
		TR.Warning << "In InterfaceAnalyzerMover, calculator " << InterfaceNeighborDefinition_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( InterfaceNeighborDefinition_, utility::pointer::make_shared< core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator >( chain1_, chain2_ ) );
	}

	InterfaceDeltaEnergetics_ = "InterfaceDeltaEnergetics_" + ijump;
	if ( CalculatorFactory::Instance().check_calculator_exists( InterfaceDeltaEnergetics_ ) ) {
		TR.Warning << "In InterfaceAnalyzerMover, calculator " << InterfaceDeltaEnergetics_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( InterfaceDeltaEnergetics_, utility::pointer::make_shared< core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator >( InterfaceNeighborDefinition_ ) );
	}

	NumberHBonds_ = "NumberHBonds_" + ijump;
	if ( CalculatorFactory::Instance().check_calculator_exists( NumberHBonds_ ) ) {
		TR.Warning << "In InterfaceAnalyzerMover, calculator " << NumberHBonds_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( NumberHBonds_, utility::pointer::make_shared< NumberHBondsCalculator >() );
	}

	BuriedUnsatisfiedPolars_ = "BuriedUnsatisfiedPolars_" + ijump;
	if ( CalculatorFactory::Instance().check_calculator_exists( BuriedUnsatisfiedPolars_ ) ) {
		TR.Warning << "In InterfaceAnalyzerMover, calculator " << BuriedUnsatisfiedPolars_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( BuriedUnsatisfiedPolars_, utility::pointer::make_shared< BuriedUnsatisfiedPolarsCalculator >( Sasa_, NumberHBonds_) );
	}
	return;
}//register_calculators


/// @details actual function to separate the chains based on the chosen jump and score
void InterfaceAnalyzerMover::score_separated_chains( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose ) {
	//using namespace core::pack::task;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;

	// TR<< "Initial complex_score: " << complex_score << "  sep_score: "
	//   << separated_score << std::endl;

	//setup mover to pack interface
	//core::Size ndruns = option[packing::ndruns];

	protocols::minimization_packing::PackRotamersMoverOP repacker( new protocols::minimization_packing::PackRotamersMover( sf_ ) );
	repacker->nloop( pack_rounds_ ); //JAB - increasing as an nloop of 1 is not enough.

	if ( pack_input_ ) {
		core::pack::task::PackerTaskOP task = setup_task( complexed_pose );
		repacker->task( task );
		repacker->apply( complexed_pose );
	}
	if ( pack_separated_ ) {
		//This is to correctly deal with out-of-date disulfides repacking in
		// the separated pose. Cannot generate task and apply to both together
		// and separated pose for this reason
		core::pack::task::PackerTaskOP task = setup_task( separated_pose );
		repacker->task( task );
		repacker->apply( separated_pose );
	}

	( *sf_ )( complexed_pose );
	( *sf_ )( separated_pose );
	//separated_pose.dump_pdb("IAM_test.pdb");
	for ( core::Size i = 1; i <= complexed_pose.size(); ++i ) {
		if ( !include_residue_[ i ] ) continue;
		included_nres_ += 1;

		InterfaceRegion region;
		if ( upstream_chains_.count( complexed_pose.chain( i ) ) ) {
			region = side1;
		} else {
			region = side2;
		}
		core::Real complexed_energy = complexed_pose.energies().residue_total_energy( i );
		core::Real separated_energy = separated_pose.energies().residue_total_energy( i );


		data_.complex_total_energy[ total ] += complexed_energy;
		data_.complex_total_energy[ region ] += complexed_energy;

		data_.separated_total_energy[ total ] += separated_energy;
		data_.separated_total_energy[ region ] += separated_energy;

		per_residue_data_.complexed_energy[ i ] = complexed_energy;
		per_residue_data_.separated_energy[ i ] = separated_energy;
	}

	//debugging step to make sure the jump is right
	//copy.dump_pdb("IAM_test.pdb");
	//  TR<< "Post Packing complex_score: " << complex_score << "  sep_score: "
	//   << separated_score << std::endl;

	TR << "included_nres: " << included_nres_ << std::endl;

	for ( core::Size i = 1; i <= 3; ++i ) {
		data_.dG[ i ] = data_.complex_total_energy[ i ] - data_.separated_total_energy[ i ];
	}
}


/// @details computes the SASA by finding difference between complex and separated SASA
/// also does the same thing for hydrophobic/polar SASA
void InterfaceAnalyzerMover::compute_separated_sasa(core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose){
	using namespace core;

	TR << "Calculating dSASA" << std::endl;

	//Speedup here can come from using SasaCalc directly instead of calculator to limit accessing all the vectors.

	//Complex Pose values.
	basic::MetricValue< Real > mv_complex_sasa;

	basic::MetricValue< vector1< core::Real > > mv_complex_res_sasa;
	basic::MetricValue< vector1< core::Real > > mv_complex_res_sasa_sc;
	basic::MetricValue< vector1< core::Real > > mv_complex_res_hsasa;
	basic::MetricValue< vector1< core::Real > > mv_complex_res_hsasa_sc;
	basic::MetricValue< vector1< core::Real > > mv_complex_res_rel_hsasa;

	complexed_pose.metric( Sasa_, "total_sasa", mv_complex_sasa );
	complexed_pose.metric( Sasa_, "residue_sasa", mv_complex_res_sasa );
	complexed_pose.metric( Sasa_, "residue_sasa_sc", mv_complex_res_sasa_sc );
	complexed_pose.metric( Sasa_, "residue_hsasa", mv_complex_res_hsasa );
	complexed_pose.metric( Sasa_, "residue_hsasa_sc", mv_complex_res_hsasa_sc );
	complexed_pose.metric( Sasa_, "residue_rel_hsasa", mv_complex_res_rel_hsasa );


	//Separated Pose values.
	basic::MetricValue< core::Real > mv_sep_sasa;

	basic::MetricValue< vector1< core::Real > > mv_sep_res_sasa;
	basic::MetricValue< vector1< core::Real > > mv_sep_res_sasa_sc;
	basic::MetricValue< vector1< core::Real > > mv_sep_res_hsasa;
	basic::MetricValue< vector1< core::Real > > mv_sep_res_hsasa_sc;
	basic::MetricValue< vector1< core::Real > > mv_sep_res_rel_hsasa;

	separated_pose.metric( Sasa_, "total_sasa", mv_sep_sasa );
	separated_pose.metric( Sasa_, "residue_sasa", mv_sep_res_sasa );
	separated_pose.metric( Sasa_, "residue_sasa_sc", mv_sep_res_sasa_sc );
	separated_pose.metric( Sasa_, "residue_hsasa", mv_sep_res_hsasa );
	separated_pose.metric( Sasa_, "residue_hsasa_sc", mv_sep_res_hsasa_sc );
	separated_pose.metric( Sasa_, "residue_rel_hsasa", mv_sep_res_rel_hsasa );

	data_.complexed_SASA = mv_complex_sasa.value();
	data_.separated_SASA = mv_sep_sasa.value();


	per_residue_data_.separated_sasa = mv_sep_res_sasa.value();
	per_residue_data_.complexed_sasa = mv_complex_res_sasa.value();




	if ( compute_separated_sasa_ ) {
		TR << "Calculating per-res dSASA data" << std::endl;
		per_residue_data_.dSASA = calc_per_residue_dSASA( complexed_pose, mv_sep_res_sasa.value(), mv_complex_res_sasa.value() );
		for ( core::Size i = 1; i <= complexed_pose.size(); ++i ) {

			per_residue_data_.dSASA_sc[ i ] = calc_per_residue_dSASA_general( i, complexed_pose,
				mv_sep_res_sasa_sc.value()[ i ],
				mv_complex_res_sasa_sc.value()[ i ],
				data_.dSASA_sc );

			per_residue_data_.dhSASA[ i ] = calc_per_residue_dSASA_general( i, complexed_pose,
				mv_sep_res_hsasa.value()[ i ],
				mv_complex_res_hsasa.value()[ i ],
				data_.dhSASA );


			per_residue_data_.dhSASA_sc[ i ] = calc_per_residue_dSASA_general( i, complexed_pose,
				mv_sep_res_hsasa_sc.value()[ i ],
				mv_complex_res_hsasa_sc.value()[ i ],
				data_.dhSASA_sc );

			per_residue_data_.dhSASA_rel_by_charge[ i ] = calc_per_residue_dSASA_general( i, complexed_pose,
				mv_sep_res_rel_hsasa.value()[ i ],
				mv_complex_res_rel_hsasa.value()[ i ],
				data_.dhSASA_rel_by_charge );
		}

		calc_interface_to_surface_fraction( separated_pose, mv_sep_res_sasa.value() );
		return;
	}
}

Real
InterfaceAnalyzerMover::calc_per_residue_dSASA_general( const core::Size i, const core::pose::Pose & complexed_pose, const Real separated_sasa, const Real complexed_sasa, vector1<Real>& regional ) {

	//Real dSASA = 0.0;
	Real dSASA  = separated_sasa - complexed_sasa;
	regional[ total ] += dSASA;

	if ( per_residue_data_.interface_residues[ i ] ) {
		InterfaceRegion region;

		if ( upstream_chains_.count( complexed_pose.chain( i ) ) ) {
			region = side1;
		} else {
			region = side2;
		}
		regional[ region ] += dSASA;
	}
	return dSASA;
}

vector1< core::Real >
InterfaceAnalyzerMover::calc_per_residue_dSASA( const core::pose::Pose & complexed_pose, const vector1<core::Real> & separated_sasa,  const vector1<core::Real> & complexed_sasa ){

	vector1< core::Real > dSASA;
	if ( separated_sasa.size() != complexed_sasa.size() ) {
		utility_exit_with_message( "Separated sasa and complexed sasa must be the same size to calculate res dSASA" );
	}

	for ( core::Size i = 1; i<= separated_sasa.size(); ++i ) {

		dSASA.push_back( separated_sasa[ i ] - complexed_sasa[ i ] );
		data_.dSASA[ total ] += dSASA[ i ];

		per_residue_data_.dSASA_fraction[ i ] = dSASA[ i ] / separated_sasa[ i ];
		if ( per_residue_data_.interface_residues[ i ] ) {
			InterfaceRegion region;
			if ( upstream_chains_.count( complexed_pose.chain( i ) ) ) {
				region = side1;
			} else {
				region = side2;
			}
			per_residue_data_.regional_avg_per_residue_dSASA[ total ] += ( dSASA[ i ] );
			per_residue_data_.regional_avg_per_residue_SASA_sep[ total ] += separated_sasa[ i ];
			per_residue_data_.regional_avg_per_residue_SASA_int[ total ] += complexed_sasa[ i ];

			per_residue_data_.regional_avg_per_residue_dSASA[ region ] += ( dSASA[ i ] );
			per_residue_data_.regional_avg_per_residue_SASA_sep[ region ] += separated_sasa[ i ];
			per_residue_data_.regional_avg_per_residue_SASA_int[ region ] += complexed_sasa[ i ];
		}
	}

	for ( core::Size i = 1; i <= 3; ++i ) {
		per_residue_data_.regional_avg_per_residue_SASA_int[ i ] = per_residue_data_.regional_avg_per_residue_SASA_int[ i ] / interface_set_.size();
		per_residue_data_.regional_avg_per_residue_SASA_sep[ i ] = per_residue_data_.regional_avg_per_residue_SASA_sep[ i ] / interface_set_.size();
		per_residue_data_.regional_avg_per_residue_dSASA[ i ] = per_residue_data_.regional_avg_per_residue_dSASA[ i ] / interface_set_.size();
	}
	return dSASA;
}

vector1< core::Real >
InterfaceAnalyzerMover::calc_per_residue_dG( core::pose::Pose & complexed_pose, const vector1< core::Real > & separated_energy, const vector1< core::Real > & complexed_energy ) {
	vector1< core::Real > dG;
	if ( separated_energy.size() != complexed_energy.size() ) {
		utility_exit_with_message("Separated energy and complexed sasa must be the same size to calculate res dSASA");
	}
	for ( core::Size i = 1; i <= separated_energy.size(); ++i ) {
		dG.push_back( complexed_energy[ i ] - separated_energy[ i ] );
		if ( per_residue_data_.interface_residues[ i ] ) {
			InterfaceRegion region;
			if ( upstream_chains_.count( complexed_pose.chain( i ) ) ) {
				region = side1;
			} else {
				region = side2;
			}
			per_residue_data_.regional_avg_per_residue_dG[ total ] += dG[ i ];
			per_residue_data_.regional_avg_per_residue_energy_int[ total ] += complexed_energy[ i ];
			per_residue_data_.regional_avg_per_residue_energy_sep[ total ] += separated_energy[ i ];

			per_residue_data_.regional_avg_per_residue_dG[ region ] += dG[ i ];
			per_residue_data_.regional_avg_per_residue_energy_int[ region ] += complexed_energy[ i ];
			per_residue_data_.regional_avg_per_residue_energy_sep[ region ] += separated_energy[ i ];
		}
	}

	for ( core::Size i = 1; i<=3; ++i ) {
		per_residue_data_.regional_avg_per_residue_dG[ i ] = per_residue_data_.regional_avg_per_residue_dG[ i ] / interface_set_.size();
		per_residue_data_.regional_avg_per_residue_energy_int[ i ] = per_residue_data_.regional_avg_per_residue_energy_int[ i ] / interface_set_.size();
		per_residue_data_.regional_avg_per_residue_energy_sep[ i ] = per_residue_data_.regional_avg_per_residue_energy_sep[ i ] / interface_set_.size();
	}
	return dG;
}

/// @details computes the interface energy of the interface
void InterfaceAnalyzerMover::compute_interface_energy( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose )
{
	//separated interface energy and ratio
	score_separated_chains( complexed_pose, separated_pose );  //sets data_.dG
	calc_per_residue_and_regional_data( complexed_pose, separated_pose );

	//crossterm interface energy and ratio
	// InterfaceDeltaEnergetics only works for two chains.
	// So, instead of giving odd values, we only calculate it if we only have two chains.
	if ( complexed_pose.conformation().num_chains() == 2 ) {
		basic::MetricValue< core::Real > mv_delta_total;
		complexed_pose.metric( InterfaceDeltaEnergetics_, "weighted_total", mv_delta_total );
		data_.crossterm_interface_energy = mv_delta_total.value();
	}

	if ( data_.dSASA[ total ] > 0 ) {
		data_.crossterm_interface_energy_dSASA_ratio = data_.crossterm_interface_energy / data_.dSASA[ total ];
		data_.dG_dSASA_ratio = data_.dG[ total ] / data_.dSASA[ total ];
	}
	return;
}

void InterfaceAnalyzerMover::compute_delta_metric( core::pose::Pose const & complexed_pose, core::pose::Pose const & separated_pose ) {

	for ( auto const & metric : delta_metrics_ ) {
		core::Real complexed_value = metric->calculate(complexed_pose);
		core::Real separated_value = metric->calculate(separated_pose);
		data_.delta_metrics.emplace_back( complexed_value - separated_value );
	}

}

/// @details calculate the average energy per residue in the interface as well as other data
void InterfaceAnalyzerMover::calc_per_residue_and_regional_data( core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose ) {
	using namespace core::scoring;
	using namespace core;

	( *sf_ ) ( complexed_pose ); //segfault prevention
	( *sf_ ) ( separated_pose );
	//itterate over all residues calc total energy for each then take avg

	core::Real complexed_residue_score;
	core::Real separated_residue_score;

	core::scoring::dssp::Dssp dssp( complexed_pose );
	vector1< core::Real > aromatic_dSASA( 3, 0.0 );
	vector1< core::Real > aromatic_dG( 3, 0.0 );

	if ( data_.interface_nres[ total ] == 0 ) {
		TR << "No interface detected.  Skipping per residue data calculation." << std::endl;
		return;
	}
	for ( core::Size it : interface_set_ ) {
		complexed_residue_score = complexed_pose.energies().residue_total_energies( it )[ ScoreType( total_score ) ];
		separated_residue_score = separated_pose.energies().residue_total_energies( it )[ ScoreType( total_score ) ];

		InterfaceRegion region;
		core::Size chain = complexed_pose.residue( it ).chain();
		if ( upstream_chains_.count( chain ) ) {
			region = side1;
		} else {
			region = side2;
		}

		//Calculate side 1 and side2 interface score and other values.
		data_.complexed_interface_score[ region ] += complexed_residue_score;
		data_.complexed_interface_score[ total ] += complexed_residue_score;
		data_.separated_interface_score[ region ] += separated_residue_score;
		data_.separated_interface_score[ total ] += separated_residue_score;

		data_.dSASA[ region ] += per_residue_data_.dSASA[ it ];
		data_.dhSASA[ region ] += per_residue_data_.dhSASA[ it ];
		data_.dhSASA_rel_by_charge[ region ] += per_residue_data_.dhSASA_rel_by_charge[ it ];

		if ( complexed_pose.residue( it ).is_aromatic() ) {
			data_.aromatic_nres[ region ] += 1;
			core::Real dG = complexed_residue_score- separated_residue_score;
			aromatic_dG[ total ] += dG;
			aromatic_dG[ region ] += dG;
			if ( data_.dSASA[ total ] > 0.0 ) {
				aromatic_dSASA[ total ]+= per_residue_data_.dSASA[ it ];
				aromatic_dSASA[ region ] += per_residue_data_.dSASA[ it ];
			}
		}
		//SS
		//skip non-protein residues for SS calculation
		if ( complexed_pose.residue( it ).is_protein() ) {
			if ( dssp.get_dssp_secstruct( it ) =='E' ) {
				data_.ss_sheet_nres[ region ] += 1;
				data_.ss_sheet_nres[ total ] += 1;
			} else if ( dssp.get_dssp_secstruct( it ) == 'H' ) {
				data_.ss_helix_nres[ region ] += 1;
				data_.ss_helix_nres[ total ] += 1;
			} else {
				data_.ss_loop_nres[ region ] += 1;
				data_.ss_loop_nres[ total ] += 1;
			}
		}
	}

	per_residue_data_.dG = calc_per_residue_dG( complexed_pose, per_residue_data_.separated_energy, per_residue_data_.complexed_energy );
	per_residue_data_.regional_avg_per_residue_energy_int[ side1 ] = data_.complexed_interface_score[ side1 ] / data_.interface_nres[ side1 ];
	per_residue_data_.regional_avg_per_residue_energy_int[ side2 ] = data_.complexed_interface_score[ side2 ] / data_.interface_nres[ side2 ];
	data_.aromatic_nres[ total ] = data_.aromatic_nres[ side1 ] + data_.aromatic_nres[ side2 ];
	for ( core::Size i = 1; i <= 3; ++i ) {
		if ( data_.interface_nres[ i ] > 0 ) {
			data_.aromatic_dG_fraction[ i ] = aromatic_dG[ i ] / data_.dG[ i ];
		}
	}

	//Loop over all regions, and avoid dividing by zero
	for ( core::Size i = 1; i <= 3; ++i ) {
		if ( data_.interface_nres[ i ] > 0 ) {
			data_.aromatic_dSASA_fraction[ i ] = aromatic_dSASA[ i ] / data_.dSASA[ i ];
		}
	}
}


//@details Averaged packstat of interface residues only
void InterfaceAnalyzerMover::compute_interface_packstat( core::pose::Pose & pose )
{
	TR << "Computing interface packstats." << std::endl;

	//calculates the packstat scores for the interface set
	utility::vector1< core::Real > interface_pack_scores;
	core::Real interface_pack_score_sum( 0.0 );
	core::Size interface_num_res( 0 );
	interface_pack_scores = core::scoring::packstat::compute_residue_packing_scores( pose, basic::options::option[ basic::options::OptionKeys::packstat::oversample ]() );

	//JAB - Still bugs in packstat  - so segfault prevention!:
	if ( interface_pack_scores.size() != pose.size() ) {
		data_.packstat = 0;
		return;
	}

	for ( core::Size it : interface_set_ ) {
		interface_pack_score_sum += interface_pack_scores[ it ];
		++interface_num_res;
	}
	//fills the selection for pymol output, doesn't print anything here
	print_pymol_selection_of_packing( pose, interface_pack_scores );
	data_.packstat = interface_pack_score_sum / interface_num_res;
}

/// @details If a polar atom at the interface is also "buried unsat" in the monomer, we don't count this one
void InterfaceAnalyzerMover::compute_interface_delta_hbond_unsat( core::pose::Pose & complexed_pose,
	core::pose::Pose & separated_pose)
{

	//do the calculations
	TR << "Computing delta unsat polar residues..." << std::endl;

	basic::MetricValue< utility::vector1< core::Size > > mv_complexed_unsat_res;
	basic::MetricValue< utility::vector1< core::Size > > mv_separated_unsat_res;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_complexed_unsat_map;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_separated_unsat_map;

	complexed_pose.metric( BuriedUnsatisfiedPolars_, "residue_bur_unsat_polars", mv_complexed_unsat_res );
	separated_pose.metric( BuriedUnsatisfiedPolars_, "residue_bur_unsat_polars", mv_separated_unsat_res );
	utility::vector1< core::Size > const complexed_unsat_res( mv_complexed_unsat_res.value() );
	utility::vector1< core::Size > const separated_unsat_res( mv_separated_unsat_res.value() );

	complexed_pose.metric( BuriedUnsatisfiedPolars_, "atom_bur_unsat", mv_complexed_unsat_map );
	separated_pose.metric( BuriedUnsatisfiedPolars_, "atom_bur_unsat", mv_separated_unsat_map );
	core::id::AtomID_Map< bool > const complexed_unsat_map( mv_complexed_unsat_map.value() );
	core::id::AtomID_Map< bool > const separated_unsat_map( mv_separated_unsat_map.value() );

	//loop over the interface set and figure out what's burried/unsat
	core::Size delta_unsat_hbond_counter( 0 );
	utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector;
	for ( core::Size it : interface_set_ ) {
		//TR << "UnsatHbond res " << *it << std::endl;
		//iterate over all its atoms, check if they're in the map of missing, and print their string name
		core::chemical::ResidueType const & res( complexed_pose.residue_type( it ) );
		for ( core::Size i=1 ; i <= res.natoms(); ++i ) {
			core::id::AtomID const atid(i, it);

			//TR << "UnsatHbond atom " << i << std::endl;

			//if atom is buried unsat in complexed pose but not separated pose, count it as delta unsat
			bool unsat_complex( complexed_unsat_map.has( atid ) && complexed_unsat_map( atid ) );
			bool unsat_separated( separated_unsat_map.has( atid ) && separated_unsat_map( atid ) );

			if ( unsat_complex && !unsat_separated ) {
				//results << " " << res.atom_name(i) ;
				delta_unsat_hbond_atid_vector.push_back( atid );
				++delta_unsat_hbond_counter;
			}
		}//end for each atom
	}//end loop of interface_set_
	data_.delta_unsat_hbonds = delta_unsat_hbond_counter; //sets the total
	//print_pymol_selection_of_interface_residues( pose, interface_set_);
	//calculated here but need to call get the selection to output it
	print_pymol_selection_of_hbond_unsat( complexed_pose, delta_unsat_hbond_atid_vector );
	return;
}

void
InterfaceAnalyzerMover::calc_interface_to_surface_fraction(core::pose::Pose const & separated_pose, const vector1<core::Real> & separated_sasa){

	// Cutoff sasa value taken from LayerDesign operations.
	// A percentage buried should be the more correct way to do this.
	// But what is the maximal for each residue?

	vector1< core::Size > surface_nres( 3, 0 );

	for ( core::Size i = 1; i <= separated_sasa.size(); ++i ) {
		if ( !include_residue_[ i ] || separated_sasa[ i ] < 40.0 ) { continue; }

		core::Size chain = separated_pose.chain( i );
		InterfaceRegion region;
		if ( upstream_chains_.count( chain ) ) {
			region = side1;
		} else {
			region = side2;
		}
		surface_nres[ region ] += 1;
	}

	//Would surface_nres ever be 0? Don't think so.
	surface_nres[ total ] = surface_nres[ side1 ] + surface_nres[ side2 ];
	for ( core::Size i = 1; i <= 3; ++i ) {
		if ( surface_nres[ i ] > 0 ) {
			data_.interface_to_surface_fraction[ i ] = data_.interface_nres[ i ] / ( core::Real ) surface_nres[ i ];
		}
	}

}


/// @details calculate the hbond energy and dampen it by exposure
void InterfaceAnalyzerMover::calc_hbond_sasaE( core::pose::Pose & pose ) {
	using namespace core::scoring::hbonds;
	using namespace core;
	using namespace core::chemical;

	( *sf_ ) ( pose ); //segfault prevention
	//calculate the per atom sasa
	// an atomID map is needed for the calc_per_atom_sasa method;
	// it stores the actual calculated sasa for every atom
	//core::id::AtomID_Map< core::Real > atom_sasa;
	//core::pose::initialize_atomid_map( atom_sasa, pose, (core::Real)0.0 ); // initialize to 0.0 for "not computed"

	//utility::vector1< Real > rsd_sasa( pose.size(), 0.0 );
	//core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];

	// create an atom_subset mask such that only the heavy atoms will
	// have their sasa computed (ignore H's to make it faster)
	core::id::AtomID_Map< bool > atom_subset;
	atom_subset.clear();
	atom_subset.resize( pose.size() );
	for ( core::Size ii=1; ii <= pose.size(); ++ii ) {
		atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
		for ( core::Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
			atom_subset[ ii ][ jj ] = true;
		}
	}
	//now actually calculate the SASA for heavy atoms only
	//   core::Real total_sasa = 0.0;
	//   total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset );
	//  //some tracer output to make sure things work
	//  for ( core::Size ii=1; ii <= pose.size(); ++ii ) {
	//   core::conformation::Residue const & rsd = pose.residue( ii );
	//   TR << "residue " << rsd.name3() << pose.pdb_info()->number(ii) << " atom_sasas: [ ";
	//   for ( core::Size at=1; at <= rsd.nheavyatoms(); ++at ) {
	//    core::id::AtomID atid( at, ii );
	//    TR << utility::trim( rsd.atom_name( at ) ) << ":" << atom_sasa[ atid ] << ", ";
	//   }
	//   TR << "], total residue SASA: " << rsd_sasa[ ii ] << std::endl;
	//  }

	//setup a vector of the sasa radii
	//Real const four_pi = 4.0f * Real( numeric::constants::d::pi );
	//AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	//utility::vector1< Real > radii( atom_type_set.n_atomtypes(), 0.0 );
	//core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "SASA_RADIUS_LEGACY" );

	// TR << "radii: [ ";
	// for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
	//  radii[ ii ] = atom_type_set[ ii ].extra_parameter( SASA_RADIUS_INDEX );
	//  TR << radii[ ii ] << ", ";
	// }
	// TR << "]" << std::endl;

	//EM options for bb-bb hbond output
	TR << "Making hbond ScoreFunction with bb-bb energy method options" << std::endl;
	core::scoring::ScoreFunctionOP new_sf = scoring::get_score_function();
	scoring::methods::EnergyMethodOptions energymethodoptions( new_sf->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies( true );
	new_sf->set_energy_method_options( energymethodoptions );

	//make an HbondSet
	//figure out energy statistics
	core::scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	//fill_hbond_set( pose,
	//       false /*calc_deriv*/,
	//       hbond_set,
	//       false /*bb only*/ );
	hbond_set.setup_for_residue_pair_energies( pose, false, false );

	//itterate through all hbonds and figure out which ones are bb-bb betas
	core::Size n_crosschain_hbonds( 0 );

	//Real total_ratios( 0 );
	//total_hb_sasa_ = 0 ;
	data_.total_hb_E = 0;
	for ( core::Size ii = 1; ii <= hbond_set.nhbonds(); ++ii ) {
		core::scoring::hbonds::HBond hbond ( hbond_set.hbond( ii ) );
		core::Size don_resnum = hbond.don_res(); core::Size acc_resnum = hbond.acc_res();
		//need to take special consideration for multichain constructor
		if ( explicit_constructor_ ) {
			//TR<< "Do multichain hbonds eval..." << std::endl;
			//if they are different chains
			if ( pose.chain( don_resnum ) != pose.chain( acc_resnum ) ) {
				//if the acceptor or donor is in a fixed chain
				if ( fixed_chains_.count( pose.chain( don_resnum ) ) || fixed_chains_.count( pose.chain( acc_resnum ) ) ) {
					//if the acceptor and donor are not both fixed chains
					if ( fixed_chains_.count( pose.chain( don_resnum ) ) != fixed_chains_.count( pose.chain( acc_resnum ) ) ) {
						TR << "Found Hbond between chains: "<< pose.chain( don_resnum ) << " and " << pose.chain( acc_resnum ) << std::endl;
						//now copy the same stuff as for the normal case:
						n_crosschain_hbonds += 1;
						data_.total_hb_E += hbond.energy();
					}
				}
			}//end if different chains
		} else { //not multichain constructor
			if ( pose.chain( hbond.don_res() ) != pose.chain( hbond.acc_res() ) ) {
				n_crosschain_hbonds += 1;
				data_.total_hb_E += hbond.energy();
			} //end if chains not equal
		} //end if not multichain
	}//end loop over all hbonds
	//get the avg exposure of the hbonds
	//hbond_exposure_ratio_ = total_ratios / n_crosschain_hbonds;

	data_.hbond_E_fraction = data_.total_hb_E/data_.dG[total];
	data_.interface_hbonds = n_crosschain_hbonds;

}//end function def

// //helper function for above to calculate hbond stats based on an hbond input
// void hbond_info_calculate( core::scoring::hbonds::HBond hbond ){
//  //does nothing for now, figure out later
// }

void
InterfaceAnalyzerMover::compute_interface_sc( core::Size &, core::pose::Pose const & complexed_pose){

	core::scoring::sc::ShapeComplementarityCalculator sc_calc;
	// Split PDB into two surfaces
	for ( core::Size i = 1; i <= complexed_pose.size(); i++ ) {
		if ( upstream_chains_.count( complexed_pose.chain( i ) ) && include_residue_[ i ] ) {
			sc_calc.AddResidue( 0, complexed_pose.residue( i ) );
		} else if ( downstream_chains_.count( complexed_pose.chain( i ) ) && include_residue_[ i ] ) {
			sc_calc.AddResidue( 1, complexed_pose.residue( i ) );
		} else {
			continue;
		}
	}
	//now calculate and print results
	TR << "Computing Shape Complementarity Score..." << std::endl;
	TR << "Upstream chain(s) numbers: ";
	for ( core::Size upstream_chain : upstream_chains_ ) {
		TR << upstream_chain << ", ";
	}
	TR << std::endl;

	TR << "Downstream chain(s) numbers: ";
	for ( core::Size downstream_chain : downstream_chains_ ) {
		TR << downstream_chain << ", ";
	}
	TR << std::endl;

	//actual calculate function
	try{
		sc_calc.Calc();
		core::scoring::sc::RESULTS const results = sc_calc.GetResults();
		data_.sc_value = results.sc;
	}
catch (utility::excn::Exception& excn){
	TR << "SC calculation failed.  Setting to 0" << std::endl;
	data_.sc_value = 0;
}

}//end compute_interface_sc


/// @details  Mutate all residues to GLY rescore complex energy and separated energy
void InterfaceAnalyzerMover::mut_to_gly( core::pose::Pose complex_pose, core::pose::Pose separated_pose ) {
	using namespace core;
	//need a copy of the pose to avoid screwing up the good one
	pose::Pose copy_complex( complex_pose );
	pose::Pose copy_separate( separated_pose );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	//setupt task info
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( ( copy_complex ) ) );
	utility::vector1_bool packable( copy_complex.size(), false ); //false = nobody is packable
	utility::vector1< bool > allowed_aa( chemical::num_canonical_aas, false ); //no allowed residues
	allowed_aa[ core::chemical::aa_from_oneletter_code( 'G' ) ] = true; //allow gly only
	//allow all interface residues to be mutated to Gly
	for ( core::Size it : interface_set_ ) {
		task->nonconst_residue_task( it ).restrict_absent_canonical_aas( allowed_aa );
		packable[ it ] = true;
	}
	task->restrict_to_residues( packable );  //prevents non interface res from changing

#ifndef NDEBUG
	TR<< "GLY Packer Task: " << *(task) << std::endl;
#endif

	//apply mutations
	protocols::minimization_packing::PackRotamersMoverOP packrot_mover( new protocols::minimization_packing::PackRotamersMover( sf_ , task ) );
	packrot_mover->apply( copy_complex );
	packrot_mover->apply( copy_separate );
	data_.gly_dG =(*sf_) ( copy_complex ) - (*sf_) ( copy_separate );
}

/// @details
void InterfaceAnalyzerMover::calc_centroid_dG ( core::pose::Pose const & complex_pose, core::pose::Pose const & separated_pose ){
	if ( !use_centroid_ ) {
		data_.centroid_dG = 0;
		return;
	}
	core::pose::Pose copy_complex( complex_pose );
	core::pose::Pose copy_separated( separated_pose );
	core::util::switch_to_residue_type_set( copy_complex , core::chemical::CENTROID_t );
	core::util::switch_to_residue_type_set( copy_separated , core::chemical::CENTROID_t );
	// use score3 but turn of RG - JAB why score3?
	core::scoring::ScoreFunctionOP scorefxn  = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scorefxn->set_weight( core::scoring::rg, 0.0 );
	//Debugging:
	TR << "Centroid score of complex: " << (*scorefxn) ( copy_complex ) << std::endl;
	TR << "Centroid score of separated: " << (*scorefxn) ( copy_separated ) << std::endl;

	data_.centroid_dG = ( *scorefxn ) ( copy_complex ) - ( *scorefxn ) ( copy_separated );
}


/// @details  sets up the packer task for the interface
core::pack::task::PackerTaskOP
InterfaceAnalyzerMover::setup_task( core::pose::Pose & pose ) {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//set up the task to match this calculation
	//set up a packer task
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >() );
	//force include current to prevent wonky results
	tf->push_back( utility::pointer::make_shared< IncludeCurrent >() );
	tf->push_back( utility::pointer::make_shared< RestrictToRepacking >() );

	if ( use_resfile_ ) tf->push_back( utility::pointer::make_shared< ReadResfile >() );
	core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	task->restrict_to_residues(data_.interface_residues[total]); //Does not turn on packing.
	return task;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
InterfaceAnalyzerMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap
)
{
	sf_ = protocols::rosetta_scripts::parse_score_function( tag, datamap )->clone();

	set_pack_separated(tag->getOption< bool >( "pack_separated", false ) );
	set_use_resfile(tag->getOption< bool >( "resfile", false ) );
	set_compute_packstat(tag->getOption< bool > ( "packstat", false ) );
	set_compute_interface_sc(tag->getOption< bool > ("interface_sc", false ) );
	set_pack_input(tag->getOption< bool > ( "pack_input", false ) );
	set_use_tracer(tag->getOption< bool > ( "tracer", false ) );
	set_use_jobname(tag->getOption< bool > ( "use_jobname", false ) );


	if ( tag->hasOption( "jump" ) && tag->hasOption( "fixedchains" ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Jump and fixedchains are mutually exclusive. Use either jump or fixedchains" );
	}

	if ( ( tag->hasOption( "jump" ) || tag->hasOption( "fixedchains" ) ) && tag->hasOption( "ligandchain" ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "if you specify Jump or fixedchains you cannot also specify ligandchain" );
	}

	if ( tag->hasOption( "fixedchains" ) ) {
		set_interface_jump( 0 );
		std::string chains_string = tag->getOption<std::string>( "fixedchains" );
		fixed_chain_strings_ = utility::string_split( chains_string, ',' );
		TR << "Fixed chains are: " ;
		for ( auto const & this_chain: fixed_chain_strings_ ) {
			TR << this_chain << ", ";
		}
		TR << "these will be moved together." << std::endl;

		explicit_constructor_ = true;
	} else if ( tag->hasOption( "interface" ) ) {
		set_interface_jump( 0 );
		explicit_constructor_ = true;
		dock_chains_ = tag->getOption< std::string >( "interface" );
	} else if ( tag->hasOption( "ligandchain" ) ) {
		ligand_chain_ = tag->getOption< std::string >( "ligandchain" );
		explicit_constructor_ = true;
	} else {
		set_interface_jump(tag->getOption< core::Size >( "jump", 1 ) );
	}

	if ( tag->hasOption( "scorefile_reporting_prefix" ) ) {
		set_scorefile_reporting_prefix( tag->getOption<std::string>( "scorefile_reporting_prefix" ));
	}

	if ( tag->hasOption( "delta_metrics" ) ) {
		utility::vector1<core::simple_metrics::SimpleMetricCOP> delta_metrics = core::simple_metrics::get_metrics_from_datamap_and_subtags(tag, datamap, "delta_metrics");
		utility::vector1<core::simple_metrics::RealMetricCOP> delta_real_metrics;
		for ( auto const & metric : delta_metrics ) {
			if ( metric->simple_metric_type() != "RealMetric" ) { utility_exit_with_message("InterfaceAnalyzerMover: Only RealMetrics can be specified for the metrics tag!");}
			core::simple_metrics::RealMetricCOP delta_real_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::RealMetric const >( metric );
			delta_real_metrics.push_back( delta_real_metric );
		}
		set_metrics( delta_real_metrics );
	}

	//      tracer_(false), //output to tracer
	//      calcs_ready_(false), //calculators are not ready
	//      use_jobname_(false), //use the pose name
	//Having set_defaults here overrides several user-set values!
	//Default ctor exists to do this.  SML July 26 2016
	//set_defaults();
}

std::string
InterfaceAnalyzerMover::scorefile_column_name(std::string const & base_name) const
{
	return scorefile_reporting_prefix_ == "" ? base_name :
		scorefile_reporting_prefix_ + "_" + base_name;
}

void InterfaceAnalyzerMover::setup_score_data() {

	score_data_["dSASA_int"] = data_.dSASA[ total ];
	score_data_["dSASA_polar"] = data_.dSASA[ total ] - data_.dhSASA[ total ];
	score_data_["dSASA_hphobic"] = data_.dhSASA[ total ];
	score_data_["dG_separated"] = data_.dG[ total ];
	score_data_["dG_separated/dSASAx100"] = data_.dG_dSASA_ratio * 100.0;
	score_data_["delta_unsatHbonds"] = data_.delta_unsat_hbonds;
	score_data_["packstat"] = data_.packstat;
	score_data_["dG_cross"] = data_.crossterm_interface_energy;
	score_data_["dG_cross/dSASAx100"] = data_.crossterm_interface_energy_dSASA_ratio * 100.0;

	if ( use_centroid_ ) {
		score_data_["cen_dG"] = data_.centroid_dG;
	}
	score_data_["nres_int"] = data_.interface_nres[ total ];
	score_data_["per_residue_energy_int"] = per_residue_data_.regional_avg_per_residue_energy_int[ total ];
	score_data_["side1_score"] = data_.complexed_interface_score[ side1 ];
	score_data_["side2_score"] = data_.complexed_interface_score[ side2 ];
	score_data_["nres_all"] = included_nres_;
	score_data_["side1_normalized"] =  data_.complexed_interface_score[ side1 ] / ( core::Real( data_.interface_nres[ side1 ] ) );
	score_data_["side2_normalized"] = data_.complexed_interface_score[ side2 ] / ( core::Real( data_.interface_nres[ side2 ] ) );
	score_data_["complex_normalized"] = data_.complex_total_energy[ total ] / ( core::Real( included_nres_ ) );
	score_data_["hbond_E_fraction"] = data_.hbond_E_fraction;
	score_data_["sc_value"] = data_.sc_value;
	score_data_["hbonds_int"] = data_.interface_hbonds;

	if ( !data_.delta_metrics.empty() ) {
		for ( core::Size index = 1; index <= delta_metrics_.size(); ++index ) {
			score_data_["delta_" + delta_metrics_[index]->get_final_sm_type()] = data_.delta_metrics[index];
		}
	}

}

/// @details reports all the cool stuff we calculate to tracer output OR puts it into the job object.
void InterfaceAnalyzerMover::report_data(){
	basic::TracerImpl my_tr(posename_base_);
	//make output
	//std::ostringstream results_oss;
	//std::ostream & results = which_ostream(my_tr, results_oss, tracer_); //easy swap between tracer/job output

	//std::ostringstream interface_sele, missingHbond;

	//set up what data is reported in a string stream
	//report to job
	if ( tracer_ ) {
		//TR<<"Debugging print interface info:" << std::endl;
		my_tr << "TOTAL SASA: " << data_.complexed_SASA << std::endl;
		my_tr << "NUMBER OF RESIDUES: " << data_.interface_nres[ total ] << std::endl;
		my_tr << "AVG RESIDUE ENERGY: " << per_residue_data_.regional_avg_per_residue_energy_int[ total ] << std::endl;
		my_tr << "INTERFACE DELTA SASA: " << data_.dSASA[ total ] << std::endl;
		my_tr << "INTERFACE HYDROPHOBIC SASA: " << data_.dhSASA[ total ] << std::endl;
		my_tr << "INTERFACE POLAR SASA: " << data_.dSASA[ total ] - data_.dhSASA[ total ] << std::endl;
		my_tr << "CROSS-INTERFACE ENERGY SUMS: " << data_.crossterm_interface_energy << std::endl;
		my_tr << "SEPARATED INTERFACE ENERGY DIFFERENCE: " << data_.dG[ total ] << std::endl;
		my_tr << "CROSS-INTERFACE ENERGY/INTERFACE DELTA SASA: " << data_.crossterm_interface_energy_dSASA_ratio << std::endl;
		my_tr << "SEPARATED INTERFACE ENERGY/INTERFACE DELTA SASA: " << data_.dG_dSASA_ratio << std::endl;
		my_tr << "DELTA UNSTAT HBONDS: " << data_.delta_unsat_hbonds << std::endl;
		//my_tr << "ALL Gly INTERFACE ENERGY:  " << data_.gly_dG << std::endl; does not help
		if ( use_centroid_ ) {
			my_tr << "CENTROID dG: " << data_.centroid_dG << std::endl;
		}
		//my_tr << "AVG HBOND EXPOSURE RATIO: " << hbond_exposure_ratio_ << std::endl; does not help
		//my_tr << "HBOND SASA / INTERFACE dSASA: " << total_hb_sasa_ / interface_delta_sasa_ << std::endl;
		my_tr << "CROSS INTERFACE HBONDS: " << data_.interface_hbonds << std::endl;
		my_tr << "HBOND ENERGY: " << data_.total_hb_E << std::endl;
		my_tr << "HBOND ENERGY/ SEPARATED INTERFACE ENERGY: " << data_.total_hb_E / data_.dG[ total ] << std::endl;
		if ( compute_packstat_ ) {
			my_tr << "INTERFACE PACK STAT: " << data_.packstat << std::endl;
		}
		if ( compute_interface_sc_ ) {
			my_tr << "SHAPE COMPLEMENTARITY VALUE: " << data_.sc_value << std::endl;
		}
		if ( !data_.delta_metrics.empty() ) {
			for ( core::Size index = 1; index <= delta_metrics_.size(); ++index ) {
				my_tr << "Delta " + delta_metrics_[index]->get_final_sm_type() + ": " << data_.delta_metrics[index] << std::endl;
			}
		}

	} else {
		//or report to job
		for ( auto & it : score_data_ ) {
			protocols::jd2::add_string_real_pair_to_current_job(it.first, it.second);
		}
	}
}

void InterfaceAnalyzerMover::add_score_info_to_pose( core::pose::Pose & pose ){

	if ( score_data_.size() == 0 ) {
		TR << "No extra scores to add to pose..." << std::endl;
		return;
	}

	typedef std::map< std::string , core::Real >::const_iterator it_type;
	for ( it_type it = score_data_.begin(); it != score_data_.end(); it++ ) {
		core::pose::setPoseExtraScore(pose, scorefile_column_name(it->first), it->second);
	}
}

/// @details prints tracer output of pymol selction of interface residues, also builds a pymol selection that can be used from a file.
void InterfaceAnalyzerMover::print_pymol_selection_of_interface_residues( core::pose::Pose const & pose, std::set< core::Size > const & interface_set )
{
	//for tracer or job output
	std::ostringstream interface_oss;
	std::ostream & pymol_interface = which_ostream( TRinterface, interface_oss, tracer_ );

	//setup naming of pymol objects
	std::ostringstream interface_sele;
	interface_sele << std::endl;
	std::string pymol_obj_for_interface_sel;
	if ( compute_packstat_ ) {
		pymol_obj_for_interface_sel = posename_base_ + "_fullpose_pack";
	} else {
		pymol_obj_for_interface_sel = posename_base_;
	}
	//setup the tracer output
	pymol_interface << "pymol-style selection for interface res \n"
		<< "select " << posename_base_ << "_interface, ";
	//itterate through the interface set and build the selection syntaxt
	bool first_sel_complete( false );
	core::Size resnum;
	std::string chain_last( "" );
	for ( core::Size it : interface_set ) {
		//sets the current values
		resnum = pose.pdb_info()->number( it );
		std::string chain = pose.pdb_info()->chain( it );
		//special print if the first time through
		if ( !first_sel_complete ) {
			interface_sele << "cmd.select(\"/" << pymol_obj_for_interface_sel << "//" << chain << "/"
				<< resnum << "\")" << std::endl;
			pymol_interface << "/" << posename_base_ << "//" << chain << "/" << resnum << "+" ;
			first_sel_complete = true;
		} else if ( chain != chain_last ) {
			interface_sele << "cmd.select(\"sele + /" << pymol_obj_for_interface_sel << "//" << chain << "/"
				<< resnum << "\")" << std::endl;
			pymol_interface <<" + "<< "/" << posename_base_ << "//" << chain << "/" << resnum << "+";
		} else {
			interface_sele << "cmd.select(\"sele + /" << pymol_obj_for_interface_sel << "//" << chain << "/"
				<< resnum << "\")" << std::endl;
			pymol_interface << resnum << "+";
		}
		chain_last = chain;
	} //end itterate over interface
	//finish up
	pymol_interface << std::endl;
	interface_sele << "cmd.create(\"" << posename_base_ << "_interface_sel\", \"sele\")" << std::endl;
	data_.pymol_sel_interface =  interface_sele.str() ;
	//job output if wanted
	if ( !tracer_ ) {
		protocols::jd2::add_string_to_current_job( interface_oss.str() );
	}
	return;
}//end


/// @details This function reports a few things: a pymol sytle selection of the unstat atoms and reports to the tracer or job what these atoms are.  The app InterfaceAnalyzer gets the multi-line string to write a file or print the selection.  Unsat hbonds to be shown as Spheres
void InterfaceAnalyzerMover::print_pymol_selection_of_hbond_unsat( core::pose::Pose & pose, utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector )
{
	//for tracer or job output
	std::ostringstream results_oss, unsathbond_oss;
	basic::TracerImpl my_tr(posename_base_);
	std::ostream & results = which_ostream( my_tr, results_oss, tracer_); //easy swap between tracer/job output
	std::ostream & unsathbond = which_ostream( TRhbonds, unsathbond_oss, tracer_ );
	results << "Residues missing H-bonds:" << std::endl;
	results << "Residue \t Chain \t Atom " << std::endl;
	bool first_sel_complete( false );
	std::ostringstream missingHbond; // for pymol selection
	missingHbond << std::endl;

	//setup the tracer output
	unsathbond << "pymol-style selection for unstat hbond res \n"
		<< "select " << posename_base_ << "_unsat, ";
	//setup for looping over all unstat hbonds
	core::Size resnum;
	std::string chain_last ("");
	std::string atomname ;
	for ( core::Size i( 1 ); i <= delta_unsat_hbond_atid_vector.size(); i++ ) {
		core::id::AtomID const id ( delta_unsat_hbond_atid_vector[ i ] );
		resnum = pose.pdb_info()->number( id.rsd() );
		std::string chain = pose.pdb_info()->chain( id.rsd() );
		atomname = pose.residue(id.rsd()).atom_name(id.atomno());
		//get rid of whitespace in the atomname
		std::string temp;
		for ( char j : atomname ) {
			if ( j != ' ' ) { temp += j; }
		}
		atomname = temp;
		//do the tracer/job output
		results << resnum << " \t " << chain << " \t "<< atomname << std::endl;
		//now setup pymol output
		if ( !first_sel_complete ) {
			missingHbond << "cmd.select(\"/" << posename_base_ << "//" << chain << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond << "/" << posename_base_ << "//" << chain << "/" << resnum << "+" ;
			first_sel_complete = true;
		} else if ( chain != chain_last ) {
			missingHbond << "cmd.select(\"sele + /" << posename_base_ << "//" << chain << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond <<" + "<< "/" << posename_base_ << "//" << chain << "/" << resnum << "+";
		} else {
			missingHbond << "cmd.select(\"sele + /" << posename_base_ << "//" << chain << "/" << resnum << "/" << atomname << "\")"<< std::endl;
			unsathbond << resnum << "+";
		}
		chain_last = chain;
	} //end itterate over all unsat AtomIDs
	unsathbond << std::endl;
	//finalize output
	missingHbond << "cmd.create(\"" << posename_base_ << "_unsat_hbond\", \"sele\")" << std::endl;
	missingHbond << "cmd.show(\"spheres\", \"" << posename_base_ << "_unsat_hbond\")" << std::endl;
	data_.pymol_sel_hbond_unsat = missingHbond.str() ;

	//if we're not doing tracer output report this stuff to the job
	if ( !tracer_ ) {
		protocols::jd2::add_string_to_current_job( results_oss.str() );
		protocols::jd2::add_string_to_current_job( unsathbond_oss.str() );
	}
	return;
}


/// @details This function doesn't do the printing itself.  The app InterfaceAnalyzer gets the multi-line string to write a file or print the selection
/// @details From best packing to worse packing, colors go as Blue, Purple, Pink, Red
void InterfaceAnalyzerMover::print_pymol_selection_of_packing( core::pose::Pose const & pose, utility::vector1< core::Real >&  interface_pack_scores) {
	std::ostringstream pymol_packing;
	pymol_packing << std::endl;

	std::string pymol_object_fullpose_pack = posename_base_ + "_fullpose_pack";
	pymol_packing << "cmd.create(\"" << pymol_object_fullpose_pack << "\", \"" << posename_base_ << "\")" << std::endl;

	///////////////////////////////////////////// WHOLE POSE COLOR BY PACKING /////////////////////
	TR << "Total: "<< pose.size() <<" PackScoresTotal: "<<interface_pack_scores.size() << std::endl;

	//JAB  - This can sometimes not match - Temp fix.  I've emailed Will to help fix this bug in packstat.
	if ( pose.size() != interface_pack_scores.size() ) {
		TR << "Packscores and residues do not match. Skipping pymol selection of packing" << std::endl;
		return;
	}

	for ( core::Size i( 1 ); i <= interface_pack_scores.size(); ++i ) {
		if ( !include_residue_[ i ] ) { continue; }

		core::Size resnum = pose.pdb_info()->number( i );
		std::string chain = pose.pdb_info()->chain( i );
		core::Size color;
		if      ( interface_pack_scores[ i ] >= 0.75 ) { color = 2; }  //blue
		else if ( interface_pack_scores[ i ] >= 0.50 ) { color = 16; } //purple
		else if ( interface_pack_scores[ i ] >= 0.25 ) { color = 12; } //pink
		else if ( interface_pack_scores[ i ] >  0.00 ) { color = 4; }  //red
		else { color = 24; } //gray, something went wrong

		pymol_packing << "cmd.select(\"/" << pymol_object_fullpose_pack << "//" << chain << "/" << resnum << "\")" << std::endl;
		pymol_packing << "cmd.color(" << color << ", \"sele\")" << std::endl;
	}
	data_.pymol_sel_packing = pymol_packing.str() ;
	return;
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InterfaceAnalyzerMover::fresh_instance() const
{
	return utility::pointer::make_shared< InterfaceAnalyzerMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InterfaceAnalyzerMover::clone() const
{
	return utility::pointer::make_shared< InterfaceAnalyzerMover >( *this );
}

core::Real InterfaceAnalyzerMover::get_interface_dG() const { return data_.dG[ total ]; } //previous functionality: redundant with get_separated_interface_energy, but supports other protocols

/// @details getters
core::Real InterfaceAnalyzerMover::get_complexed_sasa() { return data_.complexed_SASA; }
core::Real InterfaceAnalyzerMover::get_interface_delta_sasa() { return data_.dSASA[ total ]; }
core::Real InterfaceAnalyzerMover::get_separated_interface_energy() { return data_.dG[ total ]; }
core::Size InterfaceAnalyzerMover::get_num_interface_residues() {return data_.interface_nres[ total ];}
core::Real InterfaceAnalyzerMover::get_complex_energy() { return data_.complex_total_energy[ total ]; }
core::Real InterfaceAnalyzerMover::get_per_residue_energy() { return per_residue_data_.regional_avg_per_residue_energy_int[ total ];}
core::Real InterfaceAnalyzerMover::get_crossterm_interface_energy() { return data_.crossterm_interface_energy; }
core::Real InterfaceAnalyzerMover::get_separated_interface_energy_ratio() { return data_.dG_dSASA_ratio; }
core::Real InterfaceAnalyzerMover::get_crossterm_interface_energy_ratio() { return data_.crossterm_interface_energy_dSASA_ratio; }
core::Real InterfaceAnalyzerMover::get_interface_packstat() { return data_.packstat; }
core::Size InterfaceAnalyzerMover::get_interface_delta_hbond_unsat() { return data_.delta_unsat_hbonds; }
std::string InterfaceAnalyzerMover::get_pymol_sel_interface() { return data_.pymol_sel_interface; }
std::string InterfaceAnalyzerMover::get_pymol_sel_hbond_unsat() { return data_.pymol_sel_hbond_unsat; }
std::string InterfaceAnalyzerMover::get_pymol_sel_packing() { return data_.pymol_sel_packing; }
bool InterfaceAnalyzerMover::get_multichain_constructor() { return explicit_constructor_; }
std::set<int> InterfaceAnalyzerMover::get_fixed_chains() { return fixed_chains_; }
std::set<core::Size> InterfaceAnalyzerMover::get_interface_set() { return interface_set_;}
InterfaceAnalyzerMover::group_set InterfaceAnalyzerMover::get_chain_groups() { return chain_groups_; }
bool InterfaceAnalyzerMover::get_pack_input() { return pack_input_; }
core::Real InterfaceAnalyzerMover::get_gly_interface_energy() { return data_.gly_dG ; }
core::Real InterfaceAnalyzerMover::get_centroid_dG() { return data_.centroid_dG ; }
//core::Real InterfaceAnalyzerMover::get_interface_Hbond_sasa() {return total_hb_sasa_; }
//core::Real InterfaceAnalyzerMover::get_Hbond_exposure_ratio() {return hbond_exposure_ratio_; }
core::Real InterfaceAnalyzerMover::get_total_Hbond_E() { return data_.total_hb_E; }

bool InterfaceAnalyzerMover::get_use_centroid_dG() const {return use_centroid_;}

utility::vector1< core::simple_metrics::RealMetricCOP > InterfaceAnalyzerMover::get_metrics() { return delta_metrics_; }

/// @details setters
void InterfaceAnalyzerMover::set_use_resfile( bool const use_resfile ) { use_resfile_ = use_resfile; }
void InterfaceAnalyzerMover::set_use_centroid_dG( bool const use_centroid ) { use_centroid_ = use_centroid; }
void InterfaceAnalyzerMover::set_compute_packstat( bool const compute_packstat ) { compute_packstat_ = compute_packstat; }
void InterfaceAnalyzerMover::set_interface_jump( core::Size const interface_jump ) {

	interface_jump_ = interface_jump;
	explicit_constructor_ = false;
}

void InterfaceAnalyzerMover::set_interface( std::string const & interface ){
	dock_chains_ = interface;
	explicit_constructor_ = true;
}

void InterfaceAnalyzerMover::set_use_tracer( bool const tracer) {tracer_ = tracer;}
//void InterfaceAnalyzerMover::set_calcs_ready(bool const calcs_ready) {calcs_ready_ = calcs_ready;}
void InterfaceAnalyzerMover::set_use_jobname( bool const use_jobname ) { use_jobname_ = use_jobname; }
void InterfaceAnalyzerMover::set_scorefile_reporting_prefix( std::string const & prefix ) { scorefile_reporting_prefix_ = prefix; }
void InterfaceAnalyzerMover::set_scorefunction( core::scoring::ScoreFunctionCOP sf ) { sf_ = sf->clone(); }

void InterfaceAnalyzerMover::set_metrics( utility::vector1<core::simple_metrics::RealMetricCOP> const & delta_metrics ) { delta_metrics_ = delta_metrics; }

// Packing the input and/or separated pose changes how results should be interpretted
// Thus, using tracers to inform user of intended IAM packing behavior
void InterfaceAnalyzerMover::set_pack_input( bool const pack_input ) {
	if ( pack_input ) { TR << "Repacking input Pose before calculating bound state values" << std::endl; }
	pack_input_ = pack_input;
}
void InterfaceAnalyzerMover::set_pack_separated( bool const pack_separated ) {
	if ( pack_separated ) { TR << "Repacking separated Pose before running calculations" << std::endl; }
	pack_separated_ = pack_separated;
}

std::string InterfaceAnalyzerMover::get_name() const {
	return mover_name();
}

std::string InterfaceAnalyzerMover::get_scorefile_reporting_prefix() const {
	return scorefile_reporting_prefix_;
}

std::string InterfaceAnalyzerMover::mover_name() {
	return "InterfaceAnalyzerMover";
}

void InterfaceAnalyzerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "pack_separated", xsct_rosetta_bool, "repack chains after separation (better evaluation of unbound state)", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "pack_input", xsct_rosetta_bool, "repack the input before calculating bound state values", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "resfile", xsct_rosetta_bool, "use a resfile for packing operations note BOOLEAN", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "packstat", xsct_rosetta_bool, "run packstat on bound state", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "interface_sc", xsct_rosetta_bool, "run interface_sc on bound state", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "tracer", xsct_rosetta_bool, "report output to a Tracer", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "use_jobname", xsct_rosetta_bool, "when reporting output (especially to tracer), label with the job name (like 1UBQ_0001)", "false" );

	std::string const exclusive_warning(" Notice that 'jump', 'interface', 'fixedchains', and 'ligandchain' are all mutually exclusive.");

	attlist + XMLSchemaAttribute( "fixedchains", xs_string, "comma separated list of PDB chains on one side of the interface, most useful when you have multiple chains on each side. "+ exclusive_warning );
	attlist + XMLSchemaAttribute( "interface", xs_string, "docking-style interface definitions, like HL_AB for chains HL vs AB for the interface. " + exclusive_warning );
	attlist + XMLSchemaAttribute( "ligandchain", xs_string, "Move ONLY this PDB chain. " + exclusive_warning );
	attlist + XMLSchemaAttribute( "jump", xsct_non_negative_integer, "Residues upstream/downstream of this Jump are on opposite sides of the interface; this Jump moves in the separation step. " + exclusive_warning );
	attlist + XMLSchemaAttribute( "scorefile_reporting_prefix", xs_string, "Prefix to add to column names for the values that are put in the Pose for eventual output into a score file. An extra underscore is added between the value specified here and the standard column names if this option is used.");

	attlist + XMLSchemaAttribute( "delta_metrics", xs_string, "A comma separated list of RealMetrics to calculate the difference between complexed and single state for.");


	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Authors: Steven Lewis, Bryan Der, Ben Stranges, Jared Adolf-Bryfogle\n"
		"This mover reports scores and statistics useful for judging the quality of protein-protein interfaces", attlist );
}

std::string InterfaceAnalyzerMoverCreator::keyname() const {
	return InterfaceAnalyzerMover::mover_name();
}

protocols::moves::MoverOP
InterfaceAnalyzerMoverCreator::create_mover() const {
	return utility::pointer::make_shared< InterfaceAnalyzerMover >();
}

void InterfaceAnalyzerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceAnalyzerMover::provide_xml_schema( xsd );
}


}//analysis
}//protocols
