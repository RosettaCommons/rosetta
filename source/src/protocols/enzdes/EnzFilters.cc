// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/EnzFilters.cc
/// @brief
/// @author Sagar Khare (khares@uw.edu)
/// @author Roland A Pache
/// @minor-corrections Per Greisen (pgreisen@gmail.com) PG

// unit headers
#include <protocols/enzdes/EnzFilters.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>

// package headers
#include <protocols/match/downstream/LigandConformerBuilder.hh> //for ResidueConformerFilter
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/enzdes/EnzdesMovers.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.hh>
#include <protocols/enzdes/enzdes_util.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/scoring/Interface.hh>
#include <core/id/AtomID.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh> //for PDB-info-to-resid functionality
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <basic/MetricValue.hh>

#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/ChargeCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/SurfaceCalculator.hh>

#include <numeric/random/random.hh>


//Objectxxxx header
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

//utility headers
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// 21-05-2013 PG
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "protocols.enzdes.EnzFilters" );

namespace protocols {
namespace enzdes {

using namespace protocols::filters;
using namespace protocols::moves;
using namespace utility::tag;

using core::pose::Pose;

LigDSasaFilter::LigDSasaFilter( core::Real const lower_threshold, core::Real const upper_threshold ) : Filter( "DSasa" ), lower_threshold_( lower_threshold ), upper_threshold_(upper_threshold) {}

bool
LigDSasaFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const dsasa( compute( pose ) );

	TR<<"dsasa is "<<dsasa<<". ";
	if ( dsasa >= lower_threshold_ && dsasa <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
LigDSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const dsasa( compute( pose ));
	out<<"DSasa= "<< dsasa<<'\n';
}

core::Real
LigDSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const dsasa( compute( pose ));
	return( dsasa );
}

core::Real
LigDSasaFilter::compute( core::pose::Pose const & pose ) const {
	basic::MetricValue< core::Real > mv_sasa;
	core::Real sasa (0.0);
	core::Size lig_chain(2), prot_chain;
	for ( core::Size i = 1;  i<= pose.size(); ++i ) {
		if ( pose.residue_type( i ).is_ligand() ) lig_chain = pose.chain( i );
	}
	prot_chain = pose.chain( 1 );   //we're making the not so wild assumption that the the first residue in the pose belongs to the protein
	//std::string lig_ch_string = utility::to_string( lig_chain );
	//std::string prot_ch_string = utility::to_string( prot_chain );
	if ( lig_chain == prot_chain ) { utility_exit_with_message( "WTF?!? ligand and residue 1 are on the same chain... " );}
	else {
		TR<<"Now calculating SaSa"<<std::endl;
		pose.metric( "sasa_interface", "frac_ch2_dsasa", mv_sasa); // if this seems like its come out of nowhere see declaration in protocols/jd2/DockDesignParser.cc
		sasa =  mv_sasa.value() ;
	}
	return( sasa );
}

void
LigDSasaFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	lower_threshold_ = tag->getOption<core::Real>( "lower_threshold", 0 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 1 );

	TR<<"LigDSasaFilter with lower threshold of "<<lower_threshold_<<" and upper threshold of "<< upper_threshold_ <<std::endl;
}

LigDSasaFilter::~LigDSasaFilter() = default;

DiffAtomSasaFilter::DiffAtomSasaFilter( core::Size resid1, core::Size resid2, std::string atomname1, std::string atomname2, std::string sample_type ) : Filter( "DiffAtomBurial" ), resid1_( resid1 ), resid2_( resid2), aname1_(std::move( atomname1 )), aname2_(std::move( atomname2 )), sample_type_(std::move( sample_type )) {}

bool
DiffAtomSasaFilter::apply( core::pose::Pose const & pose ) const {
	bool const dsasa( compute( pose ) );

	if ( dsasa ) {
		TR<<"Diff Atom Sasa filter passing." <<std::endl;
		return true;
	} else {
		TR<<"Diff Atom Sasa filter failing."<<std::endl;
		return false;
	}
}

void
DiffAtomSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	bool const dsasa( compute( pose ));
	out<<"Diff Sasa= "<< dsasa<<'\n';
}

core::Real
DiffAtomSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	bool const dsasa( compute( pose ));
	if ( dsasa ) return 1;
	else return 0;
}

bool
DiffAtomSasaFilter::compute( core::pose::Pose const & pose ) const {
	basic::MetricValue< id::AtomID_Map< core::Real > > atom_sasa;
	bool setting (false);
	pose.metric( "sasa_interface", "delta_atom_sasa", atom_sasa );
	core::id::AtomID const atomid1 (core::id::AtomID( pose.residue( resid1_ ).atom_index(aname1_), resid1_ ));
	core::id::AtomID const atomid2 (core::id::AtomID( pose.residue( resid2_ ).atom_index(aname2_), resid2_ ));
	core::Real const atom1_delta_sasa( atom_sasa.value()( atomid1 ) );
	core::Real const atom2_delta_sasa( atom_sasa.value()( atomid2 ) );
	TR<<"Atom1: Dsasa is "<< atom1_delta_sasa <<" and Atom2 Dsasa is "<< atom2_delta_sasa << std::endl;
	if ( sample_type_ == "more" ) {
		if ( atom1_delta_sasa > atom2_delta_sasa ) setting=true;
	} else if ( atom1_delta_sasa < atom2_delta_sasa ) setting=true; //sample_type == less i.e. see if a1 is less buried than a2 upon binding
	return setting;
}

void
DiffAtomSasaFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const &pose )
{
	if ( tag->hasOption("res1_res_num") ) resid1_ =  tag->getOption<core::Size>( "res1_res_num", 0 );
	if ( tag->hasOption("res2_res_num") ) resid2_ =  tag->getOption<core::Size>( "res2_res_num", 0 );
	if ( tag->hasOption("res1_pdb_num") ) resid1_ = core::pose::get_resnum(tag, pose, "res1_");
	if ( tag->hasOption("res2_pdb_num") ) resid2_ = core::pose::get_resnum(tag, pose, "res2_");
	if ( resid1_==0 || resid2_==0 ) { //ligand
		protocols::ligand_docking::LigandBaseProtocol ligdock;
		if ( resid1_==0 ) resid1_ = ligdock.get_ligand_id(pose);
		if ( resid2_==0 ) resid2_ = ligdock.get_ligand_id(pose);
	}
	aname1_  =  tag->getOption<std::string>("atomname1", "CA" );
	aname2_  =  tag->getOption<std::string>("atomname2", "CA" );
	sample_type_  =  tag->getOption<std::string>("sample_type", "more" ); //a1 is more buried than a2 i.e. dsasa is more
	runtime_assert(resid1_>0 && resid1_<=pose.size() );
	runtime_assert(resid2_>0 && resid2_<=pose.size() );
	runtime_assert (pose.residue( resid1_ ).has( aname1_ ));
	runtime_assert (pose.residue( resid2_ ).has( aname2_ ));
	runtime_assert( sample_type_=="more" || sample_type_ == "less" );
	TR<<" Defined LigDSasaFilter "<< std::endl;
}

DiffAtomSasaFilter::~DiffAtomSasaFilter() = default;

bool
LigBurialFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	TR<<"Number of interface neighbors of ligand is "<<count_neighbors<<std::endl;
	return( count_neighbors >= neighbors_ );
}

void
LigBurialFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	out<<"Number of interface neighbors of residue is "<<count_neighbors<<'\n';
}

core::Real
LigBurialFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	return( count_neighbors );
}

/// @details counts the number of residues contacting the ligand
core::Size
LigBurialFilter::compute( core::pose::Pose const & pose ) const {

	core::Size real_lig_id (lig_id_);
	if ( real_lig_id==0 ) {
		protocols::ligand_docking::LigandBaseProtocol ligdock;
		real_lig_id = ligdock.get_ligand_id(pose);
		TR<<"Calculating neighbors of Ligand resid is " << real_lig_id << std::endl;
	}
	core::Size count_neighbors( 0 );
	core::conformation::Residue const res_target( pose.conformation().residue( real_lig_id ) );
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		core::conformation::Residue const resi( pose.residue( i ) );
		core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
		if ( distance <= distance_threshold_ ) ++count_neighbors;
	}
	return( count_neighbors);
}

/// @details: this filter basically works exactly as ResidueBurialFilter, but with the advantage that it has the capability to
/// @details: figure out resid of the ligand
void
LigBurialFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const &)
{
	lig_id_ =  tag->getOption<core::Size>( "lig_id", 0 );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );
	neighbors_ = tag->getOption<core::Size>( "neighbors", 1 );

	TR<<"LigBurialFilter with distance threshold of "<<distance_threshold_<<" around residue "<<lig_id_<<" with "<<neighbors_<<" neighbors."<<std::endl;
	TR.flush();
}

LigBurialFilter::~LigBurialFilter() = default;

LigInterfaceEnergyFilter::LigInterfaceEnergyFilter( core::scoring::ScoreFunctionOP scorefxn, core::Real const threshold, bool const include_cstE, core::Size const rb_jump, core::Real const interface_distance_cutoff ) : Filter( "LigInterfaceEnergy" ), threshold_( threshold ),  include_cstE_ ( include_cstE ), rb_jump_ ( rb_jump ), interface_distance_cutoff_ ( interface_distance_cutoff )  {

	using namespace core::scoring;

	if ( scorefxn ) scorefxn_ = scorefxn->clone();
	if ( !include_cstE ) enzutil::enable_constraint_scoreterms( scorefxn_);
}

LigInterfaceEnergyFilter::LigInterfaceEnergyFilter( LigInterfaceEnergyFilter const &init ) :
	//utility::pointer::ReferenceCount(),
	Filter( init ), threshold_( init.threshold_ ),
	include_cstE_ (init.include_cstE_), rb_jump_ (init.rb_jump_), interface_distance_cutoff_ ( init.interface_distance_cutoff_){
	using namespace core::scoring;
	if ( init.scorefxn_ ) scorefxn_ = init.scorefxn_->clone();
}


bool
LigInterfaceEnergyFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;

	if ( pose.conformation().num_chains() < 2 ) {
		TR << "pose must contain at least two chains!" << std::endl;
		return false;
	} else {
		TR<<" \n \t --------- computing  ---------- \n \t-------- ligand interface energies   --------\n \t \t------------------- \n" << std::endl;
		core::Real interf_E (compute(pose));
		return (interf_E < threshold_);
	}
}


void
LigInterfaceEnergyFilter::report( std::ostream & out, core::pose::Pose const & pose ) const

{
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

	core::pose::Pose in_pose = pose;
	// TL 4/2013: If rb_jump_ is unset, determine jump number according to the number of jumps currently in the pose
	// previously, if this filter was used in RosettaScripts, rb_jump_ was initialized to jump number according to the number of jumps in the pose when the XML was parsed. This enables adding a ligand to the pose after the XML is parsed
	core::Size jump;
	if ( rb_jump_ == 0 ) {
		jump = pose.num_jump();
	} else {
		jump = rb_jump_;
	}
	FArray1D_bool prot_res( in_pose.size(), false );
	in_pose.fold_tree().partition_by_jump( jump, prot_res );
	protocols::scoring::Interface interface_obj( jump );
	in_pose.update_residue_neighbors();
	interface_obj.distance( interface_distance_cutoff_ );
	interface_obj.calculate( in_pose );
	(*scorefxn_)( in_pose );

	out<<"\n"<<A(9, "chain")<<A( 9, "res")<<A( 9, "AA")<<A( 9, "total")<<A( 9, "contact")<<A( 9, "fa_atr")<<A( 9, "fa_rep")<<A( 9, "hb_bb_sc")<<A( 9, "hb_sc")<<A( 9, "fa_sol")<<A( 9, "fa_dun")<<A( 9, "fa_pair")<<"\n";
	for ( core::Size resnum_ = 1; resnum_ <= pose.size(); ++resnum_ ) {
		//   if ( !in_pose.residue(resnum_).is_protein() ) continue;
		if ( interface_obj.is_interface( resnum_ ) ) { // in interface

			Real total=in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( total_score ) ];
			if ( !include_cstE_ ) total-= constraint_energy(in_pose, resnum_ );
			Real weighted_fa_atr=( (*scorefxn_)[ ScoreType( fa_atr) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( fa_atr ) ]);
			Real weighted_fa_rep=( (*scorefxn_)[ ScoreType( fa_rep) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( fa_rep) ]);
			Real weighted_hbond_bb_sc=( (*scorefxn_)[ ScoreType( hbond_bb_sc) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( hbond_bb_sc ) ]);
			Real weighted_hbond_sc=( (*scorefxn_)[ ScoreType( hbond_sc) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( hbond_sc ) ]);
			Real weighted_fa_sol=( (*scorefxn_)[ ScoreType( fa_sol) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( fa_sol ) ]);
			Real weighted_contact_score = weighted_fa_atr + weighted_fa_rep + weighted_hbond_bb_sc + weighted_hbond_sc + weighted_fa_sol;

			Real weighted_fa_dun=( (*scorefxn_)[ ScoreType( fa_dun) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( fa_dun ) ]);
			Real weighted_fa_pair=( (*scorefxn_)[ ScoreType( fa_pair ) ] ) * ( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( fa_pair ) ]);

			out<<A (9 , in_pose.pdb_info()->chain( resnum_) )<< I(9,0, in_pose.pdb_info()->number(resnum_))<<A (9,in_pose.residue( resnum_).name3())
				<<F (9 , 3, total) <<" "
				<<F (9 , 3, weighted_contact_score)<<" "
				<<F (9 , 3, weighted_fa_atr) <<" "
				<<F (9 , 3, weighted_fa_rep) <<" "
				<<F (9 , 3, weighted_hbond_bb_sc )<<" "
				<<F (9 , 3, weighted_hbond_sc )<<" "
				<<F (9 , 3, weighted_fa_sol )<<" "
				<<F (9 , 3, weighted_fa_dun )<<" "
				<<F (9 , 3, weighted_fa_pair )<<"\n";


		}
	}

}


core::Real
LigInterfaceEnergyFilter::report_sm( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;

	core::Real const energy( compute( pose ) );
	return( energy );
}

core::Real
LigInterfaceEnergyFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	using namespace core::pose::metrics;
	core::pose::Pose in_pose = pose;
	in_pose.update_residue_neighbors();
	(*scorefxn_)( in_pose );
	protocols::ligand_docking::LigandBaseProtocol ligdock;
	core::Size const lig_id = ligdock.get_ligand_id(pose);
	enzutil::enable_constraint_scoreterms( scorefxn_);
	basic::MetricValue< core::Real > mv_interfE;
	in_pose.metric( "liginterfE","weighted_total", mv_interfE);  // if this seems like its come out of nowhere see declaration in protocols/jd2/DockDesignParser.cc
	core::Real weighted_score = mv_interfE.value();
	TR<<"Calculated Interface Energy is before cst correction is:"<<weighted_score<<std::endl;
	if ( !include_cstE_ ) weighted_score -= 2*constraint_energy(in_pose, lig_id);
	TR<<"Calculated Interface Energy is "<<weighted_score<<std::endl;
	return( weighted_score );
}

core::Real
LigInterfaceEnergyFilter::constraint_energy( core::pose::Pose const & in_pose , int which_res ) const
{
	using namespace core::scoring;
	core::pose::Pose pose (in_pose);
	(*scorefxn_)(pose);
	EnergyMap all_weights = pose.energies().weights();
	EnergyMap scores;

	if ( which_res == -1 ) { //means we want to know stuff for the whole pose
		scores = pose.energies().total_energies();
	} else { scores = pose.energies().residue_total_energies( which_res ); }

	return scores[ coordinate_constraint ] * all_weights[ coordinate_constraint ] + scores[atom_pair_constraint] * all_weights[ atom_pair_constraint] +
		scores[ angle_constraint ] * all_weights[ angle_constraint ] + scores[ dihedral_constraint ] * all_weights[ dihedral_constraint ];

}
void
LigInterfaceEnergyFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	threshold_ = tag->getOption<core::Real>( "energy_cutoff", 0.0 );
	include_cstE_ = tag->getOption<bool>( "include_cstE" , 0 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", 0 );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff" , 8.0 );

}

LigInterfaceEnergyFilter::~LigInterfaceEnergyFilter() = default;

EnzScoreFilter::EnzScoreFilter( core::Size const resnum, std::string  cstid, core::scoring::ScoreFunctionOP scorefxn, core::scoring::ScoreType const score_type, core::Real const threshold, bool const whole_pose, bool const is_cstE ) : Filter( "EnzScore" ), resnum_( resnum ), cstid_(std::move( cstid )), score_type_( score_type ), threshold_( threshold ),  whole_pose_ ( whole_pose ), is_cstE_ ( is_cstE ) {

	using namespace core::scoring;

	if ( scorefxn ) scorefxn_ = scorefxn->clone();
	//if(score_type_ == total_score || is_cstE_ ) enzutil::enable_constraint_scoreterms(scorefxn_);
	if ( score_type_ != total_score ) {
		core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
		scorefxn_->reset();
		scorefxn_->set_weight( score_type_, old_weight );
	}
}

EnzScoreFilter::EnzScoreFilter( EnzScoreFilter const &init ) :
	//utility::pointer::ReferenceCount(),
	Filter( init ), resnum_( init.resnum_ ), cstid_( init.cstid_ ),score_type_( init.score_type_ ), threshold_( init.threshold_ ),
	whole_pose_ (init.whole_pose_), is_cstE_ ( init.is_cstE_ ) {
	using namespace core::scoring;
	if ( init.scorefxn_ ) scorefxn_ = init.scorefxn_->clone();
}

bool
EnzScoreFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	core::Real energy (compute(pose));
	TR<<"Score found as "<<energy <<" and threshold is " << threshold_<<std::endl;
	return (energy<threshold_);
}

void
EnzScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const

{
	using namespace core::scoring;
	out<<"Weighted score "<<compute( pose )<<'\n';
}

core::Real
EnzScoreFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
EnzScoreFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose;
	using namespace core::scoring;
	PoseOP in_pose( new Pose( pose ) );
	core::Size resnum = resnum_;
	if ( resnum==0 ) { //means we want ligand scores
		protocols::ligand_docking::LigandBaseProtocol ligdock;
		resnum = ligdock.get_ligand_id(pose);
	}
	core::Real weight(0.0), score(0.0);
	if ( whole_pose_ ) {
		(*scorefxn_)( *in_pose );
		if ( !is_cstE_ ) {
			weight = ( (*scorefxn_)[ ScoreType( score_type_ ) ] );
			score = ( in_pose->energies().total_energies()[ ScoreType( score_type_ ) ]);
			TR<<"Raw Score is: " << score << std::endl;
		}
		if ( score_type_  == total_score ) {
			if ( pose.constraint_set()->has_constraints() ) { //check if constraints exist, and remove their score
				enzutil::enable_constraint_scoreterms(scorefxn_);
				(*scorefxn_)( *in_pose );
				LigInterfaceEnergyFilter liginter( scorefxn_,0.0,false, 1, 0.0 );
				core::Real all_cstE = liginter.constraint_energy(*in_pose, -1);// for whole pose;
				score -= all_cstE;
			}
			return( score );
		} else if ( is_cstE_ ) {
			TR<<"Evaluating constraint score for whole pose..."<<std::endl;
			enzutil::enable_constraint_scoreterms(scorefxn_);
			LigInterfaceEnergyFilter liginter( scorefxn_,0.0,false, 1, 0.0 );
			core::Real all_cstE = liginter.constraint_energy(*in_pose, -1); // for whole pose;
			TR<<"And constraint energy is..."<< all_cstE << std::endl;
			return( all_cstE ) ;
		} else {
			core::Real const weighted_score( weight * score );
			TR<<"Weighted Score is: " << weighted_score << std::endl;
			return( weighted_score );
		}
	} else { //score for a particular residue/cstid
		if ( ! (cstid_ =="") ) resnum = enzutil::get_resnum_from_cstid( cstid_, *in_pose );
		(*scorefxn_)( *in_pose );
		TR<<"For Resid:"<< resnum << std::endl;
		runtime_assert(resnum>0);
		if ( is_cstE_ || score_type_ == total_score ) {
			enzutil::enable_constraint_scoreterms(scorefxn_);
			(*scorefxn_)( *in_pose );
			LigInterfaceEnergyFilter liginter( scorefxn_,0.0,false, 1, 0.0 );
			core::Real resnum_cstE = liginter.constraint_energy(*in_pose, resnum);
			if ( is_cstE_ ) return resnum_cstE;
			else return (in_pose->energies().residue_total_energies( resnum )[ ScoreType( score_type_ )] - resnum_cstE );
		} else return( ( (*scorefxn_)[ ScoreType( score_type_ ) ] ) * (in_pose->energies().residue_total_energies( resnum )[ ScoreType( score_type_ )]) );
	}
	return 0.0;
}

void
EnzScoreFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &, Movers_map const &, core::pose::Pose const &pose)
{
	using namespace core::scoring;
	is_cstE_ = false;
	if ( tag->hasOption( "pdb_num" ) ) resnum_ = core::pose::get_resnum(tag, pose);
	else if ( tag->hasOption( "res_num" ) ) resnum_ =  tag->getOption<core::Size>( "res_num", 0 );
	cstid_ = tag->getOption<std::string>( "cstid", "" );

	scorefxn_ = protocols::rosetta_scripts::parse_score_function(tag, data)->clone();
	std::string sco_name = tag->getOption<std::string>( "score_type", "total_score" );
	if ( sco_name == "cstE" ) {
		is_cstE_ = true;
		score_type_=atom_pair_constraint; //dummy assignment, never actually used.
	} else score_type_ = core::scoring::score_type_from_name( sco_name );
	threshold_ = tag->getOption<core::Real>( "energy_cutoff", 0.0 );
	whole_pose_ = tag->getOption<bool>( "whole_pose" , 0 );

	//check to make sure one and only one of resnum, cstid and whole_pose are specified
	runtime_assert(tag->hasOption( "res_num" )|| tag->hasOption( "pdb_num" ) || tag->hasOption( "cstid" ) || whole_pose_==1 );
	runtime_assert(!( (tag->hasOption( "res_num" )|| tag->hasOption( "pdb_num" )) &&  tag->hasOption( "cstid" )));
	runtime_assert(!( (tag->hasOption( "res_num" )|| tag->hasOption( "pdb_num" )) &&  whole_pose_==1));
	runtime_assert(!(tag->hasOption( "cstid" ) &&  whole_pose_==1));

	if ( whole_pose_==1 ) {
		TR<<"energies for whole pose will be calculated "
			<< "\n and scorefxn " << protocols::rosetta_scripts::get_score_function_name(tag) <<" will be used" <<std::endl;
	} else {
		TR<<"EnzScoreFilter for residue or cstid with cutoff "<<threshold_<<std::endl;
	}
}
EnzScoreFilter::~EnzScoreFilter() = default;

RepackWithoutLigandFilter::RepackWithoutLigandFilter( core::scoring::ScoreFunctionOP scorefxn, core::Real rms_thresh, core::Real energy_thresh, utility::vector1< core::Size > rms_target_res  ) : Filter( "RepackWithoutLigand" ), scorefxn_(std::move(scorefxn)), rms_threshold_( rms_thresh ), energy_threshold_( energy_thresh ), target_res_( rms_target_res ) {}

bool
RepackWithoutLigandFilter::apply( core::pose::Pose const & pose ) const
{

	core::Real const value( compute( pose ) );

	if ( calc_dE_ ) {
		if ( value > energy_threshold_ ) return false;
		else return true;
	} else if ( calc_rms_ ) {
		if ( value > rms_threshold_ ) return false;
		else return true;
	}
	return false;
}

void
RepackWithoutLigandFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	//Output to tracer triggers the compute function separately. Since this is an expensive filter, we don't want to unncessarily compute, hence disabled.
	// Dummy use of pose to avoid compiler warning.
	//  core::Real value( compute( pose ));
	//  if (calc_rms_) out<<"RMS= "<< value<<'\n';
	//  else out<<"dEnergy = "<< value<<'\n';
	out<<"Ligand to take out is: "<< pose.size();
	out<<'\n';
}

core::Real
RepackWithoutLigandFilter::report_sm( core::pose::Pose const & pose ) const
{
	return ( compute( pose ));
}

core::Real
RepackWithoutLigandFilter::compute( core::pose::Pose const & pose ) const
{
	core::Real value(10000.);
	core::pose::Pose rnl_pose = pose;
	enzutil::disable_constraint_scoreterms(scorefxn_);
	(*scorefxn_)(rnl_pose) ;
	core::Real wl_score = rnl_pose.energies().total_energies()[total_score];

	enzutil::remove_all_enzdes_constraints( rnl_pose );

	protocols::enzdes::RepackLigandSiteWithoutLigandMover rnl( scorefxn_, false );
	// same length
	rnl.set_separate_prt_ligand( false );

	rnl.apply( rnl_pose );

	(*scorefxn_)(rnl_pose);
	core::Real nl_score = rnl_pose.energies().total_energies()[total_score];

	if ( calc_dE_ ) {
		TR<<"Total energy with ligand is: "<<wl_score<<" and total energy without ligand is "<<nl_score<<std::endl;
		return (wl_score - nl_score);
	} else if ( calc_rms_ ) {
		ObjexxFCL::FArray1D_bool rms_seqpos( pose.size(), false );
		utility::vector1< core::Size > trg_res;
		TR << " Length of initial pose " << rnl_pose.size() << std::endl;
		if ( rms_all_rpked_ ) {
			TR<<"Getting identities of all pack residues... "<< std::endl;
			core::pack::task::PackerTaskCOP rnl_ptask = rnl.get_ptask();
			for ( core::Size i = 1; i <= rnl_ptask->total_residue(); ++i ) {

				if ( rnl_ptask->residue_task( i ).being_packed() && pose.residue( i ).is_protein() ) {
					trg_res.push_back( i );
				}
			}
		} else if ( use_cstids_ ) {
			enzutil::get_resnum_from_cstid_list(cstid_list_, pose, trg_res);
		}

		trg_res.insert( trg_res.begin(), target_res_.begin(), target_res_.end() );
		auto last = std::unique( trg_res.begin(), trg_res.end() );
		trg_res.erase( last, trg_res.end() );

		TR<<"Calculating RMS for residues ";
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( ! (pose.residue_type( i ).is_ligand()) ) {
				utility::vector1< core::Size >::const_iterator resfind = find( trg_res.begin(), trg_res.end(), i );
				if ( resfind != trg_res.end() ) {
					rms_seqpos( i ) =true;
					TR<< i<< ", ";
				}
			}
		}

		TR<< std::endl;

		core::Real rmsd( core::scoring::rmsd_no_super_subset( pose, rnl_pose, rms_seqpos, core::scoring::is_protein_sidechain_heavyatom ) );
		// PG 21-05-2013
#ifdef WIN32
		if ( _isnan(rmsd) ) {
#else
		if ( std::isnan(rmsd) ) {
			runtime_assert(!std::isnan(rmsd));
#endif
			utility_exit_with_message( "RMSD is NaN - there is something is wrong with how the interface is defined");

		}

		TR<<"Total rms of requested region is: "<< rmsd <<std::endl;
		return rmsd;
	}
	return value; //should not get here
}

void
RepackWithoutLigandFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, Filters_map const &, Movers_map const &, core::pose::Pose const &pose )
{
	TR<<" Defining RepackWithoutLigandFilter "<< std::endl;
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	calc_rms_ = false; calc_dE_ = false; rms_all_rpked_ = false; use_cstids_ = false;
	runtime_assert( tag->hasOption("energy_threshold") ||  tag->hasOption("rms_threshold") );
	if ( tag->hasOption("rms_threshold") ) {
		rms_threshold_ =  tag->getOption<core::Real>("rms_threshold", 0.5 );
		calc_rms_ = true;
		if ( tag->hasOption("target_res") ) {
			std::string target_res = tag->getOption<std::string>("target_res", "" );
			if ( target_res=="all_repacked" ) rms_all_rpked_ = true;
			else target_res_ = core::pose::get_resnum_list(tag, "target_res",pose);
		}
		if ( tag->hasOption("target_cstids") ) {
			cstid_list_ = tag->getOption<std::string>("target_cstids", "" );
			use_cstids_ = true;
		}
	} else if ( tag->hasOption("energy_threshold") ) {
		energy_threshold_  =  tag->getOption<core::Real>("energy_threshold", 0.0 );
		calc_dE_ = true;
	}
	// 21-05-2013 PG
	if ( !calc_dE_ && target_res_.empty() && cstid_list_.empty() && !rms_all_rpked_  ) {
		throw utility::excn::EXCN_RosettaScriptsOption(" NO RESIDUES HAVE BEEN SPECIFIED FOR THE RMSD CALCULATION");
	}
	TR<<" Defined RepackWithoutLigandFilter "<< std::endl;
}

RepackWithoutLigandFilter::~RepackWithoutLigandFilter() = default;

ValueEvaluator::ValueEvaluator( CompareMode const mode, core::Real const cutoff )
: mode_(mode),
	cutoff_(cutoff)
{}

ValueEvaluator::~ValueEvaluator()= default;

bool
ValueEvaluator::value_passes( core::Real const value ) const
{
	if ( mode_ == SMALLER ) return value < cutoff_;
	else if ( mode_ == LARGER ) return value > cutoff_;
	else return value == cutoff_;
}

//initializing static member variable
std::map< std::string, std::map<std::string, ValueEvaluator > > EnzdesScorefileFilter::evaluator_map_;

EnzdesScorefileFilter::EnzdesScorefileFilter()
: Filter(),
	no_packstat_calc_( basic::options::option[basic::options::OptionKeys::enzdes::no_packstat_calculation] ),
	native_comparison_(basic::options::option[basic::options::OptionKeys::enzdes::compare_native].user() ),
	repack_no_lig_(basic::options::option[basic::options::OptionKeys::enzdes::final_repack_without_ligand]),
	keep_rnl_pose_(basic::options::option[basic::options::OptionKeys::enzdes::dump_final_repack_without_ligand_pdb]),
	rnl_pose_(/* NULL */),
	sfxn_(core::scoring::ScoreFunctionFactory::create_score_function("enzdes") ),
	enzcst_io_(/* NULL */),
	native_comp_(/* NULL */),
	reqfile_name_("")
{
	if ( native_comparison_ ) native_comp_ = DesignVsNativeComparisonOP( new protocols::enzdes::DesignVsNativeComparison() );
	residue_calculators_.clear();
	native_compare_calculators_.clear();
	silent_Es_.clear();
	relevant_scoreterms_.clear();
	spec_segments_.clear();

	relevant_scoreterms_.push_back( "total_score" );
	relevant_scoreterms_.push_back( "fa_rep" );
	relevant_scoreterms_.push_back( "hbond_sc" );
	relevant_scoreterms_.push_back( "all_cst" );
}

EnzdesScorefileFilter::EnzdesScorefileFilter( EnzdesScorefileFilter const & other )
: /*utility::pointer::ReferenceCount(),*/ Filter( other ),
	no_packstat_calc_( other.no_packstat_calc_ ),
	native_comparison_(other.native_comparison_ ),
	repack_no_lig_( other.repack_no_lig_),
	keep_rnl_pose_( other.keep_rnl_pose_),
	rnl_pose_( other.rnl_pose_),
	sfxn_( other.sfxn_->clone() ),
	enzcst_io_( other.enzcst_io_),
	residue_calculators_(other.residue_calculators_),
	native_compare_calculators_(other.native_compare_calculators_),
	native_comp_(other.native_comp_),
	silent_Es_(other.silent_Es_),
	relevant_scoreterms_(other.relevant_scoreterms_),
	spec_segments_(other.spec_segments_),
	reqfile_name_(other.reqfile_name_)
{}


EnzdesScorefileFilter::~EnzdesScorefileFilter()= default;


bool
EnzdesScorefileFilter::apply( core::pose::Pose const & pose ) const{
	runtime_assert( evaluator_map_.find(reqfile_name_) != evaluator_map_.end() );

	std::map< std::string, ValueEvaluator > const & evaluators( evaluator_map_.find(reqfile_name_)->second );

	this->examine_pose( pose );

	std::set<std::string> found_evaluators; //make sure all the desired values are returned by the scorefile

	for ( utility::vector1< core::io::silent::SilentEnergy >::const_iterator sco_it(silent_Es_.begin()), sco_end(silent_Es_.end()); sco_it != sco_end; ++sco_it ) {

		//std::cerr << "iterating " << sco_it->name() << " val " << sco_it->value() << std::endl;
		auto val_it(evaluators.find(sco_it->name() ) );
		if ( val_it != evaluators.end() ) {
			found_evaluators.insert( sco_it->name() );
			if ( !val_it->second.value_passes( sco_it->value() ) ) {
				TR << "EnzdesScorefileFilter returning false for parameter " << sco_it->name() << " with cutoff " << val_it->second.cutoff_ << " and value " << sco_it->value() << "." << std::endl;
				return false;
			}
		}
	} //loop over all score terms produced by scorefile

	if ( found_evaluators.size() != evaluators.size() ) {
		std::cerr << "Not all parameters from requirement file " << reqfile_name_ << " were found in EnzdesScorefileFilter values, filering likely not working correctly." << std::endl;
		std::cerr << "Missing parameters: ";
		 for ( auto const & evaluator : evaluators ) {
			if ( found_evaluators.find( evaluator.first ) == found_evaluators.end() ) std::cerr << evaluator.first << ", ";
		}
		std::cerr << std::endl;
	}

	return true;
}

/// @details not implemented yet
void
EnzdesScorefileFilter::parse_my_tag( utility::tag::TagCOP tag , basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	if ( tag->hasOption("requirements") ) reqfile_name_ =  tag->getOption<std::string>( "requirements","" );
	else throw utility::excn::EXCN_RosettaScriptsOption("For EnzdesScorefileFilter, a requirements file needs to be specified in the tag.");

	this->initialize_value_evaluators_from_file( reqfile_name_ );
}

void
EnzdesScorefileFilter::set_cstio(
	toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io )
{
	enzcst_io_ = enzcst_io;
}

void
EnzdesScorefileFilter::examine_pose(
	core::pose::Pose const & pose ) const
{
	using namespace core::io::silent;

	silent_Es_.clear();
	//flo jan 2010 with the new cst cache, some calcs that take the
	//ligand out of the pose will crash. for now, let's make a copy
	//of the pose and remove the cst cache from it
	core::pose::Pose calc_pose = pose;
	toolbox::match_enzdes_util::get_enzdes_observer( calc_pose )->set_cst_cache( nullptr );
	(*sfxn_)( calc_pose );

	if ( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER) ) {
		spec_segments_ = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >(pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER ) )->segments();
	}

	bool separate_out_constraints = false;
	auto cstfind = find( relevant_scoreterms_.begin(), relevant_scoreterms_.end(),"all_cst");
	if ( cstfind != relevant_scoreterms_.end() ) separate_out_constraints = true;

	setup_pose_metric_calculators( pose, separate_out_constraints );

	//first write out the relevant score terms for the pose total
	 for ( auto const & relevant_scoreterm : relevant_scoreterms_ ) {
		std::string sco_name = relevant_scoreterm;
		int width = std::max( 10, (int) sco_name.length() + 3 );

		SilentEnergy new_se;
		if ( relevant_scoreterm == "all_cst" ) {
			new_se = SilentEnergy ( sco_name, enzutil::sum_constraint_scoreterms(pose, -1 ), 1, width);
		} else if ( separate_out_constraints && ( relevant_scoreterm == "total_score" ) ) {
			core::Real desired_value = pose.energies().total_energies()[ core::scoring::score_type_from_name( relevant_scoreterm ) ] - enzutil::sum_constraint_scoreterms(pose, -1 );
			new_se = SilentEnergy(sco_name, desired_value, 1, width);
		} else {
			new_se = SilentEnergy ( sco_name, pose.energies().total_energies()[ core::scoring::score_type_from_name( relevant_scoreterm ) ] * pose.energies().weights()[ core::scoring::score_type_from_name( relevant_scoreterm ) ], 1 ,width);
		}
		silent_Es_.push_back( new_se );

	}
	//pose metric calculators for pose total
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > >::const_iterator totcalc_it = residue_calculators_.find( 0 );
	if ( totcalc_it != residue_calculators_.end() ) {

		utility::vector1< std::pair< std::string, std::string > > const & tot_calculators = totcalc_it->second;
		 for ( auto const & tot_calculator : tot_calculators ) {

			std::string calc_name = "tot_" + tot_calculator.first;
			if ( tot_calculator.first == "charges_pm" ) calc_name = "tot_" + tot_calculator.second;
			int width = std::max( 10, (int) calc_name.length() + 3 );

			core::Real calc_value;

			//following lines fairly hacky, but don't know a better solution at the moment
			if ( tot_calculator.first == "hbond_pm" || tot_calculator.first == "burunsat_pm" || tot_calculator.first == "NLconts_pm" || tot_calculator.second == "total_pos_charges" || tot_calculator.second == "total_neg_charges" ) {
				basic::MetricValue< core::Size > mval_size;
				calc_pose.metric( tot_calculator.first, tot_calculator.second, mval_size );
				calc_value = mval_size.value();
			} else if ( tot_calculator.first == "seq_recovery" ) {
				if ( toolbox::match_enzdes_util::get_enzdes_observer( calc_pose ) -> get_seq_recovery_cache() ) {
					calc_value = toolbox::match_enzdes_util::get_enzdes_observer( calc_pose ) -> get_seq_recovery_cache() -> sequence_recovery( calc_pose );
				} else {
					calc_value= 0.0;
				}
			} else {
				basic::MetricValue< core::Real > mval_real;
				calc_pose.metric( tot_calculator.first, tot_calculator.second, mval_real );
				calc_value = mval_real.value();
			}

			SilentEnergy new_se( calc_name, calc_value, 1, width);
			silent_Es_.push_back( new_se );
		}

	}

	//then write out the relevant scoreterms (and potentially pose metrics) for each of the special residues
	Size spec_res_counter(0);
	utility::vector1< Size > special_res = enzutil::catalytic_res( pose );

	for ( utility::vector1<Size>::const_iterator res_it = special_res.begin(), end = special_res.end(); res_it != end; ++res_it ) {

		spec_res_counter++;
		//for convenience, the sequence number of the residue will be written out
		std::stringstream temp;
		temp << spec_res_counter;
		std::string spec_res_name = "SR_" + temp.str();
		//std::cerr << "name for res " << *res_it << " is " << spec_res_name ;
		SilentEnergy res_name(spec_res_name, *res_it, 1, 10);
		silent_Es_.push_back( res_name );

		utility::vector1< core::Size > dummy_vec;
		dummy_vec.push_back( *res_it );
		compute_metrics_for_residue_subset( spec_res_name, separate_out_constraints, calc_pose, dummy_vec);
	}

	for ( core::Size specseg(1); specseg <= spec_segments_.size(); ++specseg ) {
		std::string segname( "Seg_" + utility::to_string( specseg ) );
		//SilentEnergy segname_e( segname, specseg, 1, 10 );
		//silent_Es_.push_back( segname_e );
		utility::vector1< core::Size > segment;
		for ( core::Size i = spec_segments_[specseg].first; i < spec_segments_[specseg].second; ++i ) segment.push_back( i );
		compute_metrics_for_residue_subset( segname, separate_out_constraints, calc_pose, segment );
	}

	//if comparison to native is requested
	if ( native_comparison_ ) {
		native_comp_->compare_to_native( calc_pose, native_compare_calculators_, sfxn_, silent_Es_ );
	}

	//if repack without lig requested
	if ( repack_no_lig_ ) {

		core::pose::PoseOP rnlpose( new core::pose::Pose( pose ) );
		protocols::enzdes::RepackLigandSiteWithoutLigandMover rnl_mov( sfxn_, true );
		if ( enzcst_io_ ) rnl_mov.set_cstio( enzcst_io_ );
		rnl_mov.apply( *rnlpose );
		for ( core::Size i =1; i<= rnl_mov.silent_Es().size(); ++i ) silent_Es_.push_back( rnl_mov.silent_Es()[i] );

		if ( keep_rnl_pose_ ) rnl_pose_ = rnlpose;
	}
} //examine pose


void
EnzdesScorefileFilter::initialize_value_evaluators_from_file( std::string const & filename )
{

	auto map_it( evaluator_map_.find( filename ) );
	if ( map_it != evaluator_map_.end() ) return; //means this has already been done

	evaluator_map_.insert( std::pair< std::string, std::map< std::string, ValueEvaluator > >( filename, std::map< std::string, ValueEvaluator >() ) );
	map_it = evaluator_map_.find( filename );

	TR << "EnzdesScorefileFilter initializing filter params from file " << filename << "..." << std::endl;

	utility::io::izstream filedata( filename.c_str() );
	utility::vector1< std::string > tokens;
	std::string line("");
	if ( !filedata ) utility_exit_with_message("File " + filename + " couldn't be opened by EnzdesScorefileFilter.");

	while ( !filedata.eof() ) {
		getline(filedata,line);

		tokens.clear(); tokens.push_back(""); //weird utilvec1 copy behaviour makes this necessary
		tokens = utility::split( line );
		if ( tokens.size() < 1 ) continue;

		if ( tokens[1] == "req" ) {
			if ( tokens.size() < 5 ) utility_exit_with_message("Could not initialize filter params from line '" + line + "' because it was too short. Check your file format.");
			if ( tokens[3] != "value" ) {
				TR << "Warning: instruction '" << tokens[3] << "' in filter requirements file could not be understand, line will be ignored." << std::endl;
				continue;
			}

			std::string param_name( tokens[2] );
			core::Real param_cutoff( core::Real(atof(tokens[5].c_str()) ) );
			ValueEvaluator::CompareMode mode(ValueEvaluator::SMALLER);
			if ( tokens[4] == "<" ) mode = ValueEvaluator::SMALLER;
			else if ( tokens[4] == "=" ) mode = ValueEvaluator::EQUALS;
			else if ( tokens[4] == ">" ) mode = ValueEvaluator::LARGER;
			else utility_exit_with_message("Comparison mode " + tokens[4] + " for requirement " + param_name + " was not understood, needs to be either '>', '<' or '='");
			map_it->second.insert( std::pair<std::string, ValueEvaluator>(param_name, ValueEvaluator(mode,param_cutoff )) );
			TR << "Instantiated " << param_name << " requirement with cutoff " << param_cutoff << " and compare mode " << mode << "..." << std::endl;

		} //if we have a line in proper format

	} //file reading
	TR << " ...finished reading " << filename << "." << std::endl;
	filedata.close();
}

void
EnzdesScorefileFilter::compute_metrics_for_residue_subset(
	std::string sub_name,
	bool separate_out_constraints,
	core::pose::Pose const & calc_pose,
	utility::vector1< core::Size > const & res_subset ) const
{

	using namespace core::io::silent;

	 for ( auto const & relevant_scoreterm : relevant_scoreterms_ ) {

		std::string sco_name = sub_name + "_" +  relevant_scoreterm ;
		int width = std::max( 10, (int) sco_name.length() + 3 );

		SilentEnergy new_se;
		if ( relevant_scoreterm == "all_cst" ) {
			core::Real value(0.0);
			for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) value += enzutil::sum_constraint_scoreterms(calc_pose, res_subset[ii]);
			new_se = SilentEnergy ( sco_name, value , 1, width);
		} else if ( separate_out_constraints && ( relevant_scoreterm == "total_score" ) ) {
			core::Real desired_value(0.0);
			for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) desired_value += calc_pose.energies().residue_total_energies( res_subset[ii] )[ core::scoring::score_type_from_name( relevant_scoreterm ) ] - enzutil::sum_constraint_scoreterms(calc_pose, res_subset[ii] );
			new_se = SilentEnergy(sco_name, desired_value, 1, width);
		} else {
			core::Real value(0.0);
			for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) value += calc_pose.energies().residue_total_energies( res_subset[ii] )[ core::scoring::score_type_from_name( relevant_scoreterm ) ];
			value *= calc_pose.energies().weights()[ core::scoring::score_type_from_name( relevant_scoreterm ) ];
			new_se = SilentEnergy ( sco_name, value, 1 ,width);
		}

		silent_Es_.push_back( new_se );
	}//loop over relevant scoreterms


	//if there are calculators that need to be evaluated for this residue, let's do that now
	//note: section still under development, right now only calculators that return reals are supported
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > >::const_iterator res_calc_it = residue_calculators_.find( res_subset[1] );
	if ( res_calc_it != residue_calculators_.end() ) {

		utility::vector1< std::pair< std::string, std::string > > calculators_this_res = res_calc_it->second;
		for ( auto & calculators_this_re : calculators_this_res ) {

			std::string res_calc_name = sub_name + "_" + calculators_this_re.first;
			int width = std::max( 10, (int) res_calc_name.length() + 3 );

			core::Real calc_value(0.0);

			basic::MetricValue< core::Real > mval_real;
			basic::MetricValue< utility::vector1< core::Size > >mval_sizevec;
			basic::MetricValue< utility::vector1< core::Real > >mval_realvec;

			//following lines fairly hacky, but don't know a better solution at the moment
			if ( ( calculators_this_re.first == "hbond_pm") || ( calculators_this_re.first == "burunsat_pm") || ( calculators_this_re.first == "NLconts_pm" ) ) {
				calc_pose.metric( calculators_this_re.first, calculators_this_re.second, mval_sizevec );
				for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) calc_value += mval_sizevec.value()[ res_subset[ii] ];
			} else if ( calculators_this_re.first == "nlsurfaceE_pm" ) {
				calc_pose.metric( calculators_this_re.first, calculators_this_re.second, mval_realvec );
				for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) calc_value += mval_realvec.value()[ res_subset[ii] ];
			} else if ( (calculators_this_re.first == "pstat_pm") || (calculators_this_re.first == "nlpstat_pm" ) ) {
				calc_pose.metric( calculators_this_re.first, calculators_this_re.second, mval_realvec );
				core::Real pstat_sum(0.0);
				for ( core::Size ii =1; ii <= res_subset.size(); ++ii ) pstat_sum += mval_realvec.value()[ res_subset[ii] ];
				calc_value = pstat_sum / res_subset.size();
			} else {
				calc_pose.metric( calculators_this_re.first, calculators_this_re.second, mval_real );
				//for( core::Size ii =1; ii <= res_subset.size(); ++ii ) calc_value += mval_realvec.value()[ res_subset[ii] ];
				calc_value = mval_real.value();
			}
			//std::cerr << " hehe, just executed pose metric for " << calc_it->first << " calculator   ";

			SilentEnergy new_se( res_calc_name, calc_value, 1, width);
			silent_Es_.push_back( new_se );

		} // for calculators this res

	}// if calculators for this res
	//else std::cerr << "resi " << *res_it << " has no calcs." << std::endl;
} //compute_metrics_for_residue_subset


/// @brief function to setup residue specific calculators for all the catalytic residues and the ligand
void
EnzdesScorefileFilter::setup_pose_metric_calculators( core::pose::Pose const & pose, bool separate_out_constraints ) const
{

	using namespace core::pose::metrics;
	residue_calculators_.clear();
	native_compare_calculators_.clear();

	//calculators setup by this function ( under development )
	//1. packstat for the complete pose (with and without ligand)
	//2. SASA and interface delta for ligand
	//3. number of Hbonds for catalytic residues + ligand
	//4. number of buried polars for catalytic residues + ligand
	//5. packstat for every catalytic residue
	//6. further down the road maybe some stuff to compare native and designed structure

	std::string hbond_calc_name = "hbond_pm";
	std::string burunsat_calc_name = "burunsat_pm";
	std::string packstat_calc_name = "pstat_pm";
	std::string noligpackstat_calc_name = "nlpstat_pm";
	std::string nonlocalcontacts_calc_name = "NLconts_pm";
	std::string surface_calc_name = "nlsurfaceE_pm";
	std::string charge_calc_name = "charges_pm";

	//before starting, we make sure that all necessary calculators are instantiated
	if ( !CalculatorFactory::Instance().check_calculator_exists( hbond_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP hb_calc( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( burunsat_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP burunsat_calc( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", hbond_calc_name) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( packstat_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP pstat_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, pstat_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( noligpackstat_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP noligpstat_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator( true ) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( noligpackstat_calc_name, noligpstat_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( nonlocalcontacts_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP nlcontacts_calc( new protocols::toolbox::pose_metric_calculators::NonlocalContactsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( nonlocalcontacts_calc_name, nlcontacts_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( surface_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP surface_calc( new protocols::toolbox::pose_metric_calculators::SurfaceCalculator( true ) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( surface_calc_name, surface_calc );
	}
	if ( !CalculatorFactory::Instance().check_calculator_exists( charge_calc_name ) ) {
		core::pose::metrics::PoseMetricCalculatorOP charge_calc( new protocols::toolbox::pose_metric_calculators::ChargeCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( charge_calc_name, charge_calc );
	}


	//first general pose calculators

	utility::vector1< std::pair< std::string, std::string > > total_pose_calculators;

	if ( !no_packstat_calc_ ) {
		total_pose_calculators.push_back( std::pair< std::string, std::string > ( packstat_calc_name, "total_packstat") );
		total_pose_calculators.push_back( std::pair< std::string, std::string > ( noligpackstat_calc_name, "total_packstat") );
	}

	total_pose_calculators.push_back( std::pair< std::string, std::string > ( burunsat_calc_name, "all_bur_unsat_polars") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > ( hbond_calc_name, "all_Hbonds") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > ( nonlocalcontacts_calc_name, "total_nlcontacts") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > ( surface_calc_name, "total_surface") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > (charge_calc_name, "total_charge") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > (charge_calc_name, "total_pos_charges") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > (charge_calc_name, "total_neg_charges") );
	total_pose_calculators.push_back( std::pair< std::string, std::string > ( "seq_recovery", "seq_recovery") );

	residue_calculators_.insert( std::pair< Size, utility::vector1< std::pair< std::string, std::string > > > ( 0, total_pose_calculators ) );

	//if native compare is requested
	if ( native_comparison_ ) {
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( burunsat_calc_name, "all_bur_unsat_polars") );
		if ( !basic::options::option[basic::options::OptionKeys::enzdes::no_packstat_calculation] ) {
			native_compare_calculators_.push_back( std::pair< std::string, std::string > ( packstat_calc_name, "total_packstat") );
		}
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( hbond_calc_name, "all_Hbonds") );
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( nonlocalcontacts_calc_name, "total_nlcontacts") );
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( surface_calc_name, "total_surface") );
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( charge_calc_name, "total_charge") );
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( charge_calc_name, "total_pos_charges") );
		native_compare_calculators_.push_back( std::pair< std::string, std::string > ( charge_calc_name, "total_neg_charges") );
	}

	//general pose calculators done, now onto the residue specific calculators
	//note: all residue calculators will also be used to calculate stuff for the special segments
	//they'll be added to the map under the seqpos of the first residue in the
	//special segment
	std::set< core::Size > spec_seg_first_res;
	utility::vector1< core::Size > scoreres( enzutil::catalytic_res( pose ) );

	for ( core::Size speccount =1; speccount <= spec_segments_.size(); ++speccount ) {
		spec_seg_first_res.insert( spec_segments_[speccount].first );
		scoreres.push_back(  spec_segments_[speccount].first );
	}

	//detect all protein chains
	std::set< core::Size > protein_chains;
	for ( utility::vector1< core::Size >::const_iterator vecit( scoreres.begin() );  vecit != scoreres.end(); ++vecit ) {
		if ( pose.residue_type( *vecit ).is_protein() ) {
			protein_chains.insert(pose.chain( *vecit ));
		}
	}

	for ( utility::vector1< core::Size >::const_iterator vecit( scoreres.begin() );  vecit != scoreres.end(); ++vecit ) {

		//if( ! (enzutil::is_catalytic_seqpos( pose, i ) || (spec_seg_first_res.find( i ) != spec_seg_first_res.end()) ) ) continue;

		utility::vector1< std::pair< std::string, std::string > > calculators_this_res;

		//first a couple of ligand specific calculators ( interface SASA and interface energetics)
		if ( pose.residue_type( *vecit ).is_ligand() ) {
			Size lig_chain = pose.chain( *vecit );
			std::string lig_ch_string = utility::to_string( lig_chain );
			for ( unsigned long prot_chain : protein_chains ) {
				std::string prot_ch_string = utility::to_string( prot_chain );
				if ( lig_chain == prot_chain ) { utility_exit_with_message( "WTF?!? ligand and residue 1 are on the same chain... " );}

				std::string lig_interface_neighbor_calc_name = "neighbor_def_" + prot_ch_string + "_" + lig_ch_string;
				std::string lig_dsasa_calc_name = "dsasa_" + prot_ch_string + "_" + lig_ch_string;
				std::string lig_interface_e_calc_name = "interf_E_" + prot_ch_string + "_" + lig_ch_string;

				if ( !CalculatorFactory::Instance().check_calculator_exists( lig_interface_neighbor_calc_name ) ) {
					PoseMetricCalculatorOP lig_neighbor_calc( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( prot_chain, lig_chain ) );
					CalculatorFactory::Instance().register_calculator( lig_interface_neighbor_calc_name, lig_neighbor_calc );
				}
				if ( !CalculatorFactory::Instance().check_calculator_exists( lig_dsasa_calc_name ) ) {
					PoseMetricCalculatorOP lig_dsasa_calc( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator( prot_chain, lig_chain ) );
					CalculatorFactory::Instance().register_calculator( lig_dsasa_calc_name, lig_dsasa_calc );
				}
				if ( !CalculatorFactory::Instance().check_calculator_exists( lig_interface_e_calc_name ) ) {
					utility::vector1<ScoreType> score_types_to_ignore;
					if ( separate_out_constraints ) {
						score_types_to_ignore.push_back( ScoreType( coordinate_constraint ) );
						score_types_to_ignore.push_back( ScoreType( atom_pair_constraint ) );
						score_types_to_ignore.push_back( ScoreType( angle_constraint ) );
						score_types_to_ignore.push_back( ScoreType( dihedral_constraint ) );
					}
					PoseMetricCalculatorOP lig_interf_E_calc( new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( lig_interface_neighbor_calc_name, score_types_to_ignore ) );
					CalculatorFactory::Instance().register_calculator( lig_interface_e_calc_name, lig_interf_E_calc );
				}

				calculators_this_res.push_back( std::pair< std::string, std::string > ( lig_interface_e_calc_name, "weighted_total") );
				calculators_this_res.push_back( std::pair< std::string, std::string > ( lig_dsasa_calc_name, "frac_ch2_dsasa") );
			}
		} //ligand specific calculators set up

		//now calculators for every residue, for starters number of Hbonds and number of buried polars
		calculators_this_res.push_back( std::pair< std::string, std::string > ( hbond_calc_name, "residue_Hbonds") );
		calculators_this_res.push_back( std::pair< std::string, std::string > ( burunsat_calc_name, "residue_bur_unsat_polars") );

		//and for protein residues, we also want to know the packstat
		if ( pose.residue_type( *vecit ).is_protein() && ( !basic::options::option[basic::options::OptionKeys::enzdes::no_packstat_calculation] ) ) {
			calculators_this_res.push_back( std::pair< std::string, std::string > ( packstat_calc_name, "residue_packstat") );
			calculators_this_res.push_back( std::pair< std::string, std::string > ( noligpackstat_calc_name, "residue_packstat") );

			//for loop residues also non local contacts and surface E
			if ( spec_seg_first_res.find( *vecit ) != spec_seg_first_res.end() ) {
				calculators_this_res.push_back( std::pair< std::string, std::string > ( nonlocalcontacts_calc_name, "residue_nlcontacts") );
				calculators_this_res.push_back( std::pair< std::string, std::string > ( surface_calc_name, "residue_surface") );
			}
		}
		//debug
		//for( utility::vector1< std::pair< std::string, std::string > >::iterator calc_it = calculators_this_res.begin();
		//   calc_it != calculators_this_res.end(); calc_it++ ){
		// std::cerr << "calculator " << calc_it->first << " created for residue " << *res_it << std::endl;
		//}
		//debug over

		residue_calculators_.insert( std::pair< Size, utility::vector1< std::pair< std::string, std::string > > > ( *vecit, calculators_this_res ) );

	} // loop over all (catalytic) residues
} //setup_pose_metric_calculators


core::pose::PoseOP
EnzdesScorefileFilter::rnl_pose()
{
	return rnl_pose_;
}

/// @brief clear rnl pose to save some memory
void
EnzdesScorefileFilter::clear_rnl_pose()
{
	rnl_pose_ = nullptr;
}

ResidueConformerFilter::ResidueConformerFilter()
: Filter(),
	restype_(nullptr), seqpos_(0),
	desired_conformer_(0), max_rms_(0.5),
	lig_conformer_builder_(/* NULL */)
{
	relevant_atom_indices_.clear();
}

ResidueConformerFilter::ResidueConformerFilter( ResidueConformerFilter const & ) = default;

ResidueConformerFilter::~ResidueConformerFilter()= default;

bool
ResidueConformerFilter::apply( core::pose::Pose const & pose ) const
{
	core::Size current_conformer( this->get_current_conformer( pose ) );
	if ( desired_conformer_ == 0 ) {
		TR << "ResidueConformerFilter did not have a desired conformer set. Apply gets passed by default, conformer in pose is " << current_conformer << std::endl;
		return true;
		//throw utility::excn::EXCN_RosettaScriptsOption("For ResidueConformerFilter, desired conformer needs to be set if it is used in actual filtering.");
	}
	if ( current_conformer == desired_conformer_ ) return true;
	else return false;
}


void
ResidueConformerFilter::report( std::ostream & stream, core::pose::Pose const & pose  ) const
{
	stream << "ResidueConformerFilter detected conformation " << this->get_current_conformer( pose ) << " in pose. " << std::endl;
}

core::Real
ResidueConformerFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real to_return = this->get_current_conformer( pose );
	return to_return;
}

void
ResidueConformerFilter::parse_my_tag(utility::tag::TagCOP tag , basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	relevant_atom_indices_.clear();

	if ( tag->hasOption("restype") ) {
		std::string resname =  tag->getOption<std::string>( "restype","" );
		restype_ = & core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( resname );
	} else throw utility::excn::EXCN_RosettaScriptsOption("For ResidueConformerFilter, the desired residue type needs to be specified.");

	if ( tag->hasOption("relevant_atoms") ) {
		utility::vector1< std::string > const atom_vec( utility::string_split( tag->getOption< std::string >("relevant_atoms", "" ), ',' ) );
		for ( core::Size i(1); i <= atom_vec.size(); ++i ) {
			relevant_atom_indices_.push_back( restype_->atom_index( atom_vec[i] ) );
		}
	} else { //use all heavy atoms
		for ( core::Size i(1); i <= restype_->nheavyatoms(); ++i ) relevant_atom_indices_.push_back( i );
	}

	if ( tag->hasOption("desired_conformer") ) desired_conformer_ = tag->getOption<core::Size>( "desired_conformer",0 );
	if ( tag->hasOption("seqpos") ) seqpos_ = tag->getOption<core::Size>( "seqpos",0 );
	if ( tag->hasOption("max_rms") ) max_rms_ = tag->getOption<core::Real>( "max_rms",0.0 );

	this->initialize_internal_data();
}

void
ResidueConformerFilter::initialize_internal_data()
{
	match::downstream::LigandConformerBuilderOP lcbuilder( new match::downstream::LigandConformerBuilder() );
	lcbuilder->set_idealize_conformers( false ); //not necessary if we're not matching
	lcbuilder->set_rmsd_unique_cutoff( max_rms_ );
	core::conformation::ResidueOP dummy_res( new core::conformation::Residue( *restype_, true ) );

	lcbuilder->initialize_from_residue( 1, 2, 3, 1, 2, 3, *dummy_res ); //all atoms set to arbitrary 1,2,3, since we're not matching
	lcbuilder->determine_redundant_conformer_groups( relevant_atom_indices_ ); //maybe this isn't even necessary

	lig_conformer_builder_ = lcbuilder;
}

/// @details
/// this might actually need to be implemented differently
/// for residue types that get their rotamers/conformers from the dunbrack
/// library, as this is bb dependent. this filter however assumes that the
/// conformer library for a given residue type is constant.
core::Size
ResidueConformerFilter::get_current_conformer( core::pose::Pose const & pose ) const
{
	core::Size seqpos_to_check( seqpos_);

	if ( seqpos_to_check == 0 ) { //this means seqpos didn't get set by the users, and we have to look in the pose
		utility::vector1< core::Size > possible_seqpos;
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue( i ).name3() == restype_->name3() ) possible_seqpos.push_back( i );
		}
		if ( possible_seqpos.size() == 0 ) {
			throw utility::excn::EXCN_RosettaScriptsOption("In ResidueConformerFilter, no residue of the desired type was found in the pose.");
		}
		if ( possible_seqpos.size() > 1 ) {
			std::cerr <<"In ResidueConformerFilter, more than one possible seqpos to check has been found. Taking only the first one (" << possible_seqpos[1] <<"), but this is probably not what you wanted." << std::endl;
		}
		seqpos_to_check = possible_seqpos[1];
	}
	core::conformation::Residue const & res_to_check( pose.residue( seqpos_to_check ) );

	return lig_conformer_builder_->assign_conformer_group_to_residue( res_to_check, relevant_atom_indices_ );

}


filters::FilterOP
DiffAtomSasaFilterCreator::create_filter() const { return filters::FilterOP( new DiffAtomSasaFilter ); }

std::string
DiffAtomSasaFilterCreator::keyname() const { return "DiffAtomBurial"; }

filters::FilterOP
EnzScoreFilterCreator::create_filter() const { return filters::FilterOP( new EnzScoreFilter ); }

std::string
EnzScoreFilterCreator::keyname() const { return "EnzScore"; }

filters::FilterOP
LigBurialFilterCreator::create_filter() const { return filters::FilterOP( new LigBurialFilter ); }

std::string
LigBurialFilterCreator::keyname() const { return "LigBurial"; }

filters::FilterOP
LigDSasaFilterCreator::create_filter() const { return filters::FilterOP( new LigDSasaFilter ); }

std::string
LigDSasaFilterCreator::keyname() const { return "DSasa"; }

filters::FilterOP
LigInterfaceEnergyFilterCreator::create_filter() const { return filters::FilterOP( new LigInterfaceEnergyFilter ); }

std::string
LigInterfaceEnergyFilterCreator::keyname() const { return "LigInterfaceEnergy"; }

filters::FilterOP
RepackWithoutLigandFilterCreator::create_filter() const { return filters::FilterOP( new RepackWithoutLigandFilter ); }

std::string
RepackWithoutLigandFilterCreator::keyname() const { return "RepackWithoutLigand"; }

filters::FilterOP
EnzdesScorefileFilterCreator::create_filter() const { return filters::FilterOP( new EnzdesScorefileFilter ); }

std::string
EnzdesScorefileFilterCreator::keyname() const { return "EnzdesScorefileFilter"; }

filters::FilterOP
ResidueConformerFilterCreator::create_filter() const { return filters::FilterOP( new ResidueConformerFilter ); }

std::string
ResidueConformerFilterCreator::keyname() const { return "ResidueConformerFilter"; }


} // enzdes
} // protocols
