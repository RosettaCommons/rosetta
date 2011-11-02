// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/dock_design_filters.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/dock_design_filters.hh>
#include <protocols/protein_interface_design/dock_design_filter_creators.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/moves/RigidBodyMover.hh>

//#include <protocols/moves/ResidueMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/MetricValue.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/chemical/AtomType.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//Objectxxxx header
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility Headers

// Unit Headers
#include <protocols/protein_interface_design/movers/ddG.hh>
#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <map>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::fmt;

static numeric::random::RandomGenerator RG( 140845 ); // <- Magic number, do not change it!!!

namespace protocols {
namespace protein_interface_design {

static basic::Tracer TR( "protocols.protein_interface_design.dock_design_filters" );
using core::pose::Pose;

DdgFilter::DdgFilter() :
	Filter( "Ddg" ),
	ddg_threshold_( -15.0 ),
	scorefxn_( NULL ),
	rb_jump_( 1 ),
	repeats_( 1 ),
	symmetry_(false),
	repack_( true )
{}

void
DdgFilter::repack( bool const repack )
{
	repack_ = repack;
}

bool
DdgFilter::repack() const
{
	return repack_;
}

bool
DdgFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const pose_ddg( compute( pose ) );
	TR<<"ddg is "<<pose_ddg<<" ";
	if( pose_ddg <= ddg_threshold_ ) {
		TR<<"passing"<<std::endl;
		return true;
	}
	TR<<"failing"<<std::endl;
	TR.flush();
	return false;
}

void
DdgFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	out<<"ddg "<<pose_ddg<<'\n';
}

core::Real
DdgFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	return( pose_ddg );
}

core::Size
DdgFilter::repeats() const
{
	return( repeats_ );
}

void
DdgFilter::repeats( core::Size const repeats )
{
	repeats_ = repeats;
}

core::Real
DdgFilter::compute( core::pose::Pose const & pose ) const {
	if( repack() ){
		protocols::protein_interface_design::movers::ddG ddg( scorefxn_, rb_jump_, symmetry_ );
		core::Real average( 0.0 );
		for( core::Size i = 1; i<=repeats_; ++i ){
			ddg.calculate( pose );
			average += ddg.sum_ddG();
			ddg.report_ddG( TR );
		}
		return average / (core::Real)repeats_;
	}
	else{
		if( repeats() > 1 && !repack() )
			utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );
		using namespace protocols::moves;

		ScoreTypeFilter const stf( scorefxn_, core::scoring::total_score, 10000/*threshold*/ );
		core::pose::Pose split_pose( pose );
		RigidBodyTransMoverOP translate( new RigidBodyTransMover( split_pose, rb_jump_ ) );
		translate->step_size( 1000.0 );
		translate->apply( split_pose );
		core::Real const bound_energy( stf.compute( pose ));
		core::Real const unbound_energy( stf.compute( split_pose ));
		core::Real const dG( bound_energy - unbound_energy );
		return( dG );
	}
}

ScoreTypeFilter::ScoreTypeFilter( core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold ) : Filter( "ScoreType" ) {
	score_type_ = score_type;
	score_type_threshold_ = score_type_threshold;
	scorefxn_ = scorefxn->clone();
}

DdgFilter::DdgFilter( core::Real const ddg_threshold, core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump/*=1*/, core::Size const repeats/*=1*/, bool const symmetry /*=false*/ ) :
	Filter("Ddg" ),
	repack_( true )
{
	ddg_threshold_ = ddg_threshold;
	scorefxn_ = scorefxn->clone();
	rb_jump_ = rb_jump;
	repeats_ = repeats;
	symmetry_ = symmetry;
}

bool
ScoreTypeFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	TR<<"score "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" is "<<score<<". ";
	if( score <= score_type_threshold_ ) {
		TR<<"passing." << std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
ScoreTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"Weighted score of "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" "<<compute( pose )<<'\n';
}

core::Real
ScoreTypeFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ScoreTypeFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose;
	using namespace core::scoring;

	PoseOP in_pose = new Pose( pose );

	// make sure that scoring weights are compatible with pose's residue type set
	// check centroid case
	if( ( (*scorefxn_)[fa_rep] == 0.0 && (*scorefxn_)[fa_atr] == 0.0 ) // full atom terms are off
				&& ( (*scorefxn_)[interchain_vdw] > 0.0 || (*scorefxn_)[vdw] > 0.0)  ) // a centroid term is on
		{
			if( in_pose->is_fullatom() ) { // but pose is full atom
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::CENTROID );
		}
	}
	else { // full atom case
		if( in_pose->is_centroid() ) { // but pose is centroid
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::FA_STANDARD );
		}
	}

	(*scorefxn_)( *in_pose );
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies( *in_pose );
	core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
	core::Real const score( in_pose->energies().total_energies()[ ScoreType( score_type_ ) ]);
	if( score_type_ == total_score ) return( score );
	core::Real const weighted_score( weight * score );
	return( weighted_score );
}

void
AlaScan::scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

AlaScan::AlaScan( bool const chain1, bool const chain2, core::Size const repeats, core::Real const dist, core::scoring::ScoreFunctionCOP scorefxn, core::Size const jump=1, bool const symmetry=false ) : Filter( "AlaScan" ),
		chain1_( chain1 ),
		chain2_( chain2 ),
		repeats_( repeats ),
		distance_threshold_( dist ),
		jump_( jump ),
		symmetry_( symmetry ),
		repack_( true )
{
	if ( symmetry_ ) scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( scorefxn );
	else scorefxn_ = scorefxn->clone();
}

bool
AlaScan::repack() const
{
	return repack_;
}

void
AlaScan::repack( bool const repack )
{
	repack_ = repack;
}

core::Real
AlaScan::ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi ) const
{
	if( !const_pose.residue( resi ).is_protein() ){
		TR<<"WARNING: Non-protein residue "<< resi<<" was requested for ala-scan. Returning 0"<<std::endl;
		return 0.0;
	}
	core::Size const rb_jump( jump_ );
	core::pose::Pose pose( const_pose );

	DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn_, rb_jump, 1 /*repeats*/ );
	if( repack() )
		TR<<"Energy calculations are carried out with repacking in the bound and unbound states (ddG)\n";
	else
		TR<<"Energy calculations are carried out without repackign in the bound and unbound states (dG)\n";
	ddg_filter.repack( repack() );
  ScoreTypeFilter const energy_filter( scorefxn_, core::scoring::total_score, 0 );

	utility::vector1< bool > allowed_aas;
	allowed_aas.assign( chemical::num_canonical_aas, false );
	allowed_aas[ chemical::aa_ala ] = true;
	using namespace core::pack::task;

	PackerTaskOP task = TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line().or_include_current( true );
	for( core::Size resj=1; resj<=pose.total_residue(); ++resj ){
		if( resi == resj )
			task->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
		else
			task->nonconst_residue_task( resj ).prevent_repacking();
	}
	core::pack::pack_rotamers( pose, *scorefxn_, task );
	core::Real accumulate_ddg = 0;

	for( core::Size r=1; r<=repeats_; ++r )
		accumulate_ddg += (rb_jump==0 ? energy_filter.compute( pose ) : ddg_filter.compute( pose ) );;
	core::Real const mut_ddg( accumulate_ddg / repeats_ );

	TR.flush();
	return( mut_ddg );
}

void
AlaScan::report( std::ostream & out, core::pose::Pose const & const_pose ) const
{
	if ( symmetry_ ) {
		report_symmetry( out, const_pose );
		return;
	}

	core::Size const rb_jump( jump_ );
	core::pose::Pose pose( const_pose );

	core::kinematics::FoldTree const fold_tree = pose.conformation().fold_tree();

	core::Size upstream_jump_res, downstream_jump_res;
	upstream_jump_res = fold_tree.upstream_jump_residue( jump_ );
	downstream_jump_res = fold_tree.downstream_jump_residue( jump_ );

	core::Size const chain_begin( chain1_ ? 1 : downstream_jump_res );
	core::Size const chain_end  ( chain2_ ? pose.total_residue() : upstream_jump_res );

	protocols::scoring::Interface interface_obj;
	interface_obj.jump( rb_jump == 0 ? 1 : rb_jump ); // 0 plays badly with interface obj.
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	interface_obj.distance( distance_threshold_ );
	interface_obj.calculate( pose );

	DdgFilter const ddg_filter( 100/*ddg_threshold*/, scorefxn_, rb_jump, 1 /*repeats*/ );
  ScoreTypeFilter const energy_filter( scorefxn_, core::scoring::total_score, 0 );

	core::Real accumulate_ddg( 0 );
	for( core::Size r=1; r<=repeats_; ++r )
		accumulate_ddg += (rb_jump==0 ? energy_filter.compute( const_pose ) : ddg_filter.compute( const_pose ) );

	core::Real const wt_ddg( accumulate_ddg / repeats_ );
	for( core::Size resi = chain_begin; resi <= chain_end; ++resi ){
		if( !pose.residue( resi ).is_protein() ) continue;
		if( interface_obj.is_interface( resi ) ){
			core::Real const mut_ddg( ddG_for_single_residue( const_pose, resi ) );
			core::Real const diff_ddg( mut_ddg - wt_ddg );

			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			char const chain( pose_info->chain( resi ) );
			core::Size const number( pose_info->number( resi ) );
			std::string const res_type( const_pose.residue( resi ).name3() );
			out<<" "<<res_type<<" "<<number<<" "<<chain<<" : "<< F (9,4,diff_ddg)<<'\n';
		}
	}
	out<<std::endl;
}

void
AlaScan::report_symmetry( std::ostream & out, core::pose::Pose const & const_pose ) const
{
	core::pose::Pose pose( const_pose );

	assert( core::pose::symmetry::is_symmetric( pose ));
  core::conformation::symmetry::SymmetricConformation & symm_conf (
        dynamic_cast<core::conformation::symmetry::SymmetricConformation & > ( pose.conformation()) );

	protocols::scoring::Interface interface_obj(1);
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	interface_obj.distance( distance_threshold_ );
	interface_obj.calculate( pose );

	DdgFilter const ddg( 100/*ddg_threshold*/, scorefxn_, 1, 1 /*repeats*/, true /*symmetry*/ );
	core::Real accumulate_ddg( 0 );
	for( core::Size r=1; r<=repeats_; ++r )
		accumulate_ddg += ddg.compute( const_pose );
	core::Real const wt_ddg( accumulate_ddg / repeats_ );

	//core::Real const wt_ddg( ddg.compute( const_pose ) );
	utility::vector1< bool > allowed_aas;
	allowed_aas.assign( chemical::num_canonical_aas, false );
	allowed_aas[ chemical::aa_ala ] = true;
	for( core::Size resi = 1; resi <= pose.total_residue(); ++resi ){
		if ( !symm_conf.Symmetry_Info()->bb_is_independent(resi) ) continue;
		if( !pose.residue( resi ).is_protein() ) continue;
		if( interface_obj.is_interface( resi ) ){
			using namespace core::pack::task;

			PackerTaskOP task = TaskFactory::create_packer_task( pose );
			task->initialize_from_command_line().or_include_current( true );
			for( core::Size resj=1; resj<=pose.total_residue(); ++resj ){
				if( !pose.residue( resi ).is_protein() ) continue;
				if( resi == resj )
					task->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
				else
					task->nonconst_residue_task( resj ).prevent_repacking();
			}
			core::pack::pack_rotamers( pose, *scorefxn_, task );
			accumulate_ddg = 0;
			for( core::Size r=1; r<=repeats_; ++r ) accumulate_ddg += ddg.compute( pose );
			core::Real const mut_ddg( accumulate_ddg / repeats_ );
			//core::Real const mut_ddg( ddg.compute( const_pose ) );
			core::Real const diff_ddg( mut_ddg - wt_ddg );

			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			//char const chain( pose_info->chain( resi ) );
			//core::Size const number( pose_info->number( resi ) );
			std::string const res_type( const_pose.residue( resi ).name3() );
			//out<<" "<<res_type<<" "<<number<<" "<<chain<<" : "<< F (9,4,diff_ddg)<<'\n';
			out<<" "<<res_type<<" "<< resi <<" : "<< F (9,4,diff_ddg)<<'\n';
			pose=const_pose;
		}
	}
}

InterfaceSasaFilter::InterfaceSasaFilter() :
	Filter( "Sasa" ),
	lower_threshold_( 0.0 ),
	hydrophobic_( false ),
	polar_( false ),
	jump_( 1 )
{}

InterfaceSasaFilter::InterfaceSasaFilter( core::Real const lower_threshold, bool const hydrophobic/*=false*/, bool const polar/*=false*/ ) :
	Filter( "Sasa" ),
	lower_threshold_( lower_threshold ),
	hydrophobic_( hydrophobic ),
	polar_( polar )
{}

InterfaceSasaFilter::~InterfaceSasaFilter(){}

FilterOP
InterfaceSasaFilter::clone() const{
	return new InterfaceSasaFilter( *this );
}

FilterOP
InterfaceSasaFilter::fresh_instance() const{
	return new InterfaceSasaFilter;
}

bool
InterfaceSasaFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ) );

	TR<<"sasa is "<<sasa<<". ";
	if( sasa >= lower_threshold_ ){
		TR<<"passing." <<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
InterfaceSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	out<<"Sasa= "<< sasa<<'\n';
}

core::Real
InterfaceSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	return( sasa );
}

void
InterfaceSasaFilter::jump( core::Size const jump )
{
	jump_ = jump;
}

core::Size
InterfaceSasaFilter::jump() const
{
	return jump_;
}

core::Real
InterfaceSasaFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose::metrics;
	using basic::MetricValue;
	using namespace core;
	using namespace protocols::moves;

	core::pose::Pose split_pose( pose );
	RigidBodyTransMoverOP translate( new RigidBodyTransMover( split_pose, jump() ) );
	translate->step_size( 1000.0 );
	translate->apply( split_pose );

	runtime_assert( !hydrophobic_ || !polar_ );
	if( !hydrophobic_ && !polar_ ){
		MetricValue< core::Real > mv_sasa;

		pose.metric( "sasa", "total_sasa", mv_sasa);
		core::Real const bound_sasa( mv_sasa.value() );
		split_pose.metric( "sasa", "total_sasa", mv_sasa );
		core::Real const unbound_sasa( mv_sasa.value() );
		core::Real const buried_sasa( unbound_sasa - bound_sasa );
		return( buried_sasa );
	}
	else{
		MetricValue< id::AtomID_Map< core::Real > > atom_sasa;
		pose.metric( "sasa_interface", "delta_atom_sasa", atom_sasa );
		core::Real polar_sasa( 0.0 ), hydrophobic_sasa( 0.0 );
		for( core::Size pos(1); pos<=pose.total_residue(); ++pos ){
			for( core::Size atomi( 1 ); atomi <= atom_sasa.value().n_atom( pos ); ++atomi ){
				core::Real const atomi_delta_sasa( atom_sasa.value()( pos, atomi ) );
				core::conformation::Residue const pos_rsd( pose.residue( pos ) );
				core::chemical::AtomType const atom_type( pos_rsd.atom_type( atomi ) );
				bool const is_polar( atom_type.is_donor() || atom_type.is_acceptor() || atom_type.is_polar_hydrogen() );
				if( is_polar ) polar_sasa += atomi_delta_sasa;
				else hydrophobic_sasa += atomi_delta_sasa;
			}
		}
		if( hydrophobic_ ) return hydrophobic_sasa;
		if( polar_ ) return polar_sasa;
	}
	utility_exit_with_message( "Execution should never have reached this point." );
	return( 0 );
}

bool
NeighborTypeFilter::apply( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	if( neighbors.size() == 0 ) return false;
	TR<<"neighbours of residue "<<pose.residue( target_residue_ ).name3()<<target_residue_<<": ";
	for( std::vector< core::Size >::const_iterator n_it=neighbors.begin(); n_it!=neighbors.end(); ++n_it ) {
		TR<<pose.residue( *n_it ).name3()<<*n_it<<" ";
	}
	TR<<std::endl;
	return true;
}

void
NeighborTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	out<<"neighbours of residue "<<pose.residue( target_residue_ ).name3()<<target_residue_<<": ";
	for( std::vector< core::Size >::const_iterator n_it=neighbors.begin(); n_it!=neighbors.end(); ++n_it ) {
		out<<pose.residue( *n_it ).name3()<<*n_it<<" ";
	}
	out<<'\n';
}

core::Real
NeighborTypeFilter::report_sm( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	return( neighbors.size() != 0 );
}

std::vector< core::Size >
NeighborTypeFilter::compute( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors;

	core::Size const chain2begin( pose.conformation().chain_begin( 2 ) );
	core::Size const residue_num( pose.total_residue() );

	core::Size const start( target_residue_ < chain2begin ? chain2begin : 1 );
	core::Size const end( target_residue_ < chain2begin ? residue_num : chain2begin-1 );
	core::conformation::Residue const res_target( pose.residue( target_residue_ ) );

	runtime_assert( target_residue_ <= residue_num );
	for( core::Size res=start; res<=end; ++res ) {
		core::conformation::Residue const resi( pose.residue( res ) );
		if( !residue_types_[ resi.aa() ] ) continue;

		core::Real const distance( res_target.xyz( res_target.nbr_atom() ).distance( resi.xyz( resi.nbr_atom() )) );
		if( distance <= distance_threshold_ )
			neighbors.push_back( res );
	}
	return neighbors;
}

bool
ResidueBurialFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	TR<<"Number of interface neighbors of residue "<<pose.residue( target_residue_ ).name3()<<target_residue_<<" is "<<count_neighbors<<std::endl;
	return( count_neighbors >= neighbors_ );
}

void
ResidueBurialFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	out<<"Number of interface neighbors of residue "<<pose.residue( target_residue_ ).name3()<<target_residue_<<" is "<<count_neighbors<<'\n';
}

core::Real
ResidueBurialFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const count_neighbors( compute( pose ) );

	return( count_neighbors );
}

/// @details counts the number of residues to target_residue_ across all chains in the pose, other than the one containing target_residue_
core::Size
ResidueBurialFilter::compute( core::pose::Pose const & pose ) const {
	core::Size chain( 1 );
	for( ; chain <= pose.conformation().num_chains(); ++chain )
		if( pose.conformation().chain_begin( chain ) <= target_residue_ && pose.conformation().chain_end( chain ) >= target_residue_ ) break;

	core::Size const chain_begin( pose.conformation().chain_begin( chain ) );
	core::Size const chain_end( pose.conformation().chain_end( chain ) );

TR<<"chain span "<<chain_begin<< " "<<chain_end<<std::endl;
	core::Size count_neighbors( 0 );
	core::conformation::Residue const res_target( pose.conformation().residue( target_residue_ ) );
	for( core::Size i=1; i<=pose.total_residue(); ++i ){
		if( i>=chain_begin && i<=chain_end ) continue;
		core::conformation::Residue const resi( pose.residue( i ) );
		core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
		if( distance <= distance_threshold_ ) ++count_neighbors;
	}
	return( count_neighbors);
}

bool
ResidueDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	TR<<"Distance between residues "<<pose.residue( res1_ ).name3()<<res1_<<" and "<<pose.residue( res2_ ).name3()<<res2_<<" is "<<distance<<std::endl;
	return( distance<=distance_threshold_ );
}

void
ResidueDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	out<<"Distance between residues "<<pose.residue( res1_ ).name3()<<res1_<<" and "<<pose.residue( res2_ ).name3()<<res2_<<" is "<<distance<<'\n';
}

core::Real
ResidueDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	return( distance );
}
core::Real
ResidueDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::conformation::Residue const res_res1( pose.conformation().residue( res1_ ) );
	core::conformation::Residue const res_res2( pose.conformation().residue( res2_ ) );
	core::Real const distance( res_res1.xyz( res_res1.nbr_atom() ).distance( res_res2.xyz( res_res2.nbr_atom() ) ) );
	return( distance );
}

/// @detailed a utility for docking_filter
core::Size
ResiduesInInterfaceFilter::compute( Pose const & pose ) const
{
	Size interface_counter( 0 );

	runtime_assert( rb_jump_ <= pose.num_jump() );

	core::pose::Pose temp_pose( pose );
	core::scoring::ScoreFunctionOP scorefxn(
											ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	(*scorefxn)(temp_pose);

	protocols::scoring::Interface interface(rb_jump_);
	interface.calculate( temp_pose );

	for (Size i = 1; i <= temp_pose.total_residue(); i++) {
		if ( !temp_pose.residue(i).is_protein() ) continue;
		if( interface.is_interface( i ) ) interface_counter++;
	}
	return( interface_counter );
}

bool
ResiduesInInterfaceFilter::apply( Pose const & pose ) const {
	Size const interface_res( compute( pose ));
	TR<<"There are "<<interface_res<<" residues in the interface.";
	if( interface_res <= residues_in_interface_threshold_ )
	{
		TR<<" Breaking out."<<std::endl;
		return( false );
	}
	else TR<<std::endl;
	return( true );
}

void
ResiduesInInterfaceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	Size const interface_res( compute( pose ));
	out<<"residues_in_interface "<<interface_res<<'\n';
}

core::Real
ResiduesInInterfaceFilter::report_sm( core::pose::Pose const & pose ) const {
	Size const interface_res( compute( pose ));
	return( (core::Real) interface_res );
}

bool
HbondsToResidueFilter::apply( Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );
	TR<<"found "<<hbonded_res<< " hbond to target residue " << resnum_;
	if( hbonded_res >= partners_ ) {
		TR << ". passing." << std::endl;
		return( true );
	}
	else {
		TR << ". failing." << std::endl;
		return( false );
	}
}

void
HbondsToResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );

	out<<"Number of residues hbonded to "<<resnum_<< " is " << hbonded_res <<'\n';
}

core::Real
HbondsToResidueFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );
	return( hbonded_res );
}

core::Size
HbondsToResidueFilter::compute( Pose const & pose ) const {
	typedef core::Size Size;
	typedef core::Real Real;

	core::pose::Pose temp_pose( pose );
	core::scoring::ScoreFunctionOP scorefxn(
											ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	(*scorefxn)(temp_pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( temp_pose );

	core::Size const chain2begin( temp_pose.conformation().chain_begin( 2 ) );
	core::Size partner_begin, partner_end;
	if( resnum_ >= chain2begin ) {
		partner_begin = 1; partner_end = chain2begin-1;
	}
	else {
		partner_begin = chain2begin; partner_end = temp_pose.total_residue();
	}
	std::set<Size> binders;
	for( Size i=partner_begin; i<=partner_end; ++i ) binders.insert( i );

	std::list< Size> hbonded_res( hbonded( temp_pose, resnum_, binders, backbone_, sidechain_, energy_cutoff_ ) );

	return( hbonded_res.size() );
}


EnergyPerResidueFilter::EnergyPerResidueFilter( core::Size const resnum, core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const threshold, bool const whole_interface, core::Size const rb_jump, core::Real const interface_distance_cutoff ) : Filter( "EnergyPerResidue" ), resnum_( resnum ), score_type_( score_type ), threshold_( threshold ), 	whole_interface_ ( whole_interface ), rb_jump_ ( rb_jump ), interface_distance_cutoff_ ( interface_distance_cutoff )  {

		using namespace core::scoring;

		if( scorefxn ) scorefxn_ = new core::scoring::ScoreFunction( *scorefxn );
		if( score_type_ != total_score ) {
			core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
			scorefxn_->reset();
			scorefxn_->set_weight( score_type_, old_weight );

		}
	}

	EnergyPerResidueFilter::EnergyPerResidueFilter( EnergyPerResidueFilter const &init ) :
	//utility::pointer::ReferenceCount(),
	Filter( init ), resnum_( init.resnum_ ), score_type_( init.score_type_ ), threshold_( init.threshold_ ),
	whole_interface_ (init.whole_interface_), rb_jump_ (init.rb_jump_), interface_distance_cutoff_ ( init.interface_distance_cutoff_){
		using namespace core::scoring;
		if( init.scorefxn_ ) scorefxn_ = new core::scoring::ScoreFunction( *init.scorefxn_ );
	}


bool
EnergyPerResidueFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;

	if ( whole_interface_)
	{
		if ( pose.conformation().num_chains() < 2 ) {
			TR << "pose must contain at least two chains!" << std::endl;
			return false;

		}
		else {
			TR<<" \n \t --------- computing  ---------- \n \t-------- interface energies   --------\n \t \t------------------- \n" << std::endl;
			return true;
		}
	}
	else
	{


	core::Real const energy( compute( pose ) );
	TR<<"Scoretype "<<name_from_score_type( score_type_ )<<" for residue: " << pose.pdb_info()->number(resnum_)<<" " <<pose.residue( resnum_).name3() <<" is "<<energy<<". ";
	bool const pass( energy <= threshold_ );
	if( pass ) TR<<"passing."<<std::endl;
	else TR<<"failing."<<std::endl;

	return pass;
 	}
}


void
EnergyPerResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const

{
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

	if( whole_interface_ )
	{
		core::pose::Pose in_pose = pose;
		FArray1D_bool partner1_( in_pose.total_residue(), false );
		in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
		protocols::scoring::Interface interface_obj(rb_jump_);
		in_pose.update_residue_neighbors();
		interface_obj.distance( interface_distance_cutoff_ );
		interface_obj.calculate( in_pose );
		(*scorefxn_)( in_pose );

		out<<A(9, "chain")<<A( 9, "res")<<A( 9, "AA")<<A( 9, "total")<<A( 9, "contact")<<A( 9, "fa_atr")<<A( 9, "fa_rep")<<A( 9, "hb_bb_sc")<<A( 9, "hb_sc")<<A( 9, "fa_sol")<<A( 9, "fa_dun")<<A( 9, "fa_pair")<<"\n";
		for ( core::Size resnum_ = 1; resnum_ <= pose.total_residue(); ++resnum_) {
			if ( !in_pose.residue(resnum_).is_protein() ) continue;
			if( interface_obj.is_interface( resnum_ ) ) { // in interface

				Real total=in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( total_score ) ];
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

	else
	{
		core::Real const energy( compute( pose ) );
		out<<"Scoretype "<<name_from_score_type( score_type_ )<<" is "<<energy<<". ";
		bool const pass( energy <= threshold_ );
		if( pass ) out<<"passing."<<'\n';
		else out<<"failing."<<'\n';
	}
}


core::Real
EnergyPerResidueFilter::report_sm( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;

	core::Real const energy( compute( pose ) );
	return( energy );
}

core::Real
EnergyPerResidueFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;

	core::pose::Pose in_pose = pose;
	if( ( (*scorefxn_)[ interchain_env ] > 0.0 )
		&& ( (*scorefxn_)[ interchain_vdw ] > 0.0 )
		&& ( (*scorefxn_)[fa_rep] == 0.0 )
		&& ( (*scorefxn_)[fa_atr] == 0.0 ) )
		{
			if( in_pose.is_fullatom() ) {
			core::util::switch_to_residue_type_set( in_pose, core::chemical::CENTROID );
		}
	}
	else {
		if( in_pose.is_centroid() ) {
			core::util::switch_to_residue_type_set( in_pose, core::chemical::FA_STANDARD );
		}
	}

	in_pose.update_residue_neighbors();
	(*scorefxn_)( in_pose );
	core::Real weighted_score;
	if( score_type_ == total_score ) weighted_score = in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( score_type_ )];
	else {
		core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
		core::Real const score( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( score_type_ ) ]);
		weighted_score = weight * score ;
	}
	return( weighted_score );
}

BuriedUnsatHbondFilter::BuriedUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num ) :
	Filter( "BuriedUnsatHbonds" ),
	upper_threshold_( upper_threshold ),
	jump_num_( jump_num )
{ }

bool
BuriedUnsatHbondFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ) );

	TR<<"# unsatisfied hbonds: "<<unsat_hbonds<<". ";
	if( unsat_hbonds <= upper_threshold_ ){
		TR<<"passing." <<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
BuriedUnsatHbondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

core::Real
BuriedUnsatHbondFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	return( unsat_hbonds );
}

core::Real
BuriedUnsatHbondFilter::compute( core::pose::Pose const & pose ) const {

	runtime_assert( jump_num_ <= pose.num_jump() );

	// score the pose to populate the 10A neighborgraph
	core::pose::Pose bound( pose );

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	( *scorefxn )( bound );
	bound.update_residue_neighbors();

	core::pose::Pose unbound( bound );
	if( jump_num_ ) {
		core::Real const unbound_dist = 1000.0;
		protocols::moves::RigidBodyTransMover trans_mover( unbound, jump_num_ );
		trans_mover.trans_axis( trans_mover.trans_axis() );
		trans_mover.step_size(unbound_dist);
		trans_mover.apply( unbound );
		unbound.update_residue_neighbors();
	}

	basic::MetricValue< core::Size > mv_bound, mv_unbound;

	using namespace protocols::toolbox::pose_metric_calculators;
	// Despite the name, it's counting H-bonders, not any old polars.
	BuriedUnsatisfiedPolarsCalculator calc_bound("default", "default"), calc_unbound("default", "default");
	calc_bound.get("all_bur_unsat_polars", mv_bound, bound);

	core::Real unsat_hbonds( 0.0 );
	if( jump_num_ ) {
		calc_unbound.get("all_bur_unsat_polars", mv_unbound, unbound);
		unsat_hbonds = mv_bound.value() - mv_unbound.value();
		TR << "unbound_unsat=" << mv_unbound.value() << "    " << "bound_unsat=" << mv_bound.value() << std::endl;
	}
	else unsat_hbonds = mv_bound.value();

	return( unsat_hbonds );
}

bool
TerminusDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	TR<<"near terminus: "<<dist<<". " ;
	bool const status = (dist <= distance_) ? (false) : (true);
	if( status ) TR << "passing." << std::endl;
	else TR << "failing." << std::endl;
	return status;
}

void
TerminusDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	out<<"near terminus: "<< dist<<'\n';
}

core::Real
TerminusDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	return( dist );
}

core::Real
TerminusDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose copy_pose = pose;
	runtime_assert( copy_pose.num_jump() >= jump_num_ );

	// scoring is necessary for Interface to work reliably
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" ) );
	(*scorefxn)(copy_pose);

	protocols::scoring::Interface iface( jump_num_ );
	iface.distance( 8 );
	iface.calculate( copy_pose );
  core::Real min_dist(1000);

	for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
		core::Real dist(1000);
		if( !iface.is_interface( i ) ) continue; // keep going if we're not at the interface

		core::Size const chain = copy_pose.residue( i ).chain();
		core::Size const N_dist = i - copy_pose.conformation().chain_begin( chain );
		core::Size const C_dist = copy_pose.conformation().chain_end( chain ) - i;
		dist = ( N_dist <= C_dist ) ? ( N_dist) : ( C_dist ) ;
		if( ( N_dist < distance_ ) || ( C_dist < distance_ ) ) {
			return dist;
		}
		min_dist = ( dist < min_dist ) ? (dist) : (min_dist);
	}
	return( min_dist );
}


//// ALL DESTRUCTORS

ResiduesInInterfaceFilter::~ResiduesInInterfaceFilter() {}
AlaScan::~AlaScan() {}
ScoreTypeFilter::~ScoreTypeFilter() {}
NeighborTypeFilter::~NeighborTypeFilter() {}
ResidueBurialFilter::~ResidueBurialFilter(){}
ResidueDistanceFilter::~ResidueDistanceFilter(){}
DdgFilter::~DdgFilter() {}
HbondsToResidueFilter::~HbondsToResidueFilter() {}
EnergyPerResidueFilter::~EnergyPerResidueFilter() {}
BuriedUnsatHbondFilter::~BuriedUnsatHbondFilter(){}
TerminusDistanceFilter::~TerminusDistanceFilter(){}

protocols::filters::FilterOP
AlaScanFilterCreator::create_filter() const { return new AlaScan; }

std::string
AlaScanFilterCreator::keyname() const { return "AlaScan"; }

protocols::filters::FilterOP
BuriedUnsatHbondFilterCreator::create_filter() const { return new BuriedUnsatHbondFilter; }

std::string
BuriedUnsatHbondFilterCreator::keyname() const { return "BuriedUnsatHbonds"; }

protocols::filters::FilterOP
DdgFilterCreator::create_filter() const { return new DdgFilter; }

std::string
DdgFilterCreator::keyname() const { return "Ddg"; }

protocols::filters::FilterOP
EnergyPerResidueFilterCreator::create_filter() const { return new EnergyPerResidueFilter; }

std::string
EnergyPerResidueFilterCreator::keyname() const { return "EnergyPerResidue"; }

protocols::filters::FilterOP
HbondsToResidueFilterCreator::create_filter() const { return new HbondsToResidueFilter; }

std::string
HbondsToResidueFilterCreator::keyname() const { return "HbondsToResidue"; }

protocols::filters::FilterOP
InterfaceSasaFilterCreator::create_filter() const { return new InterfaceSasaFilter; }

std::string
InterfaceSasaFilterCreator::keyname() const { return "Sasa"; }

protocols::filters::FilterOP
NeighborTypeFilterCreator::create_filter() const { return new NeighborTypeFilter; }

std::string
NeighborTypeFilterCreator::keyname() const { return "NeighborType"; }

protocols::filters::FilterOP
ResidueBurialFilterCreator::create_filter() const { return new ResidueBurialFilter; }

std::string
ResidueBurialFilterCreator::keyname() const { return "ResidueBurial"; }

protocols::filters::FilterOP
ResidueDistanceFilterCreator::create_filter() const { return new ResidueDistanceFilter; }

std::string
ResidueDistanceFilterCreator::keyname() const { return "ResidueDistance"; }

protocols::filters::FilterOP
ResiduesInInterfaceFilterCreator::create_filter() const { return new ResiduesInInterfaceFilter; }

std::string
ResiduesInInterfaceFilterCreator::keyname() const { return "ResInInterface"; }

protocols::filters::FilterOP
ScoreTypeFilterCreator::create_filter() const { return new ScoreTypeFilter; }

std::string
ScoreTypeFilterCreator::keyname() const { return "ScoreType"; }

protocols::filters::FilterOP
TerminusDistanceFilterCreator::create_filter() const { return new TerminusDistanceFilter; }

std::string
TerminusDistanceFilterCreator::keyname() const { return "TerminusDistance"; }


} // protein_interface_design
} // devel
