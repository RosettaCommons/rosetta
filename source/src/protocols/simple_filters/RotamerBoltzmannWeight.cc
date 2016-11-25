// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeightFilterCreator.hh>

#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/util/SelectResiduesByLayer.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <utility/vector1.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/graph/Graph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/Energies.hh>

// Neil headers 110621
#include <core/pose/symmetry/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <ObjexxFCL/format.hh>

// Jacob headers 120423/120817
#include <core/pose/symmetry/util.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/AlaScan.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.RotamerBoltzmannWeight" );

/// @brief default ctor
RotamerBoltzmannWeight::RotamerBoltzmannWeight() :
	parent( "RotamerBoltzmannWeight" ),
	task_factory_( /* NULL */ ),
	rb_jump_( 1 ),
	sym_dof_names_( "" ),
	unbound_( true ),
	scorefxn_( /* NULL */ ),
	temperature_( 0.8 ),
	ddG_threshold_( 1.5 ),
	repacking_radius_( 6.0 ),
	energy_reduction_factor_( 0.5 ),
	compute_entropy_reduction_( false ),
	compute_max_( false ),
	repack_( true ),
	type_( "" ),
	skip_ala_scan_( false ),
	skip_report_( false ),
	fast_calc_(false),
	no_modified_ddG_(false),
	write2pdb_( false )
{
	threshold_probabilities_per_residue_.assign( core::chemical::num_canonical_aas, 1 );
}

void
RotamerBoltzmannWeight::repack( bool const repack ){
	repack_ = repack;
}

bool
RotamerBoltzmannWeight::repack() const
{
	return repack_;
}

void
RotamerBoltzmannWeight::threshold_probability( core::chemical::AA const aa_type, core::Real const probability )
{
	threshold_probabilities_per_residue_[ aa_type ] = probability;
}

core::Real
RotamerBoltzmannWeight::threshold_probability( core::chemical::AA const aa_type ) const
{
	return( threshold_probabilities_per_residue_[ aa_type ] );
}

void
RotamerBoltzmannWeight::energy_reduction_factor( core::Real const factor )
{
	energy_reduction_factor_ = factor;
}

core::Real
RotamerBoltzmannWeight::energy_reduction_factor() const
{
	return( energy_reduction_factor_ );
}

core::Real
RotamerBoltzmannWeight::temperature() const{
	return temperature_;
}

void
RotamerBoltzmannWeight::temperature( core::Real const temp )
{
	temperature_ = temp;
}

core::pack::task::TaskFactoryOP
RotamerBoltzmannWeight::task_factory() const
{
	return task_factory_;
}

void
RotamerBoltzmannWeight::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

void
RotamerBoltzmannWeight::rb_jump( core::Size const jump )
{
	runtime_assert( jump );
	rb_jump_ = jump;
}

core::Size
RotamerBoltzmannWeight::rb_jump() const
{
	return rb_jump_;
}

void
RotamerBoltzmannWeight::sym_dof_names( std::string const & dof_names )
{
	sym_dof_names_ = dof_names;
}

std::string
RotamerBoltzmannWeight::sym_dof_names() const
{
	return sym_dof_names_;
}

void
RotamerBoltzmannWeight::repacking_radius( core::Real const rad )
{
	repacking_radius_ = rad;
}

core::Real
RotamerBoltzmannWeight::repacking_radius() const
{
	return repacking_radius_;
}

bool
RotamerBoltzmannWeight::apply(core::pose::Pose const & ) const
{
	return( true );
}

/// @details iterate over all packable residues according to the task factory and for each
/// determine whether it contributes to binding at least ddG_threshold R.e.u. (using alanine scanning
/// This list will later be used to do the Boltzmann weight calculations.
utility::vector1< core::Size >
RotamerBoltzmannWeight::first_pass_ala_scan( core::pose::Pose const & pose ) const
{
	runtime_assert( task_factory() != nullptr );
	TR<<"----------First pass alanine scanning to identify hot spot residues------------"<<std::endl;
	utility::vector1< core::Size > hotspot_residues;
	hotspot_residues.clear();
	protocols::simple_filters::AlaScan ala_scan;
	ala_scan.repack( repack() );
	if ( repack() ) {
		ala_scan.repeats( 3 );
	} else {
		ala_scan.repeats( 1 );
	}
	core::Size const jump( rb_jump() );
	ala_scan.jump( jump );
	ala_scan.scorefxn( scorefxn() );
	if ( repack() ) {
		TR<<"All energy calculations will be computed subject to repacking in bound and unbound states (ddG).";
	} else {
		TR<<"All energy calculations will be computed subject to no repacking in the bound and unboud states (dG)";
	}
	core::Real orig_ddG(0.0);
	if ( !skip_ala_scan() ) {
		protocols::simple_filters::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn(), jump, 1 /*repeats*/ );
		ddg_filter.repack( repack() );
		orig_ddG = ddg_filter.compute( pose );
		TR<<"\nOriginal complex ddG "<<orig_ddG<<std::endl;
	}
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	for ( core::Size resi=1; resi<=pose.size(); ++resi ) {
		if ( packer_task->being_packed( resi ) && pose.residue( resi ).is_protein() ) {
			if ( skip_ala_scan() ) {
				TR<<"Adding residue "<<resi<<" to hotspot list\n";
				hotspot_residues.push_back( resi );
			} else {
				core::Real const ala_scan_ddG( ala_scan.ddG_for_single_residue( pose, resi ) );
				core::Real const ddG( ala_scan_ddG - orig_ddG );
				TR<<"ddG for resi "<<pose.residue( resi ).name3()<<resi<<" is "<<ddG<<std::endl;
				ddGs_[ resi ] = ddG;
				if ( ddG >= ddG_threshold() ) {
					hotspot_residues.push_back( resi );
				}
			}
		}
	}
	TR<<"----------Done with first pass alanine scanning to identify hot spot residues------------"<<std::endl;
	return( hotspot_residues );
}


core::Real
RotamerBoltzmannWeight::compute( core::pose::Pose const & const_pose ) const{
	ddGs_.clear();
	rotamer_probabilities_.clear();

	utility::vector1<Size> hotspot_res = core::pose::get_resnum_list_ordered(target_residues_,const_pose);
	// for(Size i = 1; i <= hotspot_res.size(); ++i) {
	//  std::cerr << "RotamerBoltzmannWeight calc res: " << hotspot_res[i] << std::endl;
	// }
	core::pose::Pose unbound_pose( const_pose );
	if ( hotspot_res.size()==0 ) {
		if ( type_ == "monomer" ) {
			core::select::util::SelectResiduesByLayer srb( true, true, false );
			utility::vector1< core::chemical::AA > select_aa_types;
			select_aa_types.push_back( core::chemical::aa_tyr );
			select_aa_types.push_back( core::chemical::aa_phe );
			select_aa_types.push_back( core::chemical::aa_trp );
			select_aa_types.push_back( core::chemical::aa_leu );
			select_aa_types.push_back( core::chemical::aa_ile );
			srb.restrict_aatypes_for_selection( select_aa_types );
			hotspot_res = srb.compute( const_pose );
		} else {
			if ( unbound() ) {
				using namespace protocols::moves;

				if ( sym_dof_names() != "" ) { // JB 120817
					utility::vector1<std::string> sym_dof_name_list = utility::string_split( sym_dof_names() , ',' ); // JB 120817
					for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) { // JB 120817
						int sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( unbound_pose, sym_dof_name_list[i] ); // JB 120817
						rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( unbound_pose, sym_aware_jump_id ) ); // JB 120817
						translate->step_size( 1000.0 ); // JB 120817
						translate->apply( unbound_pose ); // JB 120817
					} // JB 120817
				} else { // JB 120817
					int sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num(unbound_pose, rb_jump() ); // JB 120420
					rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( unbound_pose, sym_aware_jump_id ) ); // JB 120420
					translate->step_size( 1000.0 );
					translate->apply( unbound_pose );
				}
			}
			hotspot_res = first_pass_ala_scan( const_pose );
		}
	}
	//unbound_pose.dump_pdb("unbound_pose.pdb");
	if ( hotspot_res.size() == 0 ) {
		TR<<"No hot-spot residues detected in first pass alanine scan. Doing nothing"<<std::endl;
		return( 0 );
	} else {
		TR<<hotspot_res.size()<<" hot-spot residues detected."<<std::endl;
	}

	protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator rotboltz_calc( this->scorefxn(), this->temperature(), this->repacking_radius() );
	rotboltz_calc.set_lazy( true );

	for ( core::Size const hs_res : hotspot_res ) {
		core::Real const boltz_weight( fast_calc_ ? rotboltz_calc.computeBoltzWeight( unbound_pose, hs_res ) : compute_Boltzmann_weight( unbound_pose, hs_res ) );
		TR<<const_pose.residue( hs_res ).name3()<<hs_res<<" "<<boltz_weight<<'\n';
		rotamer_probabilities_[ hs_res ] = boltz_weight;
		if ( write2pdb() ) { write_to_pdb( hs_res, const_pose.residue( hs_res ).name3(), boltz_weight ); }
	}
	TR.flush();
	return -compute_boltz_probability();
}

core::Real
RotamerBoltzmannWeight::compute_boltz_probability() const
{
	if ( compute_max_ ) {
		core::Real max = 0.0;
		for ( std::map<core::Size,core::Real>::const_iterator i = rotamer_probabilities_.begin(); i != rotamer_probabilities_.end(); ++i ) {
			if ( i->second > max ) {
				max = i->second;
			}
		}
		return max;
	} else {
		core::Real avg = 0.0;
		for ( std::map<core::Size,core::Real>::const_iterator i = rotamer_probabilities_.begin(); i != rotamer_probabilities_.end(); ++i ) {
			avg += i->second;
		}
		return avg / (core::Real)rotamer_probabilities_.size();
	}
}

core::Real
RotamerBoltzmannWeight::compute_Boltzmann_weight( core::pose::Pose const & const_pose, core::Size const resi ) const
{
	using namespace core::pack::rotamer_set;
	using namespace core::pack::task;
	using namespace core::conformation;
	using core::pack::task::operation::TaskOperationCOP;

	TR<<"-----Computing Boltzmann weight for residue "<<resi<<std::endl;
	/// build a rotamer set for resi while pruning clashes against a poly-alanine background
	core::pose::Pose pose( const_pose );

	Residue const & res = pose.residue( resi );
	RotamerSetFactory rsf;
	RotamerSetOP rotset = rsf.create_rotamer_set( res );
	rotset->set_resid( resi );
	TaskFactoryOP tf( new core::pack::task::TaskFactory );
	tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ));
	tf->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
	PackerTaskOP ptask( tf->create_task_and_apply_taskoperations( pose ) );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, ptask); // NK 110621
	}
	ResidueLevelTask & restask( ptask->nonconst_residue_task( resi ) );
	restask.restrict_to_repacking();
	utility::graph::GraphOP packer_graph( new utility::graph::Graph( pose.size() ) );
	ptask->set_bump_check( true );
	rotset->build_rotamers( pose, *scorefxn_, *ptask, packer_graph, false );
	// TR << "num rotamers for resi " << resi << " is: " << rotset->num_rotamers() << std::endl;

	/// des_around will mark all residues around resi for design, the rest for packing.
	protocols::toolbox::task_operations::DesignAroundOperationOP des_around( new protocols::toolbox::task_operations::DesignAroundOperation );
	des_around->design_shell( repacking_radius() );
	des_around->include_residue( resi );
	tf->push_back( des_around );
	PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // NK 110621
	}
	protocols::simple_filters::ScoreTypeFilter const stf( scorefxn_, core::scoring::total_score, 0 );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm ); // NK 110621
	}
	mm->set_bb( false );
	if ( unbound() ) { // the complex was split, don't minimize rb dof
		mm->set_jump( false );
	} else { // minimize rb if bound
		mm->set_jump( rb_jump(), true ); // Need to modify for multicomponent symmetries (JBB)
	}
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( task->being_designed( i ) ) {
			task->nonconst_residue_task( i ).restrict_to_repacking(); // mark all des around to repacking only
			mm->set_chi( i, true );
		} else {
			task->nonconst_residue_task( i ).prevent_repacking(); /// mark all non-desaround positions for no repacking
			mm->set_chi( i, false );
		}
	}
	task->nonconst_residue_task( resi ).prevent_repacking();
	core::pack::pack_rotamers( pose, *scorefxn_, task );
	protocols::moves::MoverOP min_mover;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		min_mover = protocols::moves::MoverOP( new protocols::simple_moves::symmetry::SymMinMover( mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) ); // NK 110621
	} else {
		min_mover = protocols::moves::MoverOP( new protocols::simple_moves::MinMover( mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) ); // NK 110621
	}
	min_mover->apply( pose );
	core::pose::Pose const const_min_pose( pose );
	core::Real const init_score( stf.compute( const_min_pose ) );
	TR<<"Total score for input pose: "<<init_score<<std::endl;
	utility::vector1< core::Real > scores;
	for ( auto const & rotamer : *rotset ) {
		pose = const_min_pose;
		pose.replace_residue( resi, *rotamer, false/*orient bb*/ );
		core::pack::pack_rotamers( pose, *scorefxn(), task );
		// std::cerr << "CHI " << rotamer->chi1() << " " << rotamer->chi2() << std::endl;
		min_mover->apply( pose );
		core::Real const score( stf.compute( pose ) );
		TR<<"This rotamer has score "<<score<<std::endl;
		scores.push_back( score );
	}
	core::Real boltz_sum ( 0.0 );
	for ( core::Real const score : scores ) {
		boltz_sum += exp(( init_score - score )/temperature());
	}

	return( 1/boltz_sum );
}

core::Real
RotamerBoltzmannWeight::ddG_threshold() const
{
	return ddG_threshold_;
}

void
RotamerBoltzmannWeight::ddG_threshold( core::Real const ddG )
{
	ddG_threshold_ = ddG;
}

core::Real
RotamerBoltzmannWeight::report_sm( core::pose::Pose const & pose ) const
{
	compute( pose );
	if ( no_modified_ddG_ ) {
		return -compute_boltz_probability();
	} else {
		return compute_modified_ddG( pose, TR );
	}
}

core::Real
RotamerBoltzmannWeight::interface_interaction_energy( core::pose::Pose const & pose, core::Size const res ) const
{
	using namespace utility::graph;
	using namespace core::scoring;

	core::pose::Pose nonconst_pose( pose );
	(*scorefxn())(nonconst_pose);
	EnergyMap const curr_weights = nonconst_pose.energies().weights();

	core::Size const res_chain( pose.residue( res ).chain() );
	core::Real total_residue_energy( 0.0 );
	for ( EdgeListConstIterator egraph_it = nonconst_pose.energies().energy_graph().get_node( res )->const_edge_list_begin(); egraph_it != nonconst_pose.energies().energy_graph().get_node( res )->const_edge_list_end(); ++egraph_it ) {
		core::Size const int_resi = (*egraph_it)->get_other_ind( res );
		core::Size const int_resi_chain( pose.residue( int_resi ).chain() );
		if ( int_resi_chain != res_chain ) { // only sum over interaction energies with the residue's non-host chain
			EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);
			core::Real const intE = Eedge->dot( curr_weights );
			total_residue_energy += intE;
		}//fi
	}//for each egraph_it
	return( total_residue_energy );
}


void
RotamerBoltzmannWeight::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	if ( skip_report_ ) return;

	if ( type_ == "monomer" || no_modified_ddG_ ) {
		out<<"RotamerBoltzmannWeightFilter returns "<<compute( pose )<<std::endl;
		out<<"RotamerBoltzmannWeightFilter final report\n";
		out<<"Residue"<<'\t'<<"ddG"<<'\t'<<"RotamerProbability"<<'\n';
		for ( auto const & rot : rotamer_probabilities_ ) {
			core::Size const res( rot.first );
			core::Real const prob( rot.second );
			core::Real const ddG( 0.0 );
			out<<pose.residue( res ).name3()<< pose.pdb_info()->number( res )<<'\t'<<ddG<<'\t'<<prob<<'\t';
			out<<'\n';
		}
	} else {
		compute_modified_ddG( pose, out );
	}
}

void
RotamerBoltzmannWeight::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	using namespace utility::tag;

	TR << "RotamerBoltzmannWeightFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	repacking_radius( tag->getOption< core::Real >( "radius", 6.0 ) );
	type_ = tag->getOption< std::string >( "type", "" );
	rb_jump( tag->getOption< core::Size >( "jump", 1 ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names" , "" ) );
	unbound( tag->getOption< bool >( "unbound", 1 ) );
	ddG_threshold( tag->getOption< core::Real >( "ddG_threshold", 1.5 ) );
	temperature( tag->getOption< core::Real >( "temperature", 0.8 ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	energy_reduction_factor( tag->getOption< core::Real >( "energy_reduction_factor", 0.5 ) );
	compute_entropy_reduction( tag->getOption< bool >( "compute_entropy_reduction", 0 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	skip_report_ = tag->getOption< bool >( "skip_report", skip_report_ );
	utility::vector0< TagCOP > const & branch( tag->getTags() );
	for ( TagCOP const tag : branch ) {
		using namespace core::chemical;

		std::string const residue_type( tag->getName() );
		AA const aa( aa_from_name( residue_type ) );
		core::Real const threshold_probability_input( tag->getOption< core::Real >( "threshold_probability" ) );

		threshold_probability( aa, threshold_probability_input );
	}

	target_residues_ = tag->getOption<std::string>("target_residues","");

	fast_calc_ = tag->getOption< bool >("fast_calc",0);
	no_modified_ddG_ = tag->getOption< bool >("no_modified_ddG",0);

	compute_max_ = tag->getOption< bool >( "compute_max", 0 );
	skip_ala_scan( tag->getOption< bool >( "skip_ala_scan", 0 ) );
	write2pdb( tag->getOption< bool >( "write2pdb", 0 ) );
	TR<<"with options repacking radius: "<<repacking_radius()<<" and jump "<<rb_jump()<<" unbound "<<unbound()<<" ddG threshold "<<ddG_threshold()<<" temperature "<<temperature()<<" energy reduction factr "<<energy_reduction_factor()<<" entropy_reduction "<<compute_entropy_reduction()<<" repack "<<repack()<<" skip_ala_scan "<<skip_ala_scan()<<std::endl;
}

/// @brief Output per-residue Boltzmann weights to the output pdb file if desired
void RotamerBoltzmannWeight::write_to_pdb( core::Size const residue, std::string const & residue_name, core::Real const boltzmann_weight ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string user_name = this->get_user_defined_name();
	std::string output_string = "RotamerBoltzmannWeight " + user_name + ": " + residue_name + ObjexxFCL::string_of(residue) + " = " + ObjexxFCL::string_of(boltzmann_weight);
	job->add_string(output_string);

}

protocols::filters::FilterOP
RotamerBoltzmannWeight::fresh_instance() const{
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight() );
}

RotamerBoltzmannWeight::~RotamerBoltzmannWeight()= default;

protocols::filters::FilterOP
RotamerBoltzmannWeight::clone() const{
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight( *this ) );
}

void
RotamerBoltzmannWeight::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
RotamerBoltzmannWeight::scorefxn() const{
	return scorefxn_;
}

bool
RotamerBoltzmannWeight::unbound() const{
	return unbound_;
}

void
RotamerBoltzmannWeight::unbound( bool const u ){
	unbound_ = u;
}

bool
RotamerBoltzmannWeight::compute_entropy_reduction() const{
	return( compute_entropy_reduction_ );
}

void
RotamerBoltzmannWeight::compute_entropy_reduction( bool const cer ){
	compute_entropy_reduction_ = cer ;
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP RotamerBoltzmannWeightFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RotamerBoltzmannWeight ); }

// XRW TEMP std::string
// XRW TEMP RotamerBoltzmannWeightFilterCreator::keyname() const { return "RotamerBoltzmannWeight"; }

bool
RotamerBoltzmannWeight::skip_ala_scan() const
{
	return skip_ala_scan_;
}

void
RotamerBoltzmannWeight::skip_ala_scan( bool const s )
{
	skip_ala_scan_ = s;
}

std::string
RotamerBoltzmannWeight::type() const
{
	return type_;
}

void
RotamerBoltzmannWeight::type(std::string const & s)
{
	type_ = s;
}

void
RotamerBoltzmannWeight::no_modified_ddG( bool const no_ddg )
{
	no_modified_ddG_ = no_ddg;
}

void
RotamerBoltzmannWeight::target_residues( std::string const & target_residues_str )
{
	target_residues_ = target_residues_str;
}

bool
RotamerBoltzmannWeight::write2pdb() const
{
	return write2pdb_;
}

void
RotamerBoltzmannWeight::write2pdb( bool const write ) {
	write2pdb_ = write;
}

/// Note that compute( pose ) needs to have been run first. This merely sums over the probabilities
core::Real
RotamerBoltzmannWeight::compute_modified_ddG( core::pose::Pose const & pose, std::ostream & out ) const
{
	if ( no_modified_ddG_ ) return (0.0);
	protocols::simple_filters::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn(), rb_jump(), repack() ? 3 : 1 /*repeats*/ ); //120423
	ddg_filter.repack( repack() );

	core::Real const ddG_in( ddg_filter.compute( pose ) ); //problem
	core::Real modified_ddG( ddG_in );
	out<<"Residue"<<'\t'<<"ddG"<<'\t'<<"RotamerProbability"<<'\t'<<"EnergyReduction\n";
	for ( auto const & rot : rotamer_probabilities_ ) {
		core::Size const res( rot.first );
		core::Real const prob( rot.second );
		core::Real const ddG( ddGs_[ res ] );
		out<<pose.residue( res ).name3()<< pose.pdb_info()->number( res )<<'\t'<<ddG<<'\t'<<prob<<'\t';
		core::Real const threshold_prob( threshold_probability( pose.residue( res ).aa() ) );

		//out<<"prob: "<<prob<<"and threshold: " <<threshold_prob<<"\n";
		if ( prob <= threshold_prob ) {
			//core::Real const int_energy( interface_interaction_energy( pose, res ) );
			//core::Real const energy_reduction( energy_reduction_factor() * int_energy );
			core::Real const energy_reduction( -1 * temperature() * log( prob ) * energy_reduction_factor() );
			//out<<"after if statemnet energy reduction: "<<energy_reduction<<"this is prob: "<<prob<<"and this is threshold: " <<threshold_prob<<"\n";
			out<<energy_reduction;
			modified_ddG += energy_reduction;
		} else {
			out<<0;
		}
		out<<'\n';
	}
	out<<"ddG before, after modification: "<<ddG_in<<", "<<modified_ddG<<std::endl;
	return( modified_ddG );
}

std::string RotamerBoltzmannWeight::name() const {
	return class_name();
}

std::string RotamerBoltzmannWeight::class_name() {
	return "RotamerBoltzmannWeight";
}

std::string RotBoltzWeight_subelement_ct_name( std::string const & name ) {
	return "RotBoltzWeight_subelement_" + name + "Type";
}

void RotamerBoltzmannWeight::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "radius", xsct_real, "XRW TO DO", "6.0")
		+ XMLSchemaAttribute::attribute_w_default( "type", xs_string, "XRW TO DO", "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_non_negative_integer, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default( "sym_dof_names", xs_string, "XRW TO DO", "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default( "unbound", xsct_rosetta_bool, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default( "ddG_threshold", xsct_real, "XRW TO DO", "1.5")
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "XRW TO DO", "0.8");

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "energy_reduction_factor", xsct_real, "XRW TO DO", "0.5")
		+ XMLSchemaAttribute::attribute_w_default( "compute_entropy_reduction", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "repack", xsct_rosetta_bool, "XRW TO DO", "1")
		+ XMLSchemaAttribute( "skip_report", xsct_rosetta_bool, "XRW TO DO" );

	AttributeList RotBoltzWeight_subtag_attributes;
	RotBoltzWeight_subtag_attributes
		+ XMLSchemaAttribute( "restype", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute( "threshold_probability", xsct_real, "XRW TO DO");

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & RotBoltzWeight_subelement_ct_name );
	subelements.add_simple_subelement( "Threshold", RotBoltzWeight_subtag_attributes, "XRW TO DO");

	attlist + XMLSchemaAttribute::attribute_w_default( "target_residues", xs_string, "XRW TO DO", "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default( "fast_calc", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "no_modified_ddG", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "compute_max", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "skip_ala_scan", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "write2pdb", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default( "confidence", xsct_real, "Probability that the pose will be filtered out if the filter fails", "1.0" );
	XMLSchemaComplexTypeGenerator complex_type_generator;
	complex_type_generator
		.element_name( class_name() )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & protocols::filters::complex_type_name_for_filter )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

std::string RotamerBoltzmannWeightFilterCreator::keyname() const {
	return RotamerBoltzmannWeight::class_name();
}

protocols::filters::FilterOP
RotamerBoltzmannWeightFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight );
}

void RotamerBoltzmannWeightFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotamerBoltzmannWeight::provide_xml_schema( xsd );
}



} // simple_filters
} // protocols
