// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_insertion/InsertionSiteTestMover.cc
/// @brief  cc file for InsertionSiteTestMover
/// @author Florian Richter, flosopher@gmail.com, february 2013


// Unit headers
#include <devel/domain_insertion/InsertionSiteTestMover.hh>
#include <devel/domain_insertion/InsertionSiteTestMoverCreator.hh>

// package headers

// Project headers
#include <basic/Tracer.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/MetricValue.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <devel/enzdes/EnzdesRemodelProtocol.hh>

#include <basic/datacache/DataMap.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>

#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>


namespace devel {
namespace domain_insertion {

static thread_local basic::Tracer tr( "devel.domain_insertion.InsertionSiteTestMover" );

InsertionSiteTestMover::InsertionSiteTestMover()
:
	sfxn_(/* NULL */), flex_window_(2),
	test_insert_ss_("LLLLLLLL"), insert_allowed_score_increase_(10.0),
	insert_attempt_sasa_cutoff_(30.0),
	length_of_insert_(test_insert_ss_.size() ), num_repeats_(5), pdb_numbering_(true),
	enz_flexbb_prot_( protocols::enzdes::EnzdesFlexBBProtocolOP( new protocols::enzdes::EnzdesFlexBBProtocol() ) ),
	insert_seqmap_(/* NULL */)
	//mostly arbitrary numbers, but can be modified through RS tag
{
	insert_test_pos_.clear();
}

InsertionSiteTestMover::InsertionSiteTestMover( InsertionSiteTestMover const & other )
: parent( other ),
	sfxn_(other.sfxn_), insert_test_pos_(other.insert_test_pos_),
	flex_window_(other.flex_window_), test_insert_ss_(other.test_insert_ss_),
	insert_allowed_score_increase_(other.insert_allowed_score_increase_),
	insert_attempt_sasa_cutoff_(other.insert_attempt_sasa_cutoff_), length_of_insert_(other.length_of_insert_),
	num_repeats_(other.num_repeats_), pdb_numbering_(other.pdb_numbering_),
	enz_flexbb_prot_( protocols::enzdes::EnzdesFlexBBProtocolOP( new protocols::enzdes::EnzdesFlexBBProtocol() ) ), insert_seqmap_(other.insert_seqmap_)
{}


InsertionSiteTestMover::~InsertionSiteTestMover(){}

protocols::moves::MoverOP
InsertionSiteTestMover::clone() const{
	return protocols::moves::MoverOP( new InsertionSiteTestMover( *this ) );
}


std::string
InsertionSiteTestMover::get_name() const
{
	return "InsertionSiteTestMover";
}


/// @details
void
InsertionSiteTestMover::apply( core::pose::Pose & pose )
{

	core::pose::Pose input_pose = pose;
	core::Real input_score = (*sfxn_)(input_pose);

	tr << "instest apply called, insert_test_pos_ has size " << insert_test_pos_.size() << std::endl;

	//spit out title line of results, this should be solved better, but ok for now
	tr <<"result Insert_pos   dSco_rawins_start    dSco_relax_start    dSco_anchor    rms_neigh    rms_anchor   diff_nlc_anchor  sasa_insert   deg_bury_insert     sasa_orig    orig_dssp" << std::endl;

	for( Size i = 1; i <= insert_test_pos_.size(); ++i){

		Size insert_pos( insert_test_pos_[ i ] );
		tr << "Starting insertion test at position " << insert_pos << "..." << std::endl;
		if( pdb_numbering_ ){
			tr << "PDB res " << insert_pos << " is pose res ";
			insert_pos = input_pose.pdb_info()->pdb2pose('A', insert_pos );
			tr << insert_pos << std::endl;
		}
		core::pose::Pose lowE_rlx_pose, lowE_rawins_pose;
		Real lowE(100000.0);
		bool one_succesful_attempt(false);

		for( Size repeat = 1; repeat <= num_repeats_; ++repeat){
			pose = input_pose;

			if( create_raw_insert_pose( pose, insert_pos ) ){
				one_succesful_attempt = true;
			//pose.dump_pdb("test_insertion_pos"+utility::to_string(insert_pos)+".pdb" );

				core::pose::Pose rawinsert_pose = pose;
				(*sfxn_)(rawinsert_pose);
				relax_raw_insert_pose( pose, rawinsert_pose, insert_pos );

				Real rlx_score = (*sfxn_)(pose);
				tr << "Test insertion at pos " << insert_pos << ", repeat " << repeat << ", produced a score diff of " << rlx_score - input_score << std::endl;
				if( (repeat == 1) || (rlx_score < lowE) ){
					lowE_rlx_pose = pose;
					lowE_rawins_pose = rawinsert_pose;
					lowE = rlx_score;
				} //new low score
			} //if remodeling worked
		} //loop over repeats

		if( one_succesful_attempt ){
			evaluate_insert_pose( input_pose, lowE_rawins_pose, lowE_rlx_pose, insert_pos );
			core::Real diff_score = lowE - input_score;
			tr << "Test insertion at pos " << insert_pos << " produced a best score diff of " << diff_score << std::endl;
			if( pdb_numbering_ ) lowE_rlx_pose.dump_pdb("test_insertrlx_pos"+utility::to_string( input_pose.pdb_info()->number(insert_pos) )+".pdb" );
			else lowE_rlx_pose.dump_pdb("test_insertrlx_pos"+utility::to_string(insert_pos)+".pdb" );


		}

		else{
			if( pdb_numbering_ ) tr<<"result " << ObjexxFCL::format::I( 8, input_pose.pdb_info()->number(insert_pos) ) << "unsuccesful." << std::endl;
			else tr<<"result " << ObjexxFCL::format::I( 8, insert_pos ) << "unsuccesful." << std::endl;
		}

	} // loop over insert pos

} //apply function

void
InsertionSiteTestMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
)
{
	insert_test_pos_.clear();

	sfxn_ = protocols::rosetta_scripts::parse_score_function( tag, "sfxn", data )->clone();
	enz_flexbb_prot_->set_scorefxn( sfxn_ );

	if( tag->hasOption("test_insert_ss") ){
		test_insert_ss_ = tag->getOption<std::string>( "test_insert_ss");
	}

	if( tag->hasOption("allowed_sc_increase") ){
		insert_allowed_score_increase_ = tag->getOption<Real>( "allowed_sc_increase");
	}

	if( tag->hasOption("seqpos")){
		insert_test_pos_.push_back( tag->getOption<Size>("seqpos") );
	}

	if( tag->hasOption("repeats")){
		num_repeats_ =  tag->getOption<Size>("repeats");
	}
	if( tag->hasOption("pdb_numbering")){
		pdb_numbering_ = tag->getOption<bool>("pdb_numbering");
	}

	utility::vector0< utility::tag::TagCOP > const subtags( tag->getTags() );
	for( utility::vector0< utility::tag::TagCOP >::const_iterator it= subtags.begin(); it!=subtags.end(); ++it ) {

		utility::tag::TagCOP const subtag = *it;
		if( subtag->getName() == "span" ) {
			core::Size const begin( core::pose::get_resnum( subtag, pose, "begin_" ) );
			core::Size const end( core::pose::get_resnum( subtag, pose, "end_" ) );
			runtime_assert( end > begin );
			runtime_assert( begin>=1);
			//runtime_assert( end<=reference_pose_->total_residue() );
			for( core::Size i=begin; i<=end; ++i ) insert_test_pos_.push_back( i );
		}
	}

	tr << "InsertionSiteTestMover parse_my_tag: trying insert string " << test_insert_ss_ << " at positions: ";
	for( Size i = 1; i <= insert_test_pos_.size(); ++i ) tr << insert_test_pos_[i] << ", ";
	tr << std::endl;
}


core::pack::task::PackerTaskCOP
InsertionSiteTestMover::make_insert_task(
	core::pose::Pose const & pose,
	Size insert_pos
) const
{
	core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_( pose ) );

	utility::vector1< bool > repack_pos( pose.total_residue(), false );
	repack_pos[ insert_pos ] = true;
	for( Size i = 1; i < flex_window_; ++i ){
		repack_pos[ insert_pos - i ] = true;
		repack_pos[ insert_pos + i ] = true;
	}
	repack_pos[ insert_pos + flex_window_] = true;

	enz_flexbb_prot_->get_tenA_neighbor_residues( pose, repack_pos );

	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( repack_pos[ i ] ) task->nonconst_residue_task( i ).restrict_to_repacking();
		else  task->nonconst_residue_task( i ).prevent_repacking();
	}
	return task;
}

devel::enzdes::EnzdesRemodelMoverOP
InsertionSiteTestMover::make_enzremodel_mover(
	core::pose::Pose & pose,
	Size insert_pos
) const
{

	core::pack::task::PackerTaskCOP insert_task( make_insert_task( pose, insert_pos ));

	devel::enzdes::EnzdesRemodelMoverOP enzremodel_mover( new devel::enzdes::EnzdesRemodelMover( enz_flexbb_prot_, insert_task, enz_flexbb_prot_->enz_flexible_region( 1 ) ) );

	std::string ss_edges;
	for( core::Size i = 1; i <= flex_window_; ++i ) ss_edges = ss_edges + "L";
	utility::vector1< std::string > ss_strings;
	ss_strings.push_back( ss_edges + test_insert_ss_ + ss_edges );
	enzremodel_mover->set_user_provided_ss( ss_strings );
	enzremodel_mover->set_max_allowed_score_increase(  insert_allowed_score_increase_ );

	enzremodel_mover->set_keep_existing_aa_identities( true );

	return enzremodel_mover;
}

bool
InsertionSiteTestMover::create_raw_insert_pose(
	core::pose::Pose & pose,
	core::Size const insert_pos
)
{

	//sanity check: only attempt at positions that are above a given sasa threshold
	core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius ];
	core::id::AtomID_Map< core::Real > atom_sasa_dummy;
	utility::vector1< core::Real > full_pose_residue_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa_dummy, full_pose_residue_sasa, probe_radius);
	core::Real insertpos_start_sasa( full_pose_residue_sasa[ insert_pos ] + full_pose_residue_sasa[ insert_pos + 1] );
	if( insertpos_start_sasa < insert_attempt_sasa_cutoff_ ) return false;

	insert_seqmap_ = NULL;
	enz_flexbb_prot_->add_flexible_region( insert_pos - (flex_window_ - 1), insert_pos + flex_window_, pose, true );

	devel::enzdes::EnzdesRemodelMoverOP enzremodel_mover(make_enzremodel_mover( pose, insert_pos ) );
	enzremodel_mover->set_keep_existing_aa_identities( true );

	enzremodel_mover->apply( pose );

	//the sequnce mapping needs to be changed by hand because the immediate insert
	//residues are considered not present anymore after vlb rebuild
	//but here we sorta know what they are
	core::id::SequenceMappingOP seqmap( new core::id::SequenceMapping( *(enzremodel_mover->get_seq_mapping()) ) );
	if( !( (*seqmap)[insert_pos]) ) (*seqmap)[insert_pos] = insert_pos;
	if( !( (*seqmap)[insert_pos+1]) ){
		runtime_assert( ( (*seqmap)[insert_pos+2]) != 0 );
		(*seqmap)[insert_pos+1] = (*seqmap)[insert_pos+2] - 1;
	}

	insert_seqmap_ = seqmap;

	if( enzremodel_mover->get_last_move_status() == protocols::moves::MS_SUCCESS ) return true;
	return false;
}

void
InsertionSiteTestMover::relax_raw_insert_pose(
	core::pose::Pose & pose,
	core::pose::Pose const & raw_pose,
	core::Size const
)
{
	using core::pack::task::operation::TaskOperationCOP;
	
	Size repeats(2); //only very few repeats to iron out the worst transgressions
	protocols::relax::FastRelaxOP frelax( new protocols::relax::FastRelax( sfxn_, repeats ) );

	//for now we relax everything, if it takes too long we could do less
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi(false);
	mm->set_bb( false );
	mm->set_jump( false );

	utility::vector1<bool> flex_res(pose.total_residue(), false );
	enz_flexbb_prot_->get_tenA_neighbor_residues( pose, flex_res );
	//std::set< core::Size > flex_res = enz_flexbb_prot_->enz_flexible_region(1)->get_10A_neighbors( pose );
	for( Size i = enz_flexbb_prot_->enz_flexible_region(1)->start(); i != enz_flexbb_prot_->enz_flexible_region(1)->stop(); ++i){
		flex_res[i] = true;
	}

	utility::vector1< Size > prevent_repack;
	for( Size i= 1; i <= pose.total_residue(); ++i){

		if( flex_res[i] == true ){
			mm->set_chi( i, true );
			mm->set_bb(i, true);
		}
		else prevent_repack.push_back( i );
	}

	core::pack::task::TaskFactoryOP taskf( new core::pack::task::TaskFactory() );
	taskf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ) );
	taskf->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ) );
	taskf->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::PreventResiduesFromRepackingOperation( prevent_repack ) ) );

	frelax->set_task_factory( taskf );
	frelax->set_movemap( mm );

	frelax->apply( pose );

	//we'll superimpose the relaxed pose onto the unrelaxed pose,
	//but ignoring the loop residues
	//makes analysis later cleaner
	core::id::AtomID_Map< core::id::AtomID> rms_atom_map(  core::id::BOGUS_ATOM_ID );
	core::pose::initialize_atomid_map( rms_atom_map, pose, core::id::BOGUS_ATOM_ID );
	for( Size i = 1; i < enz_flexbb_prot_->enz_flexible_region(1)->start(); ++i){
		if( pose.residue_type(i).is_protein() )	rms_atom_map.set(  core::id::AtomID( pose.residue( i ).atom_index("CA"), i ), core::id::AtomID( raw_pose.residue( i ).atom_index("CA"), i ) );
	}

	for( Size i = enz_flexbb_prot_->enz_flexible_region(1)->stop()+1; i <= pose.total_residue() ; ++i){
		if( pose.residue_type(i).is_protein() )	rms_atom_map.set(  core::id::AtomID( pose.residue( i ).atom_index("CA"), i ), core::id::AtomID( raw_pose.residue( i ).atom_index("CA"), i ) );
	}
	core::scoring::superimpose_pose( pose, raw_pose, rms_atom_map );
}


/// @details
///we will check the following things
///score diff rawinsert - start
///score diff relax_pose - start
///score diff relax_pose - start for anchor residues only
///difference in contacts for anchor residues and neighbor residues
///rms neighbors of relax region rawinsert -> relax
///rms of anchor residues start -> relax
///average sasa of inserted residues
///buriedness of inserted residues
///sasa of orig insert pos to compare
///for reference: dssp of the original input pose
void
InsertionSiteTestMover::evaluate_insert_pose(
	core::pose::Pose const & start_pose,
	core::pose::Pose const & rawinsert_pose,
	core::pose::Pose const & relax_pose,
	core::Size const insert_pos
)
{
	utility::vector1< Size > orig_anchor_res, shifted_anchor_res;
	for( Size i = insert_pos; i > (insert_pos -  flex_window_ ); --i){
		orig_anchor_res.push_back( i );
		shifted_anchor_res.push_back(i);
	}
	for( Size i = insert_pos + 1 ; i <= insert_pos + flex_window_; ++i ){
		orig_anchor_res.push_back( i );
		shifted_anchor_res.push_back( (*insert_seqmap_)[i] );
	}
	//score stuff begin
	core::Real diffsco_raw_start = rawinsert_pose.energies().total_energy() - start_pose.energies().total_energy();	core::Real diffsco_relax_start = relax_pose.energies().total_energy() - start_pose.energies().total_energy();
	core::Real diffsco_anchor(0.0);

	int diff_contacts_anchor(0);
	protocols::toolbox::pose_metric_calculators::NonlocalContactsCalculator nlc_calc(15, -0.5);
	basic::MetricValue< utility::vector1< Size > > start_nlc, relax_nlc;
	nlc_calc.get( "residue_nlcontacts", start_nlc, start_pose );
	nlc_calc.notify_energy_change();
	nlc_calc.get( "residue_nlcontacts", relax_nlc, relax_pose );
	for( Size i = 1; i <= orig_anchor_res.size(); ++i){
		diffsco_anchor += ( relax_pose.energies().residue_total_energy( shifted_anchor_res[i]) - start_pose.energies().residue_total_energy( orig_anchor_res[i] ) );
		diff_contacts_anchor += ( (relax_nlc.value())[ shifted_anchor_res[i] ] - (start_nlc.value())[ orig_anchor_res[i] ] );
	}

	tr << "scorestuff diffsco_raw_start: " << diffsco_raw_start << ", diffsco_relax_start: " << diffsco_relax_start << ", diffsco_anchor: " << diffsco_anchor << std::endl;
	//score stuff end

	//rms stuff begin
	//first we need to know the neighbors
	utility::vector1<bool> flex_res(relax_pose.total_residue(), false );
	enz_flexbb_prot_->get_tenA_neighbor_residues( relax_pose, flex_res );
	//but rmsd code wants an FArraz
	ObjexxFCL::FArray1D_bool flex_res_farray( relax_pose.total_residue(), false );
	tr << "Debug: when evaluating insert, the following are neighbor residues: ";
	for( Size i = 1; i <= flex_res.size(); ++ i ){
		if( flex_res[i] && !enz_flexbb_prot_->is_flexible(i) ){
			flex_res_farray(i) = true;
			tr << i << " ";
		}
	}
	tr << std::endl;
	tr << "and the pose fold tree is " << relax_pose.fold_tree() << std::endl;
	//tr << "Insert seqmapping is: " << std::endl;
	//insert_seqmap_->show( tr );
	//tr << std::endl;
	core::Real rms_neighbors( core::scoring::rmsd_no_super_subset( rawinsert_pose, relax_pose, flex_res_farray, core::scoring::is_protein_backbone ) );

	for( Size i =1; i <= relax_pose.total_residue(); ++ i ) flex_res_farray(i) = false; //reset
	//for( Size i = insert_pos; i > (insert_pos -  flex_window_ ); --i) flex_res_farray(i) = true;
	//for( Size i = (*insert_seqmap_)[ insert_pos + 1 ]; i <= (*insert_seqmap_)[ insert_pos + flex_window_ ]; ++i ) flex_res_farray(i) = true;
	tr << "anchor residues are ";
	for( Size i = 1; i <= orig_anchor_res.size(); ++i ){
		flex_res_farray( orig_anchor_res[i]) = true;
		tr << orig_anchor_res[i] << "+";
	}
	tr << std::endl;

	core::Real rms_anchors( core::scoring::rmsd_no_super_subset( start_pose, relax_pose, flex_res_farray, *insert_seqmap_, core::scoring::is_protein_backbone ) );
	//rms stuff end

	//sasa stuff begin
	core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius ];
	core::id::AtomID_Map< core::Real > atom_sasa_dummy;
	utility::vector1< core::Real > full_pose_residue_sasa;
	core::scoring::calc_per_atom_sasa( relax_pose, atom_sasa_dummy, full_pose_residue_sasa, probe_radius);
	core::Real insert_sasa(0.0);
	for( Size i = enz_flexbb_prot_->enz_flexible_region(1)->start(); i <= enz_flexbb_prot_->enz_flexible_region(1)->stop() ; ++i){
		insert_sasa += full_pose_residue_sasa[i];
	}
	//to calculate the degree of buriedness of the inserted residues,
	//we create a subpose consisting of the flexible region alone
	//and measure the sasa of this. the degree of buriedness is
	//then the sasa of the partial pose divided by the sasa
	//of the corresponding residues in the full pose
	core::kinematics::FoldTree newft(enz_flexbb_prot_->enz_flexible_region(1)->positions().size());
	newft.add_edge(1, enz_flexbb_prot_->enz_flexible_region(1)->positions().size(), core::kinematics::Edge::PEPTIDE );
	core::pose::Pose subpose;
	core::pose::create_subpose( relax_pose, enz_flexbb_prot_->enz_flexible_region(1)->positions(), newft, subpose );
	//subpose.dump_pdb("fresh_subpose.pdb");
	utility::vector1< core::Real > subpose_residue_sasa;
	core::scoring::calc_per_atom_sasa( subpose, atom_sasa_dummy, subpose_residue_sasa, probe_radius);
	core::Real insert_subpose_sasa = 0.0;
	for( Size i = 1; i <= subpose.total_residue(); ++i) insert_subpose_sasa += subpose_residue_sasa[i];
	core::Real buried_degree( 1 - (insert_sasa / insert_subpose_sasa) );

	//finally we calculate the sasa of the original two residues
	core::scoring::calc_per_atom_sasa( start_pose, atom_sasa_dummy, full_pose_residue_sasa, probe_radius);
	core::Real insertpos_start_sasa( full_pose_residue_sasa[ insert_pos ] + full_pose_residue_sasa[ insert_pos + 1] );

	//tr << "insert sasa: " << insert_sasa << ", insert_subpose_sasa: " << insert_subpose_sasa << ", degree: " << buried_degree << ", rms anchors: " << rms_anchors << ", rms_neighbors: " << rms_neighbors << std::endl;
	//sasa stuff end
	core::scoring::dssp::Dssp input_dssp( start_pose );
	std::string input_dssp_str(input_dssp.get_dssp_secstruct() );
	//tr << "complete input dssp " << input_dssp_str  << std::endl;
	std::string insert_region_dssp( input_dssp_str.substr(insert_pos-4, 8 ) );
	//now output
	using namespace ObjexxFCL;
	int output_res( insert_pos );
	if( pdb_numbering_ ) output_res = start_pose.pdb_info()->number(insert_pos);

	tr<<"result " << format::I( 8, output_res ) << format::F( 16, 2, diffsco_raw_start ) << format::F( 16, 2, diffsco_relax_start ) << format::F( 16, 2, diffsco_anchor ) << format::F( 16, 2, rms_neighbors ) << format::F( 16, 2, rms_anchors ) << format::I( 16, diff_contacts_anchor ) << format::F( 16, 2, insert_sasa ) << format::F( 16, 2, buried_degree ) << format::F( 16, 2, insertpos_start_sasa ) << "    " << insert_region_dssp << std::endl;

}


std::string
InsertionSiteTestMoverCreator::keyname() const
{
	return InsertionSiteTestMoverCreator::mover_name();
}

protocols::moves::MoverOP
InsertionSiteTestMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new InsertionSiteTestMover );
}

std::string
InsertionSiteTestMoverCreator::mover_name()
{
	return "InsertionSiteTestMover";
}



}
}
