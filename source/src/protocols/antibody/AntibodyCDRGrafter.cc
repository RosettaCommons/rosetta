// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyCDRGrafter.cc
/// @brief Class to graft CDR loops from an antibody to a new antibody or from a CDR pose into a different antibody.  
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/AntibodyCDRGrafter.hh>
#include <protocols/antibody/AntibodyCDRGrafterCreator.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>

#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.AntibodyCDRGrafter" );


namespace protocols {
namespace antibody {

	using namespace protocols::grafting;
	using namespace protocols::grafting::simple_movers;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::antibody::constraints;
	using namespace protocols::antibody::design;
	
AntibodyCDRGrafter::AntibodyCDRGrafter():
	protocols::moves::Mover( "AntibodyCDRGrafter" ),
	ab_info_(/* NULL */),
	donor_structure_(/* NULL */),
	scorefxn_(/* NULL */),
	graft_mover_(/* NULL */),
	anchored_graft_mover_(/* NULL */)
{
	setup_classes();
	set_defaults();
}

AntibodyCDRGrafter::AntibodyCDRGrafter( AntibodyInfoOP ab_info ):
	protocols::moves::Mover( "AntibodyCDRGrafter"),
	ab_info_(ab_info),
	donor_structure_(/* NULL */),
	scorefxn_(/* NULL */),
	graft_mover_(/* NULL */),
	anchored_graft_mover_(/* NULL */)
{
	setup_classes();
	set_defaults();

}


AntibodyCDRGrafter::AntibodyCDRGrafter(
		AntibodyInfoOP ab_info,
		core::pose::Pose const & donor_structure,
		utility::vector1<bool> const & cdrs_to_graft,
		core::Size nter_overhang /* 3 */,
		core::Size cter_overhang /* 3 */ ):
	
	protocols::moves::Mover( "AntibodyCDRGrafter"),
	ab_info_(ab_info),
	donor_structure_(/* NULL */),
	scorefxn_(/* NULL */),
	scorefxn_low_(/* NULL */),
	graft_mover_(/* NULL */),
	anchored_graft_mover_(/* NULL */)
{
	setup_classes();
	set_defaults();
	set_donor_structure(donor_structure);
	
	set_cdrs( cdrs_to_graft );

	nter_overhang_ = nter_overhang;
	cter_overhang_ = cter_overhang;
	
}


AntibodyCDRGrafter::~AntibodyCDRGrafter(){}
		
AntibodyCDRGrafter::AntibodyCDRGrafter( AntibodyCDRGrafter const & src ):
	protocols::moves::Mover( src ),
	ab_info_(src.ab_info_),
	donor_structure_(src.donor_structure_),
	scorefxn_(src.scorefxn_),
	scorefxn_low_(src.scorefxn_low_),
	cdrs_to_graft_(src.cdrs_to_graft_),
	graft_mover_(src.graft_mover_),
	anchored_graft_mover_(src.anchored_graft_mover_),
	use_secondary_graft_mover_(src.use_secondary_graft_mover_),
	nter_overhang_(src.nter_overhang_),
	cter_overhang_(src.cter_overhang_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_),
	optimize_cdrs_(src.optimize_cdrs_),
	include_cdr4_(src.include_cdr4_),
	neighbor_cdrs_(src.neighbor_cdrs_),
	dihedral_cst_weight_(src.dihedral_cst_weight_)
	
{
	
}



void
AntibodyCDRGrafter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption("cdrs") ) {
		TR << "Setting CDRs from settings" << std::endl;
		cdrs_to_graft_ = get_cdr_bool_from_tag(tag, "cdrs", true /* include cdr4 */);
	}

	if (tag->hasOption("cdr") ) {
		cdrs_to_graft_ = get_cdr_bool_from_tag(tag, "cdr", true /* include cdr4 */);
	}
	

	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}
	
	//Donor from File
	if (tag->hasOption("donor_structure_from_pdb") && tag->hasOption("donor_structure_from_spm")){
		utility_exit_with_message("Cannot pass both donor_structure_from_pdb and donor_structure_from_spm");
	}
	else if ( tag->hasOption( "donor_structure_from_pdb" ) ){
		std::string pdb_file = tag->getOption< std::string >( "donor_structure_from_pdb" );
		donor_structure_ = core::import_pose::pose_from_pdb( pdb_file );
	}
	else if ( tag->hasOption( "donor_structure_from_spm" ) ){
		donor_structure_ = protocols::rosetta_scripts::saved_reference_pose(tag, data, "donor_structure_from_spm");
	}
	else {
		utility_exit_with_message("donor_structure_from_pdb or donor_structure_from_spm RS option required for this mover.  Cannot continue.");
	}
	
	//use_secondary
	use_secondary_graft_mover_ = tag->getOption< bool >("use_secondary_graft_mover", use_secondary_graft_mover_);
	
	//Scorefunction
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	
	//Nter and Cter overhang
	nter_overhang_ = tag->getOption< core::Size >("nter_overhang", nter_overhang_);
	cter_overhang_ = tag->getOption< core::Size >("cter_overhang", cter_overhang_);
	
	//Stop after closure
	stop_after_closure_ = tag->getOption< bool >("stop_after_closure", stop_after_closure_);
	
	optimize_cdrs_ = tag->getOption< bool >("optimize_cdrs", optimize_cdrs_);
	include_cdr4_ = tag->getOption< bool >("optimize_cdr4_if_neighbor", include_cdr4_);
	
	dihedral_cst_weight_ = tag->getOption< core::Real >("dihedral_cst_wt", dihedral_cst_weight_);
}

protocols::moves::MoverOP
AntibodyCDRGrafter::clone() const{

	return protocols::moves::MoverOP( new AntibodyCDRGrafter( *this ) );
	
}

//AntibodyCDRGrafter & AntibodyCDRGrafter::operator=( AntibodyCDRGrafter const & src){
//	return AntibodyCDRGrafter( src );
//}

moves::MoverOP
AntibodyCDRGrafter::fresh_instance() const
{

	return protocols::moves::MoverOP( new AntibodyCDRGrafter );
	
}

std::string
AntibodyCDRGrafter::get_name() const {

	return "AntibodyCDRGrafter";
	
}

void
AntibodyCDRGrafter::setup_classes(){
	
	graft_mover_ = protocols::grafting::CCDEndsGraftMoverOP( new CCDEndsGraftMover() );
	anchored_graft_mover_ = protocols::grafting::AnchoredGraftMoverOP( new AnchoredGraftMover() );
	
}


void
AntibodyCDRGrafter::set_defaults(){

	//CCDEndsGraftMover
	graft_mover_->set_cycles(100);
	graft_mover_->set_scaffold_flexibility(2, 2);
	graft_mover_->set_insert_flexibility(2, 2);
	graft_mover_->final_repack(true);
	graft_mover_->stop_at_closure(true);
	graft_mover_->idealize_insert(true);
	graft_mover_->copy_pdbinfo(true);


	//AnchoredGraftMover
	anchored_graft_mover_->set_cycles(100);
	anchored_graft_mover_->set_scaffold_flexibility(2, 2);
	anchored_graft_mover_->set_insert_flexibility(2, 2);
	anchored_graft_mover_->final_repack(true);
	anchored_graft_mover_->stop_at_closure(true);
	anchored_graft_mover_->idealize_insert(true);
	anchored_graft_mover_->copy_pdbinfo(true);
	
	use_secondary_graft_mover_ = false;
	stop_after_closure_ = true;
	optimize_cdrs_ = false;
	include_cdr4_ = true;
	
	cdrs_to_graft_.clear();
	cdrs_to_graft_.resize(8, true);
	
	nter_overhang_ = 3;
	cter_overhang_ = 3;
	
	dihedral_cst_weight_ = 2.0;
	
	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);
	
	//Neighbor CDRs.
	//L1 - L2 L3 L4
	//L2 L1 L4
	//L3 L1 H3
	//L4 L1 L2
	
	//H1 - H2 H3 H4
	//H2 H1 H4
	//H3 H1 L3
	//H4 H1 H2
	
	neighbor_cdrs_.resize( 8, utility::vector1< bool >( 8, false ) );
	neighbor_cdrs_[ l1 ][ l2 ] = true;
	neighbor_cdrs_[ l1 ][ l3 ] = true;
	neighbor_cdrs_[ l1 ][ l4 ] = true;
	
	neighbor_cdrs_[ l2 ][ l1 ] = true;
	neighbor_cdrs_[ l2 ][ l4 ] = true;
	
	neighbor_cdrs_[ l3 ][ l1 ] = true;
	neighbor_cdrs_[ l3 ][ h3 ] = true;
	
	neighbor_cdrs_[ l4 ][ l1 ] = true;
	neighbor_cdrs_[ l4 ][ l2 ] = true;
	
	neighbor_cdrs_[ h1 ][ h2 ] = true;
	neighbor_cdrs_[ h1 ][ h3 ] = true;
	neighbor_cdrs_[ h1 ][ h4 ] = true;
	
	neighbor_cdrs_[ h2 ][ h1 ] = true;
	neighbor_cdrs_[ h2 ][ h4 ] = true;
	
	neighbor_cdrs_[ h3 ][ h1 ] = true;
	neighbor_cdrs_[ h3 ][ l3 ] = true;
	
	neighbor_cdrs_[ h4 ][ h1 ] = true;
	neighbor_cdrs_[ h4 ][ h2 ] = true;
	

	
}

void
AntibodyCDRGrafter::set_cdr_only(CDRNameEnum cdr){
	TR << "Setting " << ab_info_->get_CDR_name( cdr ) << " to graft"<<std::endl;
	cdrs_to_graft_.clear();
	cdrs_to_graft_.resize(8, false);
	cdrs_to_graft_[cdr] = true;
	
}

void
AntibodyCDRGrafter::set_cdrs(const utility::vector1<bool> & cdrs){
	cdrs_to_graft_ = cdrs;
	if ( cdrs.size() < CDRNameEnum_proto_total ){
		for (core::Size i = cdrs.size() +1; i <= CDRNameEnum_proto_total; ++i){
			cdrs_to_graft_.push_back( false );
		}
	}
}

void
AntibodyCDRGrafter::set_donor_structure(const core::pose::Pose & pose){
	donor_structure_ = core::pose::PoseOP( new core::pose::Pose(pose));
}

void
AntibodyCDRGrafter::set_idealize_insert(bool idealize_insert){
	graft_mover_->idealize_insert(idealize_insert);
	anchored_graft_mover_->idealize_insert(idealize_insert);
	
}

void
AntibodyCDRGrafter::set_overhang(core::Size nter_overhang, core::Size cter_overhang) {
	nter_overhang_ = nter_overhang;
	cter_overhang_ = cter_overhang;
	
}

void
AntibodyCDRGrafter::set_stop_after_closure( bool stop_after_closure ){
	graft_mover_->stop_at_closure( stop_after_closure );
	anchored_graft_mover_->stop_at_closure( stop_after_closure );
}

void
AntibodyCDRGrafter::set_use_secondary_graft_mover_if_needed( bool use_secondary_graft_mover ){
	use_secondary_graft_mover_ = use_secondary_graft_mover;
}

void
AntibodyCDRGrafter::set_scorefxn_pack( core::scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn->clone();
}

void
AntibodyCDRGrafter::set_scorefxn_low( core::scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_low_ = scorefxn->clone();
}

protocols::grafting::CCDEndsGraftMover &
AntibodyCDRGrafter::get_primary_graft_mover(){
	return *graft_mover_;
}

protocols::grafting::AnchoredGraftMover &
AntibodyCDRGrafter::get_secondary_graft_mover(){
	return *anchored_graft_mover_;
}

void
AntibodyCDRGrafter::set_optimize_cdrs( bool optimize_cdrs ){
	optimize_cdrs_ = optimize_cdrs;
}

void
AntibodyCDRGrafter::set_include_cdr4_in_optimization( bool include_cdr4 ) {
	include_cdr4_ = include_cdr4;
}

void
AntibodyCDRGrafter::set_dihedral_constraint_weight(core::Real dih_cst_wt){
	dihedral_cst_weight_ = dih_cst_wt;
}

void
AntibodyCDRGrafter::apply( core::pose::Pose& pose ){
	if (! donor_structure_){
		utility_exit_with_message("Donor structure required for AntibodyCDRGrafter!");
	}
	
	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoOP(new AntibodyInfo( pose, numbering_scheme_, cdr_definition_ ));
	}
	
	
	//Setup modeler for secondary graft
	CDRDihedralConstraintMover dih_cst = CDRDihedralConstraintMover( ab_info_ );
	dih_cst.set_use_cluster_csts( false ); //No use of cluster constraints until the clusters are generalized.
	
	//Pass any set scorefunctions
	if (! scorefxn_ ){
		scorefxn_ = core::scoring::get_score_function();
	}
	
	graft_mover_->set_fa_scorefunction( scorefxn_ );
	anchored_graft_mover_->set_fa_scorefunction( scorefxn_ );
	
	if (scorefxn_low_){
		graft_mover_->set_cen_scorefunction( scorefxn_low_ );
		graft_mover_->set_cen_scorefunction( scorefxn_low_ );
	}
	
	
	KeepRegionMover cookie_cutter = KeepRegionMover();
	
	for (core::Size i = 1; i <=  CDRNameEnum_proto_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (! cdrs_to_graft_[ i ]) continue;
		
		TR <<"Grafting CDR: "<< ab_info_->get_CDR_name( cdr ) << std::endl;
		
		core::Size start = ab_info_->get_CDR_start( cdr, *donor_structure_ );
		core::Size end = ab_info_->get_CDR_end( cdr, *donor_structure_ );
		
		cookie_cutter.start( start - nter_overhang_);
		cookie_cutter.end( end + nter_overhang_ );

		core::pose::Pose cdr_region = *donor_structure_; //Make copy so its in memory on the stack.
		cookie_cutter.apply(cdr_region);
		
		Pose temp_pose = pose;
		
		std::pair<bool, core::Size> cb = apply_to_cdr(temp_pose, cdr_region, cdr, graft_mover_);
		if ( cb.first && use_secondary_graft_mover_ ) {
			TR << "Graft not fully closed. Using secondary graft mover" << std::endl;
			temp_pose = pose;
			cb = apply_to_cdr( temp_pose, cdr_region, cdr, anchored_graft_mover_ );

		}
		else if ( cb.first && ! use_secondary_graft_mover_) {
			TR << "Graft not fully closed. Use of secondary graft mover false.  Nothing to be done." << std::endl;
		}
		else {
			TR << "Success.  Graft closed." << std::endl;
		}

		pose = temp_pose;
		ab_info_->setup_CDR_cluster( pose, cdr );
		
		if (cdr != l4 && cdr != h4){
			check_fix_aho_cdr_numbering( ab_info_, cdr, pose); //If not AHO won't do anything.
		}
		
	}
	
	//Fix graft - as Anchored Graft Mover will close the CDR, but usually not leave it in a great shape.
	if ( optimize_cdrs_ ){
		
		GeneralAntibodyModeler modeler = GeneralAntibodyModeler( ab_info_ );
		core::scoring::ScoreFunctionOP scorefxn_cst = scorefxn_->clone();
		scorefxn_cst->set_weight_if_zero(core::scoring::dihedral_constraint, dihedral_cst_weight_);
		
		modeler.set_scorefunction( scorefxn_cst );
		modeler.set_scorefunction_min( scorefxn_cst ); //I really need to refactor the modeler class.  Refactor Trial #2 I guess.
		TR << "optimizing grafted and neighbor cdrs" << std::endl;
		
		//Setup CDRs to minimize based on graft and neighbors.
		utility::vector1< bool > cdrs_to_minimize( 8, false );
		for ( core::Size i = 1; i <= cdrs_to_graft_.size(); ++i ){
			CDRNameEnum cdr_grafted = static_cast< CDRNameEnum >( i );
			
			for ( core::Size x = 1; x <= neighbor_cdrs_[ cdr_grafted ].size(); ++x ){
				CDRNameEnum cdr_neighbor = static_cast< CDRNameEnum >( x );
				
				if (( cdr_neighbor == l4 || cdr_neighbor == h4 ) ){
					if (include_cdr4_) cdrs_to_minimize[ cdr_neighbor ] = true;
				}
				else {
					cdrs_to_minimize[ cdr_neighbor ] = true;
				}
				
			}//for neighors
			
			cdrs_to_minimize[ cdr_grafted ] = true;
			
		}//for cdrs_to_graft
		
		//Set constraints and cdr options.
		for (core::Size i = 1; i <= cdrs_to_minimize.size(); ++i){
			if (! cdrs_to_minimize[ i ]) continue;
			
			CDRNameEnum cdr = static_cast< CDRNameEnum >( i );
			dih_cst.set_cdr( cdr );
			dih_cst.apply( pose );
			modeler.cdr_overhang(cdr, 2);
		}
		
		modeler.set_cdrs( cdrs_to_minimize );
		
		//modeler.minimize_cdrs( temp_pose, true /* min_sc */, true /* min_neighbor_cdrs */, false );
		//modeler.repack_cdrs( temp_pose, true /* min_neighbor_cdrs */ );
		
		scorefxn_->score(pose); //Segfault prevention.
		
		TR << "relaxing CDRs with constraints" <<std::endl;
		modeler.relax_cdrs( pose,  true /* min_sc */ );
	}
	
}

std::pair<bool, core::Size>
AntibodyCDRGrafter::apply_to_cdr(core::pose::Pose& pose, core::pose::Pose & cdr_region, CDRNameEnum const cdr, AnchoredGraftMoverOP grafter){

	core::Size start = ab_info_->get_CDR_start( cdr, pose );
	core::Size end = ab_info_->get_CDR_end( cdr, pose );
	
	grafter->set_insert_region( start - 1, end + 1 );
	
	core::Size nter_flex = grafter->get_nterm_insert_flexibility();
	core::Size cter_flex = grafter->get_cterm_insert_flexibility();
	grafter->set_piece(cdr_region, nter_overhang_, cter_overhang_);
	if ( cdr_region.total_residue() - nter_overhang_ - cter_overhang_ <= 4 ) {
		grafter->set_insert_flexibility( 1, 1 );
	}
	grafter->apply(pose);
	grafter->set_insert_flexibility( nter_flex, cter_flex );
	protocols::loops::remove_cutpoint_variants( pose, true );

	TR << "Checking Geometry" << std::endl;
	std::pair<bool, core::Size> cb = design::check_cb( pose, grafter->get_loops() );
	return cb;

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
AntibodyCDRGrafterCreator::create_mover() const {
	return protocols::moves::MoverOP( new AntibodyCDRGrafter );
}

std::string
AntibodyCDRGrafterCreator::keyname() const {
	return AntibodyCDRGrafterCreator::mover_name();
}

std::string
AntibodyCDRGrafterCreator::mover_name(){
	return "AntibodyCDRGrafter";
}

}//protocols
}//antibody


