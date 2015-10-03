// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/seeded_abinitio/GrowPeptides.cc
/// @brief different ways of growing peptide sequences
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#include <protocols/seeded_abinitio/GrowPeptides.hh>
#include <protocols/seeded_abinitio/GrowPeptidesCreator.hh>
#include <protocols/seeded_abinitio/SeedFoldTree.fwd.hh>
#include <protocols/seeded_abinitio/SeedFoldTree.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/id/SequenceMapping.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//other
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>


// C++ headers
#include <string>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include  <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/util.hh>

//parser
#include <utility/tag/Tag.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

//loops
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loops_main.hh>

//util
#include <utility/vector1.hh>
#include <set>
#include <boost/foreach.hpp>

using namespace core;
using namespace protocols::seeded_abinitio;
static THREAD_LOCAL basic::Tracer TR( "protocols.seeded_abinitio.GrowPeptides" );


namespace protocols {
namespace seeded_abinitio {

using namespace protocols::moves;
using namespace core;

std::string
GrowPeptidesCreator::keyname() const{
	return GrowPeptidesCreator::mover_name();
}

protocols::moves::MoverOP
GrowPeptidesCreator::create_mover() const {
	return protocols::moves::MoverOP( new GrowPeptides() );
}

std::string
GrowPeptidesCreator::mover_name(){
	return "GrowPeptides";
}


GrowPeptides::GrowPeptides()
{}

GrowPeptides::~GrowPeptides() {}


protocols::moves::MoverOP
GrowPeptides::clone() const {
	return( protocols::moves::MoverOP( new GrowPeptides( *this ) ) );
}

protocols::moves::MoverOP
GrowPeptides::fresh_instance() const {
	return protocols::moves::MoverOP( new GrowPeptides );
}

bool
GrowPeptides::ddg(){
	return ddg_;
}

void
GrowPeptides::append_residues_nterminally ( Size seq_register, Size res_pos, Size stop, std::string & nat_seq , pose::Pose & target_seeds ){
	TR<<" ---- growing N-terminal stretch from residues: "<<res_pos <<" to " <<stop <<"----------" << std::endl;
	core::chemical::ResidueTypeSet const & rsd_set( target_seeds.residue(1).residue_type_set() );// this could be changed as needed
	//std::cout<<"nseq size: " << nat_seq.size() << " template sequence: "<< nat_seq << std::endl;
	//std::cout<< "sequence register : "<< seq_register << std::endl;

	for ( Size k= res_pos; k < stop ; ++k ) {
		Size resi = stop - k + res_pos - 1 ; /* so that stop doesnt get incorporated anymore*/
		const char aa = nat_seq[ stop - k - 1 + seq_register - 1];//in case this is called within the sequence, -1 because strings start counting at 0
		TR.Debug << "RES AA N-terminal extension:  " << resi << aa <<std::endl;
		// Representative type should have no/minimal variants
		core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( aa ) );
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
		target_seeds.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, res_pos, true);
		target_seeds.set_omega( res_pos, 180.0 );
	}
	//target_seeds.dump_pdb("nextended.pdb");
}

void
GrowPeptides::append_residues_cterminally ( Size seq_register, Size res_pos, Size stop, std::string & nat_seq , pose::Pose & target_seeds ){
	TR<<" ----- growing C-terminal extension from residues: "<<res_pos <<" to " <<stop <<"--------" <<std::endl;
	core::chemical::ResidueTypeSet const & rsd_set( target_seeds.residue( res_pos - 1 ).residue_type_set() );
	//std::cout<<"cseq size: " << nat_seq.size() << " sequence: "<< nat_seq << std::endl;
	//std::cout<<"seq_register: " << seq_register << std::endl;

	for ( Size j = res_pos ; j <  stop  ; ++j  ) {
		const char aa = nat_seq[ j - res_pos + seq_register /*-1*/]; // -1 for string adjustment
		Size resi =  j ;
		TR.Debug << "RES AA C-terminal extension  " << resi << aa <<std::endl;
		// Representative type should have no/minimal variants
		core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( aa ) );
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
		target_seeds.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, resi /*- 1*/ , true );// stop
		target_seeds.set_omega( resi , 180.0 );
		//target_seeds.dump_pdb( "ctermextn.pdb" );
	}
}

/*
void
insert_segment( std::pair <Size, Size> insert_type, std::string seq ,pose::pose curr_pose, )
{
//(0, num) = insert c-terminally of given position
//(num, 0) = insert n-terminally of given position
//(0,0 ) = insert n-terminally
//(0, total resi) = insert c-terminally
///replace a segment between the two given positions
///simply insert either N or C terminal of a given residue
}
*/


/*
void
GrowPeptides::process_length_change(
core::pose::Pose & pose,
core::id::SequenceMappingCOP smap
){
enz_prot_->remap_resid( pose, *smap );
core::id::combine_sequence_mappings( *start_to_current_smap_, *smap );

for( utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP >::iterator rcg_it = rcgs_.begin();
rcg_it != rcgs_.end(); ++rcg_it ){
(*rcg_it)->set_seqmap( smap );
}
}
*/

void
GrowPeptides::grow_from_verteces(
	core::pose::Pose & curr_pose,
	std::string sequence,
	protocols::loops::Loops & seeds,
	std::set< core::Size > vertex_set
){

	using namespace core;
	using namespace kinematics;

	core::pose::Pose saved_pose;
	saved_pose = curr_pose;
	core::util::switch_to_residue_type_set( curr_pose , core::chemical::CENTROID );

	core::kinematics::FoldTree grow_foldtree = curr_pose.fold_tree() ;
	TR<<"foldtree before growing: " << grow_foldtree << std::endl;

	utility::vector1< core::Size > verteces;
	BOOST_FOREACH ( const Size vertex, vertex_set ) {
		verteces.push_back( vertex );
	}

	TR<<"start new protein: " << verteces[1] <<" end: " << verteces[verteces.size()] <<" size sequence: " << sequence.size()<< std::endl;
	TR.Debug <<"number verteces: "<< verteces.size()<< std::endl;

	if ( verteces[verteces.size()] - (verteces[1] - 1) != sequence.size() ) {
		utility_exit_with_message( "chunk pieces do not agree with the length of the submitted template pdb" );
	}

	for ( Size vertex_it = 4 ; vertex_it <= verteces.size(); vertex_it = vertex_it + 4 ) {

		if ( vertex_it < verteces.size() ) {
			grow_foldtree = curr_pose.fold_tree();
			//connect seeds
			TR.Debug <<"for temp jump --- from: " <<verteces[vertex_it - 3 ]<<", to: "<< verteces[vertex_it - 3] + (seeds[vertex_it/4].stop() - seeds[vertex_it/4].start()) + 1 <<", cutpoint: "<< verteces[vertex_it - 3] + (seeds[vertex_it/4].stop() - seeds[vertex_it/4].start()) <<std::endl;
			//ensuring that the new temporariy cutpoint is unique
			Size temp_cutpoint = verteces[vertex_it - 3] + (seeds[vertex_it/4].stop() - seeds[vertex_it/4].start());
			if ( curr_pose.fold_tree().is_cutpoint( temp_cutpoint ) ) {
				++temp_cutpoint;
			}

			grow_foldtree.new_jump( verteces[vertex_it - 3 ], verteces[vertex_it - 3] + (seeds[vertex_it/4].stop() - seeds[vertex_it/4].start()) + 1 /*start of next seed before additions*/ , temp_cutpoint );
			curr_pose.fold_tree( grow_foldtree );
			TR<<"foldtree before adding new seeds, current seed: "<<vertex_it/4 << " with foldtree " << grow_foldtree << std::endl;
		}

		//need to adjust the verteces for the N-terminally extensions since the vertex container has the numbering for the complete sequences
		//and not the trunctated starting pieces.
		TR<<"appending N-terminally from: " << verteces[vertex_it - 3 ] << " to " << verteces[vertex_it - 2] <<std::endl;
		TR<<"appending C-terminally from: " << verteces[vertex_it - 1 ] << " to " << verteces[vertex_it ] << std::endl;

		//need to adjust the extensions/numbering
		Size nseq_start =  seeds[vertex_it/4 ].start() - ( verteces[vertex_it - 2] - verteces[vertex_it - 3] );
		TR.Debug <<" grow nterm: start sequence " <<nseq_start <<" position in pose "<< verteces[vertex_it - 3] <<" stop " << verteces[vertex_it - 2] <<" sequence " << sequence << std::endl ;

		//get sequence start through the seed start minus that actual length that needs to be added
		append_residues_nterminally( seeds[vertex_it /4].start() - ( verteces[vertex_it - 2] - verteces[vertex_it - 3] ) , verteces[vertex_it - 3] , verteces[vertex_it - 2], sequence , curr_pose ) ;

		TR <<"growing foldtree: " << grow_foldtree << std::endl;
		TR.Debug <<" grow cterm: start sequence " << seeds<<" position in pose "<< verteces[vertex_it - 1 ] <<" stop " << verteces[vertex_it] <<" sequence " << sequence << std::endl ;
		append_residues_cterminally( seeds[vertex_it/4].stop(), verteces[vertex_it - 1 ] , verteces[vertex_it], sequence , curr_pose ) ;

		//using simple foldtree that connects the two seeds with each other to keep consequutives seeds constants in space
		//connect cutpoint with next seed

		grow_foldtree.clear();
		TR <<"pose size: " << curr_pose.total_residue() << std::endl;
		grow_foldtree = curr_pose.fold_tree();
		TR.Debug <<"done growing, temporary foldtree: " << grow_foldtree << std::endl;
	}

	grow_foldtree = curr_pose.fold_tree();
	TR << "growing completed, temporary foldtree: " << grow_foldtree << std::endl;
}

void GrowPeptides::apply (core::pose::Pose & pose ){

	///adding a pose observer for downstream adjustment of sequence positions
	setup_cached_observers( pose );

	///if there are loops and template, then activate grow from seeds
	if ( all_seeds_.size() > 0  &&  template_presence ) {

		utility::vector1< Size > cutpoints;

		if ( !fetch_foldtree ) {
			TR<<"taking foldtree from pose" <<std::endl;
		}
		cutpoints = pose.fold_tree().cutpoints();

		if ( fetch_foldtree ) {
			TR<<"generate a foldtree through SeedFoldTree, and get cutpoints" << std::endl;
			core::pose::PoseOP tmp_seed_target_poseOP( new core::pose::Pose( pose ) );
			SeedFoldTreeOP seed_ft_generator( new SeedFoldTree() );
			seed_ft_generator->ddg_based( ddg() );
			seed_ft_generator->scorefxn( scorefxn_ );
			seed_ft_generator->anchor_specified(anchor_specified_);
			if ( anchor_specified_ ) {
				seed_ft_generator->set_anchor_res( anchors_ );
			}
			core::scoring::dssp::Dssp dssp( *template_pdb_ );
			dssp.insert_ss_into_pose( *template_pdb_ );
			std::string secstr_template = template_pdb_->secstruct();
			TR.Debug  << "sec str for template: " << secstr_template << std::endl;
			seed_foldtree_ = seed_ft_generator->set_foldtree( /**template_pdb_ ,*/ tmp_seed_target_poseOP, secstr_template, all_seeds_, true );
			verteces_ = seed_ft_generator->get_folding_verteces();
			TR.Debug<<"verteces for folding: " <<std::endl;
			cutpoints = seed_foldtree_->cutpoints();

			//debugging stuff
			BOOST_FOREACH ( core::Size const r, verteces_ ) {
				TR.Debug<< r <<"\t";
			}
		}
		std::string seq;
		if ( seq_ != "" ) {
			seq = seq_;
		} else {
			seq = template_pdb_->sequence();
		}
		if ( seq == "" ) utility_exit_with_message("no sequence specified" );

		grow_from_verteces( pose, seq, all_seeds_ , verteces_ );

		//add the new seed foldtree to the pose
		pose.fold_tree( *seed_foldtree_ );
		TR<<"set new foldtree: "<< pose.fold_tree() <<std::endl;

	}//end grow from seeds


	/*
	// ------------------ simple pose extensions based on teh parsers input -----------------------------

	if( all_seeds_.size() == 0 ){
	if( extend_nterm > 0 ){
	std::string nseq;

	if ( nsequence_.size() > 0 ){
	nseq = nsequence_;
	if( nsequence_.size() != extend_nterm ){
	TR<<"WARNING: specified sequence is not long enough for the desired length of extension, adding extra alanine residue" <<std::endl;
	for (Size i = 0 ; i < extend_nterm - nsequence_.size(); ++i )
	nseq += "A";
	}
	}

	if( all_ala_N ){
	TR<<"overwriting N-term sequence with all ALA, if not desired turn of all_ala_N" <<std::endl;
	for(Size i = 0; i < extend_nterm ; ++i){
	nseq += "A";
	std::cout<<nseq<<"\n"<< std::endl;
	}
	}

	append_residues_nterminally ( 0 , 1 , 1 , nseq , copy_pose );
	}//end nterm extension


	if ( extend_cterm > 0 ){
	std::string cseq;
	std::cout<<"cseq "<< cseq <<std::endl;

	if ( csequence_.size() > 0 ){
	cseq = csequence_;
	if( csequence_.size() != extend_cterm ){
	TR<<"WARNING: specified sequence is not long enough for the desired length of extension, adding extra alanine residue" <<std::endl;
	for (Size i = 0 ; i < extend_cterm - csequence_.size(); ++i )
	std::cout<<cseq<<std::endl;
	cseq = cseq + "A";
	}
	}

	if( all_ala_C ){
	TR<<"overwriting C-term sequence with all ALA, if not desired turn of all_ala_C" <<std::endl;
	for(Size i = 0; i < extend_cterm ; ++i){
	std::cout<<cseq<<"\n"<<std::endl;
	cseq += "A";
	}
	}
	std::cout<<"pose: "<< copy_pose.total_residue()<< "\n" <<copy_pose.total_residue() + extend_cterm <<"seq: "<<  cseq << std::endl;

	append_residues_cterminally( 0, copy_pose.total_residue(), copy_pose.total_residue() + extend_cterm, cseq , copy_pose );

	}//end cterm extension

	pose = copy_pose;

	}//end without loops
	*/

	if ( !output_centroid ) {
		TR<<"switching back to full_atom mode" <<std::endl;
		pose.update_residue_neighbors();
		core::util::switch_to_residue_type_set( pose, chemical::FA_STANDARD );
		( *scorefxn_ )( pose );
	}//end output centroid

	TR.flush();
}


/// @details putting a LengthEventCollector into the pose
void
GrowPeptides::setup_cached_observers( core::pose::Pose & pose ){
	core::pose::datacache::LengthEventCollectorOP lencollect( new core::pose::datacache::LengthEventCollector() );
	pose.observer_cache().set( core::pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR, lencollect );
}

std::string
GrowPeptides::get_name() const {
	return GrowPeptidesCreator::mover_name();
}

void
GrowPeptides::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	TR<<"GrowPeptides mover has been initiated" <<std::endl;
	//default
	template_presence = false;
	//need scorefxn to score extended/changed pose
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	TR<<"scoring with following scorefunction: " << *scorefxn_ <<std::endl;

	fetch_foldtree = tag->getOption< bool >( "SeedFoldTree", 0 );

	ddg_ = tag->getOption< bool >( "ddg_based", 0 );

	//add_chainbreakterm_ = tag->getOption< bool >( "add_chainbreakterm", 1 );

	if ( tag->hasOption( "extend_nterm" ) ) {
		extend_nterm = tag->getOption< Size >( "extend_nterm" , 0 );
		TR<<"extending peptide n-terminally by "<< extend_nterm << std::endl;
	}
	if ( tag->hasOption( "extend_cterm" ) ) {
		extend_cterm = tag->getOption< Size >( "extend_cterm" , 0 );
		TR<<"extending peptide c-terminally by " << extend_cterm << std::endl;
	}

	all_ala_N = tag->getOption< bool >("all_ala_N" , 0 );
	if ( all_ala_N ) {
		TR<<"N-terminally added amino acids are all ALA" <<std::endl;
	}

	all_ala_C = tag->getOption< bool >("all_ala_C", 0 );
	if ( all_ala_C ) {
		TR<<"C-terminally added amino acids are all ALA" <<std::endl;
	}

	if ( tag->hasOption( "nseq" ) ) {
		csequence_ = ( tag->getOption< std::string >( "nseq" ) );
	}

	if ( tag->hasOption( "cseq" ) ) {
		nsequence_ = ( tag->getOption< std::string >( "cseq" ) );
	}

	output_centroid = tag->getOption< bool >( "output_centroid", 0 );

	if ( tag->hasOption( "template_pdb" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "template_pdb" ));
		template_pdb_ = core::pose::PoseOP( new core::pose::Pose ) ;
		core::import_pose::pose_from_pdb( *template_pdb_, template_pdb_fname );
		TR<<"read in a template pdb with " <<template_pdb_->total_residue() <<"residues"<<std::endl;
		template_presence = true;
	}

	if ( tag->hasOption("sequence" ) ) {
		seq_ = tag->getOption< std::string >("sequence" );
	}

	if ( !template_presence && seq_ == "" ) {
		utility_exit_with_message("neither template pdb nor sequence for growing is specified!!" );
	}

	//parsing branch tags
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	BOOST_FOREACH ( TagCOP const btag, branch_tags ) {

		//parse the pdb of interest, which is either the template or the input pdb depending on the users specificiation
		if ( template_presence ) {
			curr_pose_ = template_pdb_;
		} else {
			curr_pose_ = core::pose::PoseOP( new pose::Pose( pose ) );
		}

		anchor_specified_ = false;

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file
			//needs some assertions to avoid bogus input
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			core::Size const begin( core::pose::parse_resnum( beginS, *curr_pose_ ) );
			core::Size const end( core::pose::parse_resnum( endS, *curr_pose_ ) );
			all_seeds_.add_loop( begin , end , 0, 0, false );
			TR <<"parsing seeds: \n"<< begin <<" and " << end <<std::endl;

			if ( btag->hasOption( "anchor" ) ) {
				Size anchor_res = btag->getOption< core::Size >("anchor", 0 );
				TR<<"anchor residue: " << anchor_res << std::endl;
				anchors_.push_back( anchor_res );
				anchor_specified_ = true;
			}
		}//end seed tags

		//not hooked in yet...
		if ( btag->getName() == "Steal_seq_span" ) {
			if ( !template_presence ) {
				utility_exit_with_message("need to specify a template pdb to steal sequenc spans");
			}
			std::string const begin_str( btag->getOption<std::string>( "begin" ) );
			std::string const end_str( btag->getOption<std::string>( "end" ) );
			core::Size const begin( core::pose::parse_resnum( begin_str, *template_pdb_ ) );
			core::Size const end( core::pose::parse_resnum( end_str, *template_pdb_ ) );
			runtime_assert( end > begin );
			runtime_assert( begin>=1);
			runtime_assert( end<=template_pdb_->total_residue() );
			std::string seq_chunk = template_pdb_->sequence();
			sequence_chunks_.push_back( seq_chunk );
		}//end steal sequence

	}//end branch tags
}//end parse my tag
}//seeded abinitio
} //end protocols


