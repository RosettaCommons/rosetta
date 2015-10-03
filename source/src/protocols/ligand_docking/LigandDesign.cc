// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandDesign.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/LigandDesign.hh>
#include <protocols/ligand_docking/LigandDesignCreator.hh>

#include <core/pose/util.hh>


#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <utility/assert.hh>

// option key includes

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ChemicalManager.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer ligand_design_tracer( "protocols.ligand_docking.LigandDesign", basic::t_debug );

std::string
LigandDesignCreator::keyname() const
{
	return LigandDesignCreator::mover_name();
}

protocols::moves::MoverOP
LigandDesignCreator::create_mover() const {
	return protocols::moves::MoverOP( new LigandDesign );
}

std::string
LigandDesignCreator::mover_name()
{
	return "LigandDesign";
}

LigandDesign::LigandDesign(): Mover("LigandDesign"){
	Mover::type( "LigandDesign" );
	set_fragments();
}

LigandDesign::LigandDesign(LigandDesign const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that ),
	option_file_(that.option_file_),
	fragments_(that.fragments_)
{}

LigandDesign::~LigandDesign() {}

void
LigandDesign::set_fragments(){
core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
core::chemical::ResidueTypeSetCOP rsd_set= cm->residue_type_set( core::chemical::FA_STANDARD );
// following may not work -- however, i could not find test cases with FRAGMENT residue_types.
core::chemical::ResidueTypeCOPs fragment_types = core::chemical::ResidueTypeFinder( *rsd_set ).base_property( core::chemical::FRAGMENT ).get_all_possible_residue_types();
	ligand_design_tracer<< fragment_types.size()<< " fragment_types"<< std::endl;

	BOOST_FOREACH ( core::chemical::ResidueTypeCOP fragment_type, fragment_types ) {
		core::conformation::ResidueOP temp( new core::conformation::Residue( *fragment_type, true) );
		fragments_.push_back(temp);
		ligand_design_tracer<< "frag_name: "<< temp->name()<< std::endl;
	}
}

protocols::moves::MoverOP LigandDesign::clone() const {
	return protocols::moves::MoverOP( new LigandDesign( *this ) );
}

protocols::moves::MoverOP LigandDesign::fresh_instance() const {
	return protocols::moves::MoverOP( new LigandDesign );
}

std::string LigandDesign::get_name() const{
	return "LigandDesign";
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
LigandDesign::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "LigandDesign" ) {
		ligand_design_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("option_file") ) {
		option_file_ = tag->getOption<std::string>("option_file");
	}
	//parse_score_function( tag, datamap, filters, movers, pose );
	//parse_task_operations( tag, datamap, filters, movers, pose );
}

bool grow(core::pose::Pose pose, core::Size start, core::Size end){
	return passes_filters(pose,start,end) && has_incomplete_connections(pose,start,end);
}

bool has_incomplete_connections(core::pose::Pose pose, core::Size start, core::Size const end){
	for ( ; start <= end; ++start ) {
		core::conformation::Residue const & res= pose.residue(start);
		ligand_design_tracer<< res.name();
		if ( res.has_incomplete_connection() ) {
			ligand_design_tracer<<" has incomplete connection"<< std::endl;
			return true;
		}
		ligand_design_tracer<<" completely connected"<< std::endl;
	}
	ligand_design_tracer<< "no more connections"<< std::endl;
	return false;
}

bool passes_filters(core::pose::Pose const & pose, core::Size start, core::Size const end){
	/// User Defined Filters
	// example
	if ( core::pose::num_heavy_atoms(start,end,pose)>20 ) {
		ligand_design_tracer<< "Reached heavy atom limit"<< std::endl;
		return false;
	}
	if ( core::pose::num_chi_angles(start,end,pose)>10 ) {
		ligand_design_tracer<< "Reached chi angle limit"<< std::endl;
		return false;
	} else return true;
}

core::Size
random_connection(core::conformation::ResidueCOP residue){
	//find connections that are incomplete
	utility::vector1<core::Size> incomplete_connections= get_incomplete_connections(residue);
	return numeric::random::rg().random_element(incomplete_connections);
}

// should be a helper function of residue class
utility::vector1<core::Size>
get_incomplete_connections(core::conformation::ResidueCOP residue){
	utility::vector1<core::Size> incomplete_connections;
	for ( core::Size i=1; i<= residue->n_residue_connections(); ++i ) {
		if ( residue->connection_incomplete(i) ) {
			incomplete_connections.push_back(i);
		}
	}
	return incomplete_connections;
}

utility::vector1<core::Size>
find_unconnected_residues(
	core::pose::Pose const & pose,
	core::Size start,
	core::Size const end
){
	utility::vector1<core::Size> unconnected_residues;
	for ( ; start <= end; ++start ) {
		if ( pose.residue(start).has_incomplete_connection() ) {
			unconnected_residues.push_back(start);
		}
	}
	return unconnected_residues;
}

void
LigandDesign::apply( core::pose::Pose & pose )
{
	assert(!fragments_.empty());
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	core::Size ligand_residue_id= pose.n_residue();
	ASSERT_ONLY(core::conformation::Residue const & ligand= pose.residue(ligand_residue_id);)
		assert(ligand.is_ligand());
	assert( ligand.n_residue_connections() > 0);
	core::Size const & chain_id= pose.chain(ligand_residue_id);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size end = pose.conformation().chain_end(chain_id);
	while ( grow(pose, start, end) ) {
		utility::vector1<core::Size> unconnected_residues = find_unconnected_residues(pose, start,end);
		core::Size const & grow_from= numeric::random::rg().random_element(unconnected_residues);
		core::conformation::ResidueCOP const growth= numeric::random::rg().random_element(fragments_);;
		core::Size grow_from_connection= random_connection(pose.residue(grow_from).get_self_ptr());
		core::Size growth_connection= random_connection(growth);
		pose.append_residue_by_bond(*growth, true, growth_connection, grow_from, grow_from_connection);
		//end = pose.conformation().chain_end(chain_id);
		++end;
	}

}

void LigandDesign::fragments_to_string() const{
	BOOST_FOREACH ( core::conformation::ResidueCOP fragment, fragments_ ) {
		core::conformation::Residue const & res= *fragment;
		std::string name= res.name();
		core::Size total= res.n_residue_connections();
		core::Size incomplete= get_incomplete_connections(res.get_self_ptr()).size();
		ligand_design_tracer<< "name, total, incomplete"<< name<< " "<< total<<" "<< incomplete<< std::endl;
	}
}

void LigandDesign::add_scores_to_job(
	core::pose::Pose & /*pose*/
)
{
	// using namespace core::scoring;
	// ScoreFunctionCOP sfxn = command_map_.get_score_function();
	//
	// core::Real const tot_score = sfxn->score( pose );
	//
	// // Which score terms to use
	// typedef utility::vector1<ScoreType> ScoreTypeVec;
	// ScoreTypeVec score_types;
	// for(int i = 1; i <= n_score_types; ++i) {
	//  ScoreType ii = ScoreType(i);
	//  if ( sfxn->has_nonzero_weight(ii) ) score_types.push_back(ii);
	// }
	//
	// protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	// for(ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii) {
	//  job->add_string_real_pair(name_from_score_type(*ii), sfxn->get_weight(*ii) * pose.energies().total_energies()[ *ii ]);
	// }
	// job->add_string_real_pair(name_from_score_type(core::scoring::total_score), tot_score);
}

} // namespace ligand_docking
} // namespace protocols
