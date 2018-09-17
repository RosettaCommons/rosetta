// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/TryRotamers.cc
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/MapHotspot.hh>
#include <protocols/protein_interface_design/movers/MapHotspotCreator.hh>

// Package headers
#include <protocols/protein_interface_design/design_utils.hh>

// Project headers
#include <utility/graph/Graph.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/score_filters/ScoreTypeFilter.hh>
#include <protocols/protein_interface_design/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <algorithm>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.MapHotspot" );

// XRW TEMP std::string
// XRW TEMP MapHotspotCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return MapHotspot::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MapHotspotCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MapHotspot );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MapHotspot::mover_name()
// XRW TEMP {
// XRW TEMP  return "MapHotspot";
// XRW TEMP }

MapHotspot::MapHotspot() : protocols::moves::Mover( MapHotspot::mover_name() ) {}

MapHotspot::~MapHotspot() = default;

/// @details function for minimizing the interface
void
MapHotspot::MinimizeHotspots( core::pose::Pose & pose ){
	using namespace core::scoring;
	using namespace utility;

	core::Size const num_jump( pose.num_jump() );
	ScoreFunctionCOP scorefxn( minimization_scorefxns_[ num_jump ] );
	if ( scorefxn == nullptr ) {
		TR<<"skipping minimization b/c no scorefxn was defined"<<std::endl;
		return;
	}
	core::kinematics::FoldTree const hotspot_ft( make_hotspot_foldtree( pose ) );
	TR<<"imposing hotspot foldtree "<<hotspot_ft;
	pose.fold_tree( hotspot_ft );
	vector1< bool > const nomin( pose.size(), false );
	vector1< bool > minsc( pose.size(), false );
	minsc[ pose.size() ] = true;
	vector1< core::Size > const empty;
	vector1< bool > minrb_last( num_jump, false );
	minrb_last[ num_jump ] = true;
	MinimizeInterface( pose, scorefxn, nomin/*minbb*/, minsc, minrb_last, false/*optimize foldtree*/, empty/*target res*/ );
	if ( num_jump == 1 ) return;
	vector1< bool > minrb_all( num_jump, true );
	minrb_all[ num_jump ] = true;
	MinimizeInterface( pose, scorefxn, nomin/*minbb*/, minsc, minrb_all, false/*optimize foldtree*/, empty/*target res*/ );

}

/// @details utility function for dumping out a pose from GenerateMap
void
MapHotspot::output_pose( core::pose::Pose const & pose ) const{
	std::string residues("");
	for ( core::Size chain=2; chain<=pose.num_jump()+1; ++chain ) {
		core::Size const residue_num( pose.conformation().chain_begin( chain ) );
		char const residue_type( pose.residue( residue_num ).name1() );
		residues += residue_type;
	}
	//std::string name( file_name_prefix_ );
	std::ostringstream strm;
	strm<<file_name_prefix_<<"_"<<residues<<".pdb";
	pose.dump_pdb( strm.str() );
}

/// @details Utility function to deal with all of the ugliness of copying a single residue from one pose to another while conserving the jump
void
copy_hotspot_to_pose( core::pose::Pose const & src, core::pose::Pose & dest, core::Size const src_resi, core::chemical::ResidueType const restype, core::Size const jump )
{
	using namespace core::kinematics;
	using namespace core::conformation;
	using namespace core::chemical;
	Jump const saved_jump( src.jump( jump ) );
	FoldTree const saved_ft( src.fold_tree() );
	FoldTree new_ft;
	new_ft.clear();

	for ( Edge const & edge : saved_ft ) {
		if ( (core::Size) edge.start() <= src_resi && ( core::Size )edge.stop() <= src_resi ) {
			new_ft.add_edge( edge );
		}
	}//foreach edge

	ResidueOP new_res = ResidueFactory::create_residue( restype );

	dest.append_residue_by_jump( *new_res, dest.size(),"","",true/*new chain*/ );
	core::pose::add_variant_type_to_pose_residue( dest, core::chemical::SHOVE_BB, dest.size() );
	dest.fold_tree( new_ft );
	using namespace core::chemical;
	core::pose::add_upper_terminus_type_to_pose_residue( dest, src_resi );
	core::pose::add_lower_terminus_type_to_pose_residue( dest, src_resi );
	// core::pose::add_lower_terminus_type_to_pose_residue( dest, src_resi-1 );
	dest.conformation().update_polymeric_connection( src_resi );
	// dest.conformation().update_polymeric_connection( src_resi-1 );
	dest.set_jump( jump, saved_jump );
	dest.update_residue_neighbors();
}

/// @details utility function for creating a rotamer set an managing all of the rotamer options
core::pack::rotamer_set::RotamerSetOP
MapHotspot::create_rotamer_set( core::pose::Pose const & pose, core::Size const hotspot_resnum, core::Size const explosion ) const
{
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;
	using namespace core::pack::task::operation;
	using namespace core::scoring;
	TaskFactory tf;
	ScoreFunctionCOP scorefxn( get_score_function() );
	RotamerExplosionOP rotamer_exp_operation( new RotamerExplosion( hotspot_resnum, EX_THREE_THIRD_STEP_STDDEVS, explosion ) );
	InitializeFromCommandlineOP init_from_commandline( new InitializeFromCommandline );
	RestrictResidueToRepackingOP restrict_to_rep_operation( new RestrictResidueToRepacking );
	restrict_to_rep_operation->include_residue( hotspot_resnum );
	tf.push_back( rotamer_exp_operation );
	tf.push_back( restrict_to_rep_operation );
	tf.push_back( init_from_commandline );
	PackerTaskCOP ptask( tf.create_task_and_apply_taskoperations( pose ) );
	RotamerSetOP rotset = RotamerSetFactory::create_rotamer_set( pose );
	rotset->set_resid( hotspot_resnum );
	utility::graph::GraphOP packer_graph( new utility::graph::Graph( pose.size() ) );
	rotset->build_rotamers( pose, *scorefxn, *ptask, packer_graph, false );
	TR<<"Created rotamer set for residue "<<hotspot_resnum<<"with explosion="<<explosion<<std::endl;
	return( rotset );
}

/// @details A function that recursively goes over all jumps. Within each function, the residue-identities
/// related to that jump are cycled, and for each of those, it iterates over the rotamer set.
/// for each residue type, the best-energy (scorefxn) rotamer, that fulfills all of the user-defined filters
/// is chosen and then we move deeper into the recursion with the next jump. Stopping condition is the
/// final iteration over the final jump at which point a decoy is output, or if no rotamer for the
/// particular residue identity meets all filters (exit with no decoy generation).
void
MapHotspot::GenerateMap( core::pose::Pose const & start_pose, core::pose::Pose & curr_pose, core::Size const jump_number )
{
	core::Size const hotspot_resnum( start_pose.conformation().chain_begin( jump_number+1 ) );
	core::pose::Pose const saved_pose_1( curr_pose );
	TR<<"Allowed residues: "<< allowed_aas_per_jump_[ jump_number ]<<std::endl;
	for ( char const residue_type1 : allowed_aas_per_jump_[ jump_number ] ) {//iterate over residue types
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;
		using namespace core::chemical;
		using namespace core::scoring;
		ResidueType const & restype( *core::pose::get_restype_for_pose(
			start_pose,
			name_from_aa( aa_from_oneletter_code( residue_type1 ) ),
			start_pose.residue_type( hotspot_resnum ).mode()) );

		copy_hotspot_to_pose( start_pose, curr_pose, hotspot_resnum, restype, jump_number );
		core::pack::rotamer_set::RotamerSetCOP rotset( create_rotamer_set( curr_pose, hotspot_resnum, explosion_[ jump_number ]) );

		core::pose::Pose const saved_pose_2( curr_pose );
		core::Size rotset_size( 0 ); //ugly but I don't know how to find the size of cacheable data in rotset
		for ( auto rot_it = rotset->begin(); rot_it!=rotset->end(); ++rot_it, ++rotset_size ) {};
		TR<<"Iterating over "<<rotset_size<<" rotamers for residue "<<residue_type1<<" in jump #"<<jump_number<<std::endl;
		core::Size curr_rotamer_num( 1 );
		core::pose::Pose best_pose;
		core::Real lowest_energy( 100000 );
		ScoreFunctionCOP scorefxn( get_score_function() );
		score_filters::ScoreTypeFilter const pose_total_score( scorefxn, total_score, 100/*threshold*/ );
		for ( auto const & rot_it : *rotset ) {
			TR<<"Current rotamer: "<<curr_rotamer_num++<<'\n';
			core::conformation::ResidueCOP rot( rot_it );
			curr_pose.replace_residue( hotspot_resnum, *rot, false );
			jump_movers_[ jump_number ]->apply( curr_pose );
			MinimizeHotspots( curr_pose );
			bool pass_filter( true );
			for ( SizeFilter_map::const_iterator jump_filter_it( jump_filters_.begin() ); jump_filter_it!=jump_filters_.end(); ++jump_filter_it ) { // all filters must be satisfied
				if ( jump_filter_it->first <= jump_number ) {
					pass_filter = jump_filter_it->second->apply( curr_pose );
					if ( !pass_filter ) break;
				}//fi
			}//for jump filter
			if ( pass_filter ) {
				core::Real const total_score( pose_total_score.compute( curr_pose ) );
				if ( total_score <= lowest_energy ) {
					TR.Debug << "Current lowE=" << total_score << "   Prev lowE=" << lowest_energy << std::endl;
					best_pose = curr_pose;
					lowest_energy = total_score;
				}//fi total_score
			}//fi pass_filter
			curr_pose = saved_pose_2;
		}//foreach rotamer
		if ( lowest_energy >= 10000 ) {
			TR<<"No optimal pose found in jump "<<jump_number<<". Consider relaxing filters."<<std::endl;
		} else {
			curr_pose = best_pose;
			if ( jump_number == start_pose.num_jump() ) { // stopping condition
				output_pose( curr_pose );
			} else {
				GenerateMap( start_pose, curr_pose, jump_number + 1 );
			}
		}// esle
		curr_pose = saved_pose_1;
	}//foreach residue type
}

void
MapHotspot::apply( core::pose::Pose & pose ){
	core::pose::Pose chainA = *pose.split_by_chain( 1 );
	core::pose::Pose const start_pose( pose );
	pose = chainA; // this will make the viewer shows the working pose
	GenerateMap( start_pose, pose, 1 );
	set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
}

// XRW TEMP std::string
// XRW TEMP MapHotspot::get_name() const {
// XRW TEMP  return MapHotspot::mover_name();
// XRW TEMP }

void
MapHotspot::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & pose)
{
	using namespace utility::tag;

	clash_check_ = tag->getOption<bool>("clash_check", false );
	file_name_prefix_ = tag->getOption< std::string >( "file_name_prefix", "map_hs" );

	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	for ( TagCOP btag : branch_tags ) {
		std::string const btag_name( btag->getName() );
		if ( btag_name == "Jumps" ) {
			utility::vector0< TagCOP > const & jump_tags( btag->getTags() );
			runtime_assert( jump_tags.size() == pose.num_jump() );
			for ( TagCOP j_tag : jump_tags ) {
				auto const jump( j_tag->getOption< core::Size >( "jump" ) );
				bool const jump_fine( jump <= pose.num_jump() );
				if ( !jump_fine ) TR.Error<<"Jump "<<jump<<" is larger than the number of jumps in pose="<<pose.num_jump()<<std::endl;
				runtime_assert( jump_fine );
				explosion_[ jump ] = j_tag->getOption<core::Size>( "explosion", 0 );
				std::string const filter_name( j_tag->getOption<std::string>( "filter_name", "true_filter" ));
				auto find_filter( filters.find( filter_name ));
				bool const filter_found( find_filter != filters.end() );
				if ( !filter_found ) TR.Error<<"Filter "<<filter_name<<" not found in MapHotspot parsing"<<std::endl;
				runtime_assert( filter_found );
				jump_filters_[ jump ] = find_filter->second;
				std::string const mover_name( j_tag->getOption< std::string >( "mover_name", "null" ) );
				auto find_mover( movers.find( mover_name ) );
				bool const mover_found( find_mover != movers.end() );
				if ( !mover_found ) TR.Error<<"Mover "<<mover_name<<" not found in MapHotspot parsing"<<std::endl;
				runtime_assert( mover_found );
				jump_movers_[ jump ] = find_mover->second;
				std::string const allowed_aas( j_tag->getOption< std::string >( "allowed_aas", "ADEFIKLMNQRSTVWY" ) );
				allowed_aas_per_jump_[ jump ] = allowed_aas;
				std::string const scorefxn_name( rosetta_scripts::get_score_function_name(j_tag, "scorefxn_minimize") );
				if ( !data.has( "scorefxns", scorefxn_name ) ) {
					TR<<"Scorefxn "<<scorefxn_name<<" not found. Will not minimize sidechain.";
					minimization_scorefxns_[ jump ] = nullptr;
				} else {
					minimization_scorefxns_[ jump ] = rosetta_scripts::parse_score_function(j_tag, "scorefxn_minimize", data);
				}
			}//foreach j_tag
		} else { //fi btag_name=="Jumps"
			TR.Warning<<"Unrecognized branch tag from MapHotspot: "<<btag_name<<std::endl;
		}
	}//foreach branch_tag
}

std::string MapHotspot::get_name() const {
	return mover_name();
}

std::string MapHotspot::mover_name() {
	return "MapHotspot";
}

std::string jumpfunc( std::string const & manglee ) {
	return "jumpfuncmangled_" + manglee + "_type";
}

std::string jumpsfunc( std::string const & manglee ) {
	return "jumpsfuncmangled_" + manglee + "_type";
}

void MapHotspot::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "clash_check", xsct_rosetta_bool, "Perform a clash check", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "file_name_prefix", xs_string, "Prefix for filenames", "map_hs" );

	AttributeList subtag_attributes;
	subtag_attributes + XMLSchemaAttribute( "jump", xsct_non_negative_integer, "Jump in question" );
	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "explosion", xsct_non_negative_integer, "Degree of extensive sampling at this jump", "0" );
	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "filter_name", xs_string, "filter with which to eliminate bad solutions", "true_filter" );
	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "mover_name", xs_string, "mover for additional sampling", "null" );
	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "allowed_aas", xs_string, "amino acids allowed in design", "ADEFIKLMNQRSTVWY" );
	rosetta_scripts::attributes_for_parse_score_function( subtag_attributes, "scorefxn_minimize" );

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen_jump;
	ct_gen_jump.complex_type_naming_func( & jumpfunc )
		.element_name( "Jump" )
		.description( "Define a jump container" )
		.add_attributes( subtag_attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	utility::tag::XMLSchemaSimpleSubelementList jump_ssl;
	jump_ssl.add_already_defined_subelement( "Jump", & jumpfunc );

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen_jumps;
	ct_gen_jumps.complex_type_naming_func( & jumpsfunc )
		.element_name( "Jumps" )
		.description( "Map a hotspot conformation from a source to a target pose" )
		.set_subelements_repeatable( jump_ssl )
		.write_complex_type_to_schema( xsd );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Jumps", subtag_attributes, "XRW TODO"/*, 0 minoccurs*/ );

	std::string doc_string = "Map out the residues that might serve as a hotspot region on a target. This requires"
		" massive user guidance. Each hot-spot residue should be roughly placed by the user (at least as backbone)"
		" against the target. Each hot-spot residue should have a different chain ID. The method iterates over all"
		" allowed residue identities and all allowed rotamers for each residue. Tests its filters and for the subset"
		" that pass selects the lowest-energy residue by score12. Once the first hot-spot residue is identified it"
		" iterates over the next and so on until all hot-spot residues are placed. The output contains one file per"
		" residue identity combination.";

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), doc_string, attlist, ssl );
}

std::string MapHotspotCreator::keyname() const {
	return MapHotspot::mover_name();
}

protocols::moves::MoverOP
MapHotspotCreator::create_mover() const {
	return protocols::moves::MoverOP( new MapHotspot );
}

void MapHotspotCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MapHotspot::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

