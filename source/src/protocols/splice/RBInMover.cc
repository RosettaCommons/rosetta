// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RBInMover.cc
/// @brief

// Unit headers
#include <protocols/splice/RBInMover.hh>
#include <protocols/splice/RBOutMover.hh>
#include <protocols/splice/RBInMoverCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.splice.RBInMover" );
#include <core/pose/extra_pose_info_util.hh>
#include <utility/tag/Tag.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <string>
//#include <sstream>
#include <core/pose/util.hh>
#include <algorithm>
#include <utility/io/izstream.hh>
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <protocols/simple_moves/CutChainMover.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace splice {

using namespace::protocols;

// XRW TEMP std::string
// XRW TEMP RBInMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RBInMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RBInMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RBInMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RBInMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "RBIn";
// XRW TEMP }

RBInMover::RBInMover(): moves::Mover("RBIn"),
	RB_dbase_(""),
	from_entry_( 1 ),
	to_entry_( 1 ),
	current_entry_( 1 ),
	randomize_( true ),
	checkpointing_file_( "" ),
	modify_foldtree_( true ),
	db_entry_("")
{
	jump_library_.clear();
}

RBInMover::~RBInMover() = default;

void
RBInMover::jump_library( utility::vector1< std::pair <std::string,core::kinematics::Jump> > j ){ jump_library_ = j; }

utility::vector1< std::pair <std::string,core::kinematics::Jump> >
RBInMover::jump_library() const{ return jump_library_; }

bool
RBInMover::checkpoint() const{
	if ( checkpointing_file_ == "" ) {
		return false;
	}
	//write to file
	std::ofstream data;
	data.open( checkpointing_file_.c_str(), std::ios::out );
	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open checkpointing file for writing: " + checkpointing_file() + "\n" );
	}
	data<<current_entry_<<'\n';
	for ( std::pair <std::string,core::kinematics::Jump> const & jump : jump_library() ) {
		data<<jump<<'\n';
	}
	data.flush();
	TR<<"Finished writing checkpointing information to "<<checkpointing_file()<<std::endl;
	return true;
}

bool
RBInMover::checkpoint_recovery(){
	/* if ( checkpointing_file_ == "" ) {
	return false;
	}

	utility::io::izstream data( checkpointing_file_ );
	if ( !data ) {
	TR<<"Checkpointing file not yet written. Not recovering."<<std::endl;
	return false;
	}
	TR<<"Checkpointing file found. Recovering."<<std::endl;
	runtime_assert( jump_library_.size() == 0 ); /// if it's not 0, someone has already initialized the library and that doesn't make too much sense...
	std::string line;
	getline( data, line );
	std::istringstream data_stream( line );
	data_stream >> current_entry_;
	TR<<"Recovered entry: "<<current_entry_<<std::endl;
	while ( getline( data, line ) ) {
	TR<<"pushing jump into library: "<<line<<std::endl;
	core::kinematics::Jump jump;
	std::istringstream jump_stream( line );
	jump_stream >> jump;
	jump_library_.push_back( jump );
	}
	TR<<"Read "<<jump_library_.size()<<" jumps into library"<<std::endl;
	*/ return false;

}

void
RBInMover::init(){
	bool const checkpoint_recover( checkpoint_recovery() );
	if ( checkpoint_recover ) {
		return;
	}
	bool const first_pass( jump_library_.size() == 0 );
	if ( !first_pass ) { // already initialized; don't repeat
		return;
	}

	TR<<"Initializing rigid-body jump library"<<std::endl;
	utility::io::izstream data( RB_dbase_ );
	if ( !data ) {
		TR << "cannot open rigid-body database " << RB_dbase_ << std::endl;
		utility_exit();
	}
	std::string line;
	core::Size count = 0;
	while ( getline( data, line ) ) {
		TR<<"pushing jump into library: "<<line<<std::endl;
		std::vector<std::string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		std::string source_name(strs[1]);
		core::kinematics::Jump jump;
		//TR<<"Jump: "<<strs[0]<<std::endl;
		//TR<<"name: "<<strs[1]<<std::endl;
		std::istringstream data_stream( strs[0] );
		data_stream >> jump;
		count++;
		if ( count <= to_entry() && count >= from_entry() ) {
			jump_library_.push_back(std::make_pair(source_name,jump));
		}
	}
	if ( randomize() ) {
		TR<<"Randomizing the order in the jump_library"<<std::endl;
		numeric::random::random_permutation(jump_library_.begin(), jump_library_.end(), numeric::random::rg());
	}
	current_entry_ = 1;/// this should already be initialized, but just in case
	checkpoint(); /// first pass should also checkpoint to keep the order of entries
	TR<<"Done making jump library with size "<<jump_library_.size()<<std::endl;
}

void
RBInMover::set_fold_tree( Pose & pose ) {
	core::Size vl_vh_cut=find_vl_vh_cut(pose);
	utility::vector1<std::array<int, 3>> fold_tree_nodes;
	utility::vector1<core::Size> cys_pos;
	find_disulfide_postions(pose, cys_pos);
	fold_tree_nodes = set_fold_tree_nodes(pose,cys_pos, vl_vh_cut);
	core::kinematics::FoldTree ft;
	ft.clear();

	for ( auto i: fold_tree_nodes ) {
		TR<<i[0]<<","<<i[1]<<","<<i[2]<<std::endl;
		ft.add_edge(i[0],i[1],i[2]);
	}

	ft.add_edge((int)cys_pos[4],(int)cys_pos[2],1);
	ft.reorder((int)cys_pos[3]);
	ft.delete_self_edges();
	ft.check_fold_tree();
	TR<<ft<<std::endl;
	pose.fold_tree( ft );
}


void
RBInMover::apply( Pose & pose ){
	//init();
	if ( modify_foldtree() ) {
		set_fold_tree( pose );
	}
	TR<<"Setting jump now"<<std::endl;
	select_entry_by_name();
	TR<<jump_library_[ current_entry_ ].second <<std::endl;
	pose.set_jump( 1, jump_library_[ current_entry_ ].second );
	core::pose::add_comment(pose, "jump source:", jump_library_[ current_entry_ ].first);
	// Update RBO comment
	//std::ostream RBMoveStream;
	//RBMoveStream<<jump_library_[ current_entry_ ];
	//std::string RBMoveString = RBMoveStream.str();
	//std::replace( RBMoveString.begin(), RBMoveString.end(), ' ', ',');
	//TR<<"Current RBO jump: "<<RBMoveString<<std::endl;
	//core::pose::add_comment(pose,"RBO ",RBMoveString);

	if ( current_entry_ == jump_library_.size() ) { // rewind to the beginning of the library
		current_entry_ = 1;
	} else {
		++current_entry_; // go up one entry
	}
	checkpoint();
}

// XRW TEMP std::string
// XRW TEMP RBInMover::get_name() const {
// XRW TEMP  return RBInMover::mover_name();
// XRW TEMP }

moves::MoverOP
RBInMover::clone() const
{
	return moves::MoverOP( new RBInMover( *this ) );
}

moves::MoverOP
RBInMover::fresh_instance() const
{
	return moves::MoverOP( new RBInMover );
}

void
RBInMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	RB_dbase( tag->getOption< std::string >( "rigid_body_dbase", "" ) );
	from_entry( tag->getOption< core::Size >( "from_entry", 1 ) );
	to_entry( tag->getOption< core::Size >( "to_entry", 999999 /*no limit on entry size*/ ) );
	randomize( tag->getOption< bool >( "randomize", true  ) );
	checkpointing_file( tag->getOption< std::string >( "checkpointing_file", "" ) );
	modify_foldtree( tag->getOption< bool >( "modify_foldtree", true ) );
	db_entry( tag->getOption< std::string >( "db_entry", "" ) );

	runtime_assert( from_entry() <= to_entry() && from_entry() >= 1 );
	init();
}

std::string RBInMover::get_name() const {
	return mover_name();
}

std::string RBInMover::mover_name() {
	return "RBIn";
}

void RBInMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// AMW TODO
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "rigid_body_dbase", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "from_entry", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "to_entry", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "randomize", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "db_entry", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "modify_foldtree", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RBInMoverCreator::keyname() const {
	return RBInMover::mover_name();
}

protocols::moves::MoverOP
RBInMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RBInMover );
}

void RBInMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RBInMover::provide_xml_schema( xsd );
}

void RBInMover::find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos) {
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(core::chemical::DISULFIDE) ) {
			cys_pos.push_back(i);
		}
	}
	runtime_assert(cys_pos.size() >3 );

}

utility::vector1<std::array<int, 3>> RBInMover::set_fold_tree_nodes(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos, core::Size vl_vh_cut){
	utility::vector1<std::array<int, 3>> fold_tree_nodes;
	fold_tree_nodes.clear();

	/*Build the following Fold tree for the Vl/Vh fragment
	___________________________
	|                         |
	|                         |
	<------*<--------*------>//<-----*<-------*------->

	*/
	fold_tree_nodes.push_back({{(int) cys_pos[1],1,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[2],(int) cys_pos[1],-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[2],(int) vl_vh_cut,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[3],(int) vl_vh_cut+1,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[4],(int) cys_pos[3],-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[4],(int) pose.conformation().chain_end(1),-1}});
	// fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[2],1}});
	return fold_tree_nodes;
}

core::Size RBInMover::find_vl_vh_cut(core::pose::Pose pose) const {
	protocols::simple_moves::CutChainMover ccm;
	return ccm.chain_cut(pose); //
}

void
RBInMover::select_entry_by_name(){
	if ( db_entry_=="" ) {
		return;
	}
	for ( core::Size i = 1; i <= jump_library_.size(); ++i ) {
		if ( jump_library_[i].first==db_entry() ) {
			current_entry_ = i;
			return;
		}
	}
	utility_exit_with_message( "Could not find entry " + db_entry() + " in database file\n" );

}

} // simple_moves
} // protocols

