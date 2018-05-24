// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/SmartAssembly.cc
/// @brief a SEWING Assembly composed of SmartSegments
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/pose/subpose_manipulation_util.hh>
//headers for transform_segment
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/Remarks.hh>

#include <utility/utility.functions.hh>

static basic::Tracer TR( "protocols.sewing.data_storage.SmartAssembly" );


namespace protocols {
namespace sewing {
namespace data_storage {
SmartAssembly::SmartAssembly():
	utility::pointer::ReferenceCount()
{
	//pre-allocate all_basis_pairs_
	all_basis_pairs_ = utility::vector1< std::pair< core::Size, core::Size > >( 10000, std::make_pair( 0, 0 ) );
}

SmartAssembly::SmartAssembly(hashing::SegmentVectorCOP segment_vector, core::Size window_width):
	utility::pointer::ReferenceCount()
{
	modifiable_terminus_ = 'B';
	start_node_vital_segments_ = "all";
	output_partner_ = true;
	segment_vector_ = segment_vector;
	window_width_ = window_width; //Default for constructor is 4 to preserve tests/etc. Default for mover is 4.
	for ( SmartSegmentCOP current_segment : *segment_vector_ ) {
		if ( current_segment->get_length() >= window_width_ ) {
			if ( !current_segment->is_n_terminus_fixed() ) {
				n_terminal_segments_.push_back( current_segment->get_segment_id() );
			} else if ( !current_segment->is_c_terminus_fixed() ) {
				c_terminal_segments_.push_back( current_segment->get_segment_id() );
			}
		}
	}
	TR << "Total n terminal segments: " << n_terminal_segments_.size() << std::endl;
	TR << "Total c terminal segments: " << c_terminal_segments_.size() << std::endl;
	//pre-allocate all_basis_pairs_
	all_basis_pairs_ = utility::vector1< std::pair< core::Size, core::Size > >( 10000, std::make_pair( 0, 0 ) );
}

SmartAssembly::~SmartAssembly(){
	//local_segments_
	for ( std::pair< core::Size, SmartSegmentOP > local: local_segments_ ) {
		local.second->set_n_terminal_neighbor( nullptr );
		local.second->set_c_terminal_neighbor( nullptr );
		//These three would have also caused SegmentVector to leak
		local.second->set_n_terminal_parent( nullptr );
		local.second->set_c_terminal_parent( nullptr );
		local.second->set_const_reference_segment( nullptr );
		//If it's a ligand segment, also delete pointers to/from ligand
		if ( std::dynamic_pointer_cast< LigandSegment >( local.second ) != nullptr ) {
			//Should we do this? Are there cases where more than one segment thinks it owns a ligand (i.e. parent and chimaera)? Does it matter if we're deleting the assembly anyway?
			/*
			for( std::pair< core::Size, LigandResidueOP > lig: std::dynamic_pointer_cast< LigandSegment( local.second )->get_owned_ligand_residues() ){
			lig.second->set_owner_segment( nullptr );
			}
			*/
			std::dynamic_pointer_cast< LigandSegment >( local.second )->get_owned_ligand_residues().clear();
		}
	}

	//pdb_segments_ (since these are also local to a given assembly although const)
	for ( std::pair< core::Size, SmartSegmentOP > pdbseg: pdb_segments_ ) {
		pdbseg.second->set_n_terminal_neighbor( nullptr );
		pdbseg.second->set_c_terminal_neighbor( nullptr );
		//These three would have also caused SegmentVector to leak
		pdbseg.second->set_n_terminal_parent( nullptr );
		pdbseg.second->set_c_terminal_parent( nullptr );
		pdbseg.second->set_const_reference_segment( nullptr );
		//If it's a ligand segment, also delete pointers to/from ligand
		if ( std::dynamic_pointer_cast< LigandSegment >( pdbseg.second ) != nullptr ) {
			//Should we do this? Are there cases where more than one segment thinks it owns a ligand (i.e. parent and chimaera)? Does it matter if we're deleting the assembly anyway?
			//These point to the const ligands which won't be owned by any chimaerae, so we should do this here
			for ( std::pair< core::Size, LigandResidueOP > lig: std::dynamic_pointer_cast< LigandSegment >( pdbseg.second )->get_owned_ligand_residues() ) {
				lig.second->set_owner_segment( nullptr );
			}
			std::dynamic_pointer_cast< LigandSegment >( pdbseg.second )->get_owned_ligand_residues().clear();
		}
	}
	//local_ligands_ definitely
	for ( std::pair< core::Size, LigandResidueOP > lig: local_ligands_ ) {
		lig.second->set_owner_segment( nullptr );
		utility::vector1< LigandContactOP > temp;
		temp.clear();
		lig.second->set_contacts( temp );
	}
	//const_ligands_ (go along with pdb_segments_) are const
	//ligand_conformers_ (to be safe)
	for ( std::pair< core::Size, utility::vector1< LigandResidueOP > > ligvec: ligand_conformers_ ) {
		for ( LigandResidueOP lig: ligvec.second ) {
			lig->set_owner_segment( nullptr );
			utility::vector1< LigandContactOP > temp;
			temp.clear();
			lig->set_contacts( temp );
		}
		ligvec.second.clear();
	}
	//last_sampled_ligand_, first_segment_, etc. should have been already deleted b/c they were in local_segments_ and local_ligands_
}

SmartAssembly::SmartAssembly( SmartAssembly const & other) {
	segment_vector_ = other.get_segment_vector();
	pdb_segments( other.const_pdb_segments() );
	window_width_ = 1;
	for ( SmartSegmentCOP current_segment : *segment_vector_ ) {
		if ( current_segment->get_length() >= window_width_ ) {
			if ( !current_segment->is_n_terminus_fixed() && current_segment->get_dssp_code() != 'L' ) {
				n_terminal_segments_.push_back(current_segment->get_segment_id());
			}
			//A segment might be both n-terminal and c-terminal
			if ( !current_segment->is_c_terminus_fixed() && current_segment->get_dssp_code() != 'L' ) {
				c_terminal_segments_.push_back(current_segment->get_segment_id());
			}
		}

	}
	set_partner( other.get_partner() );
	set_partner_ligands( other.get_partner_ligands() );
	//pre-allocate all_basis_pairs_
	all_basis_pairs_ = utility::vector1< std::pair< core::Size, core::Size > >( 10000, std::make_pair( 0, 0 ) );
}



SmartAssemblyOP
SmartAssembly::clone() const {
	return SmartAssemblyOP( new SmartAssembly( *this ) );
}

void
SmartAssembly::set_starting_segment(SmartSegmentOP start_segment, std::string start_node_vital_segments){ // for Append mode
	if ( local_segments_.count( start_segment->get_segment_id() ) == 0 ) {
		add_segment_and_neighbors_to_local_segments( start_segment );
	}
	SmartSegmentOP local_start_segment = local_segments_.at( start_segment->get_segment_id() );
	n_terminal_segment_ = SmartSegment::get_n_most_segment( local_start_segment , false);
	c_terminal_segment_ = SmartSegment::get_c_most_segment( local_start_segment, false);

	//Setting which start node segments are vital: either just terminal or all (default)
	if ( start_node_vital_segments == "terminal" ) {
		n_terminal_segment_->set_is_vital(true);
		c_terminal_segment_->set_is_vital(true);
	} else {
		for ( core::Size i=0; i<=( (c_terminal_segment_->get_segment_id()) - (n_terminal_segment_->get_segment_id()) ); ++i ) {
			core::Size j = n_terminal_segment_->get_segment_id() + i;
			local_segments_.at(j)->set_is_vital(true);
		}
	}

	length_ = 0;
	size_ = 0;
	counter_segment_ = n_terminal_segment_;
	while ( counter_segment_ ) {
		counter_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= counter_segment_->get_length();
		counter_segment_ = counter_segment_->get_c_terminal_neighbor();
	}
}
// getters and setters for the only two pointers the Assembly has to the Segment list. Everything else, we stitch together.
SmartSegmentOP
SmartAssembly::get_n_terminal_segment() const {
	return n_terminal_segment_;
}

void
SmartAssembly::set_n_terminal_segment( SmartSegmentOP new_segment) {
	n_terminal_segment_ = new_segment;
}

SmartSegmentOP
SmartAssembly::get_c_terminal_segment() const {
	return c_terminal_segment_;
}
void
SmartAssembly::set_c_terminal_segment( SmartSegmentOP new_segment) {
	c_terminal_segment_ = new_segment;
}

std::string
SmartAssembly::get_start_node_vital_segments() {
	return start_node_vital_segments_;
}

void
SmartAssembly::set_start_node_vital_segments(std::string start_node_vital_segments) {
	start_node_vital_segments_ = start_node_vital_segments;
}

bool
SmartAssembly::has_segment( core::Size seg_id ){
	SmartSegmentOP current = n_terminal_segment_;
	while ( current ) {
		if ( current->get_segment_id() == seg_id ) {
			return true;
		}
		current = current->get_c_terminal_neighbor();
	}
	return false;
}

SmartSegmentOP
SmartAssembly::get_segment( core::Size seg_id ){
	SmartSegmentOP current = n_terminal_segment_;
	while ( current ) {
		if ( current->get_segment_id() == seg_id ) {
			return current;
		}
		current = current->get_c_terminal_neighbor();
	}
	return nullptr;
}


core::Size
SmartAssembly::get_length() const{ //length in residues
	return length_;
}

void
SmartAssembly::set_length(core::Size new_length){
	length_ = new_length;
}

core::Size
SmartAssembly::get_size() const{ //for size in segments
	return size_;

}

void
SmartAssembly::set_size(core::Size new_size){
	size_ = new_size;
}

hashing::SegmentVectorCOP
SmartAssembly::get_segment_vector() const {
	return segment_vector_;
}

void
SmartAssembly::set_segment_vector(hashing::SegmentVectorCOP segment_vector ){
	segment_vector_ = segment_vector;
}

bool
SmartAssembly::get_last_change_was_n_terminal() const {
	return last_change_was_n_terminal_;
}

void
SmartAssembly::set_last_change_was_n_terminal(bool last_change_n) {
	last_change_was_n_terminal_ = last_change_n;
}

char
SmartAssembly::get_last_change() const {
	return last_change_;
}

void
SmartAssembly::set_last_change(char last_change){
	last_change_=last_change;
}
SmartSegmentOP
SmartAssembly::get_last_chimaera() const {
	return last_chimaera_;
}

BasisPair
SmartAssembly::get_last_chimaera_deleted() const {
	return last_chimaera_deleted_;
}

void
SmartAssembly::set_last_chimaera( SmartSegmentOP last_chimaera ) {
	last_chimaera_=last_chimaera;
}

void
SmartAssembly::set_last_chimaera_deleted( BasisPair last_chimaera_deleted ) {
	last_chimaera_deleted_=last_chimaera_deleted;
}
bool
SmartAssembly::get_output_partner() const {
	return output_partner_;
}
void
SmartAssembly::set_output_partner(bool output_partner){
	output_partner_ = output_partner;
}

bool
SmartAssembly::can_delete() {
	return !(SmartSegment::get_far_n_terminal_parent( n_terminal_segment_ )->is_vital() && SmartSegment::get_far_c_terminal_parent( c_terminal_segment_ )-> is_vital() );
}

bool
SmartAssembly::get_modifiable_terminus(char op) {
	if ( op == 'A' ) {
		if ( modifiable_terminus_ == 'N' ) {
			return true;
		} else if ( modifiable_terminus_ == 'C' ) {
			return false;
		} else {
			return (numeric::random::rg().uniform() < 0.5);
		}
	}
	if ( this->get_n_terminal_segment()->is_vital() ) {
		return false;
	} else if ( this->get_c_terminal_segment()->is_vital() ) {
		return true;
	} else {
		return  (numeric::random::rg().uniform() < 0.5);
	}
}

char
SmartAssembly::get_modifiable_terminus() {
	return modifiable_terminus_;
}

void
SmartAssembly::set_modifiable_terminus(char modifiable_terminus) {
	modifiable_terminus_ = modifiable_terminus;
}

core::pose::PoseOP
SmartAssembly::get_partner() const {
	return partner_;
}

void
SmartAssembly::set_partner(core::pose::PoseOP partner) {
	partner_ = partner;
}

void
SmartAssembly::pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs ){
	pdb_segments_ = pdbsegs;
	//core::Size seg_index = 1;
	for ( std::pair< core::Size, SmartSegmentOP > current_segment: pdb_segments_ ) {
		/*  if( current_segment.second->get_length() >= window_width_ ){
		if(!current_segment.second->is_n_terminus_fixed()){
		n_terminal_segments_.push_back(current_segment.first);
		}
		else if (!current_segment.second->is_c_terminus_fixed()){
		c_terminal_segments_.push_back(current_segment.first);
		}
		}
		*/  //I'm going to assume for now that all ligands are found in pdb segments
		if ( std::dynamic_pointer_cast< LigandSegment >( current_segment.second ) != nullptr ) {
			for ( std::pair< core::Size, LigandResidueOP > ligand: std::dynamic_pointer_cast< LigandSegment >( current_segment.second )->get_const_owned_ligand_residues() ) {
				const_ligands_[ ligand.first ] = LigandResidueCOP( ligand.second );
			}
		}
		//This should also add ligands to local_ligands_
		add_segment_and_neighbors_to_local_segments( current_segment.second );
		//
	}//++seg_index;
	runtime_assert( const_ligands_.size() + partner_ligands_.size() == local_ligands_.size() ); //Make sure all of the const ligands were added to local ligands. NOTE these two will not always be the same size, i.e. in repeat proteins.
}

std::map< core::Size, data_storage::SmartSegmentOP > &
SmartAssembly::pdb_segments(){
	return pdb_segments_;
}

std::map< core::Size, data_storage::SmartSegmentOP > &
SmartAssembly::local_segments(){
	return local_segments_;
}

std::map< core::Size, data_storage::SmartSegmentOP >
SmartAssembly::const_pdb_segments() const{
	return pdb_segments_;
}

core::Size
SmartAssembly::get_random_segment_id( bool n_terminus ){
	core::Size seg_ID;
	core::Size index = 0;
	if ( n_terminus ) {
		while ( index == 0 || index > n_terminal_segments_.size() ) {
			index = std::ceil(n_terminal_segments_.size()*numeric::random::rg().uniform());
		}
		seg_ID = n_terminal_segments_.at( index );
	} else {
		while ( index == 0 || index > c_terminal_segments_.size() ) {
			index = std::ceil(c_terminal_segments_.size()*numeric::random::rg().uniform());
		}
		seg_ID = c_terminal_segments_.at( index );
	}
	return seg_ID;
}

std::string
SmartAssembly::get_forward_assembly() const{
	std::string forward_assembly;
	SmartSegmentOP current_segment = this->get_n_terminal_segment();
	while ( current_segment != nullptr ) {
		forward_assembly = forward_assembly + utility::to_string( current_segment->get_segment_id() );
		forward_assembly = forward_assembly + " - " ;
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return forward_assembly;
}

// can be used to verify equivalency of two assemblies without depending on the chimaeric segment id
std::string
SmartAssembly::get_comprehensive_forward_assembly() const{
	std::string forward_assembly;
	SmartSegmentOP current_segment = this->get_n_terminal_segment();
	while ( current_segment != nullptr ) {
		if ( current_segment->is_chimaeric() ) {
			BasisPair bp = current_segment->get_basis_pair();
			forward_assembly = forward_assembly + "[ " + utility::to_string( bp.first.segment_id() ) + " " + utility::to_string( bp.first.resnum() );
			forward_assembly = forward_assembly+ ", " + utility::to_string( bp.second.segment_id() ) + " " + utility::to_string( bp.second.resnum() ) + " ] - ";
		} else {
			forward_assembly = forward_assembly + utility::to_string( current_segment->get_segment_id() );
			if ( current_segment->is_vital() ) {
				forward_assembly = forward_assembly + "(vital)";
			}
			forward_assembly = forward_assembly + " - " ;
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return forward_assembly;
}



std::string
SmartAssembly::get_reverse_assembly() const{
	std::string reverse_assembly;
	SmartSegmentOP current_segment = this->get_c_terminal_segment();
	utility::vector1< core::Size > seg_ids;
	while ( current_segment ) {
		seg_ids.push_back(current_segment->get_segment_id());
		current_segment = current_segment->get_n_terminal_neighbor();
	}
	while ( seg_ids.size() > 0 ) {
		reverse_assembly = reverse_assembly + utility::to_string( seg_ids.back() );
		reverse_assembly = reverse_assembly + " - ";
		seg_ids.pop_back();
	}
	return reverse_assembly;
}


// can be used to verify equivalency of two assemblies without depending on the chimaeric segment id
std::string
SmartAssembly::get_comprehensive_reverse_assembly() const{
	std::string reverse_assembly;
	SmartSegmentOP current_segment = this->get_c_terminal_segment();
	utility::vector1< std::string > rev_segs;
	while ( current_segment != nullptr ) {
		std::string seg_string;
		if ( current_segment->is_chimaeric() ) {
			BasisPair bp = current_segment->get_basis_pair();
			seg_string =  "[ " + utility::to_string( bp.first.segment_id() ) + " " + utility::to_string( bp.first.resnum() );
			seg_string = seg_string + ", " + utility::to_string( bp.second.segment_id() ) + " " + utility::to_string( bp.second.resnum() ) + " ]";
		} else {
			seg_string = utility::to_string( current_segment->get_segment_id() );
			if ( current_segment->is_vital() ) {
				reverse_assembly = reverse_assembly + "(vital)";
			}
		}
		rev_segs.push_back( seg_string );
		current_segment = current_segment->get_n_terminal_neighbor();
	}
	while ( rev_segs.size() > 0 ) {
		reverse_assembly = reverse_assembly + rev_segs.back() + " - ";
		rev_segs.pop_back();
	}
	return reverse_assembly;
}
//this is pretty much only for rebuilding the lowest assembly at the end.
std::string
SmartAssembly::get_comprehensive_assembly_string() const{
	std::string comprehensive_assembly;
	SmartSegmentOP current_segment = this->get_n_terminal_segment();
	while ( current_segment != nullptr ) {
		if ( current_segment->is_chimaeric() ) {
			comprehensive_assembly = comprehensive_assembly + utility::to_string( current_segment->get_segment_id() );
			comprehensive_assembly = comprehensive_assembly+ " ";
			BasisPair bp = current_segment->get_basis_pair();
			if ( current_segment->get_n_terminal_parent()->get_segment_id() == bp.first.segment_id() ) {
				comprehensive_assembly = comprehensive_assembly + utility::to_string( bp.first.segment_id() ) + " " + utility::to_string( bp.first.resnum() );
				comprehensive_assembly = comprehensive_assembly + " " + utility::to_string( bp.second.segment_id() ) + " " + utility::to_string( bp.second.resnum()) + " ";
			} else if ( current_segment->get_c_terminal_parent()->get_segment_id() == bp.first.segment_id() ) {
				comprehensive_assembly = comprehensive_assembly + utility::to_string( bp.second.segment_id() ) + " " + utility::to_string( bp.second.resnum());
				comprehensive_assembly = comprehensive_assembly + " " + utility::to_string( bp.first.segment_id() ) + " " + utility::to_string( bp.first.resnum() ) + " ";
			}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return comprehensive_assembly;
}
std::map< core::Size, LigandResidueOP >
SmartAssembly::get_local_ligands(){
	return local_ligands_;
}

std::map< core::Size, LigandResidueOP >
SmartAssembly::get_local_ligands() const{
	return local_ligands_;
}

std::map< core::Size, LigandResidueCOP >
SmartAssembly::get_partner_ligands() const{
	return partner_ligands_;
}

void
SmartAssembly::set_partner_ligands( std::map< core::Size, LigandResidueCOP > ligs){
	partner_ligands_ = ligs;
	//Add them to local_ligands_ as well
	for ( std::pair< core::Size, LigandResidueCOP > lig: partner_ligands_ ) {
		if ( local_ligands_.count( lig.first ) == 0 ) {
			local_ligands_[ lig.first ] = lig.second->clone();
		}

	}
}
//end getters/setters


bool
SmartAssembly::add_segment( bool n_terminus ){
	TR.Debug << "Adding random segment & bp" << std::endl;
	return add_segment( n_terminus, 0, 0, 0 );
}

bool
SmartAssembly::add_segment( bool n_terminus, core::Size seg_ID_to_add, core::Size res_ID_1, core::Size res_ID_2 ){
	//runtime_assert( n_terminal_segments_.size() > 0 && c_terminal_segments_.size() > 0 );
	TR.Debug << "Adding Segment" << std::endl;
	TR.Debug << "N term: " << n_terminal_segment_->get_segment_id() << " C term: " << c_terminal_segment_->get_segment_id() << std::endl;
	if ( n_terminus ) {

		//seg1 is n terminal
		//seg2 is c terminal
		TR.Debug << "adding to n terminus" << std::endl;
		segID_1_ = n_terminal_segment_->get_segment_id();
		//We can safely assume the segment we're adding to is already in the segment vector
		runtime_assert( local_segments_.count( segID_1_ ) != 0 );
		//resID_1_ = (n_terminal_segment_->get_length()-window_width_)*(numeric::random::rg().uniform())+window_width_;
		bool chimerization_succeeded = false;
		core::Size num_attempted = 0;
		while ( !chimerization_succeeded ) {
			++num_attempted;
			if ( num_attempted > 1000 ) {
				return false;
			}
			if ( seg_ID_to_add && res_ID_1 && res_ID_2 && ( num_attempted > 1 ) ) {
				TR << "Basis Pair failed: " << segID_1_ << " " << res_ID_1 << " " << seg_ID_to_add << " " << res_ID_2 << std::endl;
				utility_exit_with_message( "provided basis pair is not valid: could not chimerize" ); // ya done messed up
			}

			SmartSegmentCOP seg1;
			SmartSegmentCOP seg2;

			//TR << "assigning seg1" << std::endl;
			//Set up seg1
			// seg1 should be in the assembly, and therefore local segments
			seg1 = local_segments_[ segID_1_ ];

			//TR << "seg1: " << segID_1_ << " is chimaera? " << seg1->is_chimaeric() << std::endl;
			if ( seg1->is_chimaeric() ) {
				//    TR << "seg1 is a chimaera" << std::endl;
				core::Size basis_res;
				if ( seg1->get_basis_pair().first.segment_id() == seg1->get_n_terminal_parent()->get_segment_id() ) {
					//first basis is n_terminal
					basis_res = seg1->get_basis_pair().first.resnum();
				} else {
					basis_res = seg1->get_basis_pair().second.resnum();
				}
				if ( !( window_width_ < basis_res ) ) {
					TR.Debug << "Cannot add to n term" << std::endl;
					return false;
				}
			}
			//seg1 is n terminal
			//seg2 is c terminal

			if ( seg_ID_to_add == 0 ) { //check if seg2 is a specific segment or random
				segID_2_ = get_random_segment_id( !n_terminus );
			} else {
				segID_2_ = seg_ID_to_add;
			}

			TR.Debug << "trying to add segment: " << segID_2_ << " " << num_attempted << std::endl;
			//Now that we have a segment id, let's find where it is.
			if ( local_segments_.count( segID_2_ ) && local_segments_.at( segID_2_ )->is_in_Assembly() ) { //check the assembly first, if it's there we'll have to make a copy
				//    TR <<"Selected segment is already in the assembly"<<std::endl;
				//    TR <<"Making a copy of the segment and all attached segments" << std::endl;
				core::Size new_segment_id;
				if ( pdb_segments_.size() > 0 ) {
					new_segment_id = std::max( ( segment_vector_->size() + pdb_segments_.rbegin()->first ), local_segments_.rbegin()->first ) + 1;
				} else {
					new_segment_id = std::max( segment_vector_->size(), local_segments_.rbegin()->first ) + 1;
				}
				runtime_assert( local_segments_.count( new_segment_id ) == 0 );
				//    TR << "Generated new segment id " << new_segment_id << std::endl;
				//core::Size new_segment_id;
				data_storage::SmartSegmentCOP reference_segment = local_segments_.at( segID_2_ )->get_const_reference_segment();

				//    TR << "Storing Copy to local segments" << std::endl;
				add_segment_and_neighbors_to_local_segments( reference_segment, new_segment_id );
				//    TR << "Done: " << local_segments_.count( new_segment_id ) << std::endl;

				// Let's find the segment we're supposed to use in the chimaera now that it has a new id
				//Set our segID_2_
				core::Size current_segid = new_segment_id;
				//TR << "Finding new_seg_id for chimaera: ";
				SmartSegmentOP local_seg = local_segments_.at( new_segment_id );
				while ( local_seg != nullptr ) {
					if ( reference_segment->get_segment_id() == segID_2_ ) {
						segID_2_ = current_segid;
						break; //we found it
					}
					local_seg = local_seg->get_c_terminal_neighbor();
					reference_segment = reference_segment->get_c_terminal_neighbor();
					++current_segid;
				}
				//    TR << segID_2_ << std::endl;
				seg2 = local_segments_.at( segID_2_ );
			} else { // end it's already in assembly //it's not in local segements, so let's find it.
				//    TR << "segment not in assembly, finding it...." << std::endl;
				if ( local_segments_.count( segID_2_ ) && local_segments_.at( segID_2_ )->get_const_reference_segment() != nullptr ) {
					seg2 = local_segments_.at( segID_2_ )->get_const_reference_segment();
				} else if (  local_segments_.count( segID_2_ ) ) {
					//We need to pick a new seg2. This one looks like a chimaera or otherwise buggy segment
					TR.Debug << "Const reference segment for seg2 not assigned!" << std::endl;
					continue;
				} else if ( segID_2_ <= segment_vector_->size() ) {
					seg2 = segment_vector_->at( segID_2_ );
				} else {
					seg2 = pdb_segments_.at( segID_2_ );
				}
				/*
				TR << "Found. Getting n most neighbor." << std::endl;
				// now that we have it, add it (and its neighbors) to local segments
				SmartSegmentCOP current_seg2 = SmartSegment::get_n_most_segment( seg2, false );
				*/
				//    TR << "calling add to local" << std::endl;
				//add_segment_and_neighbors_to_local_segments( current_seg2 );
				add_segment_and_neighbors_to_local_segments( seg2 );
				//    TR << "Done: " << local_segments_.count( segID_2_ ) << std::endl;
				//let's reset our seg2
				seg2 = local_segments_[ segID_2_ ];
			}
			//Okay, now that our segments are set up, let's figure out the basis pairs

			std::pair< bool, core::Size>  it_results = iterate_over_basis_pairs( seg1, seg2, n_terminus );
			if ( it_results.first == false ) {
				return false;
			}
			core::Size max_bp_index = it_results.second;
			TR.Debug << "max_bp_index: " << max_bp_index << std::endl;
			if ( max_bp_index == 0 && num_attempted > 1000 ) {
				return false;
			}
			if ( res_ID_1 && res_ID_2 ) { // are there specific basis pairs that we need to use?
				if ( seg1->is_chimaeric() && !( res_ID_1 <= seg1->get_basis_pair().first.resnum() ) ) {
					utility_exit_with_message( "provided basis pair is not valid: it will chimerize over original segment." );// ya done messed up
				}
				bool matched = false;
				std::pair< core::Size, core::Size > match_pair = std::make_pair( res_ID_1, res_ID_2 );
				for ( std::pair< core::Size, core::Size > basis_res_pair : all_basis_pairs_ ) { // see if the pair needed is in all_basis_pairs
					if ( basis_res_pair == match_pair ) {
						matched = true;
						break;
					}
				}
				if ( !matched ) {
					TR << "Basis Pair failed: " << segID_1_ << " " << res_ID_1 << " " << segID_2_ << " " << res_ID_2 << std::endl;
					utility_exit_with_message( "provided basis pair is not valid" ); // ya done messed up
				}
				//since we've found them, assign them.
				resID_1_ = res_ID_1;
				resID_2_ = res_ID_2;
			} else { // we're gonna use a random basis pair
				if ( max_bp_index == 0 ) {
					continue;//can't hash segments together. Find a new one.
				}
				core::Size bp_index = 0;
				while ( bp_index == 0 ) {
					bp_index = std::ceil( numeric::random::rg().uniform() * max_bp_index );
				}
				resID_1_ = all_basis_pairs_[ bp_index ].first;
				resID_2_ = all_basis_pairs_[ bp_index ].second;
			}
			if ( resID_1_ == 0 || resID_2_ == 0 ) {
				utility_exit_with_message( "Something wrong with iterate_over_basis_pairs" );
			}

			//now that we have our basis residues, lets make basis pairs and try to chimerize.
			TR.Debug << "Adding segment " << segID_2_ << " to segment " << segID_1_ << std::endl;
			TR.Debug << "Choosing residue " << resID_1_ << " of " << local_segments_.at(segID_1_)->get_length()-window_width_ << " available residues on first segment." << std::endl;
			TR.Debug << "Choosing residue " << resID_2_ << " of " << segment_vector_->at(segID_2_)->get_length()-window_width_ << " available residues on second segment." << std::endl;
			data_storage::Basis first_basis( segID_1_,resID_1_ );
			data_storage::Basis second_basis( segID_2_,resID_2_ );
			data_storage::BasisPair basis_pair = std::make_pair(first_basis,second_basis);
			chimerization_succeeded = chimerize(basis_pair, n_terminus);
		}
		//it's been chimerized, so let's add it in.
		length_-= n_terminal_segment_->get_length();
		size_--;
		//n_terminal_segment_->set_is_in_Assembly( false );//having this set to false will cause infinitely looping assemblies
		//This should replace the old N-terminal segment with a chimaera
		if ( c_terminal_segment_ == n_terminal_segment_ ) {
			c_terminal_segment_ = last_chimaera_;
		}
		if ( n_terminal_segment_->get_c_terminal_neighbor() ) {
			n_terminal_segment_ = n_terminal_segment_->get_c_terminal_neighbor()->get_n_terminal_neighbor();
		} else {
			n_terminal_segment_ = last_chimaera_;
		}
		//this looks strange, but the parent doesn't point to it's own chimaera, but the cterm neighbor points back to the chimaera.
		runtime_assert( n_terminal_segment_->is_chimaeric() );
		n_terminal_segment_->get_n_terminal_parent()->set_is_in_Assembly( true );
		//last_chimaera_ = n_terminal_segment_; This has been moved to chimerize
		first_segment_ = n_terminal_segment_;//first chimaeric segment, not the n most segment in the assembly
		while ( n_terminal_segment_->get_n_terminal_neighbor() ) {
			//   TR << "Adding segment " << n_terminal_segment_->get_segment_id() << std::endl;
			n_terminal_segment_->set_is_in_Assembly(true);
			size_++;
			length_+= n_terminal_segment_->get_length();
			//now that we added this segment, make sure everything's kosher
			n_terminal_segment_ = n_terminal_segment_->get_n_terminal_neighbor();
			if ( n_terminal_segment_->get_segment_id() == first_segment_->get_segment_id() ) {
				TR << "segment " << n_terminal_segment_->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
				utility_exit_with_message("infinite loop!");
			}
		}
		//adding the last, n most segment
		//  TR << "Adding segment " << n_terminal_segment_->get_segment_id() << std::endl;
		n_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= n_terminal_segment_->get_length();
		//TR <<" N-terminal segment is now " << n_terminal_segment_->get_segment_id();
	} else { //end n_terminal add//adding to the c-terminus
		TR.Debug << "Adding to c_terminus" << std::endl;
		segID_1_ = c_terminal_segment_->get_segment_id();
		bool chimerization_succeeded = false;
		core::Size num_attempted = 0;
		while ( !chimerization_succeeded ) {
			++num_attempted;
			if ( num_attempted > 1000 ) {
				return false;
			}
			if ( seg_ID_to_add && ( num_attempted > 1 ) ) {
				TR << "Basis Pair failed: " << segID_1_ << " " << res_ID_1 << " " << seg_ID_to_add << " " << res_ID_2 << std::endl;
				utility_exit_with_message( "provided basis pair is not valid: could not chimerize" ); // ya done messed up
			}
			SmartSegmentCOP seg1;
			SmartSegmentCOP seg2;

			// seg1 should be in the assembly, and therefore local segments
			//TR << "assigning seg1" << std::endl;
			seg1 = local_segments_[ segID_1_ ];
			//TR << "seg1: " << segID_1_ << " is chimaera? " << seg1->is_chimaeric() << std::endl;
			if ( seg1->is_chimaeric() ) {
				//    TR << "seg1 is a chimaera" << std::endl;
				core::Size basis_res;
				if ( seg1->get_basis_pair().first.segment_id() == seg1->get_n_terminal_parent()->get_segment_id() ) {
					//first basis is n_terminal
					basis_res = seg1->get_basis_pair().first.resnum();
				} else {
					basis_res = seg1->get_basis_pair().second.resnum();
				}
				if ( !( window_width_ < seg1->get_length() - basis_res ) ) {
					TR.Debug << "Cannot add to c term" << std::endl;
					return false;
				}
			}

			if ( seg_ID_to_add == 0 ) { //check if seg2 is a specific segment or random
				segID_2_ = get_random_segment_id( !n_terminus );
			} else {
				segID_2_ = seg_ID_to_add;
			}
			TR.Debug << "trying to add segment: " << segID_2_ << std::endl;
			//seg1 is c terminal
			//seg2 is n terminal




			//Add the new segment and all attached segments to the segment vector
			if ( local_segments_.count( segID_2_ ) && local_segments_.at( segID_2_ )->is_in_Assembly() == true ) {
				//    TR <<"Selected segment is already in the assembly"<<std::endl;
				//    TR <<"Making a copy of the segment" << std::endl;
				//We'll need a segment ID that isn't already taken by something else
				//Choose 1 greater than max( segvec size + pdbseg size, largest key in local_segments_ )
				//This works because std::maps are sorted by key value!
				core::Size new_segment_id;
				if ( pdb_segments_.size() > 0 ) {
					new_segment_id = std::max( ( segment_vector_->size() + pdb_segments_.rbegin()->first ), local_segments_.rbegin()->first ) + 1;
				} else {
					new_segment_id = std::max( segment_vector_->size(), local_segments_.rbegin()->first ) + 1;
				}
				if ( local_segments_.at( segID_2_ )->get_const_reference_segment() == nullptr ) {
					TR.Debug << "Const reference for seg2 not assigned!" << std::endl;
					continue;
				}
				data_storage::SmartSegmentCOP reference_segment = local_segments_.at( segID_2_ )->get_const_reference_segment();
				//    TR << "Storing Copy to local segments" << std::endl;
				add_segment_and_neighbors_to_local_segments( reference_segment, new_segment_id );
				//    TR << "Done: " << local_segments_.count( new_segment_id ) << std::endl;
				//Set our segID_2_
				core::Size current_segid = new_segment_id;
				//    TR << "Checking for it in local_segments" << std::endl;
				SmartSegmentOP local_seg = local_segments_.at( new_segment_id );
				//    TR << "It's there: " << local_seg->get_segment_id() << std::endl;
				while ( local_seg != nullptr ) {
					if ( reference_segment->get_segment_id() == segID_2_ ) {
						segID_2_ = current_segid;
						break; //we found it!
					}
					local_seg = local_seg->get_c_terminal_neighbor();
					reference_segment = reference_segment->get_c_terminal_neighbor();
					++current_segid;
				}
				seg2 = local_segments_[ segID_2_ ];
			} else { // end already in assembly //it's not in local segements, so let's find it.
				//    TR << "segment not in assembly, finding it...." << std::endl;
				if ( local_segments_.count( segID_2_ ) ) {
					seg2 = local_segments_.at( segID_2_ )->get_const_reference_segment();
				} else if ( segID_2_ <= segment_vector_->size() ) {
					seg2 = segment_vector_->at( segID_2_ );
				} else {
					seg2 = pdb_segments_.at( segID_2_ );
				}
				//    TR << "calling add to local" << std::endl;
				add_segment_and_neighbors_to_local_segments( seg2 );
				//    TR << "Done: " << local_segments_.count( segID_2_ ) << std::endl;
				//let's reset our seg2
				seg2 = local_segments_[ segID_2_ ];
			}
			//Okay, now that our segments are set up, let's figure out the basis pairs
			std::pair< bool, core::Size>  it_results = iterate_over_basis_pairs( seg1, seg2, n_terminus );
			if ( it_results.first == false ) {
				return false;
			}
			core::Size max_bp_index = it_results.second;
			TR.Debug << "max_bp_index: " << max_bp_index << std::endl;
			if ( max_bp_index == 0 /*&& num_attempted > 1000*/ ) {
				return false;
			}
			if ( res_ID_1 && res_ID_2 ) { //check if we're using specific basis residues
				if ( seg1->is_chimaeric() && !( res_ID_1 >= seg1->get_basis_pair().second.resnum() ) ) {
					utility_exit_with_message( "provided basis pair is not valid: it will chimerize over original segment." );// ya done messed up
				}
				bool matched = false;
				std::pair< core::Size, core::Size > match_pair = std::make_pair( res_ID_1, res_ID_2 );
				for ( std::pair< core::Size, core::Size > basis_res_pair : all_basis_pairs_ ) { //check if the given basis residues are valid
					if ( basis_res_pair == match_pair ) {
						matched = true;
						break;
					}
				}
				if ( !matched ) {
					utility_exit_with_message( "provided basis pair is not valid" );// ya done messed up
				}
				resID_1_ = res_ID_1;
				resID_2_ = res_ID_2;
			} else { // we're using a random basis pair
				if ( max_bp_index == 0 ) {
					continue;//can't hash segments together. Find a new one. This will continue the while(!chimerize_successful) loop
				}
				core::Size bp_index = 0;
				while ( bp_index == 0 ) {
					bp_index = std::ceil( numeric::random::rg().uniform() * max_bp_index );
				}
				resID_1_ = all_basis_pairs_[ bp_index ].first;
				resID_2_ = all_basis_pairs_[ bp_index ].second;
			}
			if ( resID_1_ == 0 || resID_2_ == 0 ) {
				utility_exit_with_message( "Something wrong with iterate_over_basis_pairs" );
			}

			//now that we have our basis residues, lets make basis pairs and try to chimerize.
			//Add the new segment and all attached segments to the segment vector
			//TR << "Adding segment " << segID_2_ << " to segment " << segID_1_ << std::endl;
			//TR << "Choosing residue " << resID_1_ << " of " << segment_vector_->at(segID_1_)->get_length()-window_width_ << " available residues on first segment." << std::endl;
			//TR << "Choosing residue " << resID_2_ << " of " << segment_vector_->at(segID_2_)->get_length()-window_width_ << " available residues on second segment." << std::endl;

			data_storage::Basis first_basis( segID_1_,resID_1_ );
			data_storage::Basis second_basis( segID_2_,resID_2_ );
			chimerization_succeeded = chimerize(std::make_pair(first_basis,second_basis), n_terminus);
		} // end while(!chimerize_successful)
		//it's been chimerized, so let's add it in.
		length_-= c_terminal_segment_->get_length();
		size_--;
		//c_terminal_segment_->set_is_in_Assembly( false );// having this set to false will cause infinitely looping assemblies.
		//This should replace the old N-terminal segment with a chimaera
		if ( c_terminal_segment_ == n_terminal_segment_ ) {
			n_terminal_segment_ = last_chimaera_;
		}
		if ( c_terminal_segment_->get_n_terminal_neighbor() ) {
			c_terminal_segment_ = c_terminal_segment_->get_n_terminal_neighbor()->get_c_terminal_neighbor();
		} else {
			c_terminal_segment_ = last_chimaera_;
		}
		//this looks strange, but the parent doesn't point to it's own chimaera, but the nterm neighbor points to the chimaera.
		runtime_assert( c_terminal_segment_->is_chimaeric() );
		c_terminal_segment_->get_c_terminal_parent()->set_is_in_Assembly( true );
		//last_chimaera_ = c_terminal_segment_; //This has been moved to Chimerize
		first_segment_ = c_terminal_segment_;
		while ( c_terminal_segment_->get_c_terminal_neighbor() ) {
			//   TR << "adding segment " << c_terminal_segment_->get_segment_id();
			c_terminal_segment_->set_is_in_Assembly(true);
			size_++;
			length_+= c_terminal_segment_->get_length();
			//now that we added this segment, make sure everything's kosher
			c_terminal_segment_ = c_terminal_segment_->get_c_terminal_neighbor();
			if ( c_terminal_segment_->get_segment_id() == first_segment_->get_segment_id() ) {
				TR << "segment " << c_terminal_segment_->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
				utility_exit_with_message("infinite loop!");
			}
		}
		// adding the last, c most segment.
		//  TR << "Adding segment " << c_terminal_segment_->get_segment_id() << std::endl;
		c_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= c_terminal_segment_->get_length();
		//TR <<" C-terminal segment is now " << c_terminal_segment_->get_segment_id();
	} // end cterm add
	// TR << "Does the Assembly have chain breaks?: " << !this->is_continuous() << std::endl;
	last_change_was_n_terminal_ = n_terminus;
	last_change_ = 'A';
	if ( n_terminal_segment_ != c_terminal_segment_ ) {
		if ( n_terminal_segment_->is_vital() ) {
			n_terminal_segment_ = c_terminal_segment_;
			while ( n_terminal_segment_->is_n_terminus_fixed() ) {
				n_terminal_segment_ = n_terminal_segment_->get_n_terminal_neighbor();
			}
		} else if ( c_terminal_segment_->is_vital() ) {
			c_terminal_segment_ = n_terminal_segment_;
			while ( c_terminal_segment_->is_c_terminus_fixed() ) {
				c_terminal_segment_ = c_terminal_segment_->get_c_terminal_neighbor();
			}
		}
	}
	return true;
}

SmartSegmentOP
SmartAssembly::recurse_revert_far_n_terminal_parent( SmartSegmentOP current_segment ){
	if ( !current_segment->is_chimaeric() ) {
		utility_exit_with_message( "recurse_revert_far_n_terminal_parent recieved non-chiameric segment" );
	}
	current_segment->set_is_in_Assembly( false );
	current_segment->get_n_terminal_parent()->set_is_in_Assembly( false );
	current_segment->get_c_terminal_parent()->set_is_in_Assembly( false );

	if ( std::dynamic_pointer_cast< LigandSegment >( current_segment ) != nullptr ) {
		unchimerize_ligand_segment( std::dynamic_pointer_cast< LigandSegment >( current_segment ) );
	}

	if ( current_segment->get_n_terminal_parent()->is_chimaeric() ) {
		SmartSegmentOP n_segment = this->recurse_revert_far_n_terminal_parent( current_segment->get_n_terminal_parent() );
		SmartSegmentOP c_segment = current_segment->get_c_terminal_parent();
		core::Size old_n_term_basis_res;
		core::Size c_basis_res;
		if ( current_segment->get_basis_pair().first.segment_id() == current_segment->get_n_terminal_parent()->get_segment_id() ) {
			//first basis is n_terminal
			old_n_term_basis_res = current_segment->get_basis_pair().first.resnum();
			c_basis_res = current_segment->get_basis_pair().second.resnum();
		} else {
			//second basis is n_terminal
			old_n_term_basis_res = current_segment->get_basis_pair().second.resnum();
			c_basis_res = current_segment->get_basis_pair().first.resnum();
		}
		core::Size old_n_parent_n_basis_res;
		core::Size old_n_parent_c_basis_res;
		//If the original segment is the far N-terminal parent
		if ( current_segment->get_n_terminal_parent()->get_basis_pair().first.segment_id() == current_segment->get_n_terminal_parent()->get_n_terminal_parent()->get_segment_id() ) {
			old_n_parent_n_basis_res = current_segment->get_n_terminal_parent()->get_basis_pair().first.resnum();
			old_n_parent_c_basis_res = current_segment->get_n_terminal_parent()->get_basis_pair().second.resnum();
		} else { //If the original segment is n_segment
			old_n_parent_n_basis_res = current_segment->get_n_terminal_parent()->get_basis_pair().second.resnum();
			old_n_parent_c_basis_res = current_segment->get_n_terminal_parent()->get_basis_pair().first.resnum();
		}
		//core::Size n_basis_res = old_n_term_basis_res - old_n_parent_n_basis_res + old_n_parent_c_basis_res;
		core::Size n_basis_res = old_n_term_basis_res + old_n_parent_c_basis_res - old_n_parent_n_basis_res;
		bool n_terminus = false; //n_term is set. add to c term

		if ( !this->add_segment( n_terminus, c_segment->get_segment_id(), n_basis_res, c_basis_res ) ) {
			TR << "Bad BasisPair!! Seg 1: " << n_segment->get_segment_id() << " res 1: " << n_basis_res << std::endl;
			TR << "\t\t\t\t   Seg 2: " << c_segment->get_segment_id() << " res 2: " << c_basis_res << std::endl;
			utility_exit_with_message( "recurse_revert_far_n_terminal_parent could not chimerize segment" );
		}
		/*
		//core::Size n_basis_res = n_segment->get_length() + current_segment->get_basis_pair().first.resnum() - current_segment->get_n_terminal_parent()->get_length();
		Basis n_basis( n_segment->get_segment_id(), n_basis_res );
		bool n_terminus = false;
		BasisPair bp;
		if( c_segment->is_vital() ){
		//c_terminal_segment_ isn't set. will this break everything?
		bp = std::make_pair( c_basis, n_basis );
		}
		else {
		bp = std::make_pair( n_basis, c_basis );
		}

		// revert current chimaeric links before rechimerizing.
		current_segment->set_is_in_Assembly( false );
		SmartSegment::link_to(current_segment->get_n_terminal_parent()->get_n_terminal_neighbor(), current_segment->get_n_terminal_parent());
		SmartSegment::link_to(current_segment->get_c_terminal_parent(), current_segment->get_c_terminal_parent()->get_c_terminal_neighbor());
		TR << "re-chimerizing" << std::endl;
		if( !chimerize( bp, n_terminus ) ){
		TR << "Bad BasisPair!! Seg 1: " << bp.first.segment_id() << " res 1: " << bp.first.resnum() << std::endl;
		TR << "       Seg 2: " << bp.second.segment_id() << " res 2: " << bp.second.resnum() << std::endl;
		utility_exit_with_message( "recurse_revert_far_n_terminal_parent could not chimerize segment" );
		}

		length_ -= n_segment->get_length();
		n_terminal_segment_ = last_chimaera_;
		c_terminal_segment_ = last_chimaera_;
		SmartSegmentOP first_segment = c_terminal_segment_;
		while(c_terminal_segment_->get_c_terminal_neighbor() != nullptr ){
		TR << "Adding segment " << c_terminal_segment_->get_segment_id() << std::endl;
		c_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= c_terminal_segment_->get_length();
		//now that we added this segment, make sure everything's kosher
		c_terminal_segment_ = c_terminal_segment_->get_c_terminal_neighbor();
		if(c_terminal_segment_->get_segment_id() == first_segment_->get_segment_id()){
		TR << "segment " << c_terminal_segment_->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
		utility_exit_with_message("infinite loop!");
		}
		}
		TR << "Adding segment " << c_terminal_segment_->get_segment_id() << std::endl;
		c_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= c_terminal_segment_->get_length();
		*/
		return n_terminal_segment_;

	} else { // found n_term parent we want to remove. Return c_term parent for re-chimerizing
		//TR << "found far n_term parent: " << current_segment->get_n_terminal_parent()->get_segment_id() << std::endl;
		//TR << "returning c_term parent: " << current_segment->get_c_terminal_parent()->get_segment_id() << std::endl;
		n_terminal_segment_ = current_segment->get_c_terminal_parent();
		c_terminal_segment_ = n_terminal_segment_;
		last_reverted_parent_chimaera_ = current_segment;

		//n_term parent is being removed. Current segment will no longer exist in assembly. Revert chimaera.
		current_segment->get_n_terminal_parent()->set_is_in_Assembly( false );

		//delete current assembly before rebuilding
		current_segment->set_is_in_Assembly( false );
		//length_ -= current_segment->get_length();
		--size_;
		length_ += n_terminal_segment_->get_length();
		size_ += 1;

		SmartSegment::link_to(current_segment->get_n_terminal_parent()->get_n_terminal_neighbor(), current_segment->get_n_terminal_parent());
		SmartSegment::link_to(current_segment->get_c_terminal_parent(), current_segment->get_c_terminal_parent()->get_c_terminal_neighbor());
		n_terminal_segment_->set_is_in_Assembly( true );

		//reset length since we're rebuilding the assembly
		//length_ = n_terminal_segment_->get_length();
		//size_ = 1;

		return n_terminal_segment_;
	}
}

SmartSegmentOP
SmartAssembly::recurse_revert_far_c_terminal_parent( SmartSegmentOP current_segment ){
	if ( !current_segment->is_chimaeric() ) {
		utility_exit_with_message( "recurse_revert_far_c_terminal_parent recieved non-chiameric segment" );
	}
	current_segment->set_is_in_Assembly( false );
	current_segment->get_n_terminal_parent()->set_is_in_Assembly( false );
	current_segment->get_c_terminal_parent()->set_is_in_Assembly( false );
	if ( std::dynamic_pointer_cast< LigandSegment >( current_segment ) != nullptr ) {
		unchimerize_ligand_segment( std::dynamic_pointer_cast< LigandSegment >( current_segment ) );
	}
	if ( current_segment->get_c_terminal_parent()->is_chimaeric() ) {
		SmartSegmentOP c_segment = this->recurse_revert_far_c_terminal_parent( current_segment->get_c_terminal_parent() );
		SmartSegmentOP n_segment = current_segment->get_n_terminal_parent();
		core::Size c_basis_res;
		core::Size n_basis_res;

		if ( current_segment->get_basis_pair().first.segment_id() == current_segment->get_n_terminal_parent()->get_segment_id() ) {
			//first basis is n_terminal
			c_basis_res = current_segment->get_basis_pair().second.resnum();
			n_basis_res = current_segment->get_basis_pair().first.resnum();
		} else {
			//second basis is n_terminal
			c_basis_res = current_segment->get_basis_pair().first.resnum();
			n_basis_res = current_segment->get_basis_pair().second.resnum();
		}
		bool n_terminus = true;
		if ( !this->add_segment( n_terminus, n_segment->get_segment_id(), c_basis_res, n_basis_res ) ) {
			TR << "Bad BasisPair!! Seg 1: " << c_segment->get_segment_id() << " res 1: " << c_basis_res << std::endl;
			TR << "\t\t\t\t   Seg 2: " << n_segment->get_segment_id() << " res 2: " << n_basis_res << std::endl;
			utility_exit_with_message( "recurse_revert_far_c_terminal_parent could not chimerize segment" );
		}

		/*
		bool n_terminus = true;
		BasisPair bp;
		if( n_segment->is_in_Assembly() ){
		//n_terminal_segment_ isn't set. will this break everything?
		bp = std::make_pair( n_basis, c_basis );
		}
		else {
		bp = std::make_pair( c_basis, n_basis );
		}

		SmartSegmentOP nterm_temp_deleted = n_segment->get_n_terminal_neighbor();
		while( nterm_temp_deleted != nullptr ){
		nterm_temp_deleted->set_is_in_Assembly( false );
		length_ -= nterm_temp_deleted->get_length();
		--size_;
		nterm_temp_deleted = nterm_temp_deleted->get_n_terminal_neighbor();
		}

		// revert current chimaeric links before rechimerizing.
		current_segment->set_is_in_Assembly( false );
		SmartSegment::link_to(current_segment->get_n_terminal_parent()->get_n_terminal_neighbor(), current_segment->get_n_terminal_parent());
		SmartSegment::link_to(current_segment->get_c_terminal_parent(), current_segment->get_c_terminal_parent()->get_c_terminal_neighbor());

		TR << "re-chimerizing" << std::endl;
		if( !chimerize( bp, n_terminus ) ){
		TR << "Bad BasisPair!! Seg 1: " << bp.first.segment_id() << " res 1: " << bp.first.resnum() << std::endl;
		TR << "       Seg 2: " << bp.second.segment_id() << " res 2: " << bp.second.resnum() << std::endl;
		utility_exit_with_message( "recurse_revert_far_c_terminal_parent could not chimerize segment" );
		}
		length_ -= c_segment->get_length();
		n_terminal_segment_ = last_chimaera_;
		c_terminal_segment_ = last_chimaera_;
		SmartSegmentOP first_segment = n_terminal_segment_;
		while(n_terminal_segment_->get_n_terminal_neighbor() != nullptr ){
		TR << "Adding segment " << n_terminal_segment_->get_segment_id() << std::endl;
		n_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= n_terminal_segment_->get_length();
		//now that we added this segment, make sure everything's kosher
		n_terminal_segment_ = n_terminal_segment_->get_n_terminal_neighbor();
		if(n_terminal_segment_->get_segment_id() == first_segment_->get_segment_id()){
		TR << "segment " << n_terminal_segment_->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
		utility_exit_with_message("infinite loop!");
		}
		}
		TR << "Adding segment " << n_terminal_segment_->get_segment_id() << std::endl;
		n_terminal_segment_->set_is_in_Assembly(true);
		size_++;
		length_+= n_terminal_segment_->get_length();
		*/
		return c_terminal_segment_;
	} else { // found c_term parent we want to remove. Return n_term parent for re-chimerizing
		TR.Debug << "found far c_term parent" << std::endl;
		c_terminal_segment_ = current_segment->get_n_terminal_parent();
		n_terminal_segment_ = c_terminal_segment_;
		last_reverted_parent_chimaera_ = current_segment;

		//c_term is being removed. Current segment will no longer exist in assembly. Revert chimaera.
		current_segment->get_c_terminal_parent()->set_is_in_Assembly( false );

		current_segment->set_is_in_Assembly( false );
		//length_ -= current_segment->get_length();
		--size_;
		length_ += c_terminal_segment_->get_length();
		size_ += 1;
		SmartSegment::link_to(current_segment->get_n_terminal_parent()->get_n_terminal_neighbor(), current_segment->get_n_terminal_parent());
		SmartSegment::link_to(current_segment->get_c_terminal_parent(), current_segment->get_c_terminal_parent()->get_c_terminal_neighbor());
		c_terminal_segment_->set_is_in_Assembly( true );

		//reset length since we're rebuilding the assembly
		//length_ = c_terminal_segment_->get_length();
		//size_ = 1;

		return c_terminal_segment_;
	}
}

bool
SmartAssembly::delete_segment(bool n_terminus){
	TR.Debug << "Deleting Segment" << std::endl;
	if ( n_terminus ) {
		TR.Debug << "Deleting segment from N terminus." << std::endl;
		if ( SmartSegment::get_far_n_terminal_parent( n_terminal_segment_ )->is_vital() ) {
			return false;
		}
		//Delete until we reach a chimaera
		while ( !n_terminal_segment_->is_chimaeric() && !n_terminal_segment_->is_vital() ) {
			//   if (n_terminal_segment_->is_vital()){
			//    TR << "non-terminal segment " << n_terminal_segment_->get_segment_id() << " is vital!" << std::endl;
			//   }
			//   TR << "removing segment " << n_terminal_segment_->get_segment_id() << std::endl;

			if ( c_terminal_segment_ == n_terminal_segment_ ) {
				TR << "This shouldn't happen! Problem with keeping vital segments." << std::endl;
				break;
			} else if ( n_terminal_segment_->get_c_terminal_neighbor() == nullptr ) {
				TR << "This shouldn't happen! Connections between segments are broken!" << std::endl;
				utility_exit_with_message( "Problem found with C-terminal connection of N-terminal segment" );
			}
			n_terminal_segment_->set_is_in_Assembly(false);
			size_--;
			length_-= n_terminal_segment_->get_length();
			n_terminal_segment_ = n_terminal_segment_->get_c_terminal_neighbor();
		}
		//  TR << "N-terminal segment is now " << n_terminal_segment_->get_segment_id();
		// current_segment is now a chimaera
		//Replace n-terminal_segment with its C-terminal parent
		if ( n_terminal_segment_->is_chimaeric() ) {
			//Reset the ligand's contacts and owner segment if applicable
			//if( n_terminal_segment_->get_name() == "LigandSegment" ){
			if ( n_terminal_segment_->get_n_terminal_parent()->is_chimaeric() ) {
				//    TR << "Calling recurse revert for n terminal parent" << std::endl;
				length_ -= n_terminal_segment_->get_length();
				//TODO
				SmartSegmentOP cterm_deleted = n_terminal_segment_->get_c_terminal_parent()->get_c_terminal_neighbor();
				while ( cterm_deleted != nullptr ) {
					cterm_deleted->set_is_in_Assembly( false );
					length_ -= cterm_deleted->get_length();
					--size_;
					cterm_deleted = cterm_deleted->get_c_terminal_neighbor();
				}
				//This includes the unchimerize_ligand_segment command
				recurse_revert_far_n_terminal_parent( n_terminal_segment_ );
				last_chimaera_deleted_ = last_reverted_parent_chimaera_->get_basis_pair();
			} else {
				if ( std::dynamic_pointer_cast< LigandSegment >( n_terminal_segment_ ) != nullptr ) {
					unchimerize_ligand_segment( std::dynamic_pointer_cast< LigandSegment >( n_terminal_segment_ ) );
				}
				n_terminal_segment_->set_is_in_Assembly( false );
				n_terminal_segment_->get_n_terminal_parent()->set_is_in_Assembly( false );
				last_chimaera_deleted_ = n_terminal_segment_->get_basis_pair();
				//    TR << "n_terminal_segment id: " << n_terminal_segment_->get_segment_id() << ". reassigned segment id: "
				//      << n_terminal_segment_->get_c_terminal_parent()->get_segment_id() << ". is it chimaeric: "
				//      << n_terminal_segment_->get_c_terminal_parent()->is_chimaeric() << std::endl;
				SmartSegment::link_to(n_terminal_segment_->get_n_terminal_parent()->get_n_terminal_neighbor(),n_terminal_segment_->get_n_terminal_parent());
				SmartSegment::link_to(n_terminal_segment_->get_c_terminal_parent(),n_terminal_segment_->get_c_terminal_parent()->get_c_terminal_neighbor());
				length_ -= n_terminal_segment_->get_length();
				if ( n_terminal_segment_->get_segment_id() == c_terminal_segment_->get_segment_id() ) { // check if we are reverting to a single segment
					c_terminal_segment_ = n_terminal_segment_->get_c_terminal_parent();
				}
				n_terminal_segment_ = n_terminal_segment_->get_c_terminal_parent();
				n_terminal_segment_->set_is_in_Assembly( true );
				length_ += n_terminal_segment_->get_length();
			}
		}
		//Update the length
		//  TR << "N-terminal segment is now " << n_terminal_segment_->get_segment_id() << std::endl;
	} else {
		TR.Debug << "Deleting segment from C terminus." << std::endl;
		if ( SmartSegment::get_far_c_terminal_parent( c_terminal_segment_ )->is_vital() ) {
			return false;
		}
		while ( !c_terminal_segment_->is_chimaeric() && !c_terminal_segment_->is_vital() ) {
			//if (c_terminal_segment_->is_vital()){
			//  TR << "non-terminal segment " << c_terminal_segment_->get_segment_id() << " is vital!" << std::endl;
			//}
			//   TR << "removing segment " << c_terminal_segment_->get_segment_id() << std::endl;
			if ( c_terminal_segment_ == n_terminal_segment_ ) {
				TR << "This shouldn't happen! Problem with keeping vital segments." << std::endl;
				break;
			} else if ( c_terminal_segment_->get_n_terminal_neighbor() == nullptr ) {
				TR << "This shouldn't happen! Connections between segments are broken!" << std::endl;
				utility_exit_with_message( "Problem found with N-terminal connection of c-terminal segment" );
			}
			c_terminal_segment_->set_is_in_Assembly(false);
			size_--;
			length_-= c_terminal_segment_->get_length();
			c_terminal_segment_ = c_terminal_segment_->get_n_terminal_neighbor();
		}
		//  TR << "C-terminal segment is now " << c_terminal_segment_->get_segment_id();
		//current_segment is now a chimaera
		if ( c_terminal_segment_->is_chimaeric() ) {
			if ( c_terminal_segment_->get_c_terminal_parent()->is_chimaeric() ) {
				//    TR << "Calling recurse revert for c terminal parent" << std::endl;
				length_ -= c_terminal_segment_->get_length();
				SmartSegmentOP nterm_deleted = c_terminal_segment_->get_n_terminal_parent()->get_n_terminal_neighbor();
				while ( nterm_deleted != nullptr ) {
					nterm_deleted->set_is_in_Assembly( false );
					length_ -= nterm_deleted->get_length();
					--size_;
					nterm_deleted = nterm_deleted->get_n_terminal_neighbor();
				}
				SmartSegmentOP new_segment = recurse_revert_far_c_terminal_parent( c_terminal_segment_ );
				Basis c_term_basis;
				core::Size c_parent_n_basis_res;
				core::Size new_c_term_basis_res;
				core::Size new_n_term_basis_res;
				if ( last_reverted_parent_chimaera_->get_basis_pair().first.segment_id() == last_reverted_parent_chimaera_->get_n_terminal_parent()->get_segment_id() ) {
					//first basis is n_terminal
					c_parent_n_basis_res = last_reverted_parent_chimaera_->get_basis_pair().first.resnum();
					c_term_basis = last_reverted_parent_chimaera_->get_basis_pair().second;
				} else {
					//second basis is n_terminal
					c_parent_n_basis_res = last_reverted_parent_chimaera_->get_basis_pair().second.resnum();
					c_term_basis = last_reverted_parent_chimaera_->get_basis_pair().first;
				}
				if ( new_segment->get_basis_pair().first.segment_id() == new_segment->get_n_terminal_parent()->get_segment_id() ) {
					new_n_term_basis_res = new_segment->get_basis_pair().first.resnum();
					new_c_term_basis_res = new_segment->get_basis_pair().second.resnum();
				} else {
					new_n_term_basis_res = new_segment->get_basis_pair().second.resnum();
					new_c_term_basis_res = new_segment->get_basis_pair().first.resnum();
				}

				//core::Size n_basis_res = c_parent_n_basis_res - new_c_term_basis_res + new_n_term_basis_res;
				core::Size n_basis_res = c_parent_n_basis_res + new_n_term_basis_res - new_c_term_basis_res;
				Basis n_basis( c_terminal_segment_->get_segment_id(), n_basis_res );
				//last_chimaera_deleted_ = last_reverted_parent_chimaera_;
				//MAKE SURE THAT THIS ORDER IS CORRECT! This is from recurse_revert far c term--we deleted from the c term, so N term will be first and C term second as shown (N term was original)
				last_chimaera_deleted_ = std::make_pair( n_basis, c_term_basis );
			} else {
				if ( std::dynamic_pointer_cast< LigandSegment >( c_terminal_segment_ ) != nullptr ) {
					unchimerize_ligand_segment( std::dynamic_pointer_cast< LigandSegment >( c_terminal_segment_ ) );
				}
				//    TR << "c_terminal_segment id: " << c_terminal_segment_->get_segment_id() << ". reassigned segment id: " << c_terminal_segment_->get_n_terminal_parent()->get_segment_id() << ". is it chimaeric: " << c_terminal_segment_->get_n_terminal_parent()->is_chimaeric() << std::endl;
				last_chimaera_deleted_ = c_terminal_segment_->get_basis_pair();
				c_terminal_segment_->set_is_in_Assembly( false );
				c_terminal_segment_->get_c_terminal_parent()->set_is_in_Assembly( false );
				SmartSegment::link_to(c_terminal_segment_->get_n_terminal_parent()->get_n_terminal_neighbor(),c_terminal_segment_->get_n_terminal_parent());
				SmartSegment::link_to(c_terminal_segment_->get_c_terminal_parent(),c_terminal_segment_->get_c_terminal_parent()->get_c_terminal_neighbor());
				//Update the length
				length_ -= c_terminal_segment_->get_length();
				if ( n_terminal_segment_->get_segment_id() == c_terminal_segment_->get_segment_id() ) { // check if we are reverting to a single segment
					n_terminal_segment_ = c_terminal_segment_->get_n_terminal_parent();
				}
				c_terminal_segment_ = c_terminal_segment_->get_n_terminal_parent();
				c_terminal_segment_->set_is_in_Assembly( true );
				length_ += c_terminal_segment_->get_length();
			}
		}
		//  TR << "C-terminal segment is now " << c_terminal_segment_->get_segment_id();
	}
	TR.Debug << "delete successful!" << std::endl;
	last_change_ = 'D';
	last_change_was_n_terminal_ = n_terminus;
	return true;
}

bool
SmartAssembly::switch_segment(bool n_terminus){
	TR.Debug << "Switching segment." << std::endl;
	//TR.Debug << this->get_comprehensive_forward_assembly() << std::endl;
	if ( this->delete_segment(n_terminus) ) {
		BasisPair switch_chimaera = last_chimaera_deleted_;
		//TR.Debug << this->get_comprehensive_forward_assembly() << std::endl;
		if ( this->add_segment(n_terminus) ) {
			//TR << this->get_comprehensive_forward_assembly() << std::endl;
			last_change_ = 'S';
			last_change_was_n_terminal_ = n_terminus;
			last_chimaera_deleted_ = switch_chimaera;
			TR.Debug << "Switch successful!" << std::endl;
			return true;
		} else {
			TR.Debug << "Switch add failed!" << std::endl;
			this->undelete_segment();
			return false;
		}
	} else {
		TR.Debug << "Switch delete failed!" << std::endl;
		return false;
	}
}


void
SmartAssembly::unchimerize_ligand_segment( LigandSegmentOP ligseg ){
	if ( !ligseg->is_chimaeric() ) {
		TR << "WARNING: Attempted to unchimerize a non-chimaeric ligand segment!" << std::endl;
		return;
	}
	//We know its N-terminal and C-terminal parents and the basis pair used to create it
	BasisPair basis_pair = ligseg->get_basis_pair();
	Basis nterm_basis;
	Basis cterm_basis;
	//Problem with the basis pair--we don't know which one was N-terminal and which was C-terminal
	//I'll just make two separate Basis objects for N- and C-terminal
	if ( basis_pair.first.segment_id() == ligseg->get_n_terminal_parent()->get_segment_id() ) {
		nterm_basis = basis_pair.first;
		cterm_basis = basis_pair.second;
	} else {
		nterm_basis = basis_pair.second;
		cterm_basis = basis_pair.first;
	}
	bool n_terminal_ligand = ( std::dynamic_pointer_cast< LigandSegment >( ligseg->get_n_terminal_parent() ) != nullptr );
	bool c_terminal_ligand = ( std::dynamic_pointer_cast< LigandSegment >( ligseg->get_c_terminal_parent() ) != nullptr );
	runtime_assert( n_terminal_ligand || c_terminal_ligand ); //One of the parents must have been a ligand segment
	//We won't need to do anything to the parents since they still think they own/contact the ligand(s)
	//Start by iterating through owned ligands
	for ( std::pair< core::Size, LigandResidueOP > ligand: ligseg->get_owned_ligand_residues() ) {
		//Fix the owner
		//Only ligand segments can handle ligands
		if ( n_terminal_ligand && c_terminal_ligand ) {
			LigandSegmentOP n_ligseg = std::dynamic_pointer_cast< LigandSegment>( ligseg->get_n_terminal_parent() );
			LigandSegmentOP c_ligseg = std::dynamic_pointer_cast< LigandSegment>( ligseg->get_c_terminal_parent() );
			//if( std::find( n_ligseg->get_owned_ligand_residues().begin(), n_ligseg->get_owned_ligand_residues().end(), ligand  ) != n_ligseg->get_owned_ligand_residues().end() ){
			if ( n_ligseg->get_owned_ligand_residues().count( ligand.first ) != 0 ) {
				//The ligand is owned by the n-terminal segment
				ligand.second->set_owner_segment( n_ligseg );
				reset_chimaera_contacts( ligand.second, ligseg, nterm_basis, cterm_basis, true );
			} else {
				//The C-terminal segment owns it
				runtime_assert( c_ligseg->get_owned_ligand_residues().count( ligand.first ) != 0 );
				ligand.second->set_owner_segment( c_ligseg );
				reset_chimaera_contacts( ligand.second, ligseg, nterm_basis, cterm_basis, false );
			}
		} else if ( n_terminal_ligand ) {
			ligand.second->set_owner_segment( std::dynamic_pointer_cast< LigandSegment >( ligseg->get_n_terminal_parent() ) );
			reset_chimaera_contacts( ligand.second, ligseg, nterm_basis, cterm_basis, true );
		} else { //The C-terminal segment must be a ligand segment
			ligand.second->set_owner_segment( std::dynamic_pointer_cast< LigandSegment >( ligseg->get_c_terminal_parent() ) );
			reset_chimaera_contacts( ligand.second, ligseg, nterm_basis, cterm_basis, false );
		}
	}
	//Then iterate through unowned ligands
	for ( core::Size ligand_id : ligseg->get_ligand_residues() ) {
		//No need to change the owner segments
		//But we need to figure out which segment this contact came from
		bool nterm = true;
		if ( n_terminal_ligand && std::dynamic_pointer_cast< LigandSegment >( ligseg->get_n_terminal_parent() )->get_ligand_residues().count( ligand_id ) != 0 ) {
			nterm = true;
		} else {
			nterm = false;
		}
		//If the ligand is not in local_ligands_, that's bad
		if ( local_ligands_.count( ligand_id ) == 0 ) {
			utility_exit_with_message( "Bad ligand ID!" );
		}
		reset_chimaera_contacts( local_ligands_.at( ligand_id ), ligseg, nterm_basis, cterm_basis, nterm );
	}
}


void
SmartAssembly::revert(){
	TR.Debug << "Reverting Last Move" << std::endl;
	if ( last_change_ == 'A' ) {
		//  this->unadd_segment();
		this->delete_segment( last_change_was_n_terminal_ );
	} else if ( last_change_ == 'D' ) {
		this->undelete_segment();
	} else if ( last_change_ == 'S' ) {
		this->unswitch_segment();
	} else if ( last_change_ == 'L' ) {
		this->unsample_ligand();
	} else {
		utility_exit_with_message("Last move is not of a revertable type!");
	}
	last_change_='R';
	// TR << "Does the assembly have chain breaks?: " << !this->is_continuous() << std::endl;
	// TR <<"N-terminal segment is: " << this->get_n_terminal_segment()->get_segment_id() << std::endl;
	// TR <<"C-terminal segment is: " << this->get_c_terminal_segment()->get_segment_id() << std::endl;
	// TR << "Current Forward Assembly is: " << this->get_forward_assembly() << std::endl;
	// TR << "Current Reverse Assembly is: " << this->get_reverse_assembly() << std::endl;
	// TR << "Assembly size is " << this->get_size() << " segments." << std::endl;
	// TR << "Assembly length is " << this->get_length() << " residues." << std::endl;
	// TR << "END REVERT" << std::endl;
}



void
SmartAssembly::reset_chimaera_contacts( LigandResidueOP ligand, LigandSegmentOP ligand_chimaera, Basis nterm_basis, Basis cterm_basis, bool nterm){
	//I think if we call this each time we unchimerize we'll be good to go
	//This part is handled outside this function
	/*
	if( nterm ){
	ligand->set_owner_segment( ligand_chimaera->get_n_terminal_parent() );
	}
	else{
	ligand->set_owner_segment( ligand_chimaera->get_c_terminal_parent() );
	}
	*/
	for ( LigandContactOP contact: ligand->get_nonconst_current_contacts() ) {
		//TR << "Old residue number: " << contact->residue_number << std::endl;
		//

		if ( contact->segment_id == 0 ) {
			TR.Debug << "Internal contacts do not change" << std::endl;
			continue;
		}
		runtime_assert( contact->residue_number != 0 );
		//Only need to deal with contacts of this chimaeric segment
		if ( contact->segment_id == ligand_chimaera->get_segment_id() ) {
			//If the residue came from the N-terminal segment, no need to change the resnum
			//if( contact->residue_number <= nterm_basis.resnum() ){
			if ( nterm ) {
				// TR << "Contact came from N terminal segment" << std::endl;
				contact->segment_id = nterm_basis.segment_id();
				//    TR << "Residue number unchanged" << std::endl;
			} else {
				//TR << "Contact came from C terminal segment" << std::endl;
				contact->segment_id = cterm_basis.segment_id();
				// TR << "Nterm basis resnum: " << nterm_basis.resnum() << " Cterm basis resnum: " << cterm_basis.resnum() << std::endl;
				//NEW = OLD + N - C, so OLD = NEW + C - N
				runtime_assert( contact->residue_number + cterm_basis.resnum() > nterm_basis.resnum() );
				if ( contact->residue_number + cterm_basis.resnum() <= nterm_basis.resnum() ) {
					utility_exit_with_message( "ERROR IN CHIMAERA RESET! Check your math because it's wrong." );
				}
				contact->residue_number = contact->residue_number + cterm_basis.resnum()  - nterm_basis.resnum();
				//TR << "New residue number: " << contact->residue_number << std::endl;
			}
		}
	}
}

bool
SmartAssembly::undelete_segment(){
	//BasisPair basis_pair = last_chimaera_deleted_->get_basis_pair();
	//return add_segment( last_change_was_n_terminal_, basis_pair.second.segment_id(), basis_pair.first.resnum(), basis_pair.second.resnum() );
	return add_segment( last_change_was_n_terminal_, last_chimaera_deleted_.second.segment_id(), last_chimaera_deleted_.first.resnum(), last_chimaera_deleted_.second.resnum() );
	/*
	if(last_change_was_n_terminal_){
	TR << "Undeleting segment from N terminus" << std::endl;
	size_--;
	length_-= n_terminal_segment_->get_length();
	//n_terminal_segment_->set_is_in_Assembly( false );
	SmartSegment::link_to(last_chimaera_,n_terminal_segment_->get_c_terminal_neighbor());
	SmartSegment::link_to(last_chimaera_->get_n_terminal_parent()->get_n_terminal_neighbor(),last_chimaera_);
	n_terminal_segment_ = last_chimaera_;
	n_terminal_segment_->get_n_terminal_parent()->set_is_in_Assembly( true );
	while(n_terminal_segment_->get_n_terminal_neighbor()){
	n_terminal_segment_->set_is_in_Assembly(true);
	size_++;
	length_+= n_terminal_segment_->get_length();
	n_terminal_segment_ = n_terminal_segment_->get_n_terminal_neighbor();
	}
	n_terminal_segment_->set_is_in_Assembly(true);
	size_++;
	length_+= n_terminal_segment_->get_length();
	} else {
	size_--;
	length_-= c_terminal_segment_->get_length();
	//c_terminal_segment_->set_is_in_Assembly( false );
	TR << "Undeleting segment from C terminus" << std::endl;
	SmartSegment::link_to(c_terminal_segment_->get_n_terminal_neighbor(), last_chimaera_ );
	runtime_assert( last_chimaera_->get_n_terminal_parent() == c_terminal_segment_ );

	SmartSegment::link_to(last_chimaera_,last_chimaera_->get_c_terminal_parent()->get_c_terminal_neighbor()); //This should be nullptr, right?

	c_terminal_segment_ = last_chimaera_;
	c_terminal_segment_->get_c_terminal_parent()->set_is_in_Assembly( true );
	while(c_terminal_segment_->get_c_terminal_neighbor()){
	c_terminal_segment_->set_is_in_Assembly(true);
	size_++;
	length_+= c_terminal_segment_->get_length();
	c_terminal_segment_ = c_terminal_segment_->get_c_terminal_neighbor();
	}
	c_terminal_segment_->set_is_in_Assembly(true);
	size_++;
	length_+= c_terminal_segment_->get_length();
	}
	*/
}
void
SmartAssembly::unswitch_segment(){
	// this->unadd_segment();
	BasisPair switched_segment = last_chimaera_deleted_;
	this->delete_segment( last_change_was_n_terminal_ );
	// TR << "Does the assembly have chain breaks?: " << !this->is_continuous() << std::endl;
	BasisPair temp = last_chimaera_deleted_;
	last_chimaera_deleted_ = switched_segment;
	if ( !this->undelete_segment() ) {
		utility_exit_with_message( "Unswitch failed!" );
		this->add_segment( last_change_was_n_terminal_, temp.second.segment_id(), temp.first.resnum(), temp.second.resnum() );
	}
	// TR << "Does the assembly have chain breaks?: " << !this->is_continuous() << std::endl;
}

void
SmartAssembly::dump_side_chains(){
}

bool
SmartAssembly::is_continuous() const {
	SmartSegmentOP current_segment = this->get_n_terminal_segment();
	while ( current_segment ) {
		for ( core::Size res = 1; res <= current_segment->get_length(); ++res ) {
			if ( res == current_segment->get_length() && current_segment->get_c_terminal_neighbor() ) {
				if ( current_segment->get_residue( res )->get_atom( 2 ).xyz().distance(
						current_segment->get_c_terminal_neighbor()->get_residue( 1 )->get_atom( 2 ).xyz() ) > 5.0// typical caplha distance is ~4, but have seen as high as 4.3
						) {
					return false;
				}
			} else {
				if ( res != current_segment->get_length() ) {
					if ( current_segment->get_residue( res )->get_atom( 2 ).xyz().distance(
							current_segment->get_residue( res + 1 )->get_atom( 2 ).xyz() ) > 5.0// typical caplha distance is ~4, but have seen as high as 4.3
							) {
						return false;
					}
				}
			}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return true;
}

void
SmartAssembly::pick_random_starting_segment(){
	core::Size starting_seg_index = std::ceil( ( segment_vector_->size() + pdb_segments_.size() )*numeric::random::rg().uniform() );
	//We'll actually need to copy over all the associated segments as well (if they're not already in the map)
	for ( std::pair< core::Size, SmartSegmentOP > pdbseg: pdb_segments_ ) {
		if ( local_segments_.count( pdbseg.first ) == 0 ) {
			TR <<"WARNING: PDB segment " << pdbseg.first << " was not added upon loading and may not be treated appropriately!" << std::endl;
			add_segment_and_neighbors_to_local_segments( pdbseg.second );
		}
	}
	//Do this once any PDB segments are in local_segments
	SmartSegmentCOP current_seg;
	if ( starting_seg_index <= segment_vector_->size() ) {
		current_seg = SmartSegment::get_n_most_segment( segment_vector_->at( starting_seg_index ), false );
	} else {
		current_seg = SmartSegment::get_n_most_segment( pdb_segments_.at( starting_seg_index ), false );
	}
	add_segment_and_neighbors_to_local_segments( current_seg );
	this->set_starting_segment( local_segments_.at( starting_seg_index ), this->get_start_node_vital_segments() );
}

core::pose::Pose // original written by Tim Jacobs
SmartAssembly::to_pose(std::string residue_type_set) {
	TR << "Generating Pose from Assembly." << std::endl;
	core::chemical::ResidueTypeSetCOP res_type_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( residue_type_set );


	//really, all the ligands should be output at the end of the pose
	utility::vector1< LigandResidueCOP > ligands_to_output;

	core::Size nres = 0;
	std::string sequence = "";
	std::set<core::Size> all_vital_residues;
	SmartSegmentCOP current_segment = this->get_n_terminal_segment();
	while ( current_segment ) {
		//first, add any vital residues to our master vital residues set
		if ( current_segment->get_vital_residues().size()>0 ) {
			for ( core::Size vital_residue : current_segment->get_vital_residues() ) {
				TR.Debug << "Found vital residue: Segment " << current_segment->get_segment_id() << " residue " << vital_residue << ", residue " << vital_residue + nres << " of pose" << std::endl;
				all_vital_residues.insert(vital_residue + nres);
			}
		}
		nres += current_segment->get_length();
		TR << "Segment " << current_segment->get_segment_id() << " is " << current_segment->get_length() << " residues long. " << std::endl;
		for ( SmartSewingResidueOP current_residue : current_segment->get_const_residue_vector() ) {
			sequence += res_type_set->name_map(current_residue->get_amino_acid_type()).name1();
			//for(core::Size atomcounter = 1; atomcounter <= current_residue->get_atom_vector().size(); atomcounter++){
			// TR << "Atom " << atomcounter << " in segment " << current_segment->get_segment_id() << " has coordinates " << current_residue->get_atom(atomcounter).xyz().at(0) << " " << current_residue->get_atom(atomcounter).xyz().at(1) << " " << current_residue->get_atom(atomcounter).xyz().at(2) << std::endl;
			//}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	TR << "Assembly is " << nres << " residues long." << std::endl;
	TR << "Assembly sequence is " << sequence.size() << "residues long." << std::endl;
	TR << "Assembly sequence is " << sequence << std::endl;

	core::pose::Pose model_pose;
	//core::pose::Pose temp_pose_to_dump;
	core::pose::make_pose_from_sequence(model_pose, sequence, *res_type_set);
	//core::pose::make_pose_from_sequence(temp_pose_to_dump, "AAA", *res_type_set);
	TR << "Made model pose" << std::endl;
	// we know it works up to this point
	core::Size counter = 1;
	current_segment=this->get_n_terminal_segment();
	while ( current_segment ) {
		//for ( SmartSewingResidueOP current_residue : current_segment->get_const_residue_vector() ) {
		for ( core::Size current_res_i = 1; current_res_i <= current_segment->get_const_residue_vector().size(); ++current_res_i ) {
			SmartSewingResidueOP current_residue = current_segment->get_const_residue_vector()[ current_res_i ];
			core::conformation::ResidueOP template_residue = model_pose.residue(counter).clone();
			//TR << "Template Residue Coordinates: " << std::endl;
			for ( core::Size template_atom = 1; template_atom <= template_residue->atoms().size(); template_atom++ ) { // should be 4
				std::string atom_name = template_residue->atom_name( template_atom );
				if ( ! ( res_type_set->name_map( current_residue->get_full_type_name() ).has( atom_name ) ) ) {
					continue;
				}
				//Also don't copy coordianates for C-terimal O
				////////??TODO!!!!//////////
				core::Size current_index = res_type_set->name_map( current_residue->get_full_type_name() ).atom_index( atom_name );
				if ( current_index > current_residue->get_atom_vector().size() ) {
					continue;
				}
				template_residue->atom( template_atom ).xyz(current_residue->get_atom( current_index ).xyz() );
				//TR << "Atom " << current_atom << " has coordinates " << template_residue->atom(current_atom).xyz().at(0) << " " << template_residue->atom(current_atom).xyz().at(1) << " " << template_residue->atom(current_atom).xyz().at(2) << std::endl;
			}
			core::conformation::ResidueOP new_residue;
			core::chemical::ResidueType const & base_type = model_pose.residue(counter).type();
			if ( counter != 1 && base_type.is_lower_terminus() ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_removed(base_type, core::chemical::LOWER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else if ( counter != model_pose.total_residue() && base_type.is_upper_terminus() ) {
				core::chemical::ResidueType const & res_type = res_type_set->get_residue_type_with_variant_removed(base_type, core::chemical::UPPER_TERMINUS_VARIANT);
				new_residue = core::conformation::ResidueFactory::create_residue(res_type);
			} else {
				new_residue = core::conformation::ResidueFactory::create_residue(base_type);
			}
			/*
			if( counter == nres ){
			temp_pose_to_dump.replace_residue( 3, *new_residue, false );
			temp_pose_to_dump.dump_pdb( "test_before_orient.pdb" );
			}
			*/
			new_residue->orient_onto_residue(*template_residue);
			/*if( counter == nres ){
			temp_pose_to_dump.replace_residue( 3, *new_residue, false );
			temp_pose_to_dump.dump_pdb( "test_after_orient.pdb" );
			}
			*/
			//TR << "New Residue Coordinates: " << std::endl;



			//To avoid misplacement of N/C terminal patch atoms, at this point we should ONLY copy atom coordinates from atoms present in the SSR
			//for(core::Size template_atom = 1; template_atom <= template_residue->atoms().size(); template_atom++) { // should be 4
			for ( core::Size current_atom = 1; current_atom <= current_residue->get_atom_vector().size(); ++current_atom ) {
				//TR << "Atom " << current_atom << " has coordinates " << new_residue->atom(current_atom).xyz().at(0) << " " << new_residue->atom(current_atom).xyz().at(1) << " " << new_residue->atom(current_atom).xyz().at(2) << std::endl;

				//std::string atom_name = template_residue->atom_name( template_atom );

				std::string atom_name = res_type_set->name_map( current_residue->get_full_type_name() ).atom_name( current_atom  );
				if ( ! ( new_residue->has( atom_name ) )  ) {
					continue;
				}
				if ( ! ( template_residue->has( atom_name ) ) ) {
					continue;
				}
				//Keep new position of c terminal oxygen
				std::string atom_type_name = res_type_set->name_map( current_residue->get_full_type_name() ).atom_type( current_atom ).atom_type_name();
				if ( counter == nres && atom_type_name == "OCbb" ) { //OCbb
					continue;
				}
				core::Size new_index = new_residue->atom_index( atom_name );
				new_residue->atom(new_index).xyz( current_residue->get_atom( current_atom ).xyz());
				//TR << "Atom " << current_atom << " has coordinates " << template_residue->atom(current_atom).xyz().at(0) << " " << template_residue->atom(current_atom).xyz().at(1) << " " << template_residue->atom(current_atom).xyz().at(2) << std::endl;
			}

			/*
			if ( current_residue->get_chis().size() == new_residue->nchi() ) {
			new_residue->set_all_chi(current_residue->get_chis());
			}
			if ( current_residue->get_chis().size() == template_residue->nchi() ) {
			template_residue->set_all_chi(current_residue->get_chis());
			}
			*/
			// Fix backbone oxygens and hydrogens
			model_pose.replace_residue(counter, *new_residue, false);
			++counter;
		}
		//Now output all of the owned ligands
		//if( current_segment->get_name() == "LigandSegment" ){
		if ( std::dynamic_pointer_cast< LigandSegment const >( current_segment ) != nullptr ) {
			LigandSegmentCOP ligseg = std::dynamic_pointer_cast< LigandSegment const >( current_segment );
			for ( std::pair< core::Size, LigandResidueOP > ligand: ligseg->get_const_owned_ligand_residues() ) {
				ligands_to_output.push_back( LigandResidueCOP( ligand.second ) );
			}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	//Now that we've gone through all the segments, add in all the ligands
	//core::chemical::ResidueTypeFinderOP ligfinder( new ResidueTypeFinder( res_type_set ) );
	for ( LigandResidueCOP ligand: ligands_to_output ) {
		//These will be a little different since
		//a) we're appending by a jump and
		//b) these residues were in no way already associated with the pose
		//Use the residue type constructor
		core::chemical::ResidueTypeCOP ligand_type = res_type_set->name_mapOP( ligand->get_full_type_name() );
		core::conformation::ResidueOP ligand_res( new core::conformation::Residue( ligand_type, true ) );
		for ( core::Size atom_index = 1; atom_index <= ligand->get_const_atom_vector().size(); ++atom_index ) {
			//ligand_res->atom( atom_index ).xyz( ligand->get_atom( atom_index ).xyz() );
			ligand_res->atom( atom_index ).xyz( ligand->get_const_atom_vector().at( atom_index ).xyz() );
			//These are ligands, so we don't have to worry about terminal variants
		}
		//Find the residue's first contact
		if ( ligand->get_current_contacts().size() == 0 ) {
			utility_exit_with_message( "Ligand found with no contacts! Something probably went wrong with cloning a ligand." );
		}
		model_pose.append_residue_by_jump( *ligand_res, model_pose.total_residue() );
	}
	runtime_assert(model_pose.total_residue() == nres + ligands_to_output.size() );

	//finally, add the partner
	if ( partner_ && output_partner_ ) {
		core::pose::Pose partner_pose_copy = *partner_->clone();
		if ( core::chemical::string_from_type_set_mode( partner_->conformation().residue_typeset_mode() ) != residue_type_set  ) {
			core::util::switch_to_residue_type_set( partner_pose_copy, residue_type_set, true );
		}
		core::pose::append_pose_to_pose(model_pose, partner_pose_copy, true);
	}
	//Now add the resfile to the remarks
	//std::string vital_residues_string = "VITAL_RESIDUES:";
	core::pose::PDBInfoOP new_info = core::pose::PDBInfoOP(new core::pose::PDBInfo(model_pose,true));
	model_pose.pdb_info(new_info);
	for ( core::Size vital_residue : all_vital_residues ) {
		TR.Debug << "Adding tag for residue " << vital_residue << std::endl;
		model_pose.pdb_info()->add_reslabel(vital_residue,"VITAL");
		//vital_residues_string.append(" ");
		//vital_residues_string.append(std::to_string(int(vital_residue)));
	}
	//core::io::RemarkInfo new_vital;
	//new_vital.value = vital_residues_string;
	//new_vital.num=0;
	//core::io::RemarksCOP new_remarks = core::io::RemarksCOP(new core::io::Remarks());
	//model_pose.pdb_info()->remarks(*new_remarks);
	//model_pose.pdb_info()->remarks().push_back(new_vital);
	//TR << model_pose.pdb_info()->remarks().at(0).value << std::endl;

	return model_pose;
}



bool
SmartAssembly::chimerize(BasisPair basis_pair, bool n_terminus ){

	TR.Debug << " Chimerizing segments " << basis_pair.first.segment_id() << " " << basis_pair.first.resnum() << " and " << basis_pair.second.segment_id() << " " << basis_pair.second.resnum() << "." << std::endl;
	//First basis should have been in assembly to start with
	SmartSegmentOP new_chimaera;
	//Add second if it's not already in there
	if ( local_segments_.count( basis_pair.second.segment_id() ) == 0 ) {
		SmartSegmentCOP current_seg;
		if ( pdb_segments_.count( basis_pair.second.segment_id() ) == 0 ) {
			current_seg = SmartSegment::get_n_most_segment( segment_vector_->at( basis_pair.second.segment_id() ), false );
		} else {
			current_seg = SmartSegment::get_n_most_segment( pdb_segments_.at( basis_pair.second.segment_id() ), false );
		}
		add_segment_and_neighbors_to_local_segments( current_seg );
	}
	if ( !local_segments_.at(basis_pair.first.segment_id())->is_in_Assembly() ) {
		utility_exit_with_message("First segment of basispair is not in Assembly!");
	}

	//Transform the segments
	if ( !transform_segments(basis_pair) ) {
		return false;
	}
	bool ligand_seg_1 = ( std::dynamic_pointer_cast< LigandSegment >( local_segments_.at( basis_pair.first.segment_id() ) ) != nullptr );
	bool ligand_seg_2 = ( std::dynamic_pointer_cast< LigandSegment >( local_segments_.at( basis_pair.second.segment_id() ) ) != nullptr );
	std::set< core::Size > ligseg1_contacts;
	std::set< core::Size > ligseg2_contacts;
	if ( ligand_seg_1 || ligand_seg_2 ) {
		LigandSegmentOP ligseg(new LigandSegment );
		if ( ligand_seg_1 ) {
			LigandSegmentCOP src = std::dynamic_pointer_cast< const LigandSegment >( local_segments_.at( basis_pair.first.segment_id() ) );
			for ( std::pair< core::Size, LigandResidueOP > owned_ligand: src->get_const_owned_ligand_residues() ) {
				ligseg->attach_ligand( local_ligands_.at( owned_ligand.second->get_ligand_id() ), true );
			}
			for ( core::Size ligand_id: src->get_ligand_residues() ) {
				ligseg->attach_unowned_ligand( local_ligands_.at( ligand_id ) );
			}
			for ( core::Size contact: src->get_ligand_contact_indices() ) {
				ligseg1_contacts.insert( contact );
			}
		}
		if ( ligand_seg_2 ) {
			LigandSegmentCOP src = std::dynamic_pointer_cast< const LigandSegment >( local_segments_.at( basis_pair.second.segment_id() ) );
			for ( std::pair< core::Size, LigandResidueOP > owned_ligand: src->get_const_owned_ligand_residues() ) {
				ligseg->attach_ligand( local_ligands_.at( owned_ligand.second->get_ligand_id() ), true );
			}
			for ( core::Size ligand_id: src->get_ligand_residues() ) {
				ligseg->attach_unowned_ligand( local_ligands_.at( ligand_id ) );
			}
			for ( core::Size contact: src->get_ligand_contact_indices() ) {
				ligseg2_contacts.insert( contact );
			}
		}
		new_chimaera = SmartSegmentOP( ligseg );
	} else {
		new_chimaera = SmartSegmentOP(new SmartSegment());
	}
	new_chimaera->set_basis_pair( std::make_pair(  basis_pair.first, basis_pair.second ) );
	new_chimaera->set_const_reference_segment( new_chimaera ); //Just in case they check for it
	//utility::vector1< std::pair< segID, LigandResidueOP > >
	new_chimaera->set_segment_id(2147483647*numeric::random::rg().uniform()+segment_vector_->size() + pdb_segments_.size() );
	local_segments_[new_chimaera->get_segment_id()] = new_chimaera;
	new_chimaera->set_is_in_Assembly(true);
	new_chimaera->set_chimaeric_status(true);
	//
	utility::vector1< SmartSewingResidueOP > residues;
	//Figure out exactly how the residues line up one to one
	//N term always starts at 1
	//Goes until it reaches the basis residue
	//C term will start at the basis residue and go to the end

	core::Size resnum = 1;
	if ( n_terminus ) { //are we adding to the N terminus of our assembly?
		//Ligseg1_contacts is C-term, Ligseg2_contacts is N-term
		core::Size nterm_resnum = 1;
		int cterm_resnum = basis_pair.first.resnum() - basis_pair.second.resnum() + 1;
		//TR << "Adding to N terminus" << std::endl;
		new_chimaera->set_c_terminal_parent(local_segments_.at(basis_pair.first.segment_id()));
		//TR << "c_terminal parent: " << new_chimaera->get_c_terminal_parent()->get_segment_id() << " basis id given: " << basis_pair.first.segment_id() << std::endl;
		new_chimaera->set_n_terminal_parent(local_segments_.at(basis_pair.second.segment_id()));
		//TR << "c_terminal parent: " << new_chimaera->get_n_terminal_parent()->get_segment_id() << " basis id given: " << basis_pair.second.segment_id() << std::endl;
		while ( resnum <= basis_pair.second.resnum() ) {
			//   TR << "Adding residue " << resnum << " of " << basis_pair.second.resnum() << std::endl;
			//If this residue is vital in cterm
			if ( cterm_resnum > 0 && core::Size( cterm_resnum ) <= new_chimaera->get_c_terminal_parent()->get_length() && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) {
				TR.Debug << "Chimerizing vital residue " << cterm_resnum << " from segment " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
				TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
				//If this residue is vital in BOTH parents, we can't chimerize because one or the other would have to be deleted
				if ( nterm_resnum <= new_chimaera->get_n_terminal_parent()->get_length() && new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
					return false;
				}
				//Otherwise take the vital residue from the C-terminal parent
				residues.push_back( new_chimaera->get_c_terminal_parent()->get_residue( core::Size( cterm_resnum ) ) );
				new_chimaera->add_vital_residue( nterm_resnum );
				if ( ligseg1_contacts.find( core::Size( cterm_resnum ) ) != ligseg1_contacts.end() ) {
					//This also sets it as a vital residue
					std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
				}
			} else {
				runtime_assert( new_chimaera->get_n_terminal_parent()->get_length() >= resnum );
				residues.push_back(new_chimaera->get_n_terminal_parent()->get_residue(resnum));
				if ( new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
					TR.Debug << "Chimerizing vital residue " << nterm_resnum << " from segment " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
					TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
					new_chimaera->add_vital_residue( nterm_resnum );
				}
			}
			//In this scenario, ligseg2_contacts will be checked first
			if ( ligseg2_contacts.find( nterm_resnum ) != ligseg2_contacts.end() ) {
				std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
			}
			resnum++;
			nterm_resnum++;
			cterm_resnum++;
		}
		//TR << "Switching to C-terminal segment" << std::endl;
		resnum = basis_pair.first.resnum()+1;
		runtime_assert( resnum == core::Size( cterm_resnum ) );
		while ( resnum <= new_chimaera->get_c_terminal_parent()->get_length() ) {
			//   TR << "Adding residue " << resnum << " of " << new_chimaera->get_c_terminal_parent()->get_length() << std::endl;
			if ( nterm_resnum <= new_chimaera->get_n_terminal_parent()->get_length() && new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
				TR.Debug << "Chimerizing vital residue " << nterm_resnum << " from segment " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
				TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
				if ( cterm_resnum > 0 && core::Size( cterm_resnum ) <= new_chimaera->get_c_terminal_parent()->get_length() && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) {
					return false;
				}
				residues.push_back(new_chimaera->get_n_terminal_parent()->get_residue(nterm_resnum));
				new_chimaera->add_vital_residue( nterm_resnum );
			} else {
				residues.push_back(new_chimaera->get_c_terminal_parent()->get_residue(resnum));
				if ( cterm_resnum > 0 && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) {
					new_chimaera->add_vital_residue( nterm_resnum );
					TR.Debug << "Chimerizing vital residue " << cterm_resnum << " from segment " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
					TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
				}
			}
			if ( cterm_resnum > 0 && ligseg1_contacts.find( core::Size( cterm_resnum ) ) != ligseg1_contacts.end() ) {
				std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
			}
			resnum++;
			++nterm_resnum;
			++cterm_resnum;
		} //End iterate over residues





		//BEGIN
		//Now fix ligand contacts
		//If the new chimaera is a ligand segment, fix up the contacts for its ligands
		//SEG1 is C TERMINAL, SEG2 is N TERMINAL

		if ( std::dynamic_pointer_cast< LigandSegment >( new_chimaera ) != nullptr ) {
			for ( std::pair< core::Size, LigandResidueOP > ligand: std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->get_owned_ligand_residues() ) {
				//For these owned ligands, we'll also need to reset their owners to the chimaera instead of the old parent
				ligand.second->set_owner_segment(std::dynamic_pointer_cast< LigandSegment >( new_chimaera ) );
				for ( LigandContactOP contact: ligand.second->get_nonconst_current_contacts() ) {
					if ( contact->segment_id == 0 ) {
						TR.Debug << "Internal contacts do not change" << std::endl;
						continue;
					}
					//TR.Debug << "Old residue number: " << contact->residue_number << std::endl;
					runtime_assert( contact->residue_number != 0 );
					//basis_pair.first is C-terminal, so resnums will need to be fixed
					if ( contact->segment_id == basis_pair.first.segment_id() ) {
						contact->segment_id = new_chimaera->get_segment_id();
						//NEW RESNUM = OLD RESNUM - CTERM BASIS + NTERM BASIS
						contact->residue_number = contact->residue_number + basis_pair.second.resnum() - basis_pair.first.resnum();
						//TR.Debug << "New residue number: " << contact->residue_number << std::endl;
					} else if ( contact->segment_id == basis_pair.second.segment_id() ) {
						//basis_pair.second is N-terminal, so resnums will match up
						contact->segment_id = new_chimaera->get_segment_id();
						TR.Debug << "New residue number: " << contact->residue_number << std::endl;
					}
				} //end for contact
			}//end for owned ligand



			//Unowned ligands. Note that changes are only made if the contacts still have the old segment IDs, so it won't overwrite any changes already made.
			for ( core::Size ligand_id: std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->get_ligand_residues() ) {
				for ( LigandContactOP contact: local_ligands_.at( ligand_id )->get_nonconst_current_contacts() ) {
					if ( contact->segment_id == 0 ) {
						TR.Debug << "Internal contacts and partner contacts do not change" << std::endl;
						continue;
					}
					//TR.Debug << "Old residue number: " << contact->residue_number << std::endl;
					runtime_assert( contact->residue_number != 0 );
					//basis_pair.first is the new C-terminal part of this chimaera so resnum must be fixed
					//TR << "C terminal basis: " << basis_pair.first.resnum() << " N terminal basis: " << basis_pair.second.resnum() << std::endl;
					//This will only run if the contact hasn't already been fixed
					if ( contact->segment_id == basis_pair.first.segment_id() ) {
						//TR << "Contact in C-terminal part of chimaera. Fixing residue numbers." << std::endl;
						contact->segment_id = new_chimaera->get_segment_id();
						//Subtract C terminal basis, add N terminal basis to get to chimaeric resnum
						//NEW RESNUM = OLD RESNUM - CTERM BASIS + NTERM BASIS
						contact->residue_number = contact->residue_number + basis_pair.second.resnum()  - basis_pair.first.resnum();
						//TR.Debug << "New residue number: " << contact->residue_number << std::endl;
					} else if ( contact->segment_id == basis_pair.second.segment_id() ) {
						//basis_pair.second is N-terminal, so resnums will match up
						contact->segment_id = new_chimaera->get_segment_id();
						//TR.Debug << "New residue number: " << contact->residue_number << std::endl;
						//TR << "Contact in N-terminal part of chimaera. Resnum unchanged." << std::endl;
					}
				}//end for contact
			}//end for unowned ligand
		} //end if

		//END
	} else { //END ADDING TO N TERMINUS
		core::Size nterm_resnum = 1; //This is the residue number from the segment that is currently in the assembly. It also tracks the current residue number in the chimaera.
		int cterm_resnum = basis_pair.second.resnum() - basis_pair.first.resnum() + 1; //This is the residue number from the new segment that is not yet in the assembly
		//TR << "Adding to C terminus" << std::endl;
		new_chimaera->set_n_terminal_parent( local_segments_.at(basis_pair.first.segment_id()));
		//TR << "n_terminal parent: " << new_chimaera->get_n_terminal_parent()->get_segment_id() << " basis id given: " << basis_pair.first.segment_id() << std::endl;
		new_chimaera->set_c_terminal_parent( local_segments_.at(basis_pair.second.segment_id()));
		//TR << "c_terminal parent: " << new_chimaera->get_c_terminal_parent()->get_segment_id() << " basis id given: " << basis_pair.second.segment_id() << std::endl;

		while ( resnum <= basis_pair.first.resnum() ) { //residues from N-terminal parent
			//   TR << "Adding residue " << resnum << " of " << basis_pair.first.resnum() << std::endl;

			//If this is a vital residue in the C-terminal segment
			if ( cterm_resnum > 0 && core::Size( cterm_resnum ) <= new_chimaera->get_c_terminal_parent()->get_length() && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) { //
				//If it's a vital residue in both segments, then we can't make this chimaera
				TR << "WARNING: This residue is vital in the other segment. Adding vital residue." << std::endl;
				TR.Debug << "Chimerizing vital residue " << cterm_resnum << " from segment " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
				TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
				if ( nterm_resnum <= new_chimaera->get_n_terminal_parent()->get_length() && new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
					return false;
				}
				//Otherwise, take this residue from the C-terminal segment instead of the N-terminal segment like we normally would
				residues.push_back( new_chimaera->get_c_terminal_parent()->get_residue( core::Size( cterm_resnum ) ) );
				new_chimaera->add_vital_residue( nterm_resnum );
				//Go ahead and take care of the ligseg2 stuff here
				if ( ligseg2_contacts.find( core::Size( cterm_resnum ) ) != ligseg2_contacts.end() ) {
					// TR << "Adding ligand contact from seg2" << std::endl;
					std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
				}
			} else {
				residues.push_back(new_chimaera->get_n_terminal_parent()->get_residue( resnum ));
				if ( new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
					//TR << "Adding vital residue!" << std::endl;
					TR.Debug << "Chimerizing vital residue " << nterm_resnum << " from segment " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
					TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
					new_chimaera->add_vital_residue( nterm_resnum );
					//TR << "Vital residue added!" << std::endl;
				}
			}
			if ( ligseg1_contacts.find( nterm_resnum ) != ligseg1_contacts.end() ) {
				std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
			}
			++resnum;
			++nterm_resnum;
			++cterm_resnum;
		}
		//TR << "Switching to C-terminal segment" << std::endl;
		resnum = basis_pair.second.resnum()+1; //Start immediately after the basis residue
		runtime_assert( resnum == core::Size( cterm_resnum ) ); //Make sure we didn't skip or repeat anything
		while ( resnum <= new_chimaera->get_c_terminal_parent()->get_length() ) {
			//   TR << "Adding residue " << resnum << " of " << new_chimaera->get_n_terminal_parent()->get_length() << std::endl;
			if ( nterm_resnum <= new_chimaera->get_n_terminal_parent()->get_length() && new_chimaera->get_n_terminal_parent()->residue_is_vital( nterm_resnum ) ) {
				TR.Debug << "Chimerizing vital residue " << nterm_resnum << " from segment " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
				TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
				if ( cterm_resnum > 0 && core::Size( cterm_resnum ) <= new_chimaera->get_c_terminal_parent()->get_length() && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) {
					return false;
				}
				residues.push_back(new_chimaera->get_n_terminal_parent()->get_residue(nterm_resnum));
				new_chimaera->add_vital_residue( nterm_resnum );
				if ( ligseg1_contacts.find( nterm_resnum ) != ligseg1_contacts.end() ) {
					//TR << "This is a ligand contact on the N-terminal segment seg1" << std::endl;
					std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
				}
			} else {
				residues.push_back(new_chimaera->get_c_terminal_parent()->get_residue(resnum));
				if ( cterm_resnum > 0 && new_chimaera->get_c_terminal_parent()->residue_is_vital( core::Size( cterm_resnum ) ) ) {
					TR.Debug << "Chimerizing vital residue " << cterm_resnum << " from segment " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
					TR.Debug << "Became residue " << nterm_resnum << " of chimaera" << new_chimaera->get_segment_id() << std::endl;
					//WAS THIS IT?
					new_chimaera->add_vital_residue( nterm_resnum );
				}
			}
			//residues.push_back(new_chimaera->get_c_terminal_parent()->get_residue(resnum));
			if ( cterm_resnum > 0 && ligseg2_contacts.find( core::Size( cterm_resnum ) ) != ligseg2_contacts.end() ) {
				std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->add_ligand_contact( nterm_resnum );
			}
			++resnum;
			++nterm_resnum;
			++cterm_resnum;
		} //end iterate over residues

		//If the new chimaera is a ligand segment, fix up the contacts for its ligands
		//TEMP
		//SEG1 is N TERMINAL, SEG2 is C TERMINAL
		if ( std::dynamic_pointer_cast< LigandSegment >( new_chimaera ) != nullptr ) {
			for ( std::pair< core::Size, LigandResidueOP > ligand: std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->get_owned_ligand_residues() ) {
				//Reset owner
				ligand.second->set_owner_segment(std::dynamic_pointer_cast< LigandSegment >( new_chimaera ) );
				for ( LigandContactOP contact: ligand.second->get_nonconst_current_contacts() ) {
					if ( contact->segment_id == 0 ) {
						TR.Debug << "Internal contacts do not change" << std::endl;
						continue;
					}
					TR.Debug << "Old residue number: " << contact->residue_number << std::endl;
					runtime_assert( contact->residue_number != 0 );
					//basis_pair.first is N-terminal, so resnums will match up
					if ( contact->segment_id == basis_pair.first.segment_id() ) {
						contact->segment_id = new_chimaera->get_segment_id();
						TR.Debug << "New residue number: " << contact->residue_number << std::endl;
					} else if ( contact->segment_id == basis_pair.second.segment_id() ) {
						//basis_pair.second is C-terminal, so resnums will need to be fixed
						//NEW RESNUM = OLD RESNUM - CTERM BASIS + NTERM BASIS
						contact->segment_id = new_chimaera->get_segment_id();
						contact->residue_number = contact->residue_number + basis_pair.first.resnum() - basis_pair.second.resnum();
						TR.Debug << "New residue number: " << contact->residue_number << std::endl;
					}
				} //end for contact
			}//end for owned ligand
			//END TEMP
			for ( core::Size ligand_id: std::dynamic_pointer_cast< LigandSegment >( new_chimaera )->get_ligand_residues() ) {
				for ( LigandContactOP contact: local_ligands_.at( ligand_id )->get_nonconst_current_contacts() ) {
					if ( contact->segment_id == 0 ) {
						TR.Debug << "Internal contacts and partner contacts do not change" << std::endl;
						continue;
					}
					//TR << "Old residue number: " << contact->residue_number << std::endl;
					runtime_assert( contact->residue_number != 0 );
					//basis_pair.first is N-terminal, so resnums will match up
					//TR << "C terminal basis: " << basis_pair.second.resnum() << " N terminal basis: " << basis_pair.first.resnum() << std::endl;
					if ( contact->segment_id == basis_pair.first.segment_id() ) {
						//TR << "Contact in N-terminal part of chimaera. Residue numbering unchanged." << std::endl;
						contact->segment_id = new_chimaera->get_segment_id();
						//TR << "New residue number: " << contact->residue_number << std::endl;
					} else if ( contact->segment_id == basis_pair.second.segment_id() ) {
						//basis_pair.second is C-terminal, so resnums will need to be fixed
						//TR << "Contact in C-terminal part of chimaera. Fixing residue numbering." << std::endl;
						//NEW RESNUM = OLD RESNUM - CTERM BASIS + NTERM BASIS
						contact->segment_id = new_chimaera->get_segment_id();
						contact->residue_number = contact->residue_number + basis_pair.first.resnum()  - basis_pair.second.resnum();
						// TR << "New residue number: " << contact->residue_number << std::endl;
					}
				}//end for contact
			}//end for owned ligand
		} //end if
	} //END IF ADDING TO C TERMINUS

	/*
	else {// are we trying to add to the middle of our assembly?
	TR << "Segment " << basis_pair.first.segment_id() << " has N-terminal neighbor " << local_segments_.at(basis_pair.first.segment_id())->get_n_terminal_neighbor()->get_segment_id() << std::endl;
	TR << "Segment " << basis_pair.first.segment_id() << " has C-terminal neighbor " << local_segments_.at(basis_pair.first.segment_id())->get_c_terminal_neighbor()->get_segment_id() << std::endl;
	utility_exit_with_message("Cannot chimerize a non-terminal segment!");
	}
	*/
	new_chimaera->set_residue_vector(residues);
	new_chimaera->set_dssp_code(new_chimaera->get_n_terminal_parent()->get_dssp_code());
	// TR << "New chimaera has n-terminal parent " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
	// TR << "New chimaera has c-terminal parent " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
	if ( !new_chimaera->get_n_terminal_parent()->get_n_terminal_neighbor() && !new_chimaera->get_c_terminal_parent()->get_c_terminal_neighbor() ) {
		TR << "Basis Pair produced isolated chimaera! Seg1: " << basis_pair.first.segment_id() << " res1: " << basis_pair.first.resnum()
			<< " Seg2: " << basis_pair.second.segment_id() << " res2: " << basis_pair.second.resnum() << std::endl;
		TR << "n_term parent" << new_chimaera->get_n_terminal_parent() << " " << new_chimaera->get_n_terminal_parent()->get_segment_id() << std::endl;
		TR << "c_term parent" << new_chimaera->get_c_terminal_parent() << " " << new_chimaera->get_c_terminal_parent()->get_segment_id() << std::endl;
		TR << "n_term neighbor" << new_chimaera->get_n_terminal_parent()->get_n_terminal_neighbor() << std::endl;
		TR << "c_term neighbor" << new_chimaera->get_c_terminal_parent()->get_c_terminal_neighbor() << std::endl;
		TR << "n_term neighbor local: " << local_segments_.at( new_chimaera->get_n_terminal_parent()->get_segment_id() )->get_n_terminal_neighbor() << std::endl;
		TR << "c_term neighbor local: " << local_segments_.at( new_chimaera->get_c_terminal_parent()->get_segment_id() )->get_c_terminal_neighbor() << std::endl;
		TR << "n_term neighbor" << new_chimaera->get_n_terminal_neighbor() /*<< " " << new_chimaera->get_n_terminal_neighbor()->get_segment_id() */<<  std::endl;
		TR << "c_term neighbor" << new_chimaera->get_c_terminal_neighbor() /*<< " " << new_chimaera->get_c_terminal_neighbor()->get_segment_id()*/ << std::endl;
		utility_exit_with_message("Isolated chimaera!");
	}
	//if(!new_chimaera->get_c_terminal_parent()->get_c_terminal_neighbor()){
	// utility_exit_with_message("Bad c-terminal neighbor");
	// }
	SmartSegment::link_to(new_chimaera, new_chimaera->get_c_terminal_parent()->get_c_terminal_neighbor());
	SmartSegment::link_to(new_chimaera->get_n_terminal_parent()->get_n_terminal_neighbor(), new_chimaera);
	if ( new_chimaera->get_n_terminal_parent()->is_vital() || new_chimaera->get_c_terminal_parent()->is_vital() ) {
		new_chimaera->set_is_vital( true );
	}
	last_chimaera_ = new_chimaera;
	return true;
}

bool
SmartAssembly::transform_segments( BasisPair basis_pair ) {
	//The first basis is the segment that is staying in place
	//The second basis is the segment that will move (note we don't know which is N/C terminal)

	//We'll only need the basis residue from the first segment to define the coordinate frame
	SmartSewingResidueOP stationary_basis_residue = local_segments_.at( basis_pair.first.segment_id() )->get_residue( basis_pair.first.resnum() );
	utility::vector1< core::conformation::Atom > & stationary_basis_atoms = stationary_basis_residue->get_atom_vector();
	//We'll define a coordinate frame for the mobile basis residue
	SmartSewingResidueOP mobile_basis_residue = local_segments_.at( basis_pair.second.segment_id() )->get_residue( basis_pair.second.resnum() );
	utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue->get_atom_vector();

	//We'll transform all of the residues in all the connected segment to the local coordinate frame of the mobile basis residue and then to the global coordinate frame of the stationary basis residue

	numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
	numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
	numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
	numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;

	SmartSegmentOP current_segment = SmartSegment::get_n_most_segment( local_segments_.at( basis_pair.second.segment_id() ) , false  ); //In theory there shouldn't be any chimaeras here!!! Should we check for that here or elsewhere?

	while ( current_segment != nullptr ) {
		for ( SmartSewingResidueOP current_res: current_segment->get_residue_vector() ) {
			for ( core::conformation::Atom & current_atom: current_res->get_atom_vector() ) {
				//we have to copy this one since the vector within the atom is constant
				numeric::xyzVector< core::Real > old_coords = current_atom.xyz();
				numeric::xyzVector< core::Real > new_coords = mobile_to_stationary_ht * old_coords;
				//current_atom.xyz( mobile_to_stationary_ht * old_coords );
				current_atom.xyz( new_coords );

			}//end for current_atom
		}//end for current_res
		//If the segment is a ligand segment, we'll also need to transform all of its owned ligands
		if ( std::dynamic_pointer_cast< LigandSegment >( current_segment ) != nullptr ) {
			for ( std::pair< core::Size, LigandResidueOP > ligand: std::dynamic_pointer_cast< LigandSegment>( current_segment )->get_owned_ligand_residues() ) {
				for ( core::conformation::Atom & current_atom: ligand.second->get_atom_vector() ) {
					numeric::xyzVector< core::Real > old_coords = current_atom.xyz();
					numeric::xyzVector< core::Real > new_coords = mobile_to_stationary_ht * old_coords;
					//current_atom.xyz( mobile_to_stationary_ht * old_coords );
					current_atom.xyz( new_coords );
				}
			}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
	}//end while current_segment
	//Now that the transformation is complete, we check our alignment

	stationary_basis_atoms = local_segments_.at( basis_pair.first.segment_id() )->get_residue( basis_pair.first.resnum() )->get_atom_vector();
	mobile_basis_atoms = local_segments_.at( basis_pair.second.segment_id() )->get_residue( basis_pair.second.resnum() )->get_atom_vector();

	//Get the edge residues
	utility::vector1< core::conformation::Atom > & stationary_edge_atoms = local_segments_.at( basis_pair.first.segment_id() )->get_residue( basis_pair.first.resnum() + (window_width_ - 1 ) )->get_atom_vector();
	utility::vector1< core::conformation::Atom > &  mobile_edge_atoms = local_segments_.at( basis_pair.second.segment_id() )->get_residue( basis_pair.second.resnum()  + ( window_width_ - 1 ) )->get_atom_vector();

	bool basis_residues_align = ( stationary_basis_atoms[ 1 ].xyz().distance( mobile_basis_atoms[ 1 ].xyz() ) <= 0.5 && stationary_basis_atoms[ 2 ].xyz().distance( mobile_basis_atoms[ 2 ].xyz() ) <= 0.5 && stationary_basis_atoms[ 3 ].xyz().distance( mobile_basis_atoms[ 3 ].xyz() ) <= 0.5 );
	bool edge_residues_align = ( stationary_edge_atoms[ 1 ].xyz().distance( mobile_edge_atoms[ 1 ].xyz() ) <= 0.5 && stationary_edge_atoms[ 2 ].xyz().distance( mobile_edge_atoms[ 2 ].xyz() ) <= 0.5 && stationary_edge_atoms[ 3 ].xyz().distance( mobile_edge_atoms[ 3 ].xyz() ) <= 0.5 );

	//In this version of sewing, the window always goes from the basis residue (N-terminal) to the ends (C-terminal)
	//return edge_resudues_align && basis_residues_align;
	// TR << "Basis distances: " << stationary_basis_atoms[ 1 ].xyz().distance( mobile_basis_atoms[ 1 ].xyz() ) << " " << stationary_basis_atoms[ 2 ].xyz().distance( mobile_basis_atoms[ 2 ].xyz() ) << " " << stationary_basis_atoms[ 3 ].xyz().distance( mobile_basis_atoms[ 3 ].xyz() ) << " " << stationary_basis_atoms[ 4 ].xyz().distance( mobile_basis_atoms[ 4 ].xyz() ) << std::endl;
	//TR << "Edge distances: " << stationary_edge_atoms[ 1 ].xyz().distance( mobile_edge_atoms[ 1 ].xyz() ) << " " << stationary_edge_atoms[ 2 ].xyz().distance( mobile_edge_atoms[ 2 ].xyz() ) << " " << stationary_edge_atoms[ 3 ].xyz().distance( mobile_edge_atoms[ 3 ].xyz() ) << " " << stationary_edge_atoms[ 4 ].xyz().distance( mobile_edge_atoms[ 4 ].xyz() ) <<  std::endl;
	//Tracers for testing purposes
	if ( edge_residues_align && basis_residues_align ) {
		return true;
	}
	if ( !basis_residues_align ) {
		TR.Debug << "Basis residues do not align!!! Something is wrong with transform_segments." << std::endl;
		TR.Debug << "Basis distances: " << stationary_basis_atoms[ 1 ].xyz().distance( mobile_basis_atoms[ 1 ].xyz() ) << " " << stationary_basis_atoms[ 2 ].xyz().distance( mobile_basis_atoms[ 2 ].xyz() ) << " " << stationary_basis_atoms[ 3 ].xyz().distance( mobile_basis_atoms[ 3 ].xyz() ) << " " << stationary_basis_atoms[ 4 ].xyz().distance( mobile_basis_atoms[ 4 ].xyz() ) << std::endl;
	}
	if ( !edge_residues_align ) {
		TR.Debug << "Edge residues do not align, probably due to irregular helices" << std::endl;
		TR.Debug << "Edge distances: " << stationary_edge_atoms[ 1 ].xyz().distance( mobile_edge_atoms[ 1 ].xyz() ) << " " << stationary_edge_atoms[ 2 ].xyz().distance( mobile_edge_atoms[ 2 ].xyz() ) << " " << stationary_edge_atoms[ 3 ].xyz().distance( mobile_edge_atoms[ 3 ].xyz() ) << " " << stationary_edge_atoms[ 4 ].xyz().distance( mobile_edge_atoms[ 4 ].xyz() ) <<  std::endl;
	}
	return false;
}


std::pair< bool, core::Size>
SmartAssembly::iterate_over_basis_pairs(
	data_storage::SmartSegmentCOP segment_1,
	data_storage::SmartSegmentCOP segment_2,
	bool n_terminus)
{
	//Fill in our vector
	//all_basis_pairs_.clear();

	//If either of the segments is not hashable, or if they are different SS types, leave the vector empty.
	if ( !( segment_1->is_hashable() && segment_2->is_hashable()
			&& ( segment_1->get_dssp_code() == segment_2->get_dssp_code() ) )
			) {
		//in this case, we will leave all_basis_pairs_ empty.
		TR.Debug << "Cannot hash segments!" << std::endl;
		if ( !segment_1->is_hashable() ) {
			return std::make_pair( false, 0 );
		} else {
			return std::make_pair(true, 0);
		}
	}//endif



	TR.Debug << "Producing basis pairs within " << utility::to_string( window_width_ ) << " residues of termini" << std::endl;
	//Find number of residues in each segment
	core::Size seg1_size = segment_1->get_length();
	core::Size seg2_size = segment_2->get_length();
	// TR << "seg1_size " << seg1_size << " seg2_size " << seg2_size << std::endl;
	//If either segment is shorter than required_num_residues, no basis pairs possible.
	if ( seg1_size < window_width_ || seg2_size < window_width_ ) {
		TR.Debug << "One of the segments is too short!" << std::endl;
		if ( seg1_size < window_width_ ) {
			//Seg1 is bad
			return std::make_pair( false, 0 );
		} else {
			//Seg1 is OK
			return std::make_pair( true, 0 );
		}
	}
	//i is index within the first segment (already in the assembly, might be N or C terminal)
	//WINDOW WIDTH INCLUDES THE BASIS RESIDUE
	core::Size i_start = 1;
	core::Size i_end  = ( seg1_size - window_width_ + 1 ); //We add 1 to account for the basis residue, which IS included in window_width_



	if ( segment_1->is_chimaeric() ) { //if we are adding to a chimaera it will always be segment_1 since it'll already be in the assembly
		//TR << "Adding to chimaera" << std::endl;
		core::Size basis_res;
		//we have to find which Basis is n_terminal
		//If the first basis is N-terminal, the numbering will be fine
		if ( segment_1->get_n_terminal_parent()->get_segment_id() == segment_1->get_basis_pair().first.segment_id() ) {
			basis_res = segment_1->get_basis_pair().first.resnum();
		} else {
			//If the first basis is C-terminal, we should go by the N-terminal numbering
			basis_res = segment_1->get_basis_pair().second.resnum();
		}
		//TR << "Chimaera has basis res of: " << basis_res << std::endl;
		if ( n_terminus ) {
			//we are adding to the n-term
			//This means that seg1 is actually C-terminal
			i_end = basis_res - window_width_;//Basis residue i is INCLUDED in window_width
		} else {
			//we are adding to c term, so this segment's basis will be closer to the end of the segment
			i_start = basis_res + window_width_; //Must be window_width_ away to avoid chimerization bug
		}
	}

	//TR << "i_start: " << i_start << " i_end: " << i_end << " window_width: " << window_width_ << std::endl;
	//Now we can check for chimerizing over required residues in seg1
	std::set< core::Size > seg1_vital_res = segment_1->get_vital_residues();
	std::set< core::Size > seg2_vital_res = segment_2->get_vital_residues();
	if ( seg1_vital_res.size() != 0 ) {
		TR.Debug << "Found vital residues in seg1 "<< std::endl;
		if ( n_terminus ) { //Seg1 is C TERMINAL! We should be closer to the BEGINNING of the segment
			//std::set is a sorted container, so this should work
			//If seg1 is C-teminal, then the new basis residue should be <= ( first vital residue - window_width )
			core::Size first_vital_residue = *(seg1_vital_res.begin() );
			i_end = std::min( i_end, first_vital_residue - window_width_ );
			TR.Debug << "Seg1 is c terminal, so we must chimerize before the first vital residue " << std::endl;
			TR.Debug << "Listing vital residues: " << std::endl;
			for ( core::Size vital: seg1_vital_res ) {
				TR.Debug << vital << " ";
			}
			TR.Debug << std::endl;
			TR.Debug << "First vital residue is  " << first_vital_residue << std::endl;
			TR.Debug << "Last possible basis residue based on vital res is " << i_end << std::endl;

		} else { //SEG1 is N TERMINAL
			core::Size last_vital_residue = *(seg1_vital_res.rbegin() ); //We should be closer to the END of the segment
			i_start = std::max( i_start, last_vital_residue + window_width_ );
			TR.Debug << "Seg1 is n terminal, so we must chimerize after the last vital residue " << std::endl;
			TR.Debug << "Listing vital residues: " << std::endl;
			for ( core::Size vital: seg1_vital_res ) {
				TR.Debug << vital << " ";
			}
			TR.Debug << std::endl;
			TR.Debug << "Last vital residue is  " << last_vital_residue << std::endl;
			TR.Debug << "First possible basis residue based on vital res is " << i_start << std::endl;
		}
	}
	//TR << "scanning seg1 residues : "<< i << " to " << i_end << std::endl;
	//Now we can check for chimerizing over required residues in seg2
	core::Size j_start = 1;
	core::Size j_end =  (seg2_size - window_width_ + 1 );
	//WE SHOULD NEVER BE CHIMERIZING TWO CHIMAERAE TO EACH OTHER. That's why we don't account for it here. If you're trying to do that, you're doing something wrong--your segments may be ordered wrong, or your basis pair was set incorrectly.
	runtime_assert(! segment_2->is_chimaeric() );
	//If seg2 is N-terminal, then the new basis residue should be >= ( last vital residue + window_width)
	//If seg2 is C-teminal, then the new basis residue should be <= ( first vital residue - window_width)
	if ( seg2_vital_res.size() != 0 ) {
		if ( n_terminus ) { //Seg2 is N terminal
			//If seg2 is N-terminal, then the new basis residue should be >= ( last vital residue + window_width_) (we count WW back toward LVR)
			core::Size last_vital_residue = *(seg2_vital_res.rbegin() );
			j_start = std::max( j_start, last_vital_residue + window_width_ ); //so this was right to start with even though the comments are all wrong
			TR.Debug << "Seg2 is n terminal, so we must chimerize after the last vital residue " << std::endl;
			TR.Debug << "Listing vital residues: " << std::endl;
			for ( core::Size vital: seg2_vital_res ) {
				TR.Debug << vital << " ";
			}
			TR.Debug << std::endl;
			TR.Debug << "Last vital residue is  " << last_vital_residue << std::endl;
			TR.Debug << "First possible basis residue based on vital res is " << j_start << std::endl;
		} else {
			core::Size first_vital_residue = *(seg2_vital_res.begin() );
			//If seg1 is C-teminal, then the new basis residue should be <= ( first vital residue - window_width )
			j_end = std::min( j_end, first_vital_residue - window_width_ );
			TR.Debug << "Seg2 is c terminal, so we must chimerize before the first vital residue " << std::endl;
			TR.Debug << "Listing vital residues: " << std::endl;
			for ( core::Size vital: seg2_vital_res ) {
				TR.Debug << vital << " ";
			}
			TR.Debug << std::endl;
			TR.Debug << "First vital residue is  " << first_vital_residue << std::endl;
			TR.Debug << "Last possible basis residue based on vital res is " << j_end << std::endl;
		}
	}
	// TR << "scanning seg2 residues : "<< j_start << " to " << j_end << std::endl;
	if ( i_start > i_end || segment_1->get_dssp_code() == 'L' ) {
		TR.Debug << "Seg1 is bad!" << std::endl;
		TR.Debug << "Start: " << i_start << " End: " << i_end << std::endl;
		TR.Debug << "Listing vital residues: " << std::endl;
		for ( core::Size vital: seg1_vital_res ) {
			TR.Debug << vital << " ";
		}
		TR.Debug << std::endl;
		return std::make_pair( false, 0);
	}
	if ( j_start > j_end || segment_2->get_dssp_code() == 'L' ) {
		return std::make_pair( true, 0);
	}
	/*
	core::Size i_first_vital_residue = 0;
	core::Size i_last_vital_residue = seg2_size + 1;
	core::Size j_first_vital_residue = 0;
	core::Size j_last_vital_residue = seg2_size + 1;
	if( seg1_vital_res.size() != 0 ){
	i_first_vital_residue = *(seg1_vital_res.begin() );
	i_last_vital_residue = *(seg1_vital_res.rbegin() );
	}
	if( seg2_vital_res.size() != 0 ){
	j_first_vital_residue = *(seg2_vital_res.begin() );
	j_last_vital_residue = *(seg2_vital_res.rbegin() );
	}
	*/
	core::Size index_of_last_entry = 0;
	//We need to check that one segment is N-terminal and one segment is C-terminal, but we don't really need to know which is which.
	if ( !( segment_1->is_c_terminus_fixed() || segment_2->is_n_terminus_fixed() ) || !( segment_2->is_c_terminus_fixed() || segment_1->is_n_terminus_fixed() ) ) {
		//TR.Debug <<" Checking segment alignments" << std::endl;
		for ( core::Size i = i_start; i <= i_end; ++i ) {
			for ( core::Size j = j_start; j <= j_end; ++j ) {
				//Only insert basis pairs that could possibly produce an overlap of the size we want
				core::Size max_overlap_size = utility::min( seg1_size - i + 1, seg2_size - j + 1 );
				//    TR << "maximum detected overlap " << max_overlap_size << std::endl;
				if ( max_overlap_size >= window_width_ ) {
					//TR << "Adding basis pair" << std::endl;
					++index_of_last_entry;
					all_basis_pairs_[ index_of_last_entry ] = std::make_pair( i, j );
				}//end if
			}//end for j
		}//end for i
	}//end if
	//this will give you some iterators, yay!
	return std::make_pair( true, index_of_last_entry );
}//end iterate_over_basis_pairs function



//BUG IS HERE

void
SmartAssembly::add_segment_and_neighbors_to_local_segments( SmartSegmentCOP ref_seg, core::Size new_id ){
	utility::vector1< LigandSegmentOP > ligsegs;
	std::map< core::Size, core::Size > copied_ligands; //Maps old ligand IDs to new ligand IDs
	core::Size max_lig_id = 0;
	if ( local_ligands_.size() > 0 ) {
		max_lig_id = local_ligands_.rbegin()->first;
	}
	SmartSegmentCOP current_seg = SmartSegment::get_n_most_segment( ref_seg, false );
	TR.Debug << "Beginning add_segment_and_neighors_to_local_segments " << std::endl;
	TR.Debug << "Listing vital residues for reference: " << std::endl;
	for ( core::Size vital: current_seg->get_vital_residues() ) {
		TR.Debug << vital << " ";
	}
	TR.Debug << std::endl;
	core::Size first_new_id = new_id;
	core::Size first_old_id = current_seg->get_segment_id();
	if ( new_id == 0 ) {
		new_id = first_old_id;
	}
	SmartSegmentOP new_seg;
	while ( current_seg != nullptr ) {
		if ( local_segments_.count( current_seg->get_segment_id() ) == 0 ) { //This segment is not already in local_segments
			if ( std::dynamic_pointer_cast< LigandSegment const >( current_seg ) != nullptr ) {
				//This makes sure the ligands are copied appropriately
				LigandSegmentCOP ref_ligseg = std::dynamic_pointer_cast< LigandSegment const >( current_seg );
				LigandSegmentOP ligseg = ref_ligseg->clone();
				for ( std::pair< core::Size, LigandResidueOP > ligand: ligseg->get_owned_ligand_residues() ) {
					//set the owner
					ligand.second->set_owner_segment( ligseg );
					//We should ONLY do this if we're also giving our segment a new ID
					if ( first_new_id != 0 ) {
						//Give the ligand a new ID
						copied_ligands[ ligand.first ] = max_lig_id + 1;
						ligand.second->set_ligand_id( max_lig_id + 1 );
						//local_ligands_[ max_lig_id + 1 ] = ligand;
						++max_lig_id;
					}
					if ( local_ligands_.count( ligand.second->get_ligand_id() ) == 0 ) {
						local_ligands_[ ligand.second->get_ligand_id() ] = ligand.second;
						ligand_conformers_[ ligand.second->get_ligand_id() ] = ligand_conformers_[ ligand.first ]; //If the ligand id didn't change, this will do nothing. If it changed, it will copy the conformers.
						for ( LigandResidueOP conformer: ligand_conformers_.at( ligand.second->get_ligand_id() ) ) {
							conformer->set_ligand_id( ligand.second->get_ligand_id() );
							conformer->set_owner_segment( ligseg );
							for ( LigandContactOP contact: conformer->get_nonconst_current_contacts() ) {
								if ( contact->segment_id == 0 ) {
									TR.Debug << "Internal contacts do not change" << std::endl;
									continue;
								}
								//TR << "Conformer contact residue number: " << contact->residue_number << std::endl;
								runtime_assert( contact->residue_number != 0 );
								if ( contact->segment_id == current_seg->get_segment_id() ) {
									contact->segment_id = new_id;
								}
							}//end iterate over contacts
						}
					}
					//Set the appropriate segment ids in all the contacts
					//This is important not just for the owned contacts but also for the nonowned contacts that are also being transformed
					for ( LigandContactOP contact: ligand.second->get_nonconst_current_contacts() ) {
						if ( contact->segment_id == 0 ) {
							TR.Debug << "Internal contacts and partner contacts do not change" << std::endl;
							continue;
						}
						//TR << "Contact residue number: " << contact->residue_number << std::endl;
						runtime_assert( contact->residue_number != 0 );
						if ( contact->segment_id == current_seg->get_segment_id() ) {
							contact->segment_id = new_id;
						}
					}//end iterate over contacts
				}//end iterate over owned ligands
				//There could still be a problem with ligands whose owners have changed IDs
				ligsegs.push_back( ligseg );
				new_seg = SmartSegmentOP( ligseg );
				TR.Debug << "Listing vital residues for ref_ligseg: " << std::endl;
				for ( core::Size vital: ref_ligseg->get_vital_residues() ) {
					TR.Debug << vital << " ";
				}
				//new_seg->set_vital_residues( ligseg->get_vital_residues() );
				TR.Debug << "Listing vital residues for new: " << std::endl;
				for ( core::Size vital: new_seg->get_vital_residues() ) {
					TR.Debug << vital << " ";
				}
				TR.Debug << std::endl;
				runtime_assert( ligseg->get_ligand_residues().size() == ref_ligseg->get_ligand_residues().size() );
			} else { // end ligand segments
				new_seg = current_seg->clone();
				//new_seg->set_vital_residues( current_seg->get_vital_residues() );
				TR.Debug << "Listing vital residues for new: " << std::endl;
				for ( core::Size vital: new_seg->get_vital_residues() ) {
					TR.Debug << vital << " ";
				}
				TR.Debug << std::endl;
			}
			new_seg->set_segment_id( new_id );


			//TODO
			if ( current_seg->get_const_reference_segment() != nullptr ) {
				new_seg->set_const_reference_segment( current_seg->get_const_reference_segment() );
			} else {
				new_seg->set_const_reference_segment( current_seg );
			}
			runtime_assert( local_segments_.count( new_id ) == 0 );
			local_segments_[ new_id ] = new_seg;

			//   TR << "Added new segment to local segments " << new_id << "(" << current_seg->get_segment_id() << ") : " << local_segments_.count( new_id ) << std::endl;

			//Possible issue!
			//TODO

			//In theory, current_seg should always be from the constant segments and should not be linked to any chimaerae
			if ( current_seg->get_const_reference_segment() != nullptr ) { ///This segment needs to worry about chimaerae
				if ( current_seg->get_const_reference_segment()->is_n_terminus_fixed() && local_segments_.count( new_id - 1 ) != 0 ) {
					SmartSegment::link_to( local_segments_.at( new_id - 1 ), local_segments_.at( new_id ) );
				}
				if ( current_seg->get_const_reference_segment()->is_c_terminus_fixed() && local_segments_.count( new_id + 1 ) != 0 ) {
					SmartSegment::link_to( local_segments_.at( new_id ), local_segments_.at( new_id + 1 ) );
				}
			} else {
				if ( current_seg->is_n_terminus_fixed() && local_segments_.count( new_id - 1 ) != 0 ) {
					SmartSegment::link_to( local_segments_.at( new_id - 1 ), local_segments_.at( new_id ) );
				}
				if ( current_seg->is_c_terminus_fixed() && local_segments_.count( new_id + 1 ) != 0 ) {
					SmartSegment::link_to( local_segments_.at( new_id ), local_segments_.at( new_id + 1 ) );
				}
			}
		} else { //end not already in local segments// it's already in local segments
			if ( first_new_id != 0 && local_segments_.count( new_id ) == 0 ) { // storing a copy? ( second half should be redundant... )
				new_seg = current_seg->clone();
				//new_seg->set_vital_residues( current_seg->get_vital_residues() );
				TR.Debug << "Listing vital residues for new segment: " << std::endl;
				for ( core::Size vital: new_seg->get_vital_residues() ) {
					TR.Debug << vital << " ";
				}
				TR.Debug << std::endl;
				new_seg->set_segment_id( new_id );
				if ( current_seg->get_const_reference_segment() != nullptr ) {
					current_seg = current_seg->get_const_reference_segment();
				}
				new_seg->set_const_reference_segment( current_seg );
				local_segments_[ new_id ] = new_seg;
				//    TR << "Added Copy segment to local segments " << new_id << "(" << current_seg->get_segment_id() << ") : " << local_segments_.count( new_id ) << std::endl;
				if ( current_seg->get_const_reference_segment() != nullptr ) { ///This segment needs to worry about chimaerae
					if ( current_seg->get_const_reference_segment()->is_n_terminus_fixed() && local_segments_.count( new_id - 1 ) != 0 ) {
						SmartSegment::link_to( local_segments_.at( new_id - 1 ), local_segments_.at( new_id ) );
					}
					if ( current_seg->get_const_reference_segment()->is_c_terminus_fixed() && local_segments_.count( new_id + 1 ) != 0 ) {
						SmartSegment::link_to( local_segments_.at( new_id ), local_segments_.at( new_id + 1 ) );
					}
				} else {
					if ( current_seg->is_n_terminus_fixed() && local_segments_.count( new_id - 1 ) != 0 ) {
						SmartSegment::link_to( local_segments_.at( new_id - 1 ), local_segments_.at( new_id ) );
					}
					if ( current_seg->is_c_terminus_fixed() && local_segments_.count( new_id + 1 ) != 0 ) {
						SmartSegment::link_to( local_segments_.at( new_id ), local_segments_.at( new_id + 1 ) );
					}
				}
			}// end store copy
		}//end already in local segments
		current_seg = current_seg->get_c_terminal_neighbor();
		++new_id;
	}// while loop
	//Now fix up ligand segment unowned contacts
	if ( first_new_id != 0 ) { //If none of the segments were renumbered, then the ligands should not have been renumbered either
		for ( LigandSegmentOP ligseg: ligsegs ) {
			TR.Debug << "Listing vital residues for ligand segment: " << std::endl;
			for ( core::Size vital: ligseg->get_vital_residues() ) {
				TR.Debug << vital << " ";
			}
			TR.Debug << std::endl;
			//Begin fix up renumbered owners for unowned ligands
			//Use our ligand renumbering map
			//First fix up any bad ligand IDs
			std::set< core::Size > & unowned_ligands = ligseg->get_nonconst_ligand_residues();
			for ( std::pair< core::Size, core::Size > renumbered_ligand: copied_ligands ) {
				if ( unowned_ligands.count( renumbered_ligand.first ) != 0 ) {
					unowned_ligands.erase( unowned_ligands.find( renumbered_ligand.first ) );
					unowned_ligands.insert( renumbered_ligand.second );
				}
			}
			//Any contact referring to the old segment ID for ligseg should take the new segment ID for ligseg instead
			//Old range:first_old_id to (first_old_id + new_id - first_new_id - 1 )
			//New range:first_new_id to new_id - 1
			//Current value of new_id is one more than the last new segment id
			core::Size num_new_segments = new_id - first_new_id; //1 to 4 = 3 new segments 1, 2, and 3

			//If a ligand isn't owned by a copied residue, then it isn't moving, and its contacts haven't changed.
			//for( core::Size unowned_ligand: unowned_ligands ){
			for ( std::pair< core::Size, LigandResidueOP > owned_ligand: ligseg->get_owned_ligand_residues() ) {
				for ( LigandContactOP contact: local_ligands_.at( owned_ligand.first )->get_nonconst_current_contacts() ) {
					if ( contact->segment_id == 0 ) {
						TR.Debug << "Internal contacts do not change" << std::endl;
						continue;
					}
					//TR << "Contact residue number: " << contact->residue_number << std::endl;
					runtime_assert( contact->residue_number != 0 );
					//Need to map all old segment ids to new segment IDs
					if ( contact->segment_id >= first_old_id && contact->segment_id < first_old_id + num_new_segments ) {
						contact->segment_id = first_new_id + ( contact->segment_id - first_old_id );
					}
				} //end for contact
				//Now fix contacts of conformers
				for ( LigandResidueOP conformer: ligand_conformers_.at( owned_ligand.first ) ) {
					for ( LigandContactOP contact: conformer->get_nonconst_current_contacts() ) {
						if ( contact->segment_id == 0 ) {
							TR.Debug << "Internal contacts do not change" << std::endl;
							continue;
						}
						//TR << "Conformer contact residue number: " << contact->residue_number << std::endl;
						runtime_assert( contact->residue_number != 0 );
						if ( contact->segment_id >= first_old_id && contact->segment_id < first_old_id + num_new_segments ) {
							contact->segment_id = first_new_id + ( contact->segment_id - first_old_id );
						}
					}
				}
			}//end for unowned ligand
		}//End for ligsegs
	}//End if segment id changed
}

bool
SmartAssembly::sample_ligand()
{
	//Pick a random ligand to sample
	core::Size index = std::floor( numeric::random::rg().uniform() * local_ligands_.size() );
	std::map< core::Size, LigandResidueOP>::iterator lig_it = local_ligands_.begin();
	std::advance( lig_it, index );
	core::Size ligand_id_to_sample = lig_it->first;
	if ( local_ligands_[ ligand_id_to_sample ]->get_partner_ligand() ) {
		return false;
	}
	//Pick a random conformer for that ligand
	LigandResidueOP conformer_to_sample = numeric::random::rg().random_element( ligand_conformers_.at( ligand_id_to_sample ) );
	last_sampled_ligand_ = std::make_pair( ligand_id_to_sample, local_ligands_.at( ligand_id_to_sample ) );
	//Make sure the new conformer has the right owner!
	conformer_to_sample->set_owner_segment( last_sampled_ligand_.second->get_owner_segment() );
	//Now we need to transform the conformer
	//We'll want to align the conformer to the copy in local_ligands_ using the specified alignment atoms
	LigandResidueCOP current_conformer = LigandResidueCOP( local_ligands_.at( ligand_id_to_sample ) );
	utility::vector1< core::Size > alignment_atoms = conformer_to_sample->get_alignment_atoms();
	utility::vector1< core::conformation::Atom > const & stationary_atoms = current_conformer->get_const_atom_vector();
	utility::vector1< core::conformation::Atom > & mobile_atoms = conformer_to_sample->get_atom_vector();
	numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_atoms.at( alignment_atoms[ 1 ] ).xyz(), stationary_atoms.at( alignment_atoms[ 2 ] ).xyz(), stationary_atoms.at( alignment_atoms[ 3 ] ).xyz() );
	numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_atoms.at( alignment_atoms[ 1 ] ).xyz(), mobile_atoms.at( alignment_atoms[ 2 ] ).xyz(), mobile_atoms.at( alignment_atoms[ 3 ] ).xyz() );
	numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
	numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;
	//Now transform the atoms
	for ( core::conformation::Atom & current_atom: mobile_atoms ) {
		numeric::xyzVector< core::Real > old_coords = current_atom.xyz();
		numeric::xyzVector< core::Real > new_coords = mobile_to_stationary_ht * old_coords;
		current_atom.xyz( new_coords );
	}
	//Change in local_ligands so that non-owner segments will know about the change
	local_ligands_[ ligand_id_to_sample ] = conformer_to_sample;
	//Change the pointer in the owner segment (it's the only one that should have a pointer to this ligand anymore)
	LigandSegmentOP owner_ligseg = conformer_to_sample->get_nonconst_owner_segment();
	owner_ligseg->get_owned_ligand_residues()[ ligand_id_to_sample ] = conformer_to_sample;
	//If the owner segment is a chimaera, we'll also need to change the pointers of its parents
	if ( owner_ligseg->is_chimaeric() ) {
		LigandSegmentOP nparent = std::dynamic_pointer_cast< LigandSegment >( owner_ligseg->get_n_terminal_parent() );
		LigandSegmentOP cparent = std::dynamic_pointer_cast< LigandSegment >( owner_ligseg->get_c_terminal_parent() );
		if ( nparent != nullptr && nparent->get_owned_ligand_residues().count( ligand_id_to_sample ) != 0 ) {
			nparent->get_owned_ligand_residues()[ ligand_id_to_sample ] = conformer_to_sample;
		}
		if ( cparent != nullptr && cparent->get_owned_ligand_residues().count( ligand_id_to_sample ) != 0 ) {
			cparent->get_owned_ligand_residues()[ ligand_id_to_sample ] = conformer_to_sample;
		}
	}
	last_change_ = 'L';
	return true;
}

void
SmartAssembly::unsample_ligand(){
	//Fix in local_ligands_
	local_ligands_[ last_sampled_ligand_.first ] = last_sampled_ligand_.second;
	//Fix in owner
	LigandSegmentOP owner_ligseg = last_sampled_ligand_.second->get_nonconst_owner_segment();

	owner_ligseg->get_owned_ligand_residues()[ last_sampled_ligand_.first ] = last_sampled_ligand_.second;

	//If owner is a chimaera, fix in parents
	if ( owner_ligseg->is_chimaeric() ) {
		LigandSegmentOP nparent = std::dynamic_pointer_cast< LigandSegment >( owner_ligseg->get_n_terminal_parent() );
		LigandSegmentOP cparent = std::dynamic_pointer_cast< LigandSegment >( owner_ligseg->get_c_terminal_parent() );
		if ( nparent != nullptr && nparent->get_owned_ligand_residues().count( last_sampled_ligand_.first ) != 0 ) {
			nparent->get_owned_ligand_residues()[ last_sampled_ligand_.first ] = last_sampled_ligand_.second;
		}
		if ( cparent != nullptr && cparent->get_owned_ligand_residues().count( last_sampled_ligand_.first ) != 0 ) {
			cparent->get_owned_ligand_residues()[ last_sampled_ligand_.first ] = last_sampled_ligand_.second;
		}
	}
}


void
SmartAssembly::load_initial_conformers( data_storage::LigandDescription ligdes ){//core::Size ligand_id, utility::vector1< core::conformation::ResidueCOP > additional_conformers ){
	//Additional conformers will be provided in the Ligand tag as a path to a file containing a list of PDBs
	//The ligand should be the only residue in the PDB
	//If the ligand is already in local_ligands, get the first conformation there. Otherwise get it from const_ligands.
	/*
	if ( ligand_conformers_.count( ligdes.ligand_id ) != 0 ) {
	utility_exit_with_message( "ERROR: load_initial_conformers should only be called once per ligand!" );
	}
	*/
	if ( local_ligands_.count( ligdes.ligand_id ) != 0 ) {
		ligand_conformers_[ ligdes.ligand_id ].push_back( local_ligands_.at( ligdes.ligand_id ) );
	} else {
		if ( const_ligands_.count( ligdes.ligand_id ) > 0 ) { //We don't need to worry about partner ligands--no conformer subs allowed
			ligand_conformers_[ ligdes.ligand_id ].push_back( const_ligands_.at( ligdes.ligand_id )->clone() );
		}
	}
	//The rest of the conformers will be cloned from this first one
	for ( core::conformation::ResidueCOP ligpose: ligdes.pdb_conformers ) {
		//I called it ligpose because I originally stored it as a pose, but then I realized that was dumb.
		LigandResidueOP new_conformer = ligand_conformers_.at( ligdes.ligand_id ).at( 1 )->clone();
		//The ligand id, amino acid type, and type will not change
		utility::vector1< core::conformation::Atom > res_atoms( ligpose->atoms().begin(), ligpose->atoms().end() );
		new_conformer->set_atom_vector( res_atoms );
		ligand_conformers_[ ligdes.ligand_id ].push_back( new_conformer );
	}
}
void
SmartAssembly::reconstitute_assembly_from_string(std::string assembly_string){
	TR << this->get_comprehensive_forward_assembly() << std::endl;
	TR << "Recovering." << std::endl;
	TR << "in_string: " << assembly_string << std::endl;

	utility::vector1<std::string> tokens = utility::string_split(assembly_string);
	if ( tokens.size()>=5 ) {
		utility::vector1<SmartSegmentOP> chimaerae;
		utility::vector1<SmartSegmentOP> n_term_segments;
		utility::vector1<core::Size> n_term_residues;
		utility::vector1<SmartSegmentOP> c_term_segments;
		utility::vector1<core::Size> c_term_residues;
		core::Size current = 1;
		while ( (tokens.size() - current) >= 4 ) {
			chimaerae.push_back(local_segments_.at(stoi(tokens[current])));
			n_term_segments.push_back(local_segments_.at(stoi(tokens[current+1])));
			n_term_residues.push_back(stoi(tokens[current+2]));
			c_term_segments.push_back(local_segments_.at(stoi(tokens[current+3])));
			c_term_residues.push_back(stoi(tokens[current+4]));
			current = current + 5;
		}
		//diagnostic tracer spam
		//  current=1;
		// while(current <= chimaerae.size());
		//delete all segments currently in assembly
		while ( !this->get_n_terminal_segment()->is_vital() ) {
			this->delete_segment(true);
		}
		while ( !this->get_c_terminal_segment()->is_vital() ) {
			this->delete_segment(false);
		}
		this->set_starting_segment(n_term_segments[1], this->get_start_node_vital_segments());
		//reconstitute

		current=1;
		BasisPair current_basis_pair;
		while ( current<=chimaerae.size() ) {
			TR << "Setting chimaera " << current << std::endl;
			current_basis_pair.first.segment_id(chimaerae[current]->get_segment_id());
			current_basis_pair.first.resnum(n_term_residues[current]);
			current_basis_pair.second.segment_id(n_term_segments[current]->get_segment_id());
			current_basis_pair.second.resnum(n_term_residues[current]);
			TR << "Transfoming n terminal segment " << current << std::endl;
			transform_segments(current_basis_pair);

			current_basis_pair.second.segment_id(c_term_segments[current]->get_segment_id());
			current_basis_pair.second.resnum(c_term_residues[current]);
			TR << "Transfoming c terminal segment " << current << std::endl;
			transform_segments(current_basis_pair);

			if ( !n_term_segments[current]->is_n_terminus_fixed() ) {
				this->set_n_terminal_segment(chimaerae[current]);
			} else {
				SmartSegment::link_to(n_term_segments[current]->get_n_terminal_neighbor(),chimaerae[current]);
			}

			if ( !c_term_segments[current]->is_c_terminus_fixed() ) {
				this->set_c_terminal_segment(chimaerae[current]);
			} else {
				SmartSegment::link_to(chimaerae[current],c_term_segments[current]->get_c_terminal_neighbor());
			}
			current = current + 1;
		}
		this->set_c_terminal_segment(SmartSegment::get_c_most_segment(this->get_n_terminal_segment(),true));
		//  SmartSegmentOP current_segment = this->get_n_terminal_segment();
		//  while(current_segment!=nullptr){
		//   TR << "C term segment is now " << current_segment->get_segment_id() << std::endl;
		//   current_segment = current_segment->get_c_terminal_neighbor();
		//  }
		//  this->set_c_terminal_segment(current_segment);
	}
	/* //restore for chimaera-destroying behavior
	current=1;
	BasisPair current_basis_pair;
	while(current<=chimaerae.size()){
	current_basis_pair.first.segment_id(n_term_segments[current]->get_segment_id());
	current_basis_pair.first.resnum(n_term_residues[current]);
	current_basis_pair.second.segment_id(c_term_segments[current]->get_segment_id());
	current_basis_pair.second.resnum(c_term_residues[current]);
	this->chimerize(current_basis_pair,false);
	current = current + 1;
	}
	*/
	//return
}
//These should really only be used in HashedSmartAssembly

void
SmartAssembly::trim_terminal_segments(char terminus, core::Size number){
	while ( number>0 ) {
		if ( !(this->get_n_terminal_segment()->is_vital()) && (terminus == 'N' || terminus == 'B') ) {
			this->set_n_terminal_segment(this->get_n_terminal_segment()->get_c_terminal_neighbor());
			this->get_n_terminal_segment()->get_n_terminal_neighbor()->isolate();
		}
		if ( !(this->get_c_terminal_segment()->is_vital()) && (terminus == 'C' || terminus == 'B') ) {
			this->set_c_terminal_segment(this->get_c_terminal_segment()->get_n_terminal_neighbor());
			this->get_c_terminal_segment()->get_c_terminal_neighbor()->isolate();
		}
		number--;
	}
}

void
SmartAssembly::set_resID_1( core::Size id){
	resID_1_ = id;
}
void
SmartAssembly::set_resID_2( core::Size id ){
	resID_2_ = id;
}
core::Size
SmartAssembly::get_resID_1() const{
	return resID_1_;
}
core::Size
SmartAssembly::get_resID_2() const{
	return resID_2_;
}

void
SmartAssembly::set_segID_1( core::Size id){
	segID_1_ = id;
}
void
SmartAssembly::set_segID_2( core::Size id ){
	segID_2_ = id;
}
core::Size
SmartAssembly::get_segID_1() const{
	return segID_1_;
}
core::Size
SmartAssembly::get_segID_2() const{
	return segID_2_;
}

core::Size
SmartAssembly::get_window_width() const{
	return window_width_;
}
void
SmartAssembly::set_window_width( core::Size width ){
	window_width_ = width;
}


} //protocols
} //sewing
} //data_storage
