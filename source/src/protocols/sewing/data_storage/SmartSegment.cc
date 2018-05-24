// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/SmartSegment.cc
/// @brief a neighbor-aware SewSegment version
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sewing.SmartSegment" );


namespace protocols {
namespace sewing {
namespace data_storage {

SmartSegment::SmartSegment():
	utility::pointer::ReferenceCount()
{ // basic constructor to make an empty SmartSegment. Intended to be used to populate the master SmartSegment vector that serves as a segment library
	segment_id_ = 0;
	n_terminal_neighbor_ = nullptr;
	c_terminal_neighbor_ = nullptr;
	dssp_code_ = 'X';
	is_a_chimaera_ = false;
	n_terminal_parent_ = nullptr;
	c_terminal_parent_ = nullptr;
	is_vital_= false;
	is_in_Assembly_ = false;
	length_ = 0;
}

SmartSegment::SmartSegment(bool is_vital):
	utility::pointer::ReferenceCount()
{ // basic constructor to make an empty SmartSegment. Intended to be used to populate the master SmartSegment vector that serves as a segment library
	segment_id_ = 0;
	n_terminal_neighbor_ = nullptr;
	c_terminal_neighbor_ = nullptr;
	dssp_code_ = 'X';
	is_a_chimaera_ = false;
	n_terminal_parent_ = nullptr;
	c_terminal_parent_ = nullptr;
	is_vital_= is_vital;
	is_in_Assembly_ = false;
	length_ = 0;
}
SmartSegment::SmartSegment(bool is_vital, core::Size max_segment_length):
	utility::pointer::ReferenceCount()
{ // basic constructor to make a SmartSegment suitable for use in an Assembly or Hasher . The important part is the preallocation of the Residue array.
	segment_id_ = 0;
	n_terminal_neighbor_ = nullptr;
	c_terminal_neighbor_ = nullptr;
	dssp_code_ = 'X';
	is_a_chimaera_ = false;
	n_terminal_parent_ = nullptr;
	c_terminal_parent_ = nullptr;
	is_vital_= is_vital;
	is_in_Assembly_ = false;
	length_ = 0;
	for ( core::Size current_residue_number = 1; current_residue_number <= (2*max_segment_length); current_residue_number++ ) {
		SmartSewingResidueOP current_residue = SmartSewingResidueOP(new SmartSewingResidue);
		residues_.push_back(current_residue);
	}


}

SmartSegment::~SmartSegment(){
	//n_terminal_neighbor_->set_c_terminal_neighbor( nullptr );
	//n_terminal_neighbor_ = nullptr;
	//c_terminal_neighbor_->set_n_terminal_neighbor( nullptr );
	//c_terminal_neighbor_ = nullptr;
	//residues_.clear();
	//None of these should have pointers to this segment anyway
	//n_terminal_parent_ = nullptr;
	//c_terminal_parent_ = nullptr;
	//const_reference_segment_ = nullptr;
}

//For use with Hasher and Assembly
SmartSegment::SmartSegment( SmartSegment const & other ) {
	segment_id_ = other.get_segment_id();
	n_terminal_neighbor_ = nullptr;
	c_terminal_neighbor_ = nullptr;
	//Deep copy the residue vector
	residues_.clear();
	for ( data_storage::SmartSewingResidueOP res: other.get_const_residue_vector() ) {
		residues_.push_back( data_storage::SmartSewingResidueOP( new data_storage::SmartSewingResidue( *res ) ) );
	}
	dssp_code_ = other.get_dssp_code();
	is_a_chimaera_ = other.is_chimaeric();
	if ( is_a_chimaera_ ) {
		n_terminal_parent_ = other.get_n_terminal_parent();
		c_terminal_parent_= other.get_c_terminal_parent();
		basis_pair_ = other.get_basis_pair();
	} else {
		n_terminal_parent_ = nullptr;
		c_terminal_parent_ = nullptr;
	}
	is_vital_ = other.is_vital();
	is_in_Assembly_ = false;
	length_ = other.get_length();
	if ( other.get_const_reference_segment() != nullptr ) {
		const_reference_segment_ = other.get_const_reference_segment();
	}
	for ( core::Size i: other.vital_residues_ ) {
		vital_residues_.insert( i );
	}
	//Otherwise we'll probably want to set the other one as our const reference segment, but not necessarily
}

SmartSegmentOP
SmartSegment::clone() const {
	return SmartSegmentOP( new SmartSegment( *this ) );
}
// this is the beginning of all the standard getters and setters
void
SmartSegment::set_segment_id(core::Size new_segment_id) {
	segment_id_ = new_segment_id;
}

core::Size
SmartSegment::get_segment_id() const {
	return segment_id_;
}

void
SmartSegment::set_n_terminal_neighbor(SmartSegmentOP new_n_terminal_neighbor) {
	n_terminal_neighbor_ = new_n_terminal_neighbor;
}

SmartSegmentOP
SmartSegment::get_n_terminal_neighbor() const {
	return n_terminal_neighbor_;
}

void
SmartSegment::set_c_terminal_neighbor(SmartSegmentOP new_c_terminal_neighbor) {
	c_terminal_neighbor_ = new_c_terminal_neighbor;
}

SmartSegmentOP
SmartSegment::get_c_terminal_neighbor() const {
	return c_terminal_neighbor_;
}

void
SmartSegment::set_dssp_code(char new_dssp_code){
	dssp_code_ = new_dssp_code;
}

char
SmartSegment::get_dssp_code() const{
	return dssp_code_;
}

void
SmartSegment::set_chimaeric_status(bool is_chimaera) {
	is_a_chimaera_ = is_chimaera;
}

bool
SmartSegment::is_chimaeric() const {
	return is_a_chimaera_;
}

void
SmartSegment::set_n_terminal_parent(SmartSegmentOP new_parent){
	n_terminal_parent_ = new_parent;
}

SmartSegmentOP
SmartSegment::get_n_terminal_parent() const { // a getter that fails gracefully if it can't get.
	if ( this->is_chimaeric() ) {
		return n_terminal_parent_;
	} else {
		TR << "Segment " << this->get_segment_id() << " is not chimaeric!" << std::endl;
		utility_exit_with_message("Tried to access parents of nonchimaeric segment");
	}
}

SmartSegmentOP
SmartSegment::get_far_n_terminal_parent( SmartSegmentOP current_segment ) { //will return itself if segment has no parent.
	if ( current_segment->is_chimaeric() ) {
		return get_far_n_terminal_parent( current_segment->get_n_terminal_parent() );
	} else {
		return current_segment;
	}
}


void
SmartSegment::set_c_terminal_parent(SmartSegmentOP new_parent){
	c_terminal_parent_ = new_parent;
}

SmartSegmentOP
SmartSegment::get_c_terminal_parent() const {
	if ( this->is_chimaeric() ) {
		return c_terminal_parent_;
	} else {
		TR << "Segment " << this->get_segment_id() << " is not chimaeric!" << std::endl;
		utility_exit_with_message("Tried to access parents of nonchimaeric segment");
	}
}

SmartSegmentOP
SmartSegment::get_far_c_terminal_parent( SmartSegmentOP current_segment ) { //will return itself if segment has no parent.
	if ( current_segment->is_chimaeric() ) {
		return get_far_c_terminal_parent( current_segment->get_c_terminal_parent() );
	} else {
		return current_segment;
	}
}

bool
SmartSegment::is_vital() const { // this one has no setter. We know whether things are vital in the beginning.
	return is_vital_;
}

void
SmartSegment::set_is_vital(bool is_vital)  { // this one has no setter. We know whether things are vital in the beginning.
	is_vital_ = is_vital;
}

bool
SmartSegment::is_in_Assembly() const {
	return is_in_Assembly_;
}

void
SmartSegment::set_is_in_Assembly(bool is_in_assembly) {
	is_in_Assembly_ = is_in_assembly;
}

void
SmartSegment::set_length(core::Size new_length){
	length_ = new_length;
}

core::Size
SmartSegment::get_length() const {
	// if (length_ > 0){
	// return length_;
	//} else {
	return residues_.size();
	//}

}
void
SmartSegment::set_basis_pair( BasisPair new_basis_pair ){
	//Apparently this just isn't working?
	basis_pair_ = std::make_pair( Basis( new_basis_pair.first ), Basis( new_basis_pair.second ) );
	//basis_pair_ = new_basis_pair;
}
//BasisPair
//SmartSegment::get_basis_pair(){
// return basis_pair_;
//}

BasisPair
SmartSegment::get_basis_pair() const {
	return basis_pair_;
}

// end of standard getters/setters
bool
SmartSegment::is_n_terminus_fixed() const {
	return (n_terminal_neighbor_ != nullptr);
}

bool
SmartSegment::is_c_terminus_fixed() const {
	return (c_terminal_neighbor_ != nullptr);
}

bool
SmartSegment::is_hashable() const { // this is purely shorthand to check if a segment is fixable
	return !(this->is_n_terminus_fixed() && this->is_c_terminus_fixed());
}
void
SmartSegment::set_residue_vector(utility::vector1<SmartSewingResidueOP> new_residue_vector){ // this should only be used when making the master list
	residues_ = new_residue_vector;
	length_ = new_residue_vector.size();
}

utility::vector1<SmartSewingResidueOP> &
SmartSegment::get_residue_vector() {
	return residues_;
}


utility::vector1<SmartSewingResidueOP> const &
SmartSegment::get_const_residue_vector() const {
	return residues_;
}

SmartSewingResidueOP
SmartSegment::get_residue(core::Size resnum) { // a slightly less cringe-inducing pass by value.
	if ( resnum > residues_.size() ) {
		TR << "resnum " << resnum << " is outside of a segment " << residues_.size() << " residues long!" << std::endl;
	}
	return residues_[resnum];
}


void
SmartSegment::link_to(SmartSegmentOP n_terminal_segment, SmartSegmentOP c_terminal_segment){
	//TR << "Linking C terminus of segment " << n_terminal_segment->get_segment_id() << " to N terminus of segment " << c_terminal_segment->get_segment_id() << std::endl;
	if ( n_terminal_segment && c_terminal_segment ) {
		n_terminal_segment->set_c_terminal_neighbor(c_terminal_segment);
		c_terminal_segment->set_n_terminal_neighbor(n_terminal_segment);
	}
}

void
SmartSegment::isolate(){
	if ( this->is_n_terminus_fixed() ) {
		this->get_n_terminal_neighbor()->set_c_terminal_neighbor(nullptr);
	}
	this->set_n_terminal_neighbor(nullptr);
	if ( this->is_c_terminus_fixed() ) {
		this->get_c_terminal_neighbor()->set_n_terminal_neighbor(nullptr);
	}
	this->set_c_terminal_neighbor(nullptr);
}

SmartSegmentOP
SmartSegment::get_n_most_segment(SmartSegmentOP start_segment, bool cross_chimaerae){
	if ( not(start_segment->is_n_terminus_fixed()) || (start_segment->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment;
	} else if ( not(start_segment->get_n_terminal_neighbor()->is_n_terminus_fixed()) || (start_segment->get_n_terminal_neighbor()->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment->get_n_terminal_neighbor();
	} else if ( not(start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor()->is_n_terminus_fixed()) || (start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor()->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor();
	} else {
		return get_n_most_segment(start_segment->get_n_terminal_neighbor() , cross_chimaerae);
	}
}
SmartSegmentCOP
SmartSegment::get_n_most_segment(SmartSegmentCOP start_segment, bool cross_chimaerae){

	if ( not(start_segment->is_n_terminus_fixed()) || (start_segment->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment;
	} else if ( not(start_segment->get_n_terminal_neighbor()->is_n_terminus_fixed()) || (start_segment->get_n_terminal_neighbor()->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment->get_n_terminal_neighbor();
	} else if ( not(start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor()->is_n_terminus_fixed()) || (start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor()->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment->get_n_terminal_neighbor()->get_n_terminal_neighbor();
	} else {
		return get_n_most_segment(start_segment->get_n_terminal_neighbor() , cross_chimaerae);
	}
	//Initial check to see if this is the n-most segment
	/* if( n_terminal_neighbor_ == NULLPTR //If this segment is N-terminal
	|| ( !cross_chimaerae && is_a_chimaera_ ) )
	{
	//We want to return a SmartSegmentOP to the current segment
	return SmartSegmentOP ( this ); //Not sure if this is allowed
	}*/
	//while(current_segment->is_n_terminus_fixed()){ // while we know what segment comes next
	// TR << "Moving Nwise one segment" << std::endl;
	// if(!cross_chimaerae && current_segment->is_chimaeric()){
	//  return current_segment;
	// }
	// current_segment = current_segment->get_n_terminal_neighbor();
	//}
	//return current_segment;
}

SmartSegmentCOP
SmartSegment::get_c_most_segment(SmartSegmentCOP start_segment, bool cross_chimaerae){
	if ( not(start_segment->is_c_terminus_fixed()) || (start_segment->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment;
	} else {
		return get_c_most_segment(start_segment->get_c_terminal_neighbor() , cross_chimaerae);
	}
}


SmartSegmentOP
SmartSegment::get_c_most_segment(SmartSegmentOP start_segment, bool cross_chimaerae){
	if ( not(start_segment->is_c_terminus_fixed()) || (start_segment->is_chimaeric() && !cross_chimaerae) ) {
		return start_segment;
	} else {
		return get_c_most_segment(start_segment->get_c_terminal_neighbor() , cross_chimaerae);
	}
}

std::set< core::Size >
SmartSegment::get_vital_residues() const
{
	return vital_residues_;
}


void
SmartSegment::set_vital_residues( std::set< core::Size > vital ){
	vital_residues_ = vital;
}

void
SmartSegment::add_vital_residue( core::Size resnum ){
	vital_residues_.insert( resnum );
}


bool
SmartSegment::residue_is_vital( core::Size resnum ){
	return ( vital_residues_.count( resnum ) != 0 );
}


void
SmartSegment::become( SmartSegmentOP changing_segment, SmartSegmentCOP src_segment ){

	changing_segment->set_segment_id( src_segment->get_segment_id() );
	changing_segment->set_n_terminal_neighbor( nullptr ); //This is an N-terminal segment, so it has no N-terminal neighbor
	changing_segment->set_c_terminal_neighbor( nullptr ); //Change this later with link_to, just did this to avoid memory leaks down the line
	changing_segment->set_is_vital( src_segment->is_vital() );
	changing_segment->set_length( src_segment->get_length() );
	changing_segment->set_dssp_code( src_segment->get_dssp_code() );
	//Deep copy residue vector
	utility::vector1< data_storage::SmartSewingResidueOP > const src_vec = src_segment->get_const_residue_vector();
	for ( core::Size i = 1; i <= changing_segment->get_residue_vector().size(); ++i ) {
		//If this is a residue in src_segment, copy it
		if ( i <= src_segment->get_length() ) {
			changing_segment->get_residue_vector().at( i )->become( src_segment->get_const_residue_vector().at( i )  );
		} else {
			changing_segment->get_residue_vector().at( i )->become( nullptr );
		}
	}
}
void
SmartSegment::set_const_reference_segment( SmartSegmentCOP ref ){
	const_reference_segment_ = ref;

}

SmartSegmentCOP
SmartSegment::get_const_reference_segment() const{
	return const_reference_segment_;
}


} //protocols
} //sewing
} //data_storage





