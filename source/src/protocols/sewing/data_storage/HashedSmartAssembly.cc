// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/HashedSmartAssembly.cc
/// @brief an Assembly that makes use of the Hasher
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/data_storage/HashedSmartAssembly.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sewing.data_storage.HashedSmartAssembly" );


namespace protocols {
namespace sewing {
namespace data_storage {

HashedSmartAssembly::HashedSmartAssembly():
	SmartAssembly()
{

}
HashedSmartAssembly::HashedSmartAssembly(hashing::SegmentVectorCOP segment_vector):
	SmartAssembly()
{
	this->set_segment_vector(segment_vector);
}
HashedSmartAssembly::~HashedSmartAssembly(){}

HashedSmartAssembly::HashedSmartAssembly( HashedSmartAssembly const & other):
	SmartAssembly( other )
{

}



HashedSmartAssemblyOP
HashedSmartAssembly::clone() const {
	return HashedSmartAssemblyOP( new HashedSmartAssembly( *this ) );
}

void
HashedSmartAssembly::set_basis_map_generator(hashing::BasisMapGeneratorOP new_bmg){
	existing_alignments_ = new_bmg;
}

//bool
//HashedSmartAssembly::add_segment(bool n_terminus){
bool
HashedSmartAssembly::add_segment(bool n_terminus, core::Size seg_ID_to_add, core::Size res_ID_1, core::Size res_ID_2 ){
	//runtime_assert( n_terminal_segments_.size() > 0 && c_terminal_segments_.size() > 0 );
	TR << "Adding segment." << std::endl;
	if ( n_terminus ) {
		TR << "adding to n terminus" << std::endl;
		this->set_segID_1( this->get_n_terminal_segment()->get_segment_id() );
		runtime_assert( local_segments().count( this->get_segID_1() ) != 0 );
		core::Size num_attempted = 0;
		bool chimerization_succeeded = false;
		while ( !chimerization_succeeded ) {
			++num_attempted;
			if ( num_attempted > 1000 ) {
				return false;
			}
			SmartSegmentCOP seg1;
			SmartSegmentCOP seg2;
			//Seg1 should already be in the assembly!
			seg1 = local_segments().at( this->get_segID_1() )->get_const_reference_segment();

			//DOUBLE CHIMAERA SPECIAL CASE!!
			//TR << "seg1: " << segID_1_ << " is chimaera? " << seg1->is_chimaeric() << std::endl;
			if ( seg1->is_chimaeric() ) {
				//TR << "seg1 is a chimaera" << std::endl;
				core::Size basis_res;
				if ( seg1->get_basis_pair().first.segment_id() == seg1->get_n_terminal_parent()->get_segment_id() ) {
					//first basis is n_terminal
					basis_res = seg1->get_basis_pair().first.resnum();
				} else {
					basis_res = seg1->get_basis_pair().second.resnum();
				}
				if ( !( this->get_window_width() < basis_res ) ) {
					TR << "Cannot add to n term" << std::endl;
					return false;
				}
			}
			BasisPair basis_pair;
			if ( seg_ID_to_add == 0 ) {
				//basis_pair = existing_alignments_->get_basis_pair_for_segment(segID_1_);
				//Since we're adding to the N terminus, this will be a C-terminal segment
				core::Size num_tries = 0;
				while ( !seg_ID_to_add ) {
					if ( num_tries >= 1000 ) {
						return false;
					}
					seg_ID_to_add = get_random_segment_id( !n_terminus );
					std::pair< bool, bool > good_segs = existing_alignments_->segs_can_chimerize( seg1, this->get_segment_vector()->at( seg_ID_to_add ), !n_terminus); //If n terminus is TRUE, the first segment is C terminal, so the bool should be false
					if ( !good_segs.first ) {
						return false;
					}
					if ( !good_segs.second ) {
						seg_ID_to_add = 0;
					}
					++num_tries;
				}
				this->set_segID_2( seg_ID_to_add );
			} else {
				this->set_segID_2( seg_ID_to_add );
			}

			if ( local_segments().count( this->get_segID_2() ) && local_segments().at( this->get_segID_2() )->is_in_Assembly() ) {
				TR <<"Selected segment is already in the assembly"<<std::endl;
				TR <<"Making a copy of the segment and all attached segments" << std::endl;
				core::Size new_segment_id;
				if ( pdb_segments().size() > 0 ) {
					new_segment_id = std::max( ( get_segment_vector()->size() + pdb_segments().rbegin()->first ), local_segments().rbegin()->first ) + 1;
				} else {
					new_segment_id = std::max( get_segment_vector()->size(), local_segments().rbegin()->first ) + 1;
				}
				runtime_assert( local_segments().count( new_segment_id ) == 0 );
				data_storage::SmartSegmentCOP reference_segment = local_segments().at( this->get_segID_2() )->get_const_reference_segment();


				add_segment_and_neighbors_to_local_segments( reference_segment, new_segment_id );
				//Set our segID_2_
				core::Size current_segid = new_segment_id;
				SmartSegmentOP local_seg = local_segments().at( new_segment_id );
				while ( local_seg != nullptr ) {
					if ( reference_segment->get_segment_id() == this->get_segID_2() ) {
						this->set_segID_2( current_segid );
					}
					local_seg = local_seg->get_c_terminal_neighbor();
					reference_segment = reference_segment->get_c_terminal_neighbor();
					++current_segid;
				}
				TR << this->get_segID_2() << std::endl;
				seg2 = local_segments().at( this->get_segID_2() ); //TODO Make sure we use seg2's const reference segment when looking for basis pairs
			} else { //Not in local segments
				if ( local_segments().count( this->get_segID_2() ) ) {
					seg2 = local_segments().at( this->get_segID_2() )->get_const_reference_segment();
				} else if ( this->get_segID_2() <= this->get_segment_vector()->size() ) {
					seg2 = this->get_segment_vector()->at( this->get_segID_2() );
				} else {
					seg2 = pdb_segments().at( this->get_segID_2() );
				}
				add_segment_and_neighbors_to_local_segments( seg2 );
				seg2 = local_segments().at( this->get_segID_2() );
			}
			//Check if these segments can even hash or not
			std::pair< bool, bool > good_segs = existing_alignments_->segs_can_chimerize( seg1, seg2, !n_terminus );
			if ( !good_segs.first ) {
				//Seg1 just can't chimerize on this terminus
				return false;
			}
			if ( !good_segs.second ) {
				//Seg1 can't chimerize specifically with seg2
				continue;
			}
			//NOW FIND THE BASIS PAIR
			if ( res_ID_1 && res_ID_2 ) {
				//We already set a specific basis pair to use
				if ( seg1->is_chimaeric() && !( res_ID_1 <= seg1->get_basis_pair().first.resnum() ) ) {
					utility_exit_with_message( "provided basis pair is not valid: it will chimerize over original segment." );// ya done messed up
				}
				basis_pair.first = Basis( this->get_segID_1(), res_ID_1  );
				basis_pair.second = Basis( this->get_segID_2(), res_ID_2 );
				//Is this basis pair a valid one?
				//if n_terminus, then seg1 is actually C-terminal
				if ( ! existing_alignments_->basis_pair_is_valid( basis_pair, !n_terminus ) ) {
					//TODO should this really be a utility exit or just a failed move? I guess it would mess up switching . . .
					utility_exit_with_message( "provided basis pair is not valid" ); // ya done messed up
				}
				this->set_resID_1( res_ID_1 );
				this->set_resID_2( res_ID_2 );
			} else {
				/*
				else if( basis_pair.first.segment_id() == this->get_segID_1() ){
				//We already found a basis pair earlier!
				this->set_resID_1( basis_pair.first.resnum() );
				this->set_resID_2( basis_pair.second.resnum() );
				}*/
				//The function supplied the segment ID but did not supply the basis residues
				//Get those basis residues from this BMG function that I WILL WRITE but haven't written yet
				std::pair< bool, data_storage::BasisPair > results = existing_alignments_->get_basis_pair_for_local_segments( seg1, seg2, !n_terminus );
				if ( results.first == false ) {
					return false;
				}
				basis_pair = results.second;
				if ( basis_pair.second.segment_id() != 0 ) {
					basis_pair.second.segment_id( this->get_segID_2() );
				}
				this->set_resID_1( basis_pair.first.resnum() );
				this->set_resID_2( basis_pair.second.resnum() );
			}
			if ( basis_pair.second.segment_id() == 0 ) {
				//Seg1 and Seg2 can't form a valid basis pair
				continue;
			}
			chimerization_succeeded = chimerize( basis_pair, n_terminus );
		}

		this->set_length(this->get_length() - this->get_n_terminal_segment()->get_length());
		this->set_size(this->get_size()-1);
		//this->get_n_terminal_segment()->set_is_in_Assembly( false ); //Apparently this was a bad idea? "having this set to false will cause infinitely looping assemblies"
		//This should replace the old N-terminal segment with a chimaera
		if ( this->get_c_terminal_segment() == this->get_n_terminal_segment() ) {
			//The assembly is only one segment. One chimaeric segment.
			this->set_c_terminal_segment( this->get_last_chimaera() );
		}
		if ( this->get_n_terminal_segment()->get_c_terminal_neighbor() ) {
			//The old n_terminal_segment ( shown ) needs to be replaced with its chimaera. If it has a neighbor, we can just do the old switcheroo like so.
			this->set_n_terminal_segment( this->get_n_terminal_segment()->get_c_terminal_neighbor()->get_n_terminal_neighbor() );
		} else {
			//Otherwise it's a chimaera
			this->set_n_terminal_segment( this->get_last_chimaera() );
		}
		runtime_assert( this->get_n_terminal_segment()->is_chimaeric() );
		this->get_n_terminal_segment()->get_n_terminal_parent()->set_is_in_Assembly( true );
		//this->set_last_chimaera(this->get_n_terminal_segment()); //This has been moved to chimerize

		data_storage::SmartSegmentOP first_segment = this->get_n_terminal_segment();
		while ( this->get_n_terminal_segment()->get_n_terminal_neighbor() ) {
			TR << "Adding segment " <<this->get_n_terminal_segment()->get_segment_id() << std::endl;
			this->get_n_terminal_segment()->set_is_in_Assembly(true);
			this->set_length(this->get_length() + this->get_n_terminal_segment()->get_length());
			this->set_size(this->get_size()+1);
			this->set_n_terminal_segment(this->get_n_terminal_segment()->get_n_terminal_neighbor());
			if ( this->get_n_terminal_segment()->get_segment_id() == first_segment->get_segment_id() ) {
				TR << "segment " << this->get_n_terminal_segment()->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
				utility_exit_with_message("infinite loop!");
			}
		}
		TR << "Adding segment " << this->get_n_terminal_segment()->get_segment_id() << std::endl;
		this->get_n_terminal_segment()->set_is_in_Assembly(true);
		this->set_length(this->get_length() + this->get_n_terminal_segment()->get_length());
		this->set_size(this->get_size()+1);
		//TR <<" N-terminal segment is now " << n_terminal_segment_->get_segment_id();
	} else { //END IF N TERMINUS //IF C TERMINUS
		TR << "adding to n terminus" << std::endl;
		this->set_segID_1( this->get_c_terminal_segment()->get_segment_id() );
		bool chimerization_succeeded = false;
		core::Size num_attempted = 0;
		while ( !chimerization_succeeded ) {
			++num_attempted;
			if ( num_attempted > 1000 ) {
				return false;
			}
			SmartSegmentCOP seg1;
			SmartSegmentCOP seg2;
			seg1 = local_segments().at( this->get_segID_1() )->get_const_reference_segment();

			//DOUBLE CHIMAERA SPECIAL CASE!!
			//TR << "seg1: " << segID_1_ << " is chimaera? " << seg1->is_chimaeric() << std::endl;
			if ( seg1->is_chimaeric() ) {
				//TR << "seg1 is a chimaera" << std::endl;
				core::Size basis_res;
				if ( seg1->get_basis_pair().first.segment_id() == seg1->get_n_terminal_parent()->get_segment_id() ) {
					//first basis is n_terminal
					basis_res = seg1->get_basis_pair().first.resnum();
				} else {
					basis_res = seg1->get_basis_pair().second.resnum();
				}
				if ( !( this->get_window_width() < seg1->get_length() - basis_res ) ) {
					//if( !( window_width_ < basis_res ) ){
					TR << "Cannot add to c term" << std::endl;
					return false;
				}
			}

			BasisPair basis_pair;
			if ( seg_ID_to_add == 0 ) {
				//basis_pair = existing_alignments_->get_basis_pair_for_segment(segID_1_);
				//Since we're adding to the C terminus, this will be an N-terminal segment
				core::Size num_tries = 0;
				while ( !seg_ID_to_add ) {
					if ( num_tries >= 1000 ) {
						return false;
					}
					seg_ID_to_add = get_random_segment_id( !n_terminus );
					std::pair< bool, bool > good_segs = existing_alignments_->segs_can_chimerize( seg1, this->get_segment_vector()->at( seg_ID_to_add ), !n_terminus); //If n terminus is TRUE, the first segment is C terminal, so the bool should be false
					if ( !good_segs.first ) {
						return false;
					}
					if ( !good_segs.second ) {
						seg_ID_to_add = 0;
					}
					++num_tries;
				}
				this->set_segID_2( seg_ID_to_add );
			} else {
				this->set_segID_2( seg_ID_to_add );
			}


			if ( local_segments().count( this->get_segID_2() ) && local_segments().at( this->get_segID_2() )->is_in_Assembly() ) {
				TR <<"Selected segment is already in the assembly"<<std::endl;
				TR <<"Making a copy of the segment and all attached segments" << std::endl;
				core::Size new_segment_id;
				if ( pdb_segments().size() > 0 ) {
					new_segment_id = std::max( ( get_segment_vector()->size() + pdb_segments().rbegin()->first ), local_segments().rbegin()->first ) + 1;
				} else {
					new_segment_id = std::max( get_segment_vector()->size(), local_segments().rbegin()->first ) + 1;
				}
				runtime_assert( local_segments().count( new_segment_id ) == 0 );
				data_storage::SmartSegmentCOP reference_segment = local_segments().at( this->get_segID_2() )->get_const_reference_segment();
				add_segment_and_neighbors_to_local_segments( reference_segment, new_segment_id );
				//Set our segID_2_
				core::Size current_segid = new_segment_id;
				SmartSegmentOP local_seg = local_segments().at( new_segment_id );
				while ( local_seg != nullptr ) {
					if ( reference_segment->get_segment_id() == this->get_segID_2() ) {
						this->set_segID_2( current_segid );
					}
					local_seg = local_seg->get_c_terminal_neighbor();
					reference_segment = reference_segment->get_c_terminal_neighbor();
					++current_segid;
				}
				TR << this->get_segID_2() << std::endl;
				seg2 = local_segments().at( this->get_segID_2() ); //TODO Make sure we use seg2's const reference segment when looking for basis pairs
			} else { //Not in local segments
				if ( local_segments().count( this->get_segID_2() ) ) {
					seg2 = local_segments().at( this->get_segID_2() )->get_const_reference_segment();
				} else if ( this->get_segID_2() <= get_segment_vector()->size() ) {
					seg2 = get_segment_vector()->at( this->get_segID_2() );
				} else {
					seg2 = pdb_segments().at( this->get_segID_2() );
				}
				add_segment_and_neighbors_to_local_segments( seg2 );
				seg2 = local_segments().at( this->get_segID_2() );
			}

			//Check if these segments can even hash or not
			std::pair< bool, bool > good_segs = existing_alignments_->segs_can_chimerize( seg1, seg2, !n_terminus );
			if ( !good_segs.first ) {
				//Seg1 just can't chimerize on this terminus
				return false;
			}
			if ( !good_segs.second ) {
				//Seg1 can't chimerize specifically with seg2
				continue;
			}


			//NOW FIND THE BASIS PAIR
			if ( res_ID_1 && res_ID_2 ) {
				//We already set the basis pair
				// seg1 is chimaeric and C terminal--we want to make sure we save at least part of the C-terminal half
				if ( seg1->is_chimaeric() && !( res_ID_1 >= seg1->get_basis_pair().first.resnum() ) ) {
					utility_exit_with_message( "provided basis pair is not valid: it will chimerize over original segment." );// ya done messed up
				}
				basis_pair.first = Basis( this->get_segID_1(), res_ID_1  );
				basis_pair.second = Basis( this->get_segID_2(), res_ID_2 );
				//Is this basis pair a valid one?
				if ( ! existing_alignments_->basis_pair_is_valid( basis_pair, !n_terminus ) ) {
					//TODO should this really be a utility exit or just a failed move? I guess it would mess up switching . . .
					utility_exit_with_message( "provided basis pair is not valid" ); // ya done messed up
				}
				this->set_resID_1( res_ID_1 );
				this->set_resID_2( res_ID_2 );
			} else {
				/*
				else if( basis_pair.first.segment_id() == this->get_segID_1() ){
				//We already found a basis pair earlier!
				this->set_resID_1(  basis_pair.first.resnum() );
				this->set_resID_2( basis_pair.second.resnum() );
				}
				*/
				//The function apparently supplied the segment ID but did not supply the basis residues
				//Get those basis residues from this BMG function that I WILL WRITE but haven't written yet
				std::pair< bool, data_storage::BasisPair > results = existing_alignments_->get_basis_pair_for_local_segments( seg1, seg2, !n_terminus );
				if ( results.first == false ) {
					return false;
				}
				basis_pair = results.second;
				if ( basis_pair.second.segment_id() != 0 ) {
					basis_pair.second.segment_id( this->get_segID_2() );
				}
				this->set_resID_1( basis_pair.first.resnum() );
				this->set_resID_2( basis_pair.second.resnum() );
			}
			if ( basis_pair.second.segment_id() == 0 ) {
				//Seg1 and Seg2 can't form a valid chimaera
				continue;
			}

			chimerization_succeeded = chimerize(basis_pair, n_terminus);
		}
		this->set_length(this->get_length() - this->get_c_terminal_segment()->get_length());
		this->set_size(this->get_size()-1);
		//this->get_c_terminal_segment()->set_is_in_Assembly( false ); //Apparently a bad idea for some reason

		if ( this->get_c_terminal_segment() == this->get_n_terminal_segment() ) {
			//The assembly is only one segment. One chimaeric segment.
			this->set_n_terminal_segment( this->get_last_chimaera() );
		}
		if ( this->get_c_terminal_segment()->get_n_terminal_neighbor() ) {
			this->set_c_terminal_segment( this->get_c_terminal_segment()->get_n_terminal_neighbor()->get_c_terminal_neighbor() );
		} else {
			this->set_c_terminal_segment( this->get_last_chimaera() );
		}
		runtime_assert( this->get_c_terminal_segment()->is_chimaeric() );
		this->get_c_terminal_segment()->get_c_terminal_parent()->set_is_in_Assembly( true );


		data_storage::SmartSegmentOP first_segment = this->get_c_terminal_segment();
		//first_segment_ = this->get_c_terminal_segment();
		while ( this->get_c_terminal_segment()->get_c_terminal_neighbor() ) {
			TR << "Adding segment " <<this->get_c_terminal_segment()->get_segment_id() << std::endl;
			this->get_c_terminal_segment()->set_is_in_Assembly(true);
			this->set_length(this->get_length() + this->get_c_terminal_segment()->get_length());
			this->set_size(this->get_size()+1);
			this->set_c_terminal_segment(this->get_c_terminal_segment()->get_c_terminal_neighbor());
			if ( this->get_c_terminal_segment()->get_segment_id() == first_segment->get_segment_id() ) {
				TR << "segment " << this->get_c_terminal_segment()->get_segment_id() << " links n-terminally to n terminal segment" << std::endl;
				utility_exit_with_message("infinite loop!");
			}
		}
		TR << "Adding segment " << this->get_c_terminal_segment()->get_segment_id() << std::endl;
		this->get_c_terminal_segment()->set_is_in_Assembly(true);
		this->set_length(this->get_length() + this->get_c_terminal_segment()->get_length());
		this->set_size(this->get_size()+1);
	} //END IF C TERMINUS
	this->set_last_change_was_n_terminal(n_terminus);
	this->set_last_change('A');


	//Cleanup taken from SmartAssembly. Not entirely sure why this is necessary, but sure.
	if ( this->get_n_terminal_segment() != this->get_c_terminal_segment() ) {
		if ( this->get_n_terminal_segment()->is_vital() ) {
			this->set_n_terminal_segment( this->get_c_terminal_segment() );
			while ( this->get_n_terminal_segment()->is_n_terminus_fixed() ) {
				this->set_n_terminal_segment( this->get_n_terminal_segment()->get_n_terminal_neighbor() );
			}
		} else if ( this->get_c_terminal_segment()->is_vital() ) {
			this->set_c_terminal_segment( this->get_n_terminal_segment() );
			while ( this->get_c_terminal_segment()->is_c_terminus_fixed() ) {
				this->set_c_terminal_segment( this->get_c_terminal_segment()->get_c_terminal_neighbor() );
			}
		}
	}
	return true;
}

core::Size
HashedSmartAssembly::get_window_width() const{
	return existing_alignments_->hasher_settings().min_hash_score / 4;
}

} //protocols
} //sewing
} //data_storage






