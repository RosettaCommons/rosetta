// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/abinitio/Template.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <core/fragment/FragSet.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/sequence/util.hh>
#include <core/sequence/DerivedSequenceMapping.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <fstream>

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.Templates" );
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace ObjexxFCL::format;

void protocols::abinitio::Template::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( templates::min_frag_size   );
	option.add_relevant( templates::max_shrink      );
	option.add_relevant( templates::shrink_step     );
	option.add_relevant( templates::shrink_pos_step );
	option.add_relevant( templates::min_padding     );
	option.add_relevant( templates::min_align_pos   );
	option.add_relevant( templates::max_align_pos   );
}

namespace protocols {
namespace abinitio {

/// @details Auto-generated virtual destructor
Template::~Template() = default;

using namespace core;
using namespace fragment;
using namespace scoring::constraints;

void
dump_movemap( kinematics::MoveMap const& mm, Size nres, std::ostream& out ) {
	for ( Size i = 1; i<=nres; i++ ) {
		if ( (i-1)%10 == 0 ) { out << i; continue; }
		//large numbers take several characters... skip appropriate
		if ( (i>=10) && (i-2)%10 == 0 ) { continue; }
		if ( (i>=100) && (i-3)%10 == 0 ) { continue; }
		if ( (i>=1000) && (i-4)%10 == 0 ) { continue; }
		out << ".";
	}
	out << std::endl;
	for ( Size i = 1; i<=nres; i++ ) {
		if ( mm.get_bb( i ) ) out << 'F'; //cuttable
		else out << '.';
	}
	out << std::endl;
}

Template::Template( std::string  name, pose::PoseCOP pose, core::sequence::DerivedSequenceMapping const& mapping )
: name_ (std::move( name ))
{
	tr.Error << "STUB ERROR: Template::Template(pose, mapping) constructure,  not really finished yet !!! " << std::endl;
	pose_ = pose;
	mapping_ = mapping;
	reverse_mapping_ = mapping;
	reverse_mapping_.reverse();
	good_ = false;
}

Template::Template(
	std::string const&  name,
	pose::PoseCOP pose,
	std::string const& map_file,
	int offset,
	Real score )
: pose_(std::move( pose )),
	name_( name ),
	score_( score )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	good_ = false;

	// read alignment
	tr.Info << "template  " << name << " read alignment file: " << map_file << std::endl;
	mapping_ = sequence::simple_mapping_from_file( map_file );

	// offset it because pdb starts at offset rather than sequence pos nr 1
	if ( offset >= 0 ) { //ignore negative offsets -- in agreement with myself
		mapping_.set_offset( offset-1 ); //offset 1 --> do nothing
	} else { //find offset automatically
		std::string pdb_seq = pose_->sequence();
		size_t found = pdb_seq.find( mapping_.seq2().substr(0,25) ); // look where query sequence starts in pdb-file
		if ( found!=std::string::npos ) { //found
			offset = mapping_.start_seq2() - found;
			mapping_.set_offset( offset-1 );
		} else {
			tr.Error << "query sequence could not be found for model " + name << std::endl;
			tr.Error << "pdb_seq:       " <<  pdb_seq << std::endl;
			tr.Error << "mapping_.seq2()" <<   mapping_.seq2() << std::endl;
			std::ofstream out( "BAD_SEQUENCES", std::ios_base::out | std::ios_base::app );
			out << name << "  "
				<< " Expecting:  "  << mapping_.seq2()
				<< " PDBseq:  " << pdb_seq << std::endl;
			return ;
		}
	}

	if ( option[ templates::min_align_pos ].user() ) {
		Size const min_align_pos( option[ templates::min_align_pos ] );
		for ( Size pos = 1; pos < min_align_pos && pos <= mapping_.size1(); ++pos ) {
			mapping_[ pos ] = 0;
		}
	}

	if ( option[ templates::max_align_pos ].user() ) {
		Size const max_align_pos( option[ templates::max_align_pos ] );
		for ( Size pos = max_align_pos+1; pos <= mapping_.size1(); ++pos ) {
			mapping_[ pos ] = 0;
		}
	}

	// get the reverse mapping: template --> target
	reverse_mapping_ = mapping_;
	reverse_mapping_.reverse();

	// get strand pairings
	strand_pairings_ = core::scoring::dssp::StrandPairingSetOP( new core::scoring::dssp::StrandPairingSet( *pose_ ) );
	tr.Info << "strand pairings \n" << *strand_pairings_;

	// for shits and giggles -- write out the pairings for the target
	core::scoring::dssp::PairingList template_pairings, target_pairings;
	strand_pairings().get_beta_pairs( template_pairings );
	map_pairings2target( template_pairings, target_pairings );
	core::scoring::dssp::StrandPairingSet target_strand_pairings( target_pairings );
	tr.Info << " strand_pairings of " << name << " aligned to target \n" << target_strand_pairings << std::endl;

	// for information purpose: write sequence where alignment starts. should be the same as seen in hhr file
	std::string seq = pose_->sequence();
	Size tpos, pos = 1; while ( (tpos = mapping_[ pos++ ])<=0 ) ; //go to first aligned residue
	tr.Info << "template sequence " << seq.substr( tpos-1 ) << std::endl; //string count from 0
	tr.Debug <<"compare with      " << mapping_.seq2() << std::endl;
	tr.Info << "first 10 aligned residues:" << std::endl;
	for ( Size i = pos-1; i<= pos + 9; i++ ) {
		tr.Info << i << " " << mapping_[ i ] << " " << ( mapping_[ i ]  ?  pose_->residue( mapping_[i] ).name1() : '_' ) << std::endl;
	}
	tr.Trace << mapping_;

	good_ = true;
}

Size
Template::pick_large_frags( FragSet& frag_set, core::fragment::SingleResidueFragDataOP srfd_type, Size ncopies ) const {
	// get vector of boolean for residues that are aligned
	// change that into vector of max_frag_length
	// [ a a _ _ a a a a a _ _ a ] --> 2 1 0 0 5 4 3 2 1 0 0 x
	typedef utility::vector1< Size > FragLengthMap;
	FrameList template_frames;
	Size const nres( std::min( pose_->size(), reverse_mapping_.size1() ) );
	FragLengthMap frag_length( nres, 0);
	for ( Size pos1 = nres; pos1 >= 2; pos1-- ) {
		//if pos1 is aligned .. assign appropriate lengths for each residue in continuous stretch, e.g., 5 4 3 2 1
		if ( reverse_mapping_[ pos1 ] ) {
			Size cont_length( 0 );
			for ( Size pos2 = pos1;
					pos2 >= 2 &&  //avoids out-of-range in reverse_mapping_[ pos2 - 1 ] ...
					reverse_mapping_[ pos2 ]; // &&
					pos2-- ) {
				cont_length++;
				frag_length[ pos2 ] = cont_length;
				if ( reverse_mapping_[ pos2 - 1 ] + 1 !=  reverse_mapping_[ pos2 ] ) break;
			}
			//forward outer loop to next unknown position
			pos1 = pos1 - cont_length + 1; //+ 1 since at end of loop there will be a pos1-- ...
		} // if pos1 aligned
	}
	if ( reverse_mapping_[ 1 ] ) {
		frag_length[ 1 ] = frag_length[ 2 ] + 1;
	}
	if ( tr.Trace.visible() ) {
		tr.Trace << "frag_lengths:\n pos length";
		for ( Size pos = 1; pos <= nres; pos++ ) {
			tr.Trace << "\n" << pos << " " << frag_length[ pos ];
		}
		tr.Trace << std::endl; //flush
		tr.Trace << "reverse mapping " << reverse_mapping_ << std::endl;
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Size const min_size( option[ templates::min_frag_size ] );
	Size const min_padding( option[ templates::min_padding ] ); //min_padding residues distance to gaps
	Size const opt_max_shrink( option[ templates::max_shrink ] );
	Size const shrink_step( option[ templates::shrink_step ]);
	Size const pos_step( option[ templates::shrink_pos_step ]);
	// for each position pick largest possible fragment and a couple of smaller ones...
	// have some option to control this.
	for ( Size pos = 1; pos<= nres; pos += std::max( (int) frag_length[ pos ], 1) ) {
		Size const max_size( frag_length[ pos ] );
		if ( max_size < min_size ) continue;
		Size const max_shrink( std::min( (int) opt_max_shrink, (int) max_size - (int) min_size ) );
		// take whole aligned region and up to max_shrink smaller fragments
		for ( Size shrink = min_padding; shrink <= max_shrink; shrink+=shrink_step ) {
			Size const size( max_size - shrink );

			// a fragment of requested size
			FragDataCOP frag_type( FragDataOP( new AnnotatedFragData( name(), pos, FragData( srfd_type, size )) ) );

			// pick this fragment for different frame positions and add to template_frames.
			for ( Size offset = min_padding; offset <= max_size - size; offset+=pos_step ) {
				FrameOP aFrame( new Frame( pos + offset, frag_type ) );
				for ( Size i = 1; i<=ncopies; i++ ) {
					aFrame->steal( *pose_ );

				}
				template_frames.push_back( aFrame );
			}
		}
	}
	tr.Trace << template_frames << std::endl;
	map2target( template_frames );
	frag_set.add( template_frames );
	return template_frames.size();
}

//@detail  pick ncopies identical fragments from template ( default ncopies = 1, higher number allow to give template frags higher weight )
// until proper weighting in fragments is set up
Size
Template::steal_frags( FrameList const& frames, FragSet &accumulator, Size ncopies ) const {
	FrameList template_frames;
	map2template( frames, template_frames );
	Size total( 0 );
	for ( auto & template_frame : template_frames ) {
		tr.Trace << name() << " pick frags at " << template_frame->start() << " " << template_frame->end() << " " << pose_->size() << std::endl;
		for ( Size ct = 1; ct <= ncopies; ct ++ ) {
			if ( template_frame->steal( *pose_ ) ) total++;
		}
	}

	//this should always work, since we created only alignable frames
	map2target( template_frames );
	accumulator.add( template_frames );
	return total;
}

bool
Template::map_pairing( core::scoring::dssp::Pairing const& in, core::scoring::dssp::Pairing &out, core::sequence::DerivedSequenceMapping const& map ) const {
	out.Pos1(map[in.Pos1()]);
	out.Pos2(map[in.Pos2()]);
	out.Orientation(in.Orientation());
	out.Pleating(in.Pleating()); // how is this different from calling out = in?
	return !( out.Pos1() <= 0 || out.Pos2() <= 0 );
}

void
Template::map_pairings2template( core::scoring::dssp::PairingsList const& in, core::scoring::dssp::PairingsList &out ) const {
	 for ( auto const & it : in ) {
		core::scoring::dssp::Pairing pairing;
		if ( map_pairing( it, pairing, mapping_ ) ) out.push_back( pairing );
	}
}

void
Template::map_pairings2target( core::scoring::dssp::PairingsList const& in, core::scoring::dssp::PairingsList &out ) const {
	 for ( auto const & it : in ) {
		core::scoring::dssp::Pairing pairing;
		if ( map_pairing( it, pairing, reverse_mapping_ ) ) out.push_back( pairing );
		tr.Trace << "template: "<<it<< "     target: " << pairing;
	}
}

void
Template::map2target( FrameList const& frames, FrameList& target_frames ) const {
	 for ( auto const & it : frames ) {
		FrameOP frame = it->clone_with_frags();
		if ( frame->align( reverse_mapping_ ) ) {
			target_frames.push_back( frame );
		}
	}
}

void
Template::map2target( FrameList &frames ) const {
	for ( FrameList::const_iterator it=frames.begin(),
			eit = frames.end(); it!=eit; ++it ) {
		FrameOP frame = (*it);
		Size start( frame->start() );
		tr.Trace << start << " " << frame->end() << std::endl;
		if ( !frame->align( reverse_mapping_ ) ) {
			std::cerr << start << " could not align frame " << *frame << std::endl;
			utility_exit_with_message("Could not align frame ");
		}
	}
}

void
Template::map2template( FrameList &frames ) const {
	for ( FrameList::const_iterator it=frames.begin(),
			eit = frames.end(); it!=eit; ++it ) {
		if ( !(*it)->align( mapping_ ) ) {
			tr.Error << "could not align frame " << **it << std::endl;
			utility_exit_with_message("Could not align frame ");
		}
	}
}

void
Template::map2template( FrameList const& target_frames, FrameList& template_frames ) const {
	 for ( auto const & target_frame : target_frames ) {
		FrameOP frame = target_frame->clone();
		*frame = *target_frame; //copy fragments
		if ( frame->align( mapping_ ) ) {
			template_frames.push_back( frame );
		}
	}
}

using namespace core::scoring::constraints;

void Template::read_constraints( std::string const& cst_file ) {
	cstfile_ = cst_file;  //lazy read
}

void Template::_read_constraints( std::string const& cst_file ) const {
	ConstraintSetOP cstset = ConstraintIO::get_instance()->read_constraints( cst_file, ConstraintSetOP( new ConstraintSet ), *pose_ );
	typedef utility::vector1< ConstraintCOP > FlatList;
	FlatList all_cst = cstset->get_all_constraints();
	for ( FlatList::const_iterator it = all_cst.begin(),
			eit = all_cst.end(); it!=eit; ++it ) {
		ConstraintCOP ptr = *it;
		AtomPairConstraintCOP ptr2 = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AtomPairConstraint const > ( ptr );
		if ( ptr2 ) {
			AtomPairConstraintOP valued_cst = utility::pointer::const_pointer_cast< AtomPairConstraint >(ptr2);
			cstset_.push_back( core::scoring::constraints::Obsolet_NamedAtomPairConstraintOP( new Obsolet_NamedAtomPairConstraint( valued_cst, *pose_ ) ) );
		} else {
			tr.Warning << "WARNING: constraint found that is not AtomPairConstraint... will be ignored by Template" << std::endl;
		}
	}
}

void
Template::map2target(
	NamedAtomPairConstraintList const& template_list,
	NamedAtomPairConstraintList& target_list ) const
{
	using namespace core::id;

	target_list.clear();
	 for ( auto const & it : template_list ) {
		Obsolet_NamedAtomPairConstraintOP new_cst = it->mapto( reverse_mapping_ );
		if ( new_cst ) {
			target_list.push_back( new_cst );
		} else {
			tr.Trace << "map2target: could not align constraint " << *it << std::endl;
		}
	}
}

void
Template::map2template(
	NamedAtomPairConstraintList const& target_list,
	NamedAtomPairConstraintList& template_list ) const
{
	using namespace core::id;
	template_list.clear();
	 for ( auto const & it : target_list ) {
		Obsolet_NamedAtomPairConstraintOP new_cst = it->mapto( mapping_ );
		if ( new_cst ) {
			template_list.push_back( new_cst );
		} else {
			tr.Trace << "map2template: could not align constraint " << *it << std::endl;
		}
	}
}

void
Template::cull_violators(
	NamedAtomPairConstraintList const& target_list,
	NamedAtomPairConstraintList& culled_list ) const
{
	using namespace core::scoring::constraints;
	 for ( auto const & it : target_list ) {
		AtomPairConstraintOP cst = it->mapto( mapping_, *pose_ );
		if ( cst ) {
			tr.Trace << "test: (target)" << cst->atom(1) << " " << cst->atom(2) <<  std::endl;
			if ( cst->score( *pose_ ) < 1.0 ) {
				culled_list.push_back( it );
			} else {
				tr.Trace <<"cull constraint with score: " << cst->score( *pose_ ) << " "
					<< cst->atom(1) << " " << cst->atom(2) << std::endl;
			} // cull this
		} else { //atoms present
			culled_list.push_back( it ); //constraint cannot be tested on template... keep it !
		}
	}
}

} //abinitio
} //protocols
