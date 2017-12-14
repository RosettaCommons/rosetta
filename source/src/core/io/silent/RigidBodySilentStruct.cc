// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/ProteinSilentStruct.tmpl.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author Oliver Lange

// Unit headers
#include <core/io/silent/RigidBodySilentStruct.hh>

#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
// C++ Headers

#include <ObjexxFCL/FArray2D.hh>

// option key includes

namespace core {
namespace io {
namespace silent {

using namespace ObjexxFCL;
//using namespace ObjexxFCL::format;


static basic::Tracer tr( "core.io.silent" );

void
RigidBodySilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	decoy_tag( tag );
	if ( tag == "empty_tag" ) set_tag_from_pose( pose );
	energies_from_pose( pose );

	// conformation information
	sequence( pose.annotated_sequence( true /* show-all-variants */ ) );
	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		jumps_.clear();
		for ( Size nr = 1; nr <= pose.num_jump(); nr++ )  {
			add_jump( pose.jump(nr) );
		}
	}
	fold_tree( pose.fold_tree() );
}


bool RigidBodySilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	bool success( false );

	utility::vector1< std::string > energy_names;
	auto iter = lines.begin();
	if ( iter->substr(0,9) != "SEQUENCE:" ) { //a full new header would change the default columns
		// get sequence and scorename data from the silent-file data object, because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			utility::pointer::static_pointer_cast< core::io::silent::EnergyNames > ( container.get_shared_silent_data( energynames ) )
		);

		SimpleSequenceDataOP seqdata = SimpleSequenceDataOP(
			utility::pointer::static_pointer_cast< core::io::silent::SimpleSequenceData > ( container.get_shared_silent_data( simplesequencedata ) )
		);

		sequence      ( seqdata->sequence()   );
		energy_names = enames ->energy_names();
	} else {
		// get sequence and scorename data from the first two lines provided, put
		// them into container for further use by other RigidBodySilentStruct
		// objects.

		// first line is SEQUENCE:
		if ( !read_sequence( *iter ) ) return false;
		++iter;
		read_score_headers( *iter, energy_names, container ); ++iter;
	} // get header information
	for ( auto end = lines.end(); iter != end; ++iter ) {
		std::string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,6) == "REMARK" ) {
			//    std::string tag;
			//    std::string comment;
			//    std::string value;
			//    runtime_assert( tag == "REMARK" );
			//    line_stream >> tag >> comment >> value;
			comment_from_line( *iter );//add_comment( comment, value );
			continue;  // don't skip comments
		}

		if ( iter->substr(0,7) == "SCORE: " ) { // SCORE: line with values from this structure.
			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}

			parse_energies( line_stream, energy_names );
		} else { // conformation lines
			// parse fold_tree and jump lines
			if ( iter->substr(0,10) == "FOLD_TREE " ) {
				kinematics::FoldTree f;
				line_stream >> f;
				fold_tree( f ); // add fold-tree to this SilentStruct
				tr.Debug << "read fold-tree " << f; //"\n" is in fold-tree output
				tr.Debug << "reading " << f.num_jump() << " jumps " << std::endl;
				continue;
			} else if ( iter->substr(0,2) == "RT" ) {
				kinematics::Jump jump;
				line_stream >> jump;
				tr.Debug << "read jump " << jump << std::endl;
				add_jump( jump );
				continue;
			} else if ( iter->substr(0,9) == "SEQUENCE:" ) {
				tr.Warning << "skipping duplicate sequence declaration " << std::endl;
				//after a SEQUENCE declaration we might find another SCORE header that should be skipped, too...
				auto iter2 = ++iter;
				if ( ( iter2 != end ) && iter2->substr(0,7) == "SCORE: " ) {
					tr.Warning << "re-reading score declaration from second line... " << std::endl;
					read_score_headers( *iter2, energy_names, container );
					++iter;
				}
				continue;
			} else if ( iter->substr(0,19) == "ANNOTATED_SEQUENCE:" ) {
				std::string annotated_seq;
				line_stream >> tag; //ANNOTATED_SEQUENCE
				line_stream >> annotated_seq;
				sequence( annotated_seq );
				continue;
			}
		} //conformation lines
	} // for ( iter ... )

	if ( fold_tree_ && fold_tree_->num_jump() != jumps_.size() ) {
		tr.Warning << "parse error:  found only " << jumps_.size()
			<< " RT-lines for a fold-tree with " << fold_tree_->num_jump()
			<< " for decoy tag " << decoy_tag() << std::endl;
		return false;
	}

	success = true;
	return success;
} // init_from_lines

void RigidBodySilentStruct::fold_tree( kinematics::FoldTree const& f ) {
	fold_tree_ = kinematics::FoldTreeOP( new kinematics::FoldTree(f) );
}

kinematics::FoldTree const& RigidBodySilentStruct::fold_tree() const {
	runtime_assert( fold_tree_ != nullptr );
	return *fold_tree_;
}

void RigidBodySilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const &
) const {
	tr.Warning << "RigidBodySilentStruct cannot regenerate the pose, and thus residue_type_set is ignored" << std::endl;
	fill_pose( pose );
}

void RigidBodySilentStruct::fill_pose(
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	tr.Debug << "fill_pose: RigidBodySilentStruct " << std::endl;

	runtime_assert( sequence() != "" );

	if ( fold_tree_ ) { // set fold_tree
		pose.fold_tree( fold_tree() );
	}

	// set jumps
	if ( pose.num_jump() != jumps_.size() ) {
		utility_exit_with_message( "RuntimeAssert failed: num_jump() == jumps in silent-struct" );
	}

	for ( Size nr = 1; nr <= pose.num_jump(); nr++ )  {
		pose.set_jump( nr, jump( nr ) );
	}

	if ( pose.size() != one_letter_sequence().length() ) {
		utility_exit_with_message( "RuntimeAssert failed: nres() == one_letter_sequence().length()" );
	}

	if ( options().keep_input_scores() ) {
		tr.Debug << "keep input scores... call energies into pose " << std::endl;
		energies_into_pose( pose );
	}
} // fill_pose

void RigidBodySilentStruct::print_conformation( std::ostream & output ) const {
	// fold tree
	// output << "REMARK RIGID_BODY_SILENTFILE" << std::endl;
	if ( write_fold_tree_ && fold_tree_ ) { //assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1
		output << "FOLD_TREE ";
		for ( auto const & it : fold_tree() ) {
			output << it;
		}
		//  output << fold_tree(); this produces a new-line --- wrong behaviour of fold_tree but I don't want to fix 1000 u-tracer unit-tests!
		output << ' ' << decoy_tag() << "\n";
	}
	for ( Size i = 1; i <= jumps_.size(); i++ ) {
		output << jumps_[ i ] << ' ' << decoy_tag() << "\n";
	}

	// output << ' ' << decoy_tag();
	// output << "\n";
} // print_conformation


RigidBodySilentStruct & RigidBodySilentStruct::operator= (
	RigidBodySilentStruct const & src
)
{
	// fold-tree and jumps
	for ( Size jj = 1; jj <= src.njumps(); ++jj ) {
		add_rt( src.jump(jj) );
	}
	if ( src.fold_tree_ ) {
		fold_tree_ = kinematics::FoldTreeOP( new kinematics::FoldTree( *src.fold_tree_ ) );
	}
	return *this;
}

ObjexxFCL::FArray2D< Real > RigidBodySilentStruct::get_CA_xyz() const {
	utility_exit_with_message( "RigidBody SilentStruct has no CA-coords" );
	ObjexxFCL::FArray2D< Real > dummy;
	return dummy;
}

} // namespace silent
} // namespace io
} // namespace core
