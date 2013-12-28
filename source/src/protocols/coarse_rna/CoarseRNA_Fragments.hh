// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_coarse_rna_CoarseRNA_Fragments_HH
#define INCLUDED_protocols_coarse_rna_CoarseRNA_Fragments_HH

#include <protocols/farna/RNA_Fragments.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>
// AUTO-REMOVED #include <protocols/farna/RNA_MatchType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <map>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Goal: to make a fragment object that can choose fragments
// "on the fly" for RNA ab inito folding.
//
// After reading in a set of torsions from, e.g., the ribosome crystal structure,
//  should be able to generate fragments of size 1, 2, or 3, with
//  exact sequence matches, partial Y/R matches, or ignoring sequence.
//
namespace protocols{
namespace coarse_rna{

	typedef std::pair< std::string, std::string > SequenceSecStructPair;

	class SourcePositions : public utility::pointer::ReferenceCount{

	public:

		SourcePositions(){}
		virtual ~SourcePositions();
		core::Size size() const { return source_positions_.size(); }
		void push_back( core::Size const & value ){ source_positions_.push_back( value ); }
		core::Size operator[]( core::Size num ){ return source_positions_[ num ]; }

	private:

		utility::vector1< core::Size > source_positions_;

	};

	typedef utility::pointer::owning_ptr< SourcePositions > SourcePositionsOP;


	/////////////////////////////////////////////////////////////////////////////////////////////////
	class CoarseRNA_Fragments : public protocols::farna::RNA_Fragments {
	public:

		//Constructor -- needs vall_torsions_file to get started.
		// CoarseRNA_Fragments();
		CoarseRNA_Fragments( std::string const & filename );
		virtual ~CoarseRNA_Fragments();

		virtual void
		apply_random_fragment(
          core::pose::Pose & pose,
					core::Size const position,
					core::Size const size,
					core::Size const type,
					protocols::toolbox::AllowInsertOP allow_insert );

		virtual bool
		is_fullatom();

	private:

		void
		insert_fragment(
										core::pose::Pose & pose,
										Size const & insert_res,
										Size const & source_res,
										Size const & frag_size,
										protocols::toolbox::AllowInsertOP allow_insert );

		void
		find_source_positions( SequenceSecStructPair const & key );

		Size
		pick_random_fragment(
												 std::string const RNA_string,
												 std::string const RNA_secstruct_string,
												 Size const type /* = MATCH_YR */);

		Size
		pick_random_fragment(
												 core::pose::Pose & pose,
												 Size const position,
												 Size const size,
												 Size const type );

		void
		initialize_frag_source_pose();

		std::string frag_source_secstruct_;

		std::string const frag_source_file_;
		core::pose::MiniPoseOP frag_source_pose_;
		std::map< SequenceSecStructPair, SourcePositionsOP > source_positions_map_;

		std::map< std::string, Size > coarse_rna_name_to_num_;
	};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////


}
}

#endif
