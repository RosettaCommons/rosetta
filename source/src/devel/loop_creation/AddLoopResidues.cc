// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddLoopResidues.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <devel/loop_creation/AddLoopResidues.hh>
#include <devel/loop_creation/AddLoopResiduesCreator.hh>

#include <numeric/random/random.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/jd2/util.hh>

//Basic
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

namespace devel {
namespace loop_creation {

using namespace std;
using namespace core;

static THREAD_LOCAL basic::Tracer TR( "devel.loop_creation.AddLoopResidues" );

/****CREATOR FUNCTIONS*****/
protocols::moves::MoverOP
AddLoopResiduesCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddLoopResidues );
}

std::string
AddLoopResiduesCreator::keyname() const
{
	return AddLoopResiduesCreator::mover_name();
}

std::string
AddLoopResiduesCreator::mover_name()
{
	return "AddLoopResidues";
}
/****END CREATOR FUNCTIONS*****/

AddLoopResidues::AddLoopResidues():
	Mover("AddLoopResidues")
{}

protocols::moves::MoverOP
AddLoopResidues::clone() const {
	return( protocols::moves::MoverOP( new AddLoopResidues( *this ) ) );
}
protocols::moves::MoverOP
AddLoopResidues::fresh_instance() const {
	return protocols::moves::MoverOP( new AddLoopResidues );
}

string
AddLoopResidues::get_name() const {
	return "AddLoopResidues";
}

bool
AddLoopResidues::reinitialize_for_new_input() const { return true; }


void
AddLoopResidues::apply( pose::Pose & pose ) {
	utility::vector1<core::Size> saved_anchors = loop_anchors_;
	core::Size saved_asym_size = asym_size_;

	protocols::loops::Loops all_loops;
	protocols::loops::Loops duplicate_loops;
	core::chemical::ResidueTypeSetCOP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::ResidueTypeCOP alanine_type = restype_set->get_representative_type_aa(core::chemical::aa_from_oneletter_code('A'));
	for ( core::Size anchor_index=1; anchor_index<=loop_anchors_.size(); ++anchor_index ) {
		//Pick a random size
		core::Size random_size_index = numeric::random::random_range(1, loop_sizes_.size());
		core::Size loop_size = loop_sizes_[random_size_index];

		TR << "Appending " << loop_size << " residues to residue: " << loop_anchors_[anchor_index] << std::endl;
		core::Size cur_anchor = loop_anchors_[anchor_index];
		for ( core::Size i=1; i<=loop_size; ++i ) {
			conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *alanine_type );
			pose.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, cur_anchor, true);
			++cur_anchor;
		}

		protocols::loops::Loop new_loop(loop_anchors_[anchor_index]+1, loop_anchors_[anchor_index]+loop_size);

		//for the loops file, extend the loop by one residue in both directions
		all_loops.push_back(protocols::loops::Loop(new_loop.start()-1, new_loop.stop()+1));

		asym_size_ += loop_size;
		core::Size dup_anchor = loop_anchors_[anchor_index] + asym_size_;
		cur_anchor = dup_anchor;

		protocols::loops::Loops just_duplicated;
		while ( dup_anchor < pose.size() ) {
			TR << "Appending " << loop_size << " residues to residue: " << dup_anchor << std::endl;
			for ( core::Size i=1; i<=loop_size; ++i ) {
				conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *alanine_type );
				pose.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, cur_anchor, true);
				++cur_anchor;
			}
			just_duplicated.push_back(protocols::loops::Loop(dup_anchor, dup_anchor+loop_size+1));
			dup_anchor+=asym_size_;
		}
		update_anchors(loop_anchors_, new_loop, anchor_index);
		duplicate_loops = update_loops(new_loop, duplicate_loops);
		for ( core::Size i=1; i<=just_duplicated.size(); ++i ) {
			duplicate_loops.push_back(just_duplicated[i]);
		}
	}

	for ( core::Size i=1; i<=duplicate_loops.size(); ++i ) {
		all_loops.push_back(duplicate_loops[i]);
	}

	std::string const job_name ( protocols::jd2::current_output_name() );
	dump_loops_file(job_name + ".loops", all_loops);

	//revert back to saved anchors.
	loop_anchors_ = saved_anchors;
	asym_size_ = saved_asym_size;
}

protocols::loops::Loops
AddLoopResidues::update_loops(
	protocols::loops::Loop const & new_loop,
	protocols::loops::Loops const & all_loops
){
	protocols::loops::Loops new_loops;
	for ( core::Size i=1; i<=all_loops.size(); ++i ) {
		if ( all_loops[i].start() > new_loop.start() ) {
			core::Size new_start = all_loops[i].start() + new_loop.size();
			core::Size new_stop = all_loops[i].stop() + new_loop.size();
			new_loops.push_back(new_start, new_stop);
		} else {
			new_loops.push_back(all_loops[i]);
		}
	}
	return new_loops;
}

/// @brief update the loops to reflect
///the position changes after an insertion
void
AddLoopResidues::update_anchors(
	utility::vector1<core::Size> & loop_anchors,
	protocols::loops::Loop const & new_loop,
	core::Size index_of_new_loop
){
	for ( core::Size i=index_of_new_loop+1; i<=loop_anchors.size(); ++i ) {
		loop_anchors_[i]+=new_loop.size();
	}
}

void
AddLoopResidues::dump_loops_file(
	std::string filename,
	protocols::loops::Loops loops
){
	utility::io::ozstream loops_file;
	loops_file.open(filename);
	for ( core::Size i = 1; i <= loops.size(); ++i ) {
		loops_file << "LOOP " << loops[i].start() << " " << loops[i].stop() << " 0 0.0 1\n";
	}
	loops_file.close();
}

void
AddLoopResidues::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/
){
	if ( tag->hasOption("asym_size") ) {
		asym_size_ = tag->getOption<core::Size>("asym_size");
	}


	if ( tag->hasOption("loop_sizes") ) {
		loop_sizes_.clear();
		utility::vector1<std::string> loop_sizes_strings =
			utility::string_split(tag->getOption< std::string >("loop_sizes"), ',');
		for ( core::Size i=1; i<=loop_sizes_strings.size(); ++i ) {
			loop_sizes_.push_back(utility::string2int(loop_sizes_strings[i]));
		}
	} else {
		utility_exit_with_message("You must specify desired loop sizes through the loop_sizes options");
	}

	if ( tag->hasOption("loop_anchors") ) {
		loop_anchors_.clear();

		string const loop_anchors_string = tag->getOption<string>("loop_anchors");
		utility::vector1<string> loop_anchor_strings=utility::string_split(loop_anchors_string, ',');
		for ( core::Size i=1; i<=loop_anchor_strings.size(); ++i ) {
			loop_anchors_.push_back(utility::string2int(loop_anchor_strings[i]));
		}
	} else {
		utility_exit_with_message("You must specify loop anchors");
	}
}

} //loop_creation
} //devel
