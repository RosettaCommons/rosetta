// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/SequenceProfileConstraint.cc
/// @brief  This is a constraint that refers to a core::sequence::SequenceProfile? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design). A note about the SequenceProfile::read_from_checkpoint function that is used to read in scores for amino acid types: the first line of the file will be ignored.
/// @author ashworth

#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/scoring/constraints/SequenceProfileConstraintCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory

#include <basic/options/option.hh>
// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

#include <boost/lexical_cast.hpp>


namespace core {
namespace scoring {
namespace constraints {

using namespace core;
using namespace chemical;
using namespace conformation;
using namespace basic::options;
using namespace scoring;
using namespace constraints;
using namespace sequence;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static thread_local basic::Tracer TR( "core.scoring.constraints.SequenceProfileConstraint" );

SequenceProfileConstraint::SequenceProfileConstraint()
	: Constraint( res_type_constraint ),
		seqpos_(0),
		sequence_profile_(/* NULL */),
		mapping_(/* NULL */),
		weight_( 1.0 )
{}

SequenceProfileConstraint::SequenceProfileConstraint(
	Pose const &,
	Size seqpos,
	SequenceProfileCOP profile /* = NULL */,
	core::id::SequenceMappingCOP mapping /* = NULL */ // current pose numbers onto profile numbers.
):
	Constraint( res_type_constraint ),
	seqpos_( seqpos ),
	sequence_profile_( profile ),
	mapping_( mapping ),
	weight_( 1.0 )
{}

SequenceProfileConstraint::SequenceProfileConstraint(
	Size seqpos,
	SequenceProfileCOP sequence_profile /* = NULL */,
	core::id::SequenceMappingCOP mapping /* = NULL */ // current pose numbers onto profile numbers.
):
	Constraint( res_type_constraint ),
	seqpos_( seqpos ),
	sequence_profile_( sequence_profile ),
	mapping_( mapping ),
	weight_( 1.0 )
{}

SequenceProfileConstraint::~SequenceProfileConstraint() {}

ConstraintOP
SequenceProfileConstraint::clone() const {
	return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( *this ) );
}

///@details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceProfileConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	func::FuncFactory const &
)
{

	int version = 0;
//	note: is >> "SequenceProfile" has already occured
	is >> version;
	if( version == -1 ) {
		std::string tmp;
		std::getline( is, tmp );
		std::stringstream tmpss(tmp);
		Size residue_index, aa_count;
		tmpss >> residue_index >> aa_count;
		utility::vector1<Real> aa_scores;
		TR(t_debug) << "Loading seqprof constraint for resi " << residue_index << " with aa_count " << aa_count << std::endl;
		for( Size i = 1; i <= aa_count; i++ ) {
			Real tmp;
			tmpss >> tmp;
			aa_scores.push_back( tmp );
		}
		if( ! sequence_profile_ )  {
			SequenceProfileOP newseqprof( new SequenceProfile );
			sequence_profile_ = newseqprof;
			mapping_ = NULL;
		}
		seqpos_ = residue_index;
		const_cast<SequenceProfile * >(sequence_profile_.get())->prof_row( aa_scores, residue_index );
	} else {
		Size residue_index(boost::lexical_cast<Size>(version));
		std::string profile_filename;

		is >> profile_filename;

		TR(t_debug) << "reading: " << residue_index << " " << profile_filename << std::endl;
		if ( residue_index < 1 || residue_index > pose.total_residue() ) {
			std::cerr << "no such residue index " << residue_index << " in pose!)" << std::endl;
			utility_exit();
		}

		seqpos_ = residue_index;

		// figure out sequence profile filename
		using namespace utility::file;
		// if specified, verify file exists
		if ( profile_filename != "none" ) {
			if ( ! file_exists( profile_filename ) ) {
				utility_exit_with_message( "no such file " + profile_filename );
			}
			// if filename not specified, load from commandline option -pssm only if sequence_profile_ is NULL
		} else if ( !sequence_profile_ && option[ OptionKeys::in::file::pssm ].user() ) {
			profile_filename = option[ OptionKeys::in::file::pssm ]().front();
			// ps. don't name your file "none"
			if ( profile_filename == "none" ) {
				utility_exit_with_message("\"none\" is not a valid value for -pssm in this context!");
			}
		}

		// if filename is not "none" by this point, read it even if sequence_profile_ is not currently NULL
		if ( profile_filename != "none" ) {
			SequenceProfileOP newseqprof( new SequenceProfile );
			newseqprof->read_from_checkpoint( FileName(profile_filename) );
			sequence_profile_ = newseqprof;
			mapping_ = NULL; // Reset mapping
		}

		// if sequence_profile_ is still NULL by this point, it is assumed that the user intended so
	}

} // read_def

void
SequenceProfileConstraint::show_def( std::ostream & os, Pose const & ) const
{
	show( os );
}

void
SequenceProfileConstraint::show( std::ostream & os ) const {
	if ( sequence_profile_ )  {
		core::Size profile_pos( seqpos_ );
		if( mapping_ ) {
			profile_pos = (*mapping_)[seqpos_];
			if( profile_pos == 0 ) return; // safety/relevance check
		}
		typedef utility::vector1<Real> RealVec;
		//if( profile_pos == 1 )
//			os << "SequenceProfile start" << std::endl;

		if( profile_pos <= sequence_profile_->size() ) {
			RealVec const & aa_scores( sequence_profile_->prof_row( profile_pos ) );
			os << "SequenceProfile -1 " << profile_pos << " " << aa_scores.size() << " ";
			for ( Size aa(1); aa <= aa_scores.size(); ++aa ) {
				os << aa_scores[aa] << " ";
			}
			if( profile_pos != sequence_profile_->size() )
				os << '\n';
		}
	}
}

void
SequenceProfileConstraint::set_sequence_profile( SequenceProfileCOP profile, core::id::SequenceMappingCOP mapping /* = NULL */ )
{
	sequence_profile_ = profile;
	mapping_ = mapping;
}

SequenceProfileCOP
SequenceProfileConstraint::sequence_profile() const { return sequence_profile_; }

core::id::SequenceMappingCOP
SequenceProfileConstraint::profile_mapping() const { return mapping_; }

utility::vector1< core::Size >
SequenceProfileConstraint::residues() const {
	utility::vector1< core::Size > pos_list(1, seqpos_); // length 1 containing "all" seqpos_ values
	return pos_list;
}

ConstraintOP
SequenceProfileConstraint::remap_resid(
	SequenceMapping const & seqmap
) const {
	Size newseqpos( seqmap[ seqpos_ ] );
	if ( newseqpos == 0 ) { return NULL; }

	TR(t_debug) << "Remapping resid " << seqpos_ << " to " << newseqpos << std::endl;

	core::id::SequenceMappingOP new_map( new core::id::SequenceMapping( seqmap ) );
	new_map->reverse(); // Go from old->new to new->old
	if ( mapping_ ) {
		new_map->downstream_combine( *mapping_ ); // Combine new->old and old->profile into new->profile
	}

	return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint(	newseqpos, sequence_profile_, new_map ) );
}

ConstraintOP
SequenceProfileConstraint::remapped_clone(
	pose::Pose const&,
	pose::Pose const&,
	id::SequenceMappingCOP map /*=NULL*/ ) const
{
	if ( ! map ) {
		return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint(	seqpos_, sequence_profile_, mapping_ ) );
	}

	// Hereafter map is valid
	Size newseqpos( (*map)[ seqpos_ ] );

	if ( newseqpos == 0 ) { return NULL; }

	TR(t_debug) << "Remapping resid " << seqpos_ << " to " << newseqpos << std::endl;

	core::id::SequenceMappingOP new_map( new core::id::SequenceMapping( *map ) );
	new_map->reverse(); // Go from old->new to new->old
	if ( mapping_ ) {
		new_map->downstream_combine( *mapping_ ); // Combine new->old and old->profile into new->profile
	}

	return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint(	newseqpos, sequence_profile_, new_map ) );
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
SequenceProfileConstraint::score(
	func::XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?
	runtime_assert( sequence_profile_ != 0 );

	core::Size profile_pos( seqpos_ );
	if( mapping_ ) {
		profile_pos = (*mapping_)[seqpos_];
		if( profile_pos == 0 ) return; // safety/relevance check
	}
	chemical::AA aa( xyz_func.residue( seqpos_ ).type().aa() );
	utility::vector1< utility::vector1< Real > > const & profile( sequence_profile_->profile() );
	if ( profile_pos > profile.size() ) return; // safety/relevance check
	utility::vector1< Real > const & position_profile( profile[ profile_pos ] );
	if ( size_t(aa) > position_profile.size() ) return; // safety/relevance check
	Real const score( position_profile[aa] );
	TR(t_trace) << "seqpos " << seqpos_ << " aa " << aa << " " << weight() * score << std::endl;

	if( sequence_profile_->negative_better() ) {
		emap[ this->score_type() ] += weight() * score;
	} else {
		emap[ this->score_type() ] -= weight() * score;
	}
}

void
SequenceProfileConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const & ,//conformation,
	Vector & ,//f1,
	Vector & ,//f2,
	EnergyMap const & //weights
) const
{
	// Do nothing, as the value of this function doesn't change with respect to
	// the torsions.
}

SequenceProfileConstraintCreator::SequenceProfileConstraintCreator() {}
SequenceProfileConstraintCreator::~SequenceProfileConstraintCreator() {}

core::scoring::constraints::ConstraintOP
SequenceProfileConstraintCreator::create_constraint() const
{
        return core::scoring::constraints::ConstraintOP( new SequenceProfileConstraint );
}

std::string
SequenceProfileConstraintCreator::keyname() const
{
        return "SequenceProfile";
}

core::Real
SequenceProfileConstraint::weight() const{
	return weight_;
}

void
SequenceProfileConstraint::weight( core::Real const w ){
	weight_ = w;
}

} // namespace constraints
} // namespace scoring
} // namespace core
