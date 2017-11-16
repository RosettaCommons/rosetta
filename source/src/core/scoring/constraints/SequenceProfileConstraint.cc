// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


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
static basic::Tracer TR( "core.scoring.constraints.SequenceProfileConstraint" );

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
	return SequenceProfileConstraintOP( new SequenceProfileConstraint(*this) );
}

bool SequenceProfileConstraint::operator == ( Constraint const & rhs ) const {
	if ( ! same_type_as_me( rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	SequenceProfileConstraint const & rhs_spc( static_cast< SequenceProfileConstraint const & > (rhs) );
	if ( seqpos_ != rhs_spc.seqpos_ ) return false;
	if ( sequence_profile_ != rhs_spc.sequence_profile_ && !( *sequence_profile_ == *rhs_spc.sequence_profile_ ) ) return false;
	if ( mapping_ != rhs_spc.mapping_ && ! ( *mapping_ == *rhs_spc.mapping_ ) ) return false;
	if ( weight_ != rhs_spc.weight_ ) return false;

	return true;
}

bool SequenceProfileConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< SequenceProfileConstraint const * > (&other);
}


/// @details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceProfileConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	func::FuncFactory const &
)
{

	int version = 0;
	// note: is >> "SequenceProfile" has already occured
	is >> version;
	if ( version == -1 ) {
		std::string tmp;
		std::getline( is, tmp );
		std::stringstream tmpss(tmp);
		Size residue_index, aa_count;
		tmpss >> residue_index >> aa_count;
		utility::vector1<Real> aa_scores;
		TR(t_debug) << "Loading seqprof constraint for resi " << residue_index << " with aa_count " << aa_count << std::endl;
		for ( Size i = 1; i <= aa_count; i++ ) {
			Real tmp;
			tmpss >> tmp;
			aa_scores.push_back( tmp );
		}
		if ( ! sequence_profile_ )  {
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
		if ( residue_index < 1 || residue_index > pose.size() ) {
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
		if ( mapping_ ) {
			profile_pos = (*mapping_)[seqpos_];
			if ( profile_pos == 0 ) return; // safety/relevance check
		}
		typedef utility::vector1<Real> RealVec;
		//if( profile_pos == 1 )
		//   os << "SequenceProfile start" << std::endl;

		if ( profile_pos <= sequence_profile_->size() ) {
			RealVec const & aa_scores( sequence_profile_->prof_row( profile_pos ) );
			os << "SequenceProfile -1 " << seqpos_ << " "<< profile_pos << " " << aa_scores.size() << " ";
			for ( Size aa(1); aa <= aa_scores.size(); ++aa ) {
				os << aa_scores[aa] << " ";
			}
			if ( profile_pos != sequence_profile_->size() ) {
				os << '\n';
			}
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

	return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( newseqpos, sequence_profile_, new_map ) );
}

ConstraintOP
SequenceProfileConstraint::remapped_clone(
	pose::Pose const&,
	pose::Pose const&,
	id::SequenceMappingCOP map /*=NULL*/ ) const
{
	if ( ! map ) {
		return clone();
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

	return ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint(
		newseqpos,
		sequence_profile_ ?
		utility::pointer::dynamic_pointer_cast< sequence::SequenceProfile > ( sequence_profile_->clone() ) :
		sequence_profile_,
		new_map ) );
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

	chemical::AA aa( xyz_func.residue( seqpos_ ).type().aa() );
	emap[ this->score_type() ] += weight() * dist( aa );

}

core::Real
SequenceProfileConstraint::dist( chemical::AA aa ) const {
	runtime_assert( sequence_profile_ != 0 );

	core::Size profile_pos( seqpos_ );
	if ( mapping_ ) {
		profile_pos = (*mapping_)[seqpos_];
		if ( profile_pos == 0 ) return 0.0; // safety/relevance check
	}
	utility::vector1< utility::vector1< Real > > const & profile( sequence_profile_->profile() );
	if ( profile_pos > profile.size() ) return 0.0; // safety/relevance check
	utility::vector1< Real > const & position_profile( profile[ profile_pos ] );
	if ( size_t(aa) > position_profile.size() ) return 0.0; // safety/relevance check
	Real const score( position_profile[aa] );
	TR(t_trace) << "seqpos " << seqpos_ << " aa " << aa << " " << weight() * score << std::endl;
	if ( sequence_profile_->negative_better() ) {
		return score;
	} else {
		return -score;
	}
}


core::Real
SequenceProfileConstraint::dist( core::scoring::func::XYZ_Func const & xyz_func ) const {
	chemical::AA aa( xyz_func.residue( seqpos_ ).type().aa() );
	return dist( aa );
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

SequenceProfileConstraint::SequenceProfileConstraint( SequenceProfileConstraint const & src ) :
	Constraint( src ),
	seqpos_( src.seqpos_ ),
	sequence_profile_( src.sequence_profile_ ?
	utility::pointer::dynamic_pointer_cast< SequenceProfile > ( src.sequence_profile_->clone() ) : src.sequence_profile_ ),
	mapping_( src.mapping_ ? core::id::SequenceMappingCOP( new core::id::SequenceMapping( *src.mapping_ )) : src.mapping_ ),
	weight_(src.weight_ )
{}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::SequenceProfileConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( CEREAL_NVP( seqpos_ ) ); // core::Size
	arc( CEREAL_NVP( sequence_profile_ ) ); // SequenceProfileCOP
	arc( CEREAL_NVP( mapping_ ) ); // core::id::SequenceMappingCOP
	arc( CEREAL_NVP( weight_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::SequenceProfileConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( seqpos_ ); // core::Size
	std::shared_ptr< core::sequence::SequenceProfile > local_sequence_profile;
	arc( local_sequence_profile ); // SequenceProfileCOP
	sequence_profile_ = local_sequence_profile; // copy the non-const pointer(s) into the const pointer(s)
	std::shared_ptr< core::id::SequenceMapping > local_mapping;
	arc( local_mapping ); // core::id::SequenceMappingCOP
	mapping_ = local_mapping; // copy the non-const pointer(s) into the const pointer(s)
	arc( weight_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::SequenceProfileConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::SequenceProfileConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_SequenceProfileConstraint )
#endif // SERIALIZATION
