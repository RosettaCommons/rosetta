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
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/vector1.hh>


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
static basic::Tracer TR("core.scoring.constraints.SequenceProfileConstraint");

SequenceProfileConstraint::SequenceProfileConstraint()
	: Constraint( res_type_constraint ),
		seqpos_(0),
		sequence_profile_(NULL)
{}

SequenceProfileConstraint::SequenceProfileConstraint(
	Pose const & pose,
	Size seqpos,
	SequenceProfileOP profile /* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos_( seqpos ),
	sequence_profile_( profile )
{
	Residue const & rsd( pose.residue(seqpos_) );
	for( Size i(1), i_end( rsd.nheavyatoms() ); i <= i_end; ++i ) {
		atom_ids_.push_back( AtomID( i, seqpos_ ) );
	}
}

SequenceProfileConstraint::SequenceProfileConstraint(
	Size seqpos,
	utility::vector1< AtomID > const & atoms_in,
	SequenceProfileOP sequence_profile /* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos_( seqpos ),
	sequence_profile_( sequence_profile )
{
	for( utility::vector1< AtomID >::const_iterator at_it( atoms_in.begin() ), end( atoms_in.end() ); at_it != end; ++at_it ) {
		atom_ids_.push_back( *at_it );
	}
}

SequenceProfileConstraint::~SequenceProfileConstraint() {}

ConstraintOP
SequenceProfileConstraint::clone() const {
	return new core::scoring::constraints::SequenceProfileConstraint( *this );
}

///@details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceProfileConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	FuncFactory const &
)
{
	Size residue_index(0);
	std::string profile_filename;

//	note: is >> "SequenceProfile" has already occured
	is >> residue_index >> profile_filename;

	TR(t_debug) << "reading: " << residue_index << " " << profile_filename << std::endl;
	if ( residue_index < 1 || residue_index > pose.total_residue() ) {
		std::cerr << "no such residue index " << residue_index << " in pose!)" << std::endl;
		utility_exit();
	}

	seqpos_ = residue_index;

	Residue const & rsd( pose.residue(seqpos_) );
	for( Size i(1), i_end( rsd.nheavyatoms() ); i <= i_end; ++i ) {
		atom_ids_.push_back( AtomID( i, seqpos_ ) );
	}

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
		sequence_profile_ = new SequenceProfile;
		sequence_profile_->read_from_checkpoint( FileName(profile_filename) );
	}

	// if sequence_profile_ is still NULL by this point, it is assumed that the user intended so

} // read_def

void
SequenceProfileConstraint::show_def( std::ostream & os, Pose const & ) const
{
	show( os );
}

void
SequenceProfileConstraint::show( std::ostream & os ) const {
	os << "SequenceProfileConstraint at seqpos " << seqpos_ << ": ";
	if ( ! sequence_profile_ ) os << "(uninitialized sequence profile)";
//	else {
//		typedef utility::vector1<Real> RealVec;
//		RealVec const & aa_scores( sequence_profile_->prof_row( seqpos_ ) );
//		runtime_assert( aa_scores.size() >= num_canonical_aas );
//		for ( Size aa(1); aa <= num_canonical_aas; ++aa ) {
//			os << aa_scores[aa] << " ";
//		}
//	}
	os << '\n';
}

void
SequenceProfileConstraint::set_sequence_profile( SequenceProfileOP profile )
{
	sequence_profile_ = profile;
}

SequenceProfileOP
SequenceProfileConstraint::sequence_profile() { return sequence_profile_; }

SequenceProfileCOP
SequenceProfileConstraint::sequence_profile() const { return sequence_profile_; }

Size
SequenceProfileConstraint::natoms() const
{
	return atom_ids_.size();
}

id::AtomID const &
SequenceProfileConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}

utility::vector1< id::AtomID > const &
SequenceProfileConstraint::atom_ids() const {
	return atom_ids_;
}

void SequenceProfileConstraint::atom_ids(utility::vector1< id::AtomID > & atomIds){
	atom_ids_ = atomIds;
}

ConstraintOP
SequenceProfileConstraint::remap_resid(
	SequenceMapping const & seqmap
) const {
	Size newseqpos( seqmap[ seqpos_ ] );
	if ( newseqpos != 0 ) {
		TR(t_debug) << "Remapping resid " << seqpos_ << " to " << newseqpos << std::endl;

		utility::vector1< AtomID > new_atomids;
		for ( utility::vector1< AtomID >::const_iterator at_it( atom_ids_.begin() ), end( atom_ids_.end() ); at_it != end; ++at_it ) {
			if ( seqmap[ at_it->rsd() ] != 0 ) {
				new_atomids.push_back( AtomID( at_it->atomno(), seqmap[ at_it->rsd() ] ) );
			}
		}
		return new core::scoring::constraints::SequenceProfileConstraint(	newseqpos, new_atomids, sequence_profile_ );
	}
	else return NULL;
}

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
SequenceProfileConstraint::score(
	XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?
	runtime_assert( sequence_profile_ );

	chemical::AA aa( xyz_func.residue( seqpos_ ).type().aa() );
	utility::vector1< utility::vector1< Real > > const & profile( sequence_profile_->profile() );
	if ( seqpos_ > profile.size() ) return; // safety/relevance check
	utility::vector1< Real > const & position_profile( profile[ seqpos_ ] );
	if ( size_t(aa) > position_profile.size() ) return; // safety/relevance check
	// this assumes that the profile contains energies and not probabilities
	Real const score( sequence_profile_->profile()[seqpos_][aa] );
	TR(t_trace) << "seqpos " << seqpos_ << " aa " << aa << " " << score << std::endl;

	if( sequence_profile_->negative_better() ) {
		emap[ this->score_type() ] += score;
	} else {
		emap[ this->score_type() ] -= score;
	}
}

void
SequenceProfileConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	XYZ_Func const & ,//conformation,
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
        return new SequenceProfileConstraint;
}

std::string
SequenceProfileConstraintCreator::keyname() const
{
        return "SequenceProfile";
}

} // namespace constraints
} // namespace scoring
} // namespace core
