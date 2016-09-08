// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/SequenceCoupling1BDConstraint.cc
/// @brief  This is a constraint that refers to a core::sequence::SequenceProfile? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design). A note about the SequenceProfile::read_from_checkpoint function that is used to read in scores for amino acid types: the first line of the file will be ignored.
/// @author ashworth

#include <protocols/constraints_additional/SequenceCoupling1BDConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory

// option key includes

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/sequence/SequenceCoupling.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace constraints_additional {

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
static THREAD_LOCAL basic::Tracer TR( "protocols.constraints_additional.SequenceCoupling1BDConstraint" );

SequenceCoupling1BDConstraint::SequenceCoupling1BDConstraint()
:core::scoring::constraints::SequenceProfileConstraint( )
{}

SequenceCoupling1BDConstraint::SequenceCoupling1BDConstraint(
	Pose const & pose,
	core::Size numpos,
	SequenceProfileCOP profile
):core::scoring::constraints::SequenceProfileConstraint(pose, numpos,profile)
{}

SequenceCoupling1BDConstraint::SequenceCoupling1BDConstraint(
	core::Size numpos,
	SequenceProfileCOP profile
):core::scoring::constraints::SequenceProfileConstraint(numpos, profile)
{}
SequenceCoupling1BDConstraint::~SequenceCoupling1BDConstraint() = default;

ConstraintOP
SequenceCoupling1BDConstraint::clone() const {
	return ConstraintOP( new SequenceCoupling1BDConstraint( *this ) );
}

bool SequenceCoupling1BDConstraint::operator == ( core::scoring::constraints::Constraint const & other ) const
{
	return SequenceProfileConstraint::operator==( other );
}

bool SequenceCoupling1BDConstraint::same_type_as_me( core::scoring::constraints::Constraint const & other ) const
{
	return dynamic_cast< SequenceCoupling1BDConstraint const * > (&other);
}



/// @details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceCoupling1BDConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	FuncFactory const &
)
{
	Size residue_index(0);
	std::string profile_filename;

	// note: is >> "SequenceProfile" has already occured
	is >> residue_index >> profile_filename;

	TR(t_debug) << "reading: " << residue_index << " " << profile_filename << std::endl;
	if ( residue_index < 1 || residue_index > pose.size() ) {
		std::cerr << "no such residue index " << residue_index << " in pose!)" << std::endl;
		utility_exit();
	}

	seqpos(residue_index);

	// figure out sequence profile filename
	using namespace utility::file;
	// if specified, verify file exists
	if ( profile_filename != "none" ) {
		if ( ! file_exists( profile_filename ) ) {
			utility_exit_with_message( "no such file " + profile_filename );
		}
		// if filename not specified, load from commandline option -pssm only if sequence_profile_ is NULL
	} else {
		utility_exit_with_message("\"none\" is not a valid value for -pssm in this context!");
	}

	// if filename is not "none" by this point, read it even if sequence_profile_ is not currently NULL
	if ( profile_filename != "none" ) {
		SequenceCouplingOP c( new SequenceCoupling );
		c->read_from_file( FileName(profile_filename) );
		set_sequence_profile(c);
	}

	// if sequence_profile_ is still NULL by this point, it is assumed that the user intended so

} // read_def

void
SequenceCoupling1BDConstraint::show( std::ostream & os ) const {
	os << "SequenceCoupling1BD Constraint at seqpos " << seqpos() << ": ";
	if ( ! sequence_profile() ) os << "(uninitialized sequence profile)";
	os << '\n';
}

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
SequenceCoupling1BDConstraint::score(
	XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?
	runtime_assert( sequence_profile() != nullptr );

	chemical::AA aa( xyz_func.residue( seqpos()).type().aa() );
	utility::vector1< utility::vector1< Real > > const & profile( sequence_profile()->profile() );

	if ( seqpos() > profile.size() ) return; // safety/relevance check
	utility::vector1< Real > const & position_profile( profile[ seqpos() ] );
	if ( size_t(aa) > position_profile.size() ) return; // safety/relevance check
	Real const score( profile[seqpos()][aa] );
	TR(t_trace) << "seqpos " << seqpos() << " aa " << aa << " " << score << std::endl;
	emap[ this->score_type() ] += score;
}

void
SequenceCoupling1BDConstraint::fill_f1_f2(
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

} // namespace constraints_additional
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::constraints_additional::SequenceCoupling1BDConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::SequenceProfileConstraint >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::constraints_additional::SequenceCoupling1BDConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::SequenceProfileConstraint >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::constraints_additional::SequenceCoupling1BDConstraint );
CEREAL_REGISTER_TYPE( protocols::constraints_additional::SequenceCoupling1BDConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_constraints_additional_SequenceCoupling1BDConstraint )
#endif // SERIALIZATION
