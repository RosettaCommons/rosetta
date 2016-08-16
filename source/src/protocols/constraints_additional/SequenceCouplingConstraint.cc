// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/SequnceCouplingConstraint.cc
/// @brief  This is a constraint that refers to a core::sequence::SequnceCopuling? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design). A note about the SequenceProfile::read_from_checkpoint function that is used to read in scores for amino acid types: the first line of the file will be ignored.
/// @author HetuKamisetty

#include <protocols/constraints_additional/SequenceCouplingConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/sequence/SequenceCoupling.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory

// option key includes

//Auto Headers
#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>

//#include <core/id/SequenceMapping.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
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
static THREAD_LOCAL basic::Tracer TR( "protocols.constraints_additional.SequenceCouplingConstraint" );

SequenceCouplingConstraint::SequenceCouplingConstraint()
: Constraint( res_type_constraint ),
	seqpos1_(0),
	seqpos2_(0),
	sequence_coupling_(/* NULL */)
{}

SequenceCouplingConstraint::SequenceCouplingConstraint(
	Pose const &,
	Size seqpos1,
	Size seqpos2,
	SequenceCouplingOP coupling/* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos1_(seqpos1),
	seqpos2_(seqpos2),
	sequence_coupling_( coupling)
{}

SequenceCouplingConstraint::SequenceCouplingConstraint(
	Size seqpos1,
	Size seqpos2,
	SequenceCouplingOP sequence_coupling /* = NULL */
):
	Constraint( res_type_constraint ),
	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 ),
	sequence_coupling_( sequence_coupling )
{}

SequenceCouplingConstraint::~SequenceCouplingConstraint() {}

ConstraintOP
SequenceCouplingConstraint::clone() const {
	return ConstraintOP( new SequenceCouplingConstraint( *this ) );
}

bool SequenceCouplingConstraint::operator == ( core::scoring::constraints::Constraint const & other ) const
{
	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( other ) ) return false;

	SequenceCouplingConstraint const & other_downcast( static_cast< SequenceCouplingConstraint const & > ( other ) );
	if ( seqpos1_ != other_downcast.seqpos1_ ) return false;
	if ( seqpos2_ != other_downcast.seqpos2_ ) return false;
	return sequence_coupling_ == other_downcast.sequence_coupling_ ||
		( sequence_coupling_ && other_downcast.sequence_coupling_ && *sequence_coupling_ == *other_downcast.sequence_coupling_ );
}

bool SequenceCouplingConstraint::same_type_as_me( core::scoring::constraints::Constraint const & other ) const
{
	return dynamic_cast< SequenceCouplingConstraint const * > (&other);
}


/// @details one line definition "SequenceProfile resindex profilefilename" (profilefilename can also be set to "none" in the constraints file, and specified by -in::file::pssm)
void
SequenceCouplingConstraint::read_def(
	std::istream & is,
	Pose const & pose,
	core::scoring::func::FuncFactory const &
)
{
	Size residue_index1(0);
	Size residue_index2(0);
	std::string coupling_filename;

	// note: is >> "SequenceProfile" has already occured
	is >> residue_index1 >> residue_index2 >> coupling_filename;

	TR(t_debug) << "reading: " << residue_index1 << " " << residue_index2 << " " << coupling_filename << std::endl;
	if ( residue_index1 < 1 || residue_index1 > pose.total_residue() ) {
		std::cerr << "no such residue index " << residue_index1 << " in pose!)" << std::endl;
		utility_exit();
	}
	if ( residue_index2 < 1 || residue_index2 > pose.total_residue() ) {
		std::cerr << "no such residue index " << residue_index2 << " in pose!)" << std::endl;
		utility_exit();
	}

	seqpos1_ = residue_index1;
	seqpos2_ = residue_index2;

	// figure out sequence profile filename
	using namespace utility::file;
	// if specified, verify file exists
	if ( coupling_filename!= "none" ) {
		if ( ! file_exists( coupling_filename ) ) {
			utility_exit_with_message( "no such file " + coupling_filename );
		}
		// if filename not specified, load from commandline option -pssm only if sequence_coupling_ is NULL
	} else {
		utility_exit_with_message("\"none\" is not a valid value for -pssm in this context!");
	}

	// if filename is not "none" by this point, read it even if sequence_coupling_ is not currently NULL
	if ( coupling_filename != "none" ) {
		sequence_coupling_ = SequenceCouplingOP( new SequenceCoupling );
		sequence_coupling_->read_from_file( FileName(coupling_filename) );
	}

	// if sequence_coupling_ is still NULL by this point, it is assumed that the user intended so

} // read_def

void
SequenceCouplingConstraint::show_def( std::ostream & os, Pose const & ) const
{
	show( os );
}

void
SequenceCouplingConstraint::show( std::ostream & os ) const {
	os << "SequenceCouplingConstraint between seqpos " << seqpos1_ << " " << seqpos2_ << ": ";
	if ( ! sequence_coupling_ ) os << "(uninitialized sequence profile)";
	// else {
	//  typedef utility::vector1<Real> RealVec;
	//  RealVec const & aa_scores( sequence_profile_->prof_row( seqpos_ ) );
	//  runtime_assert( aa_scores.size() >= num_canonical_aas );
	//  for ( Size aa(1); aa <= num_canonical_aas; ++aa ) {
	//   os << aa_scores[aa] << " ";
	//  }
	// }
	os << '\n';
}

void
SequenceCouplingConstraint::set_sequence_coupling( SequenceCouplingOP profile )
{
	sequence_coupling_ = profile;
}

SequenceCouplingOP
SequenceCouplingConstraint::sequence_coupling() { return sequence_coupling_; }

SequenceCouplingCOP
SequenceCouplingConstraint::sequence_coupling() const { return sequence_coupling_; }

utility::vector1< core::Size >
SequenceCouplingConstraint::residues() const {
	utility::vector1< core::Size > pos_list;
	pos_list.push_back(seqpos1_);
	pos_list.push_back(seqpos2_);
	return pos_list;
}

/*
* hk: how does one fix sequencecoupling on a remap_resid call?
*
ConstraintOP
SequenceCouplingConstraint::remap_resid(
SequenceMapping const & seqmap
) const {
Size newseqpos1( seqmap[ seqpos1_ ] );
if ( newseqpos != 0 ) {
TR(t_debug) << "Remapping resid " << seqpos1_ << " to " << newseqpos1 << std::endl;

return new SequenceCouplingConstraint( newseqpos1, newseqpos2, sequence_coupling_ );
}
else return NULL;
}
*/

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
SequenceCouplingConstraint::score(
	core::scoring::func::XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?
	runtime_assert( sequence_coupling_ != 0 );

	chemical::AA aa1( xyz_func.residue( seqpos1_ ).type().aa() );
	chemical::AA aa2( xyz_func.residue( seqpos2_ ).type().aa() );
	if ( seqpos1_ > sequence_coupling_->npos() || seqpos2_ > sequence_coupling_->npos() ) return; // safety/relevance check

	Size edgeId = sequence_coupling_->findEdgeId(seqpos1_, seqpos2_);
	Real score( 0);
	if ( edgeId<=0 ) { //direction important. if (i,j) edge exists, will not return if (j,i) asked
		edgeId = sequence_coupling_->findEdgeId(seqpos2_, seqpos1_);

		//epot = sequence_coupling_->edgePotBetween(edgeId);
		utility::vector1< utility::vector1 < Real > > const &  epot(sequence_coupling_->edgePotBetween(edgeId));
		//runtime_assert(epot);
		/*
		if(epot==NULL){
		std::cerr << "no such edge " << seqpos1_ << " " << seqpos2_ << " in sequence coupling !)" << std::endl;
		utility_exit();
		}
		*/
		score = epot[aa2][aa1];
	} else {
		utility::vector1< utility::vector1 < Real > > const &  epot(sequence_coupling_->edgePotBetween(edgeId));
		//runtime_assert(epot);
		score = epot[aa1][aa2];
	}

	TR(t_trace) << "seqpos1 " << seqpos1_ << " aa1 " << aa1 << " " << "seqpos2 " << seqpos2_<< " aa2 "<< aa2 << " " << score << std::endl;

	emap[ this->score_type() ] += score;//pot scores are like energies; lower better.
}

void
SequenceCouplingConstraint::fill_f1_f2(
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
protocols::constraints_additional::SequenceCouplingConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( CEREAL_NVP( seqpos1_ ) ); // core::Size
	arc( CEREAL_NVP( seqpos2_ ) ); // core::Size
	arc( CEREAL_NVP( sequence_coupling_ ) ); // SequenceCouplingOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::constraints_additional::SequenceCouplingConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( seqpos1_ ); // core::Size
	arc( seqpos2_ ); // core::Size
	arc( sequence_coupling_ ); // SequenceCouplingOP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::constraints_additional::SequenceCouplingConstraint );
CEREAL_REGISTER_TYPE( protocols::constraints_additional::SequenceCouplingConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_constraints_additional_SequenceCouplingConstraint )
#endif // SERIALIZATION
