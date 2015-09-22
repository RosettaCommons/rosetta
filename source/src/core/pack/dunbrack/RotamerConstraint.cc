// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/RotamerConstraint.cc
///
/// @brief
/// @author Ian W. Davis


#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <basic/Tracer.hh>


// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/io/pdb/file_data.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

using namespace scoring;
using namespace scoring::constraints;

static THREAD_LOCAL basic::Tracer TR( "core.pack.dunbrack.RotamerConstraint" );


void load_unboundrot(pose::Pose & pose)
{
	using namespace basic::options;
	using namespace core::pose;
	if ( !option[ OptionKeys::packing::unboundrot ].active() ) return;

	static core::pose::PoseCOPs unboundrot_poses;
	if ( unboundrot_poses.empty() ) {
		for ( Size i = 1; i <= option[ OptionKeys::packing::unboundrot ]().size(); ++i ) {
			std::string filename = option[ OptionKeys::packing::unboundrot ]()[i].name();
			TR << "Adding 'unbound' rotamers from " << filename << std::endl;
			PoseOP pose( new Pose() );
			//core::import_pose::pose_from_pdb( *pose, filename );
			core::io::pdb::build_pose_from_pdb_as_is( *pose, filename );
			unboundrot_poses.push_back( pose );
		}
	}

	load_unboundrot(pose, unboundrot_poses);
}

void load_unboundrot(pose::Pose & pose, core::pose::PoseCOPs const & unboundrot_poses)
{
	if ( unboundrot_poses.empty() ) return; // guaranteed at least one pose now

	using namespace std;
	for ( Size rsd_num = 1; rsd_num <= pose.total_residue(); ++rsd_num ) {
		// Each constraint can contain only one ResidueType, so we have to sort them out here.
		// We should get any scoring overlap since ResidueTypes are mutually exclusive.
		map< string, RotamerConstraintOP > by_res_type;
		for ( Size pose_num = 1; pose_num <= unboundrot_poses.size(); ++pose_num ) {
			if ( rsd_num > unboundrot_poses[pose_num]->total_residue() ) continue;
			conformation::Residue const & rsd = unboundrot_poses[pose_num]->residue(rsd_num);
			if ( !rsd.is_protein() ) {
				// Can't determine rotamer number for anything but protein.
				TR << "Can't use " << rsd.type().name() << " " << rsd_num << " for residue constraint -- protein only." << std::endl;
				continue;
			}
			// no point in creating constrains if the residue doesn't have a rotamer library
			if ( ! core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() ) ) continue;

			if ( by_res_type.find( rsd.type().name() ) == by_res_type.end() ) { // first one, create constraint
				TR.Debug << "Creating rotamer constraint for " << rsd.type().name() << " at " << rsd_num << std::endl;
				RotamerConstraintOP constraint( new RotamerConstraint( *unboundrot_poses[pose_num], rsd_num ) );
				pose.add_constraint( constraint );
				by_res_type[ rsd.type().name() ] = constraint;
			} else { // subsequent one, just add residue
				RotamerConstraintOP constraint = by_res_type[ rsd.type().name() ];
				constraint->add_residue( unboundrot_poses[pose_num]->residue(rsd_num) );
			}
		}
	}
}

RotamerConstraint::RotamerConstraint():
	Constraint( core::scoring::fa_dun ), // most like a Dunbrack rotamer energy, and should share the same weight
	seqpos_( 0 ),
	rsd_type_name_(""),
	atom_ids_(),
	rotlib_( 0 ),
	favored_rotamers_(),
	favored_rotamer_numbers_()
{}

RotamerConstraint::RotamerConstraint(RotamerConstraint const & other):
	Constraint( core::scoring::fa_dun ), // most like a Dunbrack rotamer energy, and should share the same weight
	seqpos_( other.seqpos_ ),
	rsd_type_name_( other.rsd_type_name_ ),
	atom_ids_( other.atom_ids_ ),
	rotlib_( other.rotlib_ ),
	favored_rotamers_( other.favored_rotamers_ ),
	favored_rotamer_numbers_( other.favored_rotamer_numbers_ )
{}

RotamerConstraint::RotamerConstraint(
	pose::Pose const & pose,
	Size seqpos
):
	Constraint( core::scoring::fa_dun ), // most like a Dunbrack rotamer energy, and should share the same weight
	seqpos_( seqpos ),
	rsd_type_name_( pose.residue_type(seqpos).name() ),
	atom_ids_(),
	rotlib_( core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue_type(seqpos) ) ), // may be NULL
	favored_rotamers_(),
	favored_rotamer_numbers_()
{
	// Depends on ~ all heavy atoms (all chis + phi and psi)
	// Could cause problems if this position mutates to something with fewer atoms?
	conformation::Residue const & rsd( pose.residue(seqpos_) );
	for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
	}

	// Note: Although the Dunbrack score depends on phi & psi, it shouldn't need to depend on atoms of other residues,
	// as the phi and psi angles are internalized within a residue to speed things.

	add_residue( rsd );
}


RotamerConstraint::~RotamerConstraint() {}


void
RotamerConstraint::add_residue( conformation::Residue const & rsd )
{
	debug_assert( rsd.type().name() == rsd_type_name_ );
	favored_rotamers_.push_back( rsd.chi() );
	pack::dunbrack::RotVector rot;
	pack::dunbrack::rotamer_from_chi( rsd, rot );
	favored_rotamer_numbers_.push_back( rot );
}


Size
RotamerConstraint::natoms() const
{
	return atom_ids_.size();
}


id::AtomID const &
RotamerConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
RotamerConstraint::score(
	func::XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	rotamers::SingleResidueRotamerLibraryCOP rotlib( rotlib_ );
	if ( ! rotlib ) return;
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );

	if ( rsd.type().name() != rsd_type_name_ ) return; // residue types must match

	pack::dunbrack::RotVector rot;
	pack::dunbrack::rotamer_from_chi( rsd, rot );
	for ( Size i = 1, i_end = favored_rotamer_numbers_.size(); i <= i_end; ++i ) {
		if ( rot == favored_rotamer_numbers_[i] ) {
			pack::dunbrack::RotamerLibraryScratchSpace scratch;
			Real const best_rotE = rotlib->best_rotamer_energy(rsd, false /* => global min */, scratch);
			Real const this_rotE = rotlib->best_rotamer_energy(rsd, true /* => local min */, scratch);
			debug_assert( best_rotE <= this_rotE );
			TR << "rotamer constraint active for " << seqpos_ << " thisE = " << this_rotE << " bestE = " << best_rotE << " dE = " << ( best_rotE - this_rotE ) << std::endl;
			emap[ this->score_type() ] +=  ( best_rotE - this_rotE );
			return; // quit once we find a match
		}
	}
	// no match found, don't adjust score for this rotamer
}


void
RotamerConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero (because it's constant
	// within each rotamer well), so we just "add zero" to F1 and F2.
}


void RotamerConstraint::show( std::ostream & out ) const
{
	out << type() << ": " << seqpos_ << " " << rsd_type_name_ << " " << favored_rotamer_numbers_.size();
	for ( Size i = 1, i_end = favored_rotamer_numbers_.size(); i <= i_end; ++i ) {
		out << ",";
		for ( Size j = 1, j_end = favored_rotamer_numbers_[i].size(); j <= j_end; ++j ) {
			out << " " << favored_rotamer_numbers_[i][j];
		}
	}
	out << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core
