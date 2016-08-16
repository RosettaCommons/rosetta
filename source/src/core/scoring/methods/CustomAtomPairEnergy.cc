// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CustomAtomPairEnergy.cc
/// @brief  Simple EnergyMethod for calculating residue-residue constraints
/// under some interaction cutoff.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/CustomAtomPairEnergy.hh>
#include <core/scoring/methods/CustomAtomPairEnergyCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <istream>

#include <core/id/AtomID.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <numeric/deriv/distance_deriv.hh>

#include <utility/vector1.hh>


// Utility headers

namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer tr( "core.scoring.methods.CustomAtomPairEnergy" );

/// @details This must return a fresh instance of the CustomAtomPairEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CustomAtomPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & opts
) const {
	return methods::EnergyMethodOP( new CustomAtomPairEnergy( opts.cst_max_seq_sep() ) );
}

ScoreTypes
CustomAtomPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( custom_atom_pair );
	return sts;
}

CustomAtomPairEnergy::CustomAtomPairEnergy( Size const max_cst_seq_sep ) :
	parent( methods::EnergyMethodCreatorOP( new CustomAtomPairEnergyCreator ) ),
	max_cst_seq_sep_( max_cst_seq_sep )
{}

/// clone
EnergyMethodOP
CustomAtomPairEnergy::clone() const
{
	return EnergyMethodOP( new CustomAtomPairEnergy(max_cst_seq_sep_) );
}


void
CustomAtomPairEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}

// make seqpos1 < seqpos2
void swap_seqpos(
	core::Size & seqpos1,
	core::Size & seqpos2
) {
	if ( seqpos1 > seqpos2 ) {
		//  Size temp = seqpos2;
		seqpos2 = seqpos1;
		seqpos1 = seqpos2;
	}
}


void
CustomAtomPairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using std::string;
	using utility::vector1;

	//string const filename( option[ score::custom_atom_pair ]() );
	string filename;
	bool success = get_comment( pose, "custom_atom_pair_fn", filename );
	if ( !success ) {
		filename = option[ score::custom_atom_pair ]();
	}

	if ( filename == fn_ ) return; // we're up to date

	tr.Debug << "reading from " << filename << std::endl;
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open input file: '" + filename + "'"
		);
	}

	Size const n_res( pose.total_residue() );
	// could save some space by only doing upper triangle (resi < resj)
	have_cst_ = vector1< vector1< bool > >( n_res, vector1< bool >( n_res, false ) );
	funcs_ = vector1< vector1< func::SOGFunc_Impl > >(
		n_res, vector1< func::SOGFunc_Impl >( n_res, func::SOGFunc_Impl() )
	);

	string line;
	while ( getline(data,line) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream line_stream(line);
		string atomi, atomj, cst_tag, func_tag;
		Size resi, resj;
		line_stream >> cst_tag >> atomi >> resi >> atomj >> resj
			>> func_tag;

		if ( cst_tag == "AtomPair" && func_tag == "SOGFUNC" ) {
			func::SOGFunc_Impl func;
			func.read_data( line_stream );
			//swap_seqpos(resi,resj);
			have_cst_[resi][resj] = true;
			have_cst_[resj][resi] = true;
			funcs_[resj][resi] = func;
			funcs_[resi][resj] = func;
			//funcs_[resj][resi].show_definition( tr.Debug );
		} else {
			tr.Error << cst_tag << "," << func_tag << std::endl;
			tr.Error << "Ignoring function on line " << line << std::endl;
		}
	}
	tr.Debug << "finished reading." << std::endl;

	tr.flush_all_channels();
}


void
CustomAtomPairEnergy::setup_for_derivatives(
	pose::Pose & pose, ScoreFunction const &
) const {
	pose.update_residue_neighbors();
	pose.update_actcoords();
}

void
CustomAtomPairEnergy::prepare_rotamers_for_packing(
	pose::Pose const & /*pose*/,
	conformation::RotamerSetBase & set
) const {
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		set.nonconst_rotamer( ii )->update_actcoord();
	}
}

void
CustomAtomPairEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const {
	pose.update_actcoord( resid );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
CustomAtomPairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size seqpos1( rsd1.seqpos() ), seqpos2( rsd2.seqpos() );
	if ( rsd1.seqpos() == rsd2.seqpos() ) return;

	// return early if we don't have a weight for this ScoreType
	if ( scorefxn.has_zero_weight(custom_atom_pair) ) return;

	// return early if we don't have constraints
	if ( !have_cst_[seqpos1][seqpos2] ) return;

	// Enforce interaction cutoff. Currently using a 10A C-alpha distance
	// cutoff. It would be nice to do this programatically from the
	// input cst file. Should be faster to use EnergyGraph code for this ...
	static std::string const atom_name( "CA" );
	if ( !( rsd1.type().has(atom_name) && rsd2.type().has(atom_name) ) ) return;
	Distance dist( rsd1.xyz(atom_name).distance( rsd2.xyz(atom_name) ) );
	if ( dist > interaction_cutoff() ) {
		emap[ custom_atom_pair ] += dist - interaction_cutoff();
	}
} // residue_pair_energy

/////////////////////////////////////////////////////////////////////////////
void
CustomAtomPairEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const & scorefxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	if ( scorefxn.has_zero_weight(custom_atom_pair) ) return;

	using core::Real;
	using core::Size;
	Size const seqpos1( atom_id.rsd() );
	if ( atom_id.atomno() != 2 ) return; // tex - fix this!

	for ( Size seqpos2 = 1; seqpos2 <= pose.total_residue(); ++seqpos2 ) {
		if ( seqpos1 == seqpos2 || ! have_cst_[seqpos1][seqpos2] ) continue;

		core::id::AtomID other_atom( atom_id.atomno(), seqpos2 );

		Real dist(0.0);
		Vector f1(0.0), f2(0.0);
		numeric::deriv::distance_f1_f2_deriv(
			pose.conformation().xyz( atom_id ), pose.conformation().xyz( other_atom ),
			dist, f1, f2 );

		if ( dist <= interaction_cutoff() ) {
			Real const wderiv(
				weights[ custom_atom_pair ] * funcs_[seqpos1][seqpos2].dfunc(dist)
			);
			F1 += f1 * wderiv;
			F2 += f2 * wderiv;
		}
	}
} // eval_atom_derivative

/// @brief CustomAtomPairEnergy distance cutoff set to be a constant
Distance
CustomAtomPairEnergy::atomic_interaction_cutoff() const
{
	return interaction_cutoff();
}

/// @details non-virtual accessor for speed; assumption: CustomAtomPairEnergy is
/// not inherited from.
Distance
CustomAtomPairEnergy::interaction_cutoff() const {
	return 10;
}

/// @brief CustomAtomPairEnergy requires no graphs (yet)
void
CustomAtomPairEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */
) const {}

/// @brief CustomAtomPairEnergy does not define intraresidue interactions
bool
CustomAtomPairEnergy::defines_intrares_energy(
	EnergyMap const & /*weights*/
) const {
	return false;
}

void
CustomAtomPairEnergy::eval_intrares_energy(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const {}
core::Size
CustomAtomPairEnergy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
