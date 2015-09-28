// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/HbondsToResidueFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
/// @author Refactored considerably by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilterCreator.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/MetricValue.hh>
#include <numeric/random/random.hh>
#include <core/chemical/AtomType.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//Objectxxxx header
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility Headers

// Unit Headers
#include <protocols/simple_moves/ddG.hh>
#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <map>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <protocols/simple_filters/DdgFilter.hh>


using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;


namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.HbondsToResidueFilter" );
using core::pose::Pose;

protocols::filters::FilterOP
HbondsToResidueFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HbondsToResidueFilter ); }

std::string
HbondsToResidueFilterCreator::keyname() const { return "HbondsToResidue"; }

/// @brief Default constructor.
///
HbondsToResidueFilter::HbondsToResidueFilter() :
	Filter( "HbondsToResidue" ),
	resnum_(""),
	partners_(0),
	energy_cutoff_(0.0),
	backbone_(true),
	sidechain_(true),
	bb_bb_(true),
	from_other_chains_(true),
	from_same_chain_(true),
	sfxn_()
{}

/// @brief Constructor
///
HbondsToResidueFilter::HbondsToResidueFilter(
	Size const resnum,
	Size const partners,
	Real const &energy_cutoff,
	bool const backbone,
	bool const sidechain,
	bool const bb_bb,
	bool const from_other_chains,
	bool const from_same_chain
) :
	Filter( "HbondsToResidue" ),
	resnum_(""),
	partners_(partners),
	energy_cutoff_(energy_cutoff),
	backbone_(backbone),
	sidechain_(sidechain),
	bb_bb_(bb_bb),
	from_other_chains_(from_other_chains),
	from_same_chain_(from_same_chain),
	sfxn_()
{
	std::stringstream resnumstr;
	resnumstr << resnum;
	resnum_ = resnumstr.str();

	runtime_assert( energy_cutoff_ <= 0 );
}

/// @brief Copy constructor
///
HbondsToResidueFilter::HbondsToResidueFilter( HbondsToResidueFilter const &src ) :
	Filter( "HbondsToResidue" ),
	resnum_(src.resnum_),
	partners_(src.partners_),
	energy_cutoff_(src.energy_cutoff_),
	backbone_(src.backbone_),
	sidechain_(src.sidechain_),
	bb_bb_(src.bb_bb_),
	from_other_chains_(src.from_other_chains_),
	from_same_chain_(src.from_same_chain_),
	sfxn_( )
{
	if( src.sfxn_ ) {
		sfxn_ = src.sfxn_->clone();
	}
}

bool
HbondsToResidueFilter::apply( Pose const & pose ) const {
	core::Size const resnum_rosetta( core::pose::parse_resnum( resnum_, pose, true ) );
	core::Size hbonded_res( compute( pose, resnum_rosetta ) ); //Number of hbonded residues
	if ( TR.visible() ) TR<<"found "<<hbonded_res<< " hbond to target residue " << resnum_rosetta;
	if ( hbonded_res >= partners_ ) {
		if ( TR.visible() ) TR << ". passing." << std::endl;
		return( true );
	} else {
		if ( TR.visible() ) TR << ". failing." << std::endl;
		return( false );
	}
}

void
HbondsToResidueFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
) {
	runtime_assert_string_msg( !tag->hasOption("res_num") && !tag->hasOption("pdb_num"), "Error in HbondsToResidueFilter::parse_my_tag():  The \"res_num\" and \"pdb_num\" options have been deprecated.  Use \"residue\" instead, and provide a Rosetta number (e.g. \"32\"), a PDB number (e.g. \"12B\"), or a reference pose number (e.g. \"refpose(snapshot1,17)+3\").");
	runtime_assert_string_msg( tag->hasOption("residue"), "Error in HbondsToResidueFilter::parse_my_tag():  No \"residue\" option found in the options for the HbondsToResidue filter." );

	set_partners( tag->getOption<core::Size>( "partners" ) );
	set_energy_cutoff( tag->getOption<core::Real>( "energy_cutoff", -0.5 ) );
	set_bb_bb( tag->getOption<bool>( "bb_bb", 0 ) );
	set_backbone( tag->getOption<bool>( "backbone", 0 ) );
	set_sidechain( tag->getOption<bool>( "sidechain", 1 ) );
	set_resnum( tag->getOption<std::string>( "residue" ) );

	if ( tag->hasOption("scorefxn") ) {
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	}

	set_from_same_chain( tag->getOption<bool>( "from_same_chain", true ) );
	set_from_other_chains( tag->getOption<bool>( "from_other_chains", true ) );
	runtime_assert_string_msg( from_same_chain() || from_other_chains(), "Error in HbondsToResidueFilter::parse_my_tag():  The user has specified that neither hydrogen bonds to other chains, nor hydrogen bonds within the same chain should be considered.  This means that there is nothing for this filter to count.  Please turn on either \"from_same_chain\" or \"from_other_chains\" (or both)." );

	if ( TR.visible() ) {
		TR << "Hbonds to residue filter for resnum " << resnum() << " with " << partners() << " hbonding partners as the cutoff threshold." << std::endl;
		if ( from_other_chains() && from_same_chain() ) {
			TR << "Hbonds to other chains and to the same chain will be considered." << std::endl;
		} else if ( from_other_chains() && !from_same_chain() ) {
			TR << "Hbonds to other chains, but not to the same chain, will be considered." << std::endl;
		} else if ( !from_other_chains() && from_same_chain() ) {
			TR << "Hbonds to the same chain, but not to other chains, will be considered." << std::endl;
		}
		TR.flush();
	}
}

void
HbondsToResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const resnum_rosetta( core::pose::parse_resnum(resnum(), pose, true) );
	core::Size hbonded_res( compute( pose, resnum_rosetta ) );

	out<<"Number of residues hbonded to "<<resnum_rosetta<< " is " << hbonded_res <<'\n';
}

core::Real
HbondsToResidueFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size hbonded_res( compute( pose, core::pose::parse_resnum(resnum(), pose, true) ) );
	return( hbonded_res );
}

core::Size
HbondsToResidueFilter::compute( Pose const & pose, core::Size const resnum_rosetta ) const {
	using core::Size;

	core::pose::Pose temp_pose( pose );
	core::scoring::ScoreFunctionOP scorefxn( sfxn_ );
	if ( !scorefxn ) {
		TR << "No scorefunction loaded.  Getting global default scorefunction." << std::endl; //DELETE ME.
		scorefxn=get_score_function();
	}
	(*scorefxn)(temp_pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( temp_pose );

	std::set<Size> binders;
	for ( Size i=1, imax=pose.n_residue(); i<=imax; ++i ) {
		if ( i == resnum_rosetta ) continue; //Don't consider hbonds of this residue to itself.
		if ( pose.chain(i) == pose.chain(resnum_rosetta) && !from_same_chain() ) continue; //Skip hbonds from same chain if the from_same_chain option is not set.
		if ( pose.chain(i) != pose.chain(resnum_rosetta) && !from_other_chains() ) continue; //Skip hbonds from different chains if the from_other_chain option is not set.
		binders.insert( i );
	}
	std::list< Size> hbonded_res( hbonded( temp_pose, resnum_rosetta, binders, backbone_, sidechain_, energy_cutoff_, bb_bb_, scorefxn) );

	return( hbonded_res.size() );
}

HbondsToResidueFilter::~HbondsToResidueFilter() {}

} // filters
} // protein_interface_design
} // protocols
