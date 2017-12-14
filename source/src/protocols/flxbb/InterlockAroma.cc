// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/flxbb/InterlockAroma.cc
/// @brief place interlocking aromatic residues
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#include <protocols/flxbb/InterlockAroma.hh>
#include <protocols/flxbb/InterlockAromaCreator.hh>

#include <protocols/fldsgn/filters/InterlockingAromaFilter.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Parser headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>


#include <utility/vector0.hh>

//Auto Headers
#include <utility/graph/Graph.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
static basic::Tracer TR( "protocols.flxbb.InterlockAroma" );

namespace protocols {
namespace flxbb {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP std::string
// XRW TEMP InterlockAromaCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return InterlockAroma::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP InterlockAromaCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new InterlockAroma );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP InterlockAroma::mover_name()
// XRW TEMP {
// XRW TEMP  return "InterlockAroma";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
InterlockAroma::InterlockAroma() :
	Mover( "InterlockAroma" ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::SOFT_REP_WTS ) ),
	input_ss_( "" ),
	max_repulsion_energy_( 30.0 ),
	min_env_energy_( -0.5 ),
	limit_aroma_chi2_( true ),
	output_pdbs_( false ),
	verbose_( false )
{}


/// @brief copy constructor
InterlockAroma::InterlockAroma( InterlockAroma const & ) = default;

/// @brief destructor
InterlockAroma::~InterlockAroma()= default;

/// @brief clone this object
InterlockAroma::MoverOP
InterlockAroma::clone() const
{
	return InterlockAroma::MoverOP( new InterlockAroma( *this ) );
}

/// @brief create this type of object
InterlockAroma::MoverOP
InterlockAroma::fresh_instance() const
{
	return InterlockAroma::MoverOP( new InterlockAroma() );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief mover apply
void
InterlockAroma::apply( Pose & pose )
{
	using core::pack::task::PackerTaskOP;
	using core::pack::task::TaskFactory;
	using core::scoring::dssp::Dssp;
	using core::pack::rotamer_set::Rotamers;
	using core::pack::rotamer_set::RotamerSet;
	using core::pack::rotamer_set::RotamerSetOP;
	using core::pack::rotamer_set::RotamerSetFactory;
	using protocols::toolbox::task_operations::LimitAromaChi2Operation;
	using protocols::fldsgn::filters::InterlockingAromaFilter;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::simple_moves::MakePolyXMover;

	// set pose to fullatom
	runtime_assert( pose.is_fullatom() );

	// set secondary structure
	String ss("");
	if ( ! input_ss_.empty() ) {
		ss = input_ss_;
		runtime_assert( input_ss_.length() == pose.size() );
	} else {
		Dssp dssp( pose );
		ss = dssp.get_dssp_secstruct();
	}
	SS_Info2_OP ssinfo( new SS_Info2( pose, ss ) );

	// create packer task
	TaskFactory taskf;
	if ( limit_aroma_chi2_ ) {  // default is true
		taskf.push_back( core::pack::task::operation::TaskOperationOP( new LimitAromaChi2Operation ) );
	}
	PackerTaskOP ptask( taskf.create_task_and_apply_taskoperations( pose ) );

	// debug
	std::ostringstream filename;

	// create poly-Ala pose
	MakePolyXMover bap( "ALA", false/*kee_pro*/, true /*keep_gly*/, false /*keep_disulfide_cys*/ );
	Pose polyala_pose( pose );
	bap.apply( polyala_pose );

	// farep score of original pose
	scorefxn_->set_weight( core::scoring::envsmooth, 1.0 );

	(*scorefxn_)( polyala_pose );
	Real farepE_ref = polyala_pose.energies().total_energies()[ core::scoring::fa_rep ];
	Real faatrE_ref = polyala_pose.energies().total_energies()[ core::scoring::fa_atr ];
	Real envE_ref   = polyala_pose.energies().total_energies()[ core::scoring::envsmooth ];

	// set search residue types to packer task
	InterlockingAromaFilter ilfilter;
	for ( Size ii=1; ii<=polyala_pose.size(); ++ii ) {

		// skip residues of strands
		if ( ss.at( ii-1 ) != 'L' && ss.at( ii-1 ) != 'H' ) continue;

		// set packer task
		utility::vector1<bool> restrict_to_aa( core::chemical::num_canonical_aas, false );
		restrict_to_aa[ core::chemical::aa_from_name( "PHE" ) ] = true;
		ptask->nonconst_residue_task( ii ).restrict_absent_canonical_aas( restrict_to_aa );
		ptask->set_bump_check( true );

		// create rotamer set
		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( polyala_pose.residue( ii ) );
		rotset->set_resid( ii );
		rotset->build_rotamers( polyala_pose, *scorefxn_, *ptask, utility::graph::GraphCOP( utility::graph::GraphOP( new utility::graph::Graph( polyala_pose.size() ) ) ), false );

		Size rotnum( 0 );
		for ( auto const & rotamer : *rotset ) {

			rotnum++;
			Pose work_pose( polyala_pose );

			// farep score after the insertion of rotamer
			work_pose.replace_residue( ii, *rotamer, false/*orient bb*/ );

			(*scorefxn_)( work_pose );
			Real farepE_mut = work_pose.energies().total_energies()[ core::scoring::fa_rep ];
			Real faatrE_mut = work_pose.energies().total_energies()[ core::scoring::fa_atr ];
			Real envE_mut = work_pose.energies().total_energies()[ core::scoring::envsmooth ];

			// evaluate clash
			Real dfarep = farepE_mut - farepE_ref;
			if ( dfarep > max_repulsion_energy_ ) continue;

			// evaluate buriedness
			Real denvE = envE_mut - envE_ref;
			if ( denvE > min_env_energy_ ) continue;

			// evaluate vdw attraction
			Real dfaatrE = faatrE_mut - faatrE_ref;

			if ( ilfilter.compute( ii, work_pose, ssinfo ) ) {

				if ( output_pdbs_ ) {
					filename.str("");
					filename << "interlock_" << ii << "." << rotnum << ".pdb";
					work_pose.dump_pdb( filename.str() );
				}

				if ( verbose_ ) {
					TR << ii << " " << rotnum << " " << dfarep << " " << denvE << " " << dfaatrE << std::endl;
				}

			}
		}
	}


} // InterlockAroma::apply

// XRW TEMP std::string
// XRW TEMP InterlockAroma::get_name() const {
// XRW TEMP  return InterlockAroma::mover_name();
// XRW TEMP }

/// @brief parse xml
void
InterlockAroma::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	String const sfxn ( tag->getOption<String>( "scorefxn", "soft_rep" ));
	scorefxn_ = data.get_ptr<ScoreFunction>( "scorefxns", sfxn );

	max_repulsion_energy_ = tag->getOption<Real>( "max_repulsion_energy", 30.0 );
	min_env_energy_ = tag->getOption<Real>( "min_env_energy", -0.5 );

	TR << "max_repulsion_energy: " << max_repulsion_energy_ << ", min_env_energy: " << min_env_energy_ <<  std::endl;

	// read secondary structure info through blueprint
	std::string const blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	if ( blueprint != "" ) {
		protocols::parser::BluePrintOP bop( new protocols::parser::BluePrint( blueprint ) );
		input_ss_ = bop->secstruct();
	}

	// Exclude aromatic chi2 rotamers, of which angles are around 0
	limit_aroma_chi2_ = tag->getOption<bool>( "limit_aroma_chi2", true );

	output_pdbs_ = tag->getOption<bool>( "output_pdbs", false );
	verbose_     = tag->getOption<bool>( "verbose", false );

}

std::string InterlockAroma::get_name() const {
	return mover_name();
}

std::string InterlockAroma::mover_name() {
	return "InterlockAroma";
}

void InterlockAroma::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "scorefxn", xs_string, "XRW TO DO", "soft_rep" )
		+ XMLSchemaAttribute::attribute_w_default( "max_repulsion_energy", xsct_real, "XRW TO DO", "30.0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_env_energy", xsct_real, "XRW TO DO", "-0.5" )
		+ XMLSchemaAttribute( "blueprint", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "limit_aroma_chi2", xsct_rosetta_bool, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "output_pdbs", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "XRW TO DO", "0" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string InterlockAromaCreator::keyname() const {
	return InterlockAroma::mover_name();
}

protocols::moves::MoverOP
InterlockAromaCreator::create_mover() const {
	return protocols::moves::MoverOP( new InterlockAroma );
}

void InterlockAromaCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterlockAroma::provide_xml_schema( xsd );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace of flxbb
} // namespace of protocols
