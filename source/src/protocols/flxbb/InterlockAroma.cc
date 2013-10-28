// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/flxbb/InterlockAroma.cc
/// @brief place interlocking aromatic residues
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#include <protocols/flxbb/InterlockAroma.hh>
#include <protocols/flxbb/InterlockAromaCreator.hh>

#include <protocols/fldsgn/filters/InterlockingAromaFilter.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>
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
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Parser headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

// AUTO-REMOVED #include <protocols/fldsgn/topology/HSSTriplet.hh>

#include <utility/vector0.hh>

//Auto Headers
#include <core/graph/Graph.hh>
static basic::Tracer TR("protocols.flxbb.InterlockAroma");

namespace protocols{
namespace flxbb{

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
InterlockAromaCreator::keyname() const
{
	return InterlockAromaCreator::mover_name();
}

protocols::moves::MoverOP
InterlockAromaCreator::create_mover() const {
	return new InterlockAroma;
}

std::string
InterlockAromaCreator::mover_name()
{
	return "InterlockAroma";
}

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
InterlockAroma::InterlockAroma( InterlockAroma const & rval ) :
  Super( rval ),
	scorefxn_( rval.scorefxn_ ),
 	input_ss_ ( rval.input_ss_ ),
	max_repulsion_energy_( rval.max_repulsion_energy_ ),
	min_env_energy_( rval.min_env_energy_ ),
	limit_aroma_chi2_( rval.limit_aroma_chi2_ ),
	output_pdbs_( rval.output_pdbs_ ),
	verbose_( rval.verbose_ )
{}

/// @brief destructor
InterlockAroma::~InterlockAroma(){}

/// @brief clone this object
InterlockAroma::MoverOP
InterlockAroma::clone() const
{
	return new InterlockAroma( *this );
}

/// @brief create this type of object
InterlockAroma::MoverOP
InterlockAroma::fresh_instance() const
{
	return new InterlockAroma();
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
	if( ! input_ss_.empty() ) {
		ss = input_ss_;
		runtime_assert( input_ss_.length() == pose.total_residue() );
	} else {
		Dssp dssp( pose );
		ss = dssp.get_dssp_secstruct();
	}
	SS_Info2_OP ssinfo = new SS_Info2( pose, ss );

	// create packer task
	TaskFactory taskf;
	if( limit_aroma_chi2_ ) {  // default is true
		taskf.push_back( new LimitAromaChi2Operation );
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
	for( Size ii=1; ii<=polyala_pose.total_residue(); ++ii ) {

		// skip residues of strands
		if( ss.at( ii-1 ) != 'L' && ss.at( ii-1 ) != 'H' ) continue;

		// set packer task
		utility::vector1<bool> restrict_to_aa( core::chemical::num_canonical_aas, false );
		restrict_to_aa[ core::chemical::aa_from_name( "PHE" ) ] = true;
		ptask->nonconst_residue_task( ii ).restrict_absent_canonical_aas( restrict_to_aa );
		ptask->set_bump_check( true );

		// create rotamer set
		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( polyala_pose.residue( ii ) );
		rotset->set_resid( ii );
		rotset->build_rotamers( polyala_pose, *scorefxn_, *ptask, new core::graph::Graph( polyala_pose.total_residue() ), false );

		Size rotnum( 0 );
		for( Rotamers::const_iterator rotamer = rotset->begin(); rotamer != rotset->end(); ++rotamer ) {

			rotnum++;
			Pose work_pose( polyala_pose );

			// farep score after the insertion of rotamer
			work_pose.replace_residue( ii, **rotamer, false/*orient bb*/ );

			(*scorefxn_)( work_pose );
			Real farepE_mut = work_pose.energies().total_energies()[ core::scoring::fa_rep ];
			Real faatrE_mut = work_pose.energies().total_energies()[ core::scoring::fa_atr ];
			Real envE_mut = work_pose.energies().total_energies()[ core::scoring::envsmooth ];

			// evaluate clash
			Real dfarep = farepE_mut - farepE_ref;
			if( dfarep > max_repulsion_energy_ ) continue;

			// evaluate buriedness
			Real denvE = envE_mut - envE_ref;
			if( denvE > min_env_energy_ ) continue;

			// evaluate vdw attraction
			Real dfaatrE = faatrE_mut - faatrE_ref;

			if( ilfilter.compute( ii, work_pose, ssinfo ) ) {

				if( output_pdbs_ ) {
					filename.str("");
					filename << "interlock_" << ii << "." << rotnum << ".pdb";
					work_pose.dump_pdb( filename.str() );
				}

				if( verbose_ ) {
					 TR << ii << " " << rotnum << " " << dfarep << " " << denvE << " " << dfaatrE << std::endl;
				}

			}
		}
	}


} // InterlockAroma::apply

std::string
InterlockAroma::get_name() const {
	return InterlockAromaCreator::mover_name();
}

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
	scorefxn_ = data.get< ScoreFunction * >( "scorefxns", sfxn );

	max_repulsion_energy_ = tag->getOption<Real>( "max_repulsion_energy", 30.0 );
	min_env_energy_ = tag->getOption<Real>( "min_env_energy", -0.5 );

	TR << "max_repulsion_energy: " << max_repulsion_energy_ << ", min_env_energy: " << min_env_energy_ <<  std::endl;

	// read secondary structure info through blueprint
 	std::string const blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	if( blueprint != "" ){
		protocols::jd2::parser::BluePrintOP bop = new protocols::jd2::parser::BluePrint( blueprint );
		input_ss_ = bop->secstruct();
	}

	// Exclude aromatic chi2 rotamers, of which angles are around 0
	limit_aroma_chi2_ = tag->getOption<bool>( "limit_aroma_chi2", 1 );

	output_pdbs_ = tag->getOption<bool>( "output_pdbs", 0 );
	verbose_     = tag->getOption<bool>( "verbose", 0 );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace of flxbb
} // namespace of protocols
