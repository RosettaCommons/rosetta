// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/ShapeComplementarityFilter.cc
/// @brief  Filter structures by shape complementarity and/or interface area
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

// Unit Headers
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>
#include <protocols/simple_filters/ShapeComplementarityFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/util.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/jd2/util.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static basic::Tracer tr( "protocols.simple_filters.ShapeComplementarityFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
ShapeComplementarityFilter::ShapeComplementarityFilter():
	Filter( "ShapeComplementarity" ),
	filtered_sc_( 0.50 ),
	filtered_area_( 250 ),
	filtered_d_median_( 1000.0f ),  // off by default
	jump_id_( 1 ),
	quick_( false ),
	verbose_( false ),
	selector1_(),
	selector2_(),
	write_int_area_( false ),
	multicomp_( false ),
	sym_dof_name_("")
{}


// @brief constructor with arguments
ShapeComplementarityFilter::ShapeComplementarityFilter( Real const & filtered_sc, Real const & filtered_area,
	Size const & jump_id, Size const & quick, Size const & verbose, Real const & filtered_median_distance /*= 1000.0f*/):
	Filter( "ShapeComplementarity" ),
	filtered_sc_( filtered_sc ),
	filtered_area_( filtered_area ),
	filtered_d_median_( filtered_median_distance ),
	jump_id_( jump_id ),
	quick_( quick ),
	verbose_( verbose ),
	selector1_(),
	selector2_(),
	sym_dof_name_("")
{}

void ShapeComplementarityFilter::filtered_sc( Real const & filtered_sc ) { filtered_sc_ = filtered_sc; }
void ShapeComplementarityFilter::filtered_area( Real const & filtered_area ) { filtered_area_ = filtered_area; }
void ShapeComplementarityFilter::filtered_median_distance( Real const & filtered_median_distance ) { filtered_d_median_ = filtered_median_distance; }
void ShapeComplementarityFilter::jump_id( Size const & jump_id ) { jump_id_ = jump_id; }
void ShapeComplementarityFilter::quick( Size const & quick ) { quick_ = quick; }
void ShapeComplementarityFilter::verbose( Size const & verbose ) { verbose_ = verbose; }

void ShapeComplementarityFilter::residues1( std::string const & res_string )
{
	using core::select::residue_selector::ResidueIndexSelector;
	using core::select::residue_selector::ResidueIndexSelectorOP;
	selector1_ = ResidueIndexSelectorOP( new ResidueIndexSelector( res_string ) );
}

void ShapeComplementarityFilter::residues2( std::string const & res_string )
{
	using core::select::residue_selector::ResidueIndexSelector;
	using core::select::residue_selector::ResidueIndexSelectorOP;
	selector2_ = ResidueIndexSelectorOP( new ResidueIndexSelector( res_string ) );
}


void ShapeComplementarityFilter::sym_dof_name( std::string const & sym_dof_name ) { sym_dof_name_ = sym_dof_name; }
std::string ShapeComplementarityFilter::sym_dof_name() const { return sym_dof_name_; }
void ShapeComplementarityFilter::write_int_area( bool write_int_area ) { write_int_area_ = write_int_area; }
bool ShapeComplementarityFilter::write_int_area() const { return write_int_area_; }
void ShapeComplementarityFilter::write_median_distance( bool write_median_distance ) { write_d_median_ = write_median_distance; }
bool ShapeComplementarityFilter::write_median_distance() const { return write_d_median_; }
void ShapeComplementarityFilter::multicomp( bool multicomp ) { multicomp_ = multicomp; }
bool ShapeComplementarityFilter::multicomp() const { return multicomp_; }


/// @brief
ShapeComplementarityFilter::ShapeComplementarityCalculatorResults
ShapeComplementarityFilter::compute( Pose const & pose ) const
{
	ShapeComplementarityCalculator scc;

	if ( !scc.Init() ) {
		throw CREATE_EXCEPTION(EXCN_InitFailed, "");
	}

	if ( quick_ ) {
		scc.settings.density = 5.0;
	}
	scc.Reset(); // this may not be needed anymore, but I'm leaving it here for safety

	bool symm = core::pose::symmetry::is_symmetric( pose );
	core::Real nsubs_scalefactor = 1.0;

	if ( selector1_ && selector2_ ) {
		// selector-based
		setup_from_selectors( pose, scc );
		if ( !scc.Calc() ) {
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}
	} else if ( !symm ) {
		// jump-based
		if ( !scc.Calc( pose, jump_id_ ) ) {
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}
	} else if ( multicomp_ ) {
		// MULTI COMPONENT SYMM
		setup_multi_component_symm( pose, scc, nsubs_scalefactor );
		if ( !scc.Calc() ) {
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}
	} else {
		// SINGLE COMPONENT SYMM
		setup_single_component_symm( pose, scc, nsubs_scalefactor );
		if ( !scc.Calc() ) {
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}
	}

	ShapeComplementarityCalculatorResults const & r = scc.GetResults();
	if ( verbose_ ) print_sc_results( tr, r, nsubs_scalefactor );

	return r;
}

/// @brief writes area value to current jd2 job
/// @param[in] pose     Pose being analyzed
/// @param[in] area_val Area to be reported, before correcting for symmetry.  If area < 0,
///                     the uncorrected value will be reported. If the pose isn't symmetric,
///                     the uncorrected value will be reported.
void
ShapeComplementarityFilter::write_area( Pose const & pose, core::Real const area_val ) const
{
	std::string column_header = this->get_user_defined_name() + "_int_area";

	if ( area_val < 0.0 ) {
		protocols::jd2::add_string_real_pair_to_current_job( column_header, area_val );
		return;
	}

	core::Real int_area = area_val;
	// symmetric scalefactor
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		if ( multicomp_ ) {
			utility::vector1<Size> subs = core::pose::symmetry::get_jump_name_to_subunits( pose, sym_dof_name() );
			int_area /= (Real) subs.size() ;
		} else {
			ObjexxFCL::FArray1D_bool is_upstream( pose.size(), false );
			utility::vector1<Size> sym_aware_jump_ids;
			if ( sym_dof_name() != "" ) {
				sym_aware_jump_ids.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() ) );
			} else {
				Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
				for ( Size j = 1; j <= nslidedofs; j++ ) sym_aware_jump_ids.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
			}
			core::pose::symmetry::partition_by_symm_jumps( sym_aware_jump_ids, pose.fold_tree(), core::pose::symmetry::symmetry_info(pose), is_upstream );
			Size ndownstream=0;
			for ( Size i=1; i<=pose.size(); ++i ) {
				if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
				if ( !is_upstream(i) ) ndownstream++;
			}
			int_area /= (Real)( ndownstream / core::pose::symmetry::symmetry_info(pose)->get_nres_subunit() );
		}
	}
	protocols::jd2::add_string_real_pair_to_current_job( column_header, int_area );
}

/// @brief writes median distance value to current jd2 job
/// @param[in] pose     Pose being analyzed
/// @param[in] d_median Median distance to be reported
void
ShapeComplementarityFilter::write_median_distance( Pose const &, core::Real const d_median ) const
{
	std::string column_header = this->get_user_defined_name() + "_median_dist";
	protocols::jd2::add_string_real_pair_to_current_job( column_header, d_median );
}


/// @brief
core::Real
ShapeComplementarityFilter::report_sm( Pose const & pose ) const
{
	ShapeComplementarityCalculatorResults r;
	try {
		r = compute( pose );
	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing shape complementarity calculator - returning -2 instead." << std::endl;
		if ( write_int_area_ ) write_area( pose, -2 );
		if ( write_d_median_ ) write_median_distance( pose, -2 );
		return -2;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running shape complementarity calculator - returning -1 instead." << std::endl;
		if ( write_int_area_ ) write_area( pose, -1 );
		if ( write_d_median_ ) write_median_distance( pose, -1 );
		return -1;
	}

	if ( write_int_area_ ) write_area( pose, r.area );
	if ( write_d_median_ ) write_median_distance( pose, r.distance );
	return r.sc;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose has high enough shape
// complementarity.
bool ShapeComplementarityFilter::apply( Pose const & pose ) const
{
	ShapeComplementarityCalculatorResults r;
	try {
		r = compute( pose );
	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing shape complementarity calculator - failing filter." << std::endl;
		return false;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running shape complementarity calculator - failing filter." << std::endl;
		return false;
	}

	Real const sc = r.sc;
	Real const area = r.area;
	Real const d_median = r.distance;

	if ( sc < filtered_sc_ ) {
		tr << "Filter failed current < threshold sc: " << sc << " < " << filtered_sc_ << std::endl;
		return false;
	}

	if ( area < filtered_area_ ) {
		tr << "Filter failed current < threshold interface area: " << area << " < " << filtered_area_ << std::endl;
		return false;
	}

	if ( d_median < filtered_d_median_ ) {
		tr << "Filter failed current > threshold median distance: " << d_median << " < " << filtered_d_median_ << std::endl;
		return false;
	}

	tr << "Successfully filtered: " << sc << std::endl;
	return true;
} // apply_filter

/// @brief parse xml
void
ShapeComplementarityFilter::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	filtered_sc_ = tag->getOption<Real>( "min_sc", 0.50 );
	filtered_area_ = tag->getOption<Real>( "min_interface", 0 );
	filtered_d_median_ = tag->getOption<Real>( "max_median_dist", 0 );
	verbose_ = tag->getOption<Size>( "verbose", false );
	quick_ = tag->getOption<Size>( "quick", false );
	jump_id_ = tag->getOption<Size>( "jump", 1 );
	write_int_area_ = tag->getOption<bool>( "write_int_area", false );
	write_d_median_ = tag->getOption<bool>( "write_median_dist", false );
	sym_dof_name(tag->getOption<std::string>( "sym_dof_name", "" ));
	multicomp( tag->getOption< bool >("multicomp", false) );

	if ( tag->hasOption("residues1") ) {
		residues1( tag->getOption< std::string >( "residues1" ) );
		if ( ! selector1_ ) {
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues1") << ". Using default." << std::endl;
		}
	}

	if ( tag->hasOption("residues2") ) {
		residues2( tag->getOption< std::string >( "residues2" ) );
		if ( ! selector2_ ) {
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues2") << ". Using default." << std::endl;
		}
	}

	std::string const selector1name = tag->getOption< std::string >( "residue_selector1", "" );
	if ( !selector1name.empty() ) selector1_ = core::select::residue_selector::get_residue_selector( selector1name, data );

	std::string const selector2name = tag->getOption< std::string >( "residue_selector2", "" );
	if ( !selector2name.empty() ) selector2_ = core::select::residue_selector::get_residue_selector( selector2name, data );

	tr.Info << "Structures with shape complementarity < " << filtered_sc_ << ", interface area < " <<
		filtered_area_ << " A^2, median distance > " << filtered_d_median_ << " will be filtered." << std::endl;

	if ( quick_ ) {
		tr.Info << "Calculating shape complementarity in quick mode with less accuracy." << std::endl;
	}

	if ( !selector1_ || !selector2_ ) {
		tr.Info << "Ignoring residue range selection since residues" << (selector2_ ? 1 : 2) << " is empty." << std::endl;
	}
	if ( jump_id_ != 1 ) {
		tr.Info << "Using Jump ID " << jump_id_ << " to define surfaces." << std::endl;
	}
}

/// @brief Uses residue selectors to set up the ShapeComplementarityCalculator
/// @param[in]  pose Pose to be analyzed
/// @param[out] scc Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
void
ShapeComplementarityFilter::setup_from_selectors( Pose const & pose, ShapeComplementarityCalculator & scc ) const
{
	using core::select::residue_selector::ResidueVector;

	ResidueVector const residues1( selector1_->apply( pose ) );
	ResidueVector const residues2( selector2_->apply( pose ) );

	// Dump information about residues
	if ( tr.Info.visible() ) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		for ( auto r = residues1.begin(); r != residues1.end(); ++r ) {
			tr.Info << (r == residues1.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		for ( auto r = residues2.begin(); r != residues2.end(); ++r ) {
			tr.Info << (r == residues2.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
	}

	for ( core::Size r : residues1 ) {
		scc.AddResidue( 0, pose.residue( r ) );
	}

	for ( core::Size r : residues2 ) {
		scc.AddResidue( 1, pose.residue( r ) );
	}
}

/// @brief Uses multi-component symmetric interfaces to set up the ShapeComplementarityCalculator
/// @param[in]  pose              Pose to be analyzed
/// @param[out] scc               Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
/// @param[out] nsubs_scalefactor Writes number of subunits, to be used as scaling factor for sc calculations
void
ShapeComplementarityFilter::setup_multi_component_symm(
	Pose const & pose,
	ShapeComplementarityCalculator & scc,
	core::Real & nsubs_scalefactor ) const
{
	runtime_assert( sym_dof_name() != "" );
	utility::vector1<std::string> sym_dof_name_list = utility::string_split( sym_dof_name() , ',' );

	Size sym_dof_index = 1;

	if ( sym_dof_name_list.size() > 1 ) {
		Size intracontact_count = 0;
		for ( core::Size i=1; i<=sym_dof_name_list.size(); i++ ) {
			if ( core::pose::symmetry::intracomponent_contact(pose, sym_dof_name_list[i], 12.0) ) {
				intracontact_count++;
				sym_dof_index = i;
			}
		}
		if ( intracontact_count > 1 ) {
			tr.Warning << "Intracontacts detected between multiple components.  Calculating sc based off of the last component with intracontacts.  Separate sc calculations are recommended for each component with intracontacts." << std::endl;
		}
	}

	utility::vector1<Size> full_intracomp_resis = core::pose::symmetry::get_full_intracomponent_resis(pose, sym_dof_name_list[sym_dof_index]);
	utility::vector1<Size> full_intracomp_neighbor_resis = core::pose::symmetry::get_full_intracomponent_neighbor_resis(pose, sym_dof_name_list[sym_dof_index], 12.0 );

	//core::pose::Pose full_intracomp_subpose = core::pose::symmetry::get_full_intracomponent_subpose(pose, sym_dof_name_list[sym_dof_index]);
	//core::pose::Pose full_intracomp_neighbor_subpose = core::pose::symmetry::get_full_intracomponent_neighbor_subpose(pose, sym_dof_name_list[sym_dof_index], 12.0);
	//full_intracomp_subpose.dump_pdb("full_intracomp_subpose_" + protocols::jd2::current_output_name() + ".pdb");
	//full_intracomp_neighbor_subpose.dump_pdb("full_intracomp_neighbor_subpose_" + protocols::jd2::current_output_name() + ".pdb");

	for ( core::Size i=1; i<=full_intracomp_resis.size(); i++ ) {
		scc.AddResidue(0,pose.residue(full_intracomp_resis[i]));
	}
	for ( core::Size i=1; i<=full_intracomp_neighbor_resis.size(); i++ ) {
		scc.AddResidue(1,pose.residue(full_intracomp_neighbor_resis[i]));
	}
	utility::vector1<Size> const subs = core::pose::symmetry::get_jump_name_to_subunits( pose, sym_dof_name_list[sym_dof_index] );
	nsubs_scalefactor = subs.size() ;
}

/// @brief Uses single-component symmetric interfaces to set up the ShapeComplementarityCalculator
/// @param[in]  pose              Pose to be analyzed
/// @param[out] scc               Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
/// @param[out] nsubs_scalefactor Writes number of subunits, to be used as scaling factor for sc calculations
void
ShapeComplementarityFilter::setup_single_component_symm(
	Pose const & pose,
	ShapeComplementarityCalculator & scc,
	core::Real & nsubs_scalefactor ) const
{
	ObjexxFCL::FArray1D_bool is_upstream ( pose.size(), false );
	utility::vector1<Size> sym_aware_jump_ids;

	if ( sym_dof_name() != "" ) {
		sym_aware_jump_ids.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() ) );
	} else {
		// all slidable jumps
		Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
		for ( Size j = 1; j <= nslidedofs; j++ ) {
			sym_aware_jump_ids.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
		}
	}

	// partition & fill residueX_ vectors
	core::pose::symmetry::partition_by_symm_jumps( sym_aware_jump_ids, pose.fold_tree(), core::pose::symmetry::symmetry_info(pose), is_upstream );
	Size ndownstream=0;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		scc.AddResidue(is_upstream(i)?1:0, pose.residue(i));
		if ( !is_upstream(i) ) ndownstream++;
	}

	// scalefactor
	nsubs_scalefactor = core::Real(ndownstream)/ core::Real(core::pose::symmetry::symmetry_info(pose)->get_nres_subunit());
}

/// @brief prints results to given tracer in a human-readable format
/// @param[out] tr std::ostream object to write to
/// @param[in]  r  ShapeComplementarityCalculatorResults object containing results
void
ShapeComplementarityFilter::print_sc_results(
	std::ostream & tr,
	ShapeComplementarityCalculatorResults const & r,
	core::Real const nsubs_scalefactor ) const
{
	// Verbose view
	tr << "==================================================" << std::endl;
	tr << std::endl;
	for ( int i = 0; i <= 2; i++ ) {
		if ( i < 2 ) {
			tr << "Molecule " << (i+1) << ":" << std::endl;
		} else {
			tr << "Total/Average for both molecules:" << std::endl;
		}

		tr << "          Total Atoms: " << r.surface[i].nAtoms << std::endl;
		tr << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << std::endl;
		tr << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << std::endl;
		tr << "           Total Dots: " << r.surface[i].nAllDots << std::endl;
		tr << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << std::endl;
		tr << "         Trimmed Area: " << r.surface[i].trimmedArea << " (avg) " << std::endl;
		tr << std::endl;
	}
	tr << std::endl;

	for ( int i = 0; i <= 2; i++ ) {
		if ( i < 2 ) {
			tr << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << std::endl;
		} else {
			tr << "Average for both molecules:" << std::endl;
		}
		tr << "      Mean Separation: " << r.surface[i].d_mean << std::endl;
		tr << "    Median Separation: " << r.surface[i].d_median << std::endl;
		tr << "    Mean Shape Compl.: " << r.surface[i].s_mean << std::endl;
		tr << "  Median Shape Compl.: " << r.surface[i].s_median << std::endl;
		tr << std::endl;
	}

	tr << "Shape complementarity: " << r.sc << std::endl;
	tr << "Interface area: " << r.area << std::endl;
	if ( nsubs_scalefactor != 1 ) {
		tr << "Area per monomer: " << ( (core::Real) r.area / nsubs_scalefactor ) << std::endl ;
	}
	tr << "Interface seperation: " << r.distance << std::endl;
}

// XRW TEMP filters::FilterOP
// XRW TEMP ShapeComplementarityFilterCreator::create_filter() const { return filters::FilterOP( new ShapeComplementarityFilter ); }

// XRW TEMP std::string
// XRW TEMP ShapeComplementarityFilterCreator::keyname() const { return "ShapeComplementarity"; }

std::string ShapeComplementarityFilter::name() const {
	return class_name();
}

std::string ShapeComplementarityFilter::class_name() {
	return "ShapeComplementarity";
}

void ShapeComplementarityFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "min_sc" , xsct_real , "The filter fails if the calculated sc is less than the given value." , "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "min_interface" , xsct_real , "The filter fails is the calculated interface area is less than the given value." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_median_dist" , xsct_real , "The filter fails is the calculated median distance between the molecular surfaces is greater than the given value." , "1000" )
		// The default value for min_interface is strange. I feel like it should be '9999' so that things do not automatically fail if the user does not specificy the filter.
		+ XMLSchemaAttribute::attribute_w_default( "verbose" , xsct_rosetta_bool , "If true, print extra calculation details to the tracer." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "quick" , xsct_rosetta_bool , "If true, do a quicker, less accurate calculation by reducing the density." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "jump" , xs_integer , "For non-symmetric poses, which jump over which to calculate the interface." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "write_int_area" , xsct_rosetta_bool , "If true, write interface area to scorefile." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "write_median_dist" , xsct_rosetta_bool , "If true, write interface median distance to scorefile." , "false" )
		+ XMLSchemaAttribute( "sym_dof_name" , xs_string , "For symmetric poses, which dof over which to calculate the interface." )
		+ XMLSchemaAttribute::attribute_w_default( "multicomp" , xsct_rosetta_bool , "If true, multiple component system. If false, single component system." , "false" )
		+ XMLSchemaAttribute( "residues1" , xs_string , "Explicitly set which residues are on each side of the interface (both symmetric and non-symmetric poses.)" )
		+ XMLSchemaAttribute( "residues2" , xs_string , "Explicitly set which residues are on each side of the interface (both symmetric and non-symmetric poses.)" )
		+ XMLSchemaAttribute( "residue_selector1" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." )
		+ XMLSchemaAttribute( "residue_selector2" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates the Lawrence and Coleman shape complementarity using a port of the original Fortran code from CCP4's sc. Symmetry aware. Can be calculated across a jump (default behavior) or the two surfaces can be specified by explicitly providing lists of the residues making up each surface.", attlist );
}

std::string ShapeComplementarityFilterCreator::keyname() const {
	return ShapeComplementarityFilter::class_name();
}

protocols::filters::FilterOP
ShapeComplementarityFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ShapeComplementarityFilter );
}

void ShapeComplementarityFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShapeComplementarityFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
