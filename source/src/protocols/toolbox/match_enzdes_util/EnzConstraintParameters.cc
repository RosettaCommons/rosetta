// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for enzyme design constraints
/// @author Florian Richter, floric@u.washington.edu

// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>

// Package headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh> //need for changing residue type sets
#include <core/chemical/ResidueTypeSet.hh> //have to include complete file
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Patch.hh> //needed for residue type base name function
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CharmmPeriodicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh> //need function in this file
#include <core/scoring/ScoreFunction.hh> //scoring ambiguous constraints
#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh> //needed for clean residue type modification
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Remarks.hh> //reading remarks
#include <core/id/SequenceMapping.hh>

// Basic headers
#include <basic/options/option.hh> //options
#include <basic/Tracer.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// External header
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <iostream>
#include <string>
#include <sstream>


static basic::Tracer tr("protocols.toolbox.match_enzdes_util.EnzConstraintParameters");

namespace protocols {
namespace toolbox{
namespace match_enzdes_util {

using namespace ObjexxFCL;

CovalentConnectionReplaceInfo::CovalentConnectionReplaceInfo(
	std::string resA_base_in,
 	std::string resB_base_in,
 	std::string resA_var_in,
 	std::string resB_var_in,
 	core::Size Apos_in,
	core::Size Bpos_in,
 	core::chemical::ResidueTypeSetCAP restype_set_in
) : ReferenceCount(),
		resA_basename_(resA_base_in), resB_basename_(resB_base_in),
		resA_varname_(resA_var_in), resB_varname_(resB_var_in),
		resA_modname_(resA_base_in + resA_var_in), resB_modname_(resB_base_in+resB_var_in),
		resA_seqpos_(Apos_in), resB_seqpos_(Bpos_in),
		restype_set_(restype_set_in)
{}


CovalentConnectionReplaceInfo::CovalentConnectionReplaceInfo( CovalentConnectionReplaceInfo const & other )
	: ReferenceCount(),
		resA_basename_(other.resA_basename_), resB_basename_(other.resB_basename_),
		resA_varname_(other.resA_varname_), resB_varname_(other.resB_varname_),
		resA_modname_(other.resA_modname_), resB_modname_(other.resB_modname_),
		resA_seqpos_(other.resA_seqpos_), resB_seqpos_(other.resB_seqpos_),
		restype_set_(other.restype_set_)
{}

CovalentConnectionReplaceInfo::~CovalentConnectionReplaceInfo(){}

void
CovalentConnectionReplaceInfo::remove_covalent_connection_from_pose(
	core::pose::Pose & pose
) const {

	//first check whether the right residue types are still at their positions
	if( residue_type_base_name( pose.residue_type(resA_seqpos_) ) != resA_modname_ ){
			utility_exit_with_message("Error when trying to remove covalent connection: for resA, pose residue basename should be "+resA_modname_+", but is "+residue_type_base_name( pose.residue_type(resA_seqpos_))+"\n.");
		}

	if( residue_type_base_name( pose.residue_type(resB_seqpos_) ) != resB_modname_ ){
			utility_exit_with_message("Error when trying to remove covalent connection: for resB, pose residue basename should be "+resB_modname_+", but is "+residue_type_base_name( pose.residue_type(resB_seqpos_))+"\n.");
		}

		core::conformation::Residue newA_res( restype_set_->name_map(resA_basename_), true);

		utility::vector1< std::string > curA_variants =
				pose.residue_type( resA_seqpos_ ).properties().get_list_of_variants();
		replace_residue_keeping_all_atom_positions( pose, newA_res, resA_seqpos_ );

		for ( core::Size var = 1; var <= curA_variants.size(); ++var ) {
			if( curA_variants[ var ] != resA_varname_ ){
				core::pose::add_variant_type_to_pose_residue(
						pose,
						core::chemical::ResidueProperties::get_variant_from_string( curA_variants[ var ] ),
						resA_seqpos_ );
			}
		}

		core::conformation::Residue newB_res( restype_set_->name_map(resB_basename_), true);

		utility::vector1< std::string > curB_variants =
				pose.residue_type( resB_seqpos_ ).properties().get_list_of_variants();
		replace_residue_keeping_all_atom_positions( pose, newB_res, resB_seqpos_ );

		for ( core::Size var = 1; var <= curB_variants.size(); ++var ) {
			if( curB_variants[ var ] != resB_varname_ ){
				core::pose::add_variant_type_to_pose_residue(
						pose,
						core::chemical::ResidueProperties::get_variant_from_string( curB_variants[ var ] ),
						resB_seqpos_ );
			}
		}
		//std::cerr << "done removing covalent connection between res " << resA_seqpos_ << " and res " << resB_seqpos_ << std::endl;
}


void
CovalentConnectionReplaceInfo::remap_resid(
	core::id::SequenceMapping const & smap
){
	resA_seqpos_ = smap[resA_seqpos_];
	resB_seqpos_ = smap[resB_seqpos_];

	if( (resA_seqpos_ == 0) || (resB_seqpos_ == 0 ) ) utility_exit_with_message("A covalently constrained residue apparently got deleted from the pose");
}

EnzConstraintParameters::EnzConstraintParameters(
	core::Size cst_block,
	core::chemical::ResidueTypeSetCAP src_restype_set,
	EnzConstraintIOCAP src_enz_io
)
: utility::pointer::ReferenceCount(),
	mcfi_(NULL),
	ndisAB_(0), nangleA_(0), nangleB_(0), ntorsionA_(0), ntorsionB_(0), ntorsionAB_(0),
	restype_set_( src_restype_set), enz_io_( src_enz_io ), cst_block_(cst_block)
{
	resA_ = new EnzCstTemplateRes(src_restype_set, this);
	resA_->set_param_index( 1 );
	resB_ = new EnzCstTemplateRes(src_restype_set, this);
	resB_->set_param_index( 2 );
	disAB_ = NULL; angleA_ = NULL; angleB_ = NULL; torsionA_ = NULL; torsionB_ = NULL; torsionAB_ = NULL;
	is_covalent_ = false;
	empty_ = true;
}
EnzConstraintParameters::EnzConstraintParameters()
: utility::pointer::ReferenceCount(),
	resA_(NULL), resB_(NULL), mcfi_(NULL),disAB_(NULL),
	angleA_(NULL), angleB_(NULL), torsionA_(NULL),
	torsionB_(NULL), torsionAB_(NULL),
	ndisAB_(0), nangleA_(0), nangleB_(0), ntorsionA_(0),
	ntorsionB_(0), ntorsionAB_(0),
	is_covalent_(false), empty_(true),
	restype_set_( NULL), enz_io_( NULL), cst_block_(0)
{}
EnzConstraintParameters::~EnzConstraintParameters(){}

/// @brief copy constructor
/// @brief WARNING: currently this probably doesn't copy the functions or active pose constraints
EnzConstraintParameters::EnzConstraintParameters( EnzConstraintParameters const & other )
:
	utility::pointer::ReferenceCount(),
	mcfi_(other.mcfi_),
	disAB_(other.disAB_), angleA_(other.angleA_), angleB_(other.angleB_),
	torsionA_(other.torsionA_), torsionB_(other.torsionB_), torsionAB_(other.torsionAB_),
	ndisAB_(other.ndisAB_), nangleA_(other.nangleA_), nangleB_(other.nangleB_),
	ntorsionA_(other.ntorsionA_), ntorsionB_(other.ntorsionB_), ntorsionAB_(other.ntorsionAB_),
	is_covalent_(other.is_covalent_), empty_(other.empty_),
	restype_set_( other.restype_set_),
	cst_block_(other.cst_block_) //debatable whether this should be kept, but it's useful
{

	resA_ = new EnzCstTemplateRes( other.resA_, this );
	resB_ = new EnzCstTemplateRes( other.resB_, this );
	// we want an independent copy of the object, so the original reference is erased
	enz_io_ = NULL;
}


void
EnzConstraintParameters::set_mcfi(
 	toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi )
{
	//if this mcfi has been previously set,
	//we don't need to do anything
	if( mcfi_ == mcfi ) return;

	mcfi_ = mcfi;

	resA_ = new EnzCstTemplateRes( mcfi_->enz_cst_template_res( 1 ), this );
	resA_->set_param_index( 1 );
	resB_ = new EnzCstTemplateRes( mcfi_->enz_cst_template_res( 2 ), this );
	resB_->set_param_index( 2 );
	resA_->identical_info_consistency_check();
	resB_->identical_info_consistency_check();

	angleA_ = convert_GeomSampleInfo_to_FuncOP( mcfi_->ang_U1D2(), nangleA_ );
	angleB_ = convert_GeomSampleInfo_to_FuncOP( mcfi_->ang_U2D1(), nangleB_ );
	torsionA_ = convert_GeomSampleInfo_to_FuncOP( mcfi_->tor_U1D3(), ntorsionA_ );
	torsionAB_ = convert_GeomSampleInfo_to_FuncOP( mcfi_->tor_U2D2(), ntorsionAB_ );
	torsionB_ = convert_GeomSampleInfo_to_FuncOP( mcfi_->tor_U3D1(), ntorsionB_ );

	//the distance is special, we'll do it explicitly
	if( mcfi_->dis_U1D1() ){
		ndisAB_ = mcfi_->dis_U1D1()->ideal_val();

		core::Real min_dis = std::max(0.0, ndisAB_ - mcfi_->dis_U1D1()->tolerance() );
		core::Real max_dis = ndisAB_ + mcfi_->dis_U1D1()->tolerance();
		core::Real force_k_dis = mcfi_->dis_U1D1()->force_const();

		disAB_ = new core::scoring::constraints::BoundFunc(
			min_dis, max_dis,	sqrt(1/ force_k_dis),	"dis");

		//if( mcfi_->dis_U1D1()->periodicity() == 1.0 ) is_covalent_ = true;
		if( mcfi->is_covalent() ) is_covalent_ = true;
		else is_covalent_ = false;
	}
}


core::scoring::func::FuncOP
EnzConstraintParameters::convert_GeomSampleInfo_to_FuncOP(
	toolbox::match_enzdes_util::GeomSampleInfoCOP gsi,
	core::Real & ideal_val)
{

	core::Real const rad_per_deg = numeric::constants::f::degrees_to_radians;
	core::Real const twopi = numeric::constants::f::pi_2;

	core::scoring::func::FuncOP to_return( NULL );

	//check if this gsi has a force constant resp. if it even exists
	//if not, return right away. this will lead to no constraints being produced
	if( !gsi ) return to_return;

	ideal_val = gsi->ideal_val() * rad_per_deg;

	if( gsi->force_const() == 0.0 ) return to_return;

	core::Real x0 = gsi->ideal_val() * rad_per_deg;
	core::Real x_sd = gsi->tolerance() * rad_per_deg;
	core::Real period_rad = gsi->periodicity() * rad_per_deg;

	if( gsi->function_tag() == "PERIODIC" ){

		if( ( gsi->tag() == "angle_A:" ) || ( gsi->tag() == "angle_B:" ) ){
			std::cerr << "Error in cstfile reading: periodic functions are not supported for angle constraints." << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		if(x_sd > 0.0){
			tr.Info << "WARNING: Constraint specified for tag " << gsi->tag() << " requests a periodic function. Standard-Deviation/tolerance is meaningless in this case. The tolerance value of " << gsi->tolerance() << " that you specified will be ignored!" << std::endl;
		}

		to_return = new core::scoring::func::CharmmPeriodicFunc(x0, gsi->force_const(), twopi/period_rad);

	}

	else{
		to_return = new core::scoring::constraints::OffsetPeriodicBoundFunc(-x_sd,x_sd, sqrt( 1/gsi->force_const() ), "offsetperiodicbound", period_rad, x0);
	}

	return to_return;

}

/// @brief process the information in template residues and the func op pointers to add
/// @brief to the constraint set in the right format
void
EnzConstraintParameters::generate_active_pose_constraints(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx) const
{

	using namespace core::scoring::constraints;

	//sequentially process the six possible constraints
	bool all_constraints_empty = true;

	EnzdesCstParamCacheOP param_cache(protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );
	runtime_assert( param_cache );
	param_cache->active_pose_constraints_.clear();
	if( param_cache->covalent_connections_.size() != 0 ) remove_covalent_connections_from_pose( pose );

	Size number_constraints_added(0);

	//std::cerr << "generating constraints for block " << cst_block_ << ". There are " << resA_->respos_map_size() << " positions in resA and " << resB_->respos_map_size() << " positions in resB." << std::endl;

	for (std::map< Size, EnzCstTemplateResAtomsOP  >::const_iterator resApos_it = param_cache->template_res_cache_[1]->seqpos_map_begin(), resApos_end = param_cache->template_res_cache_[1]->seqpos_map_end(); resApos_it != resApos_end; ++resApos_it){

		for (std::map< Size, EnzCstTemplateResAtomsOP  >::const_iterator resBpos_it = param_cache->template_res_cache_[2]->seqpos_map_begin(), resBpos_end = param_cache->template_res_cache_[2]->seqpos_map_end(); resBpos_it != resBpos_end; ++resBpos_it){

			//std::cerr << "making stuff between resA " << resApos_it->first << " and resB " << resBpos_it->first << std::endl;

			Size ambig_resA(resApos_it->second->atom1_.size()), ambig_resB(resBpos_it->second->atom1_.size());

			Size number_ambiguous_constraints = ambig_resA * ambig_resB;
			if( number_ambiguous_constraints != 1){
				tr.Info << "Constraint specified between residues " << resApos_it->first << " and " << resBpos_it->first << " is ambiguous, " << number_ambiguous_constraints << " different possibilities ( " << ambig_resA << " for " << resApos_it->first << ", " << ambig_resB << " for " << resBpos_it->first << ".)" << std::endl;
			}


			utility::vector1< ConstraintCOP > ambig_csts;

			for( Size ambig_resA_count(0); ambig_resA_count < ambig_resA; ambig_resA_count++){

				for( Size ambig_resB_count(0); ambig_resB_count < ambig_resB; ambig_resB_count++){


					utility::vector1< ConstraintCOP > this_pair_csts;

					if(disAB_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new AtomPairConstraint( resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], disAB_) );
						number_constraints_added++;

						if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] )
							{
								std::string at1name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom1_[ambig_resA_count].atomno());
								std::string at2name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom1_[ambig_resB_count].atomno());
								tr << "adding distance constraint between atom " << at1name  << "of res " << resApos_it->first << " and atom " <<  at2name  << " of resi " << resBpos_it->first << "." << std::endl;
							}
						//debug stuff over
					}

					if(angleA_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new AngleConstraint( resApos_it->second->atom2_[ambig_resA_count], resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], angleA_) );
						number_constraints_added++;

						//debug stuff
						if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] )
							{
								std::string at1name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom2_[ambig_resA_count].atomno());
								std::string at2name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom1_[ambig_resA_count].atomno());
								std::string at3name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom1_[ambig_resB_count].atomno());
								tr << "adding angle constraint between atoms " << at1name  << " and " << at2name << "of res " << resApos_it->first << " and atom " <<  at3name  << " of resi " << resBpos_it->first << "." << std::endl;
							}
						//debug stuff over
					}

					if(angleB_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new AngleConstraint( resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], resBpos_it->second->atom2_[ambig_resB_count], angleB_) );
						number_constraints_added++;

						//debug stuff
						if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] )
							{
								std::string at1name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom1_[ambig_resA_count].atomno());
								std::string at2name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom1_[ambig_resB_count].atomno());
								std::string at3name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom2_[ambig_resB_count].atomno());
								tr << "adding angle constraint between atoms " << at1name  <<  "of res " << resApos_it->first << " and atom " <<  at2name  << " and " << at3name << " of resi " << resBpos_it->first << "." << std::endl;
							}
						//debug stuff over

					}

					if(torsionA_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new DihedralConstraint( resApos_it->second->atom3_[ambig_resA_count], resApos_it->second->atom2_[ambig_resA_count], resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], torsionA_) );
						number_constraints_added++;

						//debug stuff
						if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] )
							{
								std::string at1name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom3_[ambig_resA_count].atomno());
								std::string at2name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom2_[ambig_resA_count].atomno());
								std::string at3name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom1_[ambig_resA_count].atomno());
								std::string at4name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom1_[ambig_resB_count].atomno());
								tr << "adding dihedral constraint between atoms " << at1name  <<  " and " << at2name << " and " << at3name << " of res " << resApos_it->first << " and atom " <<  at4name  << " of resi " << resBpos_it->first << "." << std::endl;
							}
						//debug stuff over

					}

					if(torsionB_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new DihedralConstraint( resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], resBpos_it->second->atom2_[ambig_resB_count], resBpos_it->second->atom3_[ambig_resB_count], torsionB_) );
						number_constraints_added++;

						//debug stuff
						if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] )
							{
								std::string at1name = pose.residue(resApos_it->first).atom_name(resApos_it->second->atom1_[ambig_resA_count].atomno());
								std::string at2name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom1_[ambig_resB_count].atomno());
								std::string at3name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom2_[ambig_resB_count].atomno());
								std::string at4name = pose.residue(resBpos_it->first).atom_name(resBpos_it->second->atom3_[ambig_resB_count].atomno());
								tr << "adding dihedral constraint between atoms " << at1name  <<  " of res " << resApos_it->first << " and " << at2name << " and " << at3name << " and " <<  at4name  << " of resi " << resBpos_it->first << "." << std::endl;
							}
						//debug stuff over
					}

					if(torsionAB_ != 0){
						all_constraints_empty = false;
						this_pair_csts.push_back( new DihedralConstraint( resApos_it->second->atom2_[ambig_resA_count], resApos_it->second->atom1_[ambig_resA_count], resBpos_it->second->atom1_[ambig_resB_count], resBpos_it->second->atom2_[ambig_resB_count], torsionAB_) );
						number_constraints_added++;
					}

					if( (this_pair_csts.size() > 0 ) && ( number_ambiguous_constraints == 1) ){

						param_cache->active_pose_constraints_.push_back( new MultiConstraint(this_pair_csts) );
					}
					else if( (this_pair_csts.size() > 0) && ( number_ambiguous_constraints > 1) ) {
						ambig_csts.push_back( new MultiConstraint(this_pair_csts) );
					}
				} //ambig_resB_count
			} //ambig_resA_count

			if( ambig_csts.size() > 0 ) {
				//cstset->show( tr );

				if( is_covalent_ ) {
					//we first have to resolve the ambiguity to make the covalent connection
					//strategy: score each constraint in succession, remember which one is the lowest
					tr.Info << "The covalent constraint specified for this block is ambiguous. ";
					tr.Info << "Ambiguity will be resolved according to the best constraint in the input pose." << std::endl;
					Size n_best_constraint = determine_best_constraint( pose, scofx, ambig_csts );

					param_cache->active_pose_constraints_.push_back( ambig_csts[ n_best_constraint ] );

					//additional feat: we have to map the information about the best constraint back to what
					//atoms it contains, in the context of where these atoms are saved in the EnzCstTemplateRes object
					Size best_resA_At(0);
					Size best_resB_At(0);
					Size n_mod_ambigA = n_best_constraint % ambig_resA;
					if( n_mod_ambigA == 0) {
						best_resA_At =  ambig_resA - 1;
						best_resB_At = ( n_best_constraint / ambig_resA ) - 1;
					}
					else{
						best_resA_At = n_mod_ambigA - 1;
						best_resB_At = n_best_constraint / ambig_resA;  //integer division
					}
					make_constraint_covalent( pose, resApos_it->first, resBpos_it->first, best_resA_At, best_resB_At);

				}
				else{
					if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
						tr.Info << "Adding an ambiguous constraint containing " << ambig_csts.size() << "constraints." << std::endl;
					}
					param_cache->active_pose_constraints_.push_back( new AmbiguousConstraint(ambig_csts) );
				}
			}
			else if( is_covalent_ ){
				make_constraint_covalent( pose, resApos_it->first, resBpos_it->first, 0, 0);
			} //if is_covalent_


		}//iterator over the residue numbers for template residue B
	}//iterator over the residue numbers for template residue A

	empty_ = all_constraints_empty;

	if(all_constraints_empty) {
		tr.Info << "Warning: no constraints were added for constraint block " << cst_block_ << "." << std::endl;
	}
	else {
		tr.Info << "for block " << cst_block_ << ", " << number_constraints_added << " newly generated constraints were added " << std::endl;
	}

} //add_constraints_to_cst_set


void
EnzConstraintParameters::make_constraint_covalent(
	core::pose::Pose & pose,
	Size resA_pos,
	Size resB_pos,
	Size resA_At,
	Size resB_At) const
{
	using namespace core::chemical;

	if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
		pose.dump_pdb("bef_resmod.pdb");
	}
	std::string resA_base = residue_type_base_name( pose.residue_type(resA_pos) );
	std::string resB_base = residue_type_base_name( pose.residue_type(resB_pos) );
	std::string resA_var, resB_var;

	make_constraint_covalent_helper( pose, resA_, resA_pos, resA_At, ntorsionA_, nangleA_, disAB_, resA_var );

	make_constraint_covalent_helper( pose, resB_, resB_pos, resB_At, ntorsionB_, nangleB_, disAB_, resB_var );

	if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
		pose.dump_pdb("after_resmod.pdb");
	}

	std::string resA_atomname = pose.residue( resA_pos ).atom_name( (resA_->get_template_atoms_at_pos(pose, resA_pos))->atom1_[resA_At].atomno());
	std::string resB_atomname = pose.residue( resB_pos ).atom_name( (resB_->get_template_atoms_at_pos(pose, resB_pos))->atom1_[resB_At].atomno());

	tr.Debug << "Adding chemical bond between " << pose.residue( resA_pos ).name() << " " << (resA_->get_template_atoms_at_pos(pose, resA_pos))->atom1_[resA_At].atomno() <<  " "<< resA_pos << " " << resA_atomname << " and "
					 << pose.residue( resB_pos ).name() << " " << resB_pos << " " << (resB_->get_template_atoms_at_pos(pose, resB_pos))->atom1_[resB_At].atomno() << " " << resB_atomname << std::endl;

	pose.conformation().declare_chemical_bond(
		resA_pos, resA_atomname,
		resB_pos, resB_atomname
	);
	EnzdesCstParamCacheOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );

	param_cache->covalent_connections_.push_back( new CovalentConnectionReplaceInfo(resA_base, resB_base, resA_var, resB_var, resA_pos, resB_pos, restype_set_ ) ); //new

} //make_constraint_covalent


/// @brief helper function so stuff doesn't need to be written twice
void
EnzConstraintParameters::make_constraint_covalent_helper(
	core::pose::Pose & pose,
	EnzCstTemplateResOP template_res,
	core::Size res_pos,
	core::Size Atpos,
	core::Real itorsion,
	core::Real iangle,
	core::Real idis,
 	std::string & res_varname
) const
{
	//std::cout << "APL DEBUG EnzConstraintParameters.cc::make_constraint_covalent_helper begin" << std::endl;

	using namespace core::chemical;
	using namespace core::pack::dunbrack;

	std::string res_atom = pose.residue(res_pos).atom_name( (template_res->get_template_atoms_at_pos(pose, res_pos) )->atom1_[Atpos].atomno() );

	ObjexxFCL::strip_whitespace( res_atom );

	std::string current_residue_type_basename( residue_type_base_name( pose.residue_type(res_pos) ) );
	std::string current_residue_type_patches_name( residue_type_all_patches_name( pose.residue_type(res_pos) ) );

	res_varname = "_connect" + res_atom;  // FIXME: variant names and patch names are not the same thing! ~Labonte
	{// scope
		// Find a name for the new residue type / variant name that will be added to the existing
		// residue so that, if the existing residue already has this variant type, then the
		// new residue type will get one with a new name.
		core::chemical::ResidueType const & currres( pose.residue_type( res_pos ));
		Size count=0;
		while ( true ) {
			if ( count > 1000 ) {
				utility_exit_with_message( "Encountered infinite loop trying to find a new variant name for residue type " + currres.name() + " in EnzConstraintParameters.  Talk to Andrew.");
			}
			++count;
			if ( count == 1 ) {
				if ( ! currres.has_variant_type( res_varname )) break;
			} else {
				res_varname = "_"+utility::to_string(count)+"connect"+res_atom;
				if ( ! currres.has_variant_type( res_varname )) break;
			}
		}
	}
	std::string res_type_mod_name( current_residue_type_basename + res_varname + current_residue_type_patches_name );

	//check whether the modified residues have already been created earlier
	if( !restype_set_->has_name(res_type_mod_name) ){

		//holy jesus, we have to change the residue type set.
		//we not only need to clone, modify and add  the residue type
		//currently in the pose, but all similar ones (same basename )
		//in the ResidueTypeSet. We also need to create a copy of the
		// SingleLigandRotamerLibrary in case it exists.
		//reminds me of open heart surgery...

		//the following line is necessary to ensure that a ligand rotamer library exists
		//if this function is called before any scoring happened
		RotamerLibrary::get_instance().get_rsd_library( pose.residue_type( res_pos ));

		SingleLigandRotamerLibraryOP new_lrots = NULL;
		if( pose.residue_type(res_pos).is_ligand() &&
				RotamerLibrary::get_instance().rsd_library_already_loaded( pose.residue_type(res_pos) ) ) {
			new_lrots = new SingleLigandRotamerLibrary();
		}

		utility::pointer::access_ptr< ResidueTypeSet > mod_restype_set = & ChemicalManager::get_instance()->nonconst_residue_type_set( restype_set_->name() );

		//first get all residue types that correspond to the type in question
		ResidueTypeCOPs res_to_modify = mod_restype_set->name3_map( pose.residue_type(res_pos).name3() );

		for ( utility::vector1< ResidueTypeCOP >::iterator res_it = res_to_modify.begin(); res_it != res_to_modify.end(); ++res_it) {
			std::string const base_name( residue_type_base_name( *(*res_it) ) );
			//std::cerr << "contemplating modification of residuetype " << (*res_it)->name() << " with basename " << base_name << std::endl;

			if( current_residue_type_basename == base_name ){

				ResidueTypeOP mod_res;
				core::Size con_res(0);
				//std::cerr << " MODIFYING" << std::endl;
				std::string patches_name( residue_type_all_patches_name( *(*res_it) ) );
				std::string new_name( base_name + res_varname + patches_name );

				mod_res = (*res_it)->clone();
				con_res = mod_res->add_residue_connection( res_atom );
				mod_res->name( new_name );
				assert( ! mod_res->has_variant_type( res_varname ) );
				mod_res->enable_custom_variant_types();
				mod_res->add_variant_type( res_varname ); //necessary to restrict the packer to only use this residue variant in packing

				mod_res->set_icoor( "CONN"+string_of( con_res ), itorsion, iangle, idis, res_atom, pose.residue(res_pos).atom_name( (template_res->get_template_atoms_at_pos( pose, res_pos) )->atom2_[Atpos].atomno() ), pose.residue(res_pos).atom_name( (template_res->get_template_atoms_at_pos(pose, res_pos) )->atom3_[Atpos].atomno() ), true );

				//new_lrots is empty at the moment, but will be filled a couple of lines down
				if( pose.residue_type( res_pos ).is_ligand() ) {
					RotamerLibrary::get_instance().add_residue_library( *mod_res, new_lrots );
				}

				//finalize again just to make sure
				mod_res->nondefault(true);
				mod_res->base_restype_name( base_name );
				mod_res->finalize();

				mod_restype_set->add_residue_type( mod_res );
			}
		}

		//and last but not least we have to regenerate the rotamer library for the ligand
		if( pose.residue_type( res_pos ).is_ligand() ) {

			SingleLigandRotamerLibraryCAP old_lrots(
				static_cast< SingleLigandRotamerLibrary const * >
				( RotamerLibrary::get_instance().get_rsd_library( pose.residue_type( res_pos ))() ));

			if( old_lrots != 0 ){

				using namespace core::conformation;
				utility::vector1< ResidueOP > new_rotamers;
				new_rotamers.clear();
				utility::vector1< ResidueOP > const old_rotamers = old_lrots->get_rotamers();

				//std::cerr << "old rotamer library has " << old_rotamers.size() << "members";
				for( utility::vector1< ResidueOP>::const_iterator oldrot_it = old_rotamers.begin(); oldrot_it != old_rotamers.end(); ++oldrot_it){
					ResidueOP new_rot_res = new Residue( restype_set_->name_map(res_type_mod_name), true);
					//set the coordinates
					//1. we go over the atoms of the NEW residue on purpose, to make sure that no atom gets skipped
					for( core::Size at_ct = 1; at_ct <= new_rot_res->natoms(); at_ct++){
						if( !(*oldrot_it)->has( new_rot_res->atom_name( at_ct ) ) ){
							std::cerr << "Unexpected ERROR: when regenerating ligand rotamer library (for covalent constraints), one atom wasn't found in a template rotamer." << std::endl;
							utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
						}
						else{
							new_rot_res->set_xyz( at_ct, (*oldrot_it)->xyz( new_rot_res->atom_name( at_ct ) ) );
						}
					}
					//2. we also set the chis, to make sure everything is properly in place
					new_rot_res->chi( (*oldrot_it)->chi() );

					new_rotamers.push_back( new_rot_res );
				}
				new_lrots->set_reference_energy( old_lrots->get_reference_energy() );
				new_lrots->set_rotamers( new_rotamers );
				//std::cerr << "new rotlibrary has " << new_lrots->get_rotamers().size() << " members." << std::endl;
			}
		}
	}

	core::conformation::Residue new_res( restype_set_->name_map(res_type_mod_name), true);

	replace_residue_keeping_all_atom_positions( pose, new_res, res_pos );

	//std::cout << "APL DEBUG EnzConstraintParameters.cc::make_constraint_covalent_helper end" << std::endl;

} //make_constraint_covalent_helper



void
EnzConstraintParameters::show_definitions() const
{
	tr.Info << "parameters residue 1:" << std::endl;
	resA_->show_params();
	tr.Info << "parameters residue 2:" << std::endl;
	resB_->show_params();
}

void
EnzConstraintParameters::generate_pose_specific_data(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx) const
{
	resA_->get_pose_data(pose);
	resB_->get_pose_data(pose);
	generate_active_pose_constraints( pose, scofx );
}



/// @brief function to determine the lowest scoring constraint in a vector of input constraints
/// @brief there might be a slightly faster and slightly more complicated way to implement this,
/// @brief but since this function will only be called once per input pose (if it is called at all),
/// @brief it probably doesn't matter.
core::Size
EnzConstraintParameters::determine_best_constraint(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP scofx,
	utility::vector1< core::scoring::constraints::ConstraintCOP > candidate_csts ) const
{

	if( candidate_csts.size() < 2 ){
		std::cerr << "Error, this function should not be called with a constraint input vector containing less than 2 constraints. " << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}


	core::scoring::ScoreFunction helper_scofx;

	//need to score each constraint in the context of the pose, but it's better to keep the original pose intact
	//so let's do a full copy
	core::pose::Pose helper_pose = pose;

	//only set the constraints, so we don't lose time scoring other stuff
	helper_scofx.set_weight( core::scoring::coordinate_constraint, (*scofx)[core::scoring::coordinate_constraint] );
	helper_scofx.set_weight( core::scoring::atom_pair_constraint, (*scofx)[core::scoring::atom_pair_constraint] );
	helper_scofx.set_weight( core::scoring::angle_constraint, (*scofx)[core::scoring::angle_constraint] );
	helper_scofx.set_weight( core::scoring::dihedral_constraint, (*scofx)[core::scoring::dihedral_constraint] );



	bool first_pass = true;
	Size best_constraint(0);
	core::Real cur_low_e(0.0);

	for( Size cur_cst(1); cur_cst <= candidate_csts.size(); cur_cst++ ){

		core::scoring::constraints::ConstraintSetOP cur_cstset = new core::scoring::constraints::ConstraintSet();
		cur_cstset->add_constraint( candidate_csts[cur_cst] );
		helper_pose.constraint_set( cur_cstset );
		helper_scofx( helper_pose );

		core::Real total_cst_e = helper_pose.energies().total_energies()[core::scoring::coordinate_constraint] + helper_pose.energies().total_energies()[core::scoring::atom_pair_constraint] + helper_pose.energies().total_energies()[core::scoring::angle_constraint] + helper_pose.energies().total_energies()[core::scoring::dihedral_constraint];

		if( first_pass || ( total_cst_e < cur_low_e ) ) {
			best_constraint = cur_cst;
			cur_low_e = total_cst_e;
		}
		first_pass = false;
	}

	return best_constraint;
} //determine best constraint function


bool
EnzConstraintParameters::missing_in_pose( core::pose::Pose const & pose ) const
{
	utility::vector1< EnzCstTemplateResCacheOP > const & template_res_cache(
		protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ )->template_res_cache_ );

	for( core::Size i = 1; i <= template_res_cache.size(); ++i ){
		if( template_res_cache[i]->not_in_pose() )return true;
	}
	return false;
}

/*
utility::vector1< core::conformation::ResidueCOP >
EnzConstraintParameters::inverse_rotamers_for_residue_missing_in_pose(
	core::pose::Pose const & pose ) const
{

	EnzCstTemplateResCOP missing_template( this->get_missing_template_res( pose ) );

	for( core::Size restype_index = 1; restype_index <= missing_template->allowed_res_types().size(); ++restype_index ){
		core::chemical::ResidueTypeCOP restype( restype_set_->name( missing_template->allowed_res_types()[ restype_index ] ) );
		utility::vector1< core::conformation::ResidueCOP > raw_rotamers;
		// for now we'll only work with one rotamer
		raw_rotamers.push_back( new core::conformation::Residue( restype, true ) );


	} //loop over all allowed restypes of missing_template
}
*/

bool
EnzConstraintParameters::update_pdb_remarks(
	core::pose::Pose & pose
) const
{

	using namespace core::pose;

	core::pose::PDBInfo & pdbinfo( *(pose.pdb_info() ) );
	Remarks & rems(pose.pdb_info()->remarks() );
	EnzdesCstParamCacheOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );

	for( std::vector< core::pose::RemarkInfo >::iterator remark_it = rems.begin(); remark_it != rems.end(); remark_it++) {

		bool remark_changed(false);
		std::string chainA(""), chainB(""), resA(""), resB("");
		core::Size cst_block(0), exgeom_id( 0 );
		int rem_pdbposA(0), rem_pdbposB(0);

		if( !split_up_remark_line( remark_it->value, chainA, resA, rem_pdbposA, chainB, resB, rem_pdbposB, cst_block, exgeom_id ) ) continue;

		if( cst_block != cst_block_ ) continue;

		if( (param_cache->template_res_cache_[1]->seqpos_map_size() > 1 ) || (param_cache->template_res_cache_[2]->seqpos_map_size() > 1 ) ){
			tr << "Error in updating remarks for cst block " << cst_block_ << ". More than one seqpos in template resA or resB." << std::endl;
			return false;
		}

		core::Size seqposA( param_cache->template_res_cache_[1]->seqpos_map_begin()->first );
		core::Size seqposB( param_cache->template_res_cache_[2]->seqpos_map_begin()->first );

		if( rem_pdbposA != pdbinfo.number( seqposA ) ){
			//stupid backwards compatibility: if a ligand is present as resA and pdbpos is 0,
			//leave untouched for now
			if( (rem_pdbposA != 0 ) || !pose.residue( seqposA ).is_ligand() ){
				remark_changed = true;
				rem_pdbposA = pdbinfo.number( seqposA );
			}
		}
		if( chainA[0] != pdbinfo.chain( seqposA ) ){
			remark_changed = true;
			chainA = pdbinfo.chain( seqposA );
		}
		if( resA != pose.residue( seqposA ).name3() ){
			remark_changed = true;
			resA = pose.residue( seqposA ).name3();
		}

		if( rem_pdbposB != pdbinfo.number( seqposB) ){
			remark_changed = true;
			rem_pdbposB = pdbinfo.number( seqposB );
		}
		if( chainB[0] != pdbinfo.chain( seqposB ) ){
			remark_changed = true;
			chainB = pdbinfo.chain( seqposB );
		}
		if( resB != pose.residue( seqposB ).name3() ){
			remark_changed = true;
			resB = pose.residue( seqposB ).name3();
		}

		if( remark_changed ){
			remark_it->value = assemble_remark_line( chainA, resA, rem_pdbposA, chainB, resB, rem_pdbposB, cst_block, exgeom_id );
		}
	}
	return true;
}


EnzCstTemplateResOP
EnzConstraintParameters::nonconst_resA()
{
	return resA_;
}

EnzCstTemplateResOP
EnzConstraintParameters::nonconst_resB()
{
	return resB_;
}


EnzCstTemplateResCOP
EnzConstraintParameters::resA() const
{
	return resA_;
}


EnzCstTemplateResCOP
EnzConstraintParameters::resB() const
{
	return resB_;
}


EnzCstTemplateResCOP
EnzConstraintParameters::get_missing_template_res( core::pose::Pose const & pose ) const
{

	utility::vector1< EnzCstTemplateResCacheOP > const &  template_res_cache =
		protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ )->template_res_cache_;

	if( template_res_cache.size() != 2 ) utility_exit_with_message( "More or less than 2 template res caches detected in enzdes cst param cache");
	if( template_res_cache[1]->not_in_pose() && template_res_cache[2]->not_in_pose() ) utility_exit_with_message("Error: Both template residues are missing in the pose. This shouldn't happen...\n");

	if( template_res_cache[1]->not_in_pose() ) return resA_;
	else if( template_res_cache[2]->not_in_pose() ) return resB_;
	else{
		utility_exit_with_message("Error: no template residue is missing in the pose, this shouldn't have happened... \n");
	}
	//unreachable
	return NULL;
}

EnzCstTemplateResCOP
EnzConstraintParameters::get_missing_template_other_res( core::pose::Pose const & pose ) const
{
	utility::vector1< EnzCstTemplateResCacheOP > const &  template_res_cache =
		protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ )->template_res_cache_;

	if( template_res_cache.size() != 2 ) utility_exit_with_message( "More or less than 2 template res caches detected in enzdes cst param cache");
	if( template_res_cache[1]->not_in_pose() && template_res_cache[2]->not_in_pose() ) utility_exit_with_message("Error: Both template residues are missing in the pose. This shouldn't happen...\n");

	if( template_res_cache[1]->not_in_pose() ) return resB_;
	else if( template_res_cache[2]->not_in_pose() ) return resA_;
	else{
		utility_exit_with_message("Error: no template residue is missing in the pose, this shouldn't have happened... \n");
	}

	//unreachable
	return NULL;
}


std::set< std::string >
EnzConstraintParameters::allowed_res_name3_at_position(
	core::pose::Pose const & pose,
	core::Size seqpos ) const
{

	EnzCstTemplateResCOP template_res;
	EnzdesCstParamCacheCOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );

	std::set< std::string > to_return;

	if( param_cache->template_res_cache(1)->contains_position( seqpos ) ) template_res = resA_;
	else if( param_cache->template_res_cache(2)->contains_position( seqpos ) ) template_res = resB_;

	else return to_return;

	for( core::Size i = 1; i <= template_res->allowed_res_types().size(); ++i ){
		to_return.insert( template_res->allowed_res_types()[i] );
	}

	return to_return;
}

void
EnzConstraintParameters::set_external_position_for_resA( core::Size pos )
{
	resA_->set_external_position( pos );
}


void
EnzConstraintParameters::set_external_position_for_resB( core::Size pos )
{
	resB_->set_external_position( pos );
}


void
EnzConstraintParameters::remove_covalent_connections_from_pose( core::pose::Pose & pose ) const {

	EnzdesCstParamCacheOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );
	for( utility::vector1< CovalentConnectionReplaceInfoCOP >::iterator cov_it = param_cache->covalent_connections_.begin();
			 cov_it != param_cache->covalent_connections_.end(); ++cov_it ){
		(*cov_it)->remove_covalent_connection_from_pose( pose );
	}
	param_cache->covalent_connections_.clear();
}

void
EnzConstraintParameters::remap_resid( core::id::SequenceMapping const & smap )
{
	resA_->remap_resid( smap );
	resB_->remap_resid( smap );
}

}
} //enzdes
} //protocols
