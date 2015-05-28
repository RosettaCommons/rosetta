// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensity.cc
/// @brief  Scoring a structure against an electron density map
/// @author Frank DiMaio

// Unit Headers
#include <core/scoring/cryst/PhenixInterface.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/cryst/util.hh>

#include <core/chemical/AtomType.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <numeric/fourier/FFT.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/statistics/functions.hh>


#include <basic/options/keys/cryst.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>

// Python interpreter
#ifdef WITH_PYTHON
#include <Python.h>
#define HANDLE_PYTHON_ERROR(msg) \
  if (PyErr_Occurred()) { \
    PyErr_PrintEx(0); \
    utility_exit_with_message(msg); \
  }
#endif

// C++ headers
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>

#ifndef WIN32
#include <unistd.h>
#endif

namespace core {
namespace scoring {
namespace cryst {

static thread_local basic::Tracer TR( "core.scoring.cryst.PhenixInterface" );

////////////////////////////////

PhenixInterface &getPhenixInterface() {
	static PhenixInterface PHENIX;
	return PHENIX;
}


/// @brief constructor
PhenixInterface::PhenixInterface() {
#ifdef WITH_PYTHON
	char filename[16] = "_ros_tempXXXXXX";
	tempdir_ = std::string( mkdtemp(filename) );

	// default options
	mtzfile_ = basic::options::option[ basic::options::OptionKeys::cryst::mtzfile ]();
	adp_strategy_ = "individual";
	target_function_ = "ml";
	twin_law_ = "";
	map_type_ = "";
	res_low_ = res_high_ = sharpen_b_ = 0;  // data from mtz file

	// set up python environment
	// use the phenix version of python...
	std::string PHENIXHOME ( getenv("ROSETTA_PHENIX_DIST") );
	std::string PHENIX_BIN ( getenv("ROSETTA_PHENIX_BIN") );
	std::string PHENIX_PYTHON_PATHS ( getenv("ROSETTA_PHENIX_MODULES") );


	std::vector< char > PHENIX_BIN_char(PHENIX_BIN.c_str(), PHENIX_BIN.c_str()+PHENIX_BIN.length()+1);
	Py_SetProgramName( &PHENIX_BIN_char[0] );

	Py_Initialize();  // initialize the python interpreter

	// set default phenix paths
	char path[] = "path";
	PyObject *sys_path = PySys_GetObject(path);
  std::stringstream path_strings(PHENIX_PYTHON_PATHS);
	std::string item;
  while (std::getline(path_strings, item, ':')) {
		std::vector< char > item_char(item.c_str(), item.c_str()+item.length()+1);
		PyObject *search_path = PyUnicode_FromString(&item_char[0]);
		PyList_Append(sys_path, search_path);
  }
  HANDLE_PYTHON_ERROR("failed setting up Phenix interface");

	phenix_home_ = PHENIXHOME;

	target_evaluator_ = NULL;
#endif
}


/// @brief score a structure
core::Real PhenixInterface::getScore (
#ifdef WITH_PYTHON
		core::pose::Pose const & pose)
#else
		core::pose::Pose const & /*pose*/)
#endif
{
#ifdef WITH_PYTHON
	if (!target_evaluator_) {
		initialize_target_evaluator( pose );
	}

	PyObject *pCoords = pose_to_pycoords( pose );
	PyObject *pMethod = PyString_FromString("update_sites_1d");
	PyObject *pValue = PyObject_CallMethodObjArgs(target_evaluator_, pMethod, pCoords, NULL);
	HANDLE_PYTHON_ERROR("PhenixInterface::getScore() : error updating sites");

	Py_XDECREF(pCoords);
	Py_XDECREF(pMethod);
	Py_XDECREF(pValue);

	pMethod = PyString_FromString("compute_target");
	PyObject *result_tuple = PyObject_CallMethodObjArgs(target_evaluator_, pMethod, Py_False, Py_False, NULL);
  HANDLE_PYTHON_ERROR("error computing X-ray target");

	Py_XDECREF(pMethod);

	// get score, r_work, and r_free
	core::Real score = PyFloat_AsDouble( PyTuple_GetItem(result_tuple, 0) );
	Py_XDECREF(result_tuple);

	//core::Real rwork = getR();
	//core::Real rfree = getRfree();

	return score;
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return 0;
}


/// @brief score a structure with derivatives
core::Real PhenixInterface::getScoreAndDerivs (
#ifdef WITH_PYTHON
		core::pose::Pose const & pose,
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > & grads)
#else
		core::pose::Pose const & /*pose*/,
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > & /*grads*/)
#endif
{
#ifdef WITH_PYTHON
	if (!target_evaluator_) {
		initialize_target_evaluator( pose );
	}

	PyObject *pCoords = pose_to_pycoords( pose );
	PyObject *pMethod = PyString_FromString("update_sites_1d");
	PyObject *pValue = PyObject_CallMethodObjArgs(target_evaluator_, pMethod, pCoords, NULL);
	HANDLE_PYTHON_ERROR("PhenixInterface::getScoreAndDerivs() : error updating sites");
	Py_XDECREF(pCoords);
	Py_XDECREF(pMethod);
	Py_XDECREF(pValue);

	pMethod = PyString_FromString("compute_functional_and_gradients_rosetta");
	PyObject *result_tuple = PyObject_CallMethodObjArgs(target_evaluator_, pMethod, NULL);
  HANDLE_PYTHON_ERROR("error computing X-ray target and gradients");

	// get score, r_work, and r_free
	core::Real score = PyFloat_AsDouble( PyTuple_GetItem(result_tuple, 0) );  // borrowed object no need to deref
	PyObject *grad_list = PyTuple_GetItem(result_tuple, 1); // borrowed object no need to deref
 debug_assert (grad_list != NULL);
	HANDLE_PYTHON_ERROR("error getting gradients list");

	// parse the python result
	pylist_to_grads(pose, grad_list, grads);
	HANDLE_PYTHON_ERROR("error converting gradients list");

	// finally we can free the results
	Py_XDECREF(result_tuple);
	Py_XDECREF(pMethod);

	//core::Real rwork = getR();
	//core::Real rfree = getRfree();

	return score;
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return 0;
}


/// @brief fit bfactors
void PhenixInterface::fitBfactors (
#ifdef WITH_PYTHON
		core::pose::Pose & pose)
#else
		core::pose::Pose & /*pose*/)
#endif
{
#ifdef WITH_PYTHON
	// this may not be the best place to put this...
	// if disulfides were formed or broken they will be updated in the asymm pose
	//   - this may cause the python-side conformation and the rosetta-side conformation to go out of sync
	//   - to fix this we will update the symmetric pose's disulfides here, when we invalidate the python-side conformation
	// Ideally this would be done when _initializing_ the python-side conformation
	//   - but, the python-side conf. may be initialized during scoring so we only have const access to the pose
	//pose.conformation().detect_disulfides();
	core::pose::initialize_disulfide_bonds( pose );

	chdir(tempdir_.c_str());

	char pdbout[256] = "outtmp_adp.pdb";
	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	fix_bfactorsH( pose_asu );
	fix_bfactorsMissing( pose_asu );

	pose_asu.dump_pdb( pdbout , "1" );

	// load module
	PyObject *pModule = PyImport_ImportModule( "phenix.refinement.command_line" );
	HANDLE_PYTHON_ERROR("importing phenix.refine failed");

	// prepare arguments
	// 'phenix.refine'
	// ['rosetta_mr.pdb' '2W9Q_4A.mtz' 'strategy=group_adp' 'group_adp_refinement_mode=two_adp_groups_per_residue' 'main.number_of_macro_cycles=1' '--overwrite' ...
	// 'refinement.input.xray_data.high_resolution=4.0' 'refinement.input.xray_data.low_resolution=4.0' ...
	//  'write_eff_file=false' 'write_geo_file=false' 'write_def_file=false' 'write_maps=false' 'write_map_coefficients=false']
	std::string mtzout = mtzfile_;
	if (mtzout[0] != '/') mtzout = "../"+mtzout;   // if path is relative, update
	std::string adp_strategy; // = "strategy=group_adp";
	std::string adp_groupstrat; // = "group_adp_refinement_mode=two_adp_groups_per_residue";

 	if (adp_strategy_ == "group" || adp_strategy_ == "group1") {
		TR << "Refining group ADPs (1 group per residue)" << std::endl;
 		adp_strategy = "strategy=group_adp";
		adp_groupstrat = "group_adp_refinement_mode=one_adp_group_per_residue";
	} else if (adp_strategy_ == "group2") {
		TR << "Refining group ADPs (2 groups per residue)" << std::endl;
 		adp_strategy = "strategy=group_adp";
		adp_groupstrat = "group_adp_refinement_mode=two_adp_groups_per_residue";
 	} else if (adp_strategy_ == "individual") {
		TR << "Refining individual ADPs" << std::endl;
 		adp_strategy = "strategy=individual_adp";
		adp_groupstrat = "group_adp_refinement_mode=one_adp_group_per_residue";  // ignored
 	} else {
		TR.Error << "Unknown ADP strategy " << adp_strategy_ << std::endl;
		utility_exit();
	}

	std::string hires_str="refinement.input.xray_data.high_resolution=None",
	            lowres_str="refinement.input.xray_data.low_resolution=None";

	if (res_high_ > 0) {
		std::ostringstream oss;
		oss << "refinement.input.xray_data.high_resolution=" << res_high_;
		hires_str=oss.str();
	}
	if (res_low_ > 0) {
		std::ostringstream oss;
		oss << "refinement.input.xray_data.low_resolution=" << res_low_;
		lowres_str=oss.str();
	}

	bool has_twin_law=false;
	std::string twinlawstr;
	if (( twin_law_.length() != 0 ) && (twin_law_ != "None")) {
		has_twin_law=true;
		std::ostringstream oss;
		oss << "twin_law=" << "\"" << twin_law_ << "\"";
		twinlawstr = oss.str();
	}

	// we need to grab stdout from the python script so set up a redirect here
	{
		// buid the input list
		Size nargs_no_cif = 13 + (has_twin_law?1:0);
		Size nargs = cif_files_.size() + nargs_no_cif;

		PyObject* inlist = PyList_New(nargs);
		PyList_SetItem(inlist, 0, PyString_FromString(pdbout));
		std::vector< char > mtzout_char(mtzout.c_str(), mtzout.c_str()+mtzout.length()+1);
		PyList_SetItem(inlist, 1, PyString_FromString(&mtzout_char[0]));
		std::vector< char > adp_strategy_char(adp_strategy.c_str(), adp_strategy.c_str()+adp_strategy.length()+1);
		PyList_SetItem(inlist, 2, PyString_FromString(&adp_strategy_char[0]));
		std::vector< char > adp_groupstrat_char(adp_groupstrat.c_str(), adp_groupstrat.c_str()+adp_groupstrat.length()+1);
		PyList_SetItem(inlist, 3, PyString_FromString(&adp_groupstrat_char[0]));
		std::vector< char > hires_str_char(hires_str.c_str(), hires_str.c_str()+hires_str.length()+1);
		PyList_SetItem(inlist, 4, PyString_FromString(&hires_str_char[0]));
		std::vector< char > lowres_str_char(lowres_str.c_str(), lowres_str.c_str()+lowres_str.length()+1);
		PyList_SetItem(inlist, 5, PyString_FromString(&lowres_str_char[0]));
		char arg6_str[] = "main.number_of_macro_cycles=1";
		PyList_SetItem(inlist, 6, PyString_FromString(arg6_str));
		char arg7_str[] = "write_eff_file=false";
		PyList_SetItem(inlist, 7, PyString_FromString(arg7_str));
		char arg8_str[] = "write_geo_file=false";
		PyList_SetItem(inlist, 8, PyString_FromString(arg8_str));
		char arg9_str[] = "write_def_file=false";
		PyList_SetItem(inlist, 9, PyString_FromString(arg9_str));
		char arg10_str[] = "write_maps=false";
		PyList_SetItem(inlist, 10, PyString_FromString(arg10_str));
		char arg11_str[] = "write_map_coefficients=false";
		PyList_SetItem(inlist, 11, PyString_FromString(arg11_str));
		char arg12_str[] = "--overwrite";
		PyList_SetItem(inlist, 12, PyString_FromString(arg12_str));
		if (has_twin_law) {
			std::vector< char > twinlawstr_char(twinlawstr.c_str(), twinlawstr.c_str()+twinlawstr.length()+1);
			PyList_SetItem(inlist, 13, PyString_FromString(&twinlawstr_char[0]));
		}

		for (int i=1; i<= (int)cif_files_.size(); ++i) {
			std::string cif_file_i = cif_files_[i];
			if (cif_file_i[0] != '/') cif_file_i = "../"+cif_file_i;   // if path is relative, update
			std::vector< char > cif_file_i_char(cif_file_i.c_str(), cif_file_i.c_str()+cif_file_i.length()+1);
			PyList_SetItem(inlist, nargs_no_cif+i-1, PyString_FromString(&cif_file_i_char[0]));
		}

		PyObject *cmd = PyString_FromString("run");
		PyObject *exe = PyString_FromString("phenix.refine");
		PyObject *pValue = PyObject_CallMethodObjArgs( pModule, cmd, exe, inlist, NULL);

		//PyObject *pValue = PyObject_CallMethod(pModule, "run", "(s[sssssssssssss])",
		//	"phenix.refine", pdbout, mtzout.c_str(), adp_strategy.c_str(), adp_groupstrat.c_str(), hires_str.c_str(), lowres_str.c_str(),
		//	"main.number_of_macro_cycles=1", "write_eff_file=false", "write_geo_file=false",
		//	"write_def_file=false", "write_maps=false", "write_map_coefficients=false", "--overwrite");
		HANDLE_PYTHON_ERROR("phenix.refine failed");

		if (!pValue) {
			utility_exit_with_message( "In PhenixInterface::fitBfactors: error!" );
		}

		// cleanup
		Py_XDECREF(inlist);
		Py_XDECREF(pValue);
		Py_XDECREF(cmd);
		Py_XDECREF(exe);
	}

	// cleanup
	chdir("..");

	// read in b factors & update
	// outfile is outtmp_adp_refine_001.pdb
	stealBfactorsFromFile( pose , tempdir_+"/outtmp_adp_refine_001.pdb" );

	// after B factor refinement the python fmodel is out of date
	// nuke the target evaluator
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}

	// get rid of ref pose
	if (ref_pose_) {
		ref_pose_ = NULL;
	}
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}


#ifdef WITH_PYTHON
PyObject* PhenixInterface::pose_to_pycoords( core::pose::Pose const & pose ) {
	using namespace core::pose::symmetry;
	using namespace core::conformation::symmetry;

	// copy coordinates from pose
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	if (is_symmetric( pose )) {
		SymmetricConformation const & symm_conf (
		      dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		for (int i=1; i<=(int)symm_info->num_total_residues_without_pseudo(); ++i) {
			if (pose.residue(i).aa() == core::chemical::aa_vrt) continue;
			if (!symm_info->bb_is_independent(i)) continue;

			core::conformation::Residue const &rsd_i = pose.residue(i);

			// disulfide changes
			bool ignore_cys_hg = ( !pose.residue(i).has_variant_type( core::chemical::DISULFIDE ) &&
			                       ref_pose_->residue(i).has_variant_type( core::chemical::DISULFIDE ) );
			bool add_cys_hg = ( pose.residue(i).has_variant_type( core::chemical::DISULFIDE ) &&
			                    !ref_pose_->residue(i).has_variant_type( core::chemical::DISULFIDE ) );

			for (int j=1; j<=(int)rsd_i.natoms(); ++j) {
				if (pose.residue_type(i).atom_name(j) == " CEN") continue;
				if (pose.residue_type(i).atom_type(j).is_virtual()) continue;

				if (ignore_cys_hg && pose.residue(i).atom_name(j) == " HG ") {
					TR.Debug << "Warning: Correcting disulfide at residue " << i << std::endl;
					continue;
				}

				coords.push_back( rsd_i.xyz( j ) );
			}

			if (add_cys_hg) { // disulfide accross symm interface
				TR.Debug << "Warning: Correcting disulfide at residue " << i << std::endl;

				// convert CYS->CYD
				core::conformation::Residue newCys = rsd_i;
				chemical::ResidueTypeSet const & residue_type_set = newCys.type().residue_type_set();
				chemical::ResidueTypeCOPs const & possible_types = residue_type_set.name3_map( "CYS" );
				utility::vector1< std::string > variant_types = newCys.type().properties().get_list_of_variants();
				variant_types.erase( std::find( variant_types.begin(), variant_types.end(), "DISULFIDE" ) );

				// Run through all possible new residue types.
				chemical::ResidueTypeCOP matchedRes;
				bool matched;
				for ( chemical::ResidueTypeCOPs::const_iterator
					type_iter = possible_types.begin(), type_end = possible_types.end();
					type_iter != type_end; ++type_iter ) {

					matched = true;
					for ( Size kk = 1; kk <= variant_types.size(); ++kk ) {
						if ( ! (*type_iter)->has_variant_type( variant_types[ kk ] ) ) {
							matched = false;
							break;
						}
					}

					if ( matched ) { // Do replacement.
						core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( **type_iter, newCys, pose.conformation()  );
						copy_residue_coordinates_and_rebuild_missing_atoms( newCys, *new_res, pose.conformation() );
						coords.push_back( new_res->xyz( new_res->atom_index(" HG ")) );
						break;
					}
				}
				if (!matched) {
					// fallback: just copy ref_pose_
					coords.push_back( ref_pose_->residue(i).xyz( ref_pose_->residue(i).atom_index(" HG ")) );
				}
			}
		}
	} else {
		for (int i=1; i<=(int)pose.total_residue(); ++i) {
			if (pose.residue(i).aa() == core::chemical::aa_vrt) continue;
			core::conformation::Residue const &rsd_i = pose.residue(i);
			for (int j=1; j<=(int)rsd_i.natoms(); ++j) {
				if (pose.residue_type(i).atom_name(j) == " CEN") continue;
				if (pose.residue_type(i).atom_type(j).is_virtual()) continue;
				coords.push_back( rsd_i.xyz( j ) );
			}
		}
	}

	//fpd everything is cast to a float before passing to python
	//fpd refolding the protein sometimes leads to small floating point changes
	//fpd which has sometimes not-so-small effects on the score (it is not clear to me why)
	//fpd casting to float should take care of this
	// convert to python obj
	PyObject *pCoords = PyList_New(coords.size()*3);
	for (int i=1; i<=(int)coords.size(); ++i) {
		//PyList_SET_ITEM( pCoords, 3*(i-1)+0, PyFloat_FromDouble( (float)coords[i][0]) );
		//PyList_SET_ITEM( pCoords, 3*(i-1)+1, PyFloat_FromDouble( (float)coords[i][1]) );
		//PyList_SET_ITEM( pCoords, 3*(i-1)+2, PyFloat_FromDouble( (float)coords[i][2]) );
		PyList_SetItem( pCoords, 3*(i-1)+0, PyFloat_FromDouble( (float)coords[i][0]) );
		PyList_SetItem( pCoords, 3*(i-1)+1, PyFloat_FromDouble( (float)coords[i][1]) );
		PyList_SetItem( pCoords, 3*(i-1)+2, PyFloat_FromDouble( (float)coords[i][2]) );
	}

	return pCoords;  // calling function must dereference this!
	return NULL;
}
#endif

#ifdef WITH_PYTHON
void PhenixInterface::pylist_to_grads(
			core::pose::Pose const & pose,
			PyObject* pygrads,
			utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > & grads ) {
	using namespace core::pose::symmetry;
	using namespace core::conformation::symmetry;

debug_assert (pygrads != NULL);

	SymmetryInfoCOP symm_info;
	int nres = (int)pose.total_residue();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}

	grads.resize( nres );
	int listCounter = 0;
	for (int i=1; i<=(int)pose.total_residue(); ++i) {
		if (pose.residue(i).aa() == core::chemical::aa_vrt) continue;
		if (symm_info && !symm_info->bb_is_independent(i)) continue;

		core::conformation::Residue const &rsd_i = pose.residue(i);

		// disulfide changes
		bool ignore_cys_hg = ( !pose.residue(i).has_variant_type( core::chemical::DISULFIDE ) &&
		                       ref_pose_->residue(i).has_variant_type( core::chemical::DISULFIDE ) );
		bool add_cys_hg = ( pose.residue(i).has_variant_type( core::chemical::DISULFIDE ) &&
		                    !ref_pose_->residue(i).has_variant_type( core::chemical::DISULFIDE ) );

		grads[i].resize( rsd_i.natoms() );
		for (int j=1; j<=(int)rsd_i.natoms(); ++j) {
			if (pose.residue_type(i).atom_type(j).is_virtual()) {
				grads[i][j][0] = grads[i][j][1] = grads[i][j][2] = 0; // placeholder
			} else if (pose.residue_type(i).atom_name(j) == " CEN") {
				grads[i][j][0] = grads[i][j][1] = grads[i][j][2] = 0; // placeholder
			} else if (ignore_cys_hg && pose.residue(i).atom_name(j) == " HG ") {
				grads[i][j][0] = grads[i][j][1] = grads[i][j][2] = 0; // placeholder
			} else {
				PyObject *gx = PyList_GetItem( pygrads, listCounter++);
				PyObject *gy = PyList_GetItem( pygrads, listCounter++);
				PyObject *gz = PyList_GetItem( pygrads, listCounter++);
				if (PyErr_Occurred()) {
					TR << "PyList_GetItem() failure at " << i << "," << j << "(" << listCounter << ")" << std::endl;
					utility_exit_with_message("PhenixInterface::pylist_to_grads failed");
				}
			debug_assert(! ((gx == NULL) || (gy == NULL) || (gz == NULL)));
				grads[i][j][0] = PyFloat_AsDouble( gx );
				grads[i][j][1] = PyFloat_AsDouble( gy );
				grads[i][j][2] = PyFloat_AsDouble( gz );
			}
		}

		// ignore CYS H
		if (add_cys_hg) {
			listCounter += 3; // placeholder
		}
	}
}
#endif


void PhenixInterface::stealBfactorsFromFile(
#ifdef WITH_PYTHON
		core::pose::Pose & pose,
		std::string filename)
#else
		core::pose::Pose & /*pose*/,
		std::string /*filename*/)
#endif
{
#ifdef WITH_PYTHON
	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	// open PDB file
	std::ifstream PDBIN( filename.c_str() , std::ifstream::in );
	std::string atmLine;
	while(!PDBIN.eof()) {
		std::getline(PDBIN,atmLine);
		if (atmLine.substr(0,6) != "ATOM  " && atmLine.substr(0,6) != "HETATM") continue;
	  if (atmLine.size() < 78) {
			TR  << "format error:" << std::endl;
			TR  << atmLine << std::endl;
			continue;
		}

		core::Size resid = (core::Size) atoi( atmLine.substr(22,4).c_str() );
		std::string chainid = atmLine.substr(21,1);
		std::string icode = atmLine.substr(26,1);
		core::Size resid_ros = pose_asu.pdb_info()->pdb2pose( chainid[0], resid, icode[0] );
		core::Size residALT_ros = pose.pdb_info()->pdb2pose( chainid[0], resid, icode[0] );

		if (resid_ros != residALT_ros) {
			TR << "FATAL ERROR!" << resid_ros << " != " << residALT_ros << std::endl;
			TR << "reading line: " << atmLine << std::endl;
			utility_exit_with_message( "aborting!" );
		}

		if (resid_ros == 0) {
			TR << "error finding res! " << resid << chainid << " first res = " << pose.pdb_info()->pose2pdb(1) << std::endl;
			utility_exit_with_message( "aborting!" );
		}

		// fpd handle H internally
		if (atmLine.substr(76,2) == " H") continue;

		core::Size atmid = pose_asu.residue_type(resid_ros).atom_index( atmLine.substr(12,4) );
		core::Real B = atof( atmLine.substr(60,6).c_str() );

		pose.pdb_info()->temperature( resid_ros, atmid, B );
	}

	// fix H
	fix_bfactorsH( pose );

#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

/// @brief update fcalc
void PhenixInterface::updateFcalc () {
	// old
	return;
}

/// @brief update mask
void PhenixInterface::updateSolventMask () {
#ifdef WITH_PYTHON
	if (!target_evaluator_) return;

	char update_fmask_str[] = "update_fmask";
	PyObject_CallMethod(target_evaluator_, update_fmask_str, NULL);
  HANDLE_PYTHON_ERROR("update_fmask failed");

	core::Real rwork = getR();
	core::Real rfree = getRfree();
	TR << "After optimizeSolventMask r/rfree = " << rwork << "/" << rfree << std::endl;

#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

/// @brief update mask
void PhenixInterface::updateSolventMask (

#ifdef WITH_PYTHON
		core::pose::Pose const & pose)
{
		if (!target_evaluator_) {
		initialize_target_evaluator( pose );
	}

	PyObject *pCoords = pose_to_pycoords( pose );
	PyObject *pMethod = PyString_FromString("update_sites_1d");
	PyObject *pValue = PyObject_CallMethodObjArgs(target_evaluator_, pMethod, pCoords, NULL);
	HANDLE_PYTHON_ERROR("PhenixInterface::getScore() : error updating sites");

	core::Real rwork = getR();
	core::Real rfree = getRfree();
	TR << "After optimizeSolventMask r/rfree = " << rwork << "/" << rfree << std::endl;

	// free memory
	Py_XDECREF(pCoords);
	Py_XDECREF(pMethod);
	Py_XDECREF(pValue);

#else
		core::pose::Pose const &)
{
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

/// @brief explicitly recompute ksol/bsol
void PhenixInterface::optimizeSolventMask () {
#ifdef WITH_PYTHON
	if (!target_evaluator_) return;

	char optimize_mask_str[] = "optimize_mask";
	PyObject *pValue = PyObject_CallMethod(target_evaluator_, optimize_mask_str, NULL);
  HANDLE_PYTHON_ERROR("optimize_mask failed");

	core::Real rwork = getR();
	core::Real rfree = getRfree();
	TR << "After optimizeSolventMask r/rfree = " << rwork << "/" << rfree << std::endl;

	// free memory
	Py_XDECREF(pValue);

#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}


/// @brief explicitly recompute ksol/bsol
void PhenixInterface::optimizeSolvParams () {
#ifdef WITH_PYTHON
	if (!target_evaluator_) return;

	char update_solvent_and_scale_str[] = "update_solvent_and_scale";
	PyObject *pValue = PyObject_CallMethod(target_evaluator_, update_solvent_and_scale_str, NULL);
  HANDLE_PYTHON_ERROR("update_solvent_and_scale failed");
	core::Real rwork = getR();
	core::Real rfree = getRfree();
	TR << "After optimizeSolvParams r/rfree = " << rwork << "/" << rfree << std::endl;

	// free memory
	Py_XDECREF(pValue);

#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

/// @brief explicitly recompute ksol/bsol and fmask
void PhenixInterface::optimizeSolvParamsAndMask () {
#ifdef WITH_PYTHON
	if (!target_evaluator_) return;

	char optimize_mask_and_update_solvent_and_scale_str[] = "optimize_mask_and_update_solvent_and_scale";
	PyObject *pValue = PyObject_CallMethod(target_evaluator_, optimize_mask_and_update_solvent_and_scale_str, NULL);
	HANDLE_PYTHON_ERROR("optimizeSolvParamsAndMask failed");

	core::Real rwork = getR();
	core::Real rfree = getRfree();
	TR << "After optimize_mask_and_update_solvent_and_scale r/rfree = " << rwork << "/" << rfree << std::endl;

	// free memory
	Py_XDECREF(pValue);
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}


/// @brief update res limits
void PhenixInterface::setResLimits(
#ifdef WITH_PYTHON
		core::Real res_high /*=0.0*/,
		core::Real res_low /*=0.0*/)
#else
		core::Real /*res_high=0.0*/,
		core::Real /*res_low=0.0*/)
#endif
{
#ifdef WITH_PYTHON
	// set member vars
	res_low_ = res_low;
	res_high_ = res_high;

	// nuke the target evaluator
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

void
PhenixInterface::setTwinLaw(
#ifdef WITH_PYTHON
		std::string twin_law)
#else
		std::string /*twin_law*/)
#endif
{
#ifdef WITH_PYTHON
	twin_law_ = twin_law;

	// nuke the target evaluator
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

void
PhenixInterface::setAlgorithm(
#ifdef WITH_PYTHON
		std::string algo)
#else
		std::string /*algo*/)
#endif
{
#ifdef WITH_PYTHON
	algo_ = algo;

	// nuke the target evaluator
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

void
PhenixInterface::set_map_type(
#ifdef WITH_PYTHON
		std::string map_type)
#else
		std::string /*map_type*/)
#endif
{
#ifdef WITH_PYTHON
	map_type_ = map_type;

	// nuke the target evaluator
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}


std::string PhenixInterface::getInfoLine() {
#ifdef WITH_PYTHON
	if (!target_evaluator_) return "0.00/0.00";

	core::Real rwork = getR();
	core::Real rfree = getRfree();

	std::ostringstream oss;
	oss << rwork << "/" << rfree;
	return oss.str();
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return "";
}

core::Real PhenixInterface::getR() {
#ifdef WITH_PYTHON
	char r_work_str[] = "r_work";
	PyObject *r_work_py = PyObject_CallMethod(target_evaluator_, r_work_str, NULL);
	core::Real rwork = 0;
	if (r_work_py != NULL) {
		rwork = PyFloat_AsDouble(r_work_py);
		Py_XDECREF(r_work_py);
	}
	return rwork;
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return 0;
}


core::Real PhenixInterface::getRfree() {
#ifdef WITH_PYTHON
	char r_free_str[] = "r_free";
	PyObject *r_free_py = PyObject_CallMethod(target_evaluator_, r_free_str, NULL);
	core::Real rfree = 0;
	if (r_free_py != NULL) {
		rfree = PyFloat_AsDouble(r_free_py);
		Py_XDECREF(r_free_py);
	}
	return rfree;
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return 0;
}


void PhenixInterface::initialize_target_evaluator(
#ifdef WITH_PYTHON
		core::pose::Pose const & pose,
		std::string eff_file /*=""*/)
#else
		core::pose::Pose const & /*pose*/,
		std::string /*eff_file=""*/)
#endif
{
#ifdef WITH_PYTHON
	// nuke the target evaluator if one exists
	if (target_evaluator_) {
		Py_XDECREF(target_evaluator_);
		target_evaluator_ = NULL;
	}

	// obvious error check 1
	//   check that the crystal refinement flag is specified
	//   if not, don't die but warn the user
	if (!basic::options::option[ basic::options::OptionKeys::cryst::crystal_refine ].user())
		TR << "[WARNING] The flag -cryst::cryst_refine should be given for crystal refinement!  Continuing..." << std::endl;

	// obvious error check 3
	//    make sure mtzfile is specified
	if (mtzfile_.length() == 0)
		utility_exit_with_message( "No MTZ file loaded!  An MTZ file must be specified with -cryst::mtzfile" );

	// obvious error check 2
	//   make sure CRYST1 line is in input structure
	//   if not die (since next step will certainly fail)
	if ( !pose.pdb_info() || pose.pdb_info()->crystinfo().A()*pose.pdb_info()->crystinfo().B()*pose.pdb_info()->crystinfo().C() == 0 )
		utility_exit_with_message( "Invalid crystal parameters!  Does the input PDB contain a valid CRYST1 line?" );

	// first time we call this we read in PDB from disk
	char buffer[256] = "outtmp.pdb";
	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	core::pose::initialize_disulfide_bonds( pose_asu ); //fpd

	// sometimes H bfactors get lost/don't exist ... fix them
	//fix_bfactorsH( pose_asu );
	//fix_bfactorsMissing( pose_asu );

	ref_pose_ = core::pose::PoseOP( new core::pose::Pose( pose_asu ) );
	//ref_pose_->conformation().detect_disulfides();

	chdir(tempdir_.c_str());

	pose_asu.dump_pdb( buffer , "1" );

	// load the module
	PyObject *pModule = PyImport_ImportModule( "phenix.rosetta.xray_target" );
	HANDLE_PYTHON_ERROR("In PhenixInterface::initialize_target_evaluator: error importing phenix.rosetta.xray_target");

	// prepare arguments
	// ['outtmp.pdb' 'file.mtz' 'high_resolution=4.0' 'low_resolution=99.0']

	// (1) pdb file
	std::string arg1(buffer);
	std::vector< char > arg1_char(arg1.c_str(), arg1.c_str()+arg1.length()+1);

	// (2) mtz file
	std::string arg2 = mtzfile_;
	if (arg2[0] != '/') arg2 = "../"+arg2;   // if path is relative, update
	std::vector< char > arg2_char(arg2.c_str(), arg2.c_str()+arg2.length()+1);

	// (3-4) res limits
	std::ostringstream oss;
	if (res_high_ > 0) {
		oss << "input.xray_data.high_resolution=" << res_high_;
	} else {
		oss << "input.xray_data.high_resolution=None";
	}
	std::string arg3 = oss.str();
	std::vector< char > arg3_char(arg3.c_str(), arg3.c_str()+arg3.length()+1);

	oss.clear(); oss.str("");
	if (res_low_ > 0) {
		oss << "input.xray_data.low_resolution=" << res_low_;
	} else {
		oss << "input.xray_data.low_resolution = None";
	}
	std::string arg4 = oss.str();
	std::vector< char > arg4_char(arg4.c_str(), arg4.c_str()+arg4.length()+1);

	// (5) target fn
	oss.clear(); oss.str("");
	oss << "options.target_name=" << target_function_;
	std::string arg5 = oss.str();
	std::vector< char > arg5_char(arg5.c_str(), arg5.c_str()+arg5.length()+1);

	// (6) if target function==lsq then use lsq bulk solvent
	oss.clear(); oss.str("");
	oss << "bss.target = ";
	if (( target_function_ == "lsq" ) ||
			(( twin_law_.length() != 0 ) && (twin_law_ != "None"))) {
		oss << "ls_wunit_k1";
	} else {
		oss << "ml";
	}
	std::string arg6 = oss.str();
	std::vector< char > arg6_char(arg6.c_str(), arg6.c_str()+arg6.length()+1);

	// (7) sf calc algo
	oss.clear(); oss.str("");
	oss << "structure_factors_accuracy.algorithm = ";
	if ( algo_.length() == 0 )
		oss << "fft";
	else
		oss << algo_;
	std::string arg7 = oss.str();
	std::vector< char > arg7_char(arg7.c_str(), arg7.c_str()+arg7.length()+1);

	// (8) twin law
	oss.clear(); oss.str("");
	oss << "options.twin_law = ";
	if ( twin_law_.length() == 0 )
		oss << "None";
	else
		oss << twin_law_;
	std::string arg8 = oss.str();
	std::vector< char > arg8_char(arg8.c_str(), arg8.c_str()+arg8.length()+1);

	// (9) map type
	oss.clear(); oss.str("");
  if (( map_type_.length() != 0 ) && (map_type_ != "Auto") && (map_type_.substr(0,2) != "%%") ) { // 3rd case bad RS tag
		if (map_type_ == "dm" || map_type_ == "density_modify") {
			oss << "density_modify=True";
		} else if (map_type_ == "prime_and_switch") {
			oss << "prime_and_switch=True";
		} else {
			oss << "map_type=" << map_type_;
		}
	} else {
		// default: prime and switch ???
		oss << "prime_and_switch=True";
	}
	std::string arg9 = oss.str();
	std::vector< char > arg9_char(arg9.c_str(), arg9.c_str()+arg9.length()+1);

	// (10) b sharpen
	oss.clear(); oss.str("");
  if (sharpen_b_ != 0.0) {
			oss << "sharpen_b=" << sharpen_b_;
	}
	std::string arg10 = oss.str();
	std::vector< char > arg10_char(arg10.c_str(), arg10.c_str()+arg10.length()+1);

	// (optional) .eff file
	std::vector< char > arg_opt;
	if (eff_file.length() > 0) {
		arg_opt = std::vector< char >(eff_file.c_str(), eff_file.c_str()+eff_file.length()+1);
	}


	{ // grab stdout from the python script
		const int MAXLEN=16000;
		char buffer[MAXLEN+1] = {0};
		//int out_pipe[2], saved_stdout;

		char run_str[] = "run";
		char s10_str[] = "([ssssssssss])";
		char s11_str[] = "([sssssssssss])";

 		if (eff_file.length() == 0)
 			target_evaluator_ = PyObject_CallMethod(pModule,
				run_str,
				s10_str,
				&arg1_char[0],  // pdb file
				&arg2_char[0],  // mtz  file
				&arg3_char[0],  // high resolution
				&arg4_char[0],  // low resolution
				&arg5_char[0],  // target function
				&arg6_char[0],  // bss target
				&arg7_char[0],  // sf calc method
				&arg8_char[0],  // twin law
				&arg9_char[0],  // map type
				&arg10_char[0]); // b sharpen
 		else
 			target_evaluator_ = PyObject_CallMethod(pModule,
				run_str,
				s11_str,
 			  &arg1_char[0],  // pdb file
				&arg2_char[0],  // mtz file
				&arg_opt[0],    // eff file
        &arg3_char[0],  // high resolution
        &arg4_char[0],  // low resolution
        &arg5_char[0],  // target function
        &arg6_char[0],  // bss target
        &arg7_char[0],  // sf calc method
        &arg8_char[0],  // twin law
        &arg9_char[0],  // map type
				&arg10_char[0]); // b sharpen

		HANDLE_PYTHON_ERROR("initialization of target evaluator failed");
		//
		if (!target_evaluator_) {
			TR << buffer << std::endl;
			utility_exit_with_message( "In PhenixInterface::initialize_target_evaluator: error initializing evaluator!" );
		}
	}

	//core::Real rwork = getR();
	//core::Real rfree = getRfree();


	Py_XDECREF(pModule);

	chdir("..");
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
}

std::string PhenixInterface::calculateDensityMap (
#ifdef WITH_PYTHON
		core::pose::Pose & pose,
		bool /*no_sidechain =false*/)
#else
		core::pose::Pose & /*pose*/,
		bool /*no_sidechain=false*/)
#endif
{
#ifdef WITH_PYTHON
	if (!target_evaluator_) {
		initialize_target_evaluator( pose );
	}

	chdir(tempdir_.c_str());

	char write_map_str[] = "write_map";
	char s1_str[] = "(s)";
	char outtmp_str[] = "outtmp.map";
	PyObject_CallMethod(target_evaluator_, write_map_str, s1_str, outtmp_str);

	HANDLE_PYTHON_ERROR("PhenixInterface::calculateDensityMap() : error writing map");

	chdir("..");
	return tempdir_+"/outtmp.map";
#else
	utility_exit_with_message( "ERROR!  To use crystal refinement compile Rosetta with extras=python." );
#endif
	return "";
}

}
}
}
