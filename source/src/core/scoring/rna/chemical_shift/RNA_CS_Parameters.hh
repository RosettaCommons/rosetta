// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_Parameters.hh
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)



#ifndef INCLUDED_core_scoring_rna_RNA_CS_Parameters_HH
#define INCLUDED_core_scoring_rna_RNA_CS_Parameters_HH

#include <core/types.hh>

#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>


#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


enum atomitem{ oshi = 1, xdir, ydir, zdir, csca, suga, rcl1, rcl2, rcl3, maca, maqx, maqw, maqy, maqz, marx, mary, marz, chco, chrg, last_atomdesc}; 
                   //1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19              20

        //const char* const atomdesc[] = {"OSHI", "XDIR", "YDIR", "ZDIR", "CSCA", "SUGA",
        //              "RCL1", "RCL2", "RCL3", "MACA", "MAQX", "MAQW",
        //              "MAQY", "MAQZ", "MARX", "MARY", "MARZ", "CHCO",
        //              "CHRG"};
        // Numerical data in the description of atoms (add 4 letter abbs
        // for new data, used by new calculations in this list)
        // and add the same abbs in the same order in the enum line above
        // (before last_atomdesc)

        // OSHI, offset of the chemical shift calculation
        // XDIR, YDIR, ZDIR x-axis, y-axis, z-axis
        // from mean of all atoms tagged 1 to mean of all atoms tagged 2
        // CSCA has to be nonzero if the chemical shift must be calculated for the atom
        // SUGA has to be true if an atom counts to the sugar part
        // Magnetic anisotropy and Ring current are not calculated for base
        // protons they are however for the sugar protons.
        // MA is data for MAGNETIC ANISOTROPY MODULE
        // MACA has to be nonzero if this atom gives no contribution to MA
        // MAQ/MAR are the Q and R tensors. MAQW = the Qxy tensor!
        // CHCO must be the atom number to which te hydrogen is attached
        //   if the csd due to the electrical field has to be calculated for
        //   this hydrogen atom
        // CHRG charge information for an atom



////////////////////////////////////////////////////////////////////////////////////////////////////

class RNA_CS_residue_parameters : public utility::pointer::ReferenceCount { 


	public:

		RNA_CS_residue_parameters( chemical::AA const & res_aa ); 

		~RNA_CS_residue_parameters(); 

		std::string const base_name() const; 

		Size num_rings() const; 

		Real ring_intensity( Size const ring_ID ) const; 

		Real ring_radius( Size const ring_ID ) const; 

		Real ring_height( Size const ring_ID ) const; 

		Size get_atomnames_size() const; 

		std::string const get_atomname( Size const count ) const; 

		//Undefinded, commenting out to fix PyRosetta build utility::vector1< std::string > const & get_ring_center_representative_atoms() const;

		Real atom_data( Size const atom, atomitem const item ) const; 

		Real ring_current_coeff() const; 

		Real magentic_anisotropy_r_coeff() const; 

		Real magentic_anisotropy_q_coeff() const; 

		chemical::AA aa() const; 

	private:


	private:

		chemical::AA const res_aa_; //chemical::AA is enum_type: 	na_rgu (25),	na_rad (26),	na_rcy (27),	na_ura (28). Usage: chemical::AA const & res_aa =  pose1.residue(moving_res_1).aa();,
		Size	 const maxatoms_; //Max number of atoms in each nucleotide.
		core::Real const RCCO_; //This constant does not depend on specific RNA base, but moved here for convenient.
		core::Real const MACQ_; //This constant does not depend on specific RNA base, but moved here for convenient.
		core::Real const MACR_; //This constant does not depend on specific RNA base, but moved here for convenient.

		std::string BASE_; 
		core::Size num_rings_; 


		utility::vector1< core::Real > ring_intensity_; //RCI ring current intensity (relative to benzene?)
		utility::vector1< core::Real > ring_radius_;    //RCR: Radius of the ring (Angstrom)
		utility::vector1< core::Real > ring_height_;     //RCH: Distance of the ring current loops to the molecular plane (Angstrom)

		utility::vector1< std::string > atomnames_; 

		utility::vector1< utility::vector1< core::Real > > realatomdata_;  //[maxdiffbases][maxatoms][last_atomdesc];

}; 

////////////////////////////////////////////////////////////////////////////////////////////////////

class RNA_CS_parameters : public utility::pointer::ReferenceCount { 


	//COMM RNA data set without charge

	public:

		RNA_CS_parameters(); 


	  ~RNA_CS_parameters(); 

	public:

		RNA_CS_residue_parameters const &
		get_RNA_CS_residue_parameters( chemical::AA const res_aa ) const; 


	private:
		std::string COMM_; 
		RNA_CS_residue_parameters const CS_RAD_params_; 
		RNA_CS_residue_parameters const CS_RGU_params_; 
		RNA_CS_residue_parameters const CS_RCY_params_; 
		RNA_CS_residue_parameters const CS_URA_params_; 

}; 




} //chemical_shift
} //rna
} //scoring
} //core


#endif
