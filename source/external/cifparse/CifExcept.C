/*
FILE:     CifExcept.C
*/
/*
VERSION:  7.105
*/
/*
DATE:     10/21/2013
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2013 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.
*/
/*
               RCSB PDB SOFTWARE LICENSE AGREEMENT

BY CLICKING THE ACCEPTANCE BUTTON OR INSTALLING OR USING 
THIS "SOFTWARE, THE INDIVIDUAL OR ENTITY LICENSING THE  
SOFTWARE ("LICENSEE") IS CONSENTING TO BE BOUND BY AND IS 
BECOMING A PARTY TO THIS AGREEMENT.  IF LICENSEE DOES NOT 
AGREE TO ALL OF THE TERMS OF THIS AGREEMENT
THE LICENSEE MUST NOT INSTALL OR USE THE SOFTWARE.

1. LICENSE AGREEMENT

This is a license between you ("Licensee") and the Protein Data Bank (PDB) 
at Rutgers, The State University of New Jersey (hereafter referred to 
as "RUTGERS").   The software is owned by RUTGERS and protected by 
copyright laws, and some elements are protected by laws governing 
trademarks, trade dress and trade secrets, and may be protected by 
patent laws. 

2. LICENSE GRANT

RUTGERS grants you, and you hereby accept, non-exclusive, royalty-free 
perpetual license to install, use, modify, prepare derivative works, 
incorporate into other computer software, and distribute in binary 
and source code format, or any derivative work thereof, together with 
any associated media, printed materials, and on-line or electronic 
documentation (if any) provided by RUTGERS (collectively, the "SOFTWARE"), 
subject to the following terms and conditions: (i) any distribution 
of the SOFTWARE shall bind the receiver to the terms and conditions 
of this Agreement; (ii) any distribution of the SOFTWARE in modified 
form shall clearly state that the SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  

2. COPYRIGHT; RETENTION OF RIGHTS.  

The above license grant is conditioned on the following: (i) you must 
reproduce all copyright notices and other proprietary notices on any 
copies of the SOFTWARE and you must not remove such notices; (ii) in 
the event you compile the SOFTWARE, you will include the copyright 
notice with the binary in such a manner as to allow it to be easily 
viewable; (iii) if you incorporate the SOFTWARE into other code, you 
must provide notice that the code contains the SOFTWARE and include 
a copy of the copyright notices and other proprietary notices.  All 
copies of the SOFTWARE shall be subject to the terms of this Agreement.  

3. NO MAINTENANCE OR SUPPORT; TREATMENT OF ENHANCEMENTS 

RUTGERS is under no obligation whatsoever to: (i) provide maintenance 
or support for the SOFTWARE; or (ii) to notify you of bug fixes, patches, 
or upgrades to the features, functionality or performance of the 
SOFTWARE ("Enhancements") (if any), whether developed by RUTGERS 
or third parties.  If, in its sole discretion, RUTGERS makes an 
Enhancement available to you and RUTGERS does not separately enter 
into a written license agreement with you relating to such bug fix, 
patch or upgrade, then it shall be deemed incorporated into the SOFTWARE 
and subject to this Agreement. You are under no obligation whatsoever 
to provide any Enhancements to RUTGERS or the public that you may 
develop over time; however, if you choose to provide your Enhancements 
to RUTGERS, or if you choose to otherwise publish or distribute your 
Enhancements, in source code form without contemporaneously requiring 
end users or RUTGERS to enter into a separate written license agreement 
for such Enhancements, then you hereby grant RUTGERS a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare
derivative works, incorporate into the SOFTWARE or other computer
software, distribute, and sublicense your Enhancements or derivative
works thereof, in binary and source code form.

4. FEES.  There is no license fee for the SOFTWARE.  If Licensee
wishes to receive the SOFTWARE on media, there may be a small charge
for the media and for shipping and handling.  Licensee is
responsible for any and all taxes.

5. TERMINATION.  Without prejudice to any other rights, Licensor
may terminate this Agreement if Licensee breaches any of its terms
and conditions.  Upon termination, Licensee shall destroy all
copies of the SOFTWARE.

6. PROPRIETARY RIGHTS.  Title, ownership rights, and intellectual
property rights in the Product shall remain with RUTGERS.  Licensee 
acknowledges such ownership and intellectual property rights and will 
not take any action to jeopardize, limit or interfere in any manner 
with RUTGERS' ownership of or rights with respect to the SOFTWARE.  
The SOFTWARE is protected by copyright and other intellectual 
property laws and by international treaties.  Title and related 
rights in the content accessed through the SOFTWARE is the property 
of the applicable content owner and is protected by applicable law.  
The license granted under this Agreement gives Licensee no rights to such
content.

7. DISCLAIMER OF WARRANTY.  THE SOFTWARE IS PROVIDED FREE OF 
CHARGE, AND, THEREFORE, ON AN "AS IS" BASIS, WITHOUT WARRANTY OF 
ANY KIND, INCLUDING WITHOUT LIMITATION THE WARRANTIES THAT IT 
IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE 
OR NON-INFRINGING.  THE ENTIRE RISK AS TO THE QUALITY AND 
PERFORMANCE OF THE SOFTWARE IS BORNE BY LICENSEE.  SHOULD THE 
SOFTWARE PROVE DEFECTIVE IN ANY RESPECT, THE LICENSEE AND NOT 
LICENSOR ASSUMES THE ENTIRE COST OF ANY SERVICE AND REPAIR.  
THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL PART OF 
THIS AGREEMENT.  NO USE OF THE PRODUCT IS AUTHORIZED HEREUNDER 
EXCEPT UNDER THIS DISCLAIMER.

8. LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY
APPLICABLE LAW,  IN NO EVENT WILL LICENSOR BE LIABLE FOR ANY 
INDIRECT, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
OUT OF THE USE OF OR INABILITY TO USE THE SOFTWARE, INCLUDING, 
WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK 
STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL 
OTHER COMMERCIAL DAMAGES OR LOSSES, EVEN IF ADVISED OF THE
POSSIBILITY THEREOF. 
*/


/*!
** \file CifExcept.C
**
** \brief Implementation file for CifExcept class.
*/


#include "CifExcept.h"


using std::string;


bool CifExcept::CanBeUnknown(const string& itemName)
{
    // Optional, parent item, has value "?"
    if (itemName == "_atom_site.pdbx_PDB_ins_code")
        return (true);

    return (false);
}


bool CifExcept::CanBeInapplicable(const string& itemName)
{
    // Non-key, mandatory, parent; It seems that this item is always set to '.'
    if (itemName == "_atom_site.label_alt_id")
        return (true);
    // Non-key, mandatory, parent and child; It is '.' only if group_PDB is
    // set to HETATM, for other group_PDB it has some numeric value.
    else if (itemName == "_atom_site.label_seq_id")
        return (true);

    // Non-key, mandatory, child; It seems that this item is always set to '.'
    else if (itemName == "_struct_conn.ptnr1_label_seq_id")
        return (true);

    // Non-key, mandatory, child; It seems that this item is randomly set to '.'
    // depending on some conditions.
    else if (itemName == "_struct_conn.ptnr2_label_seq_id")
        return (true);

    // Non-key, mandatory, child; It seems that this item is always set to '.'
    else if (itemName == "_struct_site_gen.label_atom_id")
        return (true);

    // Non-key, mandatory, child; It seems that this item is randomly set to '.'
    // depending on some conditions.
    else if (itemName == "_struct_site_gen.label_seq_id")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item is always set
    // to '.'
    else if (itemName == "_atom_site_anisotrop.pdbx_label_alt_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is randomly set to '.',  depending on some conditions.
    else if (itemName == "_chem_comp.mon_nstd_flag")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_database_status.SG_entry")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_nonpoly_scheme.pdb_ins_code")
        return (true);

    // Non-key, non-mandatory, parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_poly_seq_scheme.pdb_ins_code")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_reflns.observed_criterion_sigma_I")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ref.biol_id")
        return (true);

    // Non-key, mandatory, child; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_site_gen.label_alt_id")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item
    // is randomly set to '.',  depending on some conditions.
    else if (itemName == "_atom_site_anisotrop.pdbx_label_seq_id")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_mon_prot_cis.label_alt_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_database_related.details")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine.pdbx_ls_sigma_F")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_refine_tls_group.beg_label_seq_id")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_refine_tls_group.beg_label_asym_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_refine_tls_group.end_label_seq_id")
        return (true);

    // Non-key, non-mandatory, child; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_refine_tls_group.end_label_asym_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_software.version")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine_ls_restr_ncs.pdbx_type")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_diffrn_radiation_wavelength.wavelength")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_em_imaging.specimen_holder_model")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_em_imaging.specimen_holder_type")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine_hist.d_res_low")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_exptl_crystal_grow.pdbx_pH_range")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.beg_label_alt_id")
        return (true);

    // Non-key, non-mandatory, child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.beg_label_asym_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.beg_label_comp_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.beg_label_seq_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.end_label_alt_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.end_label_asym_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.end_label_comp_id")
        return (true);

    // Non-key, non-mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_struct_ncs_dom_lim.end_label_seq_id")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine_hist.d_res_high")
        return (true);

    // Non-key, mandatory, child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_em_sample_preparation.entity_assembly_id")
        return (true);

    // Non-key, mandatory, child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine_ls_restr_ncs.pdbx_auth_asym_id")
        return (true);

    // Non-key, mandatory, child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_remediation_atom_site_mapping.label_seq_id")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_phasing_MIR_der.d_res_low")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_phasing_set.pdbx_d_res_low")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_phasing_MIR_der.d_res_high")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_phasing_set.pdbx_d_res_high")
        return (true);

    // Non-key, mandatory, child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_phasing_MIR_der.native_set_id")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_reflns_shell.d_res_high")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_pdbx_reflns_twin.fraction")
        return (true);

    // Non-key, mandatory, non-child, non-parent; It seems that this item
    // is always set to '.'
    else if (itemName == "_refine.ls_d_res_high")
        return (true);

    return (false);
}


bool CifExcept::IsBadParentRelation(const string& itemName)
{
#ifdef BAD_RELATION
chem_comp_link                           1  '_chem_comp_link.type_comp_1'                         '_chem_comp.type'                           chem_comp
chem_comp_link                           2  '_chem_comp_link.type_comp_2'                         '_chem_comp.type'                           chem_comp
#endif

    if (itemName == "_chem_comp.type")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_ref_seq_deletion             1  '_pdbx_struct_ref_seq_deletion.asym_id'               '_pdbx_poly_seq_scheme.asym_id'             pdbx_poly_seq_scheme
#endif

    else if (itemName == "_pdbx_poly_seq_scheme.asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_sheet_hbond                  4  '_pdbx_struct_sheet_hbond.range_id_1'                 '_struct_sheet_range.id'                    struct_sheet_range
pdbx_struct_sheet_hbond                  5  '_pdbx_struct_sheet_hbond.range_id_2'                 '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_hbond                       6  '_struct_sheet_hbond.range_id_1'                      '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_hbond                       7  '_struct_sheet_hbond.range_id_2'                      '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_order                       2  '_struct_sheet_order.range_id_1'                      '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_order                       3  '_struct_sheet_order.range_id_2'                      '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_topology                    2  '_struct_sheet_topology.range_id_1'                   '_struct_sheet_range.id'                    struct_sheet_range
struct_sheet_topology                    3  '_struct_sheet_topology.range_id_2'                   '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_range.id")
        return (true);

#ifdef BAD_RELATION
atom_site                                5  '_atom_site.label_atom_id'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_angle                          1  '_chem_comp_angle.atom_id_1'                          '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_angle                          2  '_chem_comp_angle.atom_id_2'                          '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_angle                          3  '_chem_comp_angle.atom_id_3'                          '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_bond                           1  '_chem_comp_bond.atom_id_1'                           '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_bond                           2  '_chem_comp_bond.atom_id_2'                           '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_chir                           2  '_chem_comp_chir.atom_id'                             '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_chir_atom                      2  '_chem_comp_chir_atom.atom_id'                        '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_plane_atom                     2  '_chem_comp_plane_atom.atom_id'                       '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_tor                            1  '_chem_comp_tor.atom_id_1'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_tor                            2  '_chem_comp_tor.atom_id_2'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_tor                            3  '_chem_comp_tor.atom_id_3'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
chem_comp_tor                            4  '_chem_comp_tor.atom_id_4'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
#endif

    else if (itemName == "_chem_comp_atom.atom_id")
        return (true);

    return (false);
}


bool CifExcept::IsBadChildRelation(const string& itemName)
{
#ifdef BAD_RELATION
refine_ls_restr_ncs                      4  '_refine_ls_restr_ncs.pdbx_auth_asym_id'              '_atom_site.auth_asym_id'                   atom_site
#endif

    if (itemName == "_refine_ls_restr_ncs.pdbx_auth_asym_id")
        return (true);

#ifdef BAD_RELATION
chem_comp_link                           1  '_chem_comp_link.type_comp_1'                         '_chem_comp.type'                           chem_comp
#endif

    else if (itemName == "_chem_comp_link.type_comp_1")
        return (true);

#ifdef BAD_RELATION
chem_comp_link                           2  '_chem_comp_link.type_comp_2'                         '_chem_comp.type'                           chem_comp
#endif

    else if (itemName == "_chem_comp_link.type_comp_2")
        return (true);

#ifdef BAD_RELATION
pdbx_nonpoly_scheme                      1  '_pdbx_nonpoly_scheme.asym_id'                        '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_pdbx_nonpoly_scheme.asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_unobs_or_zero_occ_atoms             1  '_pdbx_unobs_or_zero_occ_atoms.PDB_model_num'         '_atom_site.pdbx_PDB_model_num'             atom_site
#endif

    else if (itemName == "_pdbx_unobs_or_zero_occ_atoms.PDB_model_num")
        return (true);

#ifdef BAD_RELATION
pdbx_unobs_or_zero_occ_residues          1  '_pdbx_unobs_or_zero_occ_residues.PDB_model_num'      '_atom_site.pdbx_PDB_model_num'             atom_site
#endif

    else if (itemName == "_pdbx_unobs_or_zero_occ_residues.PDB_model_num")
        return (true);

#ifdef BAD_RELATION
pdbx_nonpoly_scheme                      1  '_pdbx_nonpoly_scheme.mon_id'                         '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_pdbx_nonpoly_scheme.mon_id")
        return (true);

#ifdef BAD_RELATION
pdbx_nonpoly_scheme                      1  '_pdbx_nonpoly_scheme.entity_id'                      '_atom_site.label_entity_id'                atom_site
#endif

    else if (itemName == "_pdbx_nonpoly_scheme.entity_id")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_ref_seq_deletion             1  '_pdbx_struct_ref_seq_deletion.asym_id'               '_pdbx_poly_seq_scheme.asym_id'             pdbx_poly_seq_scheme
#endif

    else if (itemName == "_pdbx_struct_ref_seq_deletion.asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_sheet_hbond                  4  '_pdbx_struct_sheet_hbond.range_id_1'                 '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_pdbx_struct_sheet_hbond.range_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_sheet_hbond                  5  '_pdbx_struct_sheet_hbond.range_id_2'                 '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_pdbx_struct_sheet_hbond.range_id_2")
        return (true);

#ifdef BAD_RELATION
atom_site                                5  '_atom_site.label_atom_id'                            '_chem_comp_atom.atom_id'                   chem_comp_atom
#endif

    else if (itemName == "_atom_site.label_atom_id")
        return (true);

#ifdef BAD_RELATION
struct_site_gen                          1  '_struct_site_gen.label_asym_id'                      '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_struct_site_gen.label_asym_id")
        return (true);

#ifdef BAD_RELATION
struct_site_gen                          1  '_struct_site_gen.label_comp_id'                      '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_struct_site_gen.label_comp_id")
        return (true);

#ifdef BAD_RELATION
struct_mon_prot_cis                      1  '_struct_mon_prot_cis.label_asym_id'                  '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_struct_mon_prot_cis.label_asym_id")
        return (true);

#ifdef BAD_RELATION
struct_mon_prot_cis                      1  '_struct_mon_prot_cis.label_comp_id'                  '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_struct_mon_prot_cis.label_comp_id")
        return (true);

#ifdef BAD_RELATION
struct_conn                              1  '_struct_conn.ptnr1_label_asym_id'                    '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_struct_conn.ptnr1_label_asym_id")
        return (true);

#ifdef BAD_RELATION
struct_conn                              1  '_struct_conn.ptnr1_label_comp_id'                    '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_struct_conn.ptnr1_label_comp_id")
        return (true);

#ifdef BAD_RELATION
struct_conn                              2  '_struct_conn.ptnr2_label_asym_id'                    '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_struct_conn.ptnr2_label_asym_id")
        return (true);

#ifdef BAD_RELATION
struct_conn                              2  '_struct_conn.ptnr2_label_comp_id'                    '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_struct_conn.ptnr2_label_comp_id")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_planes_atom                1  '_pdbx_validate_planes_atom.auth_asym_id'             '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_planes_atom.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_peptide_omega              2  '_pdbx_validate_peptide_omega.auth_asym_id_2'         '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_peptide_omega.auth_asym_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_peptide_omega              1  '_pdbx_validate_peptide_omega.auth_asym_id_1'         '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_peptide_omega.auth_asym_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_close_contact              2  '_pdbx_validate_close_contact.auth_asym_id_2'         '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_close_contact.auth_asym_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_rmsd_angle                 2  '_pdbx_validate_rmsd_angle.auth_asym_id_2'            '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_rmsd_angle.auth_asym_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_rmsd_angle                 3  '_pdbx_validate_rmsd_angle.auth_asym_id_3'            '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_rmsd_angle.auth_asym_id_3")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_rmsd_bond                  2  '_pdbx_validate_rmsd_bond.auth_asym_id_2'             '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_rmsd_bond.auth_asym_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_symm_contact               2  '_pdbx_validate_symm_contact.auth_asym_id_2'          '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_symm_contact.auth_asym_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_close_contact              1  '_pdbx_validate_close_contact.auth_asym_id_1'         '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_close_contact.auth_asym_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_rmsd_angle                 1  '_pdbx_validate_rmsd_angle.auth_asym_id_1'            '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_rmsd_angle.auth_asym_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_rmsd_bond                  1  '_pdbx_validate_rmsd_bond.auth_asym_id_1'             '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_rmsd_bond.auth_asym_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_symm_contact               1  '_pdbx_validate_symm_contact.auth_asym_id_1'          '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_symm_contact.auth_asym_id_1")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_chiral                     1  '_pdbx_validate_chiral.auth_asym_id'                  '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_chiral.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_main_chain_plane           1  '_pdbx_validate_main_chain_plane.auth_asym_id'        '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_main_chain_plane.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_planes                     1  '_pdbx_validate_planes.auth_asym_id'                  '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_planes.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_validate_torsion                    1  '_pdbx_validate_torsion.auth_asym_id'                 '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_validate_torsion.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_struct_mod_residue                  1  '_pdbx_struct_mod_residue.auth_asym_id'               '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_pdbx_struct_mod_residue.auth_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_domain_range                        1  '_pdbx_domain_range.beg_label_asym_id'                '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_pdbx_domain_range.beg_label_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_domain_range                        1  '_pdbx_domain_range.beg_label_comp_id'                '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_pdbx_domain_range.beg_label_comp_id")
        return (true);

#ifdef BAD_RELATION
pdbx_domain_range                        2  '_pdbx_domain_range.end_label_asym_id'                '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_pdbx_domain_range.end_label_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_domain_range                        2  '_pdbx_domain_range.end_label_comp_id'                '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_pdbx_domain_range.end_label_comp_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_monomer                     1  '_pdbx_feature_monomer.label_asym_id'                 '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_pdbx_feature_monomer.label_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_monomer                     1  '_pdbx_feature_monomer.label_comp_id'                 '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_pdbx_feature_monomer.label_comp_id")
        return (true);


#ifdef BAD_RELATION
pdbx_refine_component                    1  '_pdbx_refine_component.label_asym_id'                '_atom_site.label_asym_id'                  atom_site
#endif

    else if (itemName == "_pdbx_refine_component.label_asym_id")
        return (true);

#ifdef BAD_RELATION
pdbx_refine_component                    1  '_pdbx_refine_component.label_comp_id'                '_atom_site.label_comp_id'                  atom_site
#endif

    else if (itemName == "_pdbx_refine_component.label_comp_id")
        return (true);

#ifdef BAD_RELATION
struct_ref_seq_dif                       2  '_struct_ref_seq_dif.mon_id'                          '_entity_poly_seq.mon_id'                   entity_poly_seq
#endif

    else if (itemName == "_struct_ref_seq_dif.mon_id")
        return (true);

#ifdef BAD_RELATION
struct_ref_seq_dif                       2  '_struct_ref_seq_dif.seq_num'                         '_entity_poly_seq.num'                      entity_poly_seq
#endif

    else if (itemName == "_struct_ref_seq_dif.seq_num")
        return (true);

#ifdef BAD_RELATION
refine_ls_restr_ncs                      3  '_refine_ls_restr_ncs.pdbx_ens_id'                    '_struct_ncs_dom.pdbx_ens_id'               struct_ncs_dom
#endif

    else if (itemName == "_refine_ls_restr_ncs.pdbx_ens_id")
        return (true);

#ifdef BAD_RELATION
atom_site                                11 '_atom_site.pdbx_ncs_dom_id'                          '_struct_ncs_dom.id'                        struct_ncs_dom
#endif

    else if (itemName == "_atom_site.pdbx_ncs_dom_id")
        return (true);

#ifdef BAD_RELATION
struct_ncs_ens_gen                       1  '_struct_ncs_ens_gen.dom_id_1'                        '_struct_ncs_dom.id'                        struct_ncs_dom
#endif

    else if (itemName == "_struct_ncs_ens_gen.dom_id_1")
        return (true);

#ifdef BAD_RELATION
struct_ncs_ens_gen                       2  '_struct_ncs_ens_gen.dom_id_2'                        '_struct_ncs_dom.id'                        struct_ncs_dom
#endif

    else if (itemName == "_struct_ncs_ens_gen.dom_id_2")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_assembly                    2  '_pdbx_feature_assembly.feature_software_id'          '_software.name'                            software
#endif

    else if (itemName == "_pdbx_feature_assembly.feature_software_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_domain                      3  '_pdbx_feature_domain.feature_software_id'            '_software.name'                            software
#endif

    else if (itemName == "_pdbx_feature_domain.feature_software_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_entry                       2  '_pdbx_feature_entry.feature_software_id'             '_software.name'                            software
#endif

    else if (itemName == "_pdbx_feature_entry.feature_software_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_monomer                     3  '_pdbx_feature_monomer.feature_software_id'           '_software.name'                            software
#endif

    else if (itemName == "_pdbx_feature_monomer.feature_software_id")
        return (true);

#ifdef BAD_RELATION
pdbx_feature_sequence_range              3  '_pdbx_feature_sequence_range.feature_software_id'    '_software.name'                            software
#endif

    else if (itemName == "_pdbx_feature_sequence_range.feature_software_id")
        return (true);

#ifdef BAD_RELATION
refine_ls_restr_ncs                      4  '_refine_ls_restr_ncs.pdbx_auth_asym_id'              '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_refine_ls_restr_ncs.pdbx_auth_asym_id")
        return (true);

#ifdef BAD_RELATION
em_sample_preparation                    4  '_em_sample_preparation.entity_assembly_id'   '_em_entity_assembly.id'      em_entity_assembly
#endif

    else if (itemName == "_em_sample_preparation.entity_assembly_id")
        return (true);

#ifdef BAD_RELATION
refine_ls_restr_ncs                      4  '_refine_ls_restr_ncs.pdbx_auth_asym_id'              '_atom_site.auth_asym_id'                   atom_site
#endif

    else if (itemName == "_refine_ls_restr_ncs.pdbx_auth_asym_id")
        return (true);

#ifdef BAD_RELATION
struct_sheet_hbond                       6  '_struct_sheet_hbond.range_id_1'                      '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_hbond.range_id_1")
        return (true);

#ifdef BAD_RELATION
struct_sheet_hbond                       7  '_struct_sheet_hbond.range_id_2'                      '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_hbond.range_id_2")
        return (true);

#ifdef BAD_RELATION
struct_sheet_order                       2  '_struct_sheet_order.range_id_1'                      '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_order.range_id_1")
        return (true);

#ifdef BAD_RELATION
struct_sheet_order                       3  '_struct_sheet_order.range_id_2'                      '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_order.range_id_2")
        return (true);

#ifdef BAD_RELATION
struct_sheet_topology                    2  '_struct_sheet_topology.range_id_1'                   '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_topology.range_id_1")
        return (true);

#ifdef BAD_RELATION
struct_sheet_topology                    3  '_struct_sheet_topology.range_id_2'                   '_struct_sheet_range.id'                    struct_sheet_range
#endif

    else if (itemName == "_struct_sheet_topology.range_id_2")
        return (true);

#ifdef BAD_RELATION
struct_ref_seq                           1  '_struct_ref_seq.seq_align_beg'                       '_entity_poly_seq.num'                      entity_poly_seq
#endif

    else if (itemName == "_struct_ref_seq.seq_align_beg")
        return (true);

#ifdef BAD_RELATION
struct_ref_seq                           2  '_struct_ref_seq.seq_align_end'                       '_entity_poly_seq.num'                      entity_poly_seq
#endif

    else if (itemName == "_struct_ref_seq.seq_align_end")
        return (true);

#ifdef BAD_RELATION
entity_link                              4  '_entity_link.entity_seq_num_1'                       '_entity_poly_seq.num'                      entity_poly_seq
#endif

    else if (itemName == "_entity_link.entity_seq_num_1")
        return (true);

#ifdef BAD_RELATION
entity_link                              5  '_entity_link.entity_seq_num_2'                       '_entity_poly_seq.num'                      entity_poly_seq
#endif

    else if (itemName == "_entity_link.entity_seq_num_2")
        return (true);

// VLAD - did not finish everything. There are a lot of items. Must filter
// per parent.

    return (false);
}
