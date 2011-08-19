#ifdef MAC
#include <OpenCL/cl_platform.h>
#else
#include <CL/cl_platform.h>
#endif


#define x1_ca     -1.5294880 // ACDEFGHIKLMNPQRSTVWY
#define x1_cb      0.0000000 // ACDEFHIKLMNQRSTVWY
#define x1_cb_p   -0.0119584 
#define x1_cd1    -0.2442786 // FY
#define x1_cd1_w  -0.0962845 
#define x1_cd2     1.9679581 // FY
#define x1_cd2_h   1.8770485 
#define x1_cd2_l   1.4893895 
#define x1_cd2_w   1.9952924 
#define x1_cd     -0.4783638 // EIKLQR
#define x1_cd_p   -1.0268733 
#define x1_ce1     0.2943790 // FY
#define x1_ce1_h   0.5944039 
#define x1_ce2     2.5065409 // FY
#define x1_ce2_w   2.0616894 
#define x1_ce3_w   3.1824268 
#define x1_ce_k   -1.8711717 
#define x1_ce_m   -2.1632737 
#define x1_cg2     0.5326465 // ITV
#define x1_cg      0.5915642 // DEFHIKLMNQRVWY
#define x1_cg_p    0.1567080 
#define x1_ch2_w   4.4375384 
#define x1_cz2_w   3.2711149 
#define x1_cz3_w   4.3952531 
#define x1_cz      1.6750406 // FY
#define x1_cz_r   -2.0808292 
#define x1_c      -2.0462234 // ACDEFGHIKLMNPQRSTVWY
#define x1_nd1_h  -0.1880577 
#define x1_nd2_n   1.8976920 
#define x1_ne1_w   0.7775292 
#define x1_ne2_h   1.8507970 
#define x1_ne2_q  -0.0574527 
#define x1_ne_r   -1.8147512 
#define x1_nh1_r  -1.0945695 
#define x1_nh2_r  -3.3333551 
#define x1_nv_p   -2.0181350 
#define x1_nz_k   -1.8249921 
#define x1_n      -2.0400342 // ACDEFGHIKLMNPQRSTVWY
#define x1_od1    -0.1517762 // DN
#define x1_od2_d   1.7939259 
#define x1_oe1    -1.6545054 // EQ
#define x1_oe2_e  -0.1345519 
#define x1_og      0.5021427 // ST
#define x1_oh_y    2.2115206 
#define x1_o      -2.3056024 // ACDEFGHIKLMNPQRSTVWY
#define x1_sd_m   -0.6700778 
#define x1_sg_c    0.7385894 

#define y1_ca      0.0000000 // ACDEFGHIKLMNPQRSTVWY
#define y1_cb     -0.0000000 // ACDEFHIKLMNQRSTVWY
#define y1_cb_p    0.1993893 
#define y1_cd1     2.5042751 // FY
#define y1_cd1_w   2.5735523 
#define y1_cd2     1.5674440 // FY
#define y1_cd2_h   1.8213088 
#define y1_cd2_l   1.5685415 
#define y1_cd2_w   1.7533594 
#define y1_cd      2.4630111 // EIKLQR
#define y1_cd_p    2.3111539 
#define y1_ce1     3.7764912 // FY
#define y1_ce1_h   3.6003141 
#define y1_ce2     2.8394833 // FY
#define y1_ce2_w   3.1491190 
#define y1_ce3_w   1.0105477 
#define y1_ce_k    1.8503952 
#define y1_ce_m    1.6985493 
#define y1_cg2    -0.7466988 // ITV
#define y1_cg      1.3971722 // DEFHIKLMNQRVWY
#define y1_cg_p    1.4764299 
#define y1_ch2_w   3.0607742 
#define y1_cz2_w   3.8259826 
#define y1_cz3_w   1.6893180 
#define y1_cz      3.9408344 // FY
#define y1_cz_r    0.5692133 
#define y1_c      -0.7812559 // ACDEFGHIKLMNPQRSTVWY
#define y1_nd1_h   2.5348379 
#define y1_nd2_n   1.4782518 
#define y1_ne1_w   3.6325178 
#define y1_ne2_h   3.1942679 
#define y1_ne2_q   3.7226568 
#define y1_ne_r    1.8899669 
#define y1_nh1_r  -0.2999949 
#define y1_nh2_r   0.1474556 
#define y1_nv_p    1.2939697 
#define y1_nz_k    0.3630040 
#define y1_n       1.3656901 // ACDEFGHIKLMNPQRSTVWY
#define y1_od1     2.3674910 // DN
#define y1_od2_d   1.5112570 
#define y1_oe1     2.1388335 // EQ
#define y1_oe2_e   3.6216168 
#define y1_og      1.3081245 // ST
#define y1_oh_y    5.2079019 
#define y1_o      -0.2101090 // ACDEFGHIKLMNPQRSTVWY
#define y1_sd_m    2.6864031 
#define y1_sg_c    1.6511372 

#define z1_ca      0.0000000 // ACDEFGHIKLMNPQRSTVWY
#define z1_cb      0.0000000 // ACDEFHIKLMNQRSTVWY
#define z1_cb_p    0.0661247 
#define z1_cd1     0.0000000 // FY
#define z1_cd1_w   0.0000028 
#define z1_cd2     0.0000000 // FY
#define z1_cd2_h  -0.0009351 
#define z1_cd2_l   1.2162614 
#define z1_cd2_w   0.0008309 
#define z1_cd      0.0000007 // EIKLQR
#define z1_cd_p    0.4052198 
#define z1_ce1     0.0000000 // FY
#define z1_ce1_h   0.0000084 
#define z1_ce2     0.0000000 // FY
#define z1_ce2_w   0.0013971 
#define z1_ce3_w   0.0019356 
#define z1_ce_k    0.0000013 
#define z1_ce_m    0.0000004 
#define z1_cg2    -1.2132623 // ITV
#define z1_cg      0.0000000 // DEFHIKLMNQRVWY
#define z1_cg_p    0.8161961 
#define z1_ch2_w   0.0029114 
#define z1_cz2_w   0.0008706 
#define z1_cz3_w   0.0027443 
#define z1_cz      0.0000000 // FY
#define z1_cz_r    0.0000010 
#define z1_c       1.2012238 // ACDEFGHIKLMNPQRSTVWY
#define z1_nd1_h   0.0000002 
#define z1_nd2_n   0.0000000 
#define z1_ne1_w  -0.0006771 
#define z1_ne2_h  -0.0009374 
#define z1_ne2_q   0.0000007 
#define z1_ne_r    0.0000011 
#define z1_nh1_r   0.0000004 
#define z1_nh2_r   0.0000015 
#define z1_nv_p    0.0057931 
#define z1_nz_k    0.0000008 
#define z1_n       0.0000000 // ACDEFGHIKLMNPQRSTVWY
#define z1_od1     0.0000000 // DN
#define z1_od2_d   0.0000000 
#define z1_oe1     0.0000011 // EQ
#define z1_oe2_e   0.0000007 
#define z1_og      0.0000000 // ST
#define z1_oh_y    0.0000000 
#define z1_o       2.2604272 // ACDEFGHIKLMNPQRSTVWY
#define z1_sd_m    0.0000009 
#define z1_sg_c    0.0000000 

#define x2_ca     -0.5963369 // ACDEFGHIKLMNPQRSTVWY
#define x2_cb     -0.0000000 // ACDEFHIKLMNQRSTVWY
#define x2_cb_p    0.1789471 
#define x2_cd1     2.2108445 // FY
#define x2_cd1_w   2.3323410 
#define x2_cd2     2.2106900 // FY
#define x2_cd2_h   2.4090189 
#define x2_cd2_l   2.0251100 
#define x2_cd2_w   2.3925496 
#define x2_cd      2.0815779 // EIKLQR
#define x2_cd_p    1.7278786 
#define x2_ce1     3.5923965 // FY
#define x2_ce1_h   3.5471396 
#define x2_ce2     3.5920500 // FY
#define x2_ce2_w   3.7037366 
#define x2_ce3_w   2.1713794 
#define x2_ce_k    0.9743982 
#define x2_ce_m    0.7206807 
#define x2_cg2    -0.4799299 // ITV
#define x2_cg      1.5172470 // DEFHIKLMNQRVWY
#define x2_cg_p    1.4206848 
#define x2_ch2_w   4.5487107 
#define x2_cz2_w   4.7985800 
#define x2_cz3_w   3.2693045 
#define x2_cz      4.2820440 // FY
#define x2_cz_r   -0.2871352 
#define x2_c      -1.5172358 // ACDEFGHIKLMNPQRSTVWY
#define x2_nd1_h   2.2609087 
#define x2_nd2_n   2.1011602 
#define x2_ne1_w   3.6481941 
#define x2_ne2_h   3.6630869 
#define x2_ne2_q   3.4056457 
#define x2_ne_r    1.0328361 
#define x2_nh1_r  -0.7030185 
#define x2_nh2_r  -1.1638662 
#define x2_nv_p    0.4047081 
#define x2_nz_k   -0.3772760 
#define x2_n       0.4622142 // ACDEFGHIKLMNPQRSTVWY
#define x2_od1     2.1209516 // DN
#define x2_od2_d   2.0910957 
#define x2_oe1     1.3244862 // EQ
#define x2_oe2_e   3.2825415 
#define x2_og      1.4003816 // ST
#define x2_oh_y    5.6580058 
#define x2_o      -1.0924196 // ACDEFGHIKLMNPQRSTVWY
#define x2_sd_m    2.2125427 
#define x2_sg_c    1.8084373 

#define y2_ca      1.4084445 // ACDEFGHIKLMNPQRSTVWY
#define y2_cb      0.0000000 // ACDEFHIKLMNQRSTVWY
#define y2_cb_p    0.0887525 
#define y2_cd1     1.2013461 // FY
#define y2_cd1_w   1.0920749 
#define y2_cd2    -1.2010785 // FY
#define y2_cd2_h  -1.0183833 
#define y2_cd2_l  -0.7599558 
#define y2_cd2_w  -1.1537624 
#define y2_cd      1.4008173 // EIKLQR
#define y2_cd_p    1.8467097 
#define y2_ce1     1.2013461 // FY
#define y2_ce1_h   0.8563749 
#define y2_ce2    -1.2010785 // FY
#define y2_ce2_w  -0.6707077 
#define y2_ce3_w  -2.5365643 
#define y2_ce_k    2.4445438 
#define y2_ce_m    2.6543252 
#define y2_cg2    -0.7816257 // ITV
#define y2_cg      0.0000000 // DEFHIKLMNQRVWY
#define y2_cg_p    0.4313438 
#define y2_ch2_w  -2.8929773 
#define y2_cz2_w  -1.5205149 
#define y2_cz3_w  -3.3887598 
#define y2_cz     -0.0059738 // FY
#define y2_cz_r    2.1380849 
#define y2_c       1.5796792 // ACDEFGHIKLMNPQRSTVWY
#define y2_nd1_h   1.1614907 
#define y2_nd2_n  -1.1711487 
#define y2_ne1_w   0.7002982 
#define y2_ne2_h  -0.4589020 
#define y2_ne2_q   1.5043442 
#define y2_ne_r    2.4080171 
#define y2_nh1_r   0.8909794 
#define y2_nh2_r   3.1270457 
#define y2_nv_p    2.3629299 
#define y2_nz_k    1.8220952 
#define y2_n       2.4110593 // ACDEFGHIKLMNPQRSTVWY
#define y2_od1     1.0628332 // DN
#define y2_od2_d  -1.0627261 
#define y2_oe1     2.3574845 // EQ
#define y2_oe2_e   1.5359470 
#define y2_og      0.0476256 // ST
#define y2_oh_y   -0.0059759 
#define y2_o       2.0412172 // ACDEFGHIKLMNPQRSTVWY
#define y2_sd_m    1.6644581 
#define y2_sg_c   -0.0363704 

#define z2_ca      0.0000000 // ACDEFGHIKLMNPQRSTVWY
#define z2_cb      0.0000000 // ACDEFHIKLMNQRSTVWY
#define z2_cb_p    0.0661247 
#define z2_cd1     0.0000000 // FY
#define z2_cd1_w   0.0000028 
#define z2_cd2     0.0000000 // FY
#define z2_cd2_h  -0.0009351 
#define z2_cd2_l   1.2162614 
#define z2_cd2_w   0.0008309 
#define z2_cd      0.0000007 // EIKLQR
#define z2_cd_p    0.4052198 
#define z2_ce1     0.0000000 // FY
#define z2_ce1_h   0.0000084 
#define z2_ce2     0.0000000 // FY
#define z2_ce2_w   0.0013971 
#define z2_ce3_w   0.0019356 
#define z2_ce_k    0.0000013 
#define z2_ce_m    0.0000004 
#define z2_cg2    -1.2132623 // ITV
#define z2_cg      0.0000000 // DEFHIKLMNQRVWY
#define z2_cg_p    0.8161961 
#define z2_ch2_w   0.0029114 
#define z2_cz2_w   0.0008706 
#define z2_cz3_w   0.0027443 
#define z2_cz      0.0000000 // FY
#define z2_cz_r    0.0000010 
#define z2_c       1.2012238 // ACDEFGHIKLMNPQRSTVWY
#define z2_nd1_h   0.0000002 
#define z2_nd2_n   0.0000000 
#define z2_ne1_w  -0.0006771 
#define z2_ne2_h  -0.0009374 
#define z2_ne2_q   0.0000007 
#define z2_ne_r    0.0000011 
#define z2_nh1_r   0.0000004 
#define z2_nh2_r   0.0000015 
#define z2_nv_p    0.0057931 
#define z2_nz_k    0.0000008 
#define z2_n       0.0000000 // ACDEFGHIKLMNPQRSTVWY
#define z2_od1     0.0000000 // DN
#define z2_od2_d   0.0000000 
#define z2_oe1     0.0000011 // EQ
#define z2_oe2_e   0.0000007 
#define z2_og      0.0000000 // ST
#define z2_oh_y    0.0000000 
#define z2_o       2.2604272 // ACDEFGHIKLMNPQRSTVWY
#define z2_sd_m    0.0000009 
#define z2_sg_c    0.0000000 

#define r_ca      2.0000000 // ACDEFGHIKLMNPQRSTVWY
#define r_cb      2.0000000 // ACDEFHIKLMNQRSTVWY
#define r_cb_p    2.0000000 
#define r_cd1     2.0000000 // FY
#define r_cd1_w   2.0000000 
#define r_cd2     2.0000000 // FY
#define r_cd2_h   2.0000000 
#define r_cd2_l   2.0000000 
#define r_cd2_w   2.0000000 
#define r_cd      2.0000000 // EIKLQR
#define r_cd_p    2.0000000 
#define r_ce1     2.0000000 // FY
#define r_ce1_h   2.0000000 
#define r_ce2     2.0000000 // FY
#define r_ce2_w   2.0000000 
#define r_ce3_w   2.0000000 
#define r_ce_k    2.0000000 
#define r_ce_m    2.0000000 
#define r_cg2     2.0000000 // ITV
#define r_cg      2.0000000 // DEFHIKLMNQRVWY
#define r_cg_p    2.0000000 
#define r_ch2_w   2.0000000 
#define r_cz2_w   2.0000000 
#define r_cz3_w   2.0000000 
#define r_cz      2.0000000 // FY
#define r_cz_r    2.0000000 
#define r_c       2.0000000 // ACDEFGHIKLMNPQRSTVWY
#define r_nd1_h   1.7500000 
#define r_nd2_n   1.7500000 
#define r_ne1_w   1.7500000 
#define r_ne2_h   1.7500000 
#define r_ne2_q   1.7500000 
#define r_ne_r    1.7500000 
#define r_nh1_r   1.7500000 
#define r_nh2_r   1.7500000 
#define r_nv_p    1.7500000 
#define r_nz_k    1.7500000 
#define r_n       1.7500000 // ACDEFGHIKLMNPQRSTVWY
#define r_od1     1.5500000 // DN
#define r_od2_d   1.5500000 
#define r_oe1     1.5500000 // EQ
#define r_oe2_e   1.5500000 
#define r_og      1.5500000 // ST
#define r_oh_y    1.5500000 
#define r_o       1.5500000 // ACDEFGHIKLMNPQRSTVWY
#define r_sd_m    1.9000000 
#define r_sg_c    1.9000000 

enum UATOM {
  CA, // ACDEFGHIKLMNPQRSTVWY
  CB, // ACDEFHIKLMNQRSTVWY
  CB_P,
  CD1, // FY
  CD1_W,
  CD2, // FY
  CD2_H,
  CD2_L,
  CD2_W,
  CD, // EIKLQR
  CD_P,
  CE1, // FY
  CE1_H,
  CE2, // FY
  CE2_W,
  CE3_W,
  CE_K,
  CE_M,
  CG2, // ITV
  CG, // DEFHIKLMNQRVWY
  CG_P,
  CH2_W,
  CZ2_W,
  CZ3_W,
  CZ, // FY
  CZ_R,
  C, // ACDEFGHIKLMNPQRSTVWY
  ND1_H,
  ND2_N,
  NE1_W,
  NE2_H,
  NE2_Q,
  NE_R,
  NH1_R,
  NH2_R,
  NV_P,
  NZ_K,
  N, // ACDEFGHIKLMNPQRSTVWY
  OD1, // DN
  OD2_D,
  OE1, // EQ
  OE2_E,
  OG, // ST
  OH_Y,
  O, // ACDEFGHIKLMNPQRSTVWY
  SD_M,
  SG_C
};

