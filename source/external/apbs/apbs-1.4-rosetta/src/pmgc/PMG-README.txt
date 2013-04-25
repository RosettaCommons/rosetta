c* ///////////////////////////////////////////////////////////////////////////
c*                           Welcome to PMG (MG/XMG).
c* ///////////////////////////////////////////////////////////////////////////
c*
c* ======================
c* The Author Information
c* ======================
c*
c* The computer codes PMG (MG/XMG, "the Code"), including FORTRAN/C/C++
c* language multilevel subroutines called "MG" and the C-language X-Window
c* interface called "XMG", were developed by:
c*
c*    Michael Holst                TELE:  (858) 534-4899
c*    Department of Mathematics    FAX:   (858) 534-5273
c*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
c*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
c*
c* ====================================================
c* Licensing, permission to use, modify, and distribute
c* ====================================================
c*
c* The code has been placed under the GNU General Public License (GPL).  
c*
c* This means that essentially you can use and modify the code as you like, 
c* under a few conditions, such as acknowledging the authors of the original 
c* code in any derived works, such as new programs or research results using 
c* the code.
c*
c* Please have a look at the GPL which should have accompanied the code.
c* If it is missing, please write to the Free Software Foundation, Inc., 
c* 675 Mass Ave, Cambridge, MA 02139, USA.
c*
c* It would have been easier for me to just keep the code completely private, 
c* or to hand out very restricted use copies.  However, in the interests of 
c* freedom of information, the advancement of science, the future of the 
c* planet, and so on, I decided to distribute the code in this very free way.
c* The GNU world functions very much on the honor system.  I want you to have 
c* all of my source code because I think it will do more toward helping you 
c* understand the ideas in there than anything else.  The GNU idea is that in 
c* return for all of my work on the source code and the ideas behind it, all I 
c* ask for is some acknowledgment in your derived works for using my original 
c* code and the original ideas that went into the code.  For GNU to be viable, 
c* we all need to abide by this very minimal requirement.
c*
c* I'm typing this on a nearly all-GNU system right now, namely Linux 1.2.8, 
c* a GNU version of UN*X available under the same conditions as PMG.
c* I.e., you get all of the source code for the entire UN*X operating system.
c* For more information about Linux, see: http://sunsite.unc.edu/mdw/linux.html
c*
c* ====================================================================
c* Pointers to Papers and other Reference Materials Related to the Code
c* ====================================================================
c*
c* The methods and techniques  implemented  in  the  Code  were  developed  by 
c* Michael Holst as a part of his PhD research, and are described in detail in 
c* the PhD thesis:
c*
c*    @phdthesis{Hols93a,
c*       author    = "M. Holst",
c*       title     = "Multilevel Methods for the {Poisson-Boltzmann} Equation",
c*       note      = "Also published as Tech. Rep. UIUCDCS-R-03-1821",
c*       school    = "Numerical Computing Group,
c*                    Department of Computer Science,
c*                    University of Illinois at Urbana-Champaign",
c*       year      = 1993 }
c*
c* Also refer to the supporting papers:
c*
c*    @manual{Hols94d,
c*       author    = "M. Holst",
c*       title     = "The {Poisson-Boltzmann} Equation:
c*                    Analysis and Multilevel Numerical Solution",
c*       note      = "Monograph (updated and extended form of the
c*                    Ph.D. thesis~\cite{Hols93a})." }
c*    
c*    @article{SHHN94,
c*       author    = "R. Sampogna and J. Hecht and M. Holst and A. Nicholls
c*                    and B. Honig",
c*       title     = "Nonlinear {Poisson-Boltzmann} Calculation of {pKa} 
c*                    Values",
c*       journal   = "Biophysics",
c*       note      = "(In Progress?)",
c*       year      = 1994 }
c*    
c*    @article{HoVa94b,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
c*       journal   = SINUM,
c*       note      = "(Accepted)",
c*       year      = 1994 }
c*    
c*    @article{Hols94e,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Numerical Solution of the Nonlinear {Poisson-Boltzmann}
c*                    Equation: Developing More Robust and Efficient Methods",
c*       journal   = JCC,
c*       volume    = "16",
c*       number    = "3",
c*       pages     = "337--364",
c*       year      = 1995 }
c*    
c*    @article{HoSa93a,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Multigrid Solution of the {Poisson-Boltzmann} Equation",
c*       journal   = JCC,
c*       volume    = "14",
c*       number    = "1",
c*       pages     = "105--113",
c*       year      = 1993 }
c*    
c*    @article{HKSS94,
c*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
c*       title     = "Protein Electrostatics: Rapid Multigrid-based {Newton}
c*                    Algorithm for Solution of the Full Nonlinear
c*                    {Poisson-Boltzmann} Equation",
c*       journal   = "J. Biomol. Struct. Dyn.",
c*       volume    = "11",
c*       pages     = "1437--1445",
c*       year      = 1994 }
c*    
c*    @article{HKSS93b,
c*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
c*       title     = "Treatment of Electrostatic Effects in Proteins:
c*                    Multigrid-based-{Newton} Iterative Method for Solution 
c*                    of the Full Nonlinear {Poisson-Boltzmann} Equation",
c*       journal   = "Proteins: Structure, Function, and Genetics",
c*       volume    = "18",
c*       number    = "3",
c*       pages     = "231--245",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94h,
c*       author    = "M. Holst",
c*       title     = "A Theoretical Analysis of the {Poisson-Boltzmann} 
c*                    Equation: A priori Estimates and Well-posedness",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       year      = 1994 }
c*    
c*    @inproceedings{HoVa95a,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: To Symmetrize or not to Symmetrize",
c*       booktitle = "Proceedings of the Seventh Copper Mountain Conference on
c*                    Multigrid Methods, April 2-7, 1995, Copper Mountain, 
c*                    Colorado",
c*       editor    = "J. Mandel and S. McCormick",
c*       publisher = "NASA Langley Research Center",
c*       year      = "1995" }
c*    
c*    @inproceedings{HoSa93b,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Multigrid and Domain Decomposition Methods for 
c*                    Electrostatics Problems",
c*       booktitle = "Domain Decomposition Methods in Science and Engineering
c*                    (Proceedings of the Seventh International Conference on
c*                    Domain Decomposition, October 27-30, 1993, 
c*                    The Pennsylvania State University)",
c*       editor    = "D. E. Keyes and J. Xu",
c*       publisher = "American Mathematical Society, Providence",
c*       year      = "1995" }
c*    
c*    @techreport{HoVa94a,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-13",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94c,
c*       author    = "M. Holst",
c*       title     = "{An} {Algebraic} {Schwarz} {Theory}",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-12",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94f,
c*       author    = "M. Holst",
c*       title     = "A Robust and Efficient Numerical Method for Nonlinear
c*                    Protein Modeling Equations",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-9",
c*       year      = 1994 }
c*    
c*    @techreport{HoSa93c,
c*       author    = "M. Holst and F. Saied",
c*       title     = "A Short Note Comparing Multigrid and Domain Decomposition
c*                    for Protein Modeling Equations",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-10",
c*       year      = 1994 }
c*
c* =====================
c* PMG Bugs and Notes
c* =====================
c*
c* The  timing  routine,  found  in  "secd.c", is machine dependent.  We  have 
c* implemented  a  number  of  routines  for  different  machines, including a
c* generic  routine  for  a  standard  unix  machine with the getrusage system
c* routine.
c*
c* =============================
c* PMG Last Modification Date
c* =============================
c*
c* MG:  10-01-95
c* XMG: 10-01-95
c*
c* ///////////////////////////////////////////////////////////////////////////
