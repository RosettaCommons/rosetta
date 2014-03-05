/// @namespace protocols::loop_modeling::samplers
///
/// @brief Algorithms for generating new loop conformations.
/// 
/// @details The classes in this namespace are broadly supposed to be 
/// responsible for generating new loop conformations.  This is in contrast to 
/// the classes in the refiners namespace, which are responsible for lowering 
/// the score of a new loop conformation.  The primary loop sampling algorithms 
/// are backrub, kinematic closure (KIC), and cyclic-coordinate descent (CCD).  
/// These algorithms are all related by the fact that they generate new 
/// backbone conformations within a given window, and do not affect any DOFs 
/// outside that window.   Backrub is distinguished from KIC and CCD by the 
/// fact that it does not allow any DOFs to be explicitly set.  KIC is much 
/// faster than CCD, but does not yet support fragment insertion.
