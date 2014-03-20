Rosetta Membrane Refactor
-----
Version 1.3.1 - Added on 3/20/2014

Author: Rebecca Alford (rfalford12@gmail.com)

Main Directory - Contians Components of static and dynamic representations of membrane proteins in Rosetta. Quick Descriptions of directory structure

 - Top Level: Contains Membrane Protein Factory Code - Responsibile
   for full initialization of poses as membrane proteins. Loads as a 
   top level resource with tag "startstruct" to work with JD2
 - Config: Register protocols with mebrane protocol config manager
   to update and maintain core support
 - Geometry: Code for representing geometry of the membrane and 
   respective membrane embedding for protein chains
 - Scoring: Refactored scoring methods with modified interface for pre-build
   search and score methods
 - io: Resource manager classes respondible for loading membrane related data
 - kinematics: Membrane protein specific fold tree code
 - util: Membrane definitions and exceptions

Last Modified: 12/5/13
