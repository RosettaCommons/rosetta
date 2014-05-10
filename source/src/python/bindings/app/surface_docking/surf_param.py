# Import python modules
from string import *
import math

# Import Rosetta modules
from rosetta import *
from rosetta.protocols.rigid import *
import rosetta.core.scoring.solid_surface

#============defining objects=====================
#Vector objects in Rosetta
BodyPosition = rosetta.numeric
#for vector, use as body.xyzVector_double

#============defining functions====================
def load(file):
    file = open(file, 'r').readlines()
    loaded_file = []
    for line in file:
      loaded_file.append(strip(line))
    return loaded_file

def get_markers(pdb):
    n=0
    for atom in pdb:
        if pdb[n][:6]=='SURFA0':
            Ax = float(pdb[n][30:38])
            Ay = float(pdb[n][38:46])
            Az = float(pdb[n][46:54])
            SURFA0 = BodyPosition.xyzVector_double(Ax,Ay,Az)
        if pdb[n][:6]=='SURFA1':
            Bx = float(pdb[n][30:38])
            By = float(pdb[n][38:46])
            Bz = float(pdb[n][46:54])
            SURFA1 = BodyPosition.xyzVector_double(Bx,By,Bz)            
        if pdb[n][:6]=='SURFA2':
            Cx = float(pdb[n][30:38])
            Cy = float(pdb[n][38:46])
            Cz = float(pdb[n][46:54])
            SURFA2 = BodyPosition.xyzVector_double(Cx,Cy,Cz)          
        n = n+1
    return (SURFA0,SURFA1, SURFA2)


           
#==============protein_surfacemove_symmetry==============
'''
moving protein to the projected center of geometry of surface by symmetry
'''
#center of geometry of protein P
#*****need to define pose and prot_length

def centerG_P(pose,prot_length):
    center_of_geo_prot = BodyPosition.xyzVector_double(0,0,0)
    for i in range(1, pose.total_residue()+1):
        if (pose.residue(i).is_protein() == True):
            ca_pose = pose.residue(i).atom('CA').xyz()
            center_of_geo_prot = center_of_geo_prot + ca_pose
    center_of_geo_prot = (center_of_geo_prot/prot_length)
    return center_of_geo_prot

#center of geometry of surface S
def centerG_S(pose,prot_length):
    center_of_geo_surf = BodyPosition.xyzVector_double(0,0,0)
    for i in range(1, pose.total_residue()+1):
        if (pose.residue(i).is_ligand() == True):
            surf_pose = pose.residue(i).xyz(pose.residue(i).nbr_atom())
            center_of_geo_surf = center_of_geo_surf + surf_pose
    center_of_geo_surf = (center_of_geo_surf/(pose.total_residue() - prot_length))
    return center_of_geo_surf

#anti parallel vector
def anti_parallel_vector(vector):
    anti_parallel_vector = vector * -1
    return anti_parallel_vector

#...................Built in rosetta vector stuff.............
def length(vector):
    length = vector.norm
    return length

def unit_vector(vector):
    unit_vector = vector.normalize()
    return unit_vector

def dot_product(vector1, vector2):
    dot_product = vector1.dot_product(vector2)
    return dot_product

def inner_product(vector1, vector2):
    inner_product = vector1.inner_product(vector2)
    return inner_product

def cross_product(vector1, vector2):
    cross_product = vector1.cross_product(vector2)
    return cross_product

def normal_vector(vector1, vector2):
    normal_vector = cross_product(vector1, vector2)
    return normal_vector

def normalto3points(pointA, pointB, pointC):
    vAB = pointB - pointA
    vAC = pointC - pointA
    normal = normal_vector(vAB, vAC)
    return normal

#calculate a point from 2 intersecting lines
#--------------------------------
#        |
#   ---B-*----------
#        |
#        C 
#        |

def point2lines2(point1, point2, vectorB, vectorC):
    CB = point2 - point1
    crossBC = cross_product(vectorB, vectorC)
    CB_C = cross_product(vectorC, crossBC)
    DD = dot_product(crossBC,crossBC)
    t = dot_product(CB, CB_C)/DD
    intersection = point1 + vectorB * t
    return intersection

#calculate a point from 2 intersected lines
def point2lines3(point1, point2, vectorB, vectorC):
    CB = point2 - point1
    cross1C = cross_product(CB, vectorC)
    crossBC = cross_product(vectorB, vectorC)
    w = dot_product(cross1C, crossBC)
    DD = inner_product(crossBC, crossBC)
    t = w/DD
    intersection = point1 + vectorB * t
    return intersection

#calculate a plane from 3 points; returns A, B, C, D in A(x1-x)+B(y1-y)+C(z1-z)+D=0
def plane3points(pointA, pointB, pointC):
    normal_to_3points = normalto3points(pointA, pointB, pointC)
    Aplane = normal_to_3points.x
    Bplane = normal_to_3points.y
    Cplane = normal_to_3points.z
    Dplane = -(Aplane*pointA.x + Bplane*pointA.y + Cplane*pointA.z) 
    return (Aplane, Bplane, Cplane, Dplane)

#calculate a point in a plane hit from a point out of plane normal
def point_intersect_plane(plane_abcd, point_out_of_plane, antiparallel_normal_vector):
    t = -(plane_abcd[0] * point_out_of_plane.x + plane_abcd[1]* point_out_of_plane.y \
	 + plane_abcd[2]*point_out_of_plane.z + plane_abcd[3])/(plane_abcd[0]*\
	plane_abcd[0] + plane_abcd[1]*plane_abcd[1] + plane_abcd[2]*plane_abcd[2])
    point_of_intersection = point_out_of_plane - t * antiparallel_normal_vector
    return point_of_intersection
    

#print info:
#d_Psurfport_cG = length(Psurface_Nprotein, centerG_surface)
def print_surf_move(tDistanceB, tDistanceC, Psurface_Nprotein, centerG_surface, CA_centerG_protein, plane0, findProteinZ):
    print "tDistanceB = "+str(tDistanceB)+"  tDistanceC= "+str(tDistanceC)
    print "distance protein is away from the surface: "+str(length(findProteinZ))     
    print "centerG_surface:"
    print centerG_surface
    print "CA_centerG_protein:"
    print CA_centerG_protein
    print "the equation plane of surface"
    print "with with pt= the intersection point of the antiparallel normal vector and the surface plane:"
    print "("+str(plane0[0])+")*("+str(Psurface_Nprotein.x)+"-x)+("+str(plane0[1])+")*("+str(Psurface_Nprotein.y)+"-y)+("+str(plane0[2])+")*("+str(Psurface_Nprotein.z)+"-z)+("+str(plane0[3])+")=0"
    print "distance protein is away from the surface: "+str(length(findProteinZ))
    
#Translating                                                
def Translate_to_center(pose, start, Translate):
    dock_jump=start-1
    Prot_trans=RigidBodyTransMover(pose,dock_jump)
    Prot_trans.trans_axis(Translate)
    Prot_trans.step_size(length(Translate))
    Prot_trans.apply(pose)

def Translate_up(pose,start,axis,mag):
    dock_jump=start-1
    Prot_trans=RigidBodyTransMover(pose,dock_jump) 
    Prot_trans.trans_axis(axis)
    Prot_trans.step_size(mag)
    Prot_trans.apply(pose)
  

#=================vector stuff; the main part=====================
def surface_symm_move(pose, start, end, SURFA0, SURFA1, SURFA2):

    #Define protein length
    prot_length = end-start+1
    
    #center of geometry of protein
    CA_centerG_protein = centerG_P(pose,prot_length)

    #center of geometry of surface
    centerG_surface = centerG_S(pose,prot_length)

    #initial ABC points
    surface3pointA0 = SURFA0
    surface3pointB0 = SURFA1
    surface3pointC0 = SURFA2
     
    #two basis vectors of the surface crystal system; unit cell dimension
    vAB0 = surface3pointB0 - surface3pointA0
    vAC0 = surface3pointC0 - surface3pointA0

    #fix the ABC in the center of surface
    surface3pointA = centerG_surface
    surface3pointB = centerG_surface + vAB0
    surface3pointC = centerG_surface + vAC0
    
    
    vAB = surface3pointB - surface3pointA
    vAC = surface3pointC - surface3pointA

    #find the normal vector the surface of the vectors AB, AC
    Nprotein_to_surface = normalto3points(surface3pointA, surface3pointB, surface3pointC)
    
    #vector antiparallel
    iNprotein_to_surface = anti_parallel_vector(Nprotein_to_surface)

    #find the equation of the plane
    plane0 =plane3points(surface3pointA, surface3pointB, surface3pointC)

    #find the intersection point of the antiparallel normal vector and the surface plane
    Psurface_Nprotein = point_intersect_plane(plane0, CA_centerG_protein, iNprotein_to_surface)

    #find the vector difference of the noraml intersection point and the center of geometry of the protein
    findProteinZ = Psurface_Nprotein - CA_centerG_protein

    #find the 2 first intersection points
    IntersectPoint1 = point2lines3(centerG_surface,Psurface_Nprotein, vAC, vAB)
    IntersectPoint2 = point2lines3(centerG_surface,Psurface_Nprotein, vAB, vAC)
    
    #check = normalto3points(centerG_surface, IntersectPoint1, IntersectPoint2)
    #print "Compare with Nprotein_to_surface:"
    #print check.normalize()
    
    #find the distances from the center of surface parallel AB, AC
    tCvector = IntersectPoint1 - Psurface_Nprotein
    tDistanceC = length(tCvector)
    tBvector = IntersectPoint2 - Psurface_Nprotein
    tDistanceB = length(tBvector)

    #begin the criteria to move or not the surface to center the protein
    ABdistance = length(vAB0)
    ACdistance = length(vAC0)

    """
    #print info: initial
    """
    #print "........................coordinates..........................."
    #print_surf_move(tDistanceB, tDistanceC, Psurface_Nprotein, centerG_surface, CA_centerG_protein, plane0, findProteinZ)
      
    """
    translation conditions
    """
    
    #Initialize translation vector
    Translate = BodyPosition.xyzVector_double(0.0,0.0,0.0)
    
    if (tDistanceC <= ACdistance and tDistanceB <= ABdistance):
        New_centerG_protein = CA_centerG_protein
        #print "No move!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        
    else:
        
        """
        Begin move criteria
        """
        
        if (tDistanceC > ACdistance):
            #print "Yes move C!!!!!!!!!!!!!!!!!!!!"
            cc = tDistanceC / ACdistance
            icc = int(cc)
            floor_ceil_icc = cc - icc
            IntersectPoint1 = (cc-icc)*(IntersectPoint1 - centerG_surface)/cc + centerG_surface;
           
        if (tDistanceB > ABdistance):
            #print "Yes move B!!!!!!!!!!!!!!!!!!!!"
            bb = tDistanceB / ABdistance
            ibb = int(bb)
            floor_ceil_ibb = bb - ibb
            IntersectPoint2 = (bb-ibb)*(IntersectPoint2 - centerG_surface)/bb + centerG_surface;

        fvAB = IntersectPoint2 - centerG_surface
        fvAC = IntersectPoint1 - centerG_surface
        tDistanceC = length( fvAC ) 
        tDistanceB = length( fvAB ) 

        #Find the final intersection point
        FinalPointPlane1 = point2lines3  ( IntersectPoint1, IntersectPoint2, fvAB, fvAC );
        New_centerG_protein = FinalPointPlane1 - findProteinZ;
        
        #Calculate new translationvector 
        Translate = New_centerG_protein - CA_centerG_protein
        Translate_to_center(pose, start, Translate)
        
        #print "Check final centerG of protein:"
        #final_centerG_protein = centerG_P(pose,prot_length) 
        #print final_centerG_protein
        #Trnaslate protein up
        #if length(final_centerG_protein-centerG_surface) < 20:
        #	print "translated up"
        #	Translate_up(pose,start,Nprotein_to_surface,80)
        
        
    """
    #print info: final
    """
    #print "..................Translated coordinates......................"
    #print "Translation x= "+str(Translate.x)+"  y= "+str(Translate.y)+"  z= "+str(Translate.z)
    #print "\n"


