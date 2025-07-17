###########################################################################################
#
#   Script creating inputs for the SMD of inactivation particle pulling out of its binding pocket (recovery) 
#  
#   Series of Steered MD simulations are performed each in different direction of pulling 
#   to roughly cover a semisphere. Each simulation data is stored in separate directory. 
#   To run it, you need to have numpy, scipy and MDAnalysis libraries installed.
#
#   Author: Katarzyna Walczewska-Szewc, Nicolaus Copernicus University in Toruń
###########################################################################################

from MDAnalysis import *
from numpy import *
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import os
import sys
from zipfile import ZipFile
##########################################################################################
#
#   The input data and parameters
#
##########################################################################################

dirname = 'smd_inact_ligname' # name of the new directory
input_pdb = str(sys.argv[1]) # step5_input.pdb
input_psf = str(sys.argv[2]) # step5_input.psf
input_vel = str(sys.argv[3]) # step7_production.vel
input_coor =str(sys.argv[4]) # step7_production.coor
input_xsc = str(sys.argv[5]) # step7_production.xsc
input_par1 = str(sys.argv[6]) # toppar.zip
# Templates for namd input file and run.bash file
template_inp = str(sys.argv[7]) # template.inp
template_run = str(sys.argv[8]) # template.sh

#Selection parameters for SMD (see MDAnalysis page for more information about selection keywords)
input_sel1= sys.argv[9] # 'protein and name CA and resid 387 972 1498 1791'   # constrained atoms in SMD; here: selectivity filter
input_sel2= sys.argv[10] # 'protein and name CA and resid 1565:1567'  # pulled atoms in SMD; here: MFM motif

###########################################################################################
#
#   End of parameters setup, do not change unless you know what you are doing ;)
#
###########################################################################################
# Function generating the set of vectors samplig hemisphere. 
# Such hemisphere is alighed to the direction of the vector joining COM of the fix and pull selections from the pdb file.

def Hedgehog(pdb,fix_selection,pull_selection):

    u = Universe(pdb)
    fix = u.select_atoms(fix_selection)
    pull = u.select_atoms(pull_selection)
    fix_COM = fix.center_of_mass()
    pull_COM = pull.center_of_mass()
    pocket = u.select_atoms("protein and name CA and resid 1411 1552 1555 1739 1740 1842 1850") # selection for the IFM/MFM binding pocket of a given Nav 
    pocket_COM = pocket.center_of_mass()
    ax_principal = pull_COM - pocket_COM #creating vector pointing the general direction of the cone
    ax_principal = ax_principal/linalg.norm(ax_principal)#creating the unit vector
    print((ax_principal))

    #Generating a hedgehog of vectors in z- direction
    hedgehog = array([[0,0,1]])
    labels = array([[0, 0]])
    for theta,resolution in zip([30, 60],[90, 90]):#[theta angles-you can put more values if you want],[resolution in phi anglesfor each theta- 90 means four directions for each theta (360/90 = 4)]):
        for phi in range(0,360,resolution):
            x = cos(deg2rad(phi))*sin(deg2rad(theta))
            y = sin(deg2rad(phi))*sin(deg2rad(theta))
            z = cos(deg2rad(theta))

            hedgehog = concatenate((hedgehog,array([[x,y,z]])),axis=0)
            labels = concatenate((labels,array([[theta,phi]])),axis=0) #labels for subdirectories names construction

    #Transforming the vectors cone to the direction given by the principal ax:
    r = R.align_vectors([ax_principal],[[0,0,1]])
    hedgehog = r[0].apply(hedgehog)
    #translation of the anchor of the cone to the center of mass od the pulled protein is not neccessary because 
    # NAMD is doing it automatically.

    return hedgehog, labels, pull_COM #set of vectors [N][x,y,z]; set of angles [N][theta,phi], the anchor point [x,y,z]

##############################################################################################
#
#   Function generating a tree of directories containing input and run  files for each SMD run
#   The template.inp and run.bash files are modified and placed in proper directories. The copies of the NAMD restart files
#   are created in the output directory, to make it self-sustainable.
#
##############################################################################################

def Gen_input(name,pdb,psf,vectors,template,sel1,sel2,par1,coor=None,vel=None,xsc=None):
    
    #   name/ - the output directory containing all generated files, copy of the NAMD input structures and the simulation results
    #       SMD_constraints.pdb - the output file containing information about fixed and pulled atoms for SMD (O and B colum, respectively)
    #       SMD_theta_i_phi_j/ - subdirectories for each SMD direction run (the simulation output will be stored here)
    #           mdrun.inp - input file for a single SMD simulation
    #           run.bash - the run.bash script (generated in next function)
    #       vmd_script.tcl - the vmd script to visualize the cone of vectors (generated below)
    
    # Creating copy of NAMD restart files in the output directory (coor,vel and xsc filea are not mandatory 
    # but without them it is necessary to adjust your template.inp files to initialize temperatures and create proper PBC box)
    os.system('mkdir '+str(name))
    os.system('cp '+str(pdb)+' '+str(name)+'/')
    os.system('cp '+str(psf)+' '+str(name)+'/')
    if str(par1)[-10:] =='toppar.zip': 
        with ZipFile(par1,'r') as z:
            z.extractall(str(name)+'/')
    if str(par1)[-10:] != 'toppar.zip': os.system('cp -r '+str(par1)+' '+str(name)+'/')
    if coor: os.system('cp '+str(coor)+' '+str(name)+'/')
    if vel: os.system('cp '+str(vel)+' '+str(name)+'/')
    if xsc: os.system('cp '+str(xsc)+' '+str(name)+'/')

    # Creating SMD_constraints.pdb file with flags in O and B column indicating fixed and pulled atoms
    u = Universe(pdb)
    all = u.select_atoms('all')
    all.occupancies = [0 for i in range(len(all.occupancies))]
    all.tempfactors = [0 for i in range(len(all.occupancies))]
    fix = u.select_atoms(sel1)
    fix.tempfactors = [1 for i in range(len(fix.tempfactors))]
    pul = u.select_atoms(sel2)
    pul.occupancies = [1 for i in range(len(pul.occupancies))]
    all.write(str(name)+'/SMD_constraints.pdb')

    # Creating subdirectories for all directions in the semisphere
    for v,l in zip(vectors[0],vectors[1]):
        os.system('mkdir '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1])))

        # Modyfying the template.inp file for each SMD run (the only difference between them is the direction of pulling)
        f = open(template,'r')
        new_f= open(str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/mdrun.inp','w')

        # changing information about the structure and parameters
        new_f.write('structure	../'+str(psf.split('/')[-1])+'\n')
        new_f.write('coordinates	../'+str(pdb.split('/')[-1])+'\n')
        if coor: new_f.write('bincoordinates	../'+str(coor.split('/')[-1])+'\n')
        if vel: new_f.write('binvelocities	../'+str(vel.split('/')[-1])+'\n')
        if xsc: new_f.write('extendedSystem	../'+str(xsc.split('/')[-1])+'\n')
        new_f.write('paratypecharmm	on\n')
        if par1[-10:] =='toppar.zip':
            for filename in os.listdir(str(name)+'/toppar'):
                new_f.write('parameters	../toppar/'+str(filename)+'\n')
        if par1[-10:]!='toppar.zip':
            new_f.write('parameters ../'+str(par1.split('/')[-1])+'\n')

        for line in f.readlines():
            if line[0:7].lower() != 'structu' and line[0:7].lower() != 'coordin' and line[0:7].lower() != 'bincoor' and line[0:7].lower() != 'binvelo' and line[0:7].lower() != 'extende' and line[0:7].lower() != 'paramet' and line[0:7].lower() != 'paratyp':
              line_new = line
              if line[0:10] == 'outputname': line_new = 'outputname      md'+'\n'
              if line[0:7] == 'SMDfile': line_new = 'SMDfile ../SMD_constraints.pdb'+'\n' #The files common for each SMD_run are stored in the output directory
              if line[0:6] == 'SMDDir': line_new = 'SMDDir	'+str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n'
              if line[0:7] == 'consref': line_new = 'consref ../SMD_constraints.pdb'+'\n'
              if line[0:10] == 'conskfile ': line_new = 'conskfile ../SMD_constraints.pdb'+'\n'
              new_f.write(line_new)
        new_f.close()
        f.close()

##############################################################################################
#
#   Function generating run.bash files inside the tree of subdirectories for each SMD run
#   This function also initiate the simulation so if you do not want it to run immediately,
#   please comment out the proper line. The template.bash file should be modified to 
#   fit your computing facility.
#
##############################################################################################

def Gen_run(name,label,template):
    master_file = open(str(name)+'/master.run','w')
    for l in label:
        # Generating run.bash files
        f = open(template,'r')
        new_f= open(str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/run.bash','w')
        for line in f.readlines():
            line_new = line
            if len(line.split())>1 and line.split()[1] == '-J': line_new = line.split()[0] + ' ' + line.split()[1] +' SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'\n'
            
            if len(line.split())>1 and line.split()[1][0:4] == 'INPF': line_new = 'set INPF=SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/mdrun.inp'+'\n'
            if line[0:4] == 'INPF': line_new = 'INPF=mdrun.inp'+'\n'
            if line[0:4] == 'OUTF': line_new = 'OUTF=mdrun.log'+'\n'

            if len(line.split())>1 and line.split()[1][0:4] == 'OUTF': line_new = 'set OUTF=mdrun.log'+'\n'
            new_f.write(line_new)
        #print('qsub '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/run.bash')
        new_f.close()
        f.close()

        # Runing the simulations
        os.system('chmod +x '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/run.bash')
        #print('Running simulation for SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1])))
        #os.system('. '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/run.bash &') # If you are using this script on a supercomputer, all jobs can run simultanously 
        master_file.write('. SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'/run.bash\n') # On a single gpu station, it could be better to run one job after one (the line without & at the end so the next job is waiting for the one to finish)
    master_file.close()
    os.system('chmod +x '+str(name)+'/master.run')
        
#########################################################################################
#
#   The tcl script that can be run in the VMD tcl console to visualize the cone of vectors
#
#########################################################################################
def Gen_vmd_script(name,vector,COM):
    f = open(name+'/vmd.tcl','w')
    #array draving subroutine (from https://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html)
    f.write('proc vmd_draw_arrow {mol start end} {\n\t# an arrow is made of a cylinder and a cone\n\t set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]\n\t graphics $mol cylinder $start $middle radius 0.15\n\t graphics $mol cone $middle $end radius 0.25\n}\n')
    s=5 #scaling the vectors
    for j in vector:
        i = j + COM # The vectors need to be translated to the center of the mass of the protein and rescaled by the scaling factor (below)
        f.write('draw arrow {'+str(COM[0])+' '+str(COM[1])+' '+str(COM[2])+'} {'+str(i[0]+(i[0]-COM[0])*s)+' '+str(i[1]+(i[1]-COM[1])*s)+' '+str(i[2]+(i[2]-COM[2])*s)+'}\n' )

#########################################################################################
#
#   The body of the script
#
#########################################################################################

hedgehog = Hedgehog(input_pdb, input_sel1, input_sel2)
Gen_input(dirname, input_pdb, input_psf, hedgehog, template_inp, input_sel1, input_sel2, input_par1, vel=input_vel, coor=input_coor, xsc=input_xsc )
Gen_run(dirname, hedgehog[1], template_run)
Gen_vmd_script(dirname, hedgehog[0], hedgehog[2])
print('Generation of SMD input finished!')
