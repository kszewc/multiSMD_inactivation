from numpy import *
from MDAnalysis import *
import os,sys
import matplotlib.pyplot as plt
import matplotlib
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def Dist(a, b):
    r = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return r


name = str(sys.argv[1]) # Name of your SMD output directory
sel_const = str(sys.argv[2]) # First selection - constrained atoms
sel_pull = str(sys.argv[3]) # Second selection - pulled atoms
n_repeats = int(sys.argv[4]) # Here is the number of simulation repeats for each pulling direction
extract_forces = True # it is very time-consuming, change it to False after force files are generated
calculate_HB = True # it is very time-consuming, change it to False after HB files are generated
red = 5 # reducing factor for too dense data

list_of_dirs = [subdir for subdir in os.listdir(name) if subdir[0:19] == 'SMD_theta_60_phi_90'] # change the direction you want to analyze

## this part takes some time, since it is analysing your SMD trajectories looking for possible H-bonds between selections
u = Universe(name+'/SMD_constraints.pdb')
fix = u.select_atoms(sel_const).center_of_mass() #
bunch = [] #here the information about pulling vectors will be stored

for dir in list_of_dirs:
    print(dir," calculations started.")
    for n in range(0,n_repeats):
        f = open(name+'/'+dir+f'/mdrun{n}.log','r')

        # SMD force vs time and SMD force vs distance
        if extract_forces:

            ft = open(name + '/' + dir + f'/smd_force_time{n}.dat', 'w')
            ft.write('# time [fs] vs force of pulling\n')
            fd = open(name + '/' + dir + f'/smd_force_dist{n}.dat', 'w')
            fd.write('# distance between COMs vs force of pulling\n')

            for line in f.readlines():
                if len(line)>5:
                    if line[:5] == 'SMD  ':
                        t = line.split()[1] #timestep
                        r = array([float(line.split()[2]),float(line.split()[3]),float(line.split()[4])]) # COM position of pulled domain
                        f = array([float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]) # Pulling force, starting from COM
                        ft.write(t+' '+str(linalg.norm(f))+'\n') #((f[0]-r[0])**2 + (f[1]-r[1])**2 + (f[2]-r[2])**2)**0.5)+'\n')
                        fd.write(str(Dist(r,fix))+' '+str(linalg.norm(f))+'\n') #((f[0]-r[0])**2 + (f[1]-r[1])**2 + (f[2]-r[2])**2)**0.5)+'\n')
            ft.close()
            fd.close()

        if calculate_HB:
            ht = open(name + '/' + dir +f'/smd_hb_time{n}.dat', 'w')
            ht.write('# time [fs] vs number of hydrogen bonds between domains\n')
            # Hydrogen bonds number vs time
            psf = [file for file in os.listdir(name) if file[-3:] == 'psf']
            print('Hydrogen bonds calculating - that may take some time...')
            u = Universe(name+'/'+psf[0] ,name+'/'+dir+f'/md{n}.dcd')
            hbonds = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel='protein and resid 1566',acceptors_sel='protein and resid 1411 1552 1555 1739 1740 1842 1850') # donor and acceptor residues for Hbond analysis
            hbonds.run()
            out = hbonds.count_by_time()
            hbonds_rev = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel='protein and resid 1411 1552 1555 1739 1740 1842 1850',acceptors_sel='protein and resid 1566') # donor and acceptor residues for Hbond analysis
            hbonds_rev.run()
            out2 = hbonds_rev.count_by_time()
            t=0
            for i in (out+out2):
                ht.write(str(t)+' '+str(i)+'\n')
                t+= 0.05 #[ns]
            ht.close()

            df=open(name+'/'+dir+f'/smd_distance_force{n}.dat', 'w') # distance based on pocket residues
            df.write('# MFM-pocket distance vs force\n')
            plik_z_silami=loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')[::50,1] #### adjust sampling of log files and dcd frames
            print(len(u.trajectory))
            print(len(plik_z_silami))
                
            print('Distance vs force calculation - that may take some time...')
            for t, force in zip( u.trajectory[:], plik_z_silami) : 
                pocket=u.select_atoms('protein and name CA and resid 1411 1552 1739 1850').center_of_mass() ##### binding pocket
                mfm=u.select_atoms('protein and name CA and resid 1565:1567').center_of_mass() ### MFM
                distance=(Dist(pocket,mfm))
            
                df.write(str(distance)+' '+str(force)+'\n')
            df.close()
    
    #Writing down used pulling vectors
    theta = double(dir.split('_')[2])
    phi = double(dir.split('_')[4])
    x = cos(deg2rad(phi))*sin(deg2rad(theta))
    y = sin(deg2rad(phi))*sin(deg2rad(theta))
    z = cos(deg2rad(theta))
    bunch.append([x,y,z])
savetxt(name+'/bunch_of_vectors.dat',bunch)
    

# Plots - this script is adjusted for 1 pulling direction
# Multiple repeats results are plotted as average value + outline based on calculated standard deviation value.

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(10,3))

sorted=['SMD_theta_60_phi_90'
]

# looking for a maximum force value and HB numbers in all files to rescale figures
max_force = max( [max([max(loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')[:,1]) for n in range(0,n_repeats)]) for dir in sorted[:]] )
max_HB = max( [max([max(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:,1]) for n in range(0,n_repeats)]) for dir in sorted[:]] )
print(max_force, max_HB)
for dir in sorted[:]:
    print(dir)

    for n in range(0,n_repeats):

        if n == 0:
            ft = reshape(loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')[:,0],(len(loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')),1))
            ht = reshape(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:,0],(len(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')),1))
            fd = reshape(loadtxt(name+'/'+dir+f'/smd_force_dist{n}.dat')[:,0],(len(loadtxt(name+'/'+dir+f'/smd_force_dist{n}.dat')),1))
            df = reshape(loadtxt(name+'/'+dir+f'/smd_distance_force{n}.dat')[:,0],(len(loadtxt(name+'/'+dir+f'/smd_distance_force{n}.dat')),1))
            print(shape(ft),shape(ht),shape(fd),shape(df))
        #if n > 0:
        ft = concatenate((ft,reshape(loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')[:,1],(len(loadtxt(name+'/'+dir+f'/smd_force_time{n}.dat')),1))),axis=1)
        ht = concatenate((ht, reshape(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:,1],(len(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')),1))),axis=1)
        fd = concatenate((fd, reshape(loadtxt(name+'/'+dir+f'/smd_force_dist{n}.dat')[:,1],(len(loadtxt(name+'/'+dir+f'/smd_force_dist{n}.dat')),1))),axis=1)
        df = concatenate((df, reshape(loadtxt(name+'/'+dir+f'/smd_distance_force{n}.dat')[:,1],(len(loadtxt(name+'/'+dir+f'/smd_distance_force{n}.dat')),1))),axis=1)
        print(shape(ft),shape(ht),shape(fd), shape(df))


    #Force vs time
    ax1.fill_between(ft[::red, 0] / 1000000, mean(ft[::red,1:],axis=1) + std(ft[::red,1:],axis=1), mean(ft[::red,1:],axis=1) - std(ft[::red,1:],axis=1), linewidth=1,alpha=0.5, color='orange')
    ax1.plot(ft[::red,0]/1000000,   mean(ft[::red,1:],axis=1)    ,linewidth=2, color='orange')
    ax1.set_xlim(0,)
    ax1.set_ylim(0, max_force)
    ax1.set_xlabel('Time [ns]').set_fontsize(16)
    ax1.set_title('SMD forces')


    #Distance vs force
    ax2.fill_between(df[:, 0], mean(df[:, 1:], axis=1) + std(df[:, 1:], axis=1),
                          mean(df[:, 1:], axis=1) - std(df[:, 1:], axis=1), linewidth=1, alpha=0.5, color='orange')

# OR :
#    ax2.scatter(df[:, 0], mean(df[:, 1:], axis=1) + std(df[:, 1:], axis=1),
#                          mean(df[:, 1:], axis=1) - std(df[:, 1:], axis=1), color='darkred')
    ax2.set_xlim(0,25)
    ax2.set_ylim(0, max_force)
    ax2.set_xlabel(r'Distance [$\AA$]').set_fontsize(16)
    ax2.set_title('MFM-pocket distance')

    # #hydrogen bonds number vs time
    ax3.fill_between(ht[:, 0], mean(ht[:, 1:], axis=1) + std(ht[:, 1:], axis=1),
                          mean(ht[:, 1:], axis=1) - std(ht[:, 1:], axis=1), linewidth=1, alpha=0.5, color='orange')
    ax3.plot(ht[:,0],mean(ht[:,1:],axis=1),linewidth=2, color='orange')
    ax3.set_xlim(0,)
    ax3.set_ylim(0,max_HB)
    ax3.set_ylabel('HB number')
    ax3.set_xlabel('Time [ns]').set_fontsize(16)
    ax3.set_title('Hydrogen bonds')

fig.tight_layout()
plt.savefig(name+'/bestdir_mean.png',dpi=600)
plt.show()

