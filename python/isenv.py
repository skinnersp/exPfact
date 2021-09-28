from isenv_functions import *


pep = 1
z = 1
pfile = 'c1'

''' STEP ONE
    Evolve the theor fully deut to approach the exp fully deut
    Select the most similar as starting point 
    (STILL TO BE DEVELOPED) '''

prospector_file = str(pep)+'.prospector.ucsf.edu'

fully_prot = pd.read_csv('pep'+str(pep)+'.exp1.txt',header=None,
                         delim_whitespace=True)
theor0 = pd.read_csv(prospector_file,header=None,
                     delim_whitespace=True,skiprows=1)

mass = list(theor0[1])
fr = np.zeros(len(mass))
for i in range(len(fully_prot)):
    for j in range(len(mass)):
        if float(fully_prot[0][i])>mass[j]-0.5 and float(fully_prot[0][i])<mass[j]+0.5:
            fr[j] += fully_prot[1][i]
fr = [fr[j]/sum(fr)*100 for j in range(len(fr))]
            
fully_prot_file = str(pep)+'.fullyprot.isot'
with open(fully_prot_file,'w+') as f:
    f.write('#comment line\n')
    for i in range(0,len(mass)):
        f.write('%d ' % i)
        f.write('%5.5f ' % mass[i])
        f.write('%5.2f ' % fr[i])
        ren = fr[i]*100/max(fr)
        if i == len(mass)-1:
            f.write('%5.2f' % ren)
        else:
            f.write('%5.2f\n' % ren)
f.close()

#os.system('python iso.py --ass moprp --pep '+str(pep)+' --kint moprp --pfact best --times moprp --pi '+fully_prot_file+' --out step1')

''' add a script that selects the best envelope '''
selected = prospector_file
#selected = str(pep)+'.fullyprot.isot'

#plt.figure()
#fwd_prot3 = pd.read_csv(selected,header=None,sep='\t',skiprows=1)
#plt.plot(mass,fr,label='experimental')
#plt.plot(mass,theor0[2],label='theoretical')
#plt.plot(fwd_prot3[1],fwd_prot3[2],label='fwd fully prot')
#plt.legend()


''' STEP TWO
    Evolve from selected initial envelope to experimentaltimes available '''
os.system('python iso.py --ass moprp --pep '+str(pep)+' --kint moprp --pfact '+str(pfile)+' --times exp --pi '+selected+' --out step2')

''' STEP THREE
    Back evolve and select the most appropriate time '''

exp_idxs = [1,2,3,4]

start = -5
finish = 5
replicate = 500
with open('back.times','w+') as f:
    a = np.logspace(start,finish,replicate)
    for i in range(0,len(a)):
        f.write('%5.10f\n' % a[i])
f.close()

for exp_idx in exp_idxs:
    #plt.clf()
    corr_env = pd.read_csv('step2.'+str(exp_idx)+'.isot',header=None,
                           delim_whitespace=True,skiprows=1)
    exp_env = pd.read_csv('pep'+str(pep)+'.exp'+str(exp_idx+1)+'.txt',header=None,
                          delim_whitespace=True)

    mass = list(corr_env[1])
    fr = np.zeros(len(mass))
    for i in range(len(exp_env)):
        for j in range(len(mass)):
            if exp_env[0][i]*z>=mass[j]-0.5 and exp_env[0][i]*z<mass[j]+0.5:
                fr[j] += exp_env[1][i]
    fr = fr/sum(fr)*100
           
    os.system('python isoback.py --ass moprp --pep '+str(pep)+' --kint quench --pfact '+str(pfile)+' --times back --pi step2.'+str(exp_idx)+'.isot --out step3.'+str(exp_idx))

    files = []
    for file in glob.glob('*step3.'+str(exp_idx)+'*'):
        files.append(file)
    sort_nicely(files)

    with open('tau.'+str(exp_idx)+'.res','w+') as f:
        for i in range(0,len(files)):
            back_env = pd.read_csv(files[i],header=None,sep='\t',skiprows=1)
            times = pd.read_csv('back.times',header=None)

            c_corr = (np.average(mass,weights=back_env[2]/sum(back_env[2])*100)-np.average(mass,weights=fr))**2
            c_corr = r2_score(fr,back_env[2]/sum(back_env[2])*100)
    
            f.write('%s ' % files[i])
            f.write('%5.10f ' % times[0][i])
            f.write('%5.10f\n' % c_corr)
            #print(files[i],round(times[0][i],5),round(c_corr,3))
    print('************************************')
    f.close()

    tau = pd.read_csv('tau.'+str(exp_idx)+'.res',header=None,sep=' ')
        
    minimum = tau[1][tau[2].idxmax()]
    min_file = tau[0][tau[2].idxmax()]
    
    plt.figure()    
    sel_env = pd.read_csv(min_file,header=None,delim_whitespace=True,skiprows=1)
    plt.plot(sel_env[1],sel_env[2]/sum(sel_env[2])*100,label='pred')
    plt.plot(mass,fr,label='exp')
    plt.ylim(0,60)
    plt.legend()
    
    sel_env_norm = list(sel_env[2]/sum(sel_env[2])*100)
    
    with open(pfile+'.'+str(pep)+'.'+str(exp_idx)+'.best.isot','w+') as f:
        for i in range(len(mass)):
            f.write('%5.5f ' % mass[i])
            f.write('%5.5f ' % fr[i])
            f.write('%5.5f\n' % sel_env_norm[i])
    f.close()
    
    plt.figure()
    plt.plot(tau[1],tau[2],'o-')
    plt.xscale('log')    
    plt.xlabel('Back Time Evolution',fontsize=15)
    plt.ylabel('Distance',fontsize=15)
    plt.pause(5)        
    
# %%    