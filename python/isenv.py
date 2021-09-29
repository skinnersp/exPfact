from isenv_functions import *

ass_file   = "moprp.ass"
seq_file   = "moprp.seq"
T_label    = 298
pH_label   = 4
T_quench   = 273
pH_quench  = 2.4
lnP_file   = "test.pfact"
times_file = "exp.times"
exp_env_pre= "pep1"
pep = 1
z   = 1

if not os.path.exists(exp_env_pre):
    os.makedirs(exp_env_pre)

predict_isotopic_envelope(ass_file, seq_file, T_label, pH_label,
                          lnP_file, times_file, pep, z,
                          exchange = 'f', out_file = exp_env_pre+'/'+str(pep))

generate_back_exchange_time_points()
times = read_time_points(times_file)
exp_idxs = [i for i in range(1, len(times))]

for exp_idx in exp_idxs:

    predict_isotopic_envelope(ass_file, seq_file, T_quench, pH_quench,
                               lnP_file, "back.times", pep, z,
                               exchange = 'b', 
                               out_file = exp_env_pre+'/'+str(pep)+'.'+str(exp_idx)+'.corr',
                               pi0_file = exp_env_pre+'/'+str(pep)+'.'+str(exp_idx)+'.isot')
    
    corr_env = pd.read_csv(exp_env_pre+'/'+str(pep)+'.'+str(exp_idx)+'.corr'+'.1.isot',header=None,delim_whitespace=True,skiprows=1)
    exp_env  = pd.read_csv(exp_env_pre+'.'+str(exp_idx)+'.txt',header=None,delim_whitespace=True)

    mass, fr = sticks_from_exp_envelope(exp_env, corr_env, z)

    res = compare_predictions(fr, prefix=exp_env_pre+'/'+str(pep)+'.'+str(exp_idx)+'.corr'+'*')        
    res.to_csv(exp_env_pre+'/'+'tau.'+str(pep)+'.'+str(exp_idx)+'.res',sep='\t')
        