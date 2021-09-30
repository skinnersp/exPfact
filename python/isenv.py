from isenv_functions import *

def corrected_isotopic_envelope_prediction(ass_file, seq_file, T_label, pH_label,
                                           T_quench, pH_quench, lnP_file, times_file,
                                           exp_env_pre, pep, z):

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
        
if __name__ == '__main__':
    
    import argparse    
    parser = argparse.ArgumentParser()

    parser.add_argument("--ass")
    parser.add_argument("--seq")
    parser.add_argument("--T_label")
    parser.add_argument("--pH_label")
    parser.add_argument("--T_quench")
    parser.add_argument("--pH_quench")
    parser.add_argument("--pfact")
    parser.add_argument("--times")
    parser.add_argument("--prefix")
    parser.add_argument("--pep")
    parser.add_argument("--z")
    
    parser.add_argument("--mode")
   
    opts = parser.parse_args()
    
    if opts.ass:
        ass_file = opts.ass
    if opts.seq:
        seq_file = opts.seq
    if opts.T_label:
        T_label = float(opts.T_label)
    if opts.pH_label:
        pH_label = float(opts.pH_label)
    if opts.T_quench:
        T_quench = float(opts.T_quench)
    if opts.pH_quench:
        pH_quench = float(opts.pH_quench)
    if opts.pfact:
        lnP_file = opts.pfact
    if opts.times:
        times_file = opts.times
    if opts.prefix:
        exp_env_pre = opts.prefix
    if opts.pep:
        pep = int(opts.pep)
    if opts.z:
        z = int(opts.z)
    
    if not os.path.exists(exp_env_pre):
        os.makedirs(exp_env_pre)
    
    if opts.mode == 'p':
        predict_isotopic_envelope(ass_file, seq_file, T_label, pH_label,
                          lnP_file, times_file, pep, z,
                          exchange = 'f', out_file = exp_env_pre+'/'+str(pep))
    elif opts.mode == 'c':
        corrected_isotopic_envelope_prediction(ass_file, seq_file, T_label, pH_label,
                          T_quench, pH_quench, lnP_file, times_file, exp_env_pre, pep, z)
    else:
        print("Usage:")
        print("* Predict mode (p)")
        print("  isenv.py -mode p  --ass <ass_file>     --seq <seq_file>     --T_label <T_label> --pH_label <pH_label>")
        print("                    --pfact <pfact_file> --times <times_file> --pep <pep_index>   --z <charge state>")
        print("                    --prefix <out_prefix>")
        print("* Comparison mode (c)")
        print("  isenv.py -mode c  --ass <ass_file>      --seq <seq_file>        --T_label <T_label> --pH_label <pH_label>")
        print("                    --T_quench <T_quench> --pH_quench <pH_quench> --pfact <pfact_file> --times <times_file>")
        print("                    --pep <pep_index>   --z <charge state>        --prefix <prefix of exp envelopes>")