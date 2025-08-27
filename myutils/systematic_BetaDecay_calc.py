import glob
import os
import re
import sys
import time

element = [
    'NA', 
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]

def generate_sh_header():
    text = "#!/bin/sh" + "\n"
    #text += "export GFORTRAN_UNBUFFERED_PRECONNECTED=y" + "\n"
    if using_O:        
        text += '#PJM -L "rscgrp=regular-o"'+ "\n"
        text += '#PJM -g "ge12"' + "\n"
        text += '#PJM -L "node=64"'+ "\n"
        text += '#PJM --mpi "max-proc-per-node=1"' + "\n"
        text += '#PJM -L "elapse=12:00:00"'  + "\n"
        text += '#PJM -j' + "\n"
        text += '#PJM -s' + "\n"
        text += "#ulimit -s unlimited" + "\n"
        text += "export OMP_NUM_THREADS=48" + "\n"
        text += "export OMP_STACKSIZE=128m" + "\n"
    else:
        text += '#PJM -L "rscgrp=regular-a"'+ "\n"
        text += '#PJM -L "node=4"'+ "\n"
        text += '#PJM -L "elapse=24:00:00"'  + "\n"
        text += "ulimit -s unlimited" + "\n"
        text += "export OMP_NUM_THREADS=36" + "\n"
        #text += "#export OMP_STACKSIZE=128m" + "\n"
    return text    

def generate_parent_candidates(nuc, jpi, ptnf, J2val, n_eigen, sntf):
    oup = open("prepare_parent_candidates.sh", "w")
    text = generate_sh_header()
    text += "echo \"preparing " + nuc + "_"+jpi+" for transitions\"" + "\n"
    text += "cat > " + nuc + ".input <<EOF" + "\n"
    text += "&input" + "\n"
    text += "  beta_cm = 10.0" + "\n"
    text += "  eff_charge = 1.5, 0.5," + "\n"
    text += "  fn_int = "+'"'+sntf+'"' + "\n"
    text += "  fn_ptn = " + '"' + ptnf + '"' + "\n"
    text += "  fn_save_wave = " + '"' + nuc + "_"+ sntf.split(".")[0]+ "_" + str(jpi) + ".wav" + '"' + "\n"
    text += "  gl = 1.0, 0.0," + "\n"
    text += "  gs = 3.91, -2.678," + "\n"
    text += "  hw_type = 2" + "\n"
    text += "  is_double_j = .true." + "\n"
    text += "  max_lanc_vec = 300" + "\n"
    text += "  maxiter = 300" + "\n"
    if using_O:
        text += "  mode_lv_hdd = 1" + "\n"
    else:
        text += "  mode_lv_hdd = 0" + "\n"
    text += "  mtot = " + str(J2val) + "\n"
    text += "  n_eigen = 1\n"
    text += "  n_restart_vec = 30" + "\n"
    text += "&end" + "\n"
    text += "EOF" + "\n"
    if using_MPI:
        text += 'mpiexec -of ' +  'log_'+nuc+'_'+sntf.split(".")[0]+'_'+ str(jpi) +'.txt  ./kshell.exe '+nuc+'.input ' + "\n"
    else:
        text += './kshell.exe '+nuc+'.input > log_'+nuc+'_'+sntf.split(".")[0]+'_'+ str(jpi) +'.txt' + "\n"
    text += 'rm '+nuc+'.input' + "\n"
    text += 'rm -f tmp_lv_'+nuc+'*  tmp_snapshot_'+nuc+'*' + "\n"
    oup.write(text)
    oup.close()
    time.sleep(5)
    os.system("bash prepare_parent_candidates.sh")

def generate_daughter_states(pnuc, dnuc, jpi_parent, jpi_daughter, ptn_parent, ptn_daughter, mtot, n_eigen, sntf, make_mwav=False):
    # generate daughter lowest states for both parity
    text = ""
    if make_mwav:
        if mtot % 2 == 0:
            m_oe = 0
        else:
            m_oe = 1
        text = "\n\n#" + "%"*50 + "\n"
        text += "echo \"preparing " + dnuc + " mwav\"" + "\n"
        for parity in ["p", "n"]:
            text += "cat > " + dnuc + ".input <<EOF" + "\n"
            text += "&input" + "\n"
            text += "  beta_cm = 10.0" + "\n"
            text += "  eff_charge = 1.5, 0.5," + "\n"
            text += "  fn_int = "+'"'+sntf+'"' + "\n"
            text += "  fn_ptn = " + '"' + dnuc+"_"+sntf.split(".")[0] +"_" + parity + '.ptn"' + "\n"
            text += "  gl = 1.0, 0.0," + "\n"
            text += "  gs = 3.91, -2.678," + "\n"
            text += "  hw_type = 2" + "\n"
            text += "  is_double_j = .false." + "\n"
            text += "  max_lanc_vec = 150" + "\n"
            text += "  maxiter = 150" + "\n"
            if using_O:
                text += "  mode_lv_hdd = 1" + "\n"
            else:
                text += "  mode_lv_hdd = 0" + "\n"
            text += "  fn_save_wave = " + '"' + dnuc + "_"+ sntf.split(".")[0]+ "_m"+str(m_oe)+parity+".wav" + '"' + "\n"
            text += "  mtot = " + str(m_oe) + "\n"
            text += "  n_eigen = 5 \n"
            text += "  n_restart_vec = 30" + "\n"
            text += "&end" + "\n"
            text += "EOF" + "\n"
            if using_MPI:
                text += 'mpiexec -of ' + 'log_'+dnuc+'_'+sntf.split(".")[0]+'_'+'m'+str(m_oe)+parity+'.txt  ./kshell.exe '+nuc+'.input ' + "\n"
            else:
                text += './kshell.exe '+dnuc+'.input > log_'+dnuc+'_'+sntf.split(".")[0]+'_'+'m'+str(m_oe)+parity+'.txt' + "\n"
            text += 'rm '+dnuc+'.input' + "\n"
            #text += 'rm tmp_lv_'+dnuc+'*' + "\n"
            text += 'rm -f tmp_lv_'+dnuc+'*  tmp_snapshot_'+dnuc+'*' + "\n"
    # generate daughter states fed by GT/FF transitions
    text += "\n\n#" + "%"*50 + "\n"
    text += "echo \"preparing " + dnuc + "("+ jpi_daughter + ") for transitions\"" + "\n"
    text += "cat > " + dnuc + ".input <<EOF" + "\n"
    text += "&input" + "\n"
    text += "  beta_cm = 10.0" + "\n"
    text += "  eff_charge = 1.5, 0.5," + "\n"
    text += "  fn_int = "+'"'+sntf+'"' + "\n"
    text += "  fn_ptn = " + '"' + ptn_daughter + '"' + "\n"
    text += "  fn_save_wave = " + '"' + dnuc + "_"+ sntf.split(".")[0]+ "_" + jpi_daughter + ".wav" + '"' + "\n"
    text += "  gl = 1.0, 0.0," + "\n"
    text += "  gs = 3.91, -2.678," + "\n"
    text += "  hw_type = 2" + "\n"
    text += "  is_double_j = .true." + "\n"
    text += "  max_lanc_vec = 300" + "\n"
    text += "  maxiter = 300" + "\n"
    if using_O:
        text += "  mode_lv_hdd = 1" + "\n"
    else:
        text += "  mode_lv_hdd = 0" + "\n"
    text += "  mtot = " + str(mtot) + "\n"
    text += "  n_eigen = " + str(n_eigen) + " \n"
    text += "  n_restart_vec = 30" + "\n"
    text += "&end" + "\n"
    text += "EOF" + "\n"
    if using_MPI:
        text += 'mpiexec -of ' + 'log_'+dnuc+'_'+sntf.split(".")[0]+'_'+ jpi_daughter+'.txt  ./kshell.exe '+dnuc+'.input ' + "\n"
    else:
        text += './kshell.exe '+dnuc+'.input > log_'+dnuc+'_'+sntf.split(".")[0]+'_'+ jpi_daughter +'.txt' + "\n"
    text += 'rm '+dnuc+'.input' + "\n"
    #text += 'rm tmp_lv_'+dnuc+'*' + "\n"
    text += 'rm -f tmp_lv_'+dnuc+'*  tmp_snapshot_'+dnuc+'*' + "\n"
    return text

def generate_transition_sh(pnuc, dnuc, ptn_parent, ptn_daughter, jpi_parent, jpi_daughter, mtot, n_eigen, sntf):
    text =  "echo \"transition from " + pnuc + "("+ jpi_parent+ ") to " + dnuc +"("+jpi_daughter +") \"\n"
    text += "cat > " + dnuc + ".input <<EOF" + "\n"
    text += "&input" + "\n"
    text += "   fn_int = "+'"'+sntf+'"' + "\n"
    text += "   fn_ptn_l = " + '"' + ptn_daughter + '"' + "\n"
    text += "   fn_ptn_r = " + '"' + ptn_parent + '"' + "\n"
    text += "   fn_load_wave_l = " + '"' + dnuc + "_"+ sntf.split(".")[0]+ "_" + jpi_daughter + ".wav" + '"' + "\n"
    text += "   fn_load_wave_r = " + '"' + pnuc + "_"+ sntf.split(".")[0]+ "_" + jpi_parent + ".wav" + '"' + "\n"
    text += "   hw_type = 2" + "\n"
    text += "   eff_charge = 1.5, 0.5" + "\n"
    text += "   is_obtd = .false." + "\n"
    text += "   is_tbtd = .false." + "\n"
    text += "   gl = 1.0, 0.0" + "\n"
    text += "   gs = 3.91, -2.678" + "\n"
    text += "&end" + "\n"
    text += "EOF" + "\n"
    text += './transit.exe '+dnuc+'.input > log_'+pnuc+'_'+dnuc+"_"+sntf.split(".")[0]+'_tr_'+ jpi_parent + "_" + jpi_daughter +'.txt' + "\n"
    return text

def generate_GT_LS_sh_body(pnuc, dnuc, ptn_parent, ptn_daughter, jpi_parent, jpi_daughter, mtot, n_eigen, sntf):
    text = "echo \"preparing " + pnuc + "("+ jpi_parent + ") =>" + dnuc + "("+ jpi_daughter + ") gs jwav\"" + "\n"
    text += "cat > " + pnuc + ".input <<EOF" + "\n"
    text += "&input" + "\n"
    text += "  beta_cm = 10.0" + "\n"
    text += "  eff_charge = 1.5, 0.5," + "\n"
    text += "  fn_int = "+'"'+sntf+'"' + "\n"
    text += "  fn_ptn_init = " + '"' + ptn_parent + '"' + "\n"
    text += "  fn_ptn = " + '"' + ptn_daughter + '"' + "\n"
    text += "  op_type_init = 'GT'"
    text += "  gl = 1.0, 0.0," + "\n"
    text += "  gs = 3.91, -2.678," + "\n"
    text += "  hw_type = 2" + "\n"
    text += "  is_double_j = .true." + "\n"
    text += "  max_lanc_vec = 300" + "\n"
    text += "  maxiter = 1" + "\n"
    text += "  mode_lv_hdd = 2" + "\n"
    text += "  fn_load_wave = " + '"' + pnuc + "_"+ sntf.split(".")[0]+ "_" + jpi_parent + ".wav" + '"' + "\n"
    text += "  mtot = " + str(mtot) + "\n"
    text += "  n_eigen = " + str(n_eigen) + "\n"
    text += "  n_restart_vec = 30" + "\n"
    text += "&end" + "\n"
    text += "EOF" + "\n"
    if using_MPI:
        text += 'mpiexec -of ' +  'log_'+pnuc+jpi_parent+'_'+sntf.split(".")[0]+'_'+dnuc+jpi_daughter+'.txt  ./kshell.exe '+pnuc+'.input ' + "\n"
    else:
        text += './kshell.exe '+pnuc+'.input > log_'+pnuc+jpi_parent+'_'+sntf.split(".")[0]+'_'+dnuc+jpi_daughter+'.txt' + "\n"
    text += 'rm '+pnuc+'.input' + "\n"
    text += 'rm -f tmp_lv_'+pnuc+'*  tmp_snapshot_'+pnuc+'*' + "\n"
    text += "#" + "%"*50 + "\n\n"
    return text

def get_possible_j(Jparent, mode):
    if mode == "FF":
        Jmax = Jparent + 2 * 2
        if Jparent % 2 == 0: #even
            Jmin = max(Jparent - 4, 0)
        else: #odd
            Jmin = max(Jparent - 4, 1)
    elif mode == "GT" or mode == "GT_LS":
        Jmax = Jparent + 1 * 2
        if Jparent % 2 == 0:
            Jmin = max(Jparent - 2, 0)
        else:
            Jmin = max(Jparent - 2, 1)
        if Jparent == 0:
            Jmin = 2            
    possible_j = list(range(Jmin, Jmax+1, 2))
    return possible_j

def get_dnuc_from_pnuc(pnuc):
    Anum = int(re.findall(r'\d+', pnuc)[0])
    el = pnuc.replace(str(Anum), "")
    Z = element.index(el)
    dnuc = element[Z+1] + str(Anum)
    return dnuc

def read_summary(pnuc, sntf, Ethreshold=0.5):
    fn = "summary_" + pnuc + "_" + sntf.split(".")[0]+ "_mwav.txt"
    f = open(fn, "r")
    lines = f.readlines()
    f.close()
    jpi_parents = [ ]
    for line in lines:
        if "log_" in line:
            tmp = line.strip().split()
            n, Jpi_raw, prty, N_Jp, T, E, Ex, logf = tmp
            prty = "p" if prty == "+" else "n"
            if "/" in Jpi_raw:
                Jpi = Jpi_raw.split("/")[0] + prty
            else:
                Jpi = str(int(Jpi_raw)*2) + prty
            if float(Ex) > Ethreshold:
                continue
            jpi_text = "j"+Jpi
            if jpi_text not in jpi_parents:
                jpi_parents.append(jpi_text)            
    return jpi_parents

def generate_jwav_transition_sh_files(sntf, mode, ZNrange):
    if mode == "GT_LS":
        n_eigen = 300
    elif mode == "GT":
        n_eigen = 5          
    elif mode == "FF":
        n_eigen = 15
    else:
        print("Error: mode not recognized")
        exit()        
    if mode == "GT_LS":
        oup1 = open("eval_GT_LanczosStrength.sh", "w")
    else:
        oup1 = open("prepare_daughter_states_"+mode+".sh", "w")    
    oup2 = open("eval_transitions_"+mode+".sh", "w")
    oup1.write(generate_sh_header())
    oup2.write(generate_sh_header())
    for (Z, N) in ZNrange:
        pnuc = element[Z] + str(Z+N)
        dnuc = element[Z+1] + str(Z+N)

        natural_parity = "p" if N % 2 == 0 else "n" 
        jpis_parent = read_summary(pnuc, sntf)
        jpi_daughter_pool = [ ]
        for jpi_parent in jpis_parent:
            Jparent = int(jpi_parent.replace("j","").replace("p","").replace("n",""))            
            dnuc = get_dnuc_from_pnuc(pnuc)              
            pi_parent = pi_daughter = jpi_parent[-1]
            ptn_parent = pnuc + "_" + sntf.split(".")[0] + "_" + pi_parent + ".ptn"
            ptn_daughter = ptn_parent.replace(pnuc, dnuc)
            wav_fn = pnuc + "_"+ sntf.split(".")[0]+ "_" + jpi_parent + ".wav"            
            if len(glob.glob(wav_fn)) == 0:
                print("Error: ", wav_fn, " not found! constructing the jwav files...")
                generate_parent_candidates(pnuc, jpi_parent, ptn_parent, Jparent, n_eigen, sntf)
            if mode == "FF":
                if pi_parent == "n":
                    pi_daughter = "p"
                else:
                    pi_daughter = "n"
            ptn_daughter = dnuc + "_" + sntf.split(".")[0] + "_" + pi_daughter + ".ptn"
            possible_j = get_possible_j(Jparent, mode)
            for mtot in possible_j:                
                jpi_daughter = "j"+str(mtot)+pi_daughter
                if mode == "GT_LS":
                    sh_body = generate_GT_LS_sh_body(pnuc, dnuc, ptn_parent, ptn_daughter, jpi_parent, jpi_daughter, mtot, n_eigen, sntf)
                    oup1.write(sh_body)
                else:                    
                    if not(jpi_daughter in jpi_daughter_pool):                    
                        sh_body = generate_daughter_states(pnuc, dnuc, jpi_parent, jpi_daughter, ptn_parent, ptn_daughter, mtot, n_eigen, sntf)
                        oup1.write(sh_body)
                        jpi_daughter_pool.append(jpi_daughter)
                    tr_sh = generate_transition_sh(pnuc, dnuc, ptn_parent, ptn_daughter, jpi_parent, jpi_daughter, mtot, n_eigen, sntf)
                    oup2.write(tr_sh) 

    oup1.close()
    oup2.close()
    if mode == "GT_LS":
        fn = "eval_transitions_"+mode+".sh"
        os.system("rm "+fn)

def generate_scriptfiles(sntf, ZNrange, n_eigen=5, run_mwav_job=False):
    """
    Fisrt, calculate a couple of low-lying states for parent and daughter nuclei.
    """
    oup = open("prepare_parent_daughter_lowlying_states.sh", "w")
    header = generate_sh_header()
    oup.write(header)
    for (Z, N) in ZNrange:
        A = Z + N
        pnuc = element[Z] + str(A)
        dnuc = element[Z+1] + str(A)
        m_oe = A % 2
        for pd in ["parent", "daughter"]:
            if pd == "parent":
                nuc = pnuc
                n_eigen = 5
            else:
                nuc = dnuc
                n_eigen = 1
            sumf = "summary_" + nuc + "_" + sntf.split(".")[0]+ "_mwav.txt"
            text = "#" + "%"*50 + "\n\n"
            text += "echo \"preparing " + nuc + " mwav\"" + "\n"
            for parity in ["p", "n"]:
                text += "cat > " + nuc + ".input <<EOF" + "\n"
                text += "&input" + "\n"
                text += "  beta_cm = 10.0" + "\n"
                text += "  eff_charge = 1.5, 0.5," + "\n"
                text += "  fn_int = "+'"'+sntf+'"' + "\n"
                text += "  fn_ptn = " + '"' + nuc+"_"+sntf.split(".")[0] +"_" + parity + '.ptn"' + "\n"
                text += "  gl = 1.0, 0.0," + "\n"
                text += "  gs = 3.91, -2.678," + "\n"
                text += "  hw_type = 2" + "\n"
                text += "  is_double_j = .false." + "\n"
                text += "  max_lanc_vec = 150" + "\n"
                text += "  maxiter = 150" + "\n"
                if using_O:
                    text += "  mode_lv_hdd = 1" + "\n"
                else:
                    text += "  mode_lv_hdd = 0" + "\n"
                text += "  mtot = " + str(m_oe) + "\n"
                text += "  n_eigen = "+str(n_eigen)+" \n"
                text += "  n_restart_vec = 50" + "\n"
                text += "&end" + "\n"
                text += "EOF" + "\n"
                if using_MPI:
                    text += 'mpiexec -of ' +  'log_'+nuc+'_'+sntf.split(".")[0]+'_'+'m'+str(m_oe)+parity+'.txt  ./kshell.exe '+nuc+'.input ' + "\n"
                else:
                    text += './kshell.exe '+nuc+'.input > log_'+nuc+'_'+sntf.split(".")[0]+'_'+'m'+str(m_oe)+parity+'.txt' + "\n"
                text += 'rm '+nuc+'.input' + "\n"
                text += "rm tmp_lv_"+nuc+"* \n\n"
            text += "python3 collect_logs.py log_"+nuc+"_"+sntf.split(".")[0]+"_m*.txt > " + sumf + "\n"
            oup.write(text)
    oup.close()

    if run_mwav_job:
        os.system("bash prepare_parent_daughter_lowlying_states.sh")


def generate_0hw1hw_ptn_files(ZNrange, sntf, hw_ofst:int=0):
    for (Z,N) in ZNrange:
        A = Z + N
        pnuc = element[Z] + str(A)
        dnuc = element[Z+1] + str(A)
        print("Z", Z, "N", N, "pnuc", pnuc, "dnuc",dnuc)
        for parity in ["+", "-"]:
            c = "p" if parity == "+" else "n"
            ptn_f = pnuc + "_" + sntf.split(".")[0] + "_"+c+".ptn"
            os.system("python3 gen_partition_local.py "+ sntf+ " "+ ptn_f + " "
                      + str(Z-Zcore)+ " " + str(N-Ncore)+ " " + parity+ "  "+str(hw_ofst))
            c_d = "p" if c == "n" else "n"
            ptn_f_d = dnuc + "_" + sntf.split(".")[0] + "_"+c_d+".ptn"
            parity_d = "+" if parity == "-" else "-"
            os.system("python3 gen_partition_local.py "+ sntf+ " "+ ptn_f_d + " "
                      + str(Z-Zcore+1)+ " " + str(N-Ncore-1)+ " " + parity_d + "  "+str(hw_ofst))

class BetaDecay:
    def __init__(self, Z, N, sntf, hw_ofst:int=0):
        self.Z = Z
        self.N = N
        self.A = Z + N
        self.sntname = sntf.split(".")[0]
        self.pnuc = element[Z] + str(self.A)
        self.dnuc = element[Z+1] + str(self.A)
        self.maindir = f"logfiles_{self.pnuc}_{self.dnuc}_{self.sntname}"
        if hw_ofst != 0:
            self.maindir = f"logfiles_{self.pnuc}_{self.dnuc}_{self.sntname}_{hw_ofst}hw"
        self.parent_parent_ptn = self.pnuc + "_" + self.sntname + "_p.ptn"
        self.parent_daughter_ptn_p = self.dnuc + "_" + self.sntname + "_p.ptn"
        self.parent_daughter_ptn_n = self.dnuc + "_" + self.sntname + "_n.ptn"
        self.dir_summary = self.maindir + "/summary"
        self.dir_parent_states = self.maindir + "/parent_states"
        self.dir_daughter_states = self.maindir + "/daughter_states"
        self.dir_transition_GT = self.maindir + "/transition_GT"
        self.dir_transition_GT_LS = self.maindir + "/transition_GT_LanczosStrength"
        self.dir_transition_FF = self.maindir + "/transition_FF"

    def create_dirs(self):
        mkdirs = [self.maindir, self.dir_summary, self.dir_parent_states, self.dir_daughter_states,
                   self.dir_transition_GT, self.dir_transition_GT_LS, self.dir_transition_FF]
        for d in mkdirs:
            os.mkdir(d)

def copy_exe_files(sntf):
    candidate = glob.glob( os.environ['HOME']+"/kshell-*")
    if len(candidate) == 0:
        print("Error: kshell not found!")
    else:
        kshell_dir = candidate[0]
        print("kshell found in ", kshell_dir)
        os.system("cp " + kshell_dir + "/bin/kshell.exe ./")
        os.system("cp " + kshell_dir + "/bin/transit.exe ./")
        os.system("cp " + kshell_dir + "/bin/collect_logs.py ./")
        if not(os.path.exists(sntf)):
            print("copying snt file from kshell")
            os.system("cp " + kshell_dir + "/snt/" + sntf + " ./")

if __name__ == "__main__":
    
    sntf = "sdpf-m.snt"
    #sntf = "SDPFSDG.snt"    
    copy_exe_files(sntf)


    Zcore = Ncore = 8

    ZNrange = [ (Z, N) for Z in [17] for N in range(24, 36)]
    ZNrange = [ (9,22) ]

    #ZNrange = [ (11, 26) ]


    T = True
    F = False
    mkdir = F
    rewrite_jwav = F
    using_MPI = F 

    # if len(sys.argv) < 2:

    #     print("usage: systematic_BetaDecay_calc.py O/A")
    #     exit()

    # if sys.argv[1].strip() == "O":
    #     using_O = True
    # else:
    #     using_O = False
    # print("using... ", sys.argv[1], "using_O", using_O)

    using_O = False

    hw_ofst = 0
    
    for (Z, N) in ZNrange:
        beta_decay = BetaDecay(Z, N, sntf, hw_ofst=hw_ofst)

        generate_0hw1hw_ptn_files(ZNrange, sntf, hw_ofst=hw_ofst)
        beta_decay.create_dirs()

        generate_scriptfiles(sntf, ZNrange, run_mwav_job=True)
            
        for mode in ["FF", "GT", "GT_LS"]:
            generate_jwav_transition_sh_files(sntf, mode, ZNrange)        


        if T:
            os.system("bash prepare_daughter_states_GT.sh")
            os.system("bash prepare_daughter_states_FF.sh")
            time.sleep(1)

            os.system("bash eval_transitions_GT.sh")
            os.system("mv log_*_tr_*.txt "+beta_decay.dir_transition_GT+"/")
            time.sleep(1)
            os.system("bash eval_GT_LanczosStrength.sh")
            os.system(f"mv log_{beta_decay.pnuc}*_{beta_decay.sntname}_{beta_decay.dnuc}*.txt "+beta_decay.dir_transition_GT_LS+"/")
            time.sleep(1)
    
            os.system("bash eval_transitions_FF.sh")
            os.system("mv log_*_tr_*.txt "+beta_decay.dir_transition_FF+"/")
            os.system(f"mv summary_{beta_decay.dnuc}_{beta_decay.sntname}_mwav.txt "+beta_decay.dir_summary+"/")
            os.system(f"mv summary_{beta_decay.pnuc}_{beta_decay.sntname}_mwav.txt "+beta_decay.dir_summary+"/")

            time.sleep(1)
            os.system(f"mv log_{beta_decay.pnuc}_{beta_decay.sntname}_*.txt "+beta_decay.dir_parent_states+"/")
            time.sleep(1)
            os.system(f"mv log_{beta_decay.dnuc}_{beta_decay.sntname}_*.txt "+beta_decay.dir_daughter_states+"/")


    os.system("rm *.exe")
    os.system("rm *.input")
    os.system("rm *.ptn")  
    os.system("rm *.wav*")
    os.system("rm *.sh")