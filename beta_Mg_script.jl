using NuclearToolkit
#include("/Users/sotair/uu_gdrive/dev_NuclearToolkit.jl/src/NuclearToolkit.jl")
#using .NuclearToolkit
using Glob

function run()

    for A = 34:41
        fns_sum_parent = [ "summary_Mg$(A)_SDPFSDG.txt" ]
        fns_sum_daughter = [ "summary_Al$(A)_SDPFSDG.txt" ]
        parentJpi = "j0p"
        if A % 2 != 0
            for parentJpi in ["j1n", "j3n", "j5n"]
                fns_GT = glob("log_Mg$(A)$(parentJpi)_SDPFSDG_Al$(A)j*n.txt")
                fns_FF = glob("log_Al$(A)_Mg$(A)_SDPFSDG_tr_j*_$(parentJpi).txt")
                eval_betadecay_from_kshell_log(fns_sum_parent, fns_sum_daughter, fns_GT, fns_FF, parentJpi)
            end
        else
            fns_GT = glob("log_Mg$(A)j0p_SDPFSDG_Al$(A)j*p.txt")
            fns_FF = glob("log_Al$(A)_Mg$(A)_SDPFSDG_tr_j*")
            eval_betadecay_from_kshell_log(fns_sum_parent, fns_sum_daughter, fns_GT, fns_FF, parentJpi)
        end
    end
end
run();exit()

