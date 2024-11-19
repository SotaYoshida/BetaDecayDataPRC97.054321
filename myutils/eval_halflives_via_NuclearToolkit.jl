using Glob

# using the latest version of NuclearToolkit 2024/11/19 in my local environment
include("src/NuclearToolkit.jl")
using .NuclearToolkit
# When you use the released version of NuclearToolkit.jl, please use the following line instead.
# using NuclearToolkit

function read_summary_mwav(pnuc, sntf, A)
    fn = "$(target_dir)/summary/summary_$(pnuc)$(A)_$(sntf)_mwav.txt"
    if isfile(fn)
        lines = readlines(fn)
        p_Jpi_candidates = String[]
        d_Jpi_candidates = Dict{String,Vector{String}}[ ]
        for line in lines
            if occursin("log_", line) && occursin(".txt", line)
                tl = split(line)
                J = tl[2]
                if occursin("/", J)
                    J = parse(Int, split(J, "/")[1])
                else
                    J = parse(Int, J) * 2
                end
                parity = tl[3]
                E = tl[6]
                Ex = tl[7]
                if parse(Float64, Ex) > 0.5
                    continue
                end
                parity = ifelse(parity == "+", "p", "n")
                Jpi = "j$(J)$(parity)"
                tdict = Dict{String,Vector{String}}() 
                tdict["GT"] = String[ ]
                tdict["FF"] = String[ ]
                for dJ = max(A%2, J-2):2:J+2
                    if J == 0 && dJ == 0
                        continue
                    end
                    push!(tdict["GT"], "j$(dJ)$(parity)")
                end
                op_parity = ifelse(parity == "p", "n", "p")
                for dJ = max(A%2, J-4):2:J+4
                    push!(tdict["FF"], "j$(dJ)$(op_parity)")
                end
                push!(p_Jpi_candidates, Jpi)
                push!(d_Jpi_candidates, tdict)
            end
        end
        return p_Jpi_candidates, d_Jpi_candidates
    else
        return error("File not found: $(fn)")
    end
end

"""
target_dir: str
    The directory where the log files are stored.
We assume that the following files are present in the directory:
- logfiles/summary/summary_Mg34_SDPFSDG_mwav.txt: summary file for the parent nucleus
- logfiles/transition_GT_LanczosStrength/log_(pnuc)(A)(p_Jpi)_(snt_name)_(dnuc)(A)j*n.txt: log files for GT transitions with Lanczos strength function method
- logfiles/transition_FF/log_(pnuc)(A)_(dnuc)(A)_(snt_name)_tr_j*_(p_Jpi).txt: log files for FF transitions with ordinary pipeline, diagonalizing low-lying states and evaluating transitions among them
"""
const target_dir = homedir()*"/Desktop/betadecay_K_to_Ca/logfiles/"; pnuc = "K"; dnuc ="Ca"; Anums = collect(41:57)

#const target_dir = homedir()*"/Desktop/betadecay_Mg_to_Al/logfiles/";pnuc = "Mg"; dnuc = "Al"; Anums = collect(34:4)
#const target_dir = homedir()*"/Desktop/betadecay_Cl_to_Ar/logfiles/" ;pnuc = "Cl"; dnuc ="Ar"; Anums = collect(41:52)


function run(;sumpath="$(target_dir)/summary/",
            GTpath="$(target_dir)/transition_GT_LanczosStrength/", 
            FFpath="$(target_dir)/transition_FF/")
    verbose= true # if turned on, all the transition details will be printed

    sntf = "sdpf-m"
    sntf = "SDPFSDG"
    
    for A = Anums
        if verbose
            println("\n*** $(A)$(pnuc) => $(A)$(dnuc)***")
        end
        fns_sum_parent = [ sumpath*"summary_$(pnuc)$(A)_$(sntf)_mwav.txt" ]
        fns_sum_daughter = [ sumpath*"summary_$(dnuc)$(A)_$(sntf)_mwav.txt" ]

        parent_Jpi_candidates, daughter_Jpi_candidates = read_summary_mwav(pnuc, sntf, A)

        num_cand = length(parent_Jpi_candidates)

        for n = 1:num_cand
            parent_Jpi = parent_Jpi_candidates[n]
            Jpis_GT = daughter_Jpi_candidates[n]["GT"]
            Jpis_FF = daughter_Jpi_candidates[n]["FF"]
            if verbose; println("Parent Jpi: $(parent_Jpi)"); end
            fns_GT = String[ ]
            fns_FF = String[ ]
            if verbose; print("--GT: "); end
            for d_Jpi in Jpis_GT
                fn_GT = "$(GTpath)log_$(pnuc)$(A)$(parent_Jpi)_$(sntf)_$(dnuc)$(A)$(d_Jpi).txt"
                push!(fns_GT, fn_GT)
                if verbose; print("  -> Jpi=$(d_Jpi) ", ifelse(isfile(fn_GT), "✓ ", "Not found: $fn_GT")); end
            end
            if verbose; 
                println()
                print("--FF: ")
            end
            for d_Jpi in Jpis_FF
                fn_FF = "$(FFpath)log_$(pnuc)$(A)_$(dnuc)$(A)_$(sntf)_tr_$(parent_Jpi)_$(d_Jpi).txt"
                push!(fns_FF, fn_FF)
                if verbose
                    print("  -> Jpi=$(d_Jpi) ", ifelse(isfile(fn_FF), "✓ ", "Not found: $fn_FF"))
                end
            end
            if verbose; println(); end        
            eval_betadecay_from_kshell_log(fns_sum_parent, fns_sum_daughter, fns_GT, fns_FF, parent_Jpi; verbose=verbose)
        end
    end
end
run()
