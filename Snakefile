from os.path import join

configfile: "./config.yaml"

result_dir = config["inputs"]["result_dir"]
intermediate_results_dir = join(result_dir, "intermediate_results/")

rule all:
    input:
        join(config["inputs"]["result_dir"], "timeseries_network.txt")

rule make_final_result:
    input:
        expand(join(result_dir, "TG", "subnetwork.{i}"), \
                i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])

    output:
        join(config["inputs"]["result_dir"], "timeseries_network.txt")
    
    shell:
        "python src/makeTGDesc.py " + result_dir

rule extract_target_gene:
    input:
        nwk = join(intermediate_results_dir, config["parameters"]["prefix"] + ".nwk.t{i}"),
        DEGli = join(intermediate_results_dir, config["parameters"]["prefix"] + ".DEGli.t{i}"),
        tfranktrim = join(intermediate_results_dir, config["parameters"]["prefix"] + ".TF_rank.t{i}.trim")
    
    params:
        tf_li_file = config["inputs"]["tf_li_file"],
        output_dir = join(result_dir, "TG")

    output:
        join(result_dir, "TG", "subnetwork.{i}")

    shell:
        "python src/Target_genes.py --nwk {input.nwk} --deg {input.DEGli} --tf {params.tf_li_file} \
        --tfrank {input.tfranktrim} --out {params.output_dir}/subnetwork.{wildcards.i}"

rule extract_optimal_tf_list:
    input:
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".kccagrn.tp{i}.tsv"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])

    params:
        tf_li_file = config["inputs"]["tf_li_file"],
        exp_file = config["inputs"]["exp_file"],
        seed_file = config["inputs"]["seed_file"],
        n_threads = config["parameters"]["nThreads"],
        prefix = config["parameters"]["prefix"],
        temp = intermediate_results_dir,
        n_sample = config["parameters"]["nTimePoints"]

    output:
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".TF_rank.t{i}.trim"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)]),
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".nwk.t{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)]),
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".DEGli.t{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])

    shell:
        "python src/main_propanet.py -TFliFile {params.tf_li_file} \
        -nwkFile {params.temp}/{params.prefix}.tp -expFile {params.exp_file} \
        -binFile {params.seed_file} -N {params.n_sample} -p {params.n_threads} \
        -cond {params.prefix} -outD {params.temp}"

rule convert_pickle_to_tsv:
    input:
        join(intermediate_results_dir, config["parameters"]["prefix"] + ".kccagrn.tp{i}.pkl")
    
    output:
        join(intermediate_results_dir, config["parameters"]["prefix"] + ".kccagrn.tp{i}.tsv")

    run:
        import pickle as pkl

        with open(input, "rb") as f:
            data = pkl.load(f)

        with open(output, "w") as f:
            for key, val in data.items():
                for start, end in val:
                    f.write(f"{start}\t{end}\n")

rule inference_gene_regulatory_network:
    input:
        dict_cluster_TF = join(intermediate_results_dir, "/TFTG_cluster/community_TF.pkl"),
        dict_cluster_TG = join(intermediate_results_dir, "/TFTG_cluster/community_TG.pkl"),
        template_network = config["inputs"]["grn_file"],
        tf_list = config["inputs"]["tf_li_file"],

        gene_exp_profile = join(intermediate_results_dir, "/timepoints/gene_exp_profile.tp{i}"),
        deg_profile = join(intermediate_results_dir, "/timepoints/deg_profile.tp{i}")
    
    params:
        nThreads = config["parameters"]["nThreads"],
        reg = config["parameters"]["reg_kernel"],
        n_comp = config["parameters"]["nComp"],
        norm = "-normalize" if config["parameters"]["norm"] else "",
        edge_thr = config["parameters"]["edge_thr"]
    
    output:
        join(intermediate_results_dir, config["parameters"]["prefix"] + ".kccagrn.tp{i}.pkl")
    
    shell:
        "python src/2_GRN_inference.py {input.template_network} \
        {input.gene_exp_profile} {input.deg_profile} \
        {output} -TFli {input.tf_list} -TFmodule {input.dict_cluster_TF} \
        -TGmodule {input.dict_cluster_TG} -nComp {params.n_comp} -thr {params.edge_thr} \
        -reg {params.reg} -nThreads {params.nThreads} {params.norm}"
 
rule convert_cluster_to_dict:
    input:
        tsv = join(intermediate_results_dir, "/TFTG_cluster/community_{type}.tsv")
    output:
        join(intermediate_results_dir, "/TFTG_cluster/community_{type}.pkl")
    run:
        import pandas as pd
        import pickle
        df_TG_community_groupBy = pd.read_csv(input.tsv,sep='\t')
        with open(output,'wb') as f:
            pickle.dump(df_TG_community_groupBy.groupby('community')['gene'].apply(list).to_dict(),f)

rule cluster_network:
    input:
        inst_nwk_TF = join(intermediate_results_dir, "TFTG_cluster/condition_specific_nwk_{type}.tsv")
    output:
        cluster = join(intermediate_results_dir, "/TFTG_cluster/community_{type}.tsv"),
        edgelist = join(intermediate_results_dir, "/TFTG_cluster/communityEdgeList_{type}.tsv")
    shell:
        "Rscript src/1_clustering.R {input.inst_nwk_TF} {output.cluster} {output.edgelist}"

rule instantiate_network:
    input:
        template_nwk = join(intermediate_results_dir, "network_ppi/template_nwk_{type}.tsv"),
        gene_exp = config["inputs"]["exp_file"]
    output:
        inst_nwk = join(intermediate_results_dir, "TFTG_cluster/condition_specific_nwk_{type}.tsv")
    params:
        corrCut = config["parameters"]["corr_thr"],
        n_threads = config["parameters"]["nThreads"]
    shell:
        "python src/0_instantiate_nwk.py {input.template_nwk} \
        {input.gene_exp} -corrCut {params.corrCut} -o {output.inst_nwk} -nThreads {params.n_threads}"  

rule split_network:
    input:
        TF_li = config["inputs"]["tf_li_file"],
        PPI = config["inputs"]["ppi_file"]
    output:
        ppi_tf = join(intermediate_results_dir, "network_ppi/template_nwk_TF.tsv"),
        ppi_tg = join(intermediate_results_dir, "network_ppi/template_nwk_TG.tsv")

    run:
        import pandas as pd
        import numpy as np

        df_TF = pd.read_csv(input.TF_li,sep='\t')
        df_nwk = pd.read_csv(input.PPI,sep='\t') 
        
        df_nwk_TF = df_nwk.loc[lambda x:np.logical_and(x.protein1.isin(df_TF.Name),x.protein2.isin(df_TF.Name))]
        df_nwk_TF.to_csv(output.PPI_TF,sep='\t',index=False)

        df_nwk_TG = df_nwk.loc[lambda x:np.logical_and(~x.protein1.isin(df_TF.Name),~x.protein2.isin(df_TF.Name))]
        df_nwk_TG.to_csv(output.PPI_TG,sep='\t',index=False)

rule split_deg:
    input:
        deg_binary = config["inputs"]["seed_file"]
    params:
        n_timepoints = config["parameters"]["nTimePoints"]
    output:
        expand(join(intermediate_results_dir, "/timepoints/deg_profile.tp{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])
    run:
        import pandas as pd
        df = pd.read_csv(input.deg_binary, sep = "\t", index_col = 0)

        for i in range(0, n_timepoints):
            df[df.columns[i]].to_csv(join(intermediate_results_dir, "/timepoints/deg_profile.tp") + str(i+1), sep = "\t")

rule split_exp:
    input:
        exp_file = config["inputs"]["exp_file"]
    params:
        n_timepoints = config["parameters"]["nTimePoints"],
        n_replicates = config["parameters"]["nReplicates"]
    output:
        expand(join(intermediate_results_dir, "/timepoints/gene_exp_profile.tp{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])
    run:
        import pandas as pd
        
        df = pd.read_csv(input.exp_file, sep = "\t", index_col = 0)
        for i, sample in enumeerate(range(params.n_replicates - 1, params.n_timepoints, params.n_replicates)):
            df[df.columns[sample:params.n_replicates]].to_csv(join(intermediate_results_dir, "/timepoints/gene_exp_profile.tp") + str(i+1), sep = "\t")