from os.path import join

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
        "python src/Target_genes.py {input.nwk} {input.DEGli} {params.tf_li_file} {input.tfranktrim} {params.output_dir} {wildcards.i}"

rule extract_optimal_tf_list:
    input:
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".tp{i}.tsv"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])

    params:
        tf_li_file = config["inputs"]["tf_li_file"],
        exp_file = config["inputs"]["exp_file"],
        seed_file = config["inputs"]["seed_file"],
        n_threads = config["parameters"]["nThreads"],
        prefix = config["parameters"]["prefix"],
        temp = intermediate_results_dir,
        n_sample = config["parameter"]["nTimePoints"]

    output:
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".TF_rank.t{i}.trim"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)]),
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".nwk.t{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)]),
        expand(join(intermediate_results_dir, config["parameters"]["prefix"] + ".DEGli.t{i}"), \
            i = [tp for tp in range(1, config["parameters"]["nTimePoints"] + 1)])

    shell:
        "python src/new_TF_adding_NP_noCtrl.py -TFliFile {params.tf_li_file} -nwkFile {params.temp}/{params.prefix}.tp -expFile {params.exp_file} -binFile {params.seed_file} -N {params.n_sample} -p {params.n_threads} -cond {params.prefix} -outD {params.temp}"
    

rule construct_weighted_template_network:
    input:
        config["inputs"]["grn_prefix"] + "{i}.tsv"

    params:
        exp_file = config["inputs"]["exp_file"],
        n_threads = config["parameters"]["nThreads"],
        prefix = config["parameters"]["prefix"],

    output:
        join(intermediate_results_dir, "{prefix}.tp{i}.tsv")
    
    shell:
        "python src/new_network_weight.py -nwk {input} -exp {params.exp_file} -o {output} -p {params.n_threads}"