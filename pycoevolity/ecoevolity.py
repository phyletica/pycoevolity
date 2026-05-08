#!/usr/bin/env python

import os
import sys
import math
import re
import subprocess
import multiprocessing
import operator
import glob
import pandas as pd

import pycoevolity
import pycoevolity.ecoevolity_config as eco_config


def get_ecoevolity_dir(dir_to_check = None):
    return pycoevolity.interop.get_exe_dir(
        exe_name = "ecoevolity",
        dir_to_check = dir_to_check,
    )

def parse_info_from_output(output_str):
    run_time_pattern = re.compile(
        r'^\s*runtime:\s+(?P<run_time>\d+)\s+seconds\.\s*$',
        re.IGNORECASE,
    )
    summary_pattern = re.compile(
        r'^\s*Summary\s+of\s+data\s+from\s+(?P<num_comparisons>\d+)\s+comparisons:\s*$',
        re.IGNORECASE,
    )
    num_var_sites_pattern = re.compile(
        r'^\s*Number\s+of\s+variable\s+sites:\s+(?P<num_var_sites>\d+)\s*$',
        re.IGNORECASE,
    )
    run_time = None
    num_comparisons = None
    num_var_sites = []
    for l in output_str.splitlines():
        line = l.strip()
        m = summary_pattern.match(line)
        if m:
            num_comparisons = int(m.group("num_comparisons"))
            continue
        m = num_var_sites_pattern.match(line)
        if m:
            num_var_sites.append(int(m.group("num_var_sites")))
            continue
        m = run_time_pattern.match(line)
        if m:
            run_time = int(m.group("run_time"))
    if num_comparisons is None:
        raise Exception(
            f"Could not find number of comparisons in this output:\n{output_str}\n"
        )
    if not num_var_sites:
        raise Exception(
            f"Could not find number of variable sites in this output:\n{output_str}\n"
        )
    if run_time is None:
        raise Exception(
            f"Could not find runtime in this output:\n{output_str}\n"
        )
    if num_comparisons != len(num_var_sites):
        raise Exception(
            f"Expected number of variable sites for {num_comparisons} pairs, "
            f"but found {len(num_var_sites)} in this output:\n{output_str}\n"
        )
    return run_time, num_var_sites

def clean_up_ecoevolity_output(state_log_path):
    operator_log_path = state_log_path.replace("state", "operator")
    if os.path.exists(state_log_path):
        os.remove(state_log_path)
    if os.path.exists(operator_log_path):
        os.remove(operator_log_path)
    return

def run_ecoevolity(
    config_path,
    seed,
    output_dir,
    working_dir = None,
    eco_exe_dir = None,
    ignore_data = False,
    relax_constant_sites = False,
    relax_missing_sites = False,
    relax_triallelic_sites = False,
    timeout = None,
    max_num_attempts = 1,
    compress_log_path = True,
    extra_returns = [],
):
    out_dir = output_dir
    if working_dir:
        out_dir = os.path.abspath(
            os.path.join(
                working_dir,
                output_dir,
            )
        )
    if (not os.path.exists(out_dir)) or (not os.path.isdir(out_dir)):
        raise Exception(
            f"Output directory is not an existing directory: {out_dir}"
        )
    eco_exe = "ecoevolity"
    if eco_exe_dir:
        eco_exe = os.path.join(eco_exe_dir, eco_exe)

    config_file_name = os.path.basename(config_path)
    config_name = os.path.splitext(config_file_name)[0]

    prefix = os.path.join(
        out_dir,
        f"run-{seed}-",
    )

    state_log_path = f"{prefix}{config_name}-state-run-1.log"

    cmd = [
        eco_exe,
        f"--seed={seed}",
        f"--prefix={prefix}",
    ]
    if ignore_data:
        cmd.append("--ignore-data")
    if relax_constant_sites:
        cmd.append("--relax-constant-sites")
    if relax_missing_sites:
        cmd.append("--relax-missing-sites")
    if relax_triallelic_sites:
        cmd.append("--relax-triallelic-sites")
    cmd.append(config_path)

    gz_state_log_path = state_log_path
    if compress_log_path:
        gz_state_log_path = f"{state_log_path}.gz"

    if os.path.exists(gz_state_log_path):
        raise Exception(
            f"Output log path already exists: '{gz_state_log_path}'\n"
        )
    
    if compress_log_path and (os.path.exists(state_log_path)):
        # We have the log path but not the gzipped log path, so we'll assume
        # the prior analysis didn't finish and needs to be re-run
        clean_up_ecoevolity_output(state_log_path)

    for attempt_idx in range(max_num_attempts):
        try:
            result = pycoevolity.interop.run_cmd(
                cmd,
                timeout = timeout,
                cwd = working_dir,
            )
            # `check_returncode` will raise CalledProcessError if return code
            # is non-zero This is likely redundant given we use `check = True`
            # in `run_cmd`, but it doesn't hurt to double check
            result.check_returncode()
            break
        except Exception as e:
            if (attempt_idx + 1) < max_num_attempts:
                sys.stderr.write(
                    f"Attempt {attempt_idx + 2} for command:\n\t{cmd}\n"
                )
                clean_up_ecoevolity_output(state_log_path)
                continue
            raise e

    run_time = None
    try:
        run_time, num_var_sites = parse_info_from_output(result.stdout)
    except Exception as e:
        raise Exception(
            f"ERROR: could not parse stdout from ecoevolity run; "
            f"here is the run's stdout and stderr:\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}\n"
        )

    if compress_log_path:
        pycoevolity.fileio.compress_file(state_log_path, gz_state_log_path)
        os.remove(state_log_path)
        state_log_path = gz_state_log_path

    if working_dir:
        state_log_path = os.path.relpath(state_log_path, working_dir)

    return run_time, num_var_sites, state_log_path, *extra_returns

def get_comparison_map(
    simcoevolity_config_comparisons,
    inference_config_comparisons,
):
    if len(simcoevolity_config_comparisons) != len(inference_config_comparisons):
        raise Exception(
            f"simcoevolity config ({len(simcoevolity_config_comparisons)}) "
            f"and inference config ({len(inference_config_comparisons)}) "
            f"have different number of comparisons."
        )
    to_simco_map = {}
    for i, i_comp in enumerate(inference_config_comparisons):
        found = False
        i_path = os.path.basename(i_comp["comparison"]["path"])
        for j, s_comp in enumerate(simcoevolity_config_comparisons):
            s_path = s_comp["comparison"]["path"]
            if s_path.endswith(i_path):
                if found is True:
                    raise Exception(
                        f"Multiple comparisons' paths in the simcoevolity "
                        f"match path in inference config ({i_path})."
                    )
                found = True
                assert not i in to_simco_map
                to_simco_map[i] = j
        if not found:
            raise Exception(
                f"Comparison path {i_path} in inference config does not match "
                f"any comparison paths in the simcoevolity config."
            )
    return to_simco_map

def create_sim_configs(sim_config_path, inference_config_paths):
    s_conf = eco_config.get_yaml_config(sim_config_path)
    s_conf_path_prefix = os.path.splitext(sim_config_path)[0]
    sim_infer_config_paths = []
    for infer_conf_path in inference_config_paths:
        i_conf = eco_config.get_yaml_config(infer_conf_path)
        i_conf_name = os.path.splitext(os.path.basename(infer_conf_path))[0]
        i_out_path = f"{s_conf_path_prefix}-{i_conf_name}.yml"
        inf_to_sim_indices = get_comparison_map(
            simcoevolity_config_comparisons = s_conf["comparisons"],
            inference_config_comparisons = i_conf["comparisons"],
        )
        for i_idx, s_idx in inf_to_sim_indices.items():
            i_conf["comparisons"][i_idx]["comparison"]["path"]                      = s_conf["comparisons"][s_idx]["comparison"]["path"]
            # Ploidy is really a modeling choice
            # i_conf["comparisons"][i_idx]["ploidy"]                    = s_conf["comparisons"][s_idx]["ploidy"]
            i_conf["comparisons"][i_idx]["comparison"]["genotypes_are_diploid"]     = s_conf["comparisons"][s_idx]["comparison"]["genotypes_are_diploid"]
            i_conf["comparisons"][i_idx]["comparison"]["markers_are_dominant"]      = s_conf["comparisons"][s_idx]["comparison"]["markers_are_dominant"]
            i_conf["comparisons"][i_idx]["comparison"]["population_name_delimiter"] = s_conf["comparisons"][s_idx]["comparison"]["population_name_delimiter"]
            i_conf["comparisons"][i_idx]["comparison"]["population_name_is_prefix"] = s_conf["comparisons"][s_idx]["comparison"]["population_name_is_prefix"]
            i_conf["comparisons"][i_idx]["comparison"]["constant_sites_removed"]    = s_conf["comparisons"][s_idx]["comparison"]["constant_sites_removed"]
        eco_config.write_yaml_config(i_conf, i_out_path)
        sim_infer_config_paths.append((i_out_path, infer_conf_path))
    return sim_infer_config_paths

def run_simcoevolity(
    sim_config_path,
    infer_config_paths,
    seed,
    output_dir,
    number_of_replicates,
    working_dir = None,
    eco_exe_dir = None,
    singleton_sample_prob = None,
    locus_size = None,
    max_one_variable_site_per_locus = False,
    charsets = False,
    relax_constant_sites = False,
    relax_missing_sites = False,
    relax_triallelic_sites = False,
    output_nexus = False,
    compress_true_values_paths = True,
    extra_returns = [],
):
    out_dir = output_dir
    if working_dir:
        out_dir = os.path.abspath(
            os.path.join(
                working_dir,
                output_dir,
            )
        )
    if (not os.path.exists(out_dir)) or (not os.path.isdir(out_dir)):
        raise Exception(
            f"Output directory is not an existing directory: {out_dir}"
        )

    eco_exe = "simcoevolity"
    if eco_exe_dir:
        eco_exe = os.path.join(eco_exe_dir, eco_exe)

    prefix = f"seed-{seed}-"

    cmd = [
        eco_exe,
        f"--seed={seed}",
        f"--number-of-replicates={number_of_replicates}",
        f"--output-directory={output_dir}",
        f"--prefix={prefix}",
    ]
    if not singleton_sample_prob is None:
        cmd.append(f"--singleton-sample-probability={singleton_sample_prob}")
    if not locus_size is None:
        cmd.append(f"--locus-size={locus_size}")
    if max_one_variable_site_per_locus:
        cmd.append("--max-one-variable-site-per-locus")
    if charsets:
        cmd.append("--charsets")
    if relax_constant_sites:
        cmd.append("--relax-constant-sites")
    if relax_missing_sites:
        cmd.append("--relax-missing-sites")
    if relax_triallelic_sites:
        cmd.append("--relax-triallelic-sites")
    if output_nexus:
        cmd.append("--nexus")
    cmd.append(sim_config_path)

    full_prefix = os.path.join(
        out_dir,
        prefix,
    )
    sim_model_out_path = f"{full_prefix}simcoevolity-model-used-for-sims.yml"

    if os.path.exists(sim_model_out_path):
        raise Exception(
            f"Output path already exists: '{sim_model_out_path}'\n"
        )

    result = pycoevolity.interop.run_cmd(cmd, cwd = working_dir)
    if result.returncode != 0:
        raise Exception(
            f"ERROR: simcoevolity run returned non-zero exit code "
            f"'{result.returncode}'; here is the stderr:\n"
            f"{result.stderr}\n"
        )
    run_time, num_var_sites = parse_info_from_output(result.stderr)

    true_val_paths = glob.glob(
        f"{full_prefix}simcoevolity-sim-[0-9]*-true-values.txt"
    )
    assert len(true_val_paths) == number_of_replicates
    sim_numbers = sorted([p.split("-")[-3] for p in true_val_paths])
    
    if working_dir:
        infer_config_paths = [
            os.path.abspath(
                os.path.join(
                    working_dir,
                    p,
                )
            ) for p in infer_config_paths
        ]
    config_paths = []
    true_conf_paths = []
    for sim_num in sim_numbers:
        true_val_path = f"{full_prefix}simcoevolity-sim-{sim_num}-true-values.txt"
        assert os.path.isfile(true_val_path)
        if compress_true_values_paths:
            gz_path = f"{true_val_path}.gz"
            pycoevolity.fileio.compress_file(true_val_path, gz_path)
            os.remove(true_val_path)
            true_val_path = gz_path
        conf_path = f"{full_prefix}simcoevolity-sim-{sim_num}-config.yml"
        assert os.path.isfile(conf_path)
        config_paths.append(conf_path)
        sim_infer_config_paths = create_sim_configs(
            sim_config_path = conf_path,
            inference_config_paths = infer_config_paths,
        )
        if working_dir:
            true_val_path = os.path.relpath(true_val_path, working_dir)
            sim_infer_config_paths = [
                (
                    os.path.relpath(c1, working_dir),
                    os.path.relpath(c2, working_dir),
                ) for c1, c2 in sim_infer_config_paths
            ]
        true_conf_paths.append(
            (true_val_path, tuple(sim_infer_config_paths))
        )

    assert len(config_paths) == number_of_replicates
    assert len(true_conf_paths) == number_of_replicates

    for path in config_paths:
        os.remove(path)

    return run_time, tuple(true_conf_paths), *extra_returns

def collect_prior_samples(
    seeds,
    config_path,
    output_dir,
    number_of_procs = 4,
    working_dir = None,
    eco_exe_dir = None,
    timeout = 300,
    max_num_attempts = 2,
):
    log_paths = []
    with multiprocessing.Pool(number_of_procs) as pool:
        workers = [
            pool.apply_async(
                run_ecoevolity,
                args = (
                    config_path,
                    seed,
                    output_dir,
                    working_dir,
                    eco_exe_dir,
                    True,   # ignore_data
                    False,  # relax_constant_sites
                    False,  # relax_missing_sites
                    False,  # relax_triallelic_sites
                    timeout,
                    False,  # compress_log_path
                    max_num_attempts,
                )
            )
            for seed in seeds
        ]
        sys.stdout.write(
            f"Loaded {len(workers)} ecoevolity workers for {number_of_procs} processors\n"
        )
        for run_time, num_var_sites, state_log_path in (w.get() for w in workers):
            log_paths.append(state_log_path)
    return log_paths

def run_sumcoevolity(
    config_path,
    input_state_log_paths,
    seed,
    output_dir,
    num_prior_draws = 1000000,
    working_dir = None,
    eco_exe_dir = None,
    burnin = 0,
    compress_results = True,
    extra_returns = [],
):
    abs_out_dir = output_dir
    if working_dir:
        abs_out_dir = os.path.abspath(
            os.path.join(
                working_dir,
                output_dir,
            )
        )
    if (not os.path.exists(abs_out_dir)) or (not os.path.isdir(abs_out_dir)):
        raise Exception(
            f"Output directory is not an existing directory: {abs_out_dir}"
        )

    eco_exe = "sumcoevolity"
    if eco_exe_dir:
        eco_exe = os.path.join(eco_exe_dir, eco_exe)

    config_file_name = os.path.basename(config_path)
    config_name = os.path.splitext(config_file_name)[0]

    prefix = os.path.join(
        abs_out_dir,
        f"{config_name}-seed-{seed}-n-{num_prior_draws}-",
    )

    nevents_results_path = f"{prefix}sumcoevolity-results-nevents.txt"
    model_results_path = f"{prefix}sumcoevolity-results-model.txt"

    input_log_paths = []
    for log_path in input_state_log_paths:
        abs_log_path = log_path
        if working_dir:
            abs_log_path = os.path.abspath(
                os.path.join(
                    working_dir,
                    log_path,
                )
            )
        if abs_log_path.endswith(".gz"):
            decomp_path = abs_log_path[:-3]
            pycoevolity.fileio.decompress_file(abs_log_path, decomp_path)
            input_log_paths.append(decomp_path)
        else:
            input_log_paths.append(abs_log_path)

    cmd = [
        eco_exe,
        f"--seed={seed}",
        f"--prefix={prefix}",
        f"--config={config_path}",
        f"--number-of-samples={num_prior_draws}",
        f"--burnin={burnin}",
        *input_log_paths,
    ]

    if os.path.exists(nevents_results_path):
        raise Exception(
            f"Output results path already exists: '{nevents_results_path}'\n"
        )
    if os.path.exists(model_results_path):
        raise Exception(
            f"Output results path already exists: '{model_results_path}'\n"
        )

    result = pycoevolity.interop.run_cmd(cmd, cwd = working_dir)

    if compress_results:
        gz_nevents_results_path = f"{nevents_results_path}.gz"
        pycoevolity.fileio.compress_file(nevents_results_path, gz_nevents_results_path)
        os.remove(nevents_results_path)
        nevents_results_path = gz_nevents_results_path

        gz_model_results_path = f"{model_results_path}.gz"
        pycoevolity.fileio.compress_file(model_results_path, gz_model_results_path)
        os.remove(model_results_path)
        model_results_path = gz_model_results_path

    if working_dir:
        nevents_results_path = os.path.relpath(nevents_results_path, working_dir)
        model_results_path = os.path.relpath(model_results_path, working_dir)

    return result, nevents_results_path, model_results_path, *extra_returns

def prepare_simulations(
    rng,
    sim_configs,
    infer_configs,
    output_dir,
    working_dir = None,
    eco_exe_dir = None,
    number_of_sims = 100,
    number_of_procs = 4,
    number_of_chains = 2,
    singleton_sample_prob = None,
    locus_size = None,
    max_one_variable_site_per_locus = False,
    charsets = False,
    relax_constant_sites = False,
    relax_missing_sites = False,
    relax_triallelic_sites = False,
    output_nexus = False,
):
    if number_of_procs > number_of_sims:
        number_of_procs = number_of_sims
    num_sims_per_proc = math.floor(number_of_sims / number_of_procs)
    remainder_sims = number_of_sims - (num_sims_per_proc * number_of_procs)

    num_sims_args = [num_sims_per_proc for _ in range(number_of_procs)]
    num_sims_args[-1] += remainder_sims

    workers = []
    seeds = []
    results = {}
    with multiprocessing.Pool(number_of_procs) as pool:
        for sim_config in sim_configs:
            results[sim_config] = {}
            for num_reps in num_sims_args:
                seed = pycoevolity.rng_utils.get_safe_seed(rng)
                seeds.append(seed)
                results[sim_config][seed] = []
                workers.append(
                    pool.apply_async(
                        run_simcoevolity,
                        args = (
                            sim_config,
                            infer_configs,
                            seed,
                            output_dir,
                            num_reps,
                            working_dir,
                            eco_exe_dir,
                            singleton_sample_prob,
                            locus_size,
                            max_one_variable_site_per_locus,
                            charsets,
                            relax_constant_sites,
                            relax_missing_sites,
                            relax_triallelic_sites,
                            output_nexus,
                            True,  # compress_true_values_paths
                            (sim_config, seed), # extra_returns
                        )
                    )
                )
        num_workers = len(workers)
        sys.stdout.write(
            f"Loaded {num_workers} simcoevolity workers for {number_of_procs} processors\n"
        )
        reporting_freq = max(num_workers // 10, 1)
        count = 0
        for run_time, true_val_config_paths, sim_config, seed in (w.get() for w in workers):
            results[sim_config][seed].extend(true_val_config_paths)
            count += 1
            if (count < num_workers) and (count % reporting_freq == 0):
                sys.stdout.write(
                    f"{count} of {num_workers} simcoevolity workers finished\n"
                )
        sys.stdout.write(
            f"{count} of {num_workers} simcoevolity workers finished\n"
        )
    ret = {}
    for sim_config in sim_configs:
        ret[sim_config] = {}
        for seed in seeds:
            true_val_config_paths = results[sim_config][seed]
            for true_val_path, inf_config_paths in true_val_config_paths:
                ret[sim_config][true_val_path] = {}
                ret[sim_config][true_val_path]["analyses"] = {}
                ret[sim_config][true_val_path]["numbers_of_variable_sites"] = None
                for rep_inf_conf, orig_inf_conf in inf_config_paths:
                    assert not orig_inf_conf in ret[sim_config][true_val_path]["analyses"]
                    ret[sim_config][true_val_path]["analyses"][orig_inf_conf] = {}
                    ret[sim_config][true_val_path]["analyses"][orig_inf_conf]["replicate_inference_config"] = rep_inf_conf
                    ret[sim_config][true_val_path]["analyses"][orig_inf_conf]["chains"] = {}
                    for chain_idx in range(number_of_chains):
                        seed = pycoevolity.rng_utils.get_safe_seed(rng)
                        assert not seed in ret[sim_config][true_val_path]["analyses"][orig_inf_conf]["chains"]
                        ret[sim_config][true_val_path]["analyses"][orig_inf_conf]["chains"][seed] = {
                            "run_time" : None,
                            "state_log_path" : None,
                        }
                    ret[sim_config][true_val_path]["analyses"][orig_inf_conf]["sumcoevolity"] = {
                        "seed" : pycoevolity.rng_utils.get_safe_seed(rng),
                        "nevents_summary_path" : None,
                        "model_summary_path" : None,
                    }
    return ret

def run_analyses_on_sims(
    sim_data,
    working_dir = None,
    eco_exe_dir = None,
    number_of_procs = 4,
    relax_constant_sites = False,
    relax_missing_sites = False,
    relax_triallelic_sites = False,
    timeout = None,
    max_num_attempts = 3,
    output_dir = None,
):
    total_num_runs = 0
    out_dir = None
    for sim_config, rep_data in sim_data.items():
        for true_val_path, results_info in rep_data.items():
            if out_dir is None:
                if working_dir:
                    out_dir = os.path.abspath(
                        os.path.join(
                            working_dir,
                            os.path.dirname(true_val_path),
                        )
                    )
                else:
                    out_dir = os.path.abspath(
                        os.path.dirname(true_val_path),
                    )
            for inf_config, inf_data in results_info["analyses"].items():
                for seed, chain_info in inf_data["chains"].items():
                    if not chain_info["state_log_path"]:
                        total_num_runs += 1

    if number_of_procs > total_num_runs:
        number_of_procs = total_num_runs

    if not output_dir:
        output_dir = out_dir

    count = 0
    if total_num_runs > 0:
        workers = []
        with multiprocessing.Pool(number_of_procs) as pool:
            for sim_config, rep_data in sim_data.items():
                for true_vals_path, results_info in rep_data.items():
                    for orig_inf_config, inf_data in results_info["analyses"].items():
                        rep_inf_config = inf_data["replicate_inference_config"]
                        for seed, chain_info in inf_data["chains"].items():
                            if not chain_info["state_log_path"]:
                                workers.append(
                                    pool.apply_async(
                                        run_ecoevolity,
                                        args = (
                                            rep_inf_config,
                                            seed,
                                            output_dir,
                                            working_dir,
                                            eco_exe_dir,
                                            False,  # ignore_data
                                            False,  # relax_constant_sites
                                            False,  # relax_missing_sites
                                            False,  # relax_triallelic_sites
                                            timeout,
                                            max_num_attempts,
                                            True,   # compress_log_path
                                            (sim_config, true_vals_path, orig_inf_config, seed), # extra_returns
                                        )
                                    )
                                )
            num_workers = len(workers)
            assert num_workers == total_num_runs
            sys.stdout.write(
                f"Loaded {num_workers} ecoevolity workers for {number_of_procs} processors\n"
            )
            reporting_freq = max(num_workers // 10, 1)
            num_var_sites = {}
            for run_time, n_var_sites, state_log_path, sim_config, true_path, inf_config, seed in (w.get() for w in workers):
                if sim_data[sim_config][true_path]["numbers_of_variable_sites"]:
                    assert sim_data[sim_config][true_path]["numbers_of_variable_sites"] == n_var_sites
                else:
                    sim_data[sim_config][true_path]["numbers_of_variable_sites"] = n_var_sites
                sim_data[sim_config][true_path]["analyses"][inf_config]["chains"][seed]["run_time"] = run_time
                assert not sim_data[sim_config][true_path]["analyses"][inf_config]["chains"][seed]["state_log_path"]
                sim_data[sim_config][true_path]["analyses"][inf_config]["chains"][seed]["state_log_path"] = state_log_path
                count += 1
                if (count < num_workers) and (count % reporting_freq == 0):
                    sys.stdout.write(
                        f"{count} of {num_workers} ecoevolity workers finished\n"
                    )
            sys.stdout.write(
                f"{count} of {num_workers} ecoevolity workers finished\n"
            )
    assert count == total_num_runs
    return count

def add_sumcoevolity_to_results(
    sim_data,
    working_dir = None,
    eco_exe_dir = None,
    output_dir = None,
    num_prior_draws = 1000000,
    burnin = 0,
    number_of_procs = 4,
):
    total_num_runs = 0
    out_dir = None
    for sim_config, rep_data in sim_data.items():
        for true_val_path, results_info in rep_data.items():
            if out_dir is None:
                out_dir = os.path.dirname(true_val_path)
            for inf_config, inf_data in results_info["analyses"].items():
                if not inf_data["sumcoevolity"]["nevents_summary_path"]:
                    total_num_runs += 1

    if number_of_procs > total_num_runs:
        number_of_procs = total_num_runs

    if not output_dir:
        output_dir = out_dir

    count = 0
    if total_num_runs > 0:
        nchains = None
        workers = []
        with multiprocessing.Pool(number_of_procs) as pool:
            for sim_config, rep_data in sim_data.items():
                for true_vals_path, infer_info in rep_data.items():
                    for config_path, results_info in infer_info["analyses"].items():
                        if not results_info["sumcoevolity"]["nevents_summary_path"]:
                            state_log_paths = [results_info["chains"][seed]["state_log_path"] for seed in results_info["chains"]]
                            if nchains is None:
                                nchains = len(state_log_paths)
                            else:
                                assert nchains == len(state_log_paths)
                            seed = results_info["sumcoevolity"]["seed"]
                            extra_returns = (sim_config, true_vals_path, config_path)
                            workers.append(
                                pool.apply_async(
                                    run_sumcoevolity,
                                    args = (
                                        config_path,
                                        state_log_paths,
                                        seed,
                                        output_dir,
                                        num_prior_draws,
                                        working_dir,
                                        eco_exe_dir,
                                        burnin,
                                        True,  # compress_results
                                        extra_returns,
                                    )
                                )
                            )
            num_workers = len(workers)
            assert num_workers == total_num_runs
            sys.stdout.write(
                f"Loaded {num_workers} sumcoevolity workers for {number_of_procs} processors\n"
            )
            reporting_freq = max(num_workers // 10, 1)
            for res, nevents_results_path, model_results_path, sim_config, true_vals_path, inf_config in (w.get() for w in workers):
                assert true_vals_path in sim_data[sim_config]
                assert inf_config in sim_data[sim_config][true_vals_path]["analyses"]
                sim_data[sim_config][true_vals_path]["analyses"][inf_config]["sumcoevolity"]["nevents_summary_path"] = nevents_results_path
                sim_data[sim_config][true_vals_path]["analyses"][inf_config]["sumcoevolity"]["model_summary_path"] = model_results_path
                count += 1
                if (count < num_workers) and (count % reporting_freq == 0):
                    sys.stdout.write(
                        f"{count} of {num_workers} sumcoevolity workers finished\n"
                    )
            sys.stdout.write(
                f"{count} of {num_workers} sumcoevolity workers finished\n"
            )
    assert count == total_num_runs
    return count

def get_result_path(rel_path, results_dir):
    return os.path.abspath(os.path.join(results_dir, rel_path))

def parse_true_values(true_values_path):
    true_values = pycoevolity.parsing.get_dict_from_spreadsheets(
        [true_values_path],
        sep = "\t",
        header = None,
    )
    for v in true_values.values():
        assert len(v) == 1
    return true_values

def parse_sim_rep_results(
    sim_id,
    sim_config_name,
    inference_config_name,
    true_values,
    state_log_paths,
    run_times,
    parameter_names,
    numbers_of_variable_sites = None,
    include_time_in_coal_units = True,
    burnin = 0,
    config_labels = None,
):
    assert len(run_times) == len(state_log_paths)
    if not config_labels:
        config_labels = {
            sim_config_name : sim_config_name,
            inference_config_name : inference_config_name,
        }
    nchains = len(state_log_paths)
    post_sample = pycoevolity.posterior.PosteriorSample(
        state_log_paths,
        burnin = burnin,
        include_time_in_coal_units = include_time_in_coal_units,
    )
    results = {
        'simulation_id': sim_id,
        'simulation_config' : config_labels.get(
            sim_config_name, sim_config_name),
        'inference_config': config_labels.get(
            inference_config_name, inference_config_name),
        'mean_run_time' : sum(run_times) / len(run_times),
        'median_run_time' : pycoevolity.stats.median(run_times),
        'min_run_time' : min(run_times),
        'sample_size' : post_sample.number_of_samples,
    }

    assert post_sample.number_of_samples % nchains == 0
    nsamples_per_chain = post_sample.number_of_samples // nchains

    if include_time_in_coal_units:
        for comp_idx in range(post_sample.number_of_comparisons):
            comp_label = post_sample.height_labels[comp_idx]
            ht_key = f"root_height_{comp_label}"
            sz_keys = [f"pop_size_{l}" for l in post_sample.tip_labels[comp_idx]]
            t = float(true_values[ht_key][0])
            pop_sizes = [float(true_values[k][0]) for k in sz_keys]
            n = sum(pop_sizes) / len(pop_sizes)
            t_coal = t / (2.0 * n)
            coal_key = "coal_root_height_{0}".format(comp_label)
            true_values[coal_key] = [t_coal]

    if numbers_of_variable_sites:
        assert len(numbers_of_variable_sites) == post_sample.number_of_comparisons
        for i, n_var_sites in enumerate(numbers_of_variable_sites):
            results[f"n_var_sites_{post_sample.height_labels[i]}"] = n_var_sites
    
    true_model = tuple(int(true_values[h][0]) for h in post_sample.height_index_keys)
    true_model_p = post_sample.get_model_probability(true_model)
    true_model_cred = post_sample.get_model_credibility_level(true_model)
    map_models = post_sample.get_map_models()
    map_model = map_models[0]
    if len(map_models) > 1:
        if true_model in map_models:
            map_model = true_model
    map_model_p = post_sample.get_model_probability(map_model)
    results["true_model"] = ",".join((str(i) for i in true_model))
    results["map_model"] = ",".join((str(i) for i in map_model))
    results["true_model_cred_level"] = true_model_cred
    results["map_model_p"] = map_model_p
    results["true_model_p"] = true_model_p
    model_dist_summary = pycoevolity.stats.get_summary(
        post_sample.distances_from(true_model))
    results["mean_model_distance"] = model_dist_summary["mean"]
    results["median_model_distance"] = model_dist_summary["median"]
    results["std_dev_model_distance"] = math.sqrt(model_dist_summary["variance"])
    results["hpdi_95_lower_model_distance"] = model_dist_summary["hpdi_95"][0]
    results["hpdi_95_upper_model_distance"] = model_dist_summary["hpdi_95"][1]
    results["eti_95_lower_model_distance"] = model_dist_summary["qi_95"][0]
    results["eti_95_upper_model_distance"] = model_dist_summary["qi_95"][1]
    map_model_distances = post_sample.get_map_model_distances_from(true_model)
    if len(map_model_distances) > 1:
        map_model_dist_summary = pycoevolity.stats.get_summary(
                map_model_distances)
        results["mean_map_model_distance"] = map_model_dist_summary["mean"]
        results["median_map_model_distance"] = map_model_dist_summary["median"]
    else:
        results["mean_map_model_distance"] = map_model_distances[0]
        results["median_map_model_distance"] = map_model_distances[0]
    
    true_nevents = int(true_values["number_of_events"][0])
    true_nevents_p = post_sample.get_number_of_events_probability(true_nevents)
    true_nevents_cred = post_sample.get_number_of_events_credibility_level(true_nevents)
    map_numbers_of_events = post_sample.get_map_numbers_of_events()
    map_nevents = map_numbers_of_events[0]
    if len(map_numbers_of_events) > 1:
        if true_nevents in map_numbers_of_events:
            map_nevents = true_nevents
    map_nevents_p = post_sample.get_number_of_events_probability(map_nevents)
    results["true_num_events"] = true_nevents
    results["map_num_events"] = map_nevents
    results["true_num_events_cred_level"] = true_nevents_cred
    results["map_num_events_p"] = map_nevents_p
    results["true_num_events_p"] = true_nevents_p
    nevents_cred_set = []
    cum_prob = 0.0
    for n, p in post_sample.get_numbers_of_events():
        nevents_cred_set.append(n)
        cum_prob += p
        if cum_prob > 0.95:
            break
    hpdi_lower_nevents = min(nevents_cred_set)
    hpdi_upper_nevents = max(nevents_cred_set)
    results["hpdi_95_lower_num_events"] = hpdi_lower_nevents
    results["hpdi_95_upper_num_events"] = hpdi_upper_nevents
    results["map_num_events_distance"] = map_nevents - true_nevents
    results["hpdi_95_lower_num_events_distance"] = hpdi_lower_nevents - true_nevents
    results["hpdi_95_upper_num_events_distance"] = hpdi_upper_nevents - true_nevents
    
    sum_of_abs_mean_error_root_height = 0.0
    sum_of_abs_mean_error_pop_size_root = 0.0
    for parameter in parameter_names:
        true_val = None
        true_val_rank = None
        post_mean = None
        post_median = None
        post_stdev = None
        hpdi_lower = None
        hpdi_upper = None
        eti_lower = None
        eti_upper = None
        ess = None
        ess_sum = None
        psrf = None
        have_true_val = bool(parameter in true_values)
        have_post = bool(parameter in post_sample.parameter_samples)
        if have_true_val:
            true_val = float(true_values[parameter][0])
        if have_post:
            if have_true_val:
                true_val_rank = post_sample.get_rank(parameter, true_val)
            ss = pycoevolity.stats.get_summary(
                    post_sample.parameter_samples[parameter])
            if parameter in post_sample.get_height_keys():
                sum_of_abs_mean_error_root_height += math.fabs(
                    true_val - ss["mean"])
            elif parameter in post_sample.get_ancestral_pop_size_keys():
                sum_of_abs_mean_error_pop_size_root += math.fabs(
                    true_val - ss["mean"])
            ess = pycoevolity.stats.effective_sample_size(
                    post_sample.parameter_samples[parameter])
            ess_sum = 0.0
            samples_by_chain = []
            for i in range(nchains):
                chain_samples = post_sample.parameter_samples[parameter][
                        i * nsamples_per_chain : (i + 1) * nsamples_per_chain]
                assert(len(chain_samples) == nsamples_per_chain)
                ess_sum += pycoevolity.stats.effective_sample_size(chain_samples)
                if nchains > 1:
                    samples_by_chain.append(chain_samples)
            if nchains > 1:
                psrf = pycoevolity.stats.potential_scale_reduction_factor(samples_by_chain)
            post_mean = ss["mean"]
            post_median = ss["median"]
            post_stdev = math.sqrt(ss["variance"])
            hpdi_lower = ss["hpdi_95"][0]
            hpdi_upper = ss["hpdi_95"][1]
            eti_lower = ss["qi_95"][0]
            eti_upper = ss["qi_95"][1]
        if nchains > 1:
            results["psrf_{0}".format(parameter)] = psrf
        results["true_{0}".format(parameter)] = true_val
        results["true_{0}_rank".format(parameter)] = true_val_rank
        results["mean_{0}".format(parameter)] = post_mean
        results["median_{0}".format(parameter)] = post_median
        results["stddev_{0}".format(parameter)] = post_stdev
        results["hpdi_95_lower_{0}".format(parameter)] = hpdi_lower
        results["hpdi_95_upper_{0}".format(parameter)] = hpdi_upper
        results["eti_95_lower_{0}".format(parameter)] = eti_lower
        results["eti_95_upper_{0}".format(parameter)] = eti_upper
        results["ess_{0}".format(parameter)] = ess
        results["ess_sum_{0}".format(parameter)] = ess_sum
    results["sum_of_abs_mean_error_root_height"] = sum_of_abs_mean_error_root_height
    results["sum_of_abs_mean_error_pop_size_root"] = sum_of_abs_mean_error_pop_size_root
    return results

def get_free_parameter_labels(
    analyses_dict,
    results_dir,
    include_time_in_coal_units = True,
):
    parameters = set()
    time_parameters = None
    for analysis_config, analysis_results in analyses_dict.items():
        chains = analysis_results["chains"]
        log_paths = [chains[seed]["state_log_path"] for seed in chains]
        log_paths = [get_result_path(p, results_dir) for p in log_paths]
        post_sample = pycoevolity.posterior.PosteriorSample(
            log_paths,
            burnin = 1,
            include_time_in_coal_units = include_time_in_coal_units,
        )
        time_params = sorted(post_sample.get_height_keys())
        if not time_parameters:
            time_parameters = time_params
            parameters.update(time_params)
        else:
            if time_params != time_parameters:
                raise Exception(
                    "Time (root_height) parameters do not match among analysis "
                    "outputs from different configs."
                )
        for comp_idx in range(post_sample.number_of_comparisons):
            tip_labels = post_sample.tip_labels[comp_idx]
            comp_label = post_sample.height_labels[comp_idx]
            tip_pop_size_keys = [f"pop_size_{l}" for l in tip_labels]
            anc_pop_size_key = f"pop_size_root_{comp_label}"
            tip_pop_sizes_last = [post_sample.parameter_samples[k][-1] for k in tip_pop_size_keys]
            anc_pop_size_last = post_sample.parameter_samples[anc_pop_size_key][-1]
            anc_pop_size_first = post_sample.parameter_samples[anc_pop_size_key][0]
            pop_sizes_constrained = math.isclose(
                anc_pop_size_last - tip_pop_sizes_last[0], 0.0, abs_tol=1e-10)
            anc_pop_size_fixed = math.isclose(
                anc_pop_size_last - anc_pop_size_first, 0.0, abs_tol=1e-10)
            if not anc_pop_size_fixed:
                parameters.add(anc_pop_size_key)
            if not pop_sizes_constrained:
                tip_pop_sizes_first = [post_sample.parameter_samples[k][0] for k in tip_pop_size_keys]
                tip_size_0_fixed = math.isclose(
                    tip_pop_sizes_last[0] - tip_pop_sizes_first[0], 0.0, abs_tol=1e-10)
                tip_size_1_fixed = math.isclose(
                    tip_pop_sizes_last[1] - tip_pop_sizes_first[1], 0.0, abs_tol=1e-10)
                if not tip_size_0_fixed:
                    parameters.add(tip_pop_size_keys[0])
                if not tip_size_1_fixed:
                    parameters.add(tip_pop_size_keys[1])

        for param in post_sample.parameter_samples.keys():
            if (
                param.startswith("root_height_")
                or param.startswith("pop_size_")
                or (param == "model")
                or (param == "number_of_events")
            ):
                continue
            if not param in parameters:
                first_val = post_sample.parameter_samples[param][0]
                last_val = post_sample.parameter_samples[param][-1]
                param_fixed = math.isclose(
                    first_val - last_val, 0.0, abs_tol = 1e-10)
                if not param_fixed:
                    parameters.add(param)
    return parameters

def parse_sim_id(true_values_path):
    true_values_file_pattern_str = (
        r"^.*seed-(?P<seed>\d+)-simcoevolity-sim-(?P<sim_num>\d+)-true-values\.txt.*$"
    )
    true_values_file_pattern = re.compile(true_values_file_pattern_str)
    m = true_values_file_pattern.match(true_values_path)
    if not m:
        raise Exception(
            f"Unexpected true values file name: {true_values_path}"
        )
    seed = m.group("seed")
    sim_num = m.group("sim_num")
    return f"{seed}-{sim_num}"

def parse_sim_results(
    results_path,
    config_labels = None,
    include_time_in_coal_units = True,
    number_of_procs = 1,
):
    results_dir = os.path.dirname(results_path)
    results = pycoevolity.fileio.load_json(results_path)
    sim_configs = results["simulation_configs"]
    if not pycoevolity.fileio.file_names_are_unique(sim_configs):
        raise Exception(
            "Simulation config file names are not unique"
        )
    inf_configs = results["inference_configs"]
    if not pycoevolity.fileio.file_names_are_unique(inf_configs):
        raise Exception(
            "Inference config file names are not unique"
        )
    nchains = results["number_of_chains"]
    burnin = results["burnin"]
    parameter_names = None

    workers = []
    with multiprocessing.Pool(number_of_procs) as pool:
        for sim_conf, sims in results["simulations"].items():
            assert sim_conf in sim_configs
            for true_vals_path, rep_data in sims.items():
                true_vals_path = get_result_path(true_vals_path, results_dir)
                true_values = parse_true_values(true_vals_path)
                numbers_of_variable_sites = rep_data[
                    "numbers_of_variable_sites"]
                sim_id = parse_sim_id(os.path.basename(true_vals_path))
                if not parameter_names:
                    parameter_names = get_free_parameter_labels(
                        rep_data["analyses"],
                        results_dir,
                        include_time_in_coal_units = include_time_in_coal_units,
                    )
                for analysis_conf, analysis_results in rep_data["analyses"].items():
                    assert analysis_conf in inf_configs
                    chains = analysis_results["chains"]
                    assert len(chains) == nchains
                    run_times = [chains[seed]["run_time"] for seed in chains]
                    log_paths = [chains[seed]["state_log_path"] for seed in chains]
                    log_paths = [get_result_path(p, results_dir) for p in log_paths]
                    workers.append(
                        pool.apply_async(
                            parse_sim_rep_results,
                            args = (
                                sim_id,
                                os.path.basename(sim_conf),      # sim_config_name
                                os.path.basename(analysis_conf), # inference_config_name
                                true_values,
                                log_paths,
                                run_times,
                                parameter_names,
                                numbers_of_variable_sites,
                                include_time_in_coal_units,
                                burnin,
                                config_labels,
                            )
                        )
                    )
        num_workers = len(workers)
        sys.stdout.write(
            f"Loaded {num_workers} result parsing workers for {number_of_procs} processors\n"
        )
        reporting_freq = 100
        if num_workers < 500:
            reporting_freq = 10
        results = []
        for i, result_dict in enumerate(w.get() for w in workers):
            results.append(result_dict)
            if ((i + 1) < num_workers) and ((i + 1) % reporting_freq == 0):
                sys.stdout.write(
                    f"{i + 1} of {num_workers} result parsing workers finished\n"
                )
        sys.stdout.write(
            f"{i + 1} of {num_workers} result parsing workers finished\n"
        )
    results.sort(key=operator.itemgetter(
        'simulation_config',
        'simulation_id',
        'inference_config',
    ))
    df = pd.DataFrame(results)
    return df

def append_results(previous_results, new_results):
    for sim_config, sim_reps in new_results.items():
        for true_vals_path, infer_info in sim_reps.items():
            assert not true_vals_path in previous_results["simulations"][sim_config]
            previous_results["simulations"][sim_config][true_vals_path] = infer_info
