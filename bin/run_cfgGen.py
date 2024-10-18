import argparse
import array
import os
import random

import ROOT
from pyamptools import atiSetup
from pyamptools.utility.general import converter, example_vps_names, example_zlm_names, load_yaml, vps_amp_name, zlm_amp_name

from types import SimpleNamespace

import sys
sys.path.append("/work/halld/home/ddarulis/jlab_utils")
from python_scripts.get_vecps_plots import parse_wave_string
from python_scripts.get_vecps_plots import m_names, l_to_num
from python_scripts import make_amptools_bins
from python_scripts import run_amptools
from python_scripts import plot_fits
from python_scripts import get_best_lik
############################################################################
# This script generates AmpTools configuration files with knobs/flags to
# append additional information to the generated file
############################################################################

help_header = """#####################################
####	THIS IS A CONFIG FILE	 ####
#####################################
##
##  Blank lines or lines beginning with a "#" are ignored.
##
##  Double colons (::) are treated like a space.
##     This is sometimes useful for grouping (for example,
##     grouping strings like "reaction::sum::amplitudeName")
##
##  All non-comment lines must begin with one of the following keywords.
##
##  (note:  <word> means necessary
##	    (word) means optional)
##
##  include	  <file>
##  define	  <word> (defn1) (defn2) (defn3) ...
##  fit 	  <fitname>
##  keyword	  <keyword> <min arguments> <max arguments>
##  reaction	  <reaction> <particle1> <particle2> (particle3) ...
##  data	  <reaction> <class> (arg1) (arg2) (arg3) ...
##  genmc	  <reaction> <class> (arg1) (arg2) (arg3) ...
##  accmc	  <reaction> <class> (arg1) (arg2) (arg3) ...
##  normintfile   <reaction> <file>
##  sum 	  <reaction> <sum> (sum2) (sum3) ...
##  amplitude	  <reaction> <sum> <amp> <class> (arg1) (arg2) ([par]) ...
##  initialize    <reaction> <sum> <amp> <"events"/"polar"/"cartesian">
##		    <value1> <value2> ("fixed"/"real")
##  scale	  <reaction> <sum> <amp> <value or [parameter]>
##  constrain	  <reaction1> <sum1> <amp1> <reaction2> <sum2> <amp2> ...
##  permute	  <reaction> <sum> <amp> <index1> <index2> ...
##  parameter	  <par> <value> ("fixed"/"bounded"/"gaussian")
##		    (lower/central) (upper/error)
##    DEPRECATED:
##  datafile	  <reaction> <file> (file2) (file3) ...
##  genmcfile	  <reaction> <file> (file2) (file3) ...
##  accmcfile	  <reaction> <file> (file2) (file3) ...
##
#####################################\n\n
"""

amptools_zlm_ampName = "Zlm"
amptools_vps_ampName = "Vec_ps_refl"


def generate_amptools_cfg(
    quantum_numbers,
    angles,
    fractions,
    datas,
    gens,
    accs,
    bkgnds,
    realAmps,
    fixedAmps,
    fitName,
    cfgFileOutputName,
    basereactName,
    particles,
    header=help_header,
    init_one_val=None,
    datareader="ROOTDataReader",
    add_amp_factor=None,
    append_to_cfg=None,
    append_to_decay=None,
):
    """
    Generate an AmpTools configuration file for a Zlm fit

    Args:
        quantum_numbers (list): List of lists of quantum numbers. For Zlm then [Reflectivity, spin, spin-projection]
        angles (list): List of polarization angles in degrees
        fractions (list): List of polarization fractions
        datas (list): List of data files
        gens (list): List of gen files
        accs (list): List of acc files
        bkgnds (list): List of bkgnd files
        realAmps (list): List of amplitude names that are real
        fixedAmps (list): List of amplitude names that are fixed
        fitName (str): A FitResults (.fit) file will be created with this name prefixed
        cfgFileOutputName (str): Name of output configuration file
        basereactName (str): Base name of reaction
        particles (list): List of particles in reaction
        header (str): Header to append to the top of the file
        datareader (str): Data reader to use
        add_amp_factor (str): Additional factor to add to the amplitude
        append_to_cfg (str): Additional string args to append to the end of the cfg file
        append_to_decay (str): Additional string args to append to the decay factor
        init_one_val (float, None): If not None, all amplitudes will be set to this value

    Returns:
        None, writes a file to cfgFileOutputName
    """

    ############## SET ENVIRONMENT VARIABLES ##############
    atiSetup.setup(globals())
    #######################################################

    cfgInfo = ConfigurationInfo(fitName)

    constraintMap = {}
    refTagMap = {1: "Pos", -1: "Neg"}  # Positive / Negative
    conjugations = {"Re": "+1", "Im": "-1"}  # Real / Imaginary

    for i in range(len(angles)):
        ####################################
        #### Create ReactionInfo object ####
        ####################################

        angle = angles[i]
        fraction = fractions[i]

        reactName = f"{basereactName}_{angle:0>3}"
        scaleName = f"parScale{angle:0>3}"
        parScale = cfgInfo.createParameter(scaleName, 1.0)
        if angle == "0" or angle == "00" or angle == "000":
            parScale.setFixed(True)

        reactionInfo = cfgInfo.createReaction(reactName, particles)  # ReactionInfo*
        # reactionInfo.setNormIntFile(f"{reactName}.ni", False)

        ###################################
        #### SET DATA, GEN, ACC, BKGND ####
        ###################################

        required_srcs = set(["data", "genmc", "accmc"])
        if len(bkgnds) > 0: 
            required_srcs.add("bkgnd")
        reader_names = {}
        reader_argss = {}
        if isinstance(datareader, str):
            reader_parts = datareader.split(" ")
            reader_name = reader_parts[0]
            reader_args = reader_parts[1:] if len(reader_parts) > 1 else []
            for required_src in required_srcs:
                reader_names[required_src] = reader_name
                reader_argss[required_src] = reader_args
        if isinstance(datareader, dict):
            assert set(datareader.keys()) == required_srcs, f"Datareader keys must be {required_srcs} instead of {set(datareader.keys())}"
            for required_src, reader_parts in datareader.items():
                reader_parts = reader_parts.split(" ")
                reader_name = reader_parts[0]
                reader_args = reader_parts[1:] if len(reader_parts) > 1 else []
                reader_names[required_src] = reader_name
                reader_argss[required_src] = reader_args

        data = datas[i]
        data_args = [data] + reader_argss["data"]  # append args
        reactionInfo.setData(reader_names["data"], data_args.copy())
        gen = gens[i]
        gen_args = [gen] + reader_argss["genmc"]
        reactionInfo.setGenMC(reader_names["genmc"], gen_args.copy())
        acc = accs[i]
        acc_args = [acc] + reader_argss["accmc"]
        reactionInfo.setAccMC(reader_names["accmc"], acc_args.copy())

        if len(bkgnds) > 0:
            bkgnd = bkgnds[i]
            bkgnd_args = [bkgnd] + reader_argss["bkgnd"]
            reactionInfo.setBkgnd(reader_names["bkgnd"], bkgnd_args.copy())

        #############################################
        #### DEFINE COHERENT SUMS AND AMPLITUDES ####
        #############################################

        for conj in conjugations.items():
            conjTag, conjVal = conj

            for quantum_number in quantum_numbers:
                ref, L, M = quantum_number[:3]
                J = quantum_number[3] if len(quantum_number) > 3 else None
                assume_zlm = J is None

                if float(fraction) == 1:  # Some terms disappear if the fraction is 1
                    if (ref == 1 and conjTag == "Im") or (ref == -1 and conjTag == "Re"):
                        continue

                refTag = refTagMap[ref]

                sumName = f"{refTag}{conjTag}"
                cfgInfo.createCoherentSum(reactName, sumName)  # returns CoherentSumInfo*

                ampName = zlm_amp_name(ref, L, M) if assume_zlm else vps_amp_name(ref, J, M, L)
                ampInfo = cfgInfo.createAmplitude(reactName, sumName, ampName)  # AmplitudeInfo*

                part = "+1" if int(ref) * int(conjVal) > 0 else "-1"
                if assume_zlm:
                    angularFactor = [f"{amptools_zlm_ampName}", f"{L}", f"{M}", conjVal, part, angle, fraction]
                else:
                    # Vec_ps_refl 1 -1 0 -1  -1  LOOPPOLANG LOOPPOLVAL omega3pi
                    angularFactor = [f"{amptools_vps_ampName}", f"{J}", f"{M}", f"{L}", conjVal, part, angle, fraction]
                if append_to_decay:
                    angularFactor.extend(append_to_decay.split())
                ampInfo.addFactor(angularFactor)
                if add_amp_factor:
                    ampInfo.addFactor(add_amp_factor.split(" "))
                ampInfo.setScale(f"[{scaleName}]")

                if ampName not in constraintMap:
                    constraintMap[ampName] = [ampInfo]
                else:
                    constraintMap[ampName].append(ampInfo)

    #####################################
    ### RANDOMLY INITALIZE AMPLITUDES ###
    #####################################

    bAoV = False
    if init_one_val is not None:
        init_one_val = float(init_one_val)
        bAoV = True

    for amp, lines in constraintMap.items():
        value = init_one_val if bAoV else random.uniform(0.0, 1.0)
        if amp in realAmps:
            lines[0].setReal(True)
        else:
            value += 1j * (init_one_val if bAoV else random.uniform(0.0, 1.0))
        if amp in fixedAmps:
            lines[0].setFixed(True)
        for line in lines[1:]:
            lines[0].addConstraint(line)
        lines[0].setValue(value)

    #########################
    ### WRITE CONFIG FILE ###
    #########################

    cfgInfo.display()
    cfgInfo.write(cfgFileOutputName)

    with open(cfgFileOutputName, "r") as original:
        data = original.read()

    # prepend help_header to top of the file
    data = header + data
    with open(cfgFileOutputName, "w") as modified:
        modified.write(data)
    if append_to_cfg:
        data = data + append_to_cfg
        with open(cfgFileOutputName, "w") as modified:
            modified.write(data)

    return data


def generate_amptools_cfg_from_dict(yaml_file):
    #####################################################
    ################ GENERAL SPECIFICATION ##############
    #####################################################

    generate_success = False
    perform_config_checks = True

    fitName = "PLACEHOLDER_FITNAME"
    cfgFileOutputName = f'{yaml_file["base_directory"]}/amptools.cfg'
    basereactName = "reaction"
    data_folder = yaml_file["data_folder"]
    particles = yaml_file["reaction"].split(" ")
    print("Particles in reaction: ", particles)

    cfgFileOutputFolder = os.path.dirname(cfgFileOutputName)
    os.system(f"mkdir -p {cfgFileOutputFolder}")

    check_map = {}
    check_map["data_folder"] = [os.path.isabs(data_folder), f"Data folder, {data_folder}, is not an absolute path"]

    #####################################################
    ################# CONSTRUCT WAVESET #################
    #####################################################

    used_quantum_numbers = []

    # If user supplies waveset string we will parse and use that
    #   otherwise we will check if buttons have been clicked.
    if yaml_file["waveset"] == "":
        check_string = "You did not select any partial waves... \n Please enter a waveset string and try again."
    else:
        waveset = yaml_file["waveset"].split("_")
        # Using waveset string:  Sp0+_Dp2+
        for wave in waveset:
            if wave == "isotropic": continue
            elif wave in converter:
                used_quantum_numbers.append(converter[wave])
            else:
                check_string = f"Wave {wave} not recognized. Please check your waveset string and try again."
                check_string += f"\n\nExample partial two pseduoscalar wave names: {example_zlm_names[:10]}"
                check_string += f"\n\nExample partial vector pseudoscalar wave names: {example_vps_names[:10]}"
                return check_string, generate_success

    if yaml_file["real_waves"] == "" or yaml_file["fixed_waves"] is None:
        realAmps = yaml_file["real_waves"].split("_")
        print(f"Using real waves: {realAmps}")

    if yaml_file["fixed_waves"] == "" or yaml_file["fixed_waves"] is None:
        fixedAmps = yaml_file["fixed_waves"].split("_")
        print(f"Using fixed waves: {fixedAmps}")

    if len(used_quantum_numbers) == 0:
        check_string = "You did not select any partial waves... \n Please go make some selections or enter a waveset string and try again."
        check_string += f'\nYour waveset string: {yaml_file["waveset"]}'
        check_string += f"\n\nExample partial wave names: {example_zlm_names}"
        return check_string, generate_success

    #####################################################
    ############### CONSTRUCT POLARIZATIONS #############
    #####################################################

    ## Polarization related, reactNames are scaled by the scales parameters
    ##   Default: parScale0 is fixed to 1.0

    used_pols = []
    used_polMags = []
    pols = yaml_file["polarizations"]
    for angle, mag in pols.items():
        used_pols.append(angle)
        used_polMags.append(f"{mag}")
        check_map[angle] = [mag >= 0 and mag <= 1, f"Polarization magnitude must be between 0 and 1 instead of {mag}"]

    #####################################################
    ################# CONSTRUCT DATASETS ################
    #####################################################

    # These are the actual locations that we can still check for
    _datas = [f"{data_folder}/data{pol:0>3}.root" for pol in used_pols]
    _bkgnds = [f"{data_folder}/bkgnd{pol:0>3}.root" for pol in used_pols]
    _gens = [f"{data_folder}/genmc{pol:0>3}.root" for pol in used_pols]
    _accs = [f"{data_folder}/accmc{pol:0>3}.root" for pol in used_pols]

    # Check if these files exist
    for ftype, sources in zip(["data", "bkgnd", "genmc", "accmc"], [_datas, _bkgnds, _gens, _accs]):
        for source in sources:
            if not os.path.isfile(source) and ftype in ["genmc", "accmc"]:
                if os.path.isfile(f"{data_folder}/{ftype}.root"):
                    if ftype == "genmc":
                        _gens = [f"{data_folder}/{ftype}.root" for pol in used_pols]
                    if ftype == "accmc":
                        _accs = [f"{data_folder}/{ftype}.root" for pol in used_pols]
                else:
                    raise FileNotFoundError(f"File {source} does not exist.")
            if not os.path.isfile(source) and ftype in ["data"]:
                raise FileNotFoundError(f"File {source} does not exist.")

    # These are placeholder locations for us to split_mass with
    datas = [f"PLACEHOLDER_DATA_{pol:0>3}" for pol in used_pols]
    gens = [f"PLACEHOLDER_GENMC_{pol:0>3}" for pol in used_pols]
    accs = [f"PLACEHOLDER_ACCMC_{pol:0>3}" for pol in used_pols]
    bkgnds = [f"PLACEHOLDER_BKGND_{pol:0>3}" for pol in used_pols]

    #####################################################
    ################# GENERATE CONFIG FILE ##############
    #####################################################

    fs_check = "Number of final state particles in {}\n Actual: {}\n Expected from Reaction: {}"
    for data, gen, acc, bkgnd in zip(_datas, _gens, _accs, _bkgnds):
        bkgnd_exists = False
        for file, _file in zip([data, gen, acc, bkgnd], ["data", "gen", "acc", "bkgnd"]):
            if _file == "bkgnd":
                bkgnd_exists = os.path.isfile(file)
            else:  # bkgnd is optional so we will not check for its existence
                check_map[file] = [os.path.isfile(file), f"File {file} does not exist."]

            if _file in ["data", "gen", "acc"] or (bkgnd_exists and _file == "bkgnd"):
                _file = ROOT.TFile.Open(file)
                _tree = _file.Get("kin")
                _branch = _tree.GetBranch("NumFinalState")
                _value = array.array("i", [0])
                _branch.SetAddress(_value)
                _tree.GetEntry(0)
                actual_num_fs_parts = _value[0]
                expect_num_fs_parts = len(particles) - 1  # ignore Beam, care only about recoil + daughters
                if "num_final_state" not in check_map:
                    check_map["num_final_state"] = [actual_num_fs_parts == expect_num_fs_parts, fs_check.format(file, actual_num_fs_parts, expect_num_fs_parts)]
                else:
                    check_map["num_final_state"][1] += "\n" + fs_check.format(file, actual_num_fs_parts, expect_num_fs_parts)

    if not bkgnd_exists:
        print("No background file found. AmpTools will assume data file contains pure signal")
        bkgnds = []

    result = ""
    any_check_failed = any([not check[0] for check in check_map.values()])
    if perform_config_checks and any_check_failed:
        result += "Config checks failed:\n\n"
        for check in check_map.values():
            if not check[0]:
                result += check[1] + "\n"
    else:
        result = generate_amptools_cfg(
            used_quantum_numbers,
            used_pols,
            used_polMags,
            datas,
            gens,
            accs,
            bkgnds,
            realAmps,
            fixedAmps,
            fitName,
            cfgFileOutputName,
            basereactName,
            particles,
            header=help_header,
            datareader=yaml_file["datareader"],
            add_amp_factor=yaml_file["add_amp_factor"].strip(),
            append_to_cfg=yaml_file["append_to_cfg"].strip(),
            append_to_decay=yaml_file["append_to_decay"].strip(),
            init_one_val=yaml_file["init_one_val"],
        )

        generate_success = True

        #### FINAL CHECK ####
        valid_keys = ["data", "bkgnd", "accmc", "genmc", "include", "define", "fit", "keyword", "reaction", "normintfile", "sum", "amplitude", "initialize", "scale", "constrain", "permute", "parameter"]
        final_check_string = ""
        for line in result.split("\n"):
            line = line.strip()
            if not line.startswith("#") and line != "":  # if not a comment nor empty line
                if line.split(" ")[0] not in valid_keys:
                    final_check_string += f"Invalid keyword found in the config file: {line}\n"
                    generate_success = False
        if final_check_string != "":
            result = final_check_string

    return result, generate_success

def generate_vecps_cfg(fit_name, reaction_name, help_header, reader_type, reader_args, data_file, acc_file, gen_file, waves, real, phase_constraint, fmt = "cartesian"):
    phase_constrained_waves = {}
    phase_constrained_waves_mag = []
    if phase_constraint == "m_con":
        phase_constrained_waves = {
            "3mpf" : "3f_posm_par",
            "3mp2f" : "3f_posm_par",
            "3mp3f" : "3f_posm_par",
            "3mmf" : "3f_negm_par",
            "3mm2f" : "3f_negm_par",
            "3mm3f" : "3f_negm_par",
            }
        phase_constrained_waves_mag = {
            "3mpf" : "3f_posm_mag_par",
            "3mp2f" : "3f_posm_mag_par",
            "3mp3f": "3f_posm_mag_par",
            "3mmf" : "3f_negm_mag_par",
            "3mm2f" : "3f_negm_mag_par",
            "3mm3f" : "3f_negm_mag_par"
                }
        phase_pars = {par for par in phase_constrained_waves.values()}
        phase_mag_pars = {par for par in phase_constrained_waves_mag.values()}

    if phase_constraint == "l_con":
        phase_constrained_waves = {
            "1pps" : "1s_par",
            "1p0s": "1s_par",
            "1pms" : "1s_par",
            "1mmp" : "1p_par",
            "1m0p" : "1p_par",
            "1mpp" : "1p_par",
            "3mm3f" : "3f_par",
            "3mm2f" : "3f_par", 
            "3mmf" : "3f_par",
            "3m0f" : "3f_par",
            "3mpf" : "3f_par",
            "3mp2f" : "3f_par",
            "3mp3f" : "3f_par"
                }
        phase_constrained_waves_mag = {
            "1pps" : "1s_mag_par",
            "1p0s": "1s_mag_par",
            "1pms" : "1s_mag_par",
            "1mmp" : "1p_mag_par",
            "1m0p" : "1p_mag_par",
            "1mpp" : "1p_mag_par",
            "3mm3f" : "3f_mag_par",
            "3mm2f" : "3f_mag_par", 
            "3mmf" : "3f_mag_par",
            "3m0f" : "3f_mag_par",
            "3mpf" : "3f_mag_par",
            "3mp2f" : "3f_mag_par",
            "3mp3f" : "3f_mag_par"
                }
        phase_pars = {par for par in phase_constrained_waves.values()}
        phase_mag_pars = {par for par in phase_constrained_waves_mag.values()}
    #print(phase_constraint)
    #print(phase_constrained_waves)
    sep = "::"
    amp_type = "Vec_ps_refl"
    signs = ["PosSign", "NegSign"]
    realities = ["Real", "Imag"]
    sum_names = [reality + sign for reality in realities for sign in signs]
    reaction = f"{fit_name} Beam Proton Eta K+ K-"
    sums = [f"sum {fit_name} {temp}" for temp in sum_names]
    TEMstring = f"define TEMstring {reader_args}"
    gen_line = f"genmc {fit_name} {reader_type} {gen_file} {reader_args}"
    acc_line = f"accmc {fit_name} {reader_type} {acc_file} {reader_args}"
    data_line = f"data {fit_name} {reader_type} {data_file} {reader_args}"
    pol_angle = 0
    pol_val = 0.4
    used_phasediff_pars = set()
    cfg = "\n".join([help_header, f"fit {fit_name}", f"reaction {reaction}", "\n".join(sums), TEMstring, gen_line, acc_line, data_line, "\n"])
    for i, wave in enumerate(waves):
        if wave in phase_constrained_waves and phase_constrained_waves[wave] not in used_phasediff_pars:
            phase_par = f"parameter {phase_constrained_waves[wave]} 3.14"
            phase_mag_par = f"parameter {phase_constrained_waves_mag[wave]} 1 fixed"
            used_phasediff_pars.add(phase_constrained_waves[wave])
            used_phasediff_pars.add(phase_constrained_waves[wave])
            cfg += "\n".join([phase_par, phase_mag_par, "\n" ])
        #print(wave)
        quantum_numbers = parse_wave_string(wave)
        for _sum in sum_names:
            wave_real = "real" if real[i] is True and (_sum == "ImagNegSign" or _sum == "RealNegSign") else ""
            reality = "+1" if "Re" in _sum else "-1"
            sign = "+1" if "PosSign" in _sum else "-1"
            pars = "100 100" if wave_real != "real" and wave not in phase_constrained_waves else "100 0"
            amplitude = f"{fit_name}{sep}{_sum}{sep}{wave}"
            phasecon_ampl_decl = ""
            if phase_constraint and wave in phase_constrained_waves:
                #print("inside if statement")
                phasecon_ampl_decl = f"amplitude {amplitude} ComplexCoeff [{phase_constrained_waves_mag[wave]}] [{phase_constrained_waves[wave]}] MagPhi"
            amplitude_decl = f"amplitude {amplitude} {amp_type} {quantum_numbers['j']} {m_names(quantum_numbers['m'])} {l_to_num(quantum_numbers['l'])} {reality} {sign} {pol_angle} {pol_val} dalitz"  
            initialize = f"initialize {amplitude} {fmt} {pars} {wave_real}"
            #print(phasecon_ampl_decl)
            cfg += "\n".join([phasecon_ampl_decl, amplitude_decl, initialize, "\n"])
            
        constrain = "\n".join([f"constrain {fit_name} {sum_names[0]} {wave} {fit_name} {sum_names[2]} {wave}", f"constrain {fit_name} {sum_names[1]} {wave} {fit_name} {sum_names[3]} {wave}", "\n"])
        cfg += constrain
    return cfg


def generate_vecps_cfg_yaml(yaml):
    fit_name = yaml["fit_name"]
    reaction_name = yaml["reaction_name"]
    reader_type = yaml["reader_type"]
    reader_args = yaml["reader_args"]
    data_file = yaml["data_file"]
    accmc_file = yaml["accmc_file"]
    genmc_file = yaml["genmc_file"]
    waves = yaml["waves"]
    real_waves = yaml["real_waves"]
    phase_constraint = yaml["phase_constraint"]
    wave_format = yaml["fmt"]
    output_directory = yaml["output_directory"]
    do_binning = yaml["do_binning"]
    do_fitting = yaml["do_fitting"]
    do_plotting = yaml["do_plotting"]

    
    #if len(fit_name) == 1:
     #   fit_name = fit_name[0]
    if len(reaction_name) == 1:
        reaction_name = reaction_name[0]
    if len(reader_type) == 1 :
        reader_type = reader_type[0]
    if len(reader_args) == 1:
        reader_args = reader_args[0]

    cfgs = []
    
    for fit, data, accmc, genmc, waveset, realset in zip(fit_name, data_file, accmc_file, genmc_file, waves, real_waves):
        cfg = generate_vecps_cfg(fit, reaction_name, help_header, reader_type, reader_args, data, accmc, genmc, waveset["waveset"], realset["realset"], phase_constraint, wave_format)
        cfgs.append(cfg)
        with open(f"{output_directory}/{fit}.cfg", 'w') as cfg_file:
            cfg_file.write(cfg)
        if do_binning:
            if os.path.isdir(f"{output_directory}/{fit}/") is False:
                os.mkdir(f"{output_directory}/{fit}/")
            make_amptools_bins._main(
                        SimpleNamespace(**{"config_path" : f"{output_directory}/{fit}.cfg", 
                            "fit_name" : fit,
                            "data_path" : data,
                            "accmc_path" : accmc,
                            "genmc_path" : genmc,
                            "mass_min" : 1.5,
                            "mass_max" : 2.5,
                            "n_bins" : 25,
                            "t_min" : 0,
                            "t_max" : 2,
                            "energy_min" : 3.0,
                            "energy_max" : 11.6,
                            "bootstrap" : False
                        }))
            if do_fitting:
                run_amptools.main(
                        SimpleNamespace(**{"bin_path" : f"{output_directory}/{fit}/bins.txt",
                            "fit_name" : fit,
                            "num_fits" : 1,
                            "fit_path": "/work/halld/home/ddarulis/etaphi/results/",
                            "num_process" : 12}
                            )
                        )
                get_best_lik.main(
                        SimpleNamespace(**{"path" : output_directory}))
            if do_plotting:
                files = [f for f in os.listdir(f"{output_directory}") if re.match(rf"{fit}_bin_(0-9){1,2}.fit", f)]
                plot_fits.main(
                        SimpleNamespace(**{"files" : files,
                            "num_process" : len(files)
                            }
                            )
                        )


    return cfgs
    

def generate_etaphi_hybrid_cfg(fit_name, reaction_name, help_header, reader_type, reader_args, data_file, acc_file, gen_file, piecewise_waves, real, bws, bw_bounds, bw_real, no_bins, fmt="cartesian"): 

    sep = "::"
    amp_type = "Vec_ps_refl"
    signs = ["PosSign", "NegSign"]
    realities = ["Real", "Imag"]
    sum_names = [reality + sign for reality in realities for sign in signs]
    reaction = f"{fit_name} Beam Proton Eta K+ K-"
    sums = [f"sum {fit_name} {temp}" for temp in sum_names]
    gen_line = f"genmc {fit_name} {reader_type} {gen_file} {reader_args}"
    acc_line = f"accmc {fit_name} {reader_type} {acc_file} {reader_args}"
    data_line = f"data {fit_name} {reader_type} {data_file} {reader_args}"
    pol_angle = 0
    pol_val = 0.4
    lo_edge = 1.5
    hi_edge = 1.98
    cfg = "\n".join([f"fit {fit_name}", f"reaction {reaction}", "\n".join(sums), gen_line, acc_line, data_line, "\n"])
    for i, wave in enumerate(piecewise_waves):
        quantum_numbers = parse_wave_string(wave)
        par_name_pos = []
        par_name_neg = []
        for j in range(no_bins):
            for r in ["Re", "Im"]:
                for s in ["Pos", "Neg"]:
                    par_name = f"pcwsBin_{j}{r}{s}{quantum_numbers['l']}{quantum_numbers['m']}"
                    if s == "Pos":
                        par_name_pos.append(par_name)
                    else:
                        par_name_neg.append(par_name)
    for par_name in par_name_pos:
        cfg += f"parameter {par_name} 1 \n"
    for par_name in par_name_neg:
        cfg += f"parameter {par_name} 1 \n"

    for i, bw in enumerate(bws):
        for _sum in sum_names:
            quantum_numbers = parse_wave_string(bw)
            bw_mass = 1680 if quantum_numbers['l'] == 'p' else 1850
            widths = {1680 : 0.15, 1850 : 0.87}
            width_bounds = {1680 : (0.1, 0.2), 1850 : (0.064, 0.115)}
            if bw_bounds[i] == "True":
                wave_bounds = f"{widths[bw_mass]} bounded {str(width_bounds[bw_mass][0])} {str(width_bounds[bw_mass][1])}"
            else:
                wave_bounds = f"{widths[bw_mass]} fixed"
            amplitude = f"amplitude {fit_name}{sep}{_sum}{sep}{bw} BreitWigner {bw_mass} {wave_bounds} 1 2 34"   
            cfg += "".join([amplitude, "\n"])
    
    for i, wave in enumerate(piecewise_waves):
        quantum_numbers = parse_wave_string(wave)
        for _sum in sum_names:
            suffix = f"Neg{quantum_numbers['l']}{quantum_numbers['m']}" if "Neg" in _sum else f"Pos{quantum_numbers['l']}{quantum_numbers['m']}"
            amplitude = f"amplitude {fit_name}{sep}{_sum}{sep}{wave} Piecewise {lo_edge} {hi_edge} {no_bins} 234 {suffix} ReIm "
            if "PosSign" in _sum:
                for par in par_name_pos:
                    amplitude += f"[{par}] "
            else:
                for par in par_name_neg:
                    amplitude += f"[{par}] "
            amplitude += "\n"

            #wave_real = "real" if real[i] is True and (_sum == "ImagNegSign" or _sum == "RealNegSign") else ""
            reality = "+1" if "Re" in _sum else "-1"
            sign = "+1" if "PosSign" in _sum else "-1"
            amplitude_ang = f"{fit_name}{sep}{_sum}{sep}{wave}"
            amplitude_ang_decl = f"amplitude {amplitude_ang} {amp_type} {quantum_numbers['j']} {m_names(quantum_numbers['m'])} {l_to_num(quantum_numbers['l'])} {reality} {sign} {pol_angle} {pol_val}  dalitz"  
            if "1p0s" in amplitude_ang:
                initialize = f"initialize {amplitude_ang} {fmt} 1 0 real fixed"
            else:
                initialize = f"initialize {amplitude_ang} {fmt} 1 1 fixed"
            cfg += "\n".join([amplitude, amplitude_ang_decl, initialize, "\n"])
    
        constrain = "\n".join([f"constrain {fit_name} {sum_names[0]} {wave} {fit_name} {sum_names[2]} {wave}", f"constrain {fit_name} {sum_names[1]} {wave} {fit_name} {sum_names[3]} {wave}", "\n"])
        cfg += constrain


    return cfg
 
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an AmpTools configuration file for a Zlm fit")
    parser.add_argument("--yaml_name", type=str, default="/work/halld/home/ddarulis/etaphi/MC/PWA/fitting/automated_configs/yaml_files/test.yaml", help="Path a configuration yaml file")
    parser.add_argument("--test", type=bool, default=False, help="Test configuration generating functions")
    args = parser.parse_args()
    if args.test:
        #result = generate_vecps_cfg("test_fit", "etaphi", help_header, "ROOTDataReaderTEM", "-1 0.0 1.5 1.54", "/work/home/halld/ddarulis/etaphi/MC/PWA/DataWeightedTree_2024_06_filtered.root", "/work/halld/home/ddarulis/etaphi/MC/PWA/phi_eta_gg_2019_mc_kkin_2024_06_amptools_tree.root", "/work/halld/home/ddarulis/etaphi/MC/PWA/phi_eta_gg_2019_mc_kin_2024_06_amptools_tree.root", ["1pps","1p0s", "1pms", "1m0p", "3mpf", "3mp2f"], [True, False, False, False, False, False], "m_con")
        result_hybrid = generate_etaphi_hybrid_cfg("test_hybrid_fit", "etaphi", help_header, "ROOTDataReaderTEM", "0.0 0.5 3.0 11.6 1.5 1.82", "/work/halld/home/ddarulis/etaphi/MC/PWA/DataWeightedTree_2024_06_filtered.root", "/work/halld/home/ddarulis/etaphi/MC/PWA/phi_eta_gg_2019_mc_kin_2024_06_amptools_tree.root", "/work/halld/home/ddarulis/etaphi/MC/PWA/phi_eta_gg_2019_mc_kin_2024_06_amptools_tree.root", ["1p0s", "1pps", "1pms"], [True, False], [], ["fixed"], [False], 8)
#        print(result)
        print(result_hybrid)
        with open("/work/halld/home/ddarulis/etaphi/MC/PWA/fitting/test_hybrid.cfg", "w") as f:
            f.write(result_hybrid)
        sys.exit(0)
    yaml_name = args.yaml_name

    cwd = os.getcwd()

    print("\n---------------------")
    print(f"Running {__file__}")
    print(f"  yaml location: {yaml_name}")
    print("---------------------\n")

    yaml_file = load_yaml(yaml_name)
    
    #result, generate_success = generate_amptools_cfg_from_dict(yaml_file)
    cfgs = generate_vecps_cfg_yaml(yaml_file)
    print(cfgs[0])
    #if generate_success:
     #   with open(f"{yaml_file['base_directory']}/amptools.cfg", "w") as f:
      #      f.write(result)
       # print("\nConfiguration file successfully generated")
    #else:
     #   print(result)
      #  print("\nConfiguration file generation failed!")
