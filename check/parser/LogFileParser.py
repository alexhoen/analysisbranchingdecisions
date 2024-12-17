import math
import os, re

import pandas as pd
import numpy
from numpy.core.defchararray import lower
from numpy.core.numeric import infty
from pandas import factorize

OPTIMAL_STATUS = "opt"
GAP_STATUS = "gap"
TIMELIMIT_STATUS = "TL"
INF_STATUS = "INF"

STATUS = 'Status'
EASY = 'Easy'
SAFE_BOUNDING = 'SF'
RECONSTRUCTION =  'Recon'
FACTORIZATION = 'factor'
FIXINT =  'fixint'
EXACT_SOPLEX = 'exact'
ERROR = 'error'

HAS_UNBOUNDED_VARS_BEFORE_PRESOLVING = 'before'
HAS_UNBOUNDED_VARS_AFTER_PRESOLVING = 'after'

def get_filename(name: str, dir: str):
    for root, dirs, files in os.walk(dir):
        for filename in files:
            if filename.endswith(".out") and filename.__contains__(name):
                return dir + "/" + filename


def cast(rational_number: str):
    if "/" in rational_number:
        return float(rational_number.split("/")[0])/float(rational_number.split("/")[1])
    else:
        return float(rational_number)


def call(filename):
    if "prodplan1" in filename:
        error = 1
        nodes = 34
        safebounding = 33
        return (nodes, safebounding, 0, 0, 0, error)
    if "ns2080781" in filename:
        error = 181 + 10 +47
        nodes = 16047 + 565277 +10 -8777 +120171
        reconstruct = 109
        factor = 984
        exactsoplex = 7540 + 10033 + 373
        safebounding = (3825 + #verifying  by soplex
                            91 + #Node (objlimit by pseudoobjective) verified
                            668741)# verifying by applying Neumaier
        return (nodes, safebounding, reconstruct, factor, exactsoplex, error)

def analyse(root, filename):
    floating_point_expr = "[-+]?[0-9]*\.?[0-9]*"

    time_expr = re.compile("Instance solved in\s+(\S+)")
    timelimit = re.compile("Time limit reached in\s+(\S+)")
    nofeasexpr = re.compile("Instance solved with no feasible solution in\s+(\S+)")
    event_time_expr = re.compile("Instance solved in \d+.\d+s spent\s+(\S+)")

    best_solution_feasible_expr = re.compile("Best solution:\s+(\S+)")
    best_solution_recons_expr = re.compile("Best solution:\s+\d+\s+feasible\s+(\S+)")
    best_solution_fixint_expr = re.compile("Best solution:\s+\d+\s+feasible\s+\d+\s+reconstr\s+(\S+)")
    best_solution_error_expr = re.compile("Best solution:\s+\d+\s+feasible\s+\d+\s+reconstr\s+\d+\s+fix-ints\s+(\S+)")

    nodefeas_feasible_expr = re.compile("Node feasible:\s+(\S+)")
    nodefeas_recons_expr = re.compile("Node feasible:\s+\d+\s+feasible\s+(\S+)")
    nodefeas_factor_expr = re.compile("Node feasible:\s+\d+\s+feasible\s+\d+\s+reconstr\s+(\S+)")
    nodefeas_exact_expr = re.compile("Node feasible:\s+\d+\s+feasible\s+\d+\s+reconstr\s+\d+\s+factor\s+(\S+)")
    nodefeas_error_expr = re.compile("Node feasible:\s+\d+\s+feasible\s+\d+\s+reconstr\s+\d+\s+factor\s+\d+\s+exact\s+(\S+)")

    nodedel_soplex_expr = re.compile("Node deletion:\s+(\S+)")
    nodedel_NeuShb_expr = re.compile("Node deletion:\s+\d+\s+soplex\s+(\S+)")
    nodedel_reconstr_expr = re.compile("Node deletion:\s+\d+\s+soplex\s+\d+\s+NeuShb\s+(\S+)")
    nodedel_factor_expr = re.compile("Node deletion:\s+\d+\s+soplex\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+(\S+)")
    nodedel_exact_expr = re.compile("Node deletion:\s+\d+\s+soplex\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+\d+\s+factor\s+(\S+)")
    nodedel_error_expr = re.compile("Node deletion:\s+\d+\s+soplex\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+\d+\s+factor\s+\d+\s+exact\s+(\S+)")

    nodedpseudo_feasible_expr = re.compile("Node del-pseudo:\s+(\S+)")
    nodedpseudo_error_expr = re.compile("Node del-pseudo:\s+\d+\s+feasible\s+(\S+)")

    nodeinfeas_infeas_expr = re.compile("Node infeasible:\s+(\S+)")
    nodeinfeas_NeuShb_expr = re.compile("Node infeasible:\s+\d+\s+infeasible\s+(\S+)")
    nodeinfeas_reconstr_expr = re.compile("Node infeasible:\s+\d+\s+infeasible\s+\d+\s+NeuShb\s+(\S+)")
    nodeinfeas_exact_expr = re.compile("Node infeasible:\s+\d+\s+infeasible\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+(\S+)")
    nodeinfeas_error_expr = re.compile("Node infeasible:\s+\d+\s+infeasible\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+\d+\s+exact\s+(\S+)")

    nodeobjlim_feasible_expr = re.compile("Node objlimit:\s+(\S+)")
    nodeobjlim_NeuShb_expr = re.compile("Node objlimit:\s+\d+\s+feasible\s+(\S+)")
    nodeobjlim_reconstr_expr = re.compile("Node objlimit:\s+\d+\s+feasible\s+\d+\s+NeuShb\s+(\S+)")
    nodeobjlim_factor_expr = re.compile("Node objlimit:\s+\d+\s+feasible\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+(\S+)")
    nodeobjlim_exact_expr = re.compile("Node objlimit:\s+\d+\s+feasible\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+\d+\s+factor\s+(\S+)")
    nodeobjlim_error_expr = re.compile("Node objlimit:\s+\d+\s+feasible\s+\d+\s+NeuShb\s+\d+\s+reconstr\s+\d+\s+factor\s+\d+\s+exact\s+(\S+)")

    nodenotsolved_pseudo_expr = re.compile("Node not solved:\s+(\S+)")
    nodenotsolved_error_expr = re.compile("Node not solved:\s+\d+\s+pseudo\s+(\S+)")

    nodeobjlimps_easy_expr = re.compile("Node feasible:\s+(\S+)")
    nodeobjlimps_NeuShb_expr = re.compile("Node feasible:\s+\d+\s+easy\s+(\S+)")
    nodeobjlimps_exact_expr = re.compile("Node feasible:\s+\d+\s+easy\s+\d+\s+NeuShb\s+(\S+)")
    nodeobjlimps_error_expr = re.compile("Node feasible:\s+\d+\s+easy\s+\d+\s+NeuShb\s+\d+\s+exact\s+(\S+)")

    recalc_value_failed = re.compile("\s+recalculating floating-point failed\s+(\S+)")
    objlimit_value_failed = re.compile("\s+ObjLimit - Safe bounding failed\s+(\S+)")



    time = False
    with open(root + "/" +filename) as f:
        status  = "unknown"
        atime  = numpy.zeros([2])
        asol  = numpy.zeros([4])
        afeas = numpy.zeros([5])
        adel = numpy.zeros([6])
        apseudo = numpy.zeros([2])
        ainfeas = numpy.zeros([5])
        aobjlim = numpy.zeros([6])
        anotsolved = numpy.zeros([2])
        aobjlimpseudo = numpy.zeros([4])
        time = False

        unresolved_nodes = []
        objlim_errors = []
        optimal_value = infty
        lowerbound = infty
        unbounded_before_presolving = False
        unbounded_after_presolving = False

        # if "prodplan1" in filename or "ns2080781" in filename:
        #     (nodes, safebounding, r, factori, exactsoplex, err) = call(filename)
        #     print("read hard coded data")
        #     return [filename, TIMELIMIT_STATUS, 0, 0, 0, 0, 0, 0, safebounding, r, factori, 0, exactsoplex, err, True, False]

        for line in f.readlines():
            ovm = recalc_value_failed.match(line)
            if ovm:
                # print(ovm.groups()[0].replace("s", ""))
                unresolved_nodes.append(cast(ovm.groups()[0]))

            ovm = recalc_value_failed.match(line)
            if line.startswith("Before presolving problem has unbounded variables"):
                unbounded_before_presolving = True
            if line.startswith("Before presolving problem has no unbounded variables"):
                unbounded_before_presolving = False
            if line.startswith("Problem has unbounded variables"):
                unbounded_after_presolving = True
            if line.startswith("Problem has no unbounded variables"):
                unbounded_after_presolving = False

            if line.startswith("FAILED:") or line.startswith("	FAILED: ") :
                line = line.replace("	FAILED", "FAILED")
                if line.startswith("FAILED: NODE_FEASIBLE not v") :
                    print(line, end="")
                elif line.startswith("FAILED: NODE_FEASIBLE infeasible ") :
                    print(line, end="")
                elif line.startswith("FAILED: ERROR in Exact-SoPlex") :
                    print(line, end="")
                elif line.startswith("FAILED: Best bound"):
                    split = line.split("Best bound")[1].split("exceeds bound")
                    objlim_errors.append(cast(split[0]))
                    print(cast(split[0]), end=" exceeds ")
                    print(cast(split[1]))
                elif line.startswith("FAILED: not verified"):
                    split = line.split("not verified ")[1].split(" ")
                    objlim_errors.append(cast(split[0]))
                    print(cast(split[0]), end=" exceeds ")
                    print(cast(split[1]))
                elif line.startswith("FAILED: ERROR in floating-point SoPlex") or line.startswith("FAILED: SoPlex returned sta"):
                    continue
                elif line.startswith("FAILED: found new best solution with solution value"):
                    print(line, end="")
                elif line.startswith("FAILED: Instead of infeasible the node has a lower bound of"):
                    print(line, end="")
                    # split = line.split("value")[1].split("could")
                    # print(cast(split[0]), end=" can not be casted to a solution\n ")
                else:
                    print("untreated Fail")
                    print(line)
            if "error" in line:
                print(line, end="")



            time_match = time_expr.match(line)
            if time_match :
                atime[0] = time_match.groups()[0].replace("s", "")
                time_match = event_time_expr.match(line)
                if time_match:
                    atime[1] = time_match.groups()[0]
                time = True
                if "infeasible" in line:
                    status = INF_STATUS
                    optimal_value = "-inf"
                    lowerbound = optimal_value
                else:
                    status = GAP_STATUS
                    optimal_value = cast(line.split(",")[1].replace("].\n", ""))
                    lowerbound = cast(line.split(",")[0].split("[")[1])
            else:
                match = timelimit.match(line)
                if match :
                    time = True
                    status = TIMELIMIT_STATUS
                    optimal_value = cast(line.split(",")[1].replace("].\n", ""))
                    lowerbound = "-inf"
                else:
                    match = nofeasexpr.match(line)
                    if match :
                        time = True
                        status = GAP_STATUS
                        lowerbound = cast(line.split(",")[0].split("[")[1])
                        optimal_value = "inf"

            if time:
                # match = best_solution_feasible_expr.match(line)
                # if match:
                #     asol[0] = match.groups()[0]
                #     asol[1] = best_solution_recons_expr.match(line).groups()[0]
                #     asol[2] = best_solution_fixint_expr.match(line).groups()[0]
                #     asol[3] = best_solution_error_expr.match(line).groups()[0]
                #     continue
                match = nodefeas_feasible_expr.match(line)
                if match:
                    afeas[0] = match.groups()[0]
                    afeas[1] =nodefeas_recons_expr.match(line).groups()[0]
                    afeas[2] =nodefeas_factor_expr.match(line).groups()[0]
                    afeas[3] =nodefeas_exact_expr.match(line).groups()[0]
                    afeas[4] =nodefeas_error_expr.match(line).groups()[0]
                    continue
                match = nodedel_soplex_expr.match(line)
                if match:
                    adel[0] = match.groups()[0]
                    adel[1] = nodedel_NeuShb_expr.match(line).groups()[0]
                    adel[2] = nodedel_reconstr_expr.match(line).groups()[0]
                    adel[3] = nodedel_factor_expr.match(line).groups()[0]
                    adel[4] = nodedel_exact_expr.match(line).groups()[0]
                    adel[5] = nodedel_error_expr.match(line).groups()[0]

                match = nodedpseudo_feasible_expr.match(line)
                if match:
                    apseudo[0] = match.groups()[0]
                    apseudo[0] = nodedpseudo_error_expr.match(line).groups()[0]
                    continue

                match = nodeinfeas_infeas_expr.match(line)
                if match:
                    ainfeas[0] = match.groups()[0]
                    ainfeas[1] = nodeinfeas_NeuShb_expr.match(line).groups()[0]
                    ainfeas[2] = nodeinfeas_reconstr_expr.match(line).groups()[0]
                    ainfeas[3] = nodeinfeas_exact_expr.match(line).groups()[0]
                    ainfeas[4] = nodeinfeas_error_expr.match(line).groups()[0]
                    continue
                match = nodeobjlim_feasible_expr.match(line)
                if match:
                    aobjlim[0] = match.groups()[0]
                    aobjlim[1] = nodeobjlim_NeuShb_expr.match(line).groups()[0]
                    aobjlim[2] = nodeobjlim_reconstr_expr.match(line).groups()[0]
                    aobjlim[3] = nodeobjlim_factor_expr.match(line).groups()[0]
                    aobjlim[4] = nodeobjlim_exact_expr.match(line).groups()[0]
                    aobjlim[5] = nodeobjlim_error_expr.match(line).groups()[0]
                    continue

                match = nodenotsolved_pseudo_expr.match(line)
                if match:
                    anotsolved[0] = match.groups()[0]
                    anotsolved[1] = nodenotsolved_error_expr.match(line).groups()[0]
                    continue

                match = nodeobjlimps_easy_expr.match(line)
                if match:
                    aobjlimpseudo[0] = match.groups()[0]
                    aobjlimpseudo[1] = nodeobjlimps_NeuShb_expr.match(line).groups()[0]
                    aobjlimpseudo[2] = nodeobjlimps_exact_expr.match(line).groups()[0]
                    aobjlimpseudo[3] = nodeobjlimps_error_expr.match(line).groups()[0]
                    continue



        easy = asol[0] + afeas[0] + adel[0] + apseudo[0] + ainfeas[0] + aobjlim[0] + anotsolved[0] + aobjlimpseudo[0]
        neu = adel[1] + ainfeas[1] + aobjlim[1] + aobjlimpseudo[1]
        reconstruct = asol[1] + afeas[1] + adel[2] + ainfeas[2] + aobjlim[2]
        factor = afeas[2] + adel[3] + aobjlim[3]
        fixint = asol[2]
        exact = afeas[3] + adel[4] + ainfeas[3] + aobjlim[4] + aobjlimpseudo[2]
        error = asol[3] + afeas[4] + adel[5] + apseudo[1] + ainfeas[4] + aobjlim[5] + anotsolved[1] + aobjlimpseudo[3]
        c = 0
        d = 0
        if optimal_value != "-inf" and optimal_value != "inf":
            for o in unresolved_nodes:
                if optimal_value + 0.0000001 < o:
                    c += 1
            for o in objlim_errors:
                if optimal_value + 0.0000001 < o:
                    d += 1


        if lowerbound == optimal_value:
            lowerbound = "-"
            status = "opt"
        if lowerbound != "-inf" and lowerbound != "inf" and lowerbound != "-":
            if math.isnan(lowerbound):
                print("Rational number could not be parsed! Assuming no gap. please check")
                status = "open"

        if not time:
            status = TIMELIMIT_STATUS

        if status == TIMELIMIT_STATUS:
            print("timelimit!")
        elif len(objlim_errors) != 0 :
            print()
            print(d, end=" of objlim errors ")
            print(len(objlim_errors), end=" could be resolved\n")

        return [filename, status, optimal_value, lowerbound, c, atime[0], atime[1], easy, neu, reconstruct, factor, fixint, exact, error, unbounded_before_presolving, unbounded_after_presolving]


def print_details(dataframe):
    all = int(dataframe[EASY].sum() + dataframe[SAFE_BOUNDING].sum() + dataframe[RECONSTRUCTION].sum() + dataframe[
        FACTORIZATION].sum() + dataframe[FIXINT].sum() + dataframe[EXACT_SOPLEX].sum() + dataframe[ERROR].sum())
    print(all, end='')
    print(" & ", end='')
    # print(int(dataframe[EASY].sum() + dataframe[SAFE_BOUNDING].sum()), end='')
    # print(" & ", end='')
    print(round((dataframe[EASY].sum() + dataframe[SAFE_BOUNDING].sum())/all*100,2), end='')
    print("\%  & ", end='')
    # print(int(dataframe[RECONSTRUCTION].sum()), end='')
    # print(" & ", end='')
    print(round(dataframe[RECONSTRUCTION].sum()/all*100,2), end='')
    print("\%  & ", end='')
    # print(int(dataframe[FACTORIZATION].sum()), end='')
    # print(" & ", end='')
    print(round(dataframe[FACTORIZATION].sum()/all*100,2), end='')
    print("\% & ", end='')
    # print(int(dataframe[FIXINT].sum() + dataframe[EXACT_SOPLEX].sum()), end='')
    # print(int(dataframe[EXACT_SOPLEX].sum()), end='')
    # print(" & ", end='')
    print(round(dataframe[EXACT_SOPLEX].sum()/all*100,2), end='')
    print("\%  & ", end='')
    # print(int(dataframe[ERROR].sum()), end='')
    print(round(dataframe[ERROR].sum()/all*100,8), end='\% ')
    print("\\\\")


if __name__ == '__main__':

    files_with_errors = []
    source_folder = ("/home/alexander/Downloads/results")
    results = []
    i = 0
    for root, dirs, files in os.walk(source_folder):
        for filename in files:
            # if i <= 1:
            #     i+=1
            #     continue
            if filename.endswith(".out") and not filename.startswith("check") and filename.startswith("ahoen.numerical.10"):
                    # and (not "alu1" in filename  and not "dfn6" in filename and not "neos-105" in filename):
                print("\n" + filename)
                try:
                    result = analyse(root, filename)
                    results.append(result)
                except RuntimeError:
                    files_with_errors.append(filename)
                    print("\n Runtime error occurred")
                # i += 1
                # if i == 1:
                #     break
    print(results)
    # return [filename, status, optimal_value, lowerbound, c,
    # atime[0], atime[1], easy, neu, reconstruct,
    # factor, fixint, exact, error]

    df = pd.DataFrame(results,
                      columns = ['Instance' , STATUS , 'Optimal', 'LB', 'RecoverableNodes',
                                 'Time', 'Eventtime' , EASY, SAFE_BOUNDING, RECONSTRUCTION,
                                 FACTORIZATION, FIXINT, EXACT_SOPLEX, ERROR, HAS_UNBOUNDED_VARS_BEFORE_PRESOLVING,
                                 HAS_UNBOUNDED_VARS_AFTER_PRESOLVING])


    df = df.sort_values(by=['Instance'])
    print(df)

    print(len(df), end='')
    print(" & ", end='')
    print(len(df[(df[STATUS] == OPTIMAL_STATUS) | (df[STATUS] == INF_STATUS)] ), end='')
    print(" & ", end='')
    print(len(df[df[STATUS] == GAP_STATUS]), end='')
    print(" & ", end='')
    print(len(df[df[STATUS] == TIMELIMIT_STATUS]), end='')
    print("\\\\")

    print("-------------------")
    print_details(dataframe=df)
    print_details(dataframe=df[df[HAS_UNBOUNDED_VARS_AFTER_PRESOLVING] == True])

    print("-------------------")

    for ind in df.index:
        print(df['Instance'][ind].replace(".exact.M640.default.out", ""). replace("ahoen.", "").replace("_","\\_"), end='')
        print(" & ", end='')
        print(df[STATUS][ind], end='')
        print(" & ", end='')
        print(df['LB'][ind], end='')
        print(" & ", end='')
        print(df['Optimal'][ind], end='')
        print(" & ", end='')
        print(int(df[EASY][ind]), end='')
        print(" & ", end='')
        print(int(df[SAFE_BOUNDING][ind]), end='')
        print(" & ", end='')
        print(int(df['RecoverableNodes'][ind]), end='')
        print(" & ", end='')
        print(int(df[RECONSTRUCTION][ind]), end='')
        print(" & ", end='')
        print(int(df[FACTORIZATION][ind]), end='')
        print(" & ", end='')
        print(int(df[FIXINT][ind]), end='')
        print(" & ", end='')
        print(int(df[EXACT_SOPLEX][ind]), end='')
        print(" & ", end='')
        print(int(df[ERROR][ind]), end='')
        print("\\\\")


    print("")
    d = df[(df[STATUS] != TIMELIMIT_STATUS) & (df[ERROR] > 0 )]
    for index, row in d.iterrows():
        print(row['Instance'], end = "")
        print(row[EASY] + row[SAFE_BOUNDING] + row[RECONSTRUCTION] + row[FACTORIZATION] + row[FIXINT] + row[ERROR] +row[EXACT_SOPLEX])


    print("")
    d = df#[(df[STATUS] == TIMELIMIT_STATUS)]
    for index, row in d.iterrows():
        print(row['Instance'], end = "")
        print(row[EASY] + row[SAFE_BOUNDING] + row[RECONSTRUCTION] + row[FACTORIZATION] + row[FIXINT] + row[ERROR] +row[EXACT_SOPLEX])
    # df2 = df[df['UnbBefore'] == True]

    # print(df2)
    # total = df['Neu'].sum() + df['Easy'].sum() +df['Recon'].sum()  + df['factor'].sum() + df['fixint'].sum() + df['exact'].sum()
    # print( str(total) + " & " + str(df['Neu'].sum() + df['Easy'].sum()) + " & " + str(df['Recon'].sum()) + " & " + str(df['factor'].sum()) + " & " +str(df['fixint'].sum()) + " & " + str(df['exact'].sum()) + " & " + str(df['error'].sum()) + " \\\\")