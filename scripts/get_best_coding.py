import argparse
from pathlib import Path
import os


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-f', '--coding_file',
                               dest='coding_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-d', '--annotation_directory',
                               dest='input_dir',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--output_directory',
                               dest='out_dir',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-op', '--output_n-proteins',
                               dest='nproteins_out_file',
                               type=lambda p: Path(p).resolve(),
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def get_best_coding(genome_id, coding_results_dict):
    sorted_callers = sorted(coding_results_dict, key=lambda line: float(line[4]), reverse=True)
    diff = 0.1 * float(sorted_callers[0][4])

    if float(sorted_callers[0][4]) - float(sorted_callers[1][4]) > diff:
        caller = sorted_callers[0][5]
        n_prot = sorted_callers[0][2]
    else:
        caller = "prodigal-11"
        for line in sorted_callers:
            if line[5] == caller:
                n_prot = line[2]

    return (genome_id, caller, n_prot)



def main():
    args = parse_args()

    # read the raw coding information file
    table_raw = [line.strip().split("\t") for line in open(args.coding_file).readlines()]

    # filter PHANOTATE annotation
    table = [line for line in table_raw if line[5] != "phanotate"]

    # put the coding results in a dictionary, genome_id as key
    coding_results_genome = dict()
    for line in table:
        if line[0] in coding_results_genome:
            coding_results_genome[line[0]].append(line)
        else:
            coding_results_genome[line[0]] = [line]

    #
    final_callers = list()
    for genome, coding_results in coding_results_genome.items():
        final_callers.append(get_best_coding(genome, coding_results))


    #
    out_nproteins = open(os.path.abspath(args.nproteins_out_file), "w")
    for tuple in final_callers:
        in_abspath = os.path.abspath(args.input_dir)
        o_abspath  = os.path.abspath(args.out_dir)
        op_abspath = os.path.abspath(args.nproteins_out_file)

        print(f"Copying {tuple[0]},{tuple[1]} to {o_abspath}")
        os.system("cp {}/{}/{}_* {}".format(in_abspath,
                                           tuple[1],
                                           tuple[0],
                                           o_abspath))

        out_nproteins.write(f'{tuple[0]}\t{tuple[1]}\t{tuple[2]}\n')













if __name__ == '__main__':
    main()
