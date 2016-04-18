__author__ = 'etseng@pacb.com'

"""
for each directory in runs/<sample>,
generate a matching smrtpipe/<sample>, then

1. make smrtpipe/<sample>/input.fofn
2. generate script to make bax2bam
3. generate script to run ccs
4. submit script


"""

#!/usr/bin/env python
import os, sys, glob

def check_file_or_run(filename):
    def func(f):
        def g(bash_f, args, **kwargs):
            if args.force or not os.path.exists(filename):
                f(bash_f, args)
            else:
                print >> sys.stderr, "{0} already exists. No need to run {1}.".format(filename, f.__name__)
        return g
    return func


def run_bax2bam():
    movies = {} # movie name --> list of .bax.h5 (use this to handle cases when some .bax.h5 missing)
    for line in open('input.fofn'):
        file = line.strip()
        movie = os.path.basename(file).split('.')[0]
        if movie not in movies: movies[movie] = []
        movies[movie].append(file)

    i = 0
    for movie, files in movies.iteritems():
        bash_f = open("cmd_" + str(i), 'w')
        if not os.path.exists(movie + '.subreads.bam'):
            bash_f.write("bax2bam {0}\n".format(" ".join(files)))

        if not os.path.exists(movie + '.ccs.report'):
            cmd = "ccs --noPolish --minLength 100 --maxLength 15000 --minPasses 1 --minPredictedAccuracy 0.80 " + \
                "--numThreads {0} --minZScore -999999 --maxDropFrac 0.8 ".format(args.cpus) + \
                "--reportFile {0}.ccs.report {0}.ccs.bam {0}.subreads.bam".format(movie)
            bash_f.write(cmd + '\n')
        bash_f.close()
        i += 1

def main(args):
    for d in os.listdir('runs'):
        d2 = os.path.join('runs', d)
        if not os.path.isdir(d2): continue
        d3 = os.path.join('smrtpipe', d)
        if not os.path.exists(d3):
            os.makedirs(d3)
        print >> sys.stderr, "processing", d3

        files = glob.glob(os.path.abspath(d2) + "/*/Analysis_Results/*.bax.h5")
        fofn = os.path.join(d3, 'input.fofn')
        with open(fofn, 'w') as f:
            f.write("\n".join(files) + '\n')

        os.chdir(d3)
        run_bax2bam()
        os.chdir('../../')


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("")
    parser.add_argument("--cpus", type=int, default=12)


    args = parser.parse_args()
    main(args)