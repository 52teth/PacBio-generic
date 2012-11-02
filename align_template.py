import os, sys

# URL = task://019807/P_Mapping/align_001of008
# Input = file://019807/input.chunk001of008.fofn
# Input = file://019807/data/filtered_regions.chunk001of008.fofn
# Output = file://019807/data/aligned_reads.chunk001of008.cmp.h5

def main(input_dir, ref, **kwargs):
    """
    """
    bestn = 1
    for key in kwargs:
        if key == 'bestn':
            bestn = kwargs[key]
        
    
    input_dir = os.path.abspath(input_dir)
    ref = os.path.abspath(ref)
    files = [line.strip() for line in open(os.path.join(input_dir, 'input.fofn'))]
    total = len(files)
    
    if not os.path.exists(os.path.join(input_dir, "data", "filtered_regions")):
        os.makedirs(os.path.join(input_dir, "data", "filtered_regions"))


    for i0, file in enumerate(files):
        i = i0 + 1
        input = os.path.join(input_dir, "input.chunk{0:03d}of{1:03d}.fofn".format(i, total))
        with open(input, 'w') as f: f.write(file + '\n')        
        rgn = os.path.join(input_dir, "data", "filtered_regions", os.path.basename(file).replace('.bas.h5','.rgn.h5'))
        filtered = os.path.join(input_dir, "filtered_regions.chunk{0:03d}of{1:03d}.fofn".format(i, total))
        with open(filtered, 'w') as f: f.write(rgn + '\n')
        output = os.path.join(input_dir, "aligned_reads.chunk{0:03d}of{1:03d}.cmp.h5".format(i, total))
 
        with open(os.path.join(input_dir,"align_{0:03d}of{1:03d}.sh".format(i, total)), 'w') as f:
            f.write(
"""
echo 'Started on' `date -u`;

compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=15  --noXML --h5mode=w --h5fn={output} --minAccuracy=0.75 --minLength=50  -x -minMatch 12 -x -bestn {bestn} -x -minPctIdentity 70.0 -x -sa {ref}/sequence/{refbase}.fasta.sa --tmpDir=/scratch --regionTable={filtered} "{input}" "{ref}" || exit $?;
echo 'Alignment Complete' || exit $?;
date || exit $?;
loadPulses {input} {output} -metrics QualityValue,InsertionQV,DeletionQV,IPD,PulseWidth -byread || exit $?;
echo 'LoadPulses Complete' || exit $?;
date || exit $?;
cmph5tools.py sort {output} || exit $?;
echo 'Sorting Complete' || exit $?;
date || exit $?;

echo 'Finished on' `date -u`;

# Success
exit 0
""".format(input=input, output=output, ref=ref, refbase=os.path.basename(ref), filtered=filtered, bestn=bestn))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(dest="input_dir", help="input directory")
    parser.add_argument(dest="ref_dir", help="reference directory")
    parser.add_argument("--misc", dest="misc", required=False, help="additional arguments, sep by comma, ex: bestn=1,minLength=20")
    
    args = parser.parse_args()
    
    extra_params = {}
    if args.misc is not None:
        extra_params = dict(stuff.split('=') for stuff in args.misc.split(','))
 
    main(args.input_dir, args.ref_dir, **extra_params)
