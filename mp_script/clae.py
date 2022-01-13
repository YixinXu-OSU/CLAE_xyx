from env_check import env_check
from blast_reader import read_blastn_and_divide
from seq_reader import read_seqs
from subread import quality_trimming, seq_extraction, seq_extraction_no_blat
from consensus_ref import consensus_finding_sparc, consensus_finding_pbdagcon, consensus_finding_sparc_lseq
from consensus_no_ref import no_ref_consensus_finding_chris, no_ref_consensus_finding_pbdagcon, no_ref_consensus_finding_sparc, no_ref_consensus_finding_sparc_lseq
from merge import merge
import logging, time, os
from datetime import datetime


import sys, os, argparse, warnings
from multiprocessing import Pool

def main():
    warnings.filterwarnings("ignore")
    
    
    
    parser = argparse.ArgumentParser('Subread Extractor and Consensus Sequence Generator')
    parser.add_argument('--blast', type=str, help='Blast File Name .csv: ', required=True)
    parser.add_argument('--comma', help='Include to use comma (,) seperated blast result instead of tab seperated blast', action='store_true')
    parser.add_argument('--seq', help='Sequence Filename, must be fastq if run phred score.', type=str, required=True)
    parser.add_argument('--llen', help='Specify Llen', type=int, default=160)
    parser.add_argument('--ref', help='Include to use reference mode', action='store_true')
    parser.add_argument('--start', help='Start Chunk ID, inclusive', type=int)
    parser.add_argument('--end', help='End Chunk ID, inclusive', type=int)
    parser.add_argument('--refseq', help='Ref Sequence Filename', type=str, required='--ref' in sys.argv)
    parser.add_argument('--algo', help='Consensus Tool Selection, s for Sparc, p for pbdagcon, c for chris (c invalid for Non-ref mode)', type=str, required=True)
    parser.add_argument('--merge', help='Include to merge autodivided chunks. RECOMMEND TO INCLUDE', action='store_true')
    # parser.add_argument('--phred', help='Include to enable phred score filtering.', action='store_true')
    parser.add_argument('--phredthreshold', help='Phred Score threshold for consensus generation.', type=float)
    # parser.add_argument('--seqfastq', help='Fastq of sequencing, required for phred score generation.', type=str, required='--phred' in sys.argv)
    parser.add_argument('--historydir', help='Directory of historical data to reuse.', type=str)
    parser.add_argument('--runlseq', help='Include to run lseq consensus generation.', action='store_true')
    parser.add_argument('--lseqref', help='Reference file for Lseq')
    parser.add_argument('--strand', help='plus/minus/both')

    args = parser.parse_args()
    
    root_dir = os.getcwd()
    mugio_path = os.path.join(root_dir, 'mugio.py')
    
    exec_dirname = ''
    
    if args.historydir is not None:
        exec_dirname = args.historydir
    else:
        exec_dirname = f'{args.seq}_{datetime.now().strftime("%Y%m%d-%H%M%S")}_files'
        if not os.path.exists(exec_dirname):
            os.makedirs(exec_dirname)
    os.chdir(exec_dirname)
    
    logging.basicConfig(filename='consensus.log', level=logging.DEBUG)
    
    logging.info("==================================================")
    logging.info(">>> Running Log <<<")
    logging.info("==================================================")
    t1_start = time.perf_counter()
    t2_start = time.process_time()

    sep = '\t'
    if args.comma:
        sep = ','

    env_check()
    chunks = read_blastn_and_divide(os.path.join(root_dir, args.blast), sep)

    start = 0
    end = chunks

    if args.start is not None:
        start = args.start
    if args.end is not None:
        end = args.end + 1

    seqs = read_seqs(os.path.join(root_dir, args.seq))
    
    pool = Pool()

    # Quality Trimming

    for i in range(start, end):
        pool.apply_async(quality_trimming, args=(args.llen, seqs, i))
        
    pool.close()
    pool.join()

    # Sequence Extraction
    
    fastq_filename = '' if '--phredthreshold' not in sys.argv else os.path.join(root_dir, args.seq)
    phred_threshold = 0 if '--phredthreshold' not in sys.argv else args.phredthreshold
    strand_filter = 'both' if '--strand' not in sys.argv else args.strand
    
    # phred
    
    phred_enabled = False if ('--phredthreshold' not in sys.argv) or (args.phredthreshold == 0) else True
    
    
    '''
    for i in range(start, end):
        if 'c' not in args.algo and not args.phred:
            seq_extraction_no_blat(i)
        else:
            seq_extraction(i, args.phred, fastq_filename, phred_threshold, mugio_path, strand_filter)
    
    '''
    
    pool = Pool()

    for i in range(start, end):
        if 'c' not in args.algo and not phred_enabled and (strand_filter != 'minus' and strand_filter != 'plus'):
            pool.apply_async(seq_extraction_no_blat, args=(i,))
        else:
            pool.apply_async(seq_extraction, args=(i, phred_enabled, fastq_filename, phred_threshold, mugio_path, strand_filter))
            
    pool.close()
    pool.join()

    
    

    # Consensus Finding
    
    pool = Pool()

    for i in range(start, end):

        # Target sequence consensus

        # Ref mode
        if args.ref:
            if 's' in args.algo:
                pool.apply_async(consensus_finding_sparc, args=(i, os.path.join(root_dir, args.refseq)))
            if 'p' in args.algo:
                pool.apply_async(consensus_finding_pbdagcon, args=(i, os.path.join(root_dir, args.refseq)))
        # No ref mode
        else:
            if 's' in args.algo:
                pool.apply_async(no_ref_consensus_finding_sparc, args=(i,))
            if 'p' in args.algo:
                pool.apply_async(no_ref_consensus_finding_pbdagcon, args=(i,))
            if 'c' in args.algo:
                pool.apply_async(no_ref_consensus_finding_chris, args=(i,))
                
                
                
        # Lseq consensus
        
        if args.runlseq:
            if '--lseqref' in sys.argv:
                pool.apply_async(consensus_finding_sparc_lseq, args=(i, os.path.join(root_dir, args.lseqref)))
            else:
                pool.apply_async(no_ref_consensus_finding_sparc_lseq, args=(i,))
                    

    pool.close()
    pool.join()
    
    '''
    for i in range(start, end):

        # Ref mode
        if args.ref:
            if 's' in args.algo:
                consensus_finding_sparc(i, args.refseq)
            if 'p' in args.algo:
                consensus_finding_pbdagcon(i, args.refseq)
        # No ref mode
        else:
            if 's' in args.algo:
                no_ref_consensus_finding_sparc(i)
            if 'p' in args.algo:
                no_ref_consensus_finding_pbdagcon(i)
            if 'c' in args.algo:
                no_ref_consensus_finding_chris(i)
    '''
    
    # clean files generated by blasr
    for f in [os.path.join(dp, fi) for dp, dn, fn in os.walk(os.path.expanduser(root_dir)) for fi in fn]:
        if 'core.' in f or '.psl' in f:
            os.remove(f)
    
    
    if args.merge:
        merge(start, end, args.ref, args.algo, args.seq, args.runlseq, '--lseqref' in sys.argv)
        
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    logging.info("--------------------------------------------------")
    logging.info("Elapsed time: %.1f [min]" % ((t1_stop-t1_start)/60))
    logging.info("CPU process time: %.1f [min]" % ((t2_stop-t2_start)/60))
    logging.info("--------------------------------------------------")
    

if __name__ == "__main__":
    main()




