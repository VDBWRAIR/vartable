import vcf
from Bio import SeqIO
# builtin libs
from typing import Any, List, Iterator, Tuple, NoReturn, Dict, Callable, Union, Sequence
from itertools import groupby, starmap, chain
import re
import argparse
import copy
import csv
import os
import sys
# for types
from Bio.SeqRecord import SeqRecord

def raise_(e: Exception) -> NoReturn:
   raise e

def sanitize(fn: str) -> str: return re.sub('[^\w\-_\. ]', '_', fn)
debug_log = sys.stderr.write

def write_subfiles(vcf_reader: vcf.parser.Reader, seqio_refs: dict, refid: str, vcf_grp: Iterator[Any]) -> Tuple[str, str, str, str]:
   refseq: Dict[str, SeqRecord] =  seqio_refs.get(refid) or raise_(ValueError(f"{refid} not found inreference fasta file"))
   basename = f"{sanitize(vcf_reader.filename)}.{sanitize(refid)}"
   vcf_name = f"{basename}.vcf"
   ref_name = f"{basename}.fasta"
   vartable_fn = f"{basename}.tsv"
   with open(vcf_name, 'w') as vcf_out_file:
       vcf_writer = vcf.Writer(vcf_out_file, vcf_reader)
       for r in vcf_grp: vcf_writer.write_record(r)
   with open(ref_name, 'w') as ref_out_file:
       SeqIO.write([refseq], ref_out_file, 'fasta')
   return (refid, ref_name, vcf_name, vartable_fn)


def make_args(in_args: argparse.Namespace, chrom: str, ref_fn: str, vcf_fn: str, vartable_out_fn: str) -> argparse.Namespace:
    ns = copy.deepcopy(in_args)
    ns.vcf_path = vcf_fn
    ns.ref = ref_fn
    ns.out = vartable_out_fn
    return ns 
    #subprocess.check_call('vartable', '--bam', args.bam, '--vcf', vcf_fn, '--ref'
def chrom_(r: vcf.model._Record) -> str:
   return r.CHROM
# def read_tsv
DEBUG = False
COMPILE_SUB_RESULTS = True
def validate_file_overwrites(arg_sets: Sequence[argparse.Namespace]) -> Union[None,NoReturn]:
    for arg_set in arg_sets:
       problems=[]
       if (not arg_set.allow_overwrite) and os.path.exists(arg_set.out):
           problems.append( (f"--allow-overwrite option is off, but file {arg_set.out} exists! Stopping.") )
    if any(problems):
       raise ValueError("\n".join(problems))
    return None
   
def dispatch(args: argparse.Namespace, vartable_dispatch: Callable[[argparse.Namespace], int], writer: Callable) -> int:
    vcf_fn, ref_fn = args.vcf_path, args.ref
    with open(ref_fn) as ref_file, open(vcf_fn) as vcf_file:
        refs = { x.id : x for x in SeqIO.parse(ref_file, 'fasta') } 
        vcf_reader = vcf_reader = vcf.Reader(vcf_file)
        vcf_groups = groupby(vcf_reader, key=chrom_) # lambda x: x.CHROM)
        chrom_files =  [write_subfiles(vcf_reader, refs, chrom, grp) for chrom, grp in vcf_groups]
    subrun_args = [make_args(args, *info) for info in chrom_files]
    validate_file_overwrites(subrun_args)
    for arg_set in subrun_args:
        debug_log(f"Running: {arg_set}")
        if not DEBUG:
            headers = vartable_dispatch(arg_set)
            # debug_log(result)
    if COMPILE_SUB_RESULTS:
        readers = [csv.DictReader(open(args.out), delimiter='\t') for args in subrun_args ]
        headers.sort()
        sub_headers = [ sorted(r.fieldnames) for r in readers ]
        assert all( (x == headers for x in sub_headers)  ), f"Header missmatch of {headers} of output csv files, see:\n {subrun_args}\n{sub_headers}"
        all_rows = chain.from_iterable( readers )
        writer(args.out, all_rows, headers)
    return 0
