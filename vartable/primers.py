from dataclasses import dataclass
import pysam
from typing import Iterator, Any, List, Optional
import vcf
# from vartable.types import find
from typing import Sequence, Iterable
def first(seq):
    '''@Warning: mutates seq'''
    if isinstance(seq, Sequence):
        ret = None if len(seq) == 0 else seq[0]
    else:
        try:
            ret = next(seq)
        except StopIteration:
            ret = None
    return ret

def find(p, seq):
    return first(filter(p, seq))

@dataclass
class Alt:
    position: int
    base: str

@dataclass
class MappingCoordinate:
    reference_start: int
    reference_end: int
    query_alignment_sequence: str

@dataclass
class BamRead:
    referenceEnd: int
    referenceStart: int
    queryStart: int
    queryEnd: int
    querySequence: str # without soft clipping
    querySequenceSoftClipped: str 
  #  def query_base_at_reference_position(position: int) -> str: ....
# drawing from https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment 
@dataclass
class PrimerResult: ...
class OutsidePrimer(PrimerResult):
  outsidePrimer = True

class PrimerCount(PrimerResult):
   outsidePrimer = False
   position: int
   alt: Alt # position and alternate where this count came from. the alt call is obtained from the VCF file.
   total: int
   totalWithinEnd : int
   supportsAltTotal: int # won't have this information for non-primer regions
   supportsAltWithinEnd: int
   primer_base: str


@dataclass
class BamRowInfo:
    queryBase: str
    withinEnd: bool

  # alnfile = pysam.AlignmentFile('testdata/fullsample.bam')
AlignedSegment = AlignmentFile = Any

def single_bam_row_info(alt_position: int, thresh: int, row: AlignedSegment) -> BamRowInfo:
  # get query position
  query_position = calc_query_position(alt_position, row.reference_start, row.query_alignment_start)
  base_at_alt_position = row.query_alignment_sequence[query_position]
  # NOTE: Not worrying if soft-clipping would interfere because it's not part of manual process

  # could check for soft clipping by looking at the first/last entry of the cigar string
  not_primer_start = thresh
  not_primer_end = len(row.query_alignment_sequence) - thresh
  is_within_end = (query_position < not_primer_start) or (query_position > not_primer_end)
  return BamRowInfo(base_at_alt_position, is_within_end)

import functools 
def primer_count_at_position(thresh: int, primerCoords: List[MappingCoordinate], sampleBamFile: str, reference_id: str, alt: Alt) -> PrimerResult: 
# , rows: Iterator[AlignedSegment]) -> PrimerCount:
# TODO1 : Could have more than one primer overlapping, which would be a problem.
  def overlaps(coord: MappingCoordinate) -> bool: 
    return (coord.reference_start <= alt.position) and (alt.position <= coord.reference_end)

  overlapped_primer: Optional[MappingCoordinate] = find(overlaps, primerCoords)
  if not overlapped_primer:
    return OutsidePrimer()
  def get_base_of_primer_at_position(coord: MappingCoordinate) -> str:
      return coord.query_alignment_sequence[ (alt.position - coord.reference_start) ]
  from collections import Counter
  primer_cnt = Counter(map(get_base_of_primer_at_position, primerCoords))
  primer_base = primer_cnt.most_common()[0] # could log here
  with pysam.AlignmentFile(sampleBamFile) as alnfile:
    bam_rows = alnfile.fetch(reference=reference_id, start=alt.position, stop=alt.position+1)
  single_info = functools.partial(single_bam_row_info, alt_position=alt.position, thresh=thresh)
  row_infos = map(single_info, bam_rows)
  # TODO: Does this assume that the base
  total = totalWithinEnd = supportsAltWithinEnd = supportsAltTotal = 0
  for bri in row_infos:
      isAlt = bri.queryBase == alt.base
      total += 1
      totalWithinEnd += int(bri.withinEnd)
      supportsAltTotal += int(isAlt)
      supportsAltWithinEnd += int(bri.withinEnd and isAlt)
  pc = PrimerCount(alt.position, alt, total, totalWithinEnd, supportsAltTotal, supportsAltWithinEnd, primer_base)
  return pc

def primer_count_all_reference(thresh: int, primerCoords: List[MappingCoordinate], sampleBamFile: str, reference_id: str, alts: Iterable[Alt]) -> Iterator[PrimerResult]:
  primer_single_position = partial(primer_count_at_position, thresh, primerCoords, reference_id)
  return map(primer_single_position, alts)

def main(thresh: int, primerFile: str, sampleBam: str, alts: Iterable[Alt], reference_file: str, reference_id: str) -> Iterator[PrimerResult]:
  primer_bam_name = f'{primerFile}.sam'
  bwa_map(primerFile, reference_file, primer_bam_name)
  primer_coords =  get_coordinates_from_bam(reference_id, primer_bam_name)
  result = primer_count_all_reference(thresh, primer_coords, sampleBam, reference_id, alts)
  return result
  
def reduce_bamRowInfo(infos: Iterable[BamRowInfo], alt: str, primer_base: str, alt_position: int, alt_primer_base: str) -> PrimerResult:
  total, altTotal, totalWithinEnd, supportsAltWithinEnd = 0, 0, 0, 0
  for x in infos:
    total += 1
    if x.queryBase == alt:
       altTotal += 1
       if x.withinEnd: 
         supportsAltWithinEnd += 1
    elif x.withinEnd:
       totalWithinEnd += 1
  return PrimerCount(alt_position, alt, total, totalWithinEnd, altTotal, supportsAltWithinEnd, primer_base)
      
      
  # query_alignment_sequence starts at query_start [query_sequence does not]

import os
import sys
import subprocess


def bwa_map(query_file: str, ref_file: str, out_sam: str) -> None:
  # @Side-effect: produce index files for ref_file if they don't exist
  # @Side-effect: produces SAM file as out_file

  index_files = '{ref_file}.amb' ,'{ref_file}.bwt' ,'{ref_file}.pac' ,'{ref_file}.ann' ,'{ref_file}.sa'
  if not all(map(os.path.exists, index_files)):
     subprocess.check_call(['bwa', 'index', ref_file])
     #index_proc_out = subprocess.check_output(['bwa', 'index', ref_file]); #sys.stderr.write(index_proc_out)
  # bwa is okay with fasta query inputs it seems 
  with open(out_sam, 'w') as out:
     subprocess.check_call(['bwa', 'mem', ref_file, query_file], stdout=out) 
  return None

#  for idx in index_files:
#    assert os.path.exists(idx), 
def get_coordinates_from_bam(refid: str, bamfile: str) -> List[MappingCoordinate]:
  with pysam.AlignmentFile(bamfile) as alnfile:
    coordinates = [ MappingCoordinate(aln.reference_start, aln.reference_end, aln.query_alignment_sequence) for aln in alnfile if aln.reference_name == refid]
  return coordinates

def calc_query_position( alt_position :int, ref_start: int, query_start: int ): 
  return alt_position - ref_start

 #  query_sequence read sequence bases, including soft clipped bases (None if not present).

def base_caller_process(vcf_path, min_percent, min_depth, out_path):
  def passes(rec):
    # VCFRow -> Bool
    #ValueError: dict contains fields not in fieldnames:
    # with base_caller, check rec.ALT[0] because we want to skip those which
    # have no alts and are jus the ref; it's rec.ALT=[None] in that case.
    return rec.ALT[0] and rec.INFO['DP'] >= min_depth and max(rec.INFO['PAC']) >= min_percent
  with open(vcf_path, 'r') as vcf_handle:
    #TODO: filter on individual alts, so after flattening (I guess as dicts var, in dicts form)
    muts = filter(passes, vcf.Reader(vcf_handle))
  #  def to_dict(rec):
  #    return dict(zip(HEADERS[:-2], [rec.CHROM, rec.POS, rec.REF, rec.ALT, rec.INFO['PRC'], rec.INFO['PAC'], 'Unavailable', 'Unavailable']))
    def to_dict_single(rec, alt, pac):
      return dict(zip(HEADERS[:-2], [rec.CHROM, rec.POS, rec.INFO['DP'], rec.REF, alt, rec.INFO['PRC'], pac, 'Unavailable', 'Unavailable']))
    def to_dict(rec):
      return map(partial(to_dict_single, rec), rec.ALT, rec.INFO['PAC'])
    #dicts = map(to_dict, muts)
    def passes2(dict):
      return (dict['Total Depth'] >= min_depth) and (dict['Alt Frequency'] >= min_percent)
    raw_dicts = itertools.chain(*map(to_dict, muts))
    dicts = filter(passes2, raw_dicts)
    return dicts
global res
# from vartable import base_caller_process
if __name__ == '__main__':
  refid = 'Den1/U88535_1/WestPac/1997/Den1_1' 
  sam = 'testdata/primers.sam'
  res = get_coordinates_from_bam(refid, sam)
  print(res)
  bam = '/vdb/instr/novaseq/01/Share/for_mike/Projects/A18x3210/A18x3210.bam'
  vcf = '/vdb/instr/novaseq/01/Share/for_mike/Projects/A18x3210/A18x3210.bam.vcf'
  ref = '/vdb/instr/novaseq/01/Share/for_mike/Projects/A18x3210/JEV_XZ0934.fasta'
  primer = '/vdb/instr/novaseq/01/Share/for_mike/Projects/primer.fasta'
  reference_id = 'JF915894.1'
  dicts = base_caller_process(vcf, 10, 10, '')
  alts = list((d for d in dicts if d['Alt Base'] != '*'))
  ref_id = listdicts[0]['Reference ID'] # TODO: multiple refs in a single vcf
  results = main(20, primer, bam, alts, ref, ref_id)
  print(results)
  # main(thresh: int, primerFile: str, sampleBam: str, alts: Iterable[Alt], reference_file: str, reference_id: str) -> Iterator[PimerResult]: 


def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead], thresh: int) -> List[PrimerCount]: 
    pass 

    with pysam.AlignmentFile(bampath) as bam:
        for ref in bam.references:
            pileupgen = pileup(bam, ref, 0, bam.get_reference_length(ref), opts)
