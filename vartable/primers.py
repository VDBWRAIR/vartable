from dataclasses import dataclass
import pysam
import os
from typing import Iterator, Any, List, Optional, NamedTuple
import vcf
from functools import partial
import itertools
# from vartable.types import find
import sys
from typing import Sequence, Iterable

debug_log = sys.stderr.write
TESTING = True
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

class Alt(NamedTuple):
    position: int
    base: str

class MappingCoordinate(NamedTuple):
    reference_start: int
    reference_end: int
    query_alignment_sequence: str
    query_alignment_start: int

class BamRead(NamedTuple):
    referenceEnd: int
    referenceStart: int
    queryStart: int
    queryEnd: int
    querySequence: str # without soft clipping
    querySequenceSoftClipped: str 
  #  def query_base_at_reference_position(position: int) -> str: ....
# drawing from https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
import abc
@dataclass
class PrimerResult(metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def outsidePrimer(self) -> bool: ...

@dataclass
class PrimerCount(PrimerResult):
   position: int
   alt: Alt # position and alternate where this count came from. the alt call is obtained from the VCF file.
   total: int
   totalWithinEnd : int
   supportsAltTotal: int # won't have this information for non-primer regions
   supportsAltWithinEnd: int
   primer_base: str

   @property
   def outsidePrimer(self) -> bool:
       return False

class OutsidePrimer(PrimerResult):
   @property
   def outsidePrimer(self) -> bool:
       return True

@dataclass
class BamRowInfo:
    queryBase: str
    withinEnd: bool

  # alnfile = pysam.AlignmentFile('testdata/fullsample.bam')
AlignedSegment = AlignmentFile = Any

def calc_query_position(row: Any, alt_position :int) -> Optional[int]: #, query_start: int ): 
  query_pos_and_ref_pos = find(lambda x: x[1] == alt_position, row.get_aligned_pairs())
  if query_pos_and_ref_pos and query_pos_and_ref_pos[0]:
      assert query_pos_and_ref_pos[1] in (alt_position, None)
      return query_pos_and_ref_pos[0] - row.query_alignment_start
def single_bam_row_info(alt_position: int, thresh: int, row: AlignedSegment) -> Optional[BamRowInfo]:
  # get query position
        
    # return alt_position - ref_start
# (Pdb) row.get_aligned_pairs()
  query_position = calc_query_position(row, alt_position) # , row.reference_start) # , row.query_alignment_start)
  if not query_position: # the query may be deleted or something at that position
      return None
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

  overlapped_primers: Optional[MappingCoordinate] = filter(overlaps, primerCoords)
  if not overlapped_primers:
    return OutsidePrimer()
  def get_base_of_primer_at_position(coord: MappingCoordinate) -> str:
      query_position = calc_query_position(coord, alt.position)
      if not query_position:
          return '*' # or should it be a? IDK
      return coord.query_alignment_sequence[ query_position ] 
  from collections import Counter
  primer_cnt = Counter(map(get_base_of_primer_at_position, overlapped_primers))
  primer_base = primer_cnt.most_common()[0] # could log here
  # import pdb; pdb.set_trace()
  with pysam.AlignmentFile(sampleBamFile) as alnfile:
    bam_rows = alnfile.fetch(reference=reference_id, start=alt.position, stop=alt.position+1)
    single_info = functools.partial(single_bam_row_info,alt.position, thresh)
    row_infos = filter(bool, map(single_info, bam_rows))
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
  primer_single_position = partial(primer_count_at_position, thresh, primerCoords, sampleBamFile, reference_id)
  return map(primer_single_position, alts)

def main(thresh: int, primerFile: str, sampleBam: str, alts: Iterable[Alt], reference_file: str, reference_id: str) -> Iterator[PrimerResult]:
  primer_bam_name = f'{primerFile}.sam'
  if not TESTING: bwa_map(primerFile, reference_file, primer_bam_name)
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
  index_files = f'{ref_file}.amb' , f'{ref_file}.bwt' , f'{ref_file}.pac' ,f'{ref_file}.ann' , f'{ref_file}.sa'
  if not all(map(os.path.exists, index_files)):
     debug_log(f"could not find all required index files: {index_files}\nRuning BWA index on {ref_file}")
     subprocess.check_call(['bwa', 'index', ref_file])
     #index_proc_out = subprocess.check_output(['bwa', 'index', ref_file]); #sys.stderr.write(index_proc_out)
  # bwa is okay with fasta query inputs it seems 
  with open(out_sam, 'w') as out:
     debug_log(f"Runing BWA mem on {ref_file}")
     debug_log(f"Running [bwa mem {ref_file} {query_file} > {out}")
     subprocess.check_call(['bwa', 'mem', ref_file, query_file], stdout=out)
     debug_log(f"bwa mem finished.\n")
  return None

#  for idx in index_files:
#    assert os.path.exists(idx), 
def get_coordinates_from_bam(refid: str, bamfile: str) -> List[MappingCoordinate]:
  with pysam.AlignmentFile(bamfile) as alnfile:
    coordinates = [ aln for aln in alnfile if aln.reference_name == refid]
    # coordinates = [ MappingCoordinate(aln.reference_start, aln.reference_end, aln.query_alignment_sequence, aln.query_alignment_start) for aln in alnfile if aln.reference_name == refid]
  return coordinates


 #  query_sequence read sequence bases, including soft clipped bases (None if not present).

def base_caller_process(vcf_path, min_percent, min_depth, out_path):
  def passes(rec):
    # VCFRow -> Bool
    #ValueError: dict contains fields not in fieldnames:
    # with base_caller, check rec.ALT[0] because we want to skip those which
    # have no alts and are jus the ref; it's rec.ALT=[None] in that case.
    return rec.ALT[0] and rec.INFO['DP'] >= min_depth and max(rec.INFO['PAC']) >= min_percent
  with open(vcf_path, 'r') as vcf_handle:
    HEADERS = ['Reference ID', 'Position', 'Total Depth', 'Ref Base', 'Alt Base', 'Ref Frequency', 'Alt Frequency', 'Codon', 'Codon Type']

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
# global res
# from vartable import base_caller_process
if __name__ == '__main__':
  refid = 'Den1/U88535_1/WestPac/1997/Den1_1' 
  primer_sam = 'testdata/primers.sam'
  # print(res)
  ref = 'testdata/Den1__WestPac__1997.fasta'
  primer_fasta = 'testdata/read1.fasta'
  # if not os.path.exists(primer_sam): bwa_map(primer_fasta, ref, primer_sam)
  res = get_coordinates_from_bam(refid, 'testdata/primers.sam')
  # primer_sa_print(res)
  bam = 'testdata/fullsample.bam'
  vcf_path = 'testdata/fullsample.bam.vcf'
  # reference_id = 'JF915894.1'
  reference_id = 'Den1/U88535_1/WestPac/1997/Den1_1'
  dicts = base_caller_process(vcf_path, 10, 10, '')
  raw_alts_ = list((d for d in dicts if d['Alt Base'] != '*'))
  ref_id_ = raw_alts_[0]['Reference ID'] # TODO: multiple refs in a single vcf
  alts_ = [Alt(d['Position'], d['Alt Base']) for d in raw_alts_]
  results = main(20, 'testdata/primers', bam, alts_, ref, ref_id_)
  print(list(results))
  # main(thresh: int, primerFile: str, sampleBam: str, alts: Iterable[Alt], reference_file: str, reference_id: str) -> Iterator[PimerResult]: 


# def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead], thresh: int) -> List[PrimerCount]: 
#     pass 
# 
#     with pysam.AlignmentFile(bampath) as bam:
#         for ref in bam.references:
#             pileupgen = pileup(bam, ref, 0, bam.get_reference_length(ref), opts)
