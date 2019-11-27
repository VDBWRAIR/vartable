# from typing import * 
from dataclasses import dataclass
# from vartable.types import find
from typing import Sequence, Iterable
def first(seq):
    '''@Warning: mutates seq'''
    if isinstance(seq, Sequence):
        return None if len(seq) == 0 else seq[0]
    else:
        try:
            return next(seq)
        except StopIteration:
            return None
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
import pysam
from typing import Iterator, Any, List
AlignedSegment = AlignmentFile = Any

def single_bam_row_info(alt_position: int, thresh: int, row: AlignedSegment) -> BamRowInfo:
  # get query position
  query_position = calc_query_position(alt_position, row.reference_start, row.query_alignment_start)
  base_at_alt_position = row.query_alignment_sequence[query_position]
  # NOTE: Not worrying if soft-clipping would interfere because it's not part of manual process
  # could check for soft clipping by looking at the first/last entry of the cigar string
  not_primer_start = thresh
  not_primer_end = len(row.query_alignment_sequence) - thersh
  is_within_end = (query_position < not_primer_start) or (query_position > not_primer_end)
  return BamRowInfo(base_at_alt_position, is_within_end)

import functools 
def primer_count_at_position(thresh: int, reference_id: str, primerCoords: List[MappingCoordinate], sampleBamFile: str, alt_position: int) -> PrimerResult: # , rows: Iterator[AlignedSegment]) -> PrimerCount:
# TODO1 : Could have more than one primer overlapping, which would be a problem. 
  def overlaps(coord: MappingCoordinate) -> bool: 
    return (coord.reference_start <= alt_position) and (alt_position <= coord.reference_end)
  overlapped_primer: Optional[MappingCoordinate] = find(overlaps, primerCoords)
  if not overlapped_primer:
    return OutsidePrimer()
  with pysam.AlignmentFile(bamfile) as alnfile:
    bam_rows = alnfile.fetch(reference=reference_id, start=alt_position, stop=alt_position+1)
  single_info = functools.partial(single_bam_row_info, alt_position=alt_position, thresh=thresh)
  row_infos = map(single_info, bam_rows)
  # TODO: Does this assume that the base 
  primer_base = coord.query_alignment_sequence[ (alt_position - coord.reference_start) ]
  return reduce_bamRowInfo(row_infos, alt, coord.base)

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
  return PrimerCount(alt_position, alt, total, totalWithinEnd, supportsAltTotal, supportsAltWithinEnd, primer_base)
      
      
        
   









   
  
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

 #   query_sequence read sequence bases, including soft clipped bases (None if not present).

global res
if __name__ == '__main__':
  refid = 'Den1/U88535_1/WestPac/1997/Den1_1' 
  sam = 'testdata/primers.sam'
  res = get_coordinates_from_bam(refid, sam)
  print(res)








#def pileup(bam: AlignmentFile, ref: str, start: int, stop: int, opts: PileupOptions) -> Iterator[PileupColumn]:
#    assert ref in bam.references
#    assert start >= 0 and stop <= bam.get_reference_length(ref) 
#    assert bam.check_index()
#    assert not bam.closed
#    assert bam.is_bam
#
#    BIG_INT = 999_999_999
#
#    for x in bam.pileup(ref, start, stop, max_depth=BIG_INT,
#               min_base_quality=opts.minBaseQual,
#               ignore_orphans=opts.ignoreOrphans,
#               ignore_overlaps=opts.ignoreOverlaps):
#        yield x 




def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead], thresh: int) -> List[PrimerCount]: 
    pass 

    with pysam.AlignmentFile(bampath) as bam:
        for ref in bam.references:
            pileupgen = pileup(bam, ref, 0, bam.get_reference_length(ref), opts)
