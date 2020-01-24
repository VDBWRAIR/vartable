# type: ignore
import csv
from toolz.dicttoolz import merge, dissoc, merge_with, valfilter, keyfilter #done
import vcf #done
from typing import Tuple, Dict, List, Iterator, Iterable, Any, Callable, NamedTuple, BinaryIO, Sequence
import itertools
from functools import partial 
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os.path
from . import translation
from .bam_rc import bam_readcount_pos, BRCRow, BRCEntry
from . import multi_segment

HEADERS = ['Reference ID', 'Position', 'Total Depth', 'Ref Base', 'Alt Base', 'Ref Frequency', 'Alt Frequency', 'Codon', 'Codon Type']
AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
VCFRow = NamedTuple("VCFRow",
                    [('ref', str),
                     ('AO', List[int]),
                     ('DP', int),
                     ('AF', float),
                     ('chrom',str),
                     ('pos', int),
                     ('alt', List[str])])
def flatten_vcf_record(rec):
    # type: (_Record) -> VCFRow
    # AO = allele support depth
    # DP = total raw depth;
    # in freebayes, DP4 gives ref depth (0, 1) then allele support (2, 3). it
    # drops the other bases that don't support the allele; that info is in DP (i
    # believe). but it's weird because the sum of the alt DP4 divided by DP
    # doesn't give the AF.
    if not hasattr(rec.INFO['AF'], '__iter__'):
      rec.INFO['AF'] = [rec.INFO['AF']]
    # round is for formatting to two dec places
    rec.INFO['PAC'] = map(lambda x: round(x*100, 2), rec.INFO['AF'])
    #rec.INFO['PRC'] = map(lambda x: 100 - (x*100), rec.INFO['AF'])
    rec.INFO['PRC'] = round((1 - sum(rec.INFO['AF'])) * 100, 2)
    return rec


def write_tsv(out_path: str, xs: List[Dict[str, Any]], headers: List[str]) -> None:
     headers.sort()
     with open(out_path, 'w') as o:
        writer = csv.DictWriter(o, headers, delimiter='\t')
        writer.writeheader()
        for row in xs:
          writer.writerow(row)

import argparse



    # map over muts

#def fb_to_dict(rec):
#  pac = rec.INFO['AF'] *100
#  prc = 100 - pac
#  return dict(zip(HEADERS[:-2], [rec.CHROM, rec.POS, rec.INFO['DP'], rec.REF, \
#                                 rec.ALT, prc, pac, 'Unavailable', 'Unavailable']))

def muts_to_dicts(muts):
  def to_dict_single(rec, alt, pac):
    return dict(zip(HEADERS[:-2], [rec.CHROM, rec.POS, rec.INFO['DP'], rec.REF, alt, rec.INFO['PRC'], pac, 'Unavailable', 'Unavailable']))
  def to_dict(rec):
    return map(partial(to_dict_single, rec), rec.ALT, rec.INFO['PAC'])
  #dicts = map(to_dict, muts)
  raw_dicts = itertools.chain(*map(to_dict, muts))
  return raw_dicts


def lofreq_process(vcf_path, min_percent, min_depth, out_path):

  def passes(rec: VCFRow) -> bool:
    # VCFRow -> Bool
    #ValueError: dict contains fields not in fieldnames:
    # with base_caller, check rec.ALT[0] because we want to skip those which
    # have no alts and are jus the ref; it's rec.ALT=[None] in that case.
    return rec.ALT[0] and rec.INFO['DP'] >= min_depth and max(rec.INFO['PAC']) >= min_percent

  with open(vcf_path, 'r') as vcf_handle:
    unfiltered_muts = map(flatten_vcf_record, vcf.Reader(vcf_handle))
    muts = filter(passes, unfiltered_muts)

    _dicts = muts_to_dicts(muts)
    def passes2(dict):
      return (dict['Total Depth'] >= min_depth) and (dict['Alt Frequency'] >= min_percent)
    dicts = filter(passes2, _dicts)
    return dicts

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
from typing import Union,Dict
Row = Dict[str, Union[str, int, float]]
RCRow = Dict[str, Union[int, float]]

    
  
def dispatch_args(args: argparse.Namespace) -> Sequence[str]:
  if args.bam:
      assert args.ref, "Reference required with Bam input!"
  if args.type == 'lofreq':
      dicts = lofreq_process(args.vcf_path, args.minp, args.mind, args.out)
  elif args.type == 'base_caller':
      dicts = base_caller_process(args.vcf_path, args.minp, args.mind, args.out)
  if args.cds_rev or args.cds_fwd or args.genbank:
      assert not ((args.cds_rev or args.cds_fwd) and args.genbank), "CDS and genank files cannot be used simultaneously!"
      assert args.ref, "Reference file required when passing CDS as an argument."
      # each seqfeature represents a CDS
      fwd_cdss = [] if not args.cds_fwd else [SeqFeature(location=FeatureLocation(start, end, strand=+1), type='CDS')
         for (start, end) in args.cds_fwd ]
      rev_cdss = [] if not args.cds_rev else [SeqFeature(location=FeatureLocation(start, end, strand=-1), type='CDS')
         for (start, end) in args.cds_rev ]
      with open(args.ref) as ref_file:
          refs = list(SeqIO.parse(ref_file, format='fasta'))
          assert len(refs) == 1, "Only one reference sequence currently supported."
          rec = refs[0]
          rec.features = sorted(fwd_cdss + rev_cdss, key=lambda x: x.location._start)
      if args.genbank:
          with open(args.genbank) as genbank_file:
              genbanks = list(SeqIO.parse(genbank_file, format='genbank'))
              assert len(genbanks) == 1, "Only one genbank record currently supported."
              rec = genbanks[0]
      variants = [ (d['Position'], d['Alt Base']) for d in dicts ]
      translation_results = translation.dispatch(rec, variants)
      # Note: dicts will have None values, but "None" is acceptable in the output TSV
      dicts = [d.update(tr.__dict__) for d, tr in zip(dicts, translation_results) ]
      dupe_fields = ('position', 'alt')
      for d in dicts:
          assert d['position'] == d['Position'], "Ordering of translation return values went bad."
          # delete fields that are now duplicated.
          for k in dupe_fields:
              del d[k]
              del d[k]
      fields = list(translation.TResult.__dataclass_fields__.keys())
      fields = [f for f in fields if not (f in dupe_fields)]
      global HEADERS
      HEADERS = HEADERS[:] + fields

  if args.bam:
      listdicts = list((d for d in dicts if d['Alt Base'] != '*'))
      ref_id = listdicts[0]['Reference ID'] # TODO: multiple refs in a single vcf
      get_info = partial(bam_readcount_pos, args.bam, args.ref, args.mind, ref_id)
      bam_info = map(get_info, (x['Position'] for x in listdicts))
      def to_dict(v_d: Dict[str, Any], brc: BRCRow) -> Dict[str, Any]:
          entry = brc.entries[ str(v_d['Alt Base']) ]
          #entry = brc.entries['='] # not reliable
          # below should actually be type safe
          return {k : getattr(entry, k) for k in entry._fields if k != 'base' }
      brc_dicts =  map(to_dict, listdicts, bam_info)
      dicts = map(merge, listdicts, brc_dicts)
      HEADERS = ['Reference ID', 'Position', 'Total Depth', 'Ref Base', 'Alt Base', 'Ref Frequency', 'Alt Frequency', 'Codon', 'Codon Type']
      HEADERS = HEADERS[:] + [ f for f in BRCEntry._fields if f != 'base'] 
  write_tsv(args.out, dicts, HEADERS[:])
  return HEADERS[:]
# PAC is a lit, because will report multiple Alleles.  i.e.
# AY487947.1    142     .       A       C,G     .       .       DP=2087;RC=2083;RAQ=37;PRC=100;AC=1,3;AAQ=37,37;PAC=0,0;CBD=2083;CB=A
# Flu is a special case again. have to groupby, map, and flatten; should be an
# exmaple of this in . . . . no, actually, just print the reference ID

def main() -> None:
  parser = argparse.ArgumentParser(description='Create a table from VCF with summary info. ')
  parser.add_argument('vcf_path', help='VCF input, from lofreq or ngs_mapper.base_caller')
  parser.add_argument('--type', dest='type', required=True, choices=['lofreq', 'base_caller'], help='Was this created with lofreq or standard ngs_mapper?')
  parser.add_argument('--out',  dest='out', required=True,  help='output TSV file')
  parser.add_argument('--mindepth', dest='mind', type=int, required=True, help='Minimum depth to filter out variants.')
  parser.add_argument('--minpercent', dest='minp', type=int, required=True, help='Minimum percentage as an integer to filter out variants.')
  parser.add_argument('--bam', dest='bam', required=False, help='Optional indexed BAM input.')
  parser.add_argument('--ref', dest='ref', required=False, help='Reference fasta file required for BAM processing.')
  parser.add_argument('--genbank', dest='genbank', required=False, help='Optional genbank file with CDS feature info.')
  parser.add_argument('--cds-rev', dest='cds_rev', action='append', metavar=('start', 'end'), nargs=2, required=False, help='Optional: zero-based coordinates for a simple cds translated on the reverse strand. Requires reference.')
  parser.add_argument('--cds-fwd', dest='cds_fwd', action='append', metavar=('start', 'end'), nargs=2, required=False, help='Optional: zero-based coordinates for a simple cds translated on the forward strand. Requires reference.')
  parser.add_argument('--multi-segment', dest='multi_segment', default=False, action='store_true')
  parser.add_argument('--allow-overwrite', dest='allow_overwrite', default=False, action='store_true')
  args = parser.parse_args()
  multi_segment.validate_file_overwrites([args])
  if args.multi_segment:
      return multi_segment.dispatch(args, dispatch_args, write_tsv)
  return dispatch_args(args)

if __name__ == '__main__':
    main()
