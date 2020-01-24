from typing import NewType, List, Union, NamedTuple

from typing_extensions import Literal, TypedDict
from Bio.Seq import Seq
PRINT_BAMRC_TOOL = True
Fraction = NewType('Fraction', float)

A = Literal['A']
C = Literal['C']
T = Literal['T']
G = Literal['G']
N = Literal['N']
E = Literal['=']
Nuc = Union[A, C, T, G]
NucN = Union[Nuc, N]
NucEN = Union[NucN, E]

# could be a typeddict or an List[5] type . . . note the = on the first thing i.e. Type
class BRCEntry(NamedTuple):
  base: Nuc
  count_: int
  avg_mapping_quality: float
  avg_basequality: float
  avg_se_mapping_quality: float
  num_plus_strand: int
  num_minus_strand: int
  avg_pos_as_fraction: Fraction
  avg_num_mismatches_as_fraction: float
  avg_sum_mismatch_qualities: float
  num_q2_containing_reads: int
  avg_distance_to_q2_start_in_q2_reads: Fraction
  avg_clipped_length: float
  avg_distance_to_effective_3p_end: Fraction

BRCEntries = TypedDict('BRCEntries', { '=' : BRCEntry, 'A' : BRCEntry, 'C' : BRCEntry, 'G' : BRCEntry, 'T' : BRCEntry, 'N' : BRCEntry })
class BRCRow(NamedTuple):
  chrom: str
  pos: int
  ref_base: Nuc
  depth: int
  entries: BRCEntries

s = '''Den1/U88535_1/WestPac/1997/Den1_1	3	T	2	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:2:60.00:34.50:0.00:1:1:0.00:0.00:34.50:1:0.99:225.00:0.50	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00'''

from typing import cast,Any,Callable,Iterator
import subprocess
import io
def bam_readcount_pos(bam: str, fasta: str, mindepth: int, ref_id: str, pos: int) -> BRCRow:
    pos_str = "{0}:{1}-{1}".format(ref_id, pos)
    cmd = ['bam-readcount', bam, '-f', fasta,  '-w', '0', '-b', str(mindepth), pos_str]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    lines = io.TextIOWrapper(proc.stdout, encoding="utf-8")
    rows = list(map(parse_line, lines))
    assert len(rows) == 1, str(cmd) + "\nShould be length 1:  " + str(rows)
    return rows[0]

def bam_readcount_all(bam: str, fasta: str, mindepth: int) -> Iterator[BRCRow]:
    cmd = ['bam-readcount', bam, '-f', fasta,  '-w', '0', '-b', str(mindepth)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    lines = io.TextIOWrapper(proc.stdout, encoding="utf-8")
    return map(parse_line, lines)
    
def parse_entry(s: str) -> BRCEntry:
  fs = s.strip().split(':')
  base = fs[0]
  def coerce_(st: str, type_: Callable[[str],Union[int,float,Fraction]]) -> Union[int,float,Fraction]:
      if type_ is Fraction: # type: ignore
          return float(st)
      else:
          return type_(st)
  fields = map(coerce_, fs[1:], [x[1] for x in list(BRCEntry._field_types.items())[1:]])
  result =  BRCEntry(base, *fields)  # type: ignore
  return cast(BRCEntry, result)

def parse_line(s: str) -> BRCRow:
  if PRINT_BAMRC_TOOL:
    print(s)
  _fs = s.split('\t')
  # TODO: What's up with below???
  fs = list(filter(lambda x: not x.startswith(('-', '+')), _fs))
  assert len(fs) == 4 + 6, fs
  entries = {e.base : e for e in map(parse_entry, fs[4:]) }
  nuc = cast(Nuc, fs[2])
  return BRCRow(fs[0], int(fs[1]), nuc, int(fs[3]), cast(BRCEntries, entries))

#TODO:
'''
assert len(fs) == 4 + 6, fs
AssertionError: ['Den1/U88535_1/WestPac/1997/Den1_1', '484', 'G', '5610', '=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00', 'A:2756:59.99:36.73:0.13:1459:1297:0.48:0.01:48.66:1459:0.48:242.22:0.46', 'C:2:60.00:13.00:0.00:2:0:0.14:0.01:27.00:2:0.06:249.50:0.06', 'G:2849:60.00:36.90:0.06:1522:1327:0.49:0.00:12.76:1522:0.47:242.13:0.47', 'T:1:60.00:16.00:0.00:0:1:0.42:0.00:16.00:0:0.00:250.00:0.79', 'N:1:60.00:15.00:60.00:1:0:0.29:0.14:403.00:1:0.59:504.00:0.41', '-G:1:60.00:0.00:0.00:0:1:0.02:0.01:33.00:0:0.00:248.00:0.01\n']
'''
