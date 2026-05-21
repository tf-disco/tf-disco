"""
Utility functions and classes for pattern parsing, matching, and analysis.

Function of interest is `compute_pattern_scores()`, which computes vagueness
score, expected matches, and Z-score for each pattern.

Script authors:
- Joseph Cijo - https://github.com/joejo-joestar
- Sreenikethan Iyer - http://github.com/SreenikethanI
"""

import math
from dataclasses import dataclass
from collections import Counter
from typing import Callable, Literal, cast

import pandas as pd
import numpy as np

from .data_loading import TFInfo



#region Constants
AA_COUNT = 20
AminoAcid = Literal["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
#endregion

#region AST Parser
@dataclass
class Node:
 def __new__(cls, *args, **kwargs):
  # disable creating instances of the base class
  if cls is Node: raise TypeError("NodeBaseClass cannot be instantiated directly")
  return super().__new__(cls)
@dataclass
class CharNode(Node):
    """A literal amino acid character, e.g. `A`, `R`, `N`, etc."""
    char: AminoAcid
    """Single uppercase letter representing the amino acid."""
@dataclass
class WildcardNode(Node):
    """A wildcard matching any amino acid, represented by `.` in regex."""
    pass
@dataclass
class AnchorNode(Node):
    """An anchor matching the start or end of the sequence, represented by `^` or `$` in regex."""
    kind: Literal["^", "$"]
    """kind" is either `^` or `$`."""
@dataclass
class CharClassNode(Node):
    """A character class, e.g. `[ACD]` or `[^ACD]`."""
    chars: frozenset[AminoAcid]
    negated: bool
@dataclass
class RepeatNode(Node):
    node: Node
    lo: int
    hi: int
@dataclass
class GroupNode(Node):
    inner: Node
@dataclass
class SequenceNode(Node):
    nodes: list[Node]
@dataclass
class BranchNode(Node):
    alts: list[SequenceNode]

class PatternParser:
    s: str
    """Pattern string being parsed."""
    pos: int
    """Current position in the pattern string."""

    def __init__(self, s: str):
        """Parse ELM regex patterns into an AST for probability computation.
        Supported syntax:
        - Literals (Amino Acids): `A`, `R`, `N`, `D`, `C`, `Q`, `E`, `G`, `H`, `I`, `L`, `K`, `M`, `F`, `P`, `S`, `T`, `W`, `Y`, `V`
        - Wildcard: `.`
        - Anchors: `^`, `$`
        - Character classes: `[ACD]`, `[^ACD]`
        - Quantifiers: `{n}`, `{n,}`, `{n,m}`
        - Grouping: `(pattern)`
        - Alternation: `pattern|pattern`
        """
        self.s = s
        self.pos = 0

    def peek(self) -> str | None:
        """Peek next character, or `None` if at end of string."""
        return self.s[self.pos] if self.pos < len(self.s) else None

    def consume(self) -> str:
        """Consume and return the next character. Will raise an error if at the
        end of the string."""
        c = self.s[self.pos]
        self.pos += 1
        return c

    def parse(self) -> Node:
        """Parse the pattern and return the root of the AST."""
        node = self._branch()
        if self.pos != len(self.s):
            raise ValueError(f"Unexpected char {self.peek()!r} at pos {self.pos}")
        return node

    # ------------------------------------------------------------------------ #
    # Helper functions

    def _digits(self) -> str:
        """Parse a sequence of digits for quantifiers."""
        s = ""
        while (self.peek() or "").isdigit():
            s += self.consume()
        return s

    def _quantifier(self, unbounded_hi_cap: int=5) -> tuple[int, int] | None:
        """Parse a quantifier, e.g. `{n}`, `{n,}`, or `{n,m}`. Returns (lo, hi)
        or None if no quantifier.

        :param unbounded_hi_cap:
            For unbounded quantifiers like `{n,}`, we set `hi` to `n +
            unbounded_hi_cap` as a heuristic, since ELM patterns are generally
            short and we want to avoid infinite loops in vagueness calculation.
            The default `unbounded_hi_cap` is set to 5.
        """
        if self.peek() != "{": return None
        self.consume() # {
        lo = int(self._digits())
        hi = lo
        if self.peek() == ",":
            self.consume()
            hi_s = self._digits()
            hi = int(hi_s) if hi_s else lo + unbounded_hi_cap
        self.consume() # }
        return lo, hi

    def _char_class(self) -> CharClassNode:
        """Parse a character class, e.g. `[ACD]` or `[^ACD]`."""
        self.consume() # [
        negated = self.peek() == "^"
        if negated: self.consume()
        chars: set[AminoAcid] = set()
        while self.peek() != "]":
            chars.add(cast(AminoAcid, self.consume()))
        self.consume() # ]
        return CharClassNode(frozenset(chars), negated)

    # ------------------------------------------------------------------------ #
    # Main functions

    def _atom(self) -> WildcardNode | CharNode | CharClassNode | AnchorNode | GroupNode:
        """Parse a single atom (literal, wildcard, char class, group, anchor)."""
        char = self.peek()
        if char == ".":
            self.consume(); return WildcardNode()
        if char == "[":
            return self._char_class()
        if char == "(":
            self.consume()
            inner = self._branch()
            if self.peek() == ")": self.consume()
            return GroupNode(inner)
        if char in ("^", "$"):
            self.consume(); return AnchorNode(char)
        if char and char.isupper():
            self.consume(); return CharNode(cast(AminoAcid, char))
        raise ValueError(f"Unexpected char {char!r} at pos {self.pos}")

    def _sequence(self) -> SequenceNode:
        """Parse a sequence of atoms."""
        nodes: list[Node] = []
        while self.peek() not in (None, ")", "|"):
            atom = self._atom()
            quant = self._quantifier()
            nodes.append(RepeatNode(atom, *quant) if quant else atom)
        return SequenceNode(nodes)

    def _branch(self) -> BranchNode|SequenceNode:
        """Parse alternation branches."""
        alts = [self._sequence()]
        while self.peek() == "|":
            self.consume()
            alts.append(self._sequence())
        return BranchNode(alts) if len(alts) > 1 else alts[0]

#endregion

#region Vagueness score computation
VaguenessPenaltyFunction = Callable[[float], float]
"""A function used to penalize wildcards and long character classes. It accepts
the number/count of possible residues in the given Node (e.g. length of
character class) and returns a penalty multiplier for the node's vagueness
score."""

def make_vagueness_penalty_func(max_penalty_at_20: float=1.0) -> VaguenessPenaltyFunction:
    """Returns a function which applies a penalty based on the number of
    possible residues.

    :param float max_penalty_at_20: Control the penalty applied to a node of
        (which has 20 possible residues). If `max_penalty_at_20` is 1.0, there
        is no penalty regardless of the number of residues.
    """

    if max_penalty_at_20 == 1.0:
        return lambda count: 1.0
    else:
        return lambda count: math.pow(count, 1/math.log(20, max_penalty_at_20))

def _calc_vagueness_recursive(node: Node, res_probs: dict[AminoAcid, float], vagueness_penalty_func: VaguenessPenaltyFunction) -> float:
    """
    Calculate the vagueness score for a given node.

    Logic for each node type:
    - `Literal`: Vagueness = probability of that residue.
    - `Wildcard`: Vagueness = vagueness penalty function called with `20`
      (number of amino acids).
    - `Anchor`: Vagueness = `1.0` (doesn't affect probability).
    - `CharClass`: Vagueness = sum of probabilities of the included residues
      (or 1 minus that, if negated), multiplied by the vagueness penalty
      function applied to the number of residues in the class.
    - `Group`: Vagueness = the vagueness of the inner node.
    - `Sequence`: Vagueness = the product of the vagueness of the nodes in the
      sequence.
    - `Branch`: Vagueness = the sum of the vagueness of the alternative
      branches (assuming they are mutually exclusive, which is a heuristic for
      ELM patterns).
    - `Repeat`: Vagueness = the sum of the vagueness of the repeated node
      raised to the power of k, for k from lo to hi (inclusive), divided by the
      number of terms (hi - lo + 1).

    :param vagueness_penalty_func:
        A function which penalizes wildcards and long character classes. A
        wildcard is thought of as a char class with length=20.
    """

    if isinstance(node, CharNode):
        return res_probs.get(node.char, 1 / AA_COUNT)

    if isinstance(node, WildcardNode):
        # A wildcard can match any of the 20 amino acids, so we apply the
        # vagueness penalty function to 20.
        return vagueness_penalty_func(AA_COUNT)

    if isinstance(node, AnchorNode):
        return 1.0

    if isinstance(node, CharClassNode):
        p = sum(
            _calc_vagueness_recursive(CharNode(ch), res_probs, vagueness_penalty_func)
            for ch in node.chars
        )
        char_count = len(node.chars)
        return (
            (1 - p) * vagueness_penalty_func(AA_COUNT - char_count)
            if node.negated else
            p * vagueness_penalty_func(char_count)
        )

    if isinstance(node, GroupNode):
        return _calc_vagueness_recursive(node.inner, res_probs, vagueness_penalty_func)

    if isinstance(node, SequenceNode):
        p = 1.0
        for n in node.nodes: p *= _calc_vagueness_recursive(n, res_probs, vagueness_penalty_func)
        return p

    # ! Assumes all branches are mutually exclusive for ELM patterns.
    if isinstance(node, BranchNode):
        ps = [_calc_vagueness_recursive(a, res_probs, vagueness_penalty_func) for a in node.alts]
        # return sum(ps) / len(node.alts)
        # return min(1.0, sum(ps))
        return sum(ps)

    if isinstance(node, RepeatNode):
        sub_p = _calc_vagueness_recursive(node.node, res_probs, vagueness_penalty_func)
        n = node.hi - node.lo + 1
        # return sum(sub_p ** k for k in range(node.lo, node.hi + 1)) / n
        # return (sub_p ** node.lo) * min(1.0, sum( sub_p ** k for k in range(n) ))
        return (sub_p ** node.lo) * sum( sub_p ** k for k in range(n) )

    return 1.0

def calc_vagueness(pattern: str, residue_probabilities: dict[AminoAcid, float], vagueness_penalty_func: VaguenessPenaltyFunction) -> float:
    """
    Probability that a random protein sequence position matches this ELM pattern.
    Lower = more specific.

    For more info, see the docstring of `_calc_vagueness_recursive()`.

    Args:
        pattern (str): ELM regex pattern
        residue_probabilities (dict[str, float]): The probability of each amino
            acid residue, computed from the FASTA sequences. If a residue is
            missing, it is assumed to have uniform probability (1/20).
        vagueness_penalty_func: A function which penalizes wildcards and long character
            classes. A wildcard is thought of as a char class with length=20.
    Returns:
        float: Probability of pattern (lower = more specific). `NaN` if invalid pattern.
    """
    try:
        parsed_pattern = PatternParser(pattern).parse()
        return _calc_vagueness_recursive(parsed_pattern, residue_probabilities, vagueness_penalty_func)
    except Exception:
        return float("nan")

#endregion

def calculate_residue_probabilities(matches_df: pd.DataFrame, sequence_dict: dict[str, TFInfo]) -> tuple[dict[AminoAcid, float], int]:
    """
    Calculate the probability of each amino acid residue across all sequences
    in the given FASTA file.

    :returns:
    Returns two things: Residue probabilities, and Total count of all residues
    """

    residue_counts: Counter[AminoAcid] = Counter()
    total_length = 0

    for genus_num in matches_df["Genus_Num"].unique():
        if genus_num not in sequence_dict: continue
        seq = sequence_dict[genus_num].Sequence
        total_length += len(seq)
        residue_counts.update(seq) # type: ignore

    residue_probabilities: dict[AminoAcid, float] = {
        residue: count/total_length
        for residue,count in residue_counts.items()
    }

    return residue_probabilities, total_length

#region Score Patterns
def compute_pattern_scores(
    matches_df: pd.DataFrame,
    sequence_dict: dict[str, TFInfo],
    vagueness_penalty_func: VaguenessPenaltyFunction=make_vagueness_penalty_func(1.0),
) -> pd.DataFrame:
    """Compute vagueness score, expected matches, and Z-score for each pattern.

    Args:
        matches_df (pd.DataFrame): DataFrame with column `ELM_Acc` (observed match rows)
        fasta_path (Path): Path to FASTA file, used to compute total sequence length.
        save_file (Path | None): Optional path to save the results as TSV.
        vagueness_penalty_func: A function which penalizes wildcards and long character classes. A wildcard is thought of as a char class with length=20.
    Returns:
        pd.DataFrame: DataFrame with columns: `ELM_Acc`, `ELM_Id`, `Regex`,
            `Vagueness`, `Expected`, `Observed`, `ZScore`, `Log2FC`
    """

    # get unique patterns
    scores = matches_df[["ELM_Acc", "ELM_Id", "Regex"]].drop_duplicates()

    # calculate vagueness and expected
    residue_probabilities, total_length = calculate_residue_probabilities(matches_df, sequence_dict)
    scores["Vagueness"] = scores["Regex"].apply(calc_vagueness, args=(residue_probabilities, vagueness_penalty_func)).astype(float)
    scores["Expected"]  = scores["Vagueness"] * total_length

    # calculate observed
    observed = matches_df.groupby("ELM_Acc").size().rename("Observed")
    scores = scores.join(observed, how="left", on="ELM_Acc")
    scores["Observed"]  = scores["Observed"].fillna(0).astype(int)

    # calculate scores
    scores["ZScore"] = (scores["Observed"] - scores["Expected"]) / np.sqrt(scores["Expected"])
    scores["Log2FC"] = np.log2((scores["Observed"] + 1) / (scores["Expected"] + 1))
    scores = scores.sort_values("Log2FC", ascending=False).reset_index(drop=True)

    return scores

#endregion
