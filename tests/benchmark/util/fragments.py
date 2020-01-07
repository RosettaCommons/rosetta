

import typing

import numpy as np

FRAG_COLUMN_FORMAT = {
    "origin_pdb_code": (1, 5),
    "origin_pdb_chain": (6),
    "origin_pdb_res_num": (8, 13),
    "origin_pdb_aa": (14),
    "origin_pdb_ss": (16),
    "phi": (18, 27),
    "psi": (27, 36),
    "omega": (36, 44),
    "x_coord": (45, 53),
    "y_coord": (54, 62),
    "z_coord": (63, 72)
}


class Fragment:
    """
    Store a single fragment

        1      -- blank
        2-5    -- PDB code for the fragment origin
        7      -- chain ID for the origin PDB
        9-13   -- PDB residue number for the origin PDB
        15     -- amino acid identity in the origin PDB
        17     -- secondary structure for the origin PDB (Helix, Loop, Extended/beta)
        19-27  -- phi
        28-36  -- psi
        37-45  -- omega
        46-54  -- C-alpha x coordinate for origin PDB (optional)
        55-63  -- C-alpha y coordinate for origin PDB (optional)
        65-73  -- C-alpha z coordinate for origin PDB (optional)
        74-79  -- unknown (unused)
        80-85  -- unknown (unused)
        86     -- Literal "P" (unused)
        87-89  -- fragment position number, pose numbered (unused)
        91     -- Literal "F"(unused)
        92-94  -- fragment number (unused)
    """
    def __init__(self, fragment_lines: typing.Optional[typing.List[str]] = None) -> None:
        self.origin_pdb_codes = []
        self.origin_pdb_chains = []
        self.origin_pdb_res_nums = []
        self.origin_pdb_aas = []
        self.origin_pdb_sss = []
        self.phis = []
        self.psis = []
        self.omegas = []
        self.xyzs = []
        if not fragment_lines:
            return
        self.parse_fragment_lines(fragment_lines)

    def parse_fragment_lines(self, fragment_lines: typing.List[str]) -> None:
        """parse a list of fragment file lines(str) into this class
        """
        for line in fragment_lines:
            if not line.strip():
                continue
            self.origin_pdb_codes.append(line[slice(*FRAG_COLUMN_FORMAT["origin_pdb_code"])])
            self.origin_pdb_chains.append(line[FRAG_COLUMN_FORMAT["origin_pdb_chain"]])
            self.origin_pdb_res_nums.append(int(line[slice(*FRAG_COLUMN_FORMAT["origin_pdb_res_num"])]))
            self.origin_pdb_aas.append(line[FRAG_COLUMN_FORMAT["origin_pdb_aa"]])
            self.origin_pdb_sss.append(line[FRAG_COLUMN_FORMAT["origin_pdb_ss"]])
            self.phis.append(float(line[slice(*FRAG_COLUMN_FORMAT["phi"])]))
            self.psis.append(float(line[slice(*FRAG_COLUMN_FORMAT["psi"])]))
            self.omegas.append(float(line[slice(*FRAG_COLUMN_FORMAT["omega"])]))
            if len(line) > 50:
                x = float(line[slice(*FRAG_COLUMN_FORMAT["x_coord"])])
                y = float(line[slice(*FRAG_COLUMN_FORMAT["y_coord"])])
                z = float(line[slice(*FRAG_COLUMN_FORMAT["z_coord"])])
                self.xyzs.append(np.array((x, y, z)))
        self.size = len(self.origin_pdb_codes)

    def gen_fragment_lines(self) -> typing.List[str]:
        """generate a list of fragment file lines(str) from this class and return it
        """
        lines = []
        for i in range(len(self.origin_pdb_codes)):
            lines.append(f" {self.origin_pdb_codes[i]} {self.origin_pdb_chains[i]} {self.origin_pdb_res_nums[i]: 5}"
                         f" {self.origin_pdb_aas[i]} {self.origin_pdb_sss[i]}"
                         f" {self.phis[i]: >8.3f} {self.psis[i]: >8.3f} {self.omegas[i]: >8.3f}"
                         f" {self.xyzs[i][0]: >8.3f} {self.xyzs[i][1]: >8.3f} {self.xyzs[i][2]: >8.3f}")
        return lines

    def gen_fragment_string(self) -> str:
        """generate the string data of a fragment file from this class and return it
        """
        return "\n".join(self.gen_fragment_lines())


class Fragments:
    """
    Store a list of 'Fragment's in this (similar to how rosetta does it)
    """

    def __init__(self):
        self.per_residue_fragments = []

    def splice_fragments(self, lines: typing.List[str]) -> typing.List[Fragment]:
        """group fragment lines by their position, convert them into Fragment objects and return them.
        """
        grouped = [[]]
        for line in lines:
            if not line.strip() and grouped[-1]:
                grouped.append([])
                continue
            if not line.strip():
                continue
            grouped[-1].append(line)
        return [Fragment(group) for group in grouped if group]

    def parse_fragment_file_lines(self, fragment_file_lines: typing.List[str]) -> None:
        """given all lines of a fragment file, parse them into List[List[Fragment]]
        so that we have a list of fragments for every position.
        """
        self.per_residue_fragments = []
        start = 1
        for i, line in enumerate(fragment_file_lines):
            if i == 0:  # Skip because file starts with position line, want to start after that line
                continue
            if line.startswith("position:"):
                end = i-1
                self.per_residue_fragments.append(self.splice_fragments(fragment_file_lines[start:end]))
                start = i+1
        self.per_residue_fragments.append(self.splice_fragments(fragment_file_lines[start:]))

    def parse_fragment_file_string(self, fragment_file_string: str) -> None:
        """parse fragment file string data
        """
        self.parse_fragment_file_lines(fragment_file_string.split('\n'))

    def parse_fragment_file(self, fragment_file_name: str) -> None:
        """parse fragment file from a filename
        """
        with open(fragment_file_name) as fh:
            self.parse_fragment_file_string(fh.read())

    def gen_fragment_file_lines(self) -> typing.List[str]:
        """generate a list of fragment file lines from this class
        """
        lines = []
        for i, fragments in enumerate(self.per_residue_fragments):
            lines.append(f"position:        {i+1: 5} neighbors:        {len(fragments): 5}")
            lines.append("")
            for fragment in fragments:
                lines += fragment.gen_fragment_lines()
                lines.append("")
        else:
            lines.append("")
        return lines

    def gen_fragment_file_str(self) -> str:
        """generate the full string of a fragment file from this class
        """
        return '\n'.join(self.gen_fragment_file_lines())
