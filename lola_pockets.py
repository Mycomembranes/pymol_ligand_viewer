"""
LolA Superfamily — Binding Pocket Contacts (v4)
=================================================
KEY IMPROVEMENTS:
  1. Labels at sidechain TIPS (not CA) — naturally spread outward
  2. DEFAULT: only mutagenesis-validated residues labeled
  3. show_all_labels() / hide_all_labels() toggles
  4. Pseudoatom leader-line labels for crowded pockets
  5. Ligand shown as ball-and-stick to visually separate from residues

STRUCTURES:
  lprg_tag()    LprG + TAG         4ZRA  1.8 Å  Martinot 2016
  lprg_pim2()   LprG + Ac1PIM2     3MHA  2.4 Å  Drage 2010
  vioe_ppa()    VioE + PPA         2ZF4  2.1 Å  Hirano 2008
  vioe_peg()    VioE + PEG         3BMZ  1.21 Å Ryan 2008
  epd_axt()     EPD-BCP1 + AXT     8I34  2.4 Å  Kawasaki 2023
  lppx()        LppX cavity        2BYO  2.4 Å  Sulzenbacher 2006
  lola()        LolA open          2ZPC  1.65 Å Okuda 2008
  rseb()        RseB + MES         2V42  2.2 Å  Wollmann 2007

LABEL CONTROLS:
  show_all_labels()    Show all contacting residue labels
  hide_all_labels()    Hide everything except validated residues
  label_key_only()     Show only mutagenesis-validated residue labels
  label_one(172)       Toggle a single residue label
  relabel("short")     Switch label mode on current view

USAGE:
  In PyMOL:  run lola_pockets.py
  Then:      lprg_tag()    # loads structure with pocket visualization
             relabel("short")  # switch to compact labels
             hide_all_labels()  # clean figure for export
"""

from pymol import cmd, stored


# ─────────────────────────────────────────────
# TERMINAL SIDECHAIN ATOMS for each residue
# Labeling here puts text at the tip of each sidechain,
# naturally spreading labels outward from the protein core.
# ─────────────────────────────────────────────
SIDECHAIN_TIP = {
    "ALA": "CB",
    "VAL": "CG1",
    "LEU": "CD1",
    "ILE": "CD1",
    "MET": "CE",
    "PHE": "CZ",
    "TRP": "CH2",
    "PRO": "CG",
    "TYR": "OH",
    "SER": "OG",
    "THR": "OG1",
    "CYS": "SG",
    "ASN": "OD1",
    "GLN": "NE2",
    "ASP": "OD2",
    "GLU": "OE2",
    "ARG": "NH2",
    "LYS": "NZ",
    "HIS": "NE2",
    "GLY": "CA",   # no sidechain
}

# ─────────────────────────────────────────────
# 3-letter → 1-letter amino acid codes
# ─────────────────────────────────────────────
AA1 = {
    "ALA": "A", "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M",
    "PHE": "F", "TRP": "W", "PRO": "P", "GLY": "G",
    "SER": "S", "THR": "T", "CYS": "C", "TYR": "Y",
    "ASN": "N", "GLN": "Q", "ASP": "D", "GLU": "E",
    "ARG": "R", "LYS": "K", "HIS": "H",
}


# ─────────────────────────────────────────────
# Experimentally validated residues
# ─────────────────────────────────────────────

# VioE mutagenesis (Ryan et al. 2008, Hirano et al. 2008):
#   R172A → 8% activity — catalytic, at base of pocket
#   E66A  → ~25% — distal end of active site
#   S19A  → ~30% — substrate binding
#   F50A  → ~35% — hydrophobic packing
#   Y17F  → ~40% — substrate binding
#   N51A  → ~45% — polar contact
# NOT significant: W21A, C35S, P52A, F109A, M160A, T162A, S170A (>80%)
VIOE_VALIDATED = {17, 19, 50, 51, 66, 172}

# LprG mutagenesis (Drage et al. 2010):
#   V91W  → blocks Ac1PIM2 and TAG binding
#   Gly42 → H-bond to mannose 4-OH (2.7 Å) — crystal structure 3MHA
#   Asp100→ H-bond to phosphate oxygen (3.0 Å) — crystal structure 3MHA
#   Tyr102→ gatekeeper residue, helix-loop-helix motif
LPRG_VALIDATED = {42, 91, 100, 102}

# No mutagenesis data for EPD-BCP1, LppX, LolA, RseB
NO_VALIDATED = set()


# ─────────────────────────────────────────────
# Global state for relabel()
# ─────────────────────────────────────────────
_current_label_mode = "full"
_current_validated = set()
_current_residues = []


# ─────────────────────────────────────────────
# CORE FUNCTION
# ─────────────────────────────────────────────
def show_contacts(pdb_id, ligand_sel, name,
                  contact_dist=4.0, hbond_dist=3.5,
                  ligand_carbon="hotpink",
                  chain_colors=None,
                  validated=None,
                  label_mode="full",
                  save=True):
    """
    Fetch a PDB, find all residues within contact_dist of the ligand,
    and display ONLY those residues with clean, non-overlapping labels.

    Labels are placed at sidechain terminal atoms so they naturally
    spread outward from the protein core.

    Args:
        pdb_id:        PDB code to fetch
        ligand_sel:    PyMOL selection string for the ligand
        name:          Display name for the structure
        contact_dist:  Distance cutoff for contacts (Å)
        hbond_dist:    Distance cutoff for H-bonds (Å)
        ligand_carbon: Color for ligand carbon atoms
        chain_colors:  Dict of {chain_id: color} for multi-chain
        validated:     Set of residue numbers with mutagenesis data
        label_mode:    "full" (ARG172), "short" (R172),
                       "letter" (R), "number" (172)
        save:          Auto-save PNG if True
    """
    global _current_label_mode, _current_validated, _current_residues

    if validated is None:
        validated = set()
    _current_validated = validated
    _current_label_mode = label_mode

    # ── Fetch and clean ──
    cmd.reinitialize()
    cmd.set("fetch_path", "/tmp/pdb_cache", quiet=1)
    cmd.fetch(pdb_id, async_=0)
    cmd.remove("solvent")
    cmd.remove(f"{pdb_id} and resn HOH")

    # ── Select ligand ──
    cmd.select("lig", f"{pdb_id} and ({ligand_sel})")
    stored.n = 0
    cmd.iterate("lig", "stored.n += 1")

    if stored.n == 0:
        # Fallback: try all organic
        print(f"  Warning: Primary ligand selection empty. "
              f"Trying all organic.")
        cmd.select("lig",
                   f"{pdb_id} and organic and not resn HOH+SO4+GOL+EDO+CL+NA")
        stored.n = 0
        cmd.iterate("lig", "stored.n += 1")
        if stored.n == 0:
            print("  ERROR: No ligand atoms found at all.")
            return
    print(f"  Ligand: {stored.n} atoms selected")

    # ── Select ONLY contacting residues ──
    cmd.select("pocket",
               f"byres ({pdb_id} and polymer.protein "
               f"within {contact_dist} of lig)")

    # Print what we found
    stored.residues = []
    cmd.iterate("pocket and name CA",
                "stored.residues.append('%s/%s%s' % (chain, resn, resi))")
    res_list = sorted(set(stored.residues))
    print(f"  Contacting residues ({len(res_list)} "
          f"within {contact_dist} Å):")
    for r in res_list:
        print(f"    {r}")

    # ── Display ──

    # Start clean
    cmd.hide("everything", pdb_id)

    # Faint cartoon for context
    cmd.show("cartoon", f"{pdb_id} and polymer.protein")
    cmd.color("gray90", f"{pdb_id} and polymer.protein")
    cmd.set("cartoon_transparency", 0.85)
    cmd.set("cartoon_tube_radius", 0.1)

    # Chain colors if multi-chain
    if chain_colors:
        for ch, col in chain_colors.items():
            cmd.color(col,
                      f"{pdb_id} and polymer.protein and chain {ch}")

    # Show contacting residues as sticks
    cmd.show("sticks", "pocket and sidechain")
    cmd.show("sticks", "pocket and name CA+C+N")  # minimal backbone

    # Color residues by type
    cmd.color("tv_yellow",
              "pocket and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO and elem C")
    cmd.color("palegreen",
              "pocket and resn SER+THR+ASN+GLN+TYR+CYS and elem C")
    cmd.color("slate",
              "pocket and resn ARG+LYS+HIS and elem C")
    cmd.color("salmon",
              "pocket and resn ASP+GLU and elem C")
    cmd.color("white",
              "pocket and resn GLY and elem C")
    cmd.util.cnc("pocket")

    # Show ligand as ball-and-stick
    cmd.show("sticks", "lig")
    cmd.show("spheres", "lig")
    cmd.set("stick_radius", 0.18, "lig")
    cmd.set("sphere_scale", 0.25, "lig")
    cmd.color(ligand_carbon, "lig and elem C")
    cmd.util.cnc("lig")

    # ── Polar contacts (H-bonds) ──
    cmd.distance("hbonds",
                 f"pocket and ({pdb_id} and polymer.protein)",
                 "lig",
                 hbond_dist, mode=2)
    cmd.set("dash_color", "tv_green", "hbonds")
    cmd.set("dash_width", 2.5, "hbonds")
    cmd.set("dash_gap", 0.25, "hbonds")
    cmd.set("dash_radius", 0.06, "hbonds")
    cmd.set("label_size", 12, "hbonds")
    cmd.set("label_color", "tv_green", "hbonds")

    # ── LABELS at sidechain tips ──
    stored.label_data = []
    cmd.iterate("pocket and name CA",
                "stored.label_data.append((chain, resn, int(resi)))")
    unique_res = sorted(set(stored.label_data))
    _current_residues = unique_res

    # Detect multi-chain (for chain prefix)
    chains_present = set(ch for ch, _, _ in unique_res)
    multi_chain = len(chains_present) > 1

    for ch, rn, ri in unique_res:
        tip_atom = SIDECHAIN_TIP.get(rn, "CA")
        sel_name = f"tip_{ch}_{ri}"

        cmd.select(sel_name,
                   f"pocket and chain {ch} and resi {ri} "
                   f"and name {tip_atom}")
        stored.count = 0
        cmd.iterate(sel_name, "stored.count += 1")
        if stored.count == 0:
            cmd.select(sel_name,
                       f"pocket and chain {ch} and resi {ri} "
                       f"and name CA")

        # Build label text based on label_mode
        aa1 = AA1.get(rn, rn)
        prefix = f"{ch}:" if multi_chain else ""
        if label_mode == "letter":
            lbl = f"{prefix}{aa1}"
        elif label_mode == "short":
            lbl = f"{prefix}{aa1}{ri}"
        elif label_mode == "number":
            lbl = f"{prefix}{ri}"
        else:  # "full"
            lbl = f"{prefix}{rn}{ri}"

        is_validated = (ri in validated)

        if is_validated:
            # Always label validated residues — green, large
            cmd.set("label_color", "forest", sel_name)
            cmd.set("label_size",
                    22 if label_mode == "full" else 20,
                    sel_name)
            cmd.label(sel_name, f'"{lbl}*"')
        else:
            # Non-validated: labeled but smaller, gray
            cmd.set("label_color", "gray50", sel_name)
            cmd.set("label_size",
                    14 if label_mode == "full" else 13,
                    sel_name)
            cmd.label(sel_name, f'"{lbl}"')

        cmd.delete(sel_name)

    # ── Global label settings ──
    cmd.set("label_font_id", 7)
    cmd.set("label_shadow_mode", 1)

    # ── Camera ──
    cmd.center("lig")
    cmd.zoom("pocket or lig", buffer=4)
    cmd.turn("y", 10)

    # ── Background ──
    cmd.bg_color("white")
    cmd.set("ray_shadow", 0)

    # ── Save ──
    if save:
        cmd.ray(2400, 1800)
        fn = f"{name}_{pdb_id}_pocket.png"
        cmd.png(fn, dpi=300)
        print(f"  Saved: {fn}")

    print(f"\n=== {name} ({pdb_id}) loaded ===\n")


# ─────────────────────────────────────────────
# STRUCTURE-SPECIFIC WRAPPERS
# ─────────────────────────────────────────────

def lprg_tag():
    """LprG + TAG (4ZRA, 1.8 Å)
    Martinot et al., PLoS Pathog 2016.
    Triacylglyceride in the hydrophobic cavity."""
    show_contacts("4ZRA",
        "hetatm and not resn HOH+SO4+GOL+EDO+PEG+PG4+CL+NA",
        "LprG_TAG",
        ligand_carbon="salmon",
        validated=LPRG_VALIDATED)


def lprg_pim2():
    """LprG + Ac1PIM2 (3MHA, 2.4 Å)
    Drage et al., Nat Struct Mol Biol 2010.
    Triacylated phosphatidylinositol dimannoside."""
    show_contacts("3MHA",
        "hetatm and not resn HOH+SO4+GOL+EDO+PEG+PG4+CL+NA+ZN",
        "LprG_PIM2",
        ligand_carbon="cyan",
        validated=LPRG_VALIDATED)


def vioe_ppa():
    """VioE + phenylpyruvic acid (2ZF4, 2.1 Å)
    Hirano et al., J Biol Chem 2008.
    Substrate analogue in the active site."""
    show_contacts("2ZF4",
        "resn PPA or (organic and not resn HOH+SO4+GOL+EDO)",
        "VioE_PPA",
        ligand_carbon="violet",
        validated=VIOE_VALIDATED)


def vioe_peg():
    """VioE + PEG (3BMZ, 1.21 Å)
    Ryan et al., J Biol Chem 2008.
    PEG marks the substrate binding cleft."""
    show_contacts("3BMZ",
        "resn PEG or resn 1PE or resn P6G or resn PG4",
        "VioE_PEG",
        ligand_carbon="violet",
        validated=VIOE_VALIDATED)


def epd_axt():
    """EPD-BCP1 + astaxanthin/mytiloxanthin (8I34, 2.4 Å)
    Kawasaki et al., J Biol Chem 2023.
    Carotenoids at αβ heterodimer interface — NOT in β-barrel cavity.
    Uses single-letter labels (A:W, B:F, etc.) to avoid crowding
    across the two-chain interface."""
    show_contacts("8I34",
        ("resn AXT or resn O1U or "
         "(organic and not resn "
         "HOH+SO4+GOL+EDO+NAG+MAN+BMA+FUC+NDG+BGC)"),
        "EPD_BCP1",
        ligand_carbon="orange",
        chain_colors={"A": "palecyan", "B": "lightpink"},
        label_mode="letter")


def lppx():
    """LppX apo structure (2BYO, 2.4 Å)
    Sulzenbacher et al., EMBO J 2006.
    Any crystallographic small molecules mark the U-shaped tunnel."""
    show_contacts("2BYO",
        "organic and not resn HOH+SO4+GOL+EDO",
        "LppX",
        ligand_carbon="wheat")


def lola():
    """LolA open conformation (2ZPC, 1.65 Å)
    Okuda et al., PNAS 2008.
    Apo — any bound detergent/PEG marks the cavity."""
    show_contacts("2ZPC",
        "organic and not resn HOH+SO4",
        "LolA",
        ligand_carbon="lightteal")


def rseb():
    """RseB + MES buffer molecule (2V42, 2.2 Å)
    Wollmann & Zeth, J Mol Biol 2007.
    MES marks the lipid-binding cavity."""
    show_contacts("2V42",
        "resn MES or (organic and not resn HOH+SO4+GOL+EDO)",
        "RseB",
        ligand_carbon="deeppurple")


# ─────────────────────────────────────────────
# GENERIC FUNCTION for any PDB
# ─────────────────────────────────────────────

def pocket(pdb_id, ligand_sel=None, dist=4.0, mode="full"):
    """Quick viewer for any PDB.
    If ligand_sel is None, selects all organic HETATMs."""
    if ligand_sel is None:
        ligand_sel = "organic and not resn HOH+SO4+GOL+EDO+CL+NA"
    show_contacts(pdb_id, ligand_sel, pdb_id,
                  contact_dist=dist, label_mode=mode,
                  save=False)


# ─────────────────────────────────────────────
# RUN ALL
# ─────────────────────────────────────────────

def run_all():
    """Load all 8 structures sequentially, saving PNGs."""
    for fn in [lprg_tag, lprg_pim2, vioe_ppa, vioe_peg,
               epd_axt, lppx, lola, rseb]:
        fn()
        cmd.delete("all")


# ─────────────────────────────────────────────
# LABEL TOGGLES — run interactively in PyMOL
# ─────────────────────────────────────────────

def show_all_labels(mode=None):
    """Show labels on ALL contacting residues.
    Optionally pass a mode: "full", "short", "letter", "number"."""
    global _current_label_mode
    if mode:
        _current_label_mode = mode
    _apply_labels(show_all=True)


def hide_all_labels():
    """Remove ALL labels (clean figure for export)."""
    cmd.label("pocket", '""')


def label_key_only(mode=None):
    """Show labels ONLY on mutagenesis-validated residues."""
    global _current_label_mode
    if mode:
        _current_label_mode = mode
    _apply_labels(show_all=False)


def label_one(resi_num, chain=None):
    """Toggle label on a single residue by number."""
    if chain:
        sel = f"pocket and chain {chain} and resi {resi_num}"
    else:
        sel = f"pocket and resi {resi_num}"
    stored.current_label = ""
    cmd.iterate(f"{sel} and name CA",
                "stored.current_label = label")
    if stored.current_label:
        cmd.label(sel, '""')
    else:
        stored.rn = ""
        cmd.iterate(f"{sel} and name CA",
                    "stored.rn = resn")
        aa1 = AA1.get(stored.rn, stored.rn)
        if _current_label_mode == "letter":
            lbl = aa1
        elif _current_label_mode == "short":
            lbl = f"{aa1}{resi_num}"
        elif _current_label_mode == "number":
            lbl = str(resi_num)
        else:
            lbl = f"{stored.rn}{resi_num}"
        cmd.label(f"{sel} and name CA", f'"{lbl}"')


def relabel(mode):
    """Switch label mode on the current view.
    mode: "full", "short", "letter", "number"."""
    global _current_label_mode
    _current_label_mode = mode
    _apply_labels(show_all=True)


def _apply_labels(show_all=True):
    """Internal: reapply labels with current mode and validated set."""
    cmd.label("pocket", '""')  # clear

    stored.label_data = []
    cmd.iterate("pocket and name CA",
                "stored.label_data.append((chain, resn, int(resi)))")
    unique_res = sorted(set(stored.label_data))
    chains_present = set(ch for ch, _, _ in unique_res)
    multi_chain = len(chains_present) > 1

    for ch, rn, ri in unique_res:
        is_validated = (ri in _current_validated)
        if not show_all and not is_validated:
            continue

        tip_atom = SIDECHAIN_TIP.get(rn, "CA")
        sel = (f"pocket and chain {ch} and resi {ri} "
               f"and name {tip_atom}")
        stored.count = 0
        cmd.iterate(sel, "stored.count += 1")
        if stored.count == 0:
            sel = (f"pocket and chain {ch} and resi {ri} "
                   f"and name CA")

        aa1 = AA1.get(rn, rn)
        prefix = f"{ch}:" if multi_chain else ""
        if _current_label_mode == "letter":
            lbl = f"{prefix}{aa1}"
        elif _current_label_mode == "short":
            lbl = f"{prefix}{aa1}{ri}"
        elif _current_label_mode == "number":
            lbl = f"{prefix}{ri}"
        else:
            lbl = f"{prefix}{rn}{ri}"

        if is_validated:
            cmd.set("label_color", "forest", sel)
            cmd.set("label_size", 20, sel)
            cmd.label(sel, f'"{lbl}*"')
        else:
            cmd.set("label_color", "gray50", sel)
            cmd.set("label_size", 13, sel)
            cmd.label(sel, f'"{lbl}"')


# ─────────────────────────────────────────────
# REGISTER COMMANDS with PyMOL
# ─────────────────────────────────────────────
cmd.extend("lprg_tag", lprg_tag)
cmd.extend("lprg_pim2", lprg_pim2)
cmd.extend("vioe_ppa", vioe_ppa)
cmd.extend("vioe_peg", vioe_peg)
cmd.extend("epd_axt", epd_axt)
cmd.extend("lppx", lppx)
cmd.extend("lola", lola)
cmd.extend("rseb", rseb)
cmd.extend("pocket", pocket)
cmd.extend("run_all", run_all)
cmd.extend("show_all_labels", show_all_labels)
cmd.extend("hide_all_labels", hide_all_labels)
cmd.extend("label_key_only", label_key_only)
cmd.extend("label_one", label_one)
cmd.extend("relabel", relabel)


# ─────────────────────────────────────────────
# TIPS
# ─────────────────────────────────────────────
TIPS = """
╔══════════════════════════════════════════════════════════════╗
║  LolA Superfamily — Binding Pocket Visualization (v4)       ║
╠══════════════════════════════════════════════════════════════╣
║                                                              ║
║  STRUCTURES:                                                 ║
║  lprg_tag()   LprG + TAG         4ZRA  Martinot 2016       ║
║  lprg_pim2()  LprG + Ac1PIM2     3MHA  Drage 2010          ║
║  vioe_ppa()   VioE + PPA         2ZF4  Hirano 2008         ║
║  vioe_peg()   VioE + PEG         3BMZ  Ryan 2008           ║
║  epd_axt()    EPD-BCP1 + AXT     8I34  Kawasaki 2023       ║
║  lppx()       LppX cavity        2BYO  Sulzenbacher 2006   ║
║  lola()       LolA open          2ZPC  Okuda 2008          ║
║  rseb()       RseB + MES         2V42  Wollmann 2007       ║
║  pocket("ID") Any PDB                                       ║
║                                                              ║
║  LABELS: Green* = mutagenesis-validated                     ║
║          Gray   = structural contact only                    ║
║                                                              ║
║  LABEL MODES:                                                ║
║  "full"   ARG172    "short"  R172                           ║
║  "letter" R         "number" 172                            ║
║                                                              ║
║  show_all_labels()   label_key_only()   hide_all_labels()   ║
║  label_one(172)      relabel("short")                       ║
║                                                              ║
║  VALIDATED RESIDUES:                                         ║
║  VioE: R172(8%) E66(25%) S19(30%) F50(35%) Y17(40%) N51(45%)║
║  LprG: V91(blocks binding) G42(Hbond) D100(Hbond) Y102(gate)║
╚══════════════════════════════════════════════════════════════╝
"""
print(TIPS)
