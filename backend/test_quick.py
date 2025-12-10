# Quick Test Script
from chem.engine.bond_engine import BondEngine
from chem.core.models import Atom

engine = BondEngine()

# Test with closer atoms (typical C-H bond length is ~1.09 Å)
atoms = {
    0: Atom(id=0, element="C"),
    1: Atom(id=1, element="H"),
    2: Atom(id=2, element="H"),
    3: Atom(id=3, element="H"),
    4: Atom(id=4, element="H"),
}

positions = {
    0: (0.0, 0.0, 0.0),
    1: (1.09, 0.0, 0.0),
    2: (-1.09, 0.0, 0.0),
    3: (0.0, 1.09, 0.0),
    4: (0.0, -1.09, 0.0),
}

bonds = engine.predict_all_bonds(atoms, positions)

print(f"✅ Found {len(bonds)} bonds")
for b in bonds:
    print(f"   Bond: {b.atom_a} - {b.atom_b}, order: {b.order}")

assert len(bonds) == 4, f"Expected 4 bonds, got {len(bonds)}"
print("✅ All tests passed!")
