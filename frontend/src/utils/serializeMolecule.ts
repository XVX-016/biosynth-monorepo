// Very small .mol (V2000-like) serializer for small molecules
export function exportAsMol(mol: { atoms: any[]; bonds: any[]; name?: string }) {
    const atoms = mol.atoms;
    const bonds = mol.bonds;
    const name = mol.name || "Mol";
    const header = `${name}\n  MolForge\n\n`;
    const counts = `${String(atoms.length).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0            999 V2000\n`;

    const atomLines = atoms.map(a =>
        `${a.x.toFixed(4).toString().padStart(10, ' ')}${a.y.toFixed(4).toString().padStart(10, ' ')}${(a.z || 0).toFixed(4).toString().padStart(10, ' ')} ${a.element.padEnd(3, ' ')}  0  0  0  0  0  0  0  0  0  0`
    ).join("\n") + "\n";

    const bondLines = bonds.map((b) => {
        const ai = atoms.findIndex(at => at.id === b.a) + 1;
        const bi = atoms.findIndex(at => at.id === b.b) + 1;
        return `${String(ai).padStart(3, ' ')}${String(bi).padStart(3, ' ')}${String(b.order || 1).padStart(3, ' ')}  0  0  0  0`;
    }).join("\n") + "\nM  END\n";

    return header + counts + atomLines + bondLines;
}
