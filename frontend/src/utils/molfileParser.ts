/**
 * Molfile Parser Utility
 * Parses V2000/V3000 molfiles into atoms and bonds
 */

export interface Atom {
  x: number;
  y: number;
  z: number;
  element: string;
}

export interface Bond {
  a1: number;
  a2: number;
  order?: number;
}

export function parseMolfile(molfile: string): { atoms: Atom[]; bonds: Bond[] } {
  const lines = molfile.split(/\r?\n/);
  
  if (lines.length < 4) {
    return { atoms: [], bonds: [] };
  }

  // Counts line is usually at line index 3 (0-based)
  const countsLine = lines[3] || '';
  const atomCount = parseInt(countsLine.slice(0, 3).trim() || '0', 10) || 0;
  const bondCount = parseInt(countsLine.slice(3, 6).trim() || '0', 10) || 0;

  const atoms: Atom[] = [];
  const bonds: Bond[] = [];

  // Parse atoms (lines 4 to 4+atomCount)
  for (let i = 0; i < atomCount; i++) {
    const line = lines[4 + i];
    if (!line) continue;
    
    const x = parseFloat(line.slice(0, 10).trim());
    const y = parseFloat(line.slice(10, 20).trim());
    const z = parseFloat(line.slice(20, 30).trim());
    const element = line.slice(31, 34).trim() || line.slice(0, 3).trim() || 'C';
    
    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
      atoms.push({ x, y, z, element });
    }
  }

  // Parse bonds (lines 4+atomCount to 4+atomCount+bondCount)
  for (let i = 0; i < bondCount; i++) {
    const line = lines[4 + atomCount + i];
    if (!line) continue;
    
    const a1 = parseInt(line.slice(0, 3).trim(), 10) - 1; // Convert to 0-indexed
    const a2 = parseInt(line.slice(3, 6).trim(), 10) - 1;
    const order = parseInt(line.slice(6, 9).trim(), 10) || 1;
    
    if (!isNaN(a1) && !isNaN(a2) && a1 >= 0 && a2 >= 0 && a1 < atoms.length && a2 < atoms.length) {
      bonds.push({ a1, a2, order });
    }
  }

  return { atoms, bonds };
}

