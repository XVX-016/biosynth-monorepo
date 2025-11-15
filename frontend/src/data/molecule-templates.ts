import { MoleculeGraph } from '@biosynth/engine';

export interface MoleculeTemplate {
  name: string;
  formula: string;
  category: 'organic' | 'inorganic' | 'biomolecule' | 'aromatic';
  atoms: Array<{
    element: 'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I';
    position: [number, number, number];
  }>;
  bonds: Array<{
    from: number; // atom index
    to: number; // atom index
    order: number;
  }>;
}

export const moleculeTemplates: MoleculeTemplate[] = [
  // Water (H2O)
  {
    name: 'Water',
    formula: 'H2O',
    category: 'inorganic',
    atoms: [
      { element: 'O', position: [0, 0, 0] },
      { element: 'H', position: [0.96, 0.26, 0] },
      { element: 'H', position: [-0.24, 0.96, 0] },
    ],
    bonds: [
      { from: 0, to: 1, order: 1 },
      { from: 0, to: 2, order: 1 },
    ],
  },
  
  // Methane (CH4)
  {
    name: 'Methane',
    formula: 'CH4',
    category: 'organic',
    atoms: [
      { element: 'C', position: [0, 0, 0] },
      { element: 'H', position: [1.09, 0, 0] },
      { element: 'H', position: [-0.36, 1.03, 0] },
      { element: 'H', position: [-0.36, -0.52, 0.89] },
      { element: 'H', position: [-0.36, -0.52, -0.89] },
    ],
    bonds: [
      { from: 0, to: 1, order: 1 },
      { from: 0, to: 2, order: 1 },
      { from: 0, to: 3, order: 1 },
      { from: 0, to: 4, order: 1 },
    ],
  },
  
  // Ethanol (C2H5OH)
  {
    name: 'Ethanol',
    formula: 'C2H5OH',
    category: 'organic',
    atoms: [
      { element: 'C', position: [0, 0, 0] },
      { element: 'C', position: [1.54, 0, 0] },
      { element: 'O', position: [2.68, 0, 0] },
      { element: 'H', position: [-0.51, 1.02, 0] },
      { element: 'H', position: [-0.51, -0.51, 0.89] },
      { element: 'H', position: [-0.51, -0.51, -0.89] },
      { element: 'H', position: [1.54, 1.02, 0] },
      { element: 'H', position: [1.54, -0.51, 0.89] },
      { element: 'H', position: [3.22, 0.96, 0] },
    ],
    bonds: [
      { from: 0, to: 1, order: 1 },
      { from: 1, to: 2, order: 1 },
      { from: 0, to: 3, order: 1 },
      { from: 0, to: 4, order: 1 },
      { from: 0, to: 5, order: 1 },
      { from: 1, to: 6, order: 1 },
      { from: 1, to: 7, order: 1 },
      { from: 2, to: 8, order: 1 },
    ],
  },
  
  // Ammonia (NH3)
  {
    name: 'Ammonia',
    formula: 'NH3',
    category: 'inorganic',
    atoms: [
      { element: 'N', position: [0, 0, 0] },
      { element: 'H', position: [1.01, 0, 0] },
      { element: 'H', position: [-0.34, 0.94, 0] },
      { element: 'H', position: [-0.34, -0.47, 0.82] },
    ],
    bonds: [
      { from: 0, to: 1, order: 1 },
      { from: 0, to: 2, order: 1 },
      { from: 0, to: 3, order: 1 },
    ],
  },
  
  // Benzene (C6H6) - simplified planar structure
  {
    name: 'Benzene',
    formula: 'C6H6',
    category: 'aromatic',
    atoms: [
      { element: 'C', position: [0, 1.39, 0] },
      { element: 'C', position: [1.20, 0.70, 0] },
      { element: 'C', position: [1.20, -0.70, 0] },
      { element: 'C', position: [0, -1.39, 0] },
      { element: 'C', position: [-1.20, -0.70, 0] },
      { element: 'C', position: [-1.20, 0.70, 0] },
      { element: 'H', position: [0, 2.48, 0] },
      { element: 'H', position: [2.14, 1.25, 0] },
      { element: 'H', position: [2.14, -1.25, 0] },
      { element: 'H', position: [0, -2.48, 0] },
      { element: 'H', position: [-2.14, -1.25, 0] },
      { element: 'H', position: [-2.14, 1.25, 0] },
    ],
    bonds: [
      { from: 0, to: 1, order: 2 }, // alternating double bonds
      { from: 1, to: 2, order: 1 },
      { from: 2, to: 3, order: 2 },
      { from: 3, to: 4, order: 1 },
      { from: 4, to: 5, order: 2 },
      { from: 5, to: 0, order: 1 },
      { from: 0, to: 6, order: 1 },
      { from: 1, to: 7, order: 1 },
      { from: 2, to: 8, order: 1 },
      { from: 3, to: 9, order: 1 },
      { from: 4, to: 10, order: 1 },
      { from: 5, to: 11, order: 1 },
    ],
  },
  
  // Carbon Dioxide (CO2)
  {
    name: 'Carbon Dioxide',
    formula: 'CO2',
    category: 'inorganic',
    atoms: [
      { element: 'C', position: [0, 0, 0] },
      { element: 'O', position: [1.16, 0, 0] },
      { element: 'O', position: [-1.16, 0, 0] },
    ],
    bonds: [
      { from: 0, to: 1, order: 2 },
      { from: 0, to: 2, order: 2 },
    ],
  },
  
  // Glucose (C6H12O6) - reduced/simplified version
  {
    name: 'Glucose',
    formula: 'C6H12O6',
    category: 'biomolecule',
    atoms: [
      { element: 'C', position: [0, 0, 0] }, // C1
      { element: 'C', position: [1.52, 0, 0] }, // C2
      { element: 'C', position: [2.28, 1.39, 0] }, // C3
      { element: 'C', position: [1.52, 2.78, 0] }, // C4
      { element: 'C', position: [0, 2.78, 0] }, // C5
      { element: 'C', position: [-0.76, 1.39, 0] }, // C6
      { element: 'O', position: [-0.76, -1.39, 0] }, // O1 (ring oxygen)
      { element: 'O', position: [0, 4.17, 0] }, // O2
      { element: 'H', position: [0.51, -0.89, 0.89] }, // H1
      { element: 'H', position: [0.51, -0.89, -0.89] }, // H2
      { element: 'H', position: [1.52, -0.89, 0] }, // H3
      { element: 'H', position: [2.79, 1.39, 0.89] }, // H4
      { element: 'H', position: [2.79, 1.39, -0.89] }, // H5
      { element: 'H', position: [1.52, 3.67, 0.89] }, // H6
      { element: 'H', position: [1.52, 3.67, -0.89] }, // H7
      { element: 'H', position: [-0.51, 3.67, 0] }, // H8
      { element: 'H', position: [-1.27, 1.39, 0.89] }, // H9
      { element: 'H', position: [-1.27, 1.39, -0.89] }, // H10
      { element: 'H', position: [-0.76, -1.39, 0.89] }, // H11
      { element: 'H', position: [-0.76, -1.39, -0.89] }, // H12
      { element: 'H', position: [0, 4.17, 0.89] }, // H13
      { element: 'H', position: [0, 4.17, -0.89] }, // H14
    ],
    bonds: [
      // Ring bonds
      { from: 0, to: 1, order: 1 },
      { from: 1, to: 2, order: 1 },
      { from: 2, to: 3, order: 1 },
      { from: 3, to: 4, order: 1 },
      { from: 4, to: 5, order: 1 },
      { from: 5, to: 0, order: 1 },
      { from: 0, to: 6, order: 1 }, // C1-O1
      // Hydroxyl groups
      { from: 4, to: 7, order: 1 }, // C5-O2
      // Hydrogens
      { from: 0, to: 8, order: 1 },
      { from: 0, to: 9, order: 1 },
      { from: 1, to: 10, order: 1 },
      { from: 2, to: 11, order: 1 },
      { from: 2, to: 12, order: 1 },
      { from: 3, to: 13, order: 1 },
      { from: 3, to: 14, order: 1 },
      { from: 4, to: 15, order: 1 },
      { from: 5, to: 16, order: 1 },
      { from: 5, to: 17, order: 1 },
      { from: 6, to: 18, order: 1 },
      { from: 6, to: 19, order: 1 },
      { from: 7, to: 20, order: 1 },
      { from: 7, to: 21, order: 1 },
    ],
  },
];

/**
 * Create a MoleculeGraph from a template
 */
export function createMoleculeFromTemplate(template: MoleculeTemplate): MoleculeGraph {
  const molecule = new MoleculeGraph();
  const atomIds: string[] = [];
  
  // Add atoms
  for (const atom of template.atoms) {
    const id = molecule.addAtom({
      element: atom.element,
      position: atom.position,
    });
    atomIds.push(id);
  }
  
  // Add bonds
  for (const bond of template.bonds) {
    const fromId = atomIds[bond.from];
    const toId = atomIds[bond.to];
    if (fromId && toId) {
      molecule.addBond(fromId, toId, bond.order);
    }
  }
  
  return molecule;
}

/**
 * Get template by name
 */
export function getTemplateByName(name: string): MoleculeTemplate | undefined {
  return moleculeTemplates.find((t) => t.name.toLowerCase() === name.toLowerCase());
}

