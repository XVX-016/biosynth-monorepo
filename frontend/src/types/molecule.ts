// canonical types for the lab
export type Vec3 = [number, number, number];

export type Atom = {
  id: string;               // unique id (nanoid)
  element: string;          // 'C', 'O', 'H', ...
  position: Vec3;           // [x,y,z]
  charge?: number;          // formal charge (optional)
  mass?: number;            // optional numeric mass
  meta?: Record<string, any>; // extensible metadata
};

export type Bond = {
  id: string;
  atom1: string;            // atom.id
  atom2: string;            // atom.id
  order: 1 | 2 | 3;         // bond order
  meta?: Record<string, any>;
};

export type Molecule = {
  id?: string;              // optional DB id
  name?: string;
  atoms: Atom[];
  bonds: Bond[];
  metadata?: Record<string, any>;
};

// Tool enum
export type ToolName = 'select' | 'add_atom' | 'bond' | 'move' | 'delete' | 'inspect';

