import { Atom, Element } from "./types";
import { MoleculeGraph } from "./MoleculeGraph";
import { getBondLength, ELEMENT_RADII } from "./types";

/**
 * Places a new atom relative to an existing atom at an appropriate distance
 */
export function placeAtomRelative(
  baseAtom: Atom,
  newElement: Element,
  direction?: [number, number, number]
): [number, number, number] {
  const bondLength = getBondLength(baseAtom.element, newElement);
  const baseRadius = ELEMENT_RADII[baseAtom.element];
  const newRadius = ELEMENT_RADII[newElement];

  // Default direction: random but normalized
  let dir: [number, number, number];
  if (direction) {
    const len = Math.sqrt(
      direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2
    );
    dir = [direction[0] / len, direction[1] / len, direction[2] / len];
  } else {
    // Generate random direction
    const theta = Math.random() * Math.PI * 2;
    const phi = Math.acos(2 * Math.random() - 1);
    dir = [
      Math.sin(phi) * Math.cos(theta),
      Math.sin(phi) * Math.sin(theta),
      Math.cos(phi),
    ];
  }

  const distance = bondLength + baseRadius + newRadius;
  return [
    baseAtom.position[0] + dir[0] * distance,
    baseAtom.position[1] + dir[1] * distance,
    baseAtom.position[2] + dir[2] * distance,
  ];
}

/**
 * Optimizes geometry using basic force field principles
 */
export function optimizeGeometry(
  molecule: MoleculeGraph,
  iterations: number = 30
): void {
  const atoms = Array.from(molecule.atoms.values());
  const velocities: Map<string, [number, number, number]> = new Map();

  // Initialize velocities
  atoms.forEach((atom) => {
    velocities.set(atom.id, [0, 0, 0]);
  });

  const k = 0.1; // Spring constant
  const damping = 0.9;

  for (let iter = 0; iter < iterations; iter++) {
    const forces: Map<string, [number, number, number]> = new Map();
    atoms.forEach((atom) => {
      forces.set(atom.id, [0, 0, 0]);
    });

    // Bond stretching forces
    molecule.bonds.forEach((bond) => {
      const a1 = molecule.atoms.get(bond.a1)!;
      const a2 = molecule.atoms.get(bond.a2)!;
      const desiredLength = getBondLength(a1.element, a2.element);

      const dx = a2.position[0] - a1.position[0];
      const dy = a2.position[1] - a1.position[1];
      const dz = a2.position[2] - a1.position[2];
      const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

      if (distance > 0.01) {
        const force = k * (distance - desiredLength);
        const fx = (force * dx) / distance;
        const fy = (force * dy) / distance;
        const fz = (force * dz) / distance;

        const f1 = forces.get(bond.a1)!;
        const f2 = forces.get(bond.a2)!;
        forces.set(bond.a1, [f1[0] - fx, f1[1] - fy, f1[2] - fz]);
        forces.set(bond.a2, [f2[0] + fx, f2[1] + fy, f2[2] + fz]);
      }
    });

    // Non-bonded repulsion (simplified Lennard-Jones)
    const epsilon = 0.1;
    const sigma = 3.0;
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const a1 = atoms[i];
        const a2 = atoms[j];

        // Skip if bonded
        const isBonded = Array.from(molecule.bonds.values()).some(
          (b) =>
            (b.a1 === a1.id && b.a2 === a2.id) ||
            (b.a1 === a2.id && b.a2 === a1.id)
        );
        if (isBonded) continue;

        const dx = a2.position[0] - a1.position[0];
        const dy = a2.position[1] - a1.position[1];
        const dz = a2.position[2] - a1.position[2];
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (distance > 0.01 && distance < 5.0) {
          const r6 = (sigma / distance) ** 6;
          const r12 = r6 * r6;
          const force = (24 * epsilon / distance) * (2 * r12 - r6);

          const fx = (force * dx) / distance;
          const fy = (force * dy) / distance;
          const fz = (force * dz) / distance;

          const f1 = forces.get(a1.id)!;
          const f2 = forces.get(a2.id)!;
          forces.set(a1.id, [f1[0] - fx, f1[1] - fy, f1[2] - fz]);
          forces.set(a2.id, [f2[0] + fx, f2[1] + fy, f2[2] + fz]);
        }
      }
    }

    // Update positions
    const dt = 0.1;
    atoms.forEach((atom) => {
      const force = forces.get(atom.id)!;
      const velocity = velocities.get(atom.id)!;

      const newVel: [number, number, number] = [
        (velocity[0] + force[0] * dt) * damping,
        (velocity[1] + force[1] * dt) * damping,
        (velocity[2] + force[2] * dt) * damping,
      ];

      velocities.set(atom.id, newVel);

      const newPos: [number, number, number] = [
        atom.position[0] + newVel[0] * dt,
        atom.position[1] + newVel[1] * dt,
        atom.position[2] + newVel[2] * dt,
      ];

      molecule.atoms.set(atom.id, { ...atom, position: newPos });
    });
  }
}

/**
 * Centers the molecule at origin
 */
export function centerMolecule(molecule: MoleculeGraph): void {
  const atoms = Array.from(molecule.atoms.values());
  if (atoms.length === 0) return;

  let sumX = 0,
    sumY = 0,
    sumZ = 0;
  atoms.forEach((atom) => {
    sumX += atom.position[0];
    sumY += atom.position[1];
    sumZ += atom.position[2];
  });

  const centerX = sumX / atoms.length;
  const centerY = sumY / atoms.length;
  const centerZ = sumZ / atoms.length;

  atoms.forEach((atom) => {
    molecule.atoms.set(atom.id, {
      ...atom,
      position: [
        atom.position[0] - centerX,
        atom.position[1] - centerY,
        atom.position[2] - centerZ,
      ],
    });
  });
}

