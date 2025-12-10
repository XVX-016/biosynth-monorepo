import React, { createContext, useContext, useReducer } from "react";

/**
 * Lightweight EditorContext used by LabV2 renderer & toolbar
 * Atoms: { id, element, x,y,z }
 * Bonds: { id, a, b, order }
 */

export type Atom = { id: string; element: string; x: number; y: number; z: number };
export type Bond = { id: string; a: string; b: string; order: number };

type Tool = "select" | "add-atom" | "add-bond" | "move" | "erase";

type State = {
  atoms: Atom[];
  bonds: Bond[];
  tool: Tool;
  selection?: { type: "atom" | "bond"; id: string } | null;
  autoBond: boolean;
  busy: boolean;
  history: { atoms: Atom[]; bonds: Bond[] }[];
  future: { atoms: Atom[]; bonds: Bond[] }[];
  name?: string;
};

const initialState: State = {
  atoms: [],
  bonds: [],
  tool: "select",
  selection: null,
  autoBond: false,
  busy: false,
  history: [],
  future: [],
  name: "untitled"
};

type Action =
  | { type: "SET_TOOL"; payload: Tool }
  | { type: "ADD_ATOM"; payload: Atom }
  | { type: "MOVE_ATOM"; payload: Atom }
  | { type: "DELETE_ATOM"; payload: string }
  | { type: "ADD_BOND"; payload: Bond }
  | { type: "SET_MOLECULE"; payload: { atoms: Atom[]; bonds: Bond[]; name?: string } }
  | { type: "AUTO_BOND" }
  | { type: "APPLY_ML"; payload: { correctedAtoms: { id: string; x: number; y: number; z: number }[] } }
  | { type: "UNDO" }
  | { type: "REDO" }
  | { type: "TOGGLE_AUTOBOND" }
  | { type: "SET_BUSY"; payload: boolean }
  | { type: "SET_SELECTION"; payload: State["selection"] }
  | { type: "CLEAR_HISTORY" }
  ;

const EditorContext = createContext<any>(null);

function snapshot(state: State) {
  return { atoms: JSON.parse(JSON.stringify(state.atoms)), bonds: JSON.parse(JSON.stringify(state.bonds)) };
}

/**
 * Simple bond engine used inline for AUTO_BOND
 * (distance-based + valence)
 */
const covalentRadii: Record<string, number> = {
  H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
  P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39,
};
const maxValence: Record<string, number> = { H: 1, C: 4, N: 3, O: 2, F: 1, S: 6, P: 5, Cl: 1, Br: 1, I: 1 };
const BOND_TOL = 0.45;
const MIN_BOND_DISTANCE = 0.3;

function dist(a: Atom, b: Atom) {
  return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2);
}

function currentValence(atomId: string, bonds: Bond[]) {
  return bonds.filter(b => b.a === atomId || b.b === atomId).reduce((s, bd) => s + (bd.order || 1), 0);
}

function canAddBond(a: Atom, b: Atom, bonds: Bond[]) {
  const r1 = covalentRadii[a.element], r2 = covalentRadii[b.element];
  if (!r1 || !r2) return false;
  const ideal = r1 + r2;
  const d = dist(a, b);
  if (d < MIN_BOND_DISTANCE) return false;
  if (d > ideal + BOND_TOL) return false;
  const maxA = maxValence[a.element] ?? 0;
  const maxB = maxValence[b.element] ?? 0;
  if (currentValence(a.id, bonds) >= maxA) return false;
  if (currentValence(b.id, bonds) >= maxB) return false;
  if (bonds.some(bd => (bd.a === a.id && bd.b === b.id) || (bd.a === b.id && bd.b === a.id))) return false;
  return true;
}

function reducer(state: State, action: Action): State {
  switch (action.type) {
    case "SET_TOOL": return { ...state, tool: action.payload };
    case "ADD_ATOM": {
      const next = { ...state, atoms: [...state.atoms, action.payload] };
      return { ...next, history: [...state.history, snapshot(state)], future: [] };
    }
    case "MOVE_ATOM": {
      const atoms = state.atoms.map(a => a.id === action.payload.id ? action.payload : a);
      return { ...state, atoms, history: [...state.history, snapshot(state)], future: [] };
    }
    case "DELETE_ATOM": {
      const atoms = state.atoms.filter(a => a.id !== action.payload);
      const bonds = state.bonds.filter(b => b.a !== action.payload && b.b !== action.payload);
      return { ...state, atoms, bonds, history: [...state.history, snapshot(state)], future: [] };
    }
    case "ADD_BOND": {
      return { ...state, bonds: [...state.bonds, action.payload], history: [...state.history, snapshot(state)], future: [] };
    }
    case "SET_MOLECULE": {
      const { atoms, bonds, name } = action.payload;
      return { ...state, atoms: atoms || [], bonds: bonds || [], name: name ?? state.name, history: [...state.history, snapshot(state)], future: [] };
    }
    case "AUTO_BOND": {
      const atoms = state.atoms;
      const bonds = state.bonds.slice();
      const newBonds: Bond[] = [];
      for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
          const a = atoms[i], b = atoms[j];
          if (canAddBond(a, b, bonds)) {
            const nb: Bond = { id: crypto.randomUUID(), a: a.id, b: b.id, order: 1 };
            newBonds.push(nb);
            bonds.push(nb);
          }
        }
      }
      if (newBonds.length === 0) return state;
      return { ...state, bonds: [...state.bonds, ...newBonds], history: [...state.history, snapshot(state)], future: [] };
    }
    case "APPLY_ML": {
      // Apply ML-corrected atom positions
      const updates = new Map(action.payload.correctedAtoms.map(a => [a.id, a]));
      const atoms = state.atoms.map(atom => {
        const update = updates.get(atom.id);
        if (update) {
          return { ...atom, x: update.x, y: update.y, z: update.z };
        }
        return atom;
      });
      return { ...state, atoms, history: [...state.history, snapshot(state)], future: [] };
    }
    case "UNDO": {
      const h = state.history.slice();
      if (h.length === 0) return state;
      const last = h[h.length - 1];
      const newHistory = h.slice(0, -1);
      const future = [snapshot(state), ...state.future];
      return { ...state, atoms: last.atoms, bonds: last.bonds, history: newHistory, future };
    }
    case "REDO": {
      if (state.future.length === 0) return state;
      const f = state.future[0];
      const newFuture = state.future.slice(1);
      return { ...state, atoms: f.atoms, bonds: f.bonds, history: [...state.history, snapshot(state)], future: newFuture };
    }
    case "TOGGLE_AUTOBOND": return { ...state, autoBond: !state.autoBond };
    case "SET_BUSY": return { ...state, busy: action.payload };
    case "SET_SELECTION": return { ...state, selection: action.payload };
    case "CLEAR_HISTORY": return { ...state, history: [], future: [] };
    default: return state;
  }
}

export function EditorProvider({ children }: any) {
  const [state, dispatch] = useReducer(reducer, initialState);
  return <EditorContext.Provider value={{ state, dispatch }}>{children}</EditorContext.Provider>;
}

export const useEditorContext = () => useContext(EditorContext);

export const useEditor = useEditorContext;
