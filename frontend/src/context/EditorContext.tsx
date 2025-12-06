import React, { createContext, useContext, useReducer } from "react";
import type { Atom, Bond } from "../types/ml";

type Tool = "select" | "add-atom" | "add-bond" | "move" | "erase";

type State = {
  atoms: Atom[];
  bonds: Bond[];
  tool: Tool;
  name?: string;
  autoBond: boolean;
  busy: boolean;
  predictions?: Record<string, any>;
  history: { atoms: Atom[]; bonds: Bond[] }[];
  future: { atoms: Atom[]; bonds: Bond[] }[];
};

const initialState: State = {
  atoms: [],
  bonds: [],
  tool: "select",
  name: "untitled",
  autoBond: false,
  busy: false,
  predictions: {},
  history: [],
  future: []
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
  | { type: "TOGGLE_AUTOBOND" }
  | { type: "SET_BUSY"; payload: boolean }
  | { type: "SET_PREDICTIONS"; payload: Record<string, any> }
  | { type: "UNDO" }
  | { type: "REDO" }
  | { type: "CLEAR_HISTORY" }
  ;

const EditorContext = createContext<any>(null);

function snapshot(state: State) {
  return { atoms: JSON.parse(JSON.stringify(state.atoms)), bonds: JSON.parse(JSON.stringify(state.bonds)) };
}

function editorReducer(state: State, action: Action): State {
  switch (action.type) {
    case "SET_TOOL": return { ...state, tool: action.payload };
    case "ADD_ATOM": {
      const next = { ...state, atoms: [...state.atoms, action.payload] };
      return { ...next, history: [...state.history, snapshot(state)], future: [] };
    }
    case "MOVE_ATOM": {
      const newAtoms = state.atoms.map(a => a.id === action.payload.id ? { ...a, ...action.payload } : a);
      return { ...state, atoms: newAtoms, history: [...state.history, snapshot(state)], future: [] };
    }
    case "DELETE_ATOM": {
      const id = action.payload;
      const newAtoms = state.atoms.filter(a => a.id !== id);
      const newBonds = state.bonds.filter(b => b.a !== id && b.b !== id);
      return { ...state, atoms: newAtoms, bonds: newBonds, history: [...state.history, snapshot(state)], future: [] };
    }
    case "ADD_BOND": {
      const newBond = action.payload;
      return { ...state, bonds: [...state.bonds, newBond], history: [...state.history, snapshot(state)], future: [] };
    }
    case "SET_MOLECULE": {
      const { atoms, bonds, name } = action.payload;
      // Reset history on new molecule load
      return { ...state, atoms, bonds, name: name ?? state.name, history: [], future: [] };
    }
    case "AUTO_BOND": {
      // NOTE: Dynamic import of BondEngine to avoid circular dependency if any.
      // But we can import at top if bundler allows. Let's assume top import is cleaner for now.
      // Or follow user advice with require.
      // "const { BondEngine } = require("../engine/bonding/BondEngine");"
      // In TS with ESM, require is iffy without proper config.
      // I will use import from top, as I will add import.
      // But user specifically said "to avoid circular imports, require dynamically".
      // I'll trust them.
      // eslint-disable-next-line @typescript-eslint/no-var-requires
      // const { BondEngine } = require("../engine/bonding/BondEngine");
      // Actually, BondEngine.ts is leaf node, no circular dep with Context unless Context uses BondEngine and BondEngine uses Context.
      // BondEngine doesn't use Context. Safe to import.
      // BUT if I paste user code EXACTLY...
      // I will just use the logic from BondEngine.
      // Let's import BondEngine at the top to be standard.
      const { BondEngine } = require("../engine/bonding/BondEngine");
      const newBonds = BondEngine.autoBond(state.atoms, state.bonds);
      if (newBonds.length === 0) return state;
      return { ...state, bonds: [...state.bonds, ...newBonds], history: [...state.history, snapshot(state)], future: [] };
    }
    case "APPLY_ML": {
      const corrected = action.payload.correctedAtoms;
      const newAtoms = state.atoms.map(a => {
        const fix = corrected.find((c: any) => c.id === a.id);
        return fix ? { ...a, x: fix.x, y: fix.y, z: fix.z } : a;
      });
      return { ...state, atoms: newAtoms, history: [...state.history, snapshot(state)], future: [] };
    }
    case "TOGGLE_AUTOBOND": return { ...state, autoBond: !state.autoBond };
    case "SET_BUSY": return { ...state, busy: action.payload };
    case "SET_PREDICTIONS": return { ...state, predictions: action.payload };
    case "UNDO": {
      const h = state.history.slice();
      if (h.length === 0) return state;
      const last = h[h.length - 1];
      const newHistory = h.slice(0, -1);
      return { ...state, atoms: last.atoms, bonds: last.bonds, history: newHistory, future: [snapshot(state), ...state.future] };
    }
    case "REDO": {
      const f = state.future.slice();
      if (f.length === 0) return state;
      const nextSnap = f[0];
      const newFuture = f.slice(1);
      return { ...state, atoms: nextSnap.atoms, bonds: nextSnap.bonds, history: [...state.history, snapshot(state)], future: newFuture };
    }
    case "CLEAR_HISTORY": return { ...state, history: [], future: [] };
    default: return state;
  }
}

export function EditorProvider({ children }: { children: React.ReactNode }) {
  const [state, dispatch] = useReducer(editorReducer, initialState);
  return <EditorContext.Provider value={{ state, dispatch }}>{children}</EditorContext.Provider>;
}

export const useEditorContext = () => {
  const context = useContext(EditorContext);
  if (!context) throw new Error("useEditorContext must be used within EditorProvider");
  return context;
};

// Also export useEditor for compatibility if used elsewhere (user's Lab.tsx doesn't use it, but maybe older components do)
export const useEditor = useEditorContext;
