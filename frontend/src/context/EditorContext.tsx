import React, { createContext, useContext, useReducer, useRef, useEffect } from "react";
import UndoRedo, { type Snapshot } from "../lab/editor/undoRedo";

export type Atom = { id: string; element: string; position: [number, number, number] };
export type Bond = { id: string; a: string; b: string; order: number };
export type Tool = "select" | "add-atom" | "add-bond" | "move" | "erase";

type State = {
  atoms: Atom[];
  bonds: Bond[];
  tool: Tool;
  name?: string;
  autoBond: boolean;
  busy: boolean;
  selection?: string;
  predictions?: Record<string, any>;
};

const initialState: State = {
  atoms: [],
  bonds: [],
  tool: "select",
  name: "untitled",
  autoBond: false,
  busy: false,
  selection: undefined,
  predictions: {}
};

type Action =
  | { type: "SET_TOOL"; payload: Tool }
  | { type: "ADD_ATOM"; payload: Atom }
  | { type: "MOVE_ATOM"; payload: Atom } // Payload has full atom data
  | { type: "DELETE_ATOM"; payload: string }
  | { type: "ADD_BOND"; payload: Bond }
  | { type: "SET_SELECTION"; payload?: string }
  | { type: "LOAD_MOLECULE"; payload: { atoms: Atom[]; bonds: Bond[]; name?: string } }
  | { type: "AUTO_BOND" } // Handled via side effect or external helper? Ideally in component or thunk. Here we just update state if passed data.
  // Actually, for AUTO_BOND to work in reducer it needs logic. 
  // We'll trust the BondEngine usage in components for now, or assume this action might be used by a saga/effect.
  // For now, let's keep it simple: components call BondEngine and dispatch LOAD_MOLECULE or ADD_BOND.
  | { type: "TOGGLE_AUTOBOND" }
  | { type: "SET_BUSY"; payload: boolean }
  | { type: "SET_PREDICTIONS"; payload: Record<string, any> }
  | { type: "APPLY_ML"; payload: { correctedAtoms: { id: string; x: number; y: number; z: number }[] } };

function reducer(state: State, action: Action): State {
  switch (action.type) {
    case "SET_TOOL": return { ...state, tool: action.payload };
    case "ADD_ATOM": return { ...state, atoms: [...state.atoms, action.payload] };
    case "MOVE_ATOM": return {
      ...state,
      atoms: state.atoms.map(a => a.id === action.payload.id ? action.payload : a)
    };
    case "DELETE_ATOM": {
      const id = action.payload;
      return {
        ...state,
        atoms: state.atoms.filter(a => a.id !== id),
        bonds: state.bonds.filter(b => b.a !== id && b.b !== id)
      };
    }
    case "ADD_BOND": return { ...state, bonds: [...state.bonds, action.payload] };
    case "SET_SELECTION": return { ...state, selection: action.payload };
    case "LOAD_MOLECULE":
      return {
        ...state,
        atoms: action.payload.atoms,
        bonds: action.payload.bonds,
        name: action.payload.name ?? state.name
      };
    case "TOGGLE_AUTOBOND": return { ...state, autoBond: !state.autoBond };
    case "SET_BUSY": return { ...state, busy: action.payload };
    case "SET_PREDICTIONS": return { ...state, predictions: action.payload };
    case "APPLY_ML": {
      const updates = new Map(action.payload.correctedAtoms.map((a) => [a.id, a]));
      return {
        ...state,
        atoms: state.atoms.map((atom) => {
          const update = updates.get(atom.id);
          if (update) {
            return { ...atom, position: [update.x, update.y, update.z] };
          }
          return atom;
        }),
      };
    }
    default: return state;
  }
}

const EditorContext = createContext<{
  state: State;
  dispatch: React.Dispatch<Action>;
  undo: () => void;
  redo: () => void;
  canUndo: () => boolean;
  canRedo: () => boolean;
}>({ state: initialState, dispatch: () => null, undo: () => { }, redo: () => { }, canUndo: () => false, canRedo: () => false });

export const EditorProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [state, dispatch] = useReducer(reducer, initialState);
  const history = useRef(new UndoRedo(100));
  const suppressPush = useRef(false);

  useEffect(() => {
    if (suppressPush.current) {
      suppressPush.current = false;
      return;
    }
    history.current.push({ atoms: state.atoms, bonds: state.bonds } as Snapshot);
  }, [state.atoms, state.bonds]);

  const undo = () => {
    if (!history.current.canUndo()) return;
    const snap = history.current.undo();
    if (!snap) return;
    suppressPush.current = true;
    dispatch({ type: 'LOAD_MOLECULE', payload: { atoms: snap.atoms, bonds: snap.bonds } });
  };

  const redo = () => {
    if (!history.current.canRedo()) return;
    const snap = history.current.redo();
    if (!snap) return;
    suppressPush.current = true;
    dispatch({ type: 'LOAD_MOLECULE', payload: { atoms: snap.atoms, bonds: snap.bonds } });
  };

  return (
    <EditorContext.Provider value={{ state, dispatch, undo, redo, canUndo: () => history.current.canUndo(), canRedo: () => history.current.canRedo() }}>
      {children}
    </EditorContext.Provider>
  );
};

export const useEditor = () => useContext(EditorContext);
