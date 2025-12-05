import { createContext, useContext, useReducer } from 'react';
import type { Dispatch, ReactNode } from 'react';
type Action = 
  | { type: 'SET_MOLECULE'; payload: any }       // replace 'any' with real molecule type
  | { type: 'SET_SELECTION'; payload: any }     // replace 'any' with real selection type
  | { type: 'SET_TOOL'; payload: any }          // replace 'any' with real tool type
  | { type: 'SET_STATE'; payload: State }       // ✅ REQUIRED for Undo/Redo + restore
  | { type: 'RUN_PREDICTION'; payload?: any }
  | { type: 'UNDO' }
  | { type: 'REDO' };

type State = {
  molecule: any;        // replace with real molecule type
  stableSmiles: string | null;
  selection: any;       // replace with real selection type
  currentTool: any;     // replace with real tool type
  predictions: any[];  // replace with real prediction type
  history: { action: Action; state: State }[];
};

const initialState: State = {
  molecule: null,
  stableSmiles: '',
  selection: null,
  currentTool: null,
  predictions: [],
  history: []
};
const EditorContext = createContext<{
  state: State;
  dispatch: Dispatch<Action>;
} | undefined>(undefined);

function editorReducer(state: State, action: Action): State {
  switch (action.type) {

    case 'SET_MOLECULE':
      return { ...state, molecule: action.payload };

    case 'SET_SELECTION':
      return { ...state, selection: action.payload };

    case 'SET_TOOL':
      return { ...state, currentTool: action.payload };

    case 'SET_STATE':
      return action.payload; 

    case 'RUN_PREDICTION': {
      const newPredictions = [...state.predictions, action.payload];
      return { ...state, predictions: newPredictions };
    }

    case 'UNDO': {
      if (state.history.length === 0) return state;
      const previous = state.history[state.history.length - 1];
      const newHistory = state.history.slice(0, -1);
      return { ...previous.state, history: newHistory };
    }

    case 'REDO': {
      if (state.history.length === 0) return state;
      const next = state.history[0];
      const newHistory = state.history.slice(1);
      return { ...next.state, history: newHistory };
    }

    default:
      throw new Error(`Unhandled action type`);
  }
}

export function EditorProvider({ children }: { children: ReactNode }) {
  const [state, dispatch] = useReducer(editorReducer, initialState);

  return (
    <EditorContext.Provider value={{ state, dispatch }}>
      {children}
    </EditorContext.Provider>
  );
}
export function useEditor() {
  const context = useContext(EditorContext);

  if (!context) {
    throw new Error('useEditor must be used within an EditorProvider');
  }

  return context;
}
