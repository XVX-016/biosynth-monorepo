import { createContext, useContext, useReducer } from 'react';
import type { Dispatch, ReactNode } from 'react';

type Action = 
  | { type: 'SET_MOLECULE', payload: any } // replace 'any' with the actual molecule type
  | { type: 'SET_SELECTION', payload: any } 
  | { type: 'RUN_PREDICTION', payload?: any }  
  | { type: 'UNDO' }
  | { type: 'REDO' };

type State = {
  molecule: any; // replace 'any' with the actual molecule type
  stableSmiles: string | null;
  selection: any; // replace 'any' with the actual selection type
  currentTool: any; // replace 'any' with the actual tool type
  predictions: any[]; // replace 'any' with the actual prediction type
  history: { action: Action, state: State }[];
};

const initialState: State = {
  molecule: null, // set initial value for molecule
  stableSmiles: '', // set initial value for stable smiles
  selection: null, // set initial value for selection
  currentTool: null, // set initial value for tool
  predictions: [], // set initial value for predictions
  history: []
};

const EditorContext = createContext<{ state: State; dispatch: Dispatch<Action> } | undefined>(undefined);

function editorReducer(state: State, action: Action): State {
  switch (action.type) {
    case 'SET_MOLECULE':
      return { ...state, molecule: action.payload };
    case 'SET_SELECTION':
      return { ...state, selection: action.payload };
    case 'RUN_PREDICTION':
      // logic for running a prediction; replace this with actual implementation
      const newPredictions = [...state.predictions, action.payload];
      return { ...state, predictions: newPredictions };
    case 'UNDO':
      if (state.history.length === 0) {
        return state;
      }
      
      const prevState = state.history[state.history.length - 1].state;
      const newHistory = state.history.slice(0, state.history.length - 1);
      return { ...prevState, history: newHistory };
    case 'REDO':
      if (state.history.length === 0) {
        return state;
      }
      
      const nextState = state.history[0].state;
      const nextAction = state.history[0].action;
      const newRedoHistory = [...state.history];
      newRedoHistory.unshift({ action: nextAction, state: nextState });
      return { ...nextState, history: newRedoHistory };
    default:
      throw new Error();
  }
}

export function EditorProvider({ children }: { children: ReactNode }) {
  const [state, dispatch] = useReducer(editorReducer, initialState);
  
  return <EditorContext.Provider value={{ state, dispatch }}>{children}</EditorContext.Provider>;
}

export function useEditor() {
  const context = useContext(EditorContext);
  
  if (context === undefined) {
    throw new Error('useEditor must be used within a EditorProvider');
  }
  
  return context;
}