import { useEffect } from 'react';
import { useEditor } from '../context/EditorContext';
import { PredictionService } from '../lib/molecule/prediction/PredictionService.ts'; // replace with actual service

export default function useMoleculeEditor() {
  const { state, dispatch } = useEditor();
  
  const setTool = (tool) => dispatch({ type: 'SET_TOOL', payload: tool });

  const selectMolecule = (molecule) => {
    dispatch({ type: 'SET_MOLECULE', payload: molecule });
    localStorage.setItem('molecule', JSON.stringify(molecule)); // save to local storage
  };
  
  const selectSelection = (selection) => dispatch({ type: 'SET_SELECTION', payload: selection });

  const runPrediction = async () => {
    try {
      const predictionResult = await PredictionService.run(state.molecule); // replace with actual service call
      dispatch({ type: 'RUN_PREDICTION', payload: predictionResult });
    } catch (error) {
      console.error('Error running prediction: ', error);
    }
  };
  
  const undo = () => dispatch({ type: 'UNDO' });
  const redo = () => dispatch({ type: 'REDO' });

  // Autosave every 10 seconds
  useEffect(() => {
    const intervalId = setInterval(() => {
      try {
        localStorage.setItem('state', JSON.stringify(state));
      } catch (error) {
        console.error('Error saving state to local storage: ', error);
      }
    }, 10000);
    
    return () => clearInterval(intervalId); // cleanup on unmount
  }, [state]);
  
  // Restore molecule from LocalStorage on load
  useEffect(() => {
    const savedState = localStorage.getItem('state');
    if (savedState) {
      dispatch({ type: 'SET_STATE', payload: JSON.parse(savedState) });
    }
  }, []); // run only on mount and unmount
  
  return {
    ...state,
    setTool,
    selectMolecule,
    selectSelection,
    runPrediction,
    undo,
    redo,
  };
}