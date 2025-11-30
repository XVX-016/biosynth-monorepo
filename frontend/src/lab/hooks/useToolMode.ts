/**
 * useToolMode - Hook for accessing tool mode
 */

import { useLabStore } from '../state/LabStore';

export function useToolMode() {
  const toolManager = useLabStore((s) => s.toolManager);
  const setTool = useLabStore((s) => s.setTool);
  const activeTool = toolManager.getActive();
  
  return {
    toolManager,
    activeTool,
    setTool,
    is: (toolId: string) => activeTool?.id === toolId,
  };
}

