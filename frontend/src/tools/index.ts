import selectTool from './selectTool'
import addAtomTool from './addAtomTool'
import bondTool from './bondTool'
import deleteTool from './deleteTool'
import type { Tool } from './toolInterface'
import type { ToolName } from '../types/molecule'

export const tools: Record<ToolName, Tool> = {
  select: selectTool,
  add_atom: addAtomTool,
  bond: bondTool,
  delete: deleteTool,
  move: selectTool, // For now, move uses select behavior
  inspect: selectTool, // For now, inspect uses select behavior
}

export function getTool(name: ToolName): Tool {
  return tools[name] || tools.select
}

