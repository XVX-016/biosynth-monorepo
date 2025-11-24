import type { PointerEvent as R3FPointerEvent } from '@react-three/fiber'
import type { useLabStore as _ } from '../store/labStore'
import type { ToolName } from '../types/molecule'

export interface Tool {
  name: ToolName
  onPointerDown?: (ev: any, store: any) => void
  onPointerMove?: (ev: any, store: any) => void
  onPointerUp?: (ev: any, store: any) => void
}

