import type { Tool } from './toolInterface'

let pendingAtomId: string | null = null

const bondTool: Tool = {
  name: 'bond',
  onPointerDown: (ev: any, store: any) => {
    const pickedAtomId = ev.object?.userData?.atomId
    if (!pickedAtomId) {
      pendingAtomId = null
      return
    }
    
    if (!pendingAtomId) {
      pendingAtomId = pickedAtomId
      store.setSelectedAtomId(pickedAtomId)
      // highlight visually via store or userData
    } else {
      if (pendingAtomId !== pickedAtomId) {
        store.addBond(pendingAtomId, pickedAtomId)
      }
      pendingAtomId = null
      store.setSelectedAtomId(null)
    }
  },
  onPointerUp: () => {
    // Reset on pointer up if needed
  }
}

export default bondTool

