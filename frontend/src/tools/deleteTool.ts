import type { Tool } from './toolInterface'

const deleteTool: Tool = {
  name: 'delete',
  onPointerDown: (ev: any, store: any) => {
    const picked = ev.object
    if (!picked) return
    
    const atomId = picked.userData?.atomId
    const bondId = picked.userData?.bondId
    
    if (atomId) {
      store.deleteAtom(atomId)
    } else if (bondId) {
      store.deleteBond(bondId)
    }
  }
}

export default deleteTool

