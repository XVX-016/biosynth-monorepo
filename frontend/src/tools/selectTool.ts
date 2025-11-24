import type { Tool } from './toolInterface'

const selectTool: Tool = {
  name: 'select',
  onPointerDown: (ev: any, store: any) => {
    // expecting ev.object is the clicked mesh and ev.instanceId if instanced
    const picked = ev.object
    if (!picked) {
      store.setSelectedAtomId(null)
      store.setSelectedBondId(null)
      return
    }
    // if you use instancing, map instanceId to atom id map
    const atomId = (picked.userData && picked.userData.atomId) ? picked.userData.atomId : null
    const bondId = (picked.userData && picked.userData.bondId) ? picked.userData.bondId : null
    
    if (atomId) {
      store.setSelectedAtomId(atomId)
      store.setSelectedBondId(null)
    } else if (bondId) {
      store.setSelectedBondId(bondId)
      store.setSelectedAtomId(null)
    } else {
      store.setSelectedAtomId(null)
      store.setSelectedBondId(null)
    }
  }
}

export default selectTool

