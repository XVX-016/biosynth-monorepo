import { describe, it, expect, beforeEach } from 'vitest'
import { selectionManager } from '../SelectionManager'

describe('SelectionManager', () => {
  beforeEach(() => {
    selectionManager.reset()
  })

  it('tracks hovered atom', () => {
    selectionManager.onHover('atom_1')
    expect(selectionManager.getHoveredAtomId()).toBe('atom_1')
    expect(selectionManager.isHovered('atom_1')).toBe(true)
  })

  it('tracks selected atom', () => {
    selectionManager.onSelect('atom_1')
    expect(selectionManager.getSelectedAtomId()).toBe('atom_1')
    expect(selectionManager.isSelected('atom_1')).toBe(true)
  })

  it('tracks dragging atom', () => {
    selectionManager.startDrag('atom_1')
    expect(selectionManager.getDraggingAtomId()).toBe('atom_1')
    expect(selectionManager.isDragging('atom_1')).toBe(true)
  })

  it('ends drag correctly', () => {
    selectionManager.startDrag('atom_1')
    selectionManager.endDrag()
    expect(selectionManager.getDraggingAtomId()).toBe(null)
    expect(selectionManager.isDragging('atom_1')).toBe(false)
  })

  it('emits events on state change', () => {
    let hoverCalled = false
    let selectCalled = false
    let dragCalled = false

    selectionManager.on('hover', () => {
      hoverCalled = true
    })
    selectionManager.on('select', () => {
      selectCalled = true
    })
    selectionManager.on('drag', () => {
      dragCalled = true
    })

    selectionManager.onHover('atom_1')
    selectionManager.onSelect('atom_1')
    selectionManager.startDrag('atom_1')

    expect(hoverCalled).toBe(true)
    expect(selectCalled).toBe(true)
    expect(dragCalled).toBe(true)
  })

  it('resets all state', () => {
    selectionManager.onHover('atom_1')
    selectionManager.onSelect('atom_1')
    selectionManager.startDrag('atom_1')

    selectionManager.reset()

    expect(selectionManager.getHoveredAtomId()).toBe(null)
    expect(selectionManager.getSelectedAtomId()).toBe(null)
    expect(selectionManager.getDraggingAtomId()).toBe(null)
  })
})

