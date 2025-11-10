import { describe, it, expect, beforeEach, vi } from 'vitest'
import React from 'react'
import { renderToString } from 'react-dom/server'
vi.mock('react-router-dom', () => {
  const React = require('react');
  return {
    Link: (props: any) => React.createElement('a', props),
    useLocation: () => ({ pathname: '/lab' }),
  };
})
vi.mock('../../store/moleculeStore', () => {
  const base: any = {
    tool: 'select',
    setTool: (t: any) => { base.tool = t },
    currentBondOrder: 1,
    setBondOrder: (o: any) => { base.currentBondOrder = o },
    reset: () => {},
  };
  const hook: any = (selector?: any) => (selector ? selector(base) : base);
  hook.getState = () => base;
  hook.setState = () => {};
  return { useMoleculeStore: hook };
})
vi.mock('../../store/historyStore', () => {
  const base = { canUndo: false, canRedo: false };
  const hook: any = (selector?: any) => (selector ? selector(base) : base);
  return { useHistoryStore: hook, undo: () => {}, redo: () => {} };
})
import ToolPanel from '../ToolPanel'
import { useMoleculeStore } from '../../store/moleculeStore'

describe('ToolPanel', () => {
  beforeEach(() => {
    // Reset store
    useMoleculeStore.getState().reset()
  })

  it('renders tool panel', () => {
    const html = renderToString(<ToolPanel />)
    expect(typeof html).toBe('string')
    expect(html.length).toBeGreaterThan(0)
  })

  it('highlights active tool', () => {
    useMoleculeStore.getState().setTool('add-atom')
    renderToString(<ToolPanel />)
    expect(useMoleculeStore.getState().tool).toBe('add-atom')
  })
})

