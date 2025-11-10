import { describe, it, expect, vi } from 'vitest'
import React from 'react'
import { renderToString } from 'react-dom/server'
vi.mock('react-router-dom', () => {
  const React = require('react');
  return {
    Link: (props: any) => React.createElement('a', props),
    useLocation: () => ({ pathname: '/' }),
  };
})
vi.mock('../../store/moleculeStore', () => {
  const base = {
    currentMolecule: null,
    fetchPredictions: () => {},
    tool: 'select',
    selectedAtomId: null,
    selectedBondId: null,
  };
  const hook: any = (selector?: any) => (selector ? selector(base) : base);
  hook.getState = () => base;
  hook.setState = () => {};
  return { useMoleculeStore: hook };
})
import MoleculeViewer from '../MoleculeViewer'

describe('MoleculeViewer', () => {
  it('renders without crashing', () => {
    const html = renderToString(<MoleculeViewer />)
    expect(typeof html).toBe('string')
    expect(html.length).toBeGreaterThan(0)
  })
})

