import { describe, it, expect, vi } from 'vitest';
import React from 'react';
import { render } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import { Suspense } from 'react';
vi.mock('react-router-dom', () => {
  const React = require('react');
  return {
    Link: (props: any) => React.createElement('a', props),
    useLocation: () => ({ pathname: '/' }),
  };
});
vi.mock('../../store/moleculeStore', () => {
  const listeners: any[] = [];
  const base: any = {
    selectedAtomId: null,
    tool: 'select',
    subscribe: (fn: any) => {
      listeners.push(fn);
      return () => {};
    },
    getState: () => base,
  };
  const hook: any = (selector?: any) => (selector ? selector(base) : base);
  hook.getState = () => base;
  hook.setState = (partial: any) => Object.assign(base, typeof partial === 'function' ? partial(base) : partial);
  hook.subscribe = base.subscribe;
  return { useMoleculeStore: hook };
});
import AtomMesh from '../AtomMesh';

describe('AtomMesh', () => {
  it('renders without throwing', () => {
    const { container } = render(
      <MemoryRouter>
        <Suspense fallback={null}>
          <AtomMesh id="a1" position={[0,0,0]} element="C" />
        </Suspense>
      </MemoryRouter>
    );
    expect(container).toBeTruthy();
  });
});


