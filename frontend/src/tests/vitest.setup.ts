/**
 * vitest.setup.ts
 * Global test setup: robust mocks for React/Three/React-Three-Fiber/Drei and Zustand store stubbing.
 *
 * This file:
 *  - mocks React.useSyncExternalStore to avoid react-dom client quirks in jsdom
 *  - provides a mockStoreFactory and attempts to auto-mock common store module paths
 *  - mocks @react-three/fiber, @react-three/drei, and three primitives
 *  - attaches helpers to globalThis for convenience
 */

import { vi } from 'vitest';
import React from 'react';
import '@testing-library/jest-dom';

// -------------------- React SSR / useSyncExternalStore guard --------------------
/**
 * Harden useSyncExternalStore usage in tests that may call renderers expecting client environment.
 * We fallback to returning the snapshot directly.
 * 
 * NOTE: We don't mock React itself as it can break react-dom. Instead, we only provide
 * a shim for useSyncExternalStore if needed, but let the actual React be used.
 */
// Don't mock React - let it work normally with react-dom
// If useSyncExternalStore issues arise, we can add a polyfill, but for now let React work as-is

// -------------------- THREE / R3F / DREI lightweight mocks --------------------
/**
 * These mocks replace Three.js and R3F primitives with simple DOM-friendly components.
 * Keeps the shape of the API used by the app while avoiding WebGL.
 */
vi.mock('@react-three/fiber', async () => {
  const React = await import('react');
  // minimal wrapper: Canvas renders children into a div
  return {
    Canvas: ({ children, ...rest }: any) => {
      // NOTE: React elements from three primitives will be plain objects — fine for tests
      return React.createElement('div', { 'data-testid': 'canvas', ...rest }, children);
    },
    useFrame: vi.fn((callback) => {
      // Simulate frame updates
      if (typeof callback === 'function') {
        callback({ clock: { elapsedTime: 0 } });
      }
    }),
    useThree: () => ({ camera: {}, gl: {} }),
    extend: () => {}
  };
});

vi.mock('@react-three/drei', async () => {
  const React = await import('react');
  return {
    OrbitControls: (props: any) => React.createElement('div', { 'data-testid': 'orbit-controls', ...props }, null),
    Html: (props: any) => React.createElement('div', { 'data-testid': 'html', ...props }, null),
    Outlines: ({ color, thickness }: { color: number; thickness: number }) => (
      React.createElement('div', { 'data-testid': 'outlines', 'data-color': color, 'data-thickness': thickness })
    ),
    ContactShadows: (props: any) => React.createElement('div', { 'data-testid': 'contact-shadows', ...props }, null),
    Environment: (props: any) => React.createElement('div', { 'data-testid': 'environment', ...props }, null),
    EffectComposer: ({ children }: { children?: React.ReactNode }) => React.createElement('div', { 'data-testid': 'effect-composer' }, children),
    Bloom: (props: any) => React.createElement('div', { 'data-testid': 'bloom', ...props }, null),
    ChromaticAberration: (props: any) => React.createElement('div', { 'data-testid': 'chromatic-aberration', ...props }, null),
  };
});

// Mock 'three' low-level classes to avoid WebGL calls
vi.mock('three', async () => {
  class MockVector3 {
    x: number; y: number; z: number;
    constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
    set(x: number, y: number, z: number) { this.x = x; this.y = y; this.z = z; return this; }
    clone() { return new MockVector3(this.x, this.y, this.z); }
    subVectors(a: MockVector3, b: MockVector3) { 
      this.x = a.x - b.x; 
      this.y = a.y - b.y; 
      this.z = a.z - b.z; 
      return this; 
    }
    addVectors(a: MockVector3, b: MockVector3) { 
      this.x = a.x + b.x; 
      this.y = a.y + b.y; 
      this.z = a.z + b.z; 
      return this; 
    }
    multiplyScalar(s: number) { 
      this.x *= s; 
      this.y *= s; 
      this.z *= s; 
      return this; 
    }
    length() { return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z); }
    normalize() { 
      const len = this.length(); 
      if (len > 0) { 
        this.x /= len; 
        this.y /= len; 
        this.z /= len; 
      } 
      return this; 
    }
  }
  class MockQuaternion {
    x: number; y: number; z: number; w: number;
    constructor(x = 0, y = 0, z = 0, w = 1) { this.x = x; this.y = y; this.z = z; this.w = w; }
    setFromUnitVectors(v1: MockVector3, v2: MockVector3) { return this; }
  }
  class MockColor { constructor() {} }
  class MockScene {}
  class MockGroup {}
  class MockMesh {}
  class MockPerspectiveCamera {}
  class MockWebGLRenderer {}
  return {
    Vector3: MockVector3,
    Quaternion: MockQuaternion,
    Color: MockColor,
    Scene: MockScene,
    Group: MockGroup,
    Mesh: MockMesh,
    PerspectiveCamera: MockPerspectiveCamera,
    WebGLRenderer: MockWebGLRenderer,
    // materials/geometries - simple stubs
    MeshStandardMaterial: function() { return {}; },
    MeshPhysicalMaterial: function() { return {}; },
    SphereGeometry: function() { return {}; },
    CylinderGeometry: function() { return {}; }
  };
});

// -------------------- Utility: mockStoreFactory --------------------
/**
 * Creates a fake Zustand-like store object and exposes helpers to auto-mock module paths.
 *
 * The fake store is shaped to satisfy common test usage:
 *  - store.getState() -> returns the store object
 *  - store.setState(...) -> updates store shallowly
 *  - named methods are jest/vitest mocks (vi.fn())
 *
 * Use mockPaths([...]) to attempt auto-mocking of typical store import paths.
 */

type AnyStore = Record<string, any> & {
  getState: () => any;
  setState: (patch: any) => void;
};

function makeFakeStore(template: Record<string, any>): AnyStore {
  const store: AnyStore = { ...template } as AnyStore;

  // Attach getState
  store.getState = () => store;

  // Attach shallow setState
  store.setState = (patch: any) => {
    if (typeof patch === 'function') {
      const p = patch(store.getState());
      Object.assign(store, p);
    } else {
      Object.assign(store, patch);
    }
  };

  // Ensure functions are vi.fn()
  Object.keys(store).forEach((k) => {
    if (typeof store[k] === 'function' && !store[k]._isMockFunction) {
      store[k] = vi.fn(store[k]);
    }
  });

  return store;
}

// Note: mockPaths function removed - using direct vi.mock calls instead for better reliability

// -------------------- Auto-mock common stores used in BIOSYNTH tests --------------------
// Create store templates
const moleculeTemplate: any = {
  currentMolecule: null,
  atoms: [],
  bonds: [],
  autoBond: true,
  selectedAtomId: null,
  selectedBondId: null,
  tool: 'select' as const,
  loadingState: 'idle' as const,
  backendPredictions: null,
  error: null,
  atomToAdd: null,
  currentBondOrder: 1,
  addAtom: vi.fn(),
  removeAtom: vi.fn(),
  addBond: vi.fn(),
  removeBond: vi.fn(),
  selectAtom: vi.fn(),
  selectBond: vi.fn(),
  setMolecule: vi.fn((molecule: any) => {
    moleculeTemplate.currentMolecule = molecule;
  }),
  reset: vi.fn(() => {
    moleculeTemplate.currentMolecule = null;
    moleculeTemplate.selectedAtomId = null;
    moleculeTemplate.selectedBondId = null;
    moleculeTemplate.atoms = [];
    moleculeTemplate.bonds = [];
  }),
  setTool: vi.fn(),
  setAtomToAdd: vi.fn(),
  setBondOrder: vi.fn(),
};

const historyTemplate: any = {
  undoStack: {
    push: vi.fn(),
    undo: vi.fn(),
    redo: vi.fn(),
    clear: vi.fn(),
    canUndo: false,
    canRedo: false,
  },
  canUndo: false,
  canRedo: false,
};

const profileTemplate: any = {
  name: 'Researcher',
  avatarUrl: null,
  setName: vi.fn(),
  setAvatarUrl: vi.fn(),
  loadFromStorage: vi.fn(),
};

// create fake stores
const fakeMoleculeStore = makeFakeStore(moleculeTemplate);
const fakeHistoryStore = makeFakeStore(historyTemplate);
const fakeProfileStore = makeFakeStore(profileTemplate);

// Don't mock stores globally - let tests use real stores
// Tests that need mocks can create their own in their test files
// This allows logic tests (engineAdapter, historyStore) to use real implementations
// while component tests can mock if needed

// Provide globals for manual access if needed
(globalThis as any).__FAKE_STORES__ = {
  moleculeStore: fakeMoleculeStore,
  historyStore: fakeHistoryStore,
  profileStore: fakeProfileStore,
};

// -------------------- Test render helper --------------------
/**
 * renderWithProviders - helpful wrapper for tests that need routing + suspense
 * Usage:
 *   import { renderWithProviders } from './tests/test-utils'
 *   renderWithProviders(<MyComponent />)
 * 
 * NOTE: Moved to test-utils.tsx to avoid importing react-dom in setup file
 */

