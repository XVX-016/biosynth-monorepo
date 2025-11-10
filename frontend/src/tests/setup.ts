import '@testing-library/jest-dom';
import '../test/r3f-mock';
import { afterEach, vi } from 'vitest';

// Optional: add lightweight THREE fallbacks if tests import three types directly
vi.mock('three', async (importOriginal) => {
  try {
    const actual = await importOriginal<any>();
    return { ...actual };
  } catch {
    class V3 {
      x: number; y: number; z: number;
      constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
      set(x: number, y: number, z: number) { this.x = x; this.y = y; this.z = z; return this; }
      clone() { return new V3(this.x, this.y, this.z); }
      subVectors(a: V3, b: V3) { this.x = a.x - b.x; this.y = a.y - b.y; this.z = a.z - b.z; return this; }
      length() { return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z); }
      addVectors(a: V3, b: V3) { this.x = a.x + b.x; this.y = a.y + b.y; this.z = a.z + b.z; return this; }
      multiplyScalar(s: number) { this.x*=s; this.y*=s; this.z*=s; return this; }
      normalize() { const l = this.length() || 1; return this.multiplyScalar(1/l); }
    }
    class Q {
      x=0; y=0; z=0; w=1;
      setFromUnitVectors() { return this; }
    }
    return { Vector3: V3, Quaternion: Q };
  }
});

// Reset Zustand stores (best-effort)
afterEach(() => {
  try {
    // eslint-disable-next-line @typescript-eslint/no-var-requires
    const { useStore } = require('zustand');
    useStore?.setState?.({});
  } catch {}
});

// ========== Mock Zustand store for AtomMesh ==========
vi.mock('../store/moleculeStore', () => ({
  useMoleculeStore: vi.fn((selector?: any) => {
    const base = {
      selectedAtomId: null,
      tool: 'select',
      atoms: [],
      bonds: [],
      addAtom: vi.fn(),
      removeAtom: vi.fn(),
      reset: vi.fn(),
    };
    return selector ? selector(base) : base;
  }),
}));

// ========== Mock useSyncExternalStore to prevent SSR/client renderer errors ==========
vi.mock('react', async () => {
  const actual = await vi.importActual<any>('react');
  return {
    ...actual,
    useSyncExternalStore: (subscribe: any, getSnapshot: any) => {
      try {
        return getSnapshot();
      } catch {
        return {};
      }
    },
  };
});


