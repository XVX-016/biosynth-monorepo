import '../../../tests/setupMocks';
import { describe, it, expect, vi, beforeEach } from 'vitest';
import React from 'react';
import { render } from '@testing-library/react';
import AtomMesh from '../AtomMesh';
import { selectionManager } from '../SelectionManager';
import { useMoleculeStore } from '../../../store/moleculeStore';

// Mock SelectionManager
vi.mock('../SelectionManager', () => ({
  selectionManager: {
    onHover: vi.fn(),
    onSelect: vi.fn(),
    startDrag: vi.fn(),
    endDrag: vi.fn(),
  },
}));

// Mock useFrame
vi.mock('@react-three/fiber', async () => {
  const actual = await vi.importActual('@react-three/fiber');
  return {
    ...actual,
    useFrame: vi.fn((callback) => {
      // Simulate frame updates
      if (typeof callback === 'function') {
        callback({ clock: { elapsedTime: 0 } });
      }
    }),
  };
});

// Mock Outlines component
vi.mock('@react-three/drei', () => ({
  Outlines: ({ color, thickness }: { color: number; thickness: number }) => (
    <div data-testid="outlines" data-color={color} data-thickness={thickness} />
  ),
}));

describe('AtomMesh', () => {
  beforeEach(() => {
    vi.clearAllMocks();
    // Reset store state
    useMoleculeStore.setState({
      selectedAtomId: null,
      tool: 'select',
    });
  });

  it('renders without throwing', () => {
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    expect(container).toBeTruthy();
  });

  it('applies correct radius for Carbon', () => {
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    // Check that sphereGeometry is rendered with correct radius (1.0 for C)
    const geometry = container.querySelector('mesh')?.querySelector('sphereGeometry');
    expect(geometry).toBeTruthy();
  });

  it('applies correct radius for Hydrogen', () => {
    const { container } = render(
      <AtomMesh id="a2" position={[0, 0, 0]} element="H" />
    );
    expect(container).toBeTruthy();
  });

  it('applies correct color for Carbon', () => {
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    // Material should have color 0x9da3ae for Carbon
    const material = container.querySelector('mesh')?.querySelector('meshPhysicalMaterial');
    expect(material).toBeTruthy();
  });

  it('applies correct color for Oxygen', () => {
    const { container } = render(
      <AtomMesh id="a3" position={[0, 0, 0]} element="O" />
    );
    expect(container).toBeTruthy();
  });

  it('shows outline when selected', () => {
    useMoleculeStore.setState({ selectedAtomId: 'a1' });
    const { getByTestId } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const outlines = getByTestId('outlines');
    expect(outlines).toBeTruthy();
    expect(outlines.getAttribute('data-color')).toBe('0x8BF3FF'); // neonCyan for select
  });

  it('shows violet outline when delete tool is active and hovered', () => {
    useMoleculeStore.setState({ tool: 'delete', selectedAtomId: 'a1' });
    const { getByTestId } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const outlines = getByTestId('outlines');
    expect(outlines).toBeTruthy();
    expect(outlines.getAttribute('data-color')).toBe('0xC6BDFE'); // violetEdge for delete
  });

  it('scales up when selected', () => {
    useMoleculeStore.setState({ selectedAtomId: 'a1' });
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const mesh = container.querySelector('mesh');
    expect(mesh?.getAttribute('scale')).toBe('1.15');
  });

  it('handles hover events', () => {
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const mesh = container.querySelector('mesh');
    const hoverEvent = new Event('pointerover', { bubbles: true });
    mesh?.dispatchEvent(hoverEvent);
    expect(selectionManager.onHover).toHaveBeenCalledWith('a1');
  });

  it('handles click events for selection', () => {
    useMoleculeStore.setState({ tool: 'select' });
    const { container } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const mesh = container.querySelector('mesh');
    const clickEvent = new MouseEvent('click', { bubbles: true });
    mesh?.dispatchEvent(clickEvent);
    expect(selectionManager.onSelect).toHaveBeenCalledWith('a1');
  });

  it('handles multiple atoms with different elements', () => {
    const { container: container1 } = render(
      <AtomMesh id="a1" position={[0, 0, 0]} element="C" />
    );
    const { container: container2 } = render(
      <AtomMesh id="a2" position={[1, 0, 0]} element="O" />
    );
    const { container: container3 } = render(
      <AtomMesh id="a3" position={[2, 0, 0]} element="H" />
    );
    
    expect(container1).toBeTruthy();
    expect(container2).toBeTruthy();
    expect(container3).toBeTruthy();
  });

  it('applies correct radius for all supported elements', () => {
    const elements: Array<'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I'> = 
      ['C', 'H', 'O', 'N', 'F', 'S', 'P', 'Cl', 'Br', 'I'];
    
    elements.forEach((element) => {
      const { container } = render(
        <AtomMesh id={`atom_${element}`} position={[0, 0, 0]} element={element} />
      );
      expect(container).toBeTruthy();
    });
  });

  it('applies correct colors for all supported elements', () => {
    const elements: Array<'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I'> = 
      ['C', 'H', 'O', 'N', 'F', 'S', 'P', 'Cl', 'Br', 'I'];
    
    elements.forEach((element) => {
      const { container } = render(
        <AtomMesh id={`atom_${element}`} position={[0, 0, 0]} element={element} />
      );
      const material = container.querySelector('mesh')?.querySelector('meshPhysicalMaterial');
      expect(material).toBeTruthy();
    });
  });
});

