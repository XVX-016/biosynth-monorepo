import { describe, it, expect, beforeEach } from 'vitest';
import { KernelSelectionManager } from '../selectionManager';

describe('KernelSelectionManager', () => {
  let sm: KernelSelectionManager;

  beforeEach(() => {
    sm = new KernelSelectionManager();
  });

  it('selects and deselects atoms', () => {
    const events: any[] = [];
    sm.on('change', (p) => events.push(p));
    
    sm.selectAtom('a1');
    expect(sm.getSelectedAtomId()).toBe('a1');
    
    sm.deselect();
    expect(sm.getSelectedAtomId()).toBeNull();
    expect(events.length).toBeGreaterThanOrEqual(2);
  });

  it('emits select event when atom is selected', () => {
    const selectEvents: any[] = [];
    sm.on('select', (id) => selectEvents.push(id));
    
    sm.selectAtom('a1');
    expect(selectEvents).toContain('a1');
  });

  it('emits deselect event when atom is deselected', () => {
    const deselectEvents: any[] = [];
    sm.on('deselect', (id) => deselectEvents.push(id));
    
    sm.selectAtom('a1');
    sm.deselect();
    expect(deselectEvents).toContain('a1');
  });

  it('does not emit events when selecting same atom twice', () => {
    const events: any[] = [];
    sm.on('select', (id) => events.push(id));
    
    sm.selectAtom('a1');
    sm.selectAtom('a1'); // Same atom
    expect(events.length).toBe(1); // Only one event
  });

  it('unsubscribes from events', () => {
    const events: any[] = [];
    const unsubscribe = sm.on('select', (id) => events.push(id));
    
    sm.selectAtom('a1');
    unsubscribe();
    sm.selectAtom('a2');
    
    expect(events).toContain('a1');
    expect(events).not.toContain('a2');
  });

  it('resets state', () => {
    sm.selectAtom('a1');
    sm.reset();
    expect(sm.getSelectedAtomId()).toBeNull();
  });
});

