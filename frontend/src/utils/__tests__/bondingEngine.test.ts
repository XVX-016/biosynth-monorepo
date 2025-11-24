import { describe, test, expect } from 'vitest'
import { shouldAutoBond, computeAutoBonds, maxAllowedBonds } from '../bondingEngine'
import type { Atom, Molecule } from '../../types/molecule'

describe('bondingEngine', () => {
  test('shouldAutoBond C-H within range', () => {
    const a: Atom = { id: 'a', element: 'C', position: [0, 0, 0] }
    const b: Atom = { id: 'b', element: 'H', position: [1.0, 0, 0] } // approx within typical
    const mol: Molecule = { atoms: [a, b], bonds: [] }
    expect(shouldAutoBond(a, b, mol)).toBe(true)
  })

  test('shouldAutoBond rejects atoms too far apart', () => {
    const a: Atom = { id: 'a', element: 'C', position: [0, 0, 0] }
    const b: Atom = { id: 'b', element: 'H', position: [10, 0, 0] } // too far
    const mol: Molecule = { atoms: [a, b], bonds: [] }
    expect(shouldAutoBond(a, b, mol)).toBe(false)
  })

  test('computeAutoBonds finds pair', () => {
    const a: Atom = { id: 'a', element: 'C', position: [0, 0, 0] }
    const b: Atom = { id: 'b', element: 'H', position: [1.0, 0, 0] }
    const mol: Molecule = { atoms: [a, b], bonds: [] }
    const c = computeAutoBonds(mol)
    expect(c.length).toBe(1)
    expect(c[0]).toEqual({ a: 'a', b: 'b' })
  })

  test('maxAllowedBonds returns correct valence', () => {
    const c: Atom = { id: 'c', element: 'C', position: [0, 0, 0] }
    const h: Atom = { id: 'h', element: 'H', position: [0, 0, 0] }
    expect(maxAllowedBonds(c)).toBe(4)
    expect(maxAllowedBonds(h)).toBe(1)
  })
})

