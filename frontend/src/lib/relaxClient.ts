import type { Molecule } from '../types/molecule'

export async function relaxServer(mol: Molecule) {
  const r = await fetch('/api/relax', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(mol)
  })
  if (!r.ok) throw new Error('relax failed')
  return r.json()
}

