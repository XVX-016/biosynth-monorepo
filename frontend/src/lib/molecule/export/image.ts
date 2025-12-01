/**
 * Image export functions
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Exports molecule as SVG or PNG from canvas.
 */

import type { Molecule } from '../Molecule'

/**
 * Export canvas as PNG
 */
export function exportCanvasAsPNG(canvas: HTMLCanvasElement, filename: string = 'molecule.png'): void {
  const dataUrl = canvas.toDataURL('image/png')
  const link = document.createElement('a')
  link.href = dataUrl
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
}

/**
 * Export canvas as SVG
 */
export function exportCanvasAsSVG(
  canvas: HTMLCanvasElement,
  molecule: Molecule,
  filename: string = 'molecule.svg'
): void {
  const width = canvas.width
  const height = canvas.height
  
  // Get canvas image data
  const dataUrl = canvas.toDataURL('image/png')
  
  // Create SVG with embedded image
  const svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <rect width="${width}" height="${height}" fill="#ffffff"/>
  <image href="${dataUrl}" x="0" y="0" width="${width}" height="${height}" preserveAspectRatio="xMidYMid meet"/>
</svg>`
  
  const blob = new Blob([svg], { type: 'image/svg+xml' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}

/**
 * Generate SVG from molecule structure (vector-based, not canvas)
 * This is a more accurate representation than canvas export
 */
export function generateSVGFromMolecule(
  molecule: Molecule,
  width: number = 800,
  height: number = 600,
  scale: number = 1,
  offsetX: number = 0,
  offsetY: number = 0
): string {
  const atoms = molecule.getAtoms()
  const bonds = molecule.getBonds()
  
  // Calculate bounds
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity
  atoms.forEach(atom => {
    const x = atom.position[0] * scale + offsetX + width / 2
    const y = atom.position[1] * scale + offsetY + height / 2
    minX = Math.min(minX, x - 15)
    minY = Math.min(minY, y - 15)
    maxX = Math.max(maxX, x + 15)
    maxY = Math.max(maxY, y + 15)
  })
  
  const viewBoxWidth = Math.max(maxX - minX, width)
  const viewBoxHeight = Math.max(maxY - minY, height)
  const viewBoxX = minX - 20
  const viewBoxY = minY - 20
  
  // Element colors (CPK)
  const elementColors: Record<string, string> = {
    H: '#ffffff',
    C: '#909090',
    N: '#3050f8',
    O: '#ff0d0d',
    S: '#ffff30',
    P: '#ff8000',
    F: '#90e050',
    Cl: '#1ff01f',
    Br: '#a62929',
    I: '#940094',
  }
  
  // Bond widths by order
  const bondWidths: Record<number, number> = {
    1: 2,
    2: 3,
    3: 4,
    1.5: 2.5, // Aromatic
  }
  
  // Generate SVG
  let svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="${viewBoxX} ${viewBoxY} ${viewBoxWidth} ${viewBoxHeight}">
  <rect width="${viewBoxWidth}" height="${viewBoxHeight}" x="${viewBoxX}" y="${viewBoxY}" fill="#ffffff"/>
  
  <!-- Bonds -->
  <g stroke="#000000" stroke-linecap="round">
`
  
  // Draw bonds
  bonds.forEach(bond => {
    const atom1 = atoms.find(a => a.id === bond.atom1)
    const atom2 = atoms.find(a => a.id === bond.atom2)
    if (!atom1 || !atom2) return
    
    const x1 = atom1.position[0] * scale + offsetX + width / 2
    const y1 = atom1.position[1] * scale + offsetY + height / 2
    const x2 = atom2.position[0] * scale + offsetX + width / 2
    const y2 = atom2.position[1] * scale + offsetY + height / 2
    
    const width = bondWidths[bond.order] || 2
    
    if (bond.order === 2) {
      // Double bond - draw two parallel lines
      const dx = x2 - x1
      const dy = y2 - y1
      const len = Math.sqrt(dx * dx + dy * dy)
      const perpX = (-dy / len) * 2
      const perpY = (dx / len) * 2
      
      svg += `    <line x1="${x1 + perpX}" y1="${y1 + perpY}" x2="${x2 + perpX}" y2="${y2 + perpY}" stroke-width="${width}"/>\n`
      svg += `    <line x1="${x1 - perpX}" y1="${y1 - perpY}" x2="${x2 - perpX}" y2="${y2 - perpY}" stroke-width="${width}"/>\n`
    } else if (bond.order === 3) {
      // Triple bond - draw three parallel lines
      const dx = x2 - x1
      const dy = y2 - y1
      const len = Math.sqrt(dx * dx + dy * dy)
      const perpX = (-dy / len) * 3
      const perpY = (dx / len) * 3
      
      svg += `    <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke-width="${width}"/>\n`
      svg += `    <line x1="${x1 + perpX}" y1="${y1 + perpY}" x2="${x2 + perpX}" y2="${y2 + perpY}" stroke-width="${width}"/>\n`
      svg += `    <line x1="${x1 - perpX}" y1="${y1 - perpY}" x2="${x2 - perpX}" y2="${y2 - perpY}" stroke-width="${width}"/>\n`
    } else {
      // Single or aromatic bond
      svg += `    <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke-width="${width}"/>\n`
    }
  })
  
  svg += `  </g>
  
  <!-- Atoms -->
  <g>
`
  
  // Draw atoms
  atoms.forEach(atom => {
    const x = atom.position[0] * scale + offsetX + width / 2
    const y = atom.position[1] * scale + offsetY + height / 2
    const color = elementColors[atom.element] || '#909090'
    const radius = 12
    
    svg += `    <circle cx="${x}" cy="${y}" r="${radius}" fill="${color}" stroke="#000000" stroke-width="1"/>\n`
    
    // Element label
    if (atom.element !== 'C' || atoms.length === 1) {
      svg += `    <text x="${x}" y="${y}" text-anchor="middle" dominant-baseline="central" font-family="Arial, sans-serif" font-size="10" font-weight="bold" fill="#000000">${atom.element}</text>\n`
    }
  })
  
  svg += `  </g>
</svg>`
  
  return svg
}

/**
 * Export molecule as SVG file
 */
export function exportMoleculeAsSVG(
  molecule: Molecule,
  filename: string = 'molecule.svg',
  width: number = 800,
  height: number = 600,
  scale: number = 1,
  offsetX: number = 0,
  offsetY: number = 0
): void {
  const svg = generateSVGFromMolecule(molecule, width, height, scale, offsetX, offsetY)
  const blob = new Blob([svg], { type: 'image/svg+xml' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}

