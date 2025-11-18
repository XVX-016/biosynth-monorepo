/**
 * MolForge Minimal Monochrome Palette
 * Clean, modern, scientific interface
 * Grey, white, silver, and black only
 */

export const colors = {
  // Primary Colors
  white: '#FFFFFF',
  offwhite: '#F7F7F7',
  lightGrey: '#E5E5E5',
  midGrey: '#B5B5B5',
  darkGrey: '#3A3A3A',
  black: '#0F0F0F',
  border: 'rgba(0,0,0,0.06)',
  
  // Legacy aliases for backward compatibility (mapped to new colors)
  chrome: '#B5B5B5',
  ivory: '#F7F7F7',
  spaceGrey: '#3A3A3A',
  ionBlack: '#0F0F0F',
  frostedGlass: 'rgba(0,0,0,0.02)',
  
  // Removed neon accents - using greys instead
  neonCyan: '#B5B5B5',
  plasmaTeal: '#3A3A3A',
  violetEdge: '#E5E5E5',
  
  // Gradients - minimal
  chromeGradient: 'linear-gradient(135deg, #F7F7F7, #E5E5E5)',
  ivoryGradient: 'linear-gradient(135deg, #FFFFFF, #F7F7F7)',
  plasmaNeon: 'linear-gradient(135deg, #3A3A3A, #0F0F0F)',
} as const;

export type ColorName = keyof typeof colors;

/**
 * CSS variable names for use in stylesheets
 */
export const cssVars = {
  white: 'var(--color-white)',
  offwhite: 'var(--color-offwhite)',
  lightGrey: 'var(--color-lightGrey)',
  midGrey: 'var(--color-midGrey)',
  darkGrey: 'var(--color-darkGrey)',
  black: 'var(--color-black)',
  border: 'var(--color-border)',
  // Legacy aliases
  chrome: 'var(--color-chrome)',
  ivory: 'var(--color-ivory)',
  spaceGrey: 'var(--color-spaceGrey)',
  ionBlack: 'var(--color-ionBlack)',
  neonCyan: 'var(--color-neonCyan)',
  plasmaTeal: 'var(--color-plasmaTeal)',
  violetEdge: 'var(--color-violetEdge)',
  frostedGlass: 'var(--color-frostedGlass)',
} as const;

