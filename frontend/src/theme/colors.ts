/**
 * Ivory & Chrome Synthetic Metallic Palette
 * Futuristic synthetic-metal lab aesthetic
 */

export const colors = {
  // Primary Colors
  chrome: '#C0C5D2',
  ivory: '#F6F7F8',
  spaceGrey: '#2B2E33',
  ionBlack: '#0F1115',
  frostedGlass: 'rgba(255,255,255,0.06)',
  
  // Accents
  neonCyan: '#8BF3FF',
  plasmaTeal: '#3BC7C9',
  violetEdge: '#C6BDFE',
  
  // Gradients
  chromeGradient: 'linear-gradient(135deg, #E3E6EB, #C0C5D2)',
  ivoryGradient: 'linear-gradient(135deg, #FFFFFF, #F6F7F8)',
  plasmaNeon: 'linear-gradient(135deg, #3BC7C9, #8BF3FF)',
} as const;

export type ColorName = keyof typeof colors;

