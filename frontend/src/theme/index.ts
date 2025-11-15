/**
 * Design tokens for Ivory & Chrome theme
 */
export { colors } from './colors';

export const theme = {
  colors: {
    // Primary
    chrome: '#C0C5D2',
    ivory: '#F6F7F8',
    spaceGrey: '#2B2E33',
    ionBlack: '#0F1115',
    frostedGlass: 'rgba(255,255,255,0.06)',
    
    // Accents
    neonCyan: '#8BF3FF',
    plasmaTeal: '#3BC7C9',
    violetEdge: '#C6BDFE',
  },
  
  gradients: {
    chrome: 'linear-gradient(135deg, #E3E6EB, #C0C5D2)',
    ivory: 'linear-gradient(135deg, #FFFFFF, #F6F7F8)',
    plasmaNeon: 'linear-gradient(135deg, #3BC7C9, #8BF3FF)',
  },
  
  shadows: {
    chrome: '0 4px 6px -1px rgba(192, 197, 210, 0.1), 0 2px 4px -1px rgba(192, 197, 210, 0.06)',
    glass: '0 8px 32px 0 rgba(15, 17, 21, 0.37)',
    neon: '0 0 20px rgba(139, 243, 255, 0.3)',
  },
} as const;

