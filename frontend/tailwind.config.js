/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
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
        
        // Legacy aliases for gradual migration
        aluminum: {
          light: '#F6F7F8',
          DEFAULT: '#C0C5D2',
          dark: '#2B2E33',
        },
        panel: {
          DEFAULT: 'rgba(255,255,255,0.06)',
          soft: 'rgba(255,255,255,0.03)',
        },
        accent: {
          cyan: '#8BF3FF',
          teal: '#3BC7C9',
          violet: '#C6BDFE',
        },
        text: {
          primary: '#F6F7F8',
          secondary: '#C0C5D2',
          tertiary: '#2B2E33',
        },
      },
      backgroundImage: {
        'chrome-gradient': 'linear-gradient(135deg, #E3E6EB, #C0C5D2)',
        'ivory-gradient': 'linear-gradient(135deg, #FFFFFF, #F6F7F8)',
        'plasma-neon': 'linear-gradient(135deg, #3BC7C9, #8BF3FF)',
      },
      borderRadius: {
        xl: '22px',
        lg: '18px',
      },
      boxShadow: {
        chrome: '0 4px 6px -1px rgba(192, 197, 210, 0.1), 0 2px 4px -1px rgba(192, 197, 210, 0.06)',
        glass: '0 8px 32px 0 rgba(15, 17, 21, 0.37)',
        neon: '0 0 20px rgba(139, 243, 255, 0.3)',
        'neon-sm': '0 0 10px rgba(139, 243, 255, 0.2)',
      },
    },
  },
  plugins: [],
}
