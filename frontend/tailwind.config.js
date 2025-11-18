/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // MolForge Minimal Monochrome Palette
        white: {
          DEFAULT: "#FFFFFF",
        },
        offwhite: {
          DEFAULT: "#F7F7F7",
        },
        lightGrey: {
          DEFAULT: "#E5E5E5",
        },
        midGrey: {
          DEFAULT: "#B5B5B5",
        },
        darkGrey: {
          DEFAULT: "#3A3A3A",
        },
        black: {
          DEFAULT: "#0F0F0F",
        },
        // Legacy aliases for backward compatibility
        ivory: {
          DEFAULT: "#F7F7F7",
          50: "#FFFFFF",
          100: "#F7F7F7",
          200: "#E5E5E5",
          300: "#B5B5B5",
          400: "#3A3A3A",
          500: "#0F0F0F",
        },
        chrome: {
          DEFAULT: "#B5B5B5",
          50: "#F7F7F7",
          100: "#E5E5E5",
          200: "#B5B5B5",
          300: "#3A3A3A",
          400: "#0F0F0F",
        },
        spaceGrey: {
          DEFAULT: "#3A3A3A",
          50: "#B5B5B5",
          100: "#3A3A3A",
          200: "#0F0F0F",
        },
        ionBlack: {
          DEFAULT: "#0F0F0F",
          50: "#3A3A3A",
          100: "#0F0F0F",
        },
        neonCyan: {
          DEFAULT: "#B5B5B5",
        },
        plasmaTeal: {
          DEFAULT: "#3A3A3A",
        },
        violetEdge: {
          DEFAULT: "#E5E5E5",
        },
        frostedGlass: 'rgba(0, 0, 0, 0.02)',
      },
      backgroundImage: {
        'chrome-gradient': 'linear-gradient(135deg, #F7F7F7, #E5E5E5)',
        'ivory-gradient': 'linear-gradient(135deg, #FFFFFF, #F7F7F7)',
        'plasma-neon': 'linear-gradient(135deg, #3A3A3A, #0F0F0F)',
      },
      borderRadius: {
        xl: '22px',
        lg: '18px',
        xl2: '1.25rem',
      },
      boxShadow: {
        chrome: '0 4px 12px rgba(0, 0, 0, 0.06)',
        'chrome-glow': '0 4px 16px rgba(0, 0, 0, 0.08)',
        'ivory-soft': '0 4px 12px rgba(0, 0, 0, 0.06)',
        glass: '0 4px 12px rgba(0, 0, 0, 0.06)',
        neon: '0 4px 12px rgba(0, 0, 0, 0.06)',
        'neon-sm': '0 4px 12px rgba(0, 0, 0, 0.06)',
        'neon-glow': '0 4px 16px rgba(0, 0, 0, 0.08)',
        'neon-hover': '0 4px 16px rgba(0, 0, 0, 0.08)',
      },
      backdropBlur: {
        xs: "2px",
      },
    },
  },
  plugins: [],
}
