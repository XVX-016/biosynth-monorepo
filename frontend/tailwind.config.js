/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // Ivory & Chrome Palette - Full color scale with opacity support
        ivory: {
          DEFAULT: "#F6F7F8",
          50: "#FBFBFB",
          100: "#F6F7F8",
          200: "#EFEFF1",
          300: "#E7E8EA",
          400: "#D9DBDE",
          500: "#CCCFD2",
          600: "#A6A9AD",
          700: "#808388",
          800: "#575A5D",
          900: "#2E3032"
        },
        chrome: {
          DEFAULT: "#C0C5D2",
          50: "#F5F6F8",
          100: "#E3E6EB",
          200: "#D1D5DC",
          300: "#C0C5D2",
          400: "#A6ADB8",
          500: "#8C959E",
          600: "#737C84",
          700: "#595F66",
          800: "#404348",
          900: "#26272A"
        },
        // Dark UI neutrals
        spaceGrey: {
          DEFAULT: "#2B2E33",
          50: "#3A3D42",
          100: "#2B2E33",
          200: "#1F2124",
          300: "#141518",
        },
        ionBlack: {
          DEFAULT: "#0F1115",
          50: "#181A1F",
          100: "#0F1115",
          200: "#0A0B0E",
        },
        // Accents with opacity support
        neonCyan: {
          DEFAULT: "#8BF3FF",
          50: "#E5FCFF",
          100: "#CBF9FF",
          200: "#8BF3FF",
          300: "#4DE6FF",
          400: "#0FD9FF",
        },
        plasmaTeal: {
          DEFAULT: "#3BC7C9",
          50: "#E0F5F6",
          100: "#C1EBEC",
          200: "#3BC7C9",
          300: "#2A9FA1",
        },
        violetEdge: {
          DEFAULT: "#C6BDFE",
          50: "#F3F1FF",
          100: "#E7E0FF",
          200: "#C6BDFE",
          300: "#A59AFD",
        },
        // Frosted glass (for backward compatibility)
        frostedGlass: 'rgba(255, 255, 255, 0.06)',
        
        // Legacy aliases (deprecated - use direct colors)
        aluminum: {
          light: "#F6F7F8",
          DEFAULT: "#C0C5D2",
          dark: "#2B2E33",
        },
        panel: {
          DEFAULT: 'rgba(255, 255, 255, 0.06)',
          soft: 'rgba(255,255,255,0.03)',
        },
        accent: {
          cyan: "#8BF3FF",
          teal: "#3BC7C9",
          violet: "#C6BDFE",
        },
        text: {
          primary: "#F6F7F8",
          secondary: "#C0C5D2",
          tertiary: "#2B2E33",
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
        xl2: '1.25rem',
      },
      boxShadow: {
        chrome: '0 4px 6px -1px rgba(192, 197, 210, 0.1), 0 2px 4px -1px rgba(192, 197, 210, 0.06)',
        'chrome-glow': '0 0 20px rgba(192, 197, 210, 0.35), 0 0 40px rgba(192, 197, 210, 0.15)',
        'ivory-soft': '0 0 20px rgba(246, 247, 248, 0.25), 0 0 40px rgba(246, 247, 248, 0.15)',
        glass: '0 8px 32px 0 rgba(15, 17, 21, 0.37)',
        neon: '0 0 20px rgba(139, 243, 255, 0.3)',
        'neon-sm': '0 0 15px rgba(139, 243, 255, 0.4)',
        'neon-glow': '0 0 20px rgba(140, 255, 255, 0.4), 0 0 40px rgba(140, 255, 255, 0.2)',
      },
      backdropBlur: {
        xs: "2px",
      },
    },
  },
  plugins: [],
}
