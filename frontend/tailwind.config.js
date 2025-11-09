/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        aluminum: {
          light: '#F5F5F7',
          DEFAULT: '#E5E7EA',
          dark: '#C9CCD1',
        },
        panel: {
          DEFAULT: '#FFFFFF',
          dark: '#D9D9DB',
        },
        accent: {
          blue: '#4EA7FF',
          green: '#6EE787',
          red: '#FF6B6B',
        },
        text: {
          primary: '#1A1A1C',
          secondary: '#5C5C5F',
          tertiary: '#9A9A9E',
        },
      },
      borderRadius: {
        lg: '18px',
      },
      boxShadow: {
        'elev-1': '0 2px 4px rgba(0,0,0,0.06), 0 8px 20px rgba(0,0,0,0.04)',
        'inset-1': 'inset 0 2px 4px rgba(0,0,0,0.06)',
      },
    },
  },
  plugins: [],
}
