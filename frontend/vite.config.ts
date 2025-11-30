import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  resolve: {
    dedupe: ['react', 'react-dom'],
    preserveSymlinks: true,
    alias: {
      // Resolve @biosynth/engine - try node_modules first, then direct path
      '@biosynth/engine': path.resolve(__dirname, '../node_modules/@biosynth/engine/dist/index.js'),
    },
  },
  optimizeDeps: {
    // Don't pre-bundle the local package
    exclude: ['@biosynth/engine'],
  },
  server: {
    fs: {
      // Allow serving files from parent directories
      allow: ['..'],
    },
  },
})
