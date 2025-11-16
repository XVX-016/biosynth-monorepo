// vitest.config.ts
import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  plugins: [react() as any], // Vitest uses a bundled Vite instance → type mismatch
  test: {
    globals: true,
    environment: "jsdom",
    setupFiles: "./src/tests/vitest.setup.ts",
    include: ["src/**/*.test.{ts,tsx}"],
    alias: {
      "@biosynth/engine": path.resolve(__dirname, "../packages/engine/src"),
      "@": path.resolve(__dirname, "./src"),
    },
  },
});
