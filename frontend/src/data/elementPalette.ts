// Element palette for metallic/ivory/chrome theme
export interface ElementPaletteEntry {
  base: number      // base color (hex number for Three.js)
  accent: number    // accent color for highlights
  radius: number    // atomic radius for rendering
}

const elementPalette: Record<string, ElementPaletteEntry> = {
  H: { base: 0xF6F7F8, accent: 0xFFFFFF, radius: 0.55 },  // Soft white/ivory
  C: { base: 0x9DA3AE, accent: 0x2B2E33, radius: 1.0 },    // Charcoal grey
  O: { base: 0x6B8FA3, accent: 0x4A6FA5, radius: 0.9 },    // Soft desaturated blue
  N: { base: 0x8B95A1, accent: 0x5A6B7A, radius: 0.95 },    // Colder grey
  F: { base: 0xA8B5B8, accent: 0x7A8B8F, radius: 0.9 },    // Soft beige-grey
  S: { base: 0xC4B5A0, accent: 0xB5A085, radius: 1.2 },    // Soft beige tint
  P: { base: 0xA8A39D, accent: 0x8B7A6B, radius: 1.15 },   // Warm grey
  Cl: { base: 0x9FA8B3, accent: 0x6B7A8B, radius: 1.1 },  // Light steel grey
  Br: { base: 0x8B8B8B, accent: 0x6B6B6B, radius: 1.25 },  // Medium grey
  I: { base: 0x7A7A7A, accent: 0x5A5A5A, radius: 1.4 },    // Darker grey
}

export default elementPalette

