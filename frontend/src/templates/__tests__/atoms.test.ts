import { describe, it, expect } from "vitest";
import { atomTemplates, AtomTemplateName } from "../atoms";

describe("atoms", () => {
  it("exports atom templates", () => {
    expect(atomTemplates).toBeDefined();
    expect(Object.keys(atomTemplates).length).toBeGreaterThan(0);
  });

  it("has common elements", () => {
    expect(atomTemplates.C).toBeDefined();
    expect(atomTemplates.H).toBeDefined();
    expect(atomTemplates.O).toBeDefined();
    expect(atomTemplates.N).toBeDefined();
  });

  it("each template has correct structure", () => {
    const template = atomTemplates.C;
    expect(template).toHaveProperty("element");
    expect(template).toHaveProperty("position");
    expect(Array.isArray(template.position)).toBe(true);
    expect(template.position.length).toBe(3);
  });
});

