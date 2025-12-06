export interface Molecule {
  id: string;
  name: string;
  formula: string;
  previewImage: string;   // base64 OR URL served by backend
  createdAt: string;
  updatedAt: string;
  atoms: any[];           // JSON from editor
  bonds: any[];
  isValid: boolean;
  qualityScore: number;
}
