/**
 * Central registry for molecule templates.
 * Each template contains atoms (element + coords) and bonds (aIndex, bIndex, order).
 */

import water from "./water.json";
import methane from "./methane.json";
import ethanol from "./ethanol.json";
import benzene from "./benzene.json";

export const TEMPLATES = [
  { id: "water", name: "Water (H₂O)", data: water },
  { id: "methane", name: "Methane (CH₄)", data: methane },
  { id: "ethanol", name: "Ethanol (C₂H₅OH)", data: ethanol },
  { id: "benzene", name: "Benzene (C₆H₆)", data: benzene },
];

export default TEMPLATES;

