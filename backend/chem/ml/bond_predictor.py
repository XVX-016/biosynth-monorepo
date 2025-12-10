import os
from typing import List, Tuple

MODEL_PATH = os.environ.get('BOND_MODEL_PATH','./backend/chem/ml/model.pth')

class BondPredictor:
    def __init__(self):
        self.model = None
        try:
            import torch
            if os.path.exists(MODEL_PATH):
                # assume the saved object exposes a `predict` method or state_dict for a known wrapper
                # We'll attempt to load a state_dict into a lightweight wrapper if available.
                print('Loading bond model from', MODEL_PATH)
                # Attempt to import a wrapper class if present
                try:
                    from backend.chem.ml.model_wrapper import ModelWrapper
                    self.model = ModelWrapper()
                    self.model.load_state_dict(torch.load(MODEL_PATH, map_location='cpu'))
                    self.model.eval()
                    print('ModelWrapper loaded.')
                except Exception as ex:
                    print('No ModelWrapper available or failed to load wrapper:', ex)
                    # last-resort: load a generic object
                    try:
                        self.model = torch.load(MODEL_PATH, map_location='cpu')
                        print('Loaded raw torch object; expecting a callable predict.')
                    except Exception as e2:
                        print('Failed to load torch object:', e2)
                        self.model = None
            else:
                print('No bond model found at', MODEL_PATH)
        except Exception as e:
            print('Torch not available or import failed:', e)
            self.model = None

    def predict(self, atoms: List[Tuple[str,str,Tuple[float,float,float]]]):
        # Prefer model prediction if model is present and exposes predict()
        try:
            if self.model is not None:
                if hasattr(self.model, 'predict'):
                    return self.model.predict(atoms)
                # if it's a raw torch module, you would need a wrapper to featurize atoms
                print('Model present but no predict() — falling back to heuristic')
        except Exception as e:
            print('Model inference error, falling back to heuristic:', e)

        # Heuristic fallback (fast, deterministic)
        bonds = []
        thresholds = { 'H': 1.2, 'C': 1.8, 'N': 1.7, 'O': 1.6, 'S': 2.1 }
        n = len(atoms)
        for i in range(n):
            id_i, el_i, pos_i = atoms[i]
            for j in range(i+1, n):
                id_j, el_j, pos_j = atoms[j]
                dx = pos_i[0]-pos_j[0]; dy = pos_i[1]-pos_j[1]; dz = pos_i[2]-pos_j[2]
                d2 = dx*dx + dy*dy + dz*dz
                import math
                d = math.sqrt(d2)
                thr = thresholds.get(el_i,1.8) + thresholds.get(el_j,1.8)
                thr *= 0.6
                if d <= thr:
                    bond = { 'id': f'b_{id_i}_{id_j}', 'a': id_i, 'b': id_j, 'order': 1 }
                    bonds.append(bond)
        return bonds
