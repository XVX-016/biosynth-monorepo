from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any
import math, json
from backend import crud

router = APIRouter(tags=["ml"])


# --- Schemas ---
class AtomSchema(BaseModel):
    id: str
    element: str
    x: float
    y: float
    z: float = 0.0

class BondSchema(BaseModel):
    a: str
    b: str
    order: int

class MLRequest(BaseModel):
    atoms: List[AtomSchema]
    bonds: List[BondSchema]

class OptimizeResponse(BaseModel):
    correctedAtoms: List[Dict[str, Any]]
    score: float
    iterations: int

class BondOrderResponse(BaseModel):
    bondOrders: List[Dict[str, Any]]  # {a,b,order,confidence}

class AutoBondResponse(BaseModel):
    bonds: List[Dict[str,Any]]

class PropertyPredictionResponse(BaseModel):
    properties: Dict[str, float]

# --- Utilities (simple physics-ish stubs) ---
def simple_relax(atoms, bonds, steps=30, step_size=0.1):
    # VERY SIMPLE: apply a small attractive force for bonded atoms, slight repulsive for all others
    # This is just a placeholder -- replace with real ML/force-field
    pos = {a['id']: [a['x'], a['y'], a.get('z', 0.0)] for a in atoms}
    neighbors = {}
    for b in bonds:
        neighbors.setdefault(b['a'], []).append((b['b'], b['order']))
        neighbors.setdefault(b['b'], []).append((b['a'], b['order']))

    def distance(p, q):
        return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)

    iterations_run = 0
    for it in range(steps):
        iterations_run += 1
        forces = {aid: [0.0,0.0,0.0] for aid in pos.keys()}

        # bonded attraction toward ideal distance
        for b in bonds:
            a1, a2 = b['a'], b['b']
            # Safety check if atoms missing
            if a1 not in pos or a2 not in pos: continue
            
            p1, p2 = pos[a1], pos[a2]
            d = distance(p1, p2) + 1e-9
            # ideal distance approximate by element radii sum (very rough)
            # fallback constant
            ideal = 1.4  # default bond length
            # simple spring: f = k*(d - ideal)
            k = 0.05 * b.get('order',1)
            mag = k * (d - ideal)
            # direction
            dx = (p2[0]-p1[0]) / d
            dy = (p2[1]-p1[1]) / d
            dz = (p2[2]-p1[2]) / d
            # apply opposite forces
            forces[a1][0] += mag * dx
            forces[a1][1] += mag * dy
            forces[a1][2] += mag * dz
            forces[a2][0] -= mag * dx
            forces[a2][1] -= mag * dy
            forces[a2][2] -= mag * dz

        # weak repulsion for all pairs (prevent collapse)
        ids = list(pos.keys())
        for i in range(len(ids)):
            for j in range(i+1, len(ids)):
                a1, a2 = ids[i], ids[j]
                p1, p2 = pos[a1], pos[a2]
                d = distance(p1, p2) + 1e-9
                if d < 0.6:
                    rep = 0.02 * (0.6 - d)
                    dx = (p1[0]-p2[0]) / d
                    dy = (p1[1]-p2[1]) / d
                    dz = (p1[2]-p2[2]) / d
                    forces[a1][0] += rep * dx
                    forces[a1][1] += rep * dy
                    forces[a1][2] += rep * dz
                    forces[a2][0] -= rep * dx
                    forces[a2][1] -= rep * dy
                    forces[a2][2] -= rep * dz

        # integrate small step
        max_move = 0.0
        for aid, f in forces.items():
            pos[aid][0] += step_size * f[0]
            pos[aid][1] += step_size * f[1]
            pos[aid][2] += step_size * f[2]
            move = math.sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]) * step_size
            if move > max_move: max_move = move

        # convergence check
        if max_move < 1e-4:
            break

    corrected = []
    for aid, p in pos.items():
        corrected.append({"id": aid, "x": p[0], "y": p[1], "z": p[2]})
    # compute fake score (lower is better)
    score = 1.0 / (1.0 + len(atoms))
    return corrected, score, iterations_run

# --- End utilities ---


@router.post("/optimize", response_model=OptimizeResponse)
def optimize_structure(req: MLRequest):
    # production: call ML model here (GNN energy minimizer)
    corrected, score, iterations = simple_relax([a.model_dump() for a in req.atoms], [b.model_dump() for b in req.bonds], steps=80)
    return {"correctedAtoms": corrected, "score": score, "iterations": iterations}


@router.post("/bond-order", response_model=BondOrderResponse)
def predict_bond_order(req: MLRequest):
    # stub: return existing orders or 1 with low confidence for new bonds
    out = []
    for b in req.bonds:
        out.append({"a": b.a, "b": b.b, "order": b.order if b.order else 1, "confidence": 0.7})
    # optionally add predicted double bonds (very naive heuristic: if both atoms are C and distance < 1.35)
    # we don't have distances in req but we can compute from coords
    pos = {a.id: (a.x, a.y, a.z) for a in req.atoms}
    for b in req.bonds:
        # Check if atoms exist in pos map
        if b.a not in pos or b.b not in pos: continue
        p1 = pos[b.a]; p2 = pos[b.b]
        d = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
        if d < 1.35:
            # bump confidence for double
            out.append({"a": b.a, "b": b.b, "order": 2, "confidence": 0.4})
    return {"bondOrders": out}


@router.post("/autobond", response_model=AutoBondResponse)
def autobond(req: MLRequest):
    # simple deterministic autobond based on covalent radii (reuse same logic from frontend engine ideally)
    # fallback: create single bonds for pairs within 1.6 Angstroms
    pairs = []
    pos = {a.id: (a.x, a.y, a.z) for a in req.atoms}
    for i in range(len(req.atoms)):
        for j in range(i+1, len(req.atoms)):
            a = req.atoms[i]; b = req.atoms[j]
            p1, p2 = pos.get(a.id), pos.get(b.id)
            if not p1 or not p2: continue
            
            d = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
            if d < 1.6:
                pairs.append({"a": a.id, "b": b.id, "order": 1})
    return {"bonds": pairs}


@router.post("/predict-properties", response_model=PropertyPredictionResponse)
def predict_properties(req: MLRequest):
    # stub: fake properties like energy (negative), dipole magnitude (abs sum), homo-lumo gap
    # production: call chemprop / custom GNN
    energy = -0.1 * len(req.atoms)
    dipole = 0.05 * math.sqrt(sum((a.x)**2 + (a.y)**2 + (a.z)**2 for a in req.atoms))
    gap = 5.0 - 0.05 * len(req.atoms)
    return {"properties": {"energy": energy, "dipole": dipole, "homo_lumo_gap": gap}}
