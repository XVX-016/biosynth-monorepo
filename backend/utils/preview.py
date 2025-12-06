from io import BytesIO
from PIL import Image, ImageDraw, ImageFont
import math, json

# Creates a simple SVG-like PNG preview (white background)
# Expect molecule JSON format: { atoms: [{id, element, x, y, z}], bonds: [{a,b,order}] }
# or { atoms: [{id, element, position: [x,y,z]}], bonds: [{a,b,order}] }
def render_preview_png(mol_json, size=512) -> bytes:
    atoms = mol_json.get("atoms", [])
    bonds = mol_json.get("bonds", [])

    # Normalize atom coordinates
    # Handle both {x,y,z} and {position:[x,y,z]} formats
    normalized_atoms = []
    for a in atoms:
        pos = a.get("position")
        if pos and isinstance(pos, list) and len(pos) >= 2:
            x, y = pos[0], pos[1]
        else:
            x = a.get("x", 0)
            y = a.get("y", 0)
        
        normalized_atoms.append({
            "id": a.get("id"),
            "element": a.get("element", "C"),
            "x": x,
            "y": y
        })
    
    xs = [a["x"] for a in normalized_atoms] or [0]
    ys = [a["y"] for a in normalized_atoms] or [0]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    
    # Add padding
    if maxx - minx == 0: maxx += 1
    if maxy - miny == 0: maxy += 1
    
    # map coords to image space
    def project(a):
        # Invert Y for image coords
        norm_x = (a["x"] - minx) / (maxx-minx)
        norm_y = (a["y"] - miny) / (maxy-miny)
        
        px = int(norm_x * (size*0.8) + size*0.1)
        py = int((1 - norm_y) * (size*0.8) + size*0.1)
        return px, py

    im = Image.new("RGBA", (size, size), (255,255,255,255))
    draw = ImageDraw.Draw(im)

    # draw bonds
    for bond in bonds:
        a_id = bond.get("a") or bond.get("a1") or bond.get("from")
        b_id = bond.get("b") or bond.get("a2") or bond.get("to")
        
        a = next((x for x in normalized_atoms if x["id"] == a_id), None)
        b = next((x for x in normalized_atoms if x["id"] == b_id), None)
        
        if not a or not b: continue
        ax, ay = project(a)
        bx, by = project(b)
        draw.line((ax,ay,bx,by), fill=(80,80,80,255), width=4)

    # draw atoms
    for atom in normalized_atoms:
        x,y = project(atom)
        el = atom.get("element","C")
        # color scheme minimal
        colors = {"C":(40,40,40),"H":(200,200,200),"O":(200,40,40),"N":(80,120,200),"S":(200,160,40)}
        color = colors.get(el, (120,120,120))
        r = 18
        draw.ellipse((x-r,y-r,x+r,y+r), fill=color+(255,))
        # element text
        try:
            # Basic centering
            draw.text((x-6,y-8), el, fill=(255,255,255,255))
        except Exception:
            pass

    buf = BytesIO()
    im.save(buf, format="PNG")
    buf.seek(0)
    return buf.read()
