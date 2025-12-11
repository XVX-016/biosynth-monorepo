import sqlite3
import datetime

DB_PATH = "biosynth.db"

ITEMS = [
    ('Water', 'H2O', 'O'),
    ('Benzene', 'C6H6', 'c1ccccc1'),
    ('Methane', 'CH4', 'C'),
    ('Ethanol', 'C2H6O', 'CCO'),
    ('Aspirin', 'C9H8O4', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
    ('Caffeine', 'C8H10N4O2', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'),
    ('Glucose', 'C6H12O6', 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'),
    ('Ammonia', 'NH3', 'N')
]

def seed():
    print(f"Connecting to {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    for name, formula, smiles in ITEMS:
        # Check if exists
        cursor.execute("SELECT id FROM molecule WHERE name=?", (name,))
        if cursor.fetchone():
            print(f"Skipping {name}")
            continue
            
        print(f"Inserting {name}")
        try:
            cursor.execute("""
                INSERT INTO molecule (name, formula, smiles, created_at, updated_at, owner_id, molfile, json_graph, thumbnail_b64, properties)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                name, formula, smiles, 
                datetime.datetime.utcnow(), datetime.datetime.utcnow(), 
                1, # owner_id default 1
                None, 
                "{}", # json_graph
                None, 
                "{}" # properties
            ))
        except Exception as e:
            print(f"Error inserting {name}: {e}")
    
    conn.commit()
    conn.close()
    print("Done")

if __name__ == "__main__":
    seed()
