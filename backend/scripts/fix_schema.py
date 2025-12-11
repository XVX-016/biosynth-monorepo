import sqlite3
import os

# Point to the DB in root
DB_PATH = "biosynth.db"

def fix():
    print(f"Connecting to {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # Check columns
    cursor.execute("PRAGMA table_info(molecule)")
    cols = [row[1] for row in cursor.fetchall()]
    print("Columns:", cols)
    
    if "json_graph" not in cols:
        print("Adding json_graph")
        cursor.execute("ALTER TABLE molecule ADD COLUMN json_graph TEXT")
        
    if "thumbnail_b64" not in cols:
        print("Adding thumbnail_b64")
        cursor.execute("ALTER TABLE molecule ADD COLUMN thumbnail_b64 TEXT")

    if "formula" not in cols:
        print("Adding formula")
        cursor.execute("ALTER TABLE molecule ADD COLUMN formula TEXT")

    if "updated_at" not in cols:
        print("Adding updated_at")
        cursor.execute("ALTER TABLE molecule ADD COLUMN updated_at DATETIME")

    if "properties" not in cols:
        print("Adding properties")
        cursor.execute("ALTER TABLE molecule ADD COLUMN properties TEXT")

    conn.commit()
    conn.close()
    print("Schema fix done")

if __name__ == "__main__":
    fix()
