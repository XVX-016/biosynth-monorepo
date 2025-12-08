import React, { useState, useRef, useEffect } from "react";
import { useEditor } from "../../context/EditorContext";

/**
 * FloatingToolbar
 * - draggable container
 * - buttons wired to EditorContext dispatch
 */
export default function FloatingToolbar() {
    const { state, dispatch, undo, redo, canUndo, canRedo } = useEditor();
    const elRef = useRef<HTMLDivElement | null>(null);
    const [dragging, setDragging] = useState(false);
    const [pos, setPos] = useState({ x: 50, y: 12 });

    useEffect(() => {
        const onMove = (e: MouseEvent) => {
            if (!dragging || !elRef.current) return;
            setPos(p => ({ x: p.x + e.movementX, y: p.y + e.movementY }));
        };
        const up = () => setDragging(false);
        window.addEventListener("mousemove", onMove);
        window.addEventListener("mouseup", up);
        return () => {
            window.removeEventListener("mousemove", onMove);
            window.removeEventListener("mouseup", up);
        };
    }, [dragging]);

    const style: React.CSSProperties = {
        position: "absolute",
        left: pos.x,
        top: pos.y,
        zIndex: 70,
        background: "rgba(255,255,255,0.95)",
        padding: 10,
        borderRadius: 12,
        boxShadow: "0 6px 20px rgba(0,0,0,0.35)",
        display: "flex",
        gap: 8,
        alignItems: "center",
        cursor: dragging ? "grabbing" : "grab"
    };

    const startDrag = (e: React.MouseEvent) => {
        // Only drag if clicking on the container itself, not buttons
        if ((e.target as HTMLElement).closest("button")) return;
        setDragging(true);
    };

    return (
        <div
            ref={elRef}
            style={style}
            onMouseDown={startDrag}
        >
            <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
                <button
                    className={`lab-btn ${state.tool === "select" ? "lab-btn-active" : ""}`}
                    onClick={() => dispatch({ type: "SET_TOOL", payload: "select" })}
                >
                    Select
                </button>

                <button
                    className={`lab-btn ${state.tool === "add-atom" ? "lab-btn-active" : ""}`}
                    onClick={() => dispatch({ type: "SET_TOOL", payload: "add-atom" })}
                >
                    Add Atom
                </button>

                <button
                    className={`lab-btn ${state.tool === "add-bond" ? "lab-btn-active" : ""}`}
                    onClick={() => dispatch({ type: "SET_TOOL", payload: "add-bond" })}
                >
                    Add Bond
                </button>

                <button
                    className={`lab-btn ${state.autoBond ? "lab-btn-active" : ""}`}
                    onClick={() => dispatch({ type: "TOGGLE_AUTOBOND" })}
                >
                    AutoBond
                </button>

                <button
                    className="lab-btn"
                    onClick={undo}
                    disabled={!canUndo()}
                >
                    Undo
                </button>

                <button
                    className="lab-btn"
                    onClick={redo}
                    disabled={!canRedo()}
                >
                    Redo
                </button>

                <div style={{ width: 1, height: 28, background: "#eee", marginLeft: 8 }} />

                <button
                    className="lab-btn lab-btn-primary"
                    onClick={() => dispatch({ type: "SET_BUSY", payload: true })}
                    disabled={state.busy}
                >
                    Optimize
                </button>

                <button
                    className="lab-btn"
                    onClick={() => dispatch({ type: "SET_BUSY", payload: true })}
                    disabled={state.busy}
                >
                    Predict
                </button>

                <button
                    className="lab-btn lab-btn-outline"
                    onClick={() => dispatch({ type: "LOAD_MOLECULE", payload: { atoms: [], bonds: [] } })}
                >
                    Clear
                </button>
            </div>
        </div>
    );
}
