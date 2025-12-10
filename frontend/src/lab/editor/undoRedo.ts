// Undo/Redo command history (immutable snapshots)
export type Snapshot = { atoms: any[]; bonds: any[] };

export default class UndoRedo {
    private stack: Snapshot[] = [];
    private pointer = -1; // points to current snapshot
    private maxSize: number;

    constructor(maxSize = 50) { this.maxSize = maxSize }

    push(s: Snapshot) {
        // drop forward history
        if (this.pointer < this.stack.length - 1) this.stack = this.stack.slice(0, this.pointer + 1);
        this.stack.push(JSON.parse(JSON.stringify(s)));
        if (this.stack.length > this.maxSize) this.stack.shift();
        this.pointer = this.stack.length - 1;
    }

    canUndo() { return this.pointer > 0; }
    canRedo() { return this.pointer < this.stack.length - 1; }

    undo(): Snapshot | null {
        if (!this.canUndo()) return null;
        this.pointer -= 1;
        return JSON.parse(JSON.stringify(this.stack[this.pointer]));
    }

    redo(): Snapshot | null {
        if (!this.canRedo()) return null;
        this.pointer += 1;
        return JSON.parse(JSON.stringify(this.stack[this.pointer]));
    }

    current(): Snapshot | null {
        if (this.pointer === -1) return null;
        return JSON.parse(JSON.stringify(this.stack[this.pointer]));
    }
}
