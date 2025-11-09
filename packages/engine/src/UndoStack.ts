import { MoleculeGraph } from "./MoleculeGraph";

export class UndoStack {
  private stack: MoleculeGraph[] = [];
  private pointer: number = -1;
  private maxSize: number = 50;

  push(state: MoleculeGraph): void {
    // Remove any states after current pointer
    this.stack = this.stack.slice(0, this.pointer + 1);

    // Add new state
    this.stack.push(state.clone());
    this.pointer++;

    // Limit stack size
    if (this.stack.length > this.maxSize) {
      this.stack.shift();
      this.pointer--;
    }
  }

  undo(): MoleculeGraph | null {
    if (this.pointer <= 0) return null;
    this.pointer--;
    return this.stack[this.pointer].clone();
  }

  redo(): MoleculeGraph | null {
    if (this.pointer >= this.stack.length - 1) return null;
    this.pointer++;
    return this.stack[this.pointer].clone();
  }

  canUndo(): boolean {
    return this.pointer > 0;
  }

  canRedo(): boolean {
    return this.pointer < this.stack.length - 1;
  }

  clear(): void {
    this.stack = [];
    this.pointer = -1;
  }

  getCurrentState(): MoleculeGraph | null {
    if (this.pointer < 0 || this.pointer >= this.stack.length) return null;
    return this.stack[this.pointer].clone();
  }
}

