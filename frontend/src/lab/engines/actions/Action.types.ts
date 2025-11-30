/**
 * Action System - Types
 * 
 * All molecule editing operations are represented as Actions.
 * Actions are reversible, serializable, and loggable.
 */

export interface IAction {
  id: string;
  timestamp: number;
  label: string;
  do(state: any): any;
  undo(state: any): any;
  toJSON(): any;
}

export interface ActionPayload {
  [key: string]: any;
}

