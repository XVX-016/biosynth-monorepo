/**
 * Cursor Tracking Hook
 * Tracks mouse position and detects card hover with state machine
 */
import { useState, useEffect, useCallback, useRef } from 'react';

export type CursorState = 'idle' | 'moving' | 'hover_card' | 'click_action' | 'waiting_response';

interface UseCursorTrackingOptions {
  debounceMs?: number;
}

interface CursorStateMachine {
  state: CursorState;
  hoveredCardId: string | null;
  cursorPosition: { x: number; y: number } | null;
  transition: (newState: CursorState, data?: any) => void;
}

export function useCursorTracking(options: UseCursorTrackingOptions = {}) {
  const { debounceMs = 100 } = options;
  const [hoveredCardId, setHoveredCardId] = useState<string | null>(null);
  const [cursorPosition, setCursorPosition] = useState<{ x: number; y: number } | null>(null);
  const [state, setState] = useState<CursorState>('idle');
  const lastMoveTime = useRef<number>(0);
  const moveTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  const handleCardHover = useCallback((cardId: string | null) => {
    setHoveredCardId(cardId);
    if (cardId) {
      setState('hover_card');
    } else if (state === 'hover_card') {
      setState('idle');
    }
  }, [state]);

  const handleMouseMove = useCallback((e: MouseEvent) => {
    const now = Date.now();
    const newPosition = {
      x: e.clientX / window.innerWidth,
      y: e.clientY / window.innerHeight,
    };
    
    setCursorPosition(newPosition);
    
    // Update state to 'moving' if not already in a specific state
    if (state === 'idle') {
      setState('moving');
    }
    
    // Clear existing timeout
    if (moveTimeoutRef.current) {
      clearTimeout(moveTimeoutRef.current);
    }
    
    // Set timeout to return to idle after no movement
    moveTimeoutRef.current = setTimeout(() => {
      if (state === 'moving' && !hoveredCardId) {
        setState('idle');
      }
    }, 500);
    
    lastMoveTime.current = now;
  }, [state, hoveredCardId]);

  const transition = useCallback((newState: CursorState, data?: any) => {
    setState(newState);
    if (data?.cardId !== undefined) {
      setHoveredCardId(data.cardId);
    }
  }, []);

  useEffect(() => {
    window.addEventListener('mousemove', handleMouseMove);
    return () => {
      window.removeEventListener('mousemove', handleMouseMove);
      if (moveTimeoutRef.current) {
        clearTimeout(moveTimeoutRef.current);
      }
    };
  }, [handleMouseMove]);

  return {
    state,
    hoveredCardId,
    cursorPosition,
    setHoveredCardId: handleCardHover,
    transition,
  };
}

