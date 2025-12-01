/**
 * Editor2D - 2D molecule editor with Konva.js
 * 
 * Note: Requires react-konva package
 * Install: npm install konva react-konva
 */

import React, { useRef, useMemo, useCallback } from 'react';
import { Stage, Layer, Circle, Line, Text, Group } from 'react-konva';
import { useLab } from '../hooks/useLab';
import { useSelection } from '../hooks/useSelection';
import { useToolMode } from '../hooks/useToolMode';
import { useAttentionWeights } from './AttentionVisualizer';
import type { Atom, Bond } from '../engines/MoleculeStateEngine';

// Grid component
function Grid({ width, height, gridSize = 30 }: { width: number; height: number; gridSize?: number }) {
  const lines = useMemo(() => {
    const gridLines: JSX.Element[] = [];
    for (let x = 0; x < width; x += gridSize) {
      gridLines.push(
        <Line
          key={`v-${x}`}
          points={[x, 0, x, height]}
          stroke="#333"
          strokeWidth={1}
          listening={false}
        />
      );
    }
    for (let y = 0; y < height; y += gridSize) {
      gridLines.push(
        <Line
          key={`h-${y}`}
          points={[0, y, width, y]}
          stroke="#333"
          strokeWidth={1}
          listening={false}
        />
      );
    }
    return gridLines;
  }, [width, height, gridSize]);

  return <Group listening={false}>{lines}</Group>;
}

// Atom node component
function AtomNode({
  atom,
  selected,
  invalid,
  importance,
  onSelect,
}: {
  atom: Atom;
  selected: boolean;
  invalid: boolean;
  importance?: number; // Attention importance score (0-1)
  onSelect: () => void;
}) {
  // Color based on importance if available
  let fillColor = invalid ? '#ef4444' : selected ? '#3b82f6' : '#ffffff';
  if (importance !== undefined && !invalid && !selected) {
    // Heat color: blue (low) -> white -> red (high)
    const intensity = importance;
    if (intensity < 0.5) {
      const r = Math.floor(255 * intensity * 2);
      const g = Math.floor(255 * intensity * 2);
      const b = 255;
      fillColor = `rgb(${r}, ${g}, ${b})`;
    } else {
      const r = 255;
      const g = Math.floor(255 * (1 - (intensity - 0.5) * 2));
      const b = Math.floor(255 * (1 - (intensity - 0.5) * 2));
      fillColor = `rgb(${r}, ${g}, ${b})`;
    }
  }
  
  const strokeColor = invalid ? '#dc2626' : selected ? '#2563eb' : '#000000';
  const strokeWidth = selected || invalid ? 2 : importance !== undefined ? 1 + importance * 2 : 1;

  return (
    <Group>
      <Circle
        x={atom.x}
        y={atom.y}
        radius={14}
        fill={fillColor}
        stroke={strokeColor}
        strokeWidth={strokeWidth}
        onClick={onSelect}
        onTap={onSelect}
        draggable={false}
      />
      <Text
        x={atom.x}
        y={atom.y}
        text={atom.element}
        fontSize={12}
        fontFamily="Arial"
        fill="#000000"
        align="center"
        verticalAlign="middle"
        offsetX={atom.element.length * 3}
        offsetY={6}
        listening={false}
      />
      {selected && (
        <Circle
          x={atom.x}
          y={atom.y}
          radius={18}
          stroke="#3b82f6"
          strokeWidth={2}
          dash={[5, 5]}
          listening={false}
        />
      )}
    </Group>
  );
}

// Bond node component
function BondNode({
  bond,
  atom1,
  atom2,
  selected,
  invalid,
  attention,
  onSelect,
}: {
  bond: Bond;
  atom1: Atom;
  atom2: Atom;
  selected: boolean;
  invalid: boolean;
  attention?: number; // Attention weight (0-1)
  onSelect: () => void;
}) {
  // Color based on attention if available
  let strokeColor = invalid ? '#ef4444' : selected ? '#3b82f6' : '#666666';
  if (attention !== undefined && !invalid && !selected) {
    // Heat color: blue (low) -> white -> red (high)
    const intensity = attention;
    if (intensity < 0.5) {
      const r = Math.floor(255 * intensity * 2);
      const g = Math.floor(255 * intensity * 2);
      const b = 255;
      strokeColor = `rgb(${r}, ${g}, ${b})`;
    } else {
      const r = 255;
      const g = Math.floor(255 * (1 - (intensity - 0.5) * 2));
      const b = Math.floor(255 * (1 - (intensity - 0.5) * 2));
      strokeColor = `rgb(${r}, ${g}, ${b})`;
    }
  }
  
  // Stroke width based on attention
  const baseWidth = bond.order === 2 ? 3 : bond.order === 3 ? 4 : 2;
  const strokeWidth = attention !== undefined
    ? baseWidth + attention * 4 // 2-6px based on attention
    : baseWidth;

  return (
    <Line
      points={[atom1.x, atom1.y, atom2.x, atom2.y]}
      stroke={strokeColor}
      strokeWidth={strokeWidth}
      lineCap="round"
      lineJoin="round"
      onClick={onSelect}
      onTap={onSelect}
      hitStrokeWidth={10}
    />
  );
}

export default function Editor2D() {
  const stageRef = useRef<any>(null);
  const { currentMolecule, moleculeEngine, toolManager, currentElement, validationResult, attentionWeights } = useLab();
  const { selectedAtomId, selectedBondId, selectAtom, selectBond } = useSelection();
  const { activeTool } = useToolMode();

  // Get stage dimensions
  const [stageSize, setStageSize] = React.useState({ width: 800, height: 600 });

  React.useEffect(() => {
    const updateSize = () => {
      const container = stageRef.current?.container();
      if (container) {
        setStageSize({
          width: container.offsetWidth || 800,
          height: container.offsetHeight || 600,
        });
      }
    };

    updateSize();
    window.addEventListener('resize', updateSize);
    return () => window.removeEventListener('resize', updateSize);
  }, []);

  // Get invalid atom/bond IDs from validation
  const invalidAtomIds = useMemo(() => {
    if (!validationResult) return new Set<string>();
    return new Set(
      validationResult.errors
        .filter((e) => e.atomId)
        .map((e) => e.atomId!)
    );
  }, [validationResult]);

  const invalidBondIds = useMemo(() => {
    if (!validationResult) return new Set<string>();
    return new Set(
      validationResult.errors
        .filter((e) => e.bondId)
        .map((e) => e.bondId!)
    );
  }, [validationResult]);

  // Handle stage events
  const handleStageClick = useCallback(
    (e: any) => {
      if (!activeTool || !currentMolecule) return;

      const stage = e.target.getStage();
      const pointerPos = stage.getPointerPosition();

      const toolContext = {
        toMoleculeCoords: () => pointerPos,
        pickAtom: () => {
          // Find closest atom within hit radius
          let closest: { id: string; dist: number } | null = null;
          const hitRadius = 20;

          currentMolecule.atoms.forEach((atom) => {
            const dist = Math.sqrt(
              Math.pow(pointerPos.x - atom.x, 2) + Math.pow(pointerPos.y - atom.y, 2)
            );
            if (dist < hitRadius && (!closest || dist < closest.dist)) {
              closest = { id: atom.id, dist };
            }
          });

          return closest ? closest.id : null;
        },
        pickBond: () => {
          // Find closest bond
          let closest: { id: string; dist: number } | null = null;
          const hitRadius = 10;

          currentMolecule.bonds.forEach((bond) => {
            const atom1 = currentMolecule.atoms.get(bond.atoms[0]);
            const atom2 = currentMolecule.atoms.get(bond.atoms[1]);
            if (!atom1 || !atom2) return;

            // Distance from point to line segment
            const A = pointerPos.x - atom1.x;
            const B = pointerPos.y - atom1.y;
            const C = atom2.x - atom1.x;
            const D = atom2.y - atom1.y;

            const dot = A * C + B * D;
            const lenSq = C * C + D * D;
            let param = -1;
            if (lenSq !== 0) param = dot / lenSq;

            let xx, yy;
            if (param < 0) {
              xx = atom1.x;
              yy = atom1.y;
            } else if (param > 1) {
              xx = atom2.x;
              yy = atom2.y;
            } else {
              xx = atom1.x + param * C;
              yy = atom1.y + param * D;
            }

            const dx = pointerPos.x - xx;
            const dy = pointerPos.y - yy;
            const dist = Math.sqrt(dx * dx + dy * dy);

            if (dist < hitRadius && (!closest || dist < closest.dist)) {
              closest = { id: bond.id, dist };
            }
          });

          return closest ? closest.id : null;
        },
        currentElement,
        stateEngine: moleculeEngine,
      };

      activeTool.onPointerDown?.(e.evt, toolContext);
    },
    [activeTool, currentMolecule, currentElement, moleculeEngine]
  );

  const handleStageMouseMove = useCallback(
    (e: any) => {
      if (!activeTool || !currentMolecule) return;

      const stage = e.target.getStage();
      const pointerPos = stage.getPointerPosition();

      const toolContext = {
        toMoleculeCoords: () => pointerPos,
        pickAtom: () => null,
        pickBond: () => null,
        currentElement,
        stateEngine: moleculeEngine,
      };

      activeTool.onPointerMove?.(e.evt, toolContext);
    },
    [activeTool, currentMolecule, currentElement, moleculeEngine]
  );

  const handleStageMouseUp = useCallback(
    (e: any) => {
      if (!activeTool) return;

      const stage = e.target.getStage();
      const pointerPos = stage.getPointerPosition();

      const toolContext = {
        toMoleculeCoords: () => pointerPos,
        pickAtom: () => null,
        pickBond: () => null,
        currentElement,
        stateEngine: moleculeEngine,
      };

      activeTool.onPointerUp?.(e.evt, toolContext);
    },
    [activeTool, currentElement, moleculeEngine]
  );

  if (!currentMolecule) {
    return (
      <div className="flex-1 bg-neutral-800 border-r border-neutral-700 relative overflow-hidden flex items-center justify-center">
        <div className="text-neutral-500">No molecule loaded</div>
      </div>
    );
  }

  const atoms = Array.from(currentMolecule.atoms.values());
  const bonds = Array.from(currentMolecule.bonds.values());

  return (
    <div className="flex-1 bg-neutral-800 border-r border-neutral-700 relative overflow-hidden">
      <Stage
        ref={stageRef}
        width={stageSize.width}
        height={stageSize.height}
        onMouseDown={handleStageClick}
        onMouseMove={handleStageMouseMove}
        onMouseUp={handleStageMouseUp}
        onTouchStart={handleStageClick}
        onTouchMove={handleStageMouseMove}
        onTouchEnd={handleStageMouseUp}
        style={{ cursor: activeTool?.id === 'drag' ? 'move' : 'default' }}
      >
        <Layer>
          <Grid width={stageSize.width} height={stageSize.height} />
        </Layer>

        <Layer>
          {/* Render bonds first (behind atoms) */}
          {bonds.map((bond) => {
            const atom1 = currentMolecule.atoms.get(bond.atoms[0]);
            const atom2 = currentMolecule.atoms.get(bond.atoms[1]);
            if (!atom1 || !atom2) return null;

            return (
              <BondNode
                key={bond.id}
                bond={bond}
                atom1={atom1}
                atom2={atom2}
                selected={bond.id === selectedBondId}
                invalid={invalidBondIds.has(bond.id)}
                onSelect={() => selectBond(bond.id)}
              />
            );
          })}

          {/* Render atoms on top */}
          {atoms.map((atom) => (
            <AtomNode
              key={atom.id}
              atom={atom}
              selected={atom.id === selectedAtomId}
              invalid={invalidAtomIds.has(atom.id)}
              onSelect={() => selectAtom(atom.id)}
            />
          ))}
        </Layer>
      </Stage>

      <div className="absolute top-2 left-2 text-xs opacity-60 bg-neutral-900 px-2 py-1 rounded z-10">
        {activeTool?.label || 'No tool'} | {atoms.length} atoms, {bonds.length} bonds
        {validationResult && !validationResult.valid && (
          <span className="ml-2 text-red-400">
            ({validationResult.errors.length} errors)
          </span>
        )}
      </div>
    </div>
  );
}
