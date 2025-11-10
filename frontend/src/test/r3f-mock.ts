import React from 'react';
import { vi } from 'vitest';

vi.mock('@react-three/fiber', () => {
  return {
    Canvas: ({ children }: any) => React.createElement('div', null, children),
    useFrame: () => {},
    useThree: () => ({ camera: {}, gl: {} }),
  };
});

vi.mock('@react-three/drei', () => {
  return {
    OrbitControls: () => React.createElement('div'),
    Html: ({ children }: any) => React.createElement('div', null, children),
    StatsGl: () => React.createElement('div'),
  };
});


