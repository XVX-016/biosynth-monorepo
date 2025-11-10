import '../../../tests/setupMocks';
import { describe, it, expect } from 'vitest';
import React, { Suspense } from 'react';
import { render } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import AtomMesh from '../AtomMesh';

describe.skip('AtomMesh (skipped in unit tests; covered via integration/E2E)', () => {
  it('renders without throwing', () => {
    const { container } = render(
      <MemoryRouter>
        <Suspense fallback={null}>
          <AtomMesh id="a1" position={[0,0,0]} element="C" />
        </Suspense>
      </MemoryRouter>
    );
    expect(container).toBeTruthy();
  });
});


