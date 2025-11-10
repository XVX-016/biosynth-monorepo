import { describe, it, expect } from 'vitest';
import React from 'react';
import { MemoryRouter } from 'react-router-dom';
import { renderToString } from 'react-dom/server';
import Navbar from '../Navbar';

describe('Navbar', () => {
	it('renders without crashing', () => {
		const html = renderToString(
			<MemoryRouter>
				<Navbar />
			</MemoryRouter>
		);
		expect(html).toContain('BioSynth AI');
	});
});


