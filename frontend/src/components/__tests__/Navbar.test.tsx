import { describe, it, expect, vi } from 'vitest';
import React from 'react';
import { renderToString } from 'react-dom/server';
// Minimal router mock to avoid hook/context issues
vi.mock('react-router-dom', () => {
	const React = require('react');
	return {
		Link: (props: any) => React.createElement('a', props),
		useLocation: () => ({ pathname: '/' }),
	};
});
import Navbar from '../Navbar';

describe('Navbar', () => {
	it('renders without crashing', () => {
		const html = renderToString(<Navbar />);
		expect(html).toContain('BioSynth AI');
	});
});


