import { describe, it, expect } from 'vitest';
import { paginate } from '../pagination';

describe('paginate', () => {
	it('paginates items', () => {
		const items = Array.from({ length: 25 }, (_, i) => i + 1);
		const { data, totalPages } = paginate(items, 2, 10);
		expect(totalPages).toBe(3);
		expect(data).toEqual([11, 12, 13, 14, 15, 16, 17, 18, 19, 20]);
	});

	it('clamps page within range', () => {
		const items = Array.from({ length: 5 }, (_, i) => i + 1);
		const { data, totalPages } = paginate(items, 10, 10);
		expect(totalPages).toBe(1);
		expect(data).toEqual([1, 2, 3, 4, 5]);
	});
});


