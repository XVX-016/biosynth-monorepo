import { describe, it, expect, beforeEach } from 'vitest';
import { useProfileStore } from '../profileStore';

describe('profileStore', () => {
	beforeEach(() => {
		// reset localStorage
		localStorage.clear();
		// reset store state
		useProfileStore.setState({ name: 'Researcher', avatarUrl: null });
	});

	it('persists name to localStorage', () => {
		useProfileStore.getState().setName('Alice');
		const saved = localStorage.getItem('biosynth.profile');
		expect(saved).toBeTruthy();
		const parsed = JSON.parse(saved as string);
		expect(parsed.name).toBe('Alice');
	});

	it('loads from localStorage', () => {
		localStorage.setItem('biosynth.profile', JSON.stringify({ name: 'Bob', avatarUrl: null }));
		useProfileStore.getState().loadFromStorage();
		expect(useProfileStore.getState().name).toBe('Bob');
	});
});


